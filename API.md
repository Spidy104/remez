# API reference

Everything is in `include/remez.hpp`, namespace `remez`. It needs C++20.

## A complete design

```cpp
#include <remez.hpp>

const std::array bands{
    remez::Band{.lower = 0.0, .upper = 0.2 * pi, .target = 1.0, .weight = 1.0},
    remez::Band{.lower = 0.3 * pi, .upper = pi, .target = 0.0, .weight = 10.0},
};
remez::Exchange designer(remez::Specification{.taps = 63, .bands = bands});
if (designer.status() != remez::Status::ok) return;

const remez::Result result = designer.run();
if (result.status != remez::Status::ok) return;

const std::vector<double> h(designer.taps().begin(), designer.taps().end());
```

Construction checks the spec and allocates everything. `run()` does the design and allocates nothing. If you think in Hz and dB, see [Specs in Hz and dB](#specs-in-hz-and-db).

## Linking

```powershell
xmake f -m release --toolchain=msvc
xmake -y remez
xmake install -o <prefix> remez
```

That gives you `<prefix>/include/remez.hpp` and `<prefix>/lib/remez.lib`. In another xmake project:

```lua
includes("path/to/remez")
target("app")
    set_kind("binary")
    add_deps("remez")
    add_files("main.cpp")
```

The library is built against the dynamic CRT (`/MD`), and your program has to match. With plain `clang++`, pass `-fms-runtime-lib=dll`; otherwise the linker stops with a `RuntimeLibrary` mismatch:

```powershell
clang++ -std=c++20 -fms-runtime-lib=dll -I<prefix>/include main.cpp <prefix>/lib/remez.lib
```

## Bands and specifications

Frequencies are in radians per sample, from 0 to π, where π is Nyquist.

```cpp
struct Band {
    double lower;
    double upper;
    double target;
    double weight{1.0};
};
```

`target` is the gain you want across the band. For a differentiator band, `target` is a slope s instead, and the desired response is s·ω. A target of zero is a stopband either way.

The error in a band is `weight·(target − A(ω))`. For a differentiator band it is relative, `weight·(s·ω − A(ω))/(|s|·ω)`. Bands must be in increasing order, must not touch, and must have finite edges inside [0, π] with lower < upper. Weights must be finite and positive.

```cpp
struct Specification {
    std::size_t taps;
    std::span<const Band> bands;
    Symmetry symmetry{Symmetry::symmetric};
    Response response{Response::constant};
    std::size_t density{16};
};
```

`taps` is the filter length. `density` sets the grid: the spacing is π/(density·(d+1)), with d from the table below. The default of 16 is fine for nearly everything. Raise it only if setup reports `invalid_density` for a very narrow band.

`Lowpass{taps, pass_edge, stop_edge, {pass_weight, stop_weight}, density}` is a shorthand for an odd-length lowpass with gain 1 below `pass_edge` and 0 above `stop_edge`.

## Filter types

Symmetry and the parity of the length decide the type:

| | odd length | even length |
| --- | --- | --- |
| `symmetric` | Type I, d = (N−1)/2 | Type II, d = N/2 − 1, zero at π |
| `antisymmetric` | Type III, d = (N−3)/2, zero at 0 and π | Type IV, d = N/2 − 1, zero at 0 |

Those zeros come from the symmetry itself. No design method can avoid them. They decide which filters a given length can make:

| filter needs gain at | odd length | even length |
| --- | --- | --- |
| neither end (bandpass, Hilbert) | Type I or III | Type II or IV |
| 0 only (lowpass) | Type I | Type II |
| π only (highpass) | Type I | Type IV only |
| both ends (bandstop, notch) | Type I | not possible |

A band that asks for gain where the type has a zero is rejected with `incompatible_endpoint`. The one exception is a differentiator band starting at 0.

`Symmetry::automatic` lets the library choose. It uses symmetric whenever that fits, switches to antisymmetric when only that fits, and forces antisymmetric for differentiators. `Exchange::symmetry()` tells you what it picked. The two choices have the same constant group delay of (N−1)/2 samples. The antisymmetric one also adds a fixed 90° phase shift. That doesn't matter on its own, but it does if you sum the output with other filters, as in a crossover.

With `automatic` and unit weights, 8 taps give about −14 dB for a lowpass or bandpass and −15 dB for a highpass. At 15 taps the lowpass, highpass, bandpass and bandstop all reach −21 to −25 dB, and a narrow notch at 0.5π reaches −13 dB. All of them converge in under 10 iterations.

## Specs in Hz and dB

```cpp
enum class Region { pass, stop };

struct DecibelBand {
    Region region;
    double lower_hz;
    double upper_hz;
    double decibels;
    double gain{1.0};
};

Status convert(std::span<const DecibelBand> bands, double sample_rate, std::span<Band> out) noexcept;
double achieved_decibels(const DecibelBand& band, double max_error) noexcept;
```

For a passband, `decibels` is the allowed peak-to-peak ripple. For a stopband, it is the required attenuation. `gain` only applies to passbands.

`convert` writes one `Band` per input band. `out` must be exactly as long as `bands`. Frequencies are mapped with π·(2f/fs), so a band ending at fs/2 lands exactly on π at any sample rate. Each weight is set to 1/(allowed deviation), which means a finished design meets the whole spec exactly when `result.max_error <= 1`.

`achieved_decibels(band, result.max_error)` returns what that band actually reached: passband ripple in dB, or stopband attenuation in dB. Use it to see how far off a short filter is, or to search for the shortest length that works, as the example does.

`convert` returns:
- `invalid_sample_rate` for a rate that isn't finite and positive;
- `invalid_bands` for bad edges, edges above fs/2, an unknown region, an empty list, or an `out` of the wrong size;
- `invalid_decibels` for a value that isn't finite and positive, or one so large the weight overflows;
- `invalid_target` for a passband gain that isn't finite and positive.

Band ordering is checked later, by the `Exchange` constructor. The dB helper covers flat-target bands only, not differentiators.

## The designer

```cpp
explicit Exchange(const Specification& spec, std::size_t max_points = Exchange::default_max_points);
explicit Exchange(const Lowpass& spec, std::size_t max_points = Exchange::default_max_points);
```

The constructor validates first. An invalid spec sets `status()` and allocates nothing. A valid one makes exactly seven allocations: the grid, the reference, the bands, the interpolation nodes, the peaks, the transform samples and the taps.

In bytes, with M grid points and n = d + 2, that is 32M + 32n + 48·(bands) + 56n + 24M + 8(N/2+1) + 8N. M is roughly density·(d+1)·(total band width)/π. `max_points` caps both M and the length, with a default of one million. If an allocation fails, `std::bad_alloc` propagates and nothing leaks.

| member | what it does |
| --- | --- |
| `status()` | Setup result. A moved-from object reports `invalid_specification`. |
| `symmetry()` | The symmetry in use, which tells you what `automatic` picked. |
| `run(max_iterations = 64, tolerance = 1e-6)` | Designs from scratch. It never allocates and is deterministic, so the same input gives bit-identical taps. Results from any earlier run are dropped first, so a failed run leaves nothing behind. |
| `taps()` | The N coefficients, exactly symmetric or antisymmetric. Empty unless the last run succeeded. |
| `amplitude(ω)` | The real amplitude A(ω), where H(ω) = e^{−iαω}·A(ω) (times i if antisymmetric) and α = (N−1)/2. Empty unless the last run succeeded and 0 ≤ ω ≤ π. |
| `extrema()` | The final d + 2 alternation points, in the transformed units described below. |
| `grid()` | The design grid. Empty unless setup succeeded. |

The `Point`s in `extrema()` and `grid()` hold the transformed problem: the polynomial P(cos ω) with target D/Q and weight W·Q, where A = Q·P. They are mainly useful for checking a design.

`Result` holds the following:
- `status`;
- `iterations`;
- `ripple`, the levelled error;
- `max_error`, the largest weighted error found;
- `flatness`, (max_error − smallest peak)/max_error. The run converges when flatness ≤ `tolerance`.

If every target is zero, or a Type I design has one constant target everywhere, `run()` returns the exact answer straight away: zero ripple, and all-zero taps or a single centre tap.

Copying is disabled, so a design is never duplicated behind your back. Moving is `noexcept` and allocation-free, and leaves the source empty: it exposes nothing, refuses to run, and can be assigned a new design.

## Status codes

| status | where | meaning |
| --- | --- | --- |
| `ok` | everywhere | success |
| `invalid_specification` | constructor, `run` | unknown enum value, differentiator that isn't antisymmetric, or moved-from object |
| `invalid_taps` | constructor | zero length, length above `max_points`, antisymmetric length below 2, even `Lowpass` |
| `invalid_bands` | constructor, `convert` | missing, non-finite, out-of-range, empty, overlapping or unordered bands |
| `invalid_weights` | constructor | weight not finite and positive |
| `invalid_target` | constructor, `convert` | non-finite target, or passband gain not positive |
| `incompatible_endpoint` | constructor | gain requested where the filter type has a zero |
| `invalid_density` | constructor | density 0, or too few grid points for the degree |
| `grid_limit` | constructor | grid larger than `max_points` |
| `invalid_tolerance` | `run` | tolerance outside (0, 1) |
| `numerical_failure` | constructor, `run` | double precision can't resolve the problem |
| `nonconvergence` | `run` | iteration limit reached first |
| `invalid_sample_rate` | `convert` | sample rate not finite and positive |
| `invalid_decibels` | `convert` | dB value not finite and positive, or too large |

`to_string(status)` gives the name as a `std::string_view`. A failed setup is sticky: every `run()` returns the same status.

## Lifetimes and threads

The band array you pass in is only read during construction, so you can change or free it afterwards.

`taps()`, `extrema()` and `grid()` point into the designer's own storage. They stay valid until the designer is destroyed or moved from, or until the next `run()`, which reuses and overwrites that storage. Copy the taps out if you want to keep them across runs.

The library owns its buffers with `std::vector` and holds no raw pointers.

Const members can be called from several threads at once. `run()` changes the object, so it needs exclusive access. Separate designers don't share anything.

## Precision

Everything runs in double precision without fast-math. Convergence is checked on the grid with refined peaks, not proven over the continuous band. The taps match the internal amplitude to 64·N·ε plus the evaluation's own rounding bound.

How flat the error gets depends on how small the ripple is. Down to about 1e-8 (−160 dB) you get 1e-6 at any length, tested up to 801 taps. Around 5e-10, constant-target designs still reach 1e-6, while Hilbert transformers stop at 1e-5. Around 3e-12, Type I reaches 1e-5 and Types II to IV report `nonconvergence`. These levels are far past what coefficient quantization or converters preserve: float32 is about 6e-8 and 16-bit fixed point about 3e-5.
