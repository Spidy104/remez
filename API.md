# remez API

One header, `include/remez.hpp`, and one static library, `remez`. It needs C++20, and everything lives in `namespace remez`.

```cpp
#include <remez.hpp>

const std::array bands{
    remez::Band{.lower = 0.0, .upper = 0.2 * pi, .target = 1.0, .weight = 1.0},
    remez::Band{.lower = 0.3 * pi, .upper = pi, .target = 0.0, .weight = 10.0},
};
remez::Exchange designer(remez::Specification{.taps = 63, .bands = bands});
if (designer.status() != remez::Status::ok) return;          // invalid specification
const remez::Result result = designer.run();                 // allocation-free
if (result.status != remez::Status::ok) return;              // not converged
const std::vector<double> h(designer.taps().begin(), designer.taps().end());
```

[examples/design.cpp](examples/design.cpp) is a complete, runnable version that also filters test tones with the result.

## Linking

The library builds with xmake:

```powershell
xmake f -m release --toolchain=msvc     # or --toolchain=clang
xmake -y remez                          # static library
xmake install -o <prefix> remez         # <prefix>/lib/remez.lib, <prefix>/include/remez.hpp
```

From another xmake project, include this directory and depend on the target:

```lua
includes("path/to/remez")
target("app")
    set_kind("binary")
    add_deps("remez")      -- also adds include/ to the include path
    add_files("main.cpp")
```

Any other build system needs the include directory `<prefix>/include`, the library `<prefix>/lib/remez.lib`, and C++20. The library uses the dynamic CRT (`/MD`), so consumers must use the same runtime. Plain `clang++` defaults to the static CRT and fails to link with a `RuntimeLibrary` mismatch; pass `-fms-runtime-lib=dll`:

```powershell
clang++ -std=c++20 -fms-runtime-lib=dll -I<prefix>/include main.cpp <prefix>/lib/remez.lib
```

The `memory_checks` configuration builds the library with `/MDd` for the debug-CRT tests.

## Types

### `Band`

| Field | Meaning |
| --- | --- |
| `lower`, `upper` | Band edges in radians per sample, with 0 ≤ lower < upper ≤ π. Bands must be in increasing order and must not touch. |
| `target` | For `Response::constant`, the desired amplitude. For `Response::differentiator`, the slope s, giving a desired A(ω) = s·ω. A target of 0 makes an ordinary stopband in both cases. |
| `weight` | Error weight; must be finite and positive. The error is `weight·(target − A)`. For a differentiator band with s ≠ 0 it is the relative error `weight·(s·ω − A)/(|s|·ω)`, which has a finite limit at ω = 0. |

### `Specification`

| Field | Default | Meaning |
| --- | --- | --- |
| `taps` | — | Filter length N ≥ 1, or N ≥ 2 when antisymmetric. |
| `bands` | — | A borrowed view of at least one `Band`, read only during construction. |
| `symmetry` | `symmetric` | Together with the parity of N, this selects the filter type (see the table below). `automatic` lets the engine choose (see "Choosing the type"). |
| `response` | `constant` | `differentiator` requires `antisymmetric`. |
| `density` | 16 | Grid points per extremum. The grid spacing is π/(density·(d+1)), where d is the degree in the table below. |

| | N odd | N even |
| --- | --- | --- |
| `symmetric` | Type I, d = (N−1)/2 | Type II, d = N/2 − 1, forces A(π) = 0 |
| `antisymmetric` | Type III, d = (N−3)/2, forces A(0) = A(π) = 0 | Type IV, d = N/2 − 1, forces A(0) = 0 |

A band with a nonzero target that touches a forced zero is rejected with `incompatible_endpoint`. The one exception is a differentiator band starting at 0.

#### Choosing the type

The forced zeros are properties of linear phase itself, not of this library, so they decide which filters a given length can represent:

| Filter needs gain at | Odd N | Even N |
| --- | --- | --- |
| neither 0 nor π (bandpass, Hilbert) | Type I or III | Type II or IV |
| 0 only (lowpass) | Type I | Type II |
| π only (highpass) | Type I | Type IV only |
| both 0 and π (bandstop, notch, multiband with both ends) | Type I | **impossible** |

With `Symmetry::automatic`, the engine chooses for you:
- A differentiator becomes antisymmetric.
- Otherwise it uses `symmetric` unless that hits a forced zero, and then `antisymmetric` if that fits the bands.
- `symmetric` is preferred whenever both fit.
- If neither fits, `status()` is `incompatible_endpoint`; for even N this happens exactly when the filter needs gain at both 0 and π. Use an odd length instead.

`Exchange::symmetry()` reports the choice. Both choices have the same constant group delay of (N−1)/2 samples. The antisymmetric one also adds a constant 90° phase offset (H = i·e^{−iαω}A). That doesn't distort the delay, but it does matter when the output is summed with other filters, as in a crossover.

Measured with `automatic` at weights 1:
- **8 taps:** lowpass −14 dB, bandpass −14 dB, highpass −15 dB (chosen antisymmetric); bandstop and notch are impossible.
- **15 taps:** lowpass −24 dB, highpass −25 dB, bandpass −21 dB, bandstop −21 dB, and a 0.02π-wide notch −13 dB.

All of them converge in 4–7 iterations.

### `Lowpass`

`Lowpass{taps, pass_edge, stop_edge, {pass_weight, stop_weight}, density}` is shorthand for a Type-I specification: target 1 on [0, pass_edge] and 0 on [stop_edge, π]. An even `taps` gives `invalid_taps`.

### `Result`

| Field | Meaning |
| --- | --- |
| `status` | `ok`, or why the run stopped. |
| `iterations` | Exchange iterations performed. |
| `ripple` | The levelled error \|δ\| of the final reference, in weighted units. |
| `max_error` | The largest weighted error found on the grid after peak refinement. |
| `flatness` | (max_error − smallest selected peak)/max_error. A run converges when flatness ≤ `tolerance`. |

### `Point`

A point of the transformed problem P(x) ≈ `desired` at x = cos `omega`, with weight `weight`. A(ω) = Q(ω)·P(cos ω), where Q is 1, cos(ω/2), sin ω or sin(ω/2) for Types I–IV. `desired` and `weight` are therefore in transformed units, not the physical ones.

### `Status` and `to_string(Status)`

`constexpr std::string_view to_string(Status)` returns the enumerator's name.

| Status | Returned by | Cause |
| --- | --- | --- |
| `ok` | all | Success. |
| `invalid_specification` | constructor, `run` | Unknown enum value, a differentiator that isn't antisymmetric, or a moved-from engine. |
| `invalid_taps` | constructor | N = 0, N > `max_points`, antisymmetric N < 2, or an even `Lowpass`. |
| `invalid_bands` | constructor | No bands, non-finite edges, edges outside [0, π], empty bands, or bands overlapping or out of order. |
| `invalid_weights` | constructor | A weight that is not finite and positive. |
| `invalid_target` | constructor | A non-finite target. |
| `incompatible_endpoint` | constructor | A nonzero target at a frequency the type forces to zero. |
| `invalid_density` | constructor | Density 0, or a grid with fewer than d + 2 points. |
| `grid_limit` | constructor | The grid would exceed `max_points`. |
| `invalid_tolerance` | `run` | `tolerance` is not in (0, 1). |
| `numerical_failure` | constructor, `run` | Non-finite transformed values, a degenerate fit, fewer than d + 2 error alternations, or peaks that cannot be told apart in binary64. |
| `nonconvergence` | `run` | `max_iterations` reached without flatness ≤ `tolerance`, including `max_iterations = 0`. |

A setup failure is sticky: every `run()` returns the constructor's status.

## `Exchange`

```cpp
explicit Exchange(const Specification& spec, std::size_t max_points = Exchange::default_max_points);
explicit Exchange(const Lowpass& spec, std::size_t max_points = Exchange::default_max_points);
```

The constructor validates first. An invalid specification sets `status()` and allocates nothing.

A valid one allocates exactly seven buffers:
- the grid (32 bytes per point);
- the reference (32·(d+2));
- the bands (48 per band);
- the interpolation nodes (56·(d+2));
- the peaks (24 per grid point);
- the samples (8·(N/2+1));
- the taps (8·N).

The grid has about density·(d+1)·(total band width)/π points. `max_points` caps both the grid and N (default 1,000,000), and every count is also capped so that byte sizes fit `std::ptrdiff_t`. Allocation failure throws `std::bad_alloc`, and the buffers already allocated are released.

| Member | Behaviour |
| --- | --- |
| `Status status() const noexcept` | Setup status, or `invalid_specification` after being moved from. |
| `Symmetry symmetry() const noexcept` | The symmetry actually used: the requested one, or the type chosen for `automatic`. It is never `automatic` after a valid setup. |
| `Result run(std::size_t max_iterations = 64, double tolerance = 1e-6) noexcept` | Designs from scratch and is deterministic: the same specification and arguments give bit-identical taps. It never allocates. It first withdraws the outputs of any previous run, so a failed run leaves nothing usable. |
| `std::span<const double> taps() const noexcept` | h[0..N−1] with exact (anti)symmetry; empty unless the last run succeeded. H(ω) = e^{−iαω}A(ω), or i·e^{−iαω}A(ω) when antisymmetric, with α = (N−1)/2. |
| `std::optional<double> amplitude(double omega) const noexcept` | The design's real amplitude A(ω). It is empty unless the last run succeeded and 0 ≤ ω ≤ π, and exactly 0 at forced zeros. |
| `std::span<const Point> extrema() const noexcept` | The final d + 2 alternation points; empty unless the last run succeeded, and empty for exactly solvable targets. |
| `std::span<const Point> grid() const noexcept` | The design grid; empty unless setup succeeded. Band edges are exact grid points, and forced zeros are left out. |

**Exactly solvable targets.** If every target is 0, or a Type-I design has one shared constant target, `run()` returns `ok` with zero ripple without iterating. The taps are then 0, or a single centre tap.

**Copying and moving.** Copying is deleted, so there is never a hidden allocation. Moves are `noexcept` and never allocate. The moved-from engine reports `invalid_specification`, returns empty views and an empty `amplitude`, and refuses to run. It can be assigned a new engine.

## Ownership and lifetimes

- **Band data:** `Specification::bands` is borrowed only for the constructor call. The engine copies what it needs, so the caller's array can change or go away afterwards.
- **Returned views:** `taps()`, `extrema()` and `grid()` return views into the engine's own buffers. A view is valid until the engine is destroyed, moved from, or `run()` is called again: a new run reuses the same storage and rewrites it. To keep results across runs, copy them, as the example does.
- **Raw pointers:** the library holds none. It owns its buffers through `std::vector`, the standard RAII owner for runtime-sized arrays.

## Thread safety

The `const` members (`status`, `taps`, `amplitude`, `extrema`, `grid`) may run concurrently with each other. `run()` mutates the engine, so it needs exclusive access. Separate engines are independent.

## Numerical contract and limits

- **Arithmetic:** IEEE-754 binary64 with no fast-math. The algorithm and its numerical reasoning are in [README.md](README.md).
- **Convergence:** flatness is measured on the grid after parabolic peak refinement inside each band. It is not a rigorous continuous-band bound.
- **Accuracy:** taps reproduce the internal amplitude to 64·N·eps, plus the barycentric evaluation's forward-error bound. Checked in 50-digit arithmetic, the 223-tap lowpass taps are within 9.4e-16 of the exact minimax interpolant.
- **Reachable flatness** depends on the ripple, not the length:

| Ripple | Typical design | Flatness reached |
| --- | --- | --- |
| ≥ 1e-8 | any length (tested up to 801 taps) | 1e-6 |
| about 5e-10 | 255-tap lowpass, bandpass, Type II, differentiator | 1e-6 |
| about 3e-10 | 255/256-tap Hilbert | 1e-5 |
| about 3e-12 | 319-tap Type I | 1e-5 |
| about 3e-12 | 320-tap Type II–IV | `nonconvergence`, no taps |

Ripples below about 1e-9 lie far beyond any deployed coefficient precision: float32 is about 6e-8 and 16-bit fixed point about 3e-5. A design that hits these limits normally means the filter is longer than its specification needs.
