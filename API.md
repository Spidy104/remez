# API reference

Header `include/remez.hpp`, namespace `remez`, C++20.

```cpp
const std::array bands{
    remez::Band{.lower = 0.0, .upper = 0.2 * pi, .target = 1.0, .weight = 1.0},
    remez::Band{.lower = 0.3 * pi, .upper = pi, .target = 0.0, .weight = 10.0},
};
remez::Exchange designer(remez::Specification{.taps = 63, .bands = bands});
if (designer.status() != remez::Status::ok) return;
if (designer.run().status != remez::Status::ok) return;
const std::vector<double> h(designer.taps().begin(), designer.taps().end());
```

The constructor validates and allocates; `run()` designs without allocating.

## Bands

```cpp
struct Band { double lower, upper, target, weight{1.0}; };
```

Edges are in radians per sample, in [0, π]. Bands must be ordered and must not touch, and weights must be positive. The error is `weight·(target − A(ω))`.

For a differentiator, `target` is a slope s and the error is relative: `weight·(s·ω − A(ω))/(|s|·ω)`. A zero target is a stopband in both cases.

## Specification

```cpp
struct Specification {
    std::size_t taps;
    std::span<const Band> bands;
    Symmetry symmetry{Symmetry::symmetric};
    Response response{Response::constant};
    std::size_t density{16};
};
```

`density` sets the grid spacing to π/(density·(d+1)). Raise it only if a very narrow band gives `invalid_density`. `Lowpass{taps, pass_edge, stop_edge, {pass_weight, stop_weight}}` is shorthand for an odd-length lowpass.

## Filter types

| | odd N | even N |
| --- | --- | --- |
| `symmetric` | Type I, d = (N−1)/2 | Type II, d = N/2 − 1, zero at π |
| `antisymmetric` | Type III, d = (N−3)/2, zeros at 0 and π | Type IV, d = N/2 − 1, zero at 0 |

The zeros come from the symmetry, so they limit what a given length can make:

| needs gain at | odd N | even N |
| --- | --- | --- |
| neither end (bandpass, Hilbert) | I or III | II or IV |
| 0 (lowpass) | I | II |
| π (highpass) | I | IV |
| 0 and π (bandstop, notch) | I | impossible |

A band that asks for gain at a forced zero gives `incompatible_endpoint`, except a differentiator band starting at 0.

`Symmetry::automatic` uses symmetric when it fits, otherwise antisymmetric, and always antisymmetric for differentiators. `symmetry()` reports the choice. Both have a group delay of (N−1)/2 samples; antisymmetric adds a fixed 90° phase shift, which matters if you sum outputs, as in a crossover.

## Specs in Hz and dB

```cpp
enum class Region { pass, stop };
struct DecibelBand { Region region; double lower_hz, upper_hz, decibels, gain{1.0}; };

Status convert(std::span<const DecibelBand> bands, double sample_rate, std::span<Band> out) noexcept;
double achieved_decibels(const DecibelBand& band, double max_error) noexcept;
```

`decibels` is the peak-to-peak ripple for a passband, or the attenuation for a stopband. `convert` writes one `Band` per input, and `out` must be exactly that size.

`convert` maps edges as π·(2f/fs), so Nyquist lands exactly on π. It sets each weight to 1/deviation, which makes a design meet the spec exactly when `max_error <= 1`. `achieved_decibels(band, result.max_error)` gives the ripple or attenuation that band actually reached. Only flat bands are supported, not differentiators.

## Exchange

```cpp
explicit Exchange(const Specification&, std::size_t max_points = 1'000'000);
explicit Exchange(const Lowpass&, std::size_t max_points = 1'000'000);
```

An invalid spec allocates nothing. A valid one makes seven allocations, about 56·M + 88·n + 12·N bytes, where M is the number of grid points and n = d + 2. `max_points` caps M and N. Allocation failure throws `std::bad_alloc` and leaks nothing.

| member | |
| --- | --- |
| `status()` | setup result; `invalid_specification` after a move |
| `symmetry()` | the symmetry in use |
| `run(max_iterations = 64, tolerance = 1e-6)` | designs from scratch; deterministic, allocation-free, clears earlier results first |
| `taps()` | the N coefficients, or empty unless the last run succeeded |
| `amplitude(ω)` | fitted A(ω) as `std::optional`, where H = e^{−iαω}A (times i if antisymmetric) and α = (N−1)/2; rounded taps can differ |
| `extrema()` | final alternation points, in transformed units (target D/Q, weight W·Q) |
| `grid()` | the design grid, or empty unless setup succeeded |

`Result` holds `status`, `iterations`, `ripple`, `max_error` and `flatness`. `ripple` and `flatness` describe the fitted polynomial before tap recovery. Flatness is (largest fitted error − smallest fitted peak)/largest fitted error; it controls convergence.

After recovery, `max_error` also covers the returned taps on the design grid, final reference and refined peaks, with allowances for evaluation rounding. It can exceed `ripple` substantially near the precision limit. Use this final `max_error` for the dB acceptance check. It is a sampled estimate, not a continuous-band certificate. If every target is zero, or a Type I design has one constant target everywhere, `run()` returns the exact answer immediately.

`Exchange` can't be copied. Moves are `noexcept` and leave the source empty and unusable until it is assigned a new design.

## Status

| status | from | cause |
| --- | --- | --- |
| `invalid_specification` | ctor, `run` | bad enum, differentiator not antisymmetric, moved-from object |
| `invalid_taps` | ctor | zero length, too long, antisymmetric below 2, even `Lowpass` |
| `invalid_bands` | ctor, `convert` | missing, non-finite, out of range, empty, overlapping or unordered bands; wrong `out` size |
| `invalid_weights` | ctor | weight not finite and positive |
| `invalid_target` | ctor, `convert` | non-finite target, passband gain not positive |
| `incompatible_endpoint` | ctor | gain asked for at a forced zero |
| `invalid_density` | ctor | density 0, or too few grid points |
| `grid_limit` | ctor | more than `max_points` grid points |
| `invalid_tolerance` | `run` | tolerance outside (0, 1) |
| `numerical_failure` | ctor, `run` | double precision can't resolve the problem |
| `nonconvergence` | `run` | iteration limit reached |
| `invalid_sample_rate` | `convert` | rate not finite and positive |
| `invalid_decibels` | `convert` | dB not finite and positive, or too large |

`to_string(status)` returns the name. A failed setup is sticky, so `run()` keeps returning it.

## Lifetimes and threads

The band array is read only in the constructor. `taps()`, `extrema()` and `grid()` point into the designer and stay valid until it is destroyed or moved from, or until the next `run()` overwrites them.

Const members are safe to call from several threads. `run()` needs exclusive access.

## Linking

`xmake install -o <prefix> remez` installs `include/remez.hpp` and `lib/remez.lib`. The library uses the dynamic CRT (`/MD`). With plain clang, link with `-fms-runtime-lib=dll`.

## Precision

The library works in double precision without fast-math. Convergence is checked on a dense grid with refined peaks.

Flatness reaches 1e-6 for ripples down to about 1e-9 at any length, tested up to 801 taps. Below that, Hilbert transformers stop at 1e-5. At around 3e-12, Type II to IV designs report `nonconvergence`.

The taps of a 223-tap design match a 50-digit reference to 1e-15.
