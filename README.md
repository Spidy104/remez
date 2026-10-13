# remez

A Parks-McClellan (Remez exchange) designer for linear-phase FIR filters, written in C++20.

It handles all four linear-phase types, any number of bands, and differentiators. You can give the spec in radians with weights, or in Hz with dB tolerances. Builds with xmake on MSVC and Clang.

## Quick start

```cpp
#include <remez.hpp>

using remez::Region;
const std::array spec{
    remez::DecibelBand{Region::pass, 0.0, 9600.0, 0.5},
    remez::DecibelBand{Region::stop, 12000.0, 24000.0, 60.0},
};
std::array<remez::Band, 2> bands{};
if (remez::convert(spec, 48000.0, bands) != remez::Status::ok) return;

remez::Exchange designer(remez::Specification{.taps = 61, .bands = bands});
const auto result = designer.run();
const bool met = result.status == remez::Status::ok && result.max_error <= 1.0;
```

That asks for 0.5 dB of passband ripple up to 9.6 kHz and 60 dB of attenuation above 12 kHz, at 48 kHz. When `met` is true, `designer.taps()` holds the 61 coefficients.

`convert` sets each band's weight to 1/tolerance, so `max_error <= 1` means every band meets its dB requirement. `achieved_decibels` tells you what you actually got. The example in [examples/design.cpp](examples/design.cpp) uses that to find the shortest length that meets a 48 kHz bandpass spec. It settles on 57 taps, with 0.495 dB ripple and 60.1 dB attenuation.

The full API, including statuses, lifetimes and threading, is in [API.md](API.md).

## Building

```powershell
xmake f -m release --toolchain=msvc     # or --toolchain=clang
xmake -y                                # library, tests, example
xmake run tests
xmake run example
```

To use it from another xmake project, add `includes("path/to/remez")` and `add_deps("remez")`. Anything else can link `remez.lib` and add `include/` to the include path. API.md covers the runtime flags.

| Path | What it is |
| --- | --- |
| `include/remez.hpp` | public header |
| `src/remez.cpp` | validation, grid, exchange loop, tap recovery, dB conversion |
| `src/seed.cpp` | initial reference |
| `examples/design.cpp` | dB/Hz example |
| `tools/tests.cpp` | tests (`tests`, and `memory_tests` with the debug CRT) |
| `tools/benchmark.cpp` | timings |
| `tools/verify-sanitizers.ps1` | full release check |

## How it works

Every type is reduced to the same problem. A linear-phase filter's amplitude is A(ω) = Q(ω)·P(cos ω), where Q is 1, cos(ω/2), sin ω or sin(ω/2) for Types I to IV. That leaves a polynomial P of degree d to fit. The weighted error W(D − QP) becomes WQ(D/Q − P), so each type is fitted with target D/Q and weight WQ. Differentiators use relative error, which stays finite at DC. Grid points where Q forces A to zero are dropped.

The initial reference matters more than anything else for long filters. The optimal alternation points crowd towards band edges, following the equilibrium measure of the bands in x = cos ω. `seed.cpp` computes that measure, solving for one root per gap, and places the starting points by it. Evenly spaced starting points, which the classic implementations use, make the first error estimate collapse to rounding noise above about 190 taps. The exchange then never recovers.

Each iteration then does five things:
1. It fits P through the reference by barycentric interpolation and solves for the levelled error δ.
2. It scans the grid for the largest error in each sign run.
3. It refines each peak with a few parabolic steps.
4. It keeps the strongest alternating set.
5. It stops once the peaks are flat to within the tolerance.

Two numerical details buy most of the accuracy. First, the error is evaluated as Σ ℓₖ(x)(D − Dₖ − sₖ) rather than as D − P(x). This avoids subtracting two numbers that agree to the ripple, and lowers the precision floor about 20 times. Second, the taps come from an inverse DFT of the amplitude plus one refinement pass. Samples taken inside transition bands are badly conditioned, and the refinement removes the error they leak into the passband.

An even length forces a zero. A symmetric even-length filter is zero at π, and an antisymmetric one is zero at 0. So an 8-tap highpass only exists as the antisymmetric Type IV, and an 8-tap bandstop or notch doesn't exist in linear phase at all. With `Symmetry::automatic` the library picks whichever type fits and reports the choice. Odd lengths have no such restriction.

## Accuracy and limits

Everything runs in plain double precision with no fast-math. Convergence is checked on a dense grid with refined peaks, not proven over the continuous band. Against a 50-digit mpmath reference, the 223-tap lowpass taps are within 1e-15 of the exact interpolant they represent.

How flat the error can get depends on the size of the ripple, not on the number of taps:

| Ripple | Examples | Flatness reached |
| --- | --- | --- |
| above 1e-8 | everything up to 801 taps in the tests | 1e-6 |
| about 5e-10 | 255-tap lowpass, bandpass, Type II, differentiator | 1e-6 |
| about 3e-10 | 255-tap Hilbert transformer | 1e-5 |
| about 3e-12 | 319-tap Type I | 1e-5 |
| about 3e-12 | 320-tap Type II to IV | doesn't converge, reports `nonconvergence` |

The last rows are around −190 to −230 dB. That is far below what float32 coefficients (about −145 dB) or any converter can carry, so in practice these designs are simply longer than their spec needs.

Not supported: arbitrary target shapes (only flat bands and differentiators), minimum-phase design, FFT acceleration, SIMD and threads.

## Testing

`tools/verify-sanitizers.ps1` builds and runs everything. It runs MSVC ASan, Clang ASan and UBSan, and the debug-CRT allocation tests, then debug and release on both compilers, including the example and the benchmark.

The tests check every converged design independently of the library's internals:
- tap symmetry;
- the complex frequency response computed from the taps;
- alternation and equal ripple at the reference;
- a dense search for missed peaks.

They compare against 12 SciPy reference designs, and check that no 223/224-tap design is worse than SciPy's. They also cover every invalid input, the dB conversion against independent measurements, short designs from 2 to 16 taps, and moves.

The debug-CRT build confirms that the constructor makes exactly seven allocations and that nothing after it allocates or leaks, including when allocation fails partway through setup.

When adding a test, break the code on purpose and make sure the test fails. The current suite catches all of the deliberate bugs tried so far: wrong signs, missing refinement steps, bad seeds, broken moves, wrong dB formulas, and a stray allocation in `run()`.

## Performance

These are release timings from `xmake run benchmark`, the median of five runs, in milliseconds, at tolerance 1e-6:

| Design | Taps | MSVC | Clang |
| --- | --- | --- | --- |
| lowpass 0.2π/0.3π | 15 | 0.023 | 0.020 |
| lowpass 0.2π/0.3π | 63 | 0.18 | 0.17 |
| lowpass 0.2π/0.3π | 127 | 0.72 | 0.70 |
| lowpass 0.24π/0.26π | 255 | 2.70 | 2.72 |
| lowpass 0.2π/0.3π | 223 | 1.97 | 1.96 |
| bandpass | 223 | 5.62 | 5.59 |
| Hilbert, Type III | 223 | 2.27 | 2.26 |
| differentiator, Type III | 223 | 3.80 | 3.68 |

SciPy takes about 2.1 ms for the 255-tap narrow lowpass. The earlier double-double version of this library (commit f7bd6e7) was 1.5 to 7 times slower than these figures.
