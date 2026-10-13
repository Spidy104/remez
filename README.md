# remez

Parks-McClellan FIR filter design in C++20. It supports all four linear-phase types, any number of bands, and differentiators, with specs in radians or in Hz and dB.

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

This asks for at most 0.5 dB of ripple up to 9.6 kHz and at least 60 dB of attenuation above 12 kHz, at a 48 kHz sample rate. When `met` is true, `designer.taps()` holds the 61 coefficients. [API.md](API.md) is the full reference.

## Build

```powershell
xmake f -m release --toolchain=msvc
xmake -y
xmake run tests
```

Clang works too (`--toolchain=clang`). Every build refreshes `compile_commands.json` in the project root for clangd. To use the library from another xmake project, add `includes("path/to/remez")` and `add_deps("remez")`.

## Examples

| Target | Shows |
| --- | --- |
| `example_lowpass` | the simplest design, in radians with weights |
| `example_bandpass_db` | a 48 kHz spec in dB, and the shortest length that meets it |
| `example_short_highpass` | 8- and 9-tap highpass and bandstop with `Symmetry::automatic` |
| `example_hilbert` | a Hilbert transformer turning a cosine into a sine |
| `example_differentiator` | a differentiator checked against the true derivative |

Run any of them with `xmake run example_lowpass`, and so on.

## How it works

Each filter type is rewritten as A(ω) = Q(ω)·P(cos ω), with Q equal to 1, cos(ω/2), sin ω or sin(ω/2), so that every design becomes a polynomial fit. Each iteration does three things:
1. It fits P through the current reference points by barycentric interpolation.
2. It finds the largest error in every sign run on a dense grid, and refines each peak.
3. It keeps the strongest alternating set.

The design stops once those peaks are equal to within the tolerance.

Three details make it robust:
- **Starting points.** The initial points follow the equilibrium measure of the bands (`src/seed.cpp`). The usual evenly spaced start stalls above roughly 190 taps.
- **Error formula.** The error is computed as Σ ℓₖ(D − Dₖ − sₖ) rather than D − P. This avoids cancelling two nearly equal numbers, so long filters converge about 20 times more tightly.
- **Tap recovery.** Taps come from an inverse DFT of the amplitude plus one refinement pass. The pass removes rounding that leaks in from samples taken inside the transition bands.

## Limits

An even length forces a zero at 0 or π, so an 8-tap bandstop or notch can't be linear phase; use an odd length. Precision runs out only for ripples below about 1e-9 (−180 dB), beyond what float32 coefficients can hold. Arbitrary target shapes, minimum-phase design, FFT and SIMD are out of scope.

## Testing

`tools/verify-sanitizers.ps1` runs everything on MSVC and Clang: the tests, ASan, UBSan, the debug-CRT allocation checks, all examples and the benchmark. The tests check each design independently from its taps and compare against SciPy. They also confirm that `run()` never allocates, and they fail on deliberately broken versions of the code.

On this machine, a 63-tap lowpass designs in 0.17 ms and a 255-tap one in 2.7 ms (`xmake run benchmark`).

## License

MIT, see [LICENSE](LICENSE).
