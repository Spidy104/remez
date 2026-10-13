# Remez

Parks–McClellan (Remez exchange) design of linear-phase FIR filters in C++20 binary64 arithmetic: Types I–IV, any number of ordered disjoint bands with constant targets and weights, and relative-error differentiators. It builds with xmake, MSVC (primary) and Clang. Follow [AGENTS.md](AGENTS.md).

```cpp
#include "remez.hpp"
using namespace remez;
constexpr double pi = 3.141592653589793;

const std::array bands{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, pi, 0, 10}};
Exchange design(Specification{63, bands});        // validates and allocates everything
const Result result = design.run();                 // max_iterations = 64, flatness tolerance = 1e-6
if (result.status == Status::ok) use(design.taps()); // std::span<const double>, h[0..62]
```

`Exchange(Lowpass{taps, pass_edge, stop_edge, {pass_weight, stop_weight}})` is the odd-length Type-I shortcut. Frequencies are radians/sample in [0, π].
- **`amplitude(w)`:** gives A(w), where H = e^{-iαw}A, times i for antisymmetric filters.
- **`extrema()`:** gives the final reference.
- **`grid()`:** gives the design grid.

Every output is empty or NaN unless the last `run()` converged.

## Method

1. **Transformed problem.** Each type writes A(w) = Q(w)P(cos w) with Q = 1, cos(w/2), sin w, sin(w/2) for Types I–IV, so every type is fitted by a polynomial P of degree d.
   - **Weighted error:** W(D − QP) = WQ(D/Q − P), so the fit uses target D/Q and weight WQ.
   - **Differentiator bands** with slope s use relative error W(sw − A)/(|s|w), so the fit uses r = Q/w, target s/r and weight Wr/|s|. This stays finite at DC.
   - **Forced zeros:** grid points where Q forces A = 0 are left out, and a nonzero target there is rejected (`incompatible_endpoint`).
2. **Grid.** Each band is sampled uniformly with `std::lerp` at spacing about π/(density·(d+1)), with exact band endpoints.
3. **Initial reference: equilibrium measure.** Alternation points follow the equilibrium measure of the bands in x = cos w (Kadec). Equispaced references behave like equispaced interpolation, so their first levelled error collapses toward rounding level at high degree.
   - **Density:** for k bands it is |q(x)|/(π√|R(x)|), where R is the product of the band endpoint factors and q has one root per gap with zero gap mass (Akhiezer).
   - **Roots:** solved by Gauss–Seidel; one sweep is exact for two bands.
   - **Quadrature:** 64-node Gauss–Chebyshev.
   - **Placement:** each band gets a share of points by mass, placed at equal steps of its own measure.
   - **No mirror symmetry:** the last band leaves out its upper endpoint. A mirror-symmetric reference for a mirror-symmetric specification, such as a Hilbert band centred on π/2, makes the levelled error of an even target vanish identically.
   - **Fallback:** uniform indices when there are more bands than points.
4. **Fit.** Barycentric weights wₖ = 1/∏2(xₖ − xⱼ), with `frexp` every 16 factors.
   - **Levelled error:** from Σ wₖP(xₖ) = 0, giving δ = Σ wₖDₖ / Σ(−1)ᵏwₖ/Wₖ and node values yₖ = Dₖ + sₖ with shift sₖ = −(−1)ᵏδ/Wₖ.
   - **Evaluation set:** P is evaluated on n − 1 nodes with weights wₖ(xₖ − xₘ), which is exactly degree d. All n nodes would add a rounding-level degree-(d+1) term.
5. **Scan, refine, select.**
   - **Error form:** the Lagrange basis sums to one, so the error is evaluated as W·Σℓₖ(x)(D(x) − Dₖ − sₖ), not W(D − P).
     - **Why:** D and P agree to the ripple, so subtracting them cancels away all but a few digits. In the error form, a node in the same constant band contributes its shift exactly.
     - **Effect:** at 255 taps the flatness floor drops from about 7e-6 to about 4e-7. Bands with varying transformed targets (Hilbert's −1/sin w) keep the old floor.
   - **Scan:** each sign run of the sampled error contributes its largest point.
   - **Refine:** interior peaks are refined by successive parabolic interpolation inside their grid bracket. Band-edge peaks stay on the edge, so a bracket never crosses a transition.
   - **Select:** n alternating peaks are kept by repeatedly removing the weakest (an interior one together with its weaker equal-sign neighbour).
   - **Convergence:** when (max − min)/max over the kept peaks is at most `tolerance`. Peaks that meet are reported as `numerical_failure`.
6. **Taps.** A direct real inverse DFT of A at w_m = 2πm/N, with integer phase reduction and exact mirroring.
   - **The problem:** samples inside transition bands come from an extrapolating interpolant, whose rounding is amplified there and leaks into the bands (3e-12 at 223 taps).
   - **The fix:** one step of iterative refinement. The taps' polynomial and P are both degree d, so interpolating the residuals at the d + 1 nodes, where the values are exact data, reproduces their difference exactly.

## Memory

- **Setup:** the constructor validates first and makes no allocation for an invalid specification. Otherwise it makes exactly seven allocations, all sized from the specification and capped so byte counts fit `ptrdiff_t`: the grid, the reference, the bands, the nodes, the peaks, the samples and the taps.
- **After setup:** `run`, `taps`, `amplitude`, `extrema` and `grid` never allocate, on success or failure.
- **Allocation failure:** throws `std::bad_alloc` from the constructor, and RAII releases what was already allocated.
- **Copying and moving:** the engine copies and moves like its vectors. A moved-from engine's `run()` returns `invalid_specification`.

## Checks

```powershell
.\tools\verify-sanitizers.ps1   # MSVC ASan, Clang ASan+UBSan, debug-CRT target, debug+release on both compilers
xmake run tests --benchmark     # Release timings
```

**What the tests check (`tools/tests.cpp`)**
- **Every converged design** gets the same checks:
  - finite taps with exact (anti)symmetry;
  - the full complex tap response against `amplitude` on 2,049 points per band, plus zero quadrature (linear phase);
  - alternation and equal ripple at the reference;
  - a dense physical error search for missed peaks, using only the taps;
  - deterministic restart, and withdrawal of every output after a failed run.
- **The seed's purpose:** on seven designs up to 255 taps, the first levelled error must be at least 0.1 of the final ripple (measured 0.61–0.94) and never above it (de la Vallée Poussin). Equispaced seeds give 1e-2 to 1e-8 of it.
- **The tap-to-amplitude tolerance** is 64·N·eps for tap rounding, plus Higham's forward-error bound for the second barycentric form (IMA J. Numer. Anal. 2004), built independently from the public reference. Only the extreme-weight design needs that bound, at under 1e-3 of it.
- **References:** 12 SciPy 1.18.1 references (density 512), which are the Type-I amplitudes and coefficients and nine general-form coefficient fixtures. A minimax design must also do no worse than SciPy at 223 and 224 taps, where SciPy stops 4–5% above the minimax ripple.
- **Coverage:** lengths 1–319 across all types, highpass, bandpass, bandstop, five-band designs, Hilbert transformers and differentiators. It includes narrow (0.02π) and edge-adjacent transitions, weight ratios of 1e6 and common weight scales of 1e±100, exact zero and constant targets, and every invalid input status.
- **Debug-CRT target:** checks the seven setup allocations and their byte bounds, zero allocations in all runtime calls, no leaks after success, failure or each forced allocation failure, and positive controls for both detectors.

**Mutation checks.** Each of these deliberate bugs makes the suite fail at a specific check:
- the plain W(D − P) error;
- removing tap refinement;
- removing peak refinement;
- dropping a node without reweighting;
- a wrong antisymmetric sign;
- the equispaced seed;
- a mirror-symmetric seed;
- removing the exact-constant path;
- one hidden `std::vector` in `run()`, caught by the debug-CRT target.

**Precision check (outside the suite).** I dumped the final reference and taps as exact hex floats and compared them with the true interpolant in 50-digit mpmath:
- **223 taps:** `amplitude` is within 2.0e-15 and the taps within 9.4e-16. The taps' true maximum weighted error is 8.776885e-9 against the ripple 8.776884e-9.
- **Extreme weights (1e-3 and 1e3, 63 taps):** `amplitude` is within 2.6e-12 and the taps within 4.0e-13.

`long double` is binary64 on this MSVC-ABI target, so such checks need mpmath.

## Performance and limits

`--benchmark`, Release, five-trial medians in ms, weights 1 and 10, density 16. Transitions are 0.2π–0.3π, or 0.24π–0.26π at 255 taps. Measured on the same day; run-to-run variation is about 5%.

| Taps | Double-double version (commit f7bd6e7), MSVC / Clang | This version, MSVC / Clang |
| --- | --- | --- |
| 15 | 0.054 / 0.067 | 0.027 / 0.036 |
| 63 | 0.37 / 0.25 | 0.23 / 0.22 |
| 127 | 1.72 / 1.07 | 0.94 / 0.95 |
| 255 | 6.25 / 4.19 | 3.58 / 3.73 |

This version needs one more iteration at 15/63/127 taps. Its ripples agree with the double-double version to about 8 digits, which is what flatness 1e-6 implies. Against the first lean commit, the error form costs 0–5% on MSVC and saves 2–15% on Clang.

Long designs at tolerance 1e-6 on Clang -O2, against the double-double version:

| Design | Double-double | This version |
| --- | --- | --- |
| 223-tap lowpass | 6.7 ms | 3.1 ms |
| 223-tap bandpass | 18.5 ms | 6.2 ms |
| 223-tap Type III differentiator | 22.9 ms | 5.4 ms |
| 255-tap lowpass | 17.8 ms | 3.7 ms |

SciPy 1.18.1 takes about 2.1 ms at 255 taps.

**Binary64 limits, reported rather than hidden.** These are the flatness levels reachable at each ripple:
- **Ripples near 5e-10 (255 taps):** lowpass, bandpass, Type II and differentiator designs reach 1e-6. Hilbert transformers reach only 1e-5 (flatness 7e-6 to 1e-5), because their transformed target varies within the band.
- **Ripples near 3e-12 (319 taps):** Type-I designs reach 1e-5 after about 50 iterations (110 ms). Type II–IV designs do not converge, and `run()` returns `nonconvergence`.

For context, float32 taps carry about 6e-8 relative error and 16-bit fixed point about 3e-5, so these ripples are far below any deployed precision.
