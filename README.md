# Remez

C++20, xmake, MSVC and Clang. Follow [AGENTS.md](AGENTS.md).

Phases 1-5 are implemented: dense grids, compensated barycentric interpolation, Remez exchange, tap recovery, and Types I-IV with disjoint bands and differentiators. The core remains `remez.hpp` and `remez.cpp`; all tests share `tools/tests.cpp`. The original odd-length `Lowpass` interface is preserved.

## Grid

`taps = 2r + 1`; frequencies are radians per sample. Require `0 < pass_edge < stop_edge < pi` and positive finite weights.

```text
step = pi / (density * (r + 1))
intervals[b] = max(1, ceil((upper[b] - lower[b]) / step))
omega[b,j] = lerp(lower[b], upper[b], j / intervals[b])
x[b,j] = cos(omega[b,j])
```

Both endpoints are exact. The transition is omitted. Passband gain is 1, stopband gain is 0. Each point retains its supplied weight. Frequency increases and cosine coordinates strictly decrease.

- `grid_size(spec, max_points)` validates and returns band counts without allocation. Defaults: density 16, point budget 1,000,000.
- `fill_grid(spec, output)` writes into pre-sized storage without allocation. Extra elements are untouched. Invalid inputs and short buffers leave storage untouched; resolution failure may leave partial output that must be discarded.
- `make_grid(spec, max_points)` returns a `Grid` with status, owned points, and band split. It makes one storage allocation after sizing. Error statuses return empty storage; allocation exceptions propagate with RAII cleanup.
- `Grid::passband()` and `stopband()` borrow read-only spans. Keep the owner alive and its storage/split unchanged while using them.

## Interpolation

For `n = r+2` ordered extrema, signed weighted error is `E = W*(D-A)`:

```text
gamma[k] = 1 / product(x[k] - x[j], j != k)
s[k] = (-1)^k
delta = sum(gamma[k]*D[k]) / sum(s[k]*gamma[k]/W[k])
C[k] = D[k] - s[k]*delta/W[k]

m = argmax(abs(gamma[k]))
beta[k] = gamma[k]*(x[k] - x[m]), k != m
A(x) = sum(beta[k]*C[k]/(x-x[k]), k != m) / sum(beta[k]/(x-x[k]), k != m)
```

All n nodes solve the ripple; the adjusted beta weights interpolate on n-1 nodes, representing degree r. Omitting the largest absolute gamma bounds the omitted-node evaluation condition by n-1 in exact arithmetic; always omitting the last node needlessly rejects some edge-adjacent designs. The remaining nodes stay ordered using allocation-free `std::rotate`. Reusing the original gamma weights after dropping a node is incorrect. See the [barycentric Remez reference](https://www.chebfun.org/publications/remez.pdf).

Weights use exponent-scaled products with two-component, FMA-compensated arithmetic. Normalization is amortized over eight bounded factors; tiny differences are scaled separately. Desired values and band weights are scaled, and every fit checks the weighted equations, including the omitted node. `std::ranges::lower_bound` locates the nearest node; `max_element` projections avoid copying nodes. Evaluation uses the ordinary compensated path when its bound satisfies the requested accuracy, and a two-component path otherwise. Unresolvable cases return `numerical_failure`; the guards are not certified error bounds. Arithmetic requires binary64 and the default round-to-nearest environment.

```cpp
remez::Interpolator polynomial(extrema.size()); // One setup allocation.
const auto ripple = polynomial.fit(extrema);    // No allocations.
const auto amplitude = polynomial.evaluate(x); // x = cos(omega), no allocations.
```

`extrema` is a span of GridPoint values with strictly decreasing x in [-1,1], finite desired values, and positive finite weights. Its omega field is unused here. The count must match the constructor and be at least 2. Fitting copies data into owned storage; inputs need not outlive the fit. Copying the interpolator is disabled; moving transfers ownership.

Both methods return `Value {status, value}`; inspect status before using value. `fit` returns signed delta, not its magnitude. A failed fit invalidates the previous polynomial. `evaluate` rejects nonfinite/out-of-domain x and unfitted/moved-from state. Workspace allocation exceptions propagate during construction; fit/evaluate are noexcept. Overflowing amplitudes return `numerical_failure`. Every evaluation starts from one binary64 barycentric sum with a rounding bound, `8*eps*(sum|t*C| + |A|*sum|t|)/|sum t|` (2x margin), and falls back to the two-component sum, which has its own `16(n+4)*eps^2` bound, when that is not good enough:

- `evaluate(x)`: binary64 if the bound is within 2^-42 of |A|.
- `estimate(x)` returns `Bounded {status, value, bound}`: binary64 if the bound is within 2^-46 of the node-value scale.
- `weighted_error(point, relative_accuracy = 2^-32, absolute_accuracy = 0)`: W*(D-A), binary64 if the bound is below |E| (exact sign) and within `max(relative_accuracy*|E|, absolute_accuracy)` (relative in (0, 1), absolute finite and >= 0, else `invalid_tolerance`); the two-component path keeps tiny ripples near unit amplitude exact. Validates x, desired and weight; omega is unused.

Bounds also include an upward-rounded absolute allowance for gradual underflow, so a subnormal result cannot acquire a false zero bound. The two-component evaluator keeps amplitudes normalized internally. Weighted errors subtract desired response and amplitude in a common scale, then combine that scale with the weight before rounding; a tiny amplitude can therefore produce an accurate normal weighted error. If final rescaling loses bits and the absolute underflow allowance cannot meet the requested relative accuracy, `weighted_error` returns `numerical_failure`; exact subnormal results remain usable.

Compensated sums use Knuth's branch-free two-sum; Neumaier's magnitude branch mispredicts on alternating barycentric terms. Against two-component references, observed binary64 errors stayed below 0.15 of their bounds.

The formerly failing 255-tap initialization now evaluates all 1,846 active-grid points. An independent Chebyshev-system solve at 80 and 110 decimal digits agreed to about 1.9e-68; maximum amplitude error fell from 0.032 to 1.11e-16, with no rejected queries. Fixed binary64 inputs and reference samples are regression-tested without any external test dependency. More extreme conditioning or weight ranges can still be rejected. Public amplitude and ripple results are rounded to double; these interpolation checks alone do not establish minimax optimality.

## Exchange

```cpp
remez::Exchange engine({63, 0.2 * std::numbers::pi_v<double>, 0.3 * std::numbers::pi_v<double>, {1, 10}});
const auto result = engine.run();
if (result.status == remez::Status::ok) {
    const auto amplitude = engine.evaluate(std::cos(omega));
    const auto extrema = engine.extrema();
}
```

Construction owns the grid, interpolator, fitted extrema, and peak workspace. `run(max_iterations=64, tolerance=1e-6)` restarts deterministically from evenly distributed active-grid indices. Every iteration fits the ripple and streams a full grid scan. Each sign run contributes its largest error, including band endpoints. This implements same-sign peak replacement without storing an error array or comparing neighboring samples across the transition.

A degree-r response has at most r-1 stationary points in the active bands plus four endpoints, so peak capacity is r+3. Excess sign runs indicate numerical trouble and return `numerical_failure`. One excess alternating peak is removed from the weaker outer endpoint. Every scan and refinement error has an exact sign and is certified to `(tolerance/256)*max(|E|, ripple)`: `max_error >= ripple` and convergence resolves errors only at `tolerance` of it, so this keeps `max_error`, flatness and peak choice within 1% of that resolution. Well-conditioned points and points near zero crossings stay in binary64. The scan also visits the previous extrema, filed under their nearest grid index, so refined peaks are not lost between iterations. Each peak is refined by Brent's method (parabolic steps with golden-section fallback) over its neighboring grid interval (one-sided at band endpoints), seeded with the recorded peak location and value. It stops at `sqrt(eps)*omega + 2^-20*bracket`, typically within about ten evaluations instead of 34. Brackets stay within a single band and only improvements replace the recorded candidate. Opposite-sign peaks must be separated by a zero; if refined peaks meet or cross, the error is unresolved in binary64 and the run returns `numerical_failure` rather than passing unordered nodes to the fit. Exactly r+2 ordered, alternating peaks seed the next fit.

Success requires `(max_error-min_selected_peak)/max_error <= tolerance` and `max_error <= (1+tolerance)*abs(fitted_delta)`. The global maximum includes every grid sample and refined candidate, including candidates discarded by exchange. Equal ripple at fitted nodes alone cannot trigger convergence. This is a sampled search with local refinement, not a rigorous continuous-band error certificate; sufficiently coarse grids can miss peaks.

`ExchangeResult` reports status, fit-attempt count, nonnegative fitted ripple, maximum searched weighted error, and relative peak spread (`flatness`). Exhausting the iteration budget returns `nonconvergence`. Invalid tolerances, grid limits and numerical failures are explicit statuses. Diagnostics on failure may describe the previous completed scan; they are not a usable design. Failed runs invalidate evaluation and return an empty extrema view. Extrema borrow fitted storage; discard the view on the next run, move, or destruction. Copying is disabled and moving transfers ownership.

The general engine uses six explicit setup allocations and zero runtime allocations in the verified paths. For M grid points, B bands, polynomial degree d, peak capacity K and workspace length W, payload is `32*M + 48*B + 72*(d+2) + 24*K + 8*W` bytes. The 72 bytes per fitted point cover a 40-byte interpolation node and a 32-byte extremum. Band records own four physical parameters and two grid indices. For ordinary two-band lowpass, buffer payload grew by 112 bytes: 96 bytes for band metadata and 16 bytes for selection scratch. The engine object and allocator overhead are additional.

STL alignment and allocator bookkeeping are additional. Invalid specifications allocate no buffers; allocation exceptions propagate with RAII cleanup.

### Chebyshev scan

After each fit, `A(x) = sum c_k T_k(x)` is recovered exactly from its values at `y_m = cos(pi*m/r)` by a DCT-I with a cached cosine table. Since `cos(pi*(r-m)*k/r) = (-1)^k*cos(pi*m*k/r)`, sample pairs are staged as sums for even k and differences for odd k in the same sample buffer. This halves the inner transform length; endpoint pairs are halved and the unpaired middle sample is included only for even k. The transform's rounding allowance includes pair addition (5*eps instead of 4*eps). Grid points are evaluated by Clenshaw's recurrence four at a time: one recurrence is latency-bound and no faster than barycentric, four interleaved are about 3.5x faster per node. Refinement and previous extrema use the barycentric path.

A rounding residual at Clenshaw step k is a perturbation of c_k, which reaches A through |T_k| <= 1, so `|dA| <= 2*eps*(sum_{k>=1}|c_k| + 5*sum|b_k| + |c_0| + 2|b_1| + |A|) + series_error`. The series error is the Lebesgue constant `(2/pi) ln(r+1) + 1` times the largest sample bound (from `estimate`), plus the DCT rounding, with a 2x margin. A grid point uses the series under the same rule as the scan (exact sign, bound within `(tolerance/256)*max(|E|, ripple)`). An iteration skips the series when the smaller band weight times its error exceeds `(tolerance/256)*ripple`, as for small ripples (127 taps, ripple 2.3e-5) where `W*(1-A)` cannot be certified from any binary64 amplitude. Observed Clenshaw errors stayed below 0.053 of their bounds over 199,106 accepted points.

The regression sweep covers ordinary, narrow-transition and edge-adjacent bands through 255 taps, independent dense band sampling, common weight scaling from 1e-200 to 1e200, moves, restarts and iteration exhaustion. Twenty-three designs were compared with SciPy 1.18.1 `remez` at density 512; maximum sampled weighted amplitude difference was about 1.5e-5 of the fitted ripple. Three reference fixtures are embedded in the tests, so SciPy is not a build dependency. The wide-transition 255-tap case (edges 0.2*pi and 0.3*pi, weights 1 and 10), which SciPy rejects, now converges in 13 iterations to ripple 6.47e-10, consistent with Kaiser's length estimate (about 249 taps for that ripple and transition). Its check is internal only: alternation, dense 8,192-point band sampling, and the tap response; there is no external reference. Sampled comparisons through `evaluate()` include a `16*eps*max(W)` rounding floor, because a binary64 amplitude near unity cannot resolve 1e-10 of this ripple. The test still accepts `numerical_failure` or `nonconvergence` for this case. High-order convergence is not guaranteed.

## Taps

Phase 4 reuses the final Chebyshev coefficients with a direct inverse-cosine fallback: odd Type-I taps, caller-owned output, no FFT dependency and zero recovery allocations. The acceptance checks are known cosine-basis coefficients, full complex FIR response and group delay, independent SciPy coefficients/responses, dense band error, and buffer/lifetime/range behavior. Existing owned buffers are unchanged.

```cpp
std::vector<double> h(spec.taps);                  // Caller-owned, sized before computation.
const auto status = engine.taps(h);                // After a successful run(); no allocations.
```

For `taps = N = 2r+1`, `A(w) = h[r] + 2*sum(h[r+k]*cos(k*w), k=1..r)` has degree r, so N equispaced samples `w_m = 2*pi*m/N` invert it exactly:

```text
h[r+k] = h[r-k] = (A(0) + 2*sum(A(w_m)*cos(w_m*k), m=1..r)) / N
```

`run()` prepares the finished taps once at convergence in the existing Chebyshev workspace, leaving its cosine table intact for the next run. It uses `h[r] = c[0]`, `h[r-k] = h[r+k] = c[k]/2` only when the final series error bound is finite and at most `8*N*eps*max(abs(samples))`; otherwise it performs direct recovery once. The output overlaps the coefficient buffer one position behind it, so ascending k overwrites only consumed coefficients. `Exchange::taps(output)` then copies the N finished taps into caller storage with `std::ranges::copy`, without recomputing or changing engine state. Failed runs invalidate the cache and moves transfer it. A recovery failure is stored and reported by `taps()`; the converged amplitude remains usable. This prepares taps even for callers using only the amplitude. `Interpolator::taps(output)` remains the standalone direct inverse-cosine recovery path. That fallback recovers the degree implied by its fitted node count. Sampling uses `estimate()` and rejects an unusable bound. The r+1 amplitude samples are staged in `output[r..2r]` and normalized by an exactly reversible power-of-two scale before summation. This keeps intermediates finite for large, representable gains; coefficients are rescaled after division by N. The lower half is written first, the centre tap last, and the upper half is mirrored, so symmetry is exact and no scratch storage is needed. Each tap advances `(cos(m*w), sin(m*w))` by angle addition instead of calling cosine for every term. A direct sine/cosine restart every 32 samples limits accumulated drift; its harmonic index advances by modular addition rather than forming `m*k`. The centre tap uses a plain compensated sum of samples. Sums retain their original order and use the same compensated accumulator as evaluation. Returns `invalid_nodes` without a successful run, `insufficient_storage` for buffers shorter than N (extra elements are untouched), and `numerical_failure` for nonfinite taps; discard output on failure. The fallback remains O(r^2) arithmetic, with two trigonometric calls per tap and two per 32-sample restart. Tests recover every cosine basis through degree 127 and compare each coefficient with an independent direct-cosine transform, check zero/subnormal/near-maximum constant gains, and reject unrepresentable samples. Filter checks include exact symmetry, finiteness, the full complex frequency response with delay r, the response against internal amplitude on 16,386 band samples (tolerance 16*N*eps), the sampled minimax bound, complete 15-tap SciPy reference coefficients, and three SciPy reference responses. Binary64 sampling and transform rounding impose an absolute response floor; tiny stopband values are checked with that floor, not an unconditional relative-error guarantee.

`--benchmark` measures preparation as part of `run()` and times retrieval separately. Local MSVC Release comparisons use identical specifications and prebuilt before/after executables, alternating for three runs. Entries are medians of those runs; each run takes a five-trial median, with 10,000 tap calls per trial.

| Taps | run() before / after (ms) | taps() before / after (ns) |
| --- | --- | --- |
| 15 | 0.04719 / 0.04745 | 697.42 / 5.84 |
| 63 | 0.68352 / 0.65973 | 17.25 / 7.14 |
| 127 | 5.97195 / 5.88464 | 33.48 / 10.89 |
| 255 | 10.54970 / 10.40693 | 65.89 / 18.48 |

The 15-tap case previously repeated direct recovery; it now pays that cost once during `run()`. The other benchmark cases previously remapped certified coefficients on each tap call. Every successful retrieval now uses an O(N) contiguous copy. Overall design time changes little: the 15-tap median is about 0.6% slower, while the others are about 1-4% faster in this comparison. A Clang profile of the 255-tap final polynomial measured expansion at 82.4 us before and 71.7 us after pairing; sampling still accounts for about 60 us, so halving transform work does not halve expansion time. System load affects these local timings.

No new owned buffer or allocation is introduced; the sample buffer becomes parity pairs, then the sample/coefficient region holds finished taps. Read-only tap calls do not mutate the engine. Accuracy tolerances are unchanged; sampling and recovered coefficients have no universal relative-error guarantee. Tests compare cached taps with independent direct recovery from final fitted extrema, include even-degree DCT midpoint cases (5, 65 and 129 taps), and check identical taps after restart, invalidation and ownership transfer.

## General designs

```cpp
const double pi = std::numbers::pi_v<double>;
const std::array bands{
    remez::Band{0, 0.2*pi, 0, 10},
    remez::Band{0.3*pi, 0.6*pi, 1, 1},
    remez::Band{0.7*pi, pi, 0, 10}
};
remez::Exchange engine(remez::Specification{31, bands});
const auto result = engine.run();
std::array<double, 31> h{};
if (result.status == remez::Status::ok) {
    const auto recovered = engine.taps(h); // Check this status independently.
    const auto value = engine.amplitude(0.4*pi);
}
```

Bands must be finite, ordered, strictly separated, and have positive width and weight. `target` is a signed constant gain with `Response::constant`, or a slope s in `D(w)=s*w` with `Response::differentiator`. The constructor copies the band parameters; the specification's span need not survive construction. Invalid targets, incompatible endpoints, empty bands and invalid response/symmetry combinations have explicit statuses. Touching bands and arbitrary response callbacks are unsupported.

Symmetry and tap-count parity infer the form. With `alpha=(N-1)/2`, physical response is `H=exp(-i*alpha*w)*A` for symmetry and `H=i*exp(-i*alpha*w)*A` for antisymmetry. A conventional Hilbert target of -i therefore uses target -1. `evaluate(x)` returns the transformed polynomial P(x); `amplitude(w)` returns the signed physical amplitude, in radians/sample. `extrema()` borrows fitted points with transformed desired values and weights. Physical band parameters remain in the owned grid metadata. Views expire on owner destruction, move or another run.

| Type | N | Polynomial degree d | A(w)=Q(w)*P(cos(w)) | Forced zeros |
| --- | --- | --- | --- | --- |
| I | odd | (N-1)/2 | Q=1 | none |
| II | even | N/2-1 | Q=cos(w/2) | pi |
| III | odd, at least 3 | (N-3)/2 | Q=sin(w) | 0, pi |
| IV | even | N/2-1 | Q=sin(w/2) | 0 |

The fit uses `Dtilde=D/Q`, `Wtilde=W*Q` and d+2 nodes. The scan and refinement measure the error of that same transformed problem, `Wtilde*(Dtilde-P)`, which equals the physical `W*(D-Q*P)` up to rounding Q and D/Q. Measuring physically would bury the tiny ripples of early high-order iterations (1e-16..1e-24) under that rounding and erase their alternation. Each refinement probe recomputes its transformed target; Type I reuses the band constant. Brackets stay within the owning band. Compatible ordinary forced endpoints have exactly zero error and are omitted from the fit. Nonzero targets at forced zeros are rejected as an API restriction. Differentiator DC uses analytic limits instead of an infinite weight. Half-angle products are scaled before rounding to avoid erasing a representable subnormal physical amplitude.

For differentiator bands with s != 0 the objective is `band.weight/abs(s) * (s - A(omega)/omega)`. At DC use `A'(0)` in place of `A(omega)/omega`. `Q(omega)/omega` tends to 1 for Type III and 1/2 for Type IV, giving transformed DC values s and 2s and transformed weights `weight/abs(s)` and `weight/(2*abs(s))`. Zero-slope stopbands use ordinary absolute-error weights. Relative differentiators require antisymmetry. SciPy reference generation converts desired slopes to `2*pi*s` and nonzero-band weights to `weight/(2*pi*abs(s))` with fs=2 and bands normalized by pi.

The transformed and physical errors differ by at most `f = 16*eps*max(abs(weighted_target))` per point (rounding Q and D/Q, with margin), so the physical peak spread can exceed the transformed one by `2*f/max_error`. On convergence Types II-IV add that to the reported flatness and return `numerical_failure` if it then exceeds `tolerance`. More iterations cannot reduce it, and loosening `tolerance` is the honest remedy. With f about 3.6e-15 and tolerance 1e-6, ripples down to about 7e-9 converge. Earlier versions required f below 1% of `tolerance*max_error`, which rejected every new-form design with ripple below about 3.6e-7, around 190 taps for the reference specifications.

Only run peaks need accuracy at `tolerance/256` of the ripple, so the scan is lazy. It ranks every point by a binary64 estimate with a rigorous bound: four-point Clenshaw while the series bound is below ripple/16, otherwise binary64 barycentric. It then certifies only two kinds of point:
- points whose sign the bound leaves uncertain;
- the two largest estimates of each sign run.

If any other point of a run could still exceed that run's certified peak, the iteration rescans with every point certified. The peaks and `max_error` are therefore those of a fully certified scan.

At 223 taps a 1e-8 ripple needs about 3e-17 absolute accuracy, while Clenshaw's bound is 2.3e-13. Before this change every grid point took the two-component path. Now about 2,000 precise evaluations per iteration drop to about 560 (I-lp 223), and the scan time falls from 13.3 to 1.1 ms. The certified rescan fired once in the 50 profiled iterations. A full peak buffer is a numerical failure; since the equilibrium seed, no measured design fills it. There is no continuous-band optimality certificate.

The initial reference follows the equilibrium measure of the bands in `x = cos(omega)`, where alternation points of best approximations concentrate (Kadec). For k bands with endpoints `e(0) > ... > e(2k-1)` its density is `|q(x)|/(pi*sqrt|R(x)|)`, where `R = prod(x - e(j))` and q is monic with one root per gap, chosen so each gap carries no mass (Akhiezer).
- **Roots of q:** solved by Gauss-Seidel. With the other roots fixed, each root is the weighted mean of x over its gap, so one sweep is exact for two bands.
- **Quadrature:** 64-node Gauss-Chebyshev in `x = c + h*cos(t)`, which absorbs the interval's endpoint singularities.
- **Points per band:** each band gets one point plus a cumulative-floor share of the rest by mass. It then places them at equal steps of its own trapezoid CDF in t, endpoints included, snapped to distinct grid indices.
- **Cost:** O((k + grid)·k) per run, with the k-1 roots held in existing scratch, so it adds no allocations.
- **Why not equispaced:** an equispaced start (the Parks-McClellan and SciPy choice) behaves like equispaced interpolation. Its first levelled error fell to 1e-11..1e-20 at 223+ taps, and every form stalled at rounding level. Filip (ACM TOMS 2016) reports the same failure and uses approximate Fekete points, which follow the same measure asymptotically.

Measured at density 16 and tolerance 1e-6 (Clang -O2; transition 0.2*pi to 0.3*pi, the 0/0.2/0.3/0.6/0.7 bandpass, Hilbert band 0.05*pi to 0.95*pi):

| Taps | Before | After |
| --- | --- | --- |
| 63 | I-lp 10 iterations, 5.6 ms | 5 iterations, 0.9 ms |
| 63 | Hilbert 8–9 iterations, 10.4 ms | 4–5 iterations, 0.9–1.3 ms |
| 223 | every Type I and II form fails | all converge in 4–12 iterations (I-lp 6.7 ms, I-bp 18.5 ms, III-diff 22.9 ms) |
| 255 | Type I fails | converges (I-lp 17.8 ms) |
| 319 | Type I fails | converges (5–9 iterations, ripple 3e-12; I-lp 17.9 ms, I-bp 37.5 ms) |

The equilibrium seed alone gave I-lp 223 22.9 ms and I-lp 319 75.7 ms. The lazy scan (below) gave the rest. Designs that still fail (Type II–IV at 319 taps) take longer to fail: 207–989 ms against 210–767 ms, because their noise-level runs trigger the certified rescan.

- **Remaining limits:**
  - At 223 taps the Type III and IV Hilberts need tolerance 1e-5 (flatness about 2e-6).
  - From 255 taps, Type II–IV ripples of 6e-10 and below converge in the transformed problem, but the charged rounding `2*f/max_error` exceeds 1e-5 (1.2e-5 to 2.5e-5), so they are reported as `numerical_failure`. The Type III differentiator passes at 255 taps with 1e-5. These limits come from binary64, not from the start.
- **SciPy 1.18.1 at density 16:** it fails the 255-tap lowpass and 223-tap differentiator. At 223 and 224 taps it stops at max errors of 9.13e-9 and 8.74e-9, 4–5% above this library's minimax ripples of 8.78e-9 and 8.31e-9.

Exact zero targets in every form and equal constant targets in Type I are solved analytically. They return zero ripple/flatness and an empty extrema view because alternation is undefined for exact zero error. One-tap unequal constant bands use the active pair of weighted scalar constraints. The equilibrium reference is followed by a check for distinct target levels. It falls back to uniform active-grid indices for bands with fewer than two grid points or more bands than reference points; only that fallback needs the distinct-target-level check. The earlier multiband reseeds and peak-pair dropping never fired after the equilibrium seed (0 hits over the test suite and a 49-design sweep), so they were removed.

The general selector maximizes the weakest selected peak over ordered alternating subsequences. It sorts magnitudes in preallocated sample scratch, finds the largest feasible threshold, and reconstructs d+2 alternating peaks including a global strongest peak. The common d+2 and d+3 candidate cases retain the short path. The full searched maximum remains in the convergence check. For B constant-target bands, capacities are `max(d+2,max(d-1,0)+2B)` for I, `max(d+2,d+2B)` for II/IV, and `d+1+2B` for III. Differentiators use the conservative grid bound `M+d+2`; this consumes more memory than ordinary filters. Setup checks every count and byte bound.

With L=d+1, workspace length is `2*L+max(L,K,N-L)`: a cosine table, sample/selection scratch, and coefficients. Once converged, the scratch/coefficient region holds finished taps. Trusted coefficients map to taps through cosine/sine product identities; otherwise direct inverse-DFT recovery uses the full response phase and explicit DC/Nyquist terms. `Interpolator::taps(output,symmetry,N)` exposes that allocation-free direct recovery when the fitted degree matches N. Output symmetry is written exactly; Type III's centre is exactly zero. Recovery failure is independent of fit success, and read-only tap calls only copy the finished buffer.

The embedded sweep covers 54 designs across highpass, bandpass, bandstop, five-band signed-gain designs, even-length filters, Hilbert transformers, and relative differentiators, including the minimum lengths. Nine complete SciPy 1.18.1 coefficient fixtures are embedded. An external check at density 512 found maximum coefficient difference about 1.97e-5 of ripple over 44 reference designs; 44 small cases were checked against an independent linear-programming minimax solve on 8,193 points per band, with maximum ripple difference about 2.23e-6. Those LP comparisons are discrete. SciPy rejected the three-tap bandpass and bandstop cases; both passed the LP comparison. Exhaustive subset enumeration checks the peak selector, and analytic cosine-basis tests check direct recovery for every form through degree 127. No Python dependency is required by normal tests.

Local Clang -O2 measurements (five-trial medians, 100 restarts per trial, density 16, construction excluded):

| Design | run() (ms) | Owned buffer payload |
| --- | --- | --- |
| Type II, 32-tap lowpass | 0.150 | 9,640 bytes |
| Type III, 31-tap Hilbert | 0.120 | 8,928 bytes |
| Type IV, 32-tap full-band relative differentiator | 0.187 | 18,520 bytes |

The engine object is 240 bytes on this x64 build; caller output and allocator overhead are extra. Repeated local MSVC before/after comparisons found broadly comparable lowpass design times, with identical iteration counts and reported ripple. Timing variation exceeded some of the small differences, so no general speedup is claimed for generalization.

## Checks

```powershell
.\tools\verify-sanitizers.ps1
```

Runs all five phases through MSVC/Clang debug and release, MSVC ASan, Clang ASan + UBSan, and debug-CRT tracking. Phase 1 uses a half-angle Taylor reference and analytic cosine landmarks. Phase 2 checks known polynomials through degree 127, an independent pivoted Chebyshev linear-system solve, actual grid integration, near-node queries, signed ripple, common gain/weight scales from 1e-200 to 1e200, and invalid/ill-conditioned cases. Range regressions cover fractional-subnormal amplitudes, their weighted errors, overflow rejection, and finite weighted errors despite overflowing unweighted amplitudes.

CRT checks verify first-use/repeated allocation-free sizing, filling, fitting, evaluation, exchange and tap recovery; one allocation per owner, including the 5,160-byte interpolator payload for 129 nodes (5,207 bytes requested on this MSVC debug CRT); and cleanup on success, numerical failure, and forced allocation failure. Deliberate faults validate the detectors. `MDd` with `_ITERATOR_DEBUG_LEVEL=0` excludes STL iterator proxies while retaining debug-heap tracking. Each of the six exchange setup allocations is forced to fail in turn, checking cleanup of earlier buffers. Coverage is the monitored CRT heap. Runtime paths are scoped to the script.

Quick normal build:

```powershell
xmake f -c -y -m debug --toolchain=msvc --builddir=build/msvc-debug
xmake -y
xmake run tests
```
## Local performance

MSVC Release, five-trial medians, same known quadratic data and query sequence before/after the Phase 2 precision repair. These historical timings predate the Phase 3 omitted-node selection change. These numbers describe the scalar benchmark, not every filter specification.

| Nodes | Fit before/after (us) | Evaluation before/after (ns) |
| --- | --- | --- |
| 17 | 2.01 / 3.78 | 69.9 / 45.5 |
| 65 | 23.6 / 31.5 | 260 / 158 |
| 129 | 91.5 / 108 | 523 / 308 |

Ordinary evaluation improves about 1.5-1.7x. Fitting is slower; difficult evaluations use the more expensive compensated path. The precision repair originally grew node storage from 24 to 48 bytes. Reusing the not-yet-computed value storage for exact integer exponents reduces it to 40 bytes per node: 5,160 instead of 6,192 bytes for 129 nodes (16.7% less workspace), still in one preallocated vector. The fit bounds node count so scratch exponents round-trip exactly through binary64; both value components are overwritten before evaluation. These figures cover the vector payload, excluding allocator overhead and the interpolator object. No dependency is added to production code.

Exchange, `--benchmark` five-trial medians (ms), Release, weights 1 and 10, density 16; ripples agree to about 12 significant digits. "Original" is the 32-step golden-section refinement with an always two-component scan.

| Taps | Edges (x pi) | MSVC original / current | Clang original / current |
| --- | --- | --- | --- |
| 15 | 0.2, 0.3 | 0.62 / 0.042 | 0.50 / 0.037 |
| 63 | 0.2, 0.3 | 7.94 / 0.26 | 7.91 / 0.20 |
| 127 | 0.2, 0.3 | 35.3 / 1.14 | 35.6 / 0.85 |
| 255 | 0.24, 0.26 | 120.9 / 4.5 | 122.9 / 3.2 |

The equilibrium seed and lazy scan were measured against the equispaced-start code in the same session (MSVC / Clang ms):

| Taps | Equispaced start | Current | Change |
| --- | --- | --- | --- |
| 15 | 0.041 / 0.032 | 0.042 / 0.037 | **3% / 15% slower** |
| 63 | 0.54 / 0.41 | 0.26 / 0.20 | 2.1× / 2.0× faster |
| 127 | 4.42 / 3.60 | 1.14 / 0.85 | 3.9× / 4.2× faster |
| 255 | 8.25 / 5.87 | 4.50 / 3.23 | 1.8× / 1.8× faster |

- **Iterations:** fall from 5/8/10/9 to 4/4/5/6.
- **127 taps:** this design never qualified for the certified Clenshaw scan, so the lazy scan adds about 30% on top of the seed.
- **15 taps:** slower because the seed's fixed quadrature cost exceeds the iteration it saves.
- **Lazy scan not used everywhere:** when Clenshaw certifies a whole iteration, the plain certified scan runs instead. Running the lazy one there cost 7–24% at 15/63 taps.
- **Tried and dropped, two-stage Brent:** locating peaks in binary64 and polishing in double-double was slower in every case (for example I-lp 223 at 6.4 to 6.9 ms), because the polish needs about as many precise evaluations. Brent refinement remains the largest cost for long designs: about 4.9 of 7.8 ms at 223 taps.

The other gains come from Brent refinement, the certified binary64 path, certifying to `tolerance/256` of the ripple rather than 2^-32 of each error, the branch-free two-sum, and the Chebyshev scan (which the 127-tap design reaches only through the lazy scan). `std::reduce`/`transform_reduce` cannot help: they may reorder floating-point sums. SciPy 1.18.1 `remez` (C core, uncompensated, no refinement) remains faster: about 2.1 ms at 255 taps, against 3.2 ms (Clang) here.

Targeting AVX2/FMA (x86-64-v3) was measured and not adopted: ripples stayed bit-identical, but exchange time improved only about 10% for MSVC and not at all for Clang, because the remaining cost is scalar division and double-double arithmetic.

After a Release build, use `xmake run tests --benchmark` for timings or `--probe` to export exact inputs/results of the 255-tap case. Use `--exchange-probe` to export the Phase 3 specification sweep and response samples. The optional probes support independent reference checks; they do not alter the algorithm.

