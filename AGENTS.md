# C++20 project guidelines

These rules apply to all work in this directory. Keep the implementation compact, translate the mathematics directly, and support correctness and performance claims with measurements.

## Scope and build

- Support Types I-IV with ordered disjoint bands, constant targets/weights, and relative-error differentiators. Preserve the Type-I lowpass convenience API and direct recovery fallback.
- Use standard C++20 and xmake. MSVC is the primary compiler; Clang is the independent compiler and UBSan build. No compiler extensions in the numerical core.
- Defer arbitrary callbacks, FFT dependencies, SIMD, parallel execution, and high-order guarantees. New numerical paths require independent checks before performance claims.
- Keep warnings enabled and treat them as errors. Fix warnings; do not broadly suppress them. Keep debug, sanitizer, and optimized builds separate.

## Code and standard library

- Use small functions, value semantics, aggregates, and RAII. Prefer concrete double arithmetic initially; introduce templates only for an actual second use.
- Use const for values that do not change. Use [[nodiscard]] for results or status values callers must inspect. Mark functions noexcept only when their complete call path cannot throw.
- Use std::array for small fixed storage, pre-sized std::vector for owned runtime storage, and std::span for borrowed contiguous buffers. Borrowed views must not outlive owners or survive owner reallocation.
- Use ranges algorithms and projections for recognizable operations such as validation, transforms, and extrema selection. Use direct loops for coupled sums and numerical recurrences when they express the equations more clearly.
- Use views only when lazy evaluation is useful. Cache repeated expensive calculations such as grid cosines. Do not build long pipelines that obscure arithmetic order or recompute values.
- Use std::numbers::pi_v<double>, std::lerp for band sampling, and std::midpoint for bracket bisection where appropriate. Verify exact band endpoints.
- Concepts constrain genuinely generic interfaces; constexpr applies to genuinely constant work. Do not assume std::cos is constexpr in C++20.
- Do not accidentally use C++23 facilities: std::expected, std::mdspan, std::ranges::to, or std::views::zip.
- Avoid owning raw pointers, manual new/delete, shared ownership without need, std::function in numerical kernels, speculative class hierarchies, and custom arenas without evidence.
- Keep names tied to equations. Comments explain numerical reasoning or invariants; do not narrate obvious code.

## Memory contract

- Make all setup allocations explicit and bounded by input sizes. Check size arithmetic before allocating.
- Allocate workspace and output before computation. Remez iterations and tap recovery must perform zero dynamic allocations after successful setup.
- Size buffers before writing. reserve changes capacity, not size. Never grow containers, build strings, or log inside the measured numerical region.
- Reuse scratch storage with clear lifetimes. Avoid accidental container copies and temporary materializations. std::span indexing is not guaranteed to be bounds-checked in C++20.
- Setup may report allocation failure. Runtime numerical failures use a small status enum and scalar diagnostics; never silently return invalid taps or claim convergence.
- RAII must release resources on successful, failed, and exceptional exits. Keep allocations in the test harness outside the measured region.

## Numerical contract

- Frequencies are radians/sample. Validate type-dependent tap count, ordered disjoint bands, positive widths/weights, finite targets, forced-endpoint compatibility, and checked storage arithmetic before computation.
- Keep desired response, weight, signed error, ripple, and extrema definitions consistent with the derivation. Preserve band boundaries during peak searches; never compare across the ignored transition as if it were an active interval.
- Barycentric weights must match the evaluation node set. Weights computed for r+2 nodes cannot be reused unchanged after dropping a node.
- Preserve deliberate summation order. Use compensated accumulation when justified by accuracy checks. reduce and transform_reduce may reorder floating-point sums; do not substitute them blindly.
- Handle interpolation nodes explicitly. Tolerances must reflect scale and rounding, not an arbitrary universal epsilon. Scaling weights does not cure ill-conditioning.
- Keep fast-math disabled. Use FMA deliberately, with accuracy checks. Check non-finite intermediates and provide an iteration limit and explicit numerical-failure/nonconvergence status.
- Equal ripple at selected extrema is not sufficient: independently search the active bands for larger missed peaks. State whether verification is sampled, refined, or rigorously bounded.

## Required verification

- Compare reconstructed tap responses against the internal amplitude and trustworthy reference designs. Independently check symmetry, finite taps, signed alternation, maximum weighted error, and ripple flatness.
- Cover ordinary cases, small tap counts, narrow transitions, extreme weight ratios, edge-adjacent bands, invalid inputs, and nonconvergence. Set tolerances explicitly; do not require bitwise agreement between compilers.
- Run MSVC ASan and Clang ASan + UBSan. tools/verify-sanitizers.ps1 verifies the detector setup using clean and deliberately faulty cases. Fault probes are non-default targets, never production code.
- Add a separate MSVC debug-CRT target using MDd and _DEBUG for leak snapshots and allocation hooks. Verify no additional outstanding allocations after workspace/output destruction and no allocation/reallocation events during computation. Keep hooks allocation-free and validate detectors with deliberate leaks and temporary allocations.
- CRT tracking covers allocations routed through the monitored CRT, not arbitrary Windows heaps or foreign allocators. Warm up unrelated runtime facilities outside measurement; do not hide first-use allocations in the numerical API.
- Leak/allocation instrumentation is required before claiming those contracts; it is not implemented merely because it is described here.
- Benchmark optimized unsanitized builds with equivalent specifications and verified accuracy. Report failures and slower cases as well as wins. Add optimization only after profiling identifies a bottleneck.
