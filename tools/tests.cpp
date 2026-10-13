#include "../remez.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <source_location>
#include <stdexcept>
#include <string_view>
#include <utility>

#ifdef REMEZ_CRT_CHECKS
#include <crtdbg.h>
#endif

namespace {
using remez::Lowpass;
using remez::Status;
constexpr double pi = std::numbers::pi_v<double>;

void check(bool condition, const char* message,
           std::source_location at = std::source_location::current()) {
    if (!condition) {
        std::fprintf(stderr, "%s:%u: %s\n", at.file_name(), static_cast<unsigned>(at.line()), message);
        std::exit(1);
    }
}

bool same_point(const remez::GridPoint& a, const remez::GridPoint& b) {
    return a.omega == b.omega && a.x == b.x && a.desired == b.desired && a.weight == b.weight;
}

constexpr remez::GridPoint sentinel{-7, -7, -7, -7};
constexpr double tolerance = 16 * std::numeric_limits<double>::epsilon();

// cos(w) = 2*cos(w/2)^2 - 1. Taylor remainder on [0, pi/2] is below 4e-24.
double cosine_reference(double omega) {
    const double x = omega / 2;
    double term = 1.0, sum = 1.0;
    for (int k = 1; k <= 13; ++k) {
        term *= -x * x / ((2 * k - 1) * (2 * k));
        sum += term;
    }
    return 2 * sum * sum - 1;
}

void verify(const Lowpass& spec) {
    const auto size = remez::grid_size(spec);
    check(size.status == Status::ok, "valid sizing");
    const auto grid = remez::make_grid(spec);
    check(grid.status == Status::ok, "valid generation");
    check(grid.points.size() == size.count(), "sizing matches ownership");
    check(grid.points.size() >= spec.taps / 2 + 2, "enough exchange nodes");
    const std::array views{grid.passband(), grid.stopband()};
    const std::array edges{std::pair{0.0, spec.pass_edge}, std::pair{spec.stop_edge, pi}};
    const double step = pi / (static_cast<double>(spec.density) * static_cast<double>(spec.taps / 2 + 1));
    for (std::size_t b = 0; b < 2; ++b) {
        const auto points = views[b];
        check(points.size() >= 2, "both band endpoints exist");
        check(points.front().omega == edges[b].first, "exact lower endpoint");
        check(points.back().omega == edges[b].second, "exact upper endpoint");
        for (std::size_t j = 0; j < points.size(); ++j) {
            const auto& point = points[j];
            check(point.omega >= edges[b].first && point.omega <= edges[b].second, "active band membership");
            check(point.desired == (b == 0 ? 1.0 : 0.0), "desired response");
            check(point.weight == spec.weights[b], "weight preservation");
            const double expected_omega = edges[b].first +
                (edges[b].second - edges[b].first) * static_cast<double>(j) / static_cast<double>(points.size() - 1);
            check(std::abs(point.omega - expected_omega) <= pi * tolerance, "uniform band sampling");
            check(std::abs(point.x - cosine_reference(point.omega)) <= tolerance, "independent cosine reference");
            if (j != 0) {
                const auto gap = point.omega - points[j - 1].omega;
                check(gap > 0 && gap <= step + 8 * pi * std::numeric_limits<double>::epsilon(), "positive bounded spacing");
                check(point.x < points[j - 1].x, "distinct decreasing interpolation coordinates");
            }
        }
    }
    check(grid.passband().back().omega < grid.stopband().front().omega, "transition excluded");
    auto storage = std::vector<remez::GridPoint>(size.count() + 1);
    storage.back() = sentinel;
    check(remez::fill_grid(spec, storage) == Status::ok, "caller storage fill");
    check(std::ranges::equal(grid.points, std::span{storage}.first(size.count()), same_point), "deterministic fill");
    check(same_point(storage.back(), sentinel), "extra storage untouched");
    check(remez::fill_grid(spec, std::span{storage}.first(size.count() - 1)) == Status::insufficient_storage, "short buffer rejected");
    check(remez::make_grid(spec, size.count() - 1).status == Status::grid_limit, "resource limit enforced");
    check(remez::make_grid(spec, size.count()).status == Status::ok, "exact resource limit accepted");
}

void invalid_inputs() {
    const Lowpass base{31, 0.2 * pi, 0.3 * pi};
    const auto expect = [](Lowpass spec, Status status) {
        check(remez::grid_size(spec).status == status, "invalid sizing status");
        const auto result = remez::make_grid(spec);
        check(result.status == status && result.points.empty(), "invalid result has no grid");
        std::array storage{sentinel, sentinel};
        check(remez::fill_grid(spec, storage) == status, "invalid fill status");
        check(std::ranges::all_of(storage, [](const auto& p) { return same_point(p, sentinel); }), "invalid input leaves storage untouched");
    };
    for (auto taps : {std::size_t{0}, std::size_t{2}, std::size_t{32}}) {
        auto spec = base; spec.taps = taps; expect(spec, Status::invalid_taps);
    }
    for (double edge : {0.0, -1.0, pi, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        auto spec = base; spec.pass_edge = edge; expect(spec, Status::invalid_edges);
        spec = base; spec.stop_edge = edge; expect(spec, Status::invalid_edges);
    }
    auto spec = base;
    spec.stop_edge = spec.pass_edge; expect(spec, Status::invalid_edges);
    spec = base; std::swap(spec.pass_edge, spec.stop_edge); expect(spec, Status::invalid_edges);
    for (double weight : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        for (std::size_t b = 0; b < 2; ++b) {
            spec = base; spec.weights[b] = weight; expect(spec, Status::invalid_weights);
        }
    }
    spec = base; spec.density = 0; expect(spec, Status::invalid_density);
    spec = base; spec.density = std::numeric_limits<std::size_t>::max(); expect(spec, Status::grid_limit);
    spec = base; spec.taps = std::numeric_limits<std::size_t>::max(); expect(spec, Status::grid_limit);
    check(remez::grid_size(base, 0).status == Status::grid_limit, "zero budget");
    spec = {255, 1e-4, pi - 1e-4}; expect(spec, Status::insufficient_points);
    spec = {1, std::numeric_limits<double>::denorm_min(), pi / 2};
    const auto unresolved = remez::make_grid(spec);
    check(unresolved.status == Status::numerical_resolution && unresolved.points.empty(), "coincident cosines rejected");
    spec = {1, pi / 2, std::nextafter(pi, 0.0)};
    check(remez::make_grid(spec).status == Status::numerical_resolution, "Nyquist resolution failure");
    std::array untouched{sentinel, sentinel};
    check(remez::fill_grid(base, untouched) == Status::insufficient_storage, "short storage status");
    check(std::ranges::all_of(untouched, [](const auto& p) { return same_point(p, sentinel); }), "short storage untouched");
}

double chebyshev(std::size_t degree, double x) {
    double previous = 1.0, current = x;
    if (degree == 0) return previous;
    for (std::size_t k = 1; k < degree; ++k) {
        const double next = 2 * x * current - previous;
        previous = current;
        current = next;
    }
    return current;
}

double polynomial(std::size_t degree, double x) {
    return 0.125 + 0.5 * chebyshev(degree, x) + (degree == 0 ? 0.0 : 0.25 * x);
}

std::array<double, 12> dense_reference(std::span<const remez::GridPoint> points) {
    const auto n = points.size();
    check(n >= 2 && n <= 12, "reference matrix capacity");
    std::array<std::array<double, 13>, 12> matrix{};
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j + 1 < n; ++j)
            matrix[k][j] = chebyshev(j, points[k].x);
        matrix[k][n - 1] = (k % 2 == 0 ? 1.0 : -1.0) / points[k].weight;
        matrix[k][n] = points[k].desired;
    }
    for (std::size_t k = 0; k < n; ++k) {
        auto pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(matrix[i][k]) > std::abs(matrix[pivot][k])) pivot = i;
        check(matrix[pivot][k] != 0.0, "reference system nonsingular");
        std::swap(matrix[pivot], matrix[k]);
        const double divisor = matrix[k][k];
        for (std::size_t j = k; j <= n; ++j) matrix[k][j] /= divisor;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double multiplier = matrix[i][k];
            for (std::size_t j = k; j <= n; ++j) matrix[i][j] -= multiplier * matrix[k][j];
        }
    }
    std::array<double, 12> solution{};
    for (std::size_t k = 0; k < n; ++k) solution[k] = matrix[k][n];
    return solution;
}

void interpolation() {
    for (auto n : {std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{6}, std::size_t{10},
                   std::size_t{17}, std::size_t{33}, std::size_t{65}, std::size_t{129}}) {
        std::vector<remez::GridPoint> points(n);
        remez::Interpolator interpolator(n);
        const double bound = 128 * static_cast<double>(n) * std::numeric_limits<double>::epsilon();
        for (double delta : {0.0, 0.03125, -0.03125}) {
            for (std::size_t k = 0; k < n; ++k) {
                const double x = cosine_reference(pi * static_cast<double>(k) / static_cast<double>(n - 1));
                const double weight = std::ldexp(1.0, static_cast<int>(k % 5) - 2);
                points[k] = {0.0, x, polynomial(n - 2, x) + (k % 2 == 0 ? delta : -delta) / weight, weight};
            }
            const auto fit = interpolator.fit(points);
            check(fit.status == Status::ok, "known-polynomial fit");
            check(std::abs(fit.value - delta) <= bound, "known signed ripple");
            for (int j = 0; j <= 400; ++j) {
                const double x = -1.0 + static_cast<double>(j) / 200;
                const auto value = interpolator.evaluate(x);
                check(value.status == Status::ok && std::abs(value.value - polynomial(n - 2, x)) <= bound, "polynomial recovery");
                // polynomial() itself rounds: the three-term recurrence errs by at most (n+2) eps.
                const auto bounded = interpolator.estimate(x);
                const double reference_error = static_cast<double>(n + 2) * std::numeric_limits<double>::epsilon();
                check(bounded.status == Status::ok && std::isfinite(bounded.bound) &&
                      std::abs(bounded.value - polynomial(n - 2, x)) <= bounded.bound + reference_error,
                      "estimate bound covers known polynomial");
            }
            for (std::size_t k = 0; k < n; ++k) {
                const auto value = interpolator.evaluate(points[k].x);
                check(value.status == Status::ok, "extremum evaluation");
                const double error = points[k].weight * (points[k].desired - value.value);
                check(std::abs(error - (k % 2 == 0 ? delta : -delta)) <= 4 * bound, "all extremal equations including omitted node");
                const double nearby = std::nextafter(points[k].x, 0.0);
                const auto neighbor = interpolator.evaluate(nearby);
                check(neighbor.status == Status::ok && std::abs(neighbor.value - polynomial(n - 2, nearby)) <= bound,
                      "near-node finite evaluation");
            }
        }
    }
    for (auto n : {std::size_t{2}, std::size_t{3}, std::size_t{6}, std::size_t{10}}) {
        std::vector<remez::GridPoint> points(n);
        for (std::size_t k = 0; k < n; ++k) {
            const double x = cosine_reference(pi * static_cast<double>(k) / static_cast<double>(n - 1));
            points[k] = {0.0, x, k < n / 2 ? 1.0 : 0.0, k < n / 2 ? 1.0 : 10.0};
        }
        const auto reference = dense_reference(points);
        remez::Interpolator interpolator(n);
        const auto fit = interpolator.fit(points);
        check(fit.status == Status::ok && std::abs(fit.value - reference[n - 1]) <= 2e-12, "independent linear-system ripple");
        for (int j = 0; j <= 400; ++j) {
            const double x = -1.0 + static_cast<double>(j) / 200;
            double expected = 0.0;
            for (std::size_t k = 0; k + 1 < n; ++k) expected += reference[k] * chebyshev(k, x);
            const auto value = interpolator.evaluate(x);
            check(value.status == Status::ok && std::abs(value.value - expected) <= 2e-12, "independent linear-system amplitude");
        }
        for (double factor : {1e-200, 1e200}) {
            auto scaled = points;
            for (auto& point : scaled) point.weight *= factor;
            const auto fit_scaled = interpolator.fit(scaled);
            check(fit_scaled.status == Status::ok && std::abs(fit_scaled.value / factor - fit.value) <= 2e-12, "weight scaling ripple");
            const auto value = interpolator.evaluate(0.125);
            double expected = 0.0;
            for (std::size_t k = 0; k + 1 < n; ++k) expected += reference[k] * chebyshev(k, 0.125);
            check(value.status == Status::ok && std::abs(value.value - expected) <= 2e-12, "weight scaling preserves amplitude");
            scaled = points;
            for (auto& point : scaled) point.desired *= factor;
            const auto fit_gain = interpolator.fit(scaled);
            const auto amplitude = interpolator.evaluate(0.125);
            check(fit_gain.status == Status::ok && std::abs(fit_gain.value / factor - fit.value) <= 2e-12, "gain scaling ripple");
            check(amplitude.status == Status::ok && std::abs(amplitude.value / factor - expected) <= 2e-12, "gain scaling amplitude");
        }
    }
}

void grid_interpolation() {
    for (auto taps : {std::size_t{1}, std::size_t{3}, std::size_t{15}, std::size_t{31}, std::size_t{63}, std::size_t{127}, std::size_t{255}}) {
        const auto grid = remez::make_grid({taps, 0.2 * pi, 0.3 * pi, {1.0, 10.0}});
        check(grid.status == Status::ok, "interpolation integration grid");
        const auto n = taps / 2 + 2;
        std::vector<remez::GridPoint> extrema(n);
        for (std::size_t k = 0; k < n; ++k) extrema[k] = grid.points[k * (grid.points.size() - 1) / (n - 1)];
        remez::Interpolator interpolator(n);
        const auto fit = interpolator.fit(extrema);
        check(fit.status == Status::ok, "Phase 1 grid to Phase 2 fit");
        const auto reference = n <= 12 ? dense_reference(extrema) : std::array<double, 12>{};
        if (n <= 12) check(std::abs(fit.value - reference[n - 1]) <= 2e-12, "integrated independent ripple");
        for (const auto& point : grid.points) {
            const auto value = interpolator.evaluate(point.x);
            check(value.status == Status::ok, "full active-grid evaluation");
            if (n <= 12) {
                double expected = 0.0;
                for (std::size_t k = 0; k + 1 < n; ++k) expected += reference[k] * chebyshev(k, point.x);
                check(std::abs(value.value - expected) <= 2e-12, "integrated independent amplitude");
            }
        }
        for (std::size_t k = 0; k < n; ++k) {
            const auto value = interpolator.evaluate(extrema[k].x);
            const double error = extrema[k].weight * (extrema[k].desired - value.value);
            check(value.status == Status::ok && std::abs(error - (k % 2 == 0 ? fit.value : -fit.value)) <= 2e-11,
                  "integrated extremal equations");
        }
        // A binary64 result accepted at a looser accuracy must honour that accuracy against the
        // default 2^-32 certification, including its sign.
        for (const auto& point : grid.points) {
            const auto exact = interpolator.weighted_error(point);
            for (const double accuracy : {0x1p-10, 0x1p-20, 0x1p-28}) {
                const auto loose = interpolator.weighted_error(point, accuracy);
                check(exact.status == Status::ok && loose.status == Status::ok &&
                      std::abs(loose.value - exact.value) <= (accuracy + 0x1p-31) * std::abs(exact.value),
                      "certified weighted-error accuracy");
                // An absolute floor relaxes accuracy near zero crossings but never the sign.
                const double floor = accuracy * std::abs(fit.value);
                const auto floored = interpolator.weighted_error(point, accuracy, floor);
                check(floored.status == Status::ok && std::signbit(floored.value) == std::signbit(exact.value) &&
                      std::abs(floored.value - exact.value) <=
                          std::max(accuracy * std::abs(exact.value), floor) + 0x1p-31 * std::abs(exact.value),
                      "certified weighted-error absolute floor");
            }
        }
    }
}
void range_regression() {
    const double u = std::numeric_limits<double>::denorm_min();
    const std::array<remez::GridPoint, 3> tiny{{{0,1,2*u,1}, {0,0,u,1}, {0,-1,0,1}}};
    remez::Interpolator polynomial(3);
    check(polynomial.fit(tiny).status == Status::ok, "subnormal linear polynomial fit");
    for (const double x : {-0.75, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto estimate = polynomial.estimate(x);
        const auto value = polynomial.evaluate(x);
        check(estimate.status == Status::ok && value.status == Status::ok && std::isfinite(estimate.bound) &&
              std::abs(estimate.value/u-(1+x)) <= estimate.bound/u, "subnormal estimate includes absolute rounding");
        for (const double desired : {0.0, u, 2*u}) {
            const auto error = polynomial.weighted_error({0,x,desired,1e308});
            const double expected = (desired/u-(1+x))*(u*1e308);
            check(error.status == Status::ok && std::abs(error.value-expected) <= 16*std::numeric_limits<double>::epsilon()*std::abs(expected),
                  "weighting preserves fractional subnormal amplitudes");
        }
    }
    check(polynomial.weighted_error({0,0.25,0,1}).status == Status::numerical_failure,
          "uncertifiable relative subnormal error is rejected");
    check(polynomial.weighted_error({0,0.0,0,1}).status == Status::ok &&
          polynomial.weighted_error({0,0.0,0,1}).value == -u, "exact subnormal weighted error remains usable");
    auto zero = tiny;
    for (auto& point : zero) point.desired = 0;
    check(polynomial.fit(zero).status == Status::ok, "zero polynomial fit");
    check(polynomial.weighted_error({0,0.25,u,0.25}).status == Status::numerical_failure,
          "weighted multiplication underflow cannot certify zero");
    const double large = 0.95*std::numeric_limits<double>::max();
    const std::array<remez::GridPoint, 3> huge{{{0,0.9,large,1}, {0,0,0,1}, {0,-0.9,-large,1}}};
    check(polynomial.fit(huge).status == Status::ok, "finite large linear polynomial fit");
    for (const double x : {-1.0, 1.0}) {
        check(polynomial.evaluate(x).status == Status::numerical_failure && polynomial.estimate(x).status == Status::numerical_failure,
              "overflowed amplitudes never return success");
        const auto error = polynomial.weighted_error({0,x,0,1e-308});
        const double expected = -(large*1e-308)*(x/0.9);
        check(error.status == Status::ok && std::abs(error.value-expected) <= 16*std::numeric_limits<double>::epsilon()*std::abs(expected),
              "finite weighted error survives overflowing amplitude");
    }
    check(polynomial.evaluate(0.9).status == Status::ok && polynomial.estimate(0.9).status == Status::ok,
          "representable large node amplitude succeeds");
}

void invalid_interpolation() {
    std::array<remez::GridPoint, 3> points{{{0, 1, 1, 1}, {0, 0, 1, 1}, {0, -1, 0, 1}}};
    remez::Interpolator interpolator(3);
    check(interpolator.evaluate(0).status == Status::invalid_nodes && interpolator.estimate(0).status == Status::invalid_nodes,
          "unfitted evaluator");
    check(interpolator.fit(std::span{points}.first(2)).status == Status::invalid_nodes, "wrong node count");
    for (const double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), 2.0, -2.0}) {
        auto invalid = points; invalid[1].x = bad;
        check(interpolator.fit(invalid).status == Status::invalid_nodes, "invalid coordinate");
    }
    for (const double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        auto invalid = points; invalid[1].desired = bad;
        check(interpolator.fit(invalid).status == Status::invalid_nodes, "nonfinite desired value");
    }
    for (const double bad : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        auto invalid = points; invalid[1].weight = bad;
        check(interpolator.fit(invalid).status == Status::invalid_weights, "invalid interpolation weight");
    }
    auto invalid = points; invalid[1].x = invalid[0].x;
    check(interpolator.fit(invalid).status == Status::numerical_resolution, "duplicate nodes");
    invalid = points; std::swap(invalid[0], invalid[1]);
    check(interpolator.fit(invalid).status == Status::invalid_nodes, "unordered nodes");
    check(interpolator.fit(points).status == Status::ok, "recovery after rejected fit");
    for (double bad : {-2.0, 2.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(interpolator.evaluate(bad).status == Status::invalid_nodes, "invalid evaluation coordinate");
    for (double bad : {0.0, -0.5, 1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(interpolator.weighted_error(points[1], bad).status == Status::invalid_tolerance, "invalid weighted-error accuracy");
    for (double bad : {-1e-300, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(interpolator.weighted_error(points[1], 0x1p-32, bad).status == Status::invalid_tolerance,
              "invalid weighted-error absolute accuracy");
    for (double bad : {-2.0, 2.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(interpolator.estimate(bad).status == Status::invalid_nodes, "invalid estimate coordinate");
    invalid = points; invalid[0].weight = std::numeric_limits<double>::denorm_min();
    invalid[1].weight = std::numeric_limits<double>::max();
    check(interpolator.fit(invalid).status == Status::numerical_failure, "unrepresentable weight ratio");
    check(interpolator.evaluate(0).status == Status::invalid_nodes, "failed fit invalidates previous interpolant");
    remez::Interpolator constant(2);
    const std::array<remez::GridPoint, 2> dynamic_range{{{0, 1, 1, 1}, {0, -1, 0, 1e200}}};
    const auto constant_fit = constant.fit(dynamic_range);
    check(constant_fit.status == Status::ok && constant_fit.value == 1.0 &&
          std::abs(constant.evaluate(0).value/1e-200-1) <= tolerance,
          "conditioned omission resolves extreme-weight constant exactly");
    bool oversized = false;
    try {
        remez::Interpolator impossible(std::numeric_limits<std::size_t>::max());
    } catch (const std::length_error&) {
        oversized = true;
    }
    check(oversized, "oversized workspace rejected before allocation");
    remez::Interpolator empty(0), singleton(1);
    check(empty.fit({}).status == Status::invalid_nodes && singleton.fit(std::span{points}.first(1)).status == Status::invalid_nodes, "undersized interpolators");
    check(interpolator.fit(points).status == Status::ok, "fit before move");
    auto moved = std::move(interpolator);
    check(moved.evaluate(0.25).status == Status::ok && interpolator.evaluate(0.25).status == Status::invalid_nodes, "move preserves ownership and invalidates source");
    const auto cached = moved.evaluate(0.25);
    points[0].desired = 100.0;
    check(moved.evaluate(0.25).value == cached.value, "fit owns cached values independently of input");
    std::array<remez::GridPoint, 3> zero{{{0, 1e-200, 0, 1}, {0, 0, 0, 1}, {0, -1e-200, 0, 1}}};
    check(moved.fit(zero).status == Status::ok && moved.evaluate(0.5).value == 0.0, "log scaling avoids raw product overflow");
    std::array<remez::GridPoint, 3> clustered{{{0, 1, 1, 1}, {0, std::nextafter(1.0, 0.0), 1, 1},
                                            {0, std::nextafter(std::nextafter(1.0, 0.0), 0.0), 1, 1}}};
    const auto clustered_fit = moved.fit(clustered);
    const auto constant_value = moved.evaluate(-1);
    check(clustered_fit.status == Status::ok && constant_value.status == Status::ok && constant_value.value == 1.0,
          "compensated evaluation resolves clustered constant polynomial");
    for (std::size_t k = 0; k < clustered.size(); ++k) clustered[k].x = (1.0 - static_cast<double>(k)) * 1e-100;
    const auto extreme_fit = moved.fit(clustered);
    check(extreme_fit.status == Status::numerical_failure ||
          (extreme_fit.status == Status::ok && moved.evaluate(-1).status == Status::numerical_failure), "unresolvable extrapolation rejected");
}
// Fixed binary64 inputs; Chebyshev-system references agree at 80 and 110 decimal digits.
void precision_regression() {
    constexpr std::array<double, 129> x{
        1, 0.99976985505525917, 0.99907952615442785, 0.99782959002511984,
        0.99618727421460784, 0.99391884098347871, 0.99132778145939771, 0.98828042378034853,
        0.98451064426965618, 0.98052279748316495, 0.97574932492388344, 0.97082886411043656,
        0.96546154058642264, 0.95921776695166927, 0.95293289065008324, 0.94571243525403725,
        0.93852241736050679, 0.93033953894239163, 0.9222585780890713, 0.9137531109367687,
        0.90417349017882842, 0.89480136248949882, 0.88430286752013565, 0.87408258709105247,
        0.8634599752845592, 0.85163769955901847, 0.84019744887149372, 0.82751102338556837,
        0.81527572568168771, 0.57906685948227155, 0.56143089760175446, 0.54353615422385171,
        0.5240854199785675, 0.50568099906014341, 0.48570355044505459, 0.46682584718472919,
        0.44773296892732278, 0.42704753562367997, 0.40753704124944984, 0.38642472697423014,
        0.36653519896947523, 0.34647672314222089, 0.3248085242469626, 0.3044295799858126,
        0.28243953261559879, 0.26178041978803335, 0.23951200980668191, 0.2186135823204661,
        0.19761438882649535, 0.17501443107432912, 0.15383734632652477, 0.1310691844520542,
        0.10975608666049512, 0.088392398770359995, 0.065457781612190499, 0.044020646366324731,
        0.021030032178019058, -0.00043815795719221834, -0.023439304666128584, -0.044896095812486464,
        -0.066332192909139151, -0.089265250525316878, -0.11062706613103906, -0.13345788099230665,
        -0.15470317158716923, -0.17587715453123318, -0.1984733476863228, -0.21946861745347712,
        -0.24185103782018133, -0.26262607574397867, -0.2832800608909804, -0.30526418445907533,
        -0.32563720133913104, -0.34729862568176489, -0.36735038614893739, -0.38864627526247597,
        -0.40833712646334952, -0.42783976199337964, -0.4485163699873152, -0.46760063720336398,
        -0.48780866127070355, -0.50643682056124428, -0.52483154688734379, -0.54427151255166306,
        -0.56215585408312574, -0.58102989413326989, -0.59836848414135657, -0.61663927469522961,
        -0.633397825688368, -0.64986442332006411, -0.66717460733285738, -0.68301245669818789,
        -0.69963203465178048, -0.71480981601099891, -0.7296581186127713, -0.74519362119218346,
        -0.75933783647917985, -0.7741037875137684, -0.78751594494961641, -0.80056511104259487,
        -0.81413672111603308, -0.8264150361149758, -0.83914752173361906, -0.85063070496623072,
        -0.86249888732412394, -0.87316423051088832, -0.88342710431450988, -0.8939710521435682,
        -0.90338540438973847, -0.91301002385868069, -0.92155723743244655, -0.9296796757349407,
        -0.93790666819230473, -0.94513752951894425, -0.95240136526467478, -0.95872635040793031,
        -0.9650126650267663, -0.97041926617171759, -0.97537856997007744, -0.98019315377951066,
        -0.98421887679366926, -0.98802868512626563, -0.99111286638916318, -0.99374021196179518,
        -0.99604690572122889, -0.99772433624080514, -0.99901123760477006, -0.99973543591085401,
        -1
    };
    constexpr std::array samples{
        std::pair<double, double>{1, 1},
        std::pair<double, double>{0.95339639205493054, 1.0000000000000808},
        std::pair<double, double>{0.94370242199830279, 1.0000000000005935},
        std::pair<double, double>{0.8595696069872012, 0.99999999800574513},
        std::pair<double, double>{0.83517578283731808, 0.99999993204634963},
        std::pair<double, double>{0.82144063849595106, 1.0000007841915672},
        std::pair<double, double>{0.82056570766019543, 1.000000919272841},
        std::pair<double, double>{0.80901699437494745, 0.99997924737101695},
        std::pair<double, double>{0.58778525229247314, 2.5704865596226826e-06},
        std::pair<double, double>{0.031761110554615793, -5.6422308166217013e-19},
        std::pair<double, double>{-0.77604180145198209, -2.3184696334921911e-24},
        std::pair<double, double>{-1, -3.3662658688937235e-24}
    };
    std::array<remez::GridPoint, 129> extrema{};
    for (std::size_t k = 0; k < x.size(); ++k)
        extrema[k] = {0.0, x[k], k < 29 ? 1.0 : 0.0, k < 29 ? 1.0 : 10.0};
    remez::Interpolator interpolator(extrema.size());
    const auto fit = interpolator.fit(extrema);
    constexpr double delta = 3.366265868893723709187224e-23;
    check(fit.status == Status::ok && std::abs(fit.value / delta - 1.0) <= 2e-13, "independent high-precision ripple");
    for (std::size_t k = 0; k < extrema.size(); ++k) {
        const auto error = interpolator.weighted_error(extrema[k]);
        check(error.status == Status::ok && std::abs(error.value/delta-(k%2 == 0 ? 1 : -1)) <= 2e-8,
              "tiny weighted error survives rounded unit amplitude");
    }
    for (const auto& [coordinate, expected] : samples) {
        const auto value = interpolator.evaluate(coordinate);
        check(value.status == Status::ok && std::abs(value.value - expected) <= 2e-12 * std::abs(expected) + 1e-34,
              "independent high-precision amplitude");
        // The references carry about 17 digits, so allow one rounding of `expected` beyond the bound.
        const auto bounded = interpolator.estimate(coordinate);
        check(bounded.status == Status::ok && std::isfinite(bounded.bound) &&
              std::abs(bounded.value - expected) <= bounded.bound + std::numeric_limits<double>::epsilon() * std::abs(expected),
              "estimate bound covers independent high-precision amplitude");
    }
}

void probe() {
    const auto grid = remez::make_grid({255, 0.2 * pi, 0.3 * pi, {1.0, 10.0}});
    std::array<remez::GridPoint, 129> extrema{};
    for (std::size_t k = 0; k < extrema.size(); ++k) {
        extrema[k] = grid.points[k * (grid.points.size() - 1) / (extrema.size() - 1)];
        std::printf("node %.17g %.17g %.17g\n", extrema[k].x, extrema[k].desired, extrema[k].weight);
    }
    remez::Interpolator interpolator(extrema.size());
    const auto fit = interpolator.fit(extrema);
    std::printf("ripple %d %.17g\n", static_cast<int>(fit.status), fit.value);
    for (std::size_t j = 0; j < grid.points.size(); ++j) {
        const auto value = interpolator.evaluate(grid.points[j].x);
        std::printf("query %zu %.17g %d %.17g\n", j, grid.points[j].x, static_cast<int>(value.status), value.value);
    }
}

void benchmark() {
    using Clock = std::chrono::steady_clock;
    double checksum = 0.0;
    for (auto n : {std::size_t{17}, std::size_t{65}, std::size_t{129}}) {
        std::vector<remez::GridPoint> extrema(n);
        for (std::size_t k = 0; k < n; ++k) {
            const double x = std::cos(pi * static_cast<double>(k) / static_cast<double>(n - 1));
            extrema[k] = {0, x, polynomial(2, x) + (k % 2 == 0 ? 0.03125 : -0.03125), 1};
        }
        remez::Interpolator interpolator(n);
        check(interpolator.fit(extrema).status == Status::ok, "benchmark fit");
        std::array<double, 5> fit_times{}, evaluation_times{};
        for (std::size_t trial = 0; trial < fit_times.size(); ++trial) {
            auto start = Clock::now();
            for (int repeat = 0; repeat < 100; ++repeat) checksum += interpolator.fit(extrema).value;
            fit_times[trial] = std::chrono::duration<double, std::nano>(Clock::now() - start).count() / 100;
            start = Clock::now();
            for (int repeat = 0; repeat < 20'000; ++repeat)
                checksum += interpolator.evaluate(-1.0 + static_cast<double>(repeat % 401) / 200).value;
            evaluation_times[trial] = std::chrono::duration<double, std::nano>(Clock::now() - start).count() / 20'000;
        }
        std::ranges::sort(fit_times);
        std::ranges::sort(evaluation_times);
        std::printf("nodes=%zu fit_us=%.3f evaluate_ns=%.3f\n", n, fit_times[2] / 1000, evaluation_times[2]);
    }
    for (const auto taps : {std::size_t{15}, std::size_t{63}, std::size_t{127}, std::size_t{255}}) {
        const double pass = taps == 255 ? 0.24 : 0.2, stop = taps == 255 ? 0.26 : 0.3;
        remez::Exchange engine({taps, pass*pi, stop*pi, {1,10}, 16});
        const auto warm = engine.run();
        check(warm.status == Status::ok, "exchange benchmark converges");
        std::array<double, 5> times{};
        for (auto& time : times) {
            const auto start = Clock::now();
            for (int repeat = 0; repeat < 10; ++repeat) {
                const auto result = engine.run();
                check(result.status == Status::ok, "exchange benchmark remains valid");
                checksum += result.ripple;
            }
            time = std::chrono::duration<double, std::milli>(Clock::now()-start).count()/10;
        }
        std::ranges::sort(times);
        std::printf("taps=%zu exchange_ms=%.6f iterations=%zu ripple=%.17g\n", taps, times[2], warm.iterations, warm.ripple);
        std::vector<double> h(taps);
        for (auto& time : times) {
            const auto start = Clock::now();
            for (int repeat = 0; repeat < 10'000; ++repeat) {
                check(engine.taps(h) == Status::ok, "tap benchmark recovery");
                checksum += h[taps/2];
            }
            time = std::chrono::duration<double, std::micro>(Clock::now()-start).count()/10'000;
        }
        std::ranges::sort(times);
        std::printf("taps=%zu recovery_us=%.6f\n", taps, times[2]);
    }
    std::printf("checksum=%.17g\n", checksum);
}
void exchange_probe() {
    for (const auto taps : {std::size_t{1}, std::size_t{3}, std::size_t{15}, std::size_t{31}, std::size_t{63}, std::size_t{127}, std::size_t{255}}) {
        for (const auto edges : {std::pair{0.2, 0.3}, std::pair{0.24, 0.26}, std::pair{0.01, 0.03}, std::pair{0.97, 0.99}}) {
            const Lowpass spec{taps, edges.first*pi, edges.second*pi, {1, 10}, 32};
            remez::Exchange engine(spec);
            const auto result = engine.run();
            std::printf("case %zu %.17g %.17g %d %zu %.17g %.17g %.17g\n", taps, edges.first, edges.second,
                static_cast<int>(result.status), result.iterations, result.ripple, result.max_error, result.flatness);
            if (result.status == Status::ok) {
                for (std::size_t j = 0; j <= 100; ++j) {
                    const double omega = pi * static_cast<double>(j)/100;
                    const auto amplitude = engine.evaluate(std::cos(omega));
                    std::printf("amplitude %zu %.17g %d %.17g\n", j, omega, static_cast<int>(amplitude.status), amplitude.value);
                }
            }
        }
    }
}

// Independent of Exchange: A(w) = h[r] + 2 sum_k h[r+k] cos(k w), summed in plain order.
double tap_response(std::span<const double> h, double omega) {
    const auto r = h.size() / 2;
    double sum = h[r];
    for (std::size_t k = 1; k <= r; ++k) sum += 2 * h[r + k] * std::cos(static_cast<double>(k) * omega);
    return sum;
}

// Measured worst reconstruction error over the exchange cases is about 3*N*eps, including this
// test's own uncompensated cosine sum.
double tap_tolerance(std::size_t taps) {
    return 16 * static_cast<double>(taps) * std::numeric_limits<double>::epsilon();
}

void tap_recovery() {
    for (const auto r : {std::size_t{0}, std::size_t{1}, std::size_t{7}, std::size_t{31}, std::size_t{127}}) {
        const auto N = 2*r+1;
        std::vector<remez::GridPoint> points(r+2);
        std::vector<double> h(N+1, sentinel.omega), samples(r+1);
        remez::Interpolator polynomial(points.size());
        check(polynomial.taps(h) == Status::invalid_nodes && std::ranges::all_of(h, [](double v) { return v == sentinel.omega; }),
              "unfitted tap recovery leaves output untouched");
        for (std::size_t degree = 0; degree <= r; ++degree) {
            for (std::size_t j = 0; j < points.size(); ++j) {
                const double x = std::cos(pi*static_cast<double>(j)/static_cast<double>(points.size()-1));
                points[j] = {0, x, chebyshev(degree, x), 1};
            }
            check(polynomial.fit(points).status == Status::ok, "known cosine-basis fit");
            std::ranges::fill(h, sentinel.omega);
            check(polynomial.taps(std::span{h}.first(N-1)) == Status::insufficient_storage &&
                  std::ranges::all_of(h, [](double v) { return v == sentinel.omega; }), "short tap output stays untouched");
            check(polynomial.taps(h) == Status::ok && h.back() == sentinel.omega, "cosine-basis tap recovery");
            for (std::size_t j = 0; j < N; ++j) {
                const double expected = degree == 0 ? (j == r ? 1.0 : 0.0) : (j == r-degree || j == r+degree ? 0.5 : 0.0);
                check(std::abs(h[j]-expected) <= 64*static_cast<double>(N)*std::numeric_limits<double>::epsilon(),
                      "cosine basis has the expected sparse symmetric impulse response");
            }
            for (std::size_t m = 0; m <= r; ++m) {
                const auto sample = polynomial.estimate(std::cos(2*pi*static_cast<double>(m)/static_cast<double>(N)));
                check(sample.status == Status::ok, "direct-transform reference sample");
                samples[m] = sample.value;
            }
            for (std::size_t k = 0; k <= r; ++k) {
                double reference = samples[0], magnitude = std::abs(reference);
                for (std::size_t m = 1; m <= r; ++m) {
                    const auto q = m*k % N;
                    const double term = 2*samples[m]*std::cos(2*pi*static_cast<double>(std::min(q,N-q))/static_cast<double>(N));
                    reference += term;
                    magnitude += std::abs(term);
                }
                const double bound = 128*std::numeric_limits<double>::epsilon()*magnitude/static_cast<double>(N);
                check(std::abs(h[r+k]-reference/static_cast<double>(N)) <= bound,
                      "rotation transform agrees with independent direct cosine sums");
            }
        }
        auto moved = std::move(polynomial);
        check(moved.taps(h) == Status::ok && polynomial.taps(h) == Status::invalid_nodes, "tap recovery respects moved ownership");
        points[0].weight = 0;
        check(moved.fit(points).status == Status::invalid_weights && moved.taps(h) == Status::invalid_nodes,
              "rejected fit invalidates recovered-tap access");
    }
    std::array<remez::GridPoint, 3> points{{{0,1,0,1},{0,0,0,1},{0,-1,0,1}}};
    remez::Interpolator polynomial(3);
    std::array<double, 3> h{};
    for (const double gain : {0.0, std::numeric_limits<double>::denorm_min(), 0.95*std::numeric_limits<double>::max()}) {
        for (auto& point : points) point.desired = gain;
        check(polynomial.fit(points).status == Status::ok && polynomial.taps(h) == Status::ok, "constant tap recovery across exponent range");
        if (gain == 0.0 || gain == std::numeric_limits<double>::denorm_min()) {
            check(h[0] == 0 && h[1] == gain && h[2] == 0, "zero and subnormal constant impulses");
        } else {
            check(std::abs(h[1]/gain-1) <= tolerance && std::abs(h[0]/gain) <= tolerance && h[0] == h[2],
                  "large finite constant does not overflow the inverse sum");
        }
    }
    points = {{{0,0.9,0.95*std::numeric_limits<double>::max(),1},{0,0,0,1},{0,-0.9,-0.95*std::numeric_limits<double>::max(),1}}};
    check(polynomial.fit(points).status == Status::ok && polynomial.taps(h) == Status::numerical_failure,
          "unrepresentable amplitude samples reject tap recovery");
}

void exchange_reference() {
    // SciPy 1.18.1: grid_density=512, maxiter=128, fs=2, weights={1,10}.
    struct Reference { std::size_t taps; double pass, stop; std::array<double, 6> values; };
    const std::array references{
        Reference{15, 0.20000000000000001, 0.29999999999999999, {1.369204595481571, 1.2240989852775044, 0.63079540451842875, 0.036920459548157314, -0.0067866262183710335, 0.03692045954815687}},
        Reference{63, 0.23999999999999999, 0.26000000000000001, {1.3441013601771825, 0.85340865793979082, 0.65589863982281948, 0.034410136017717363, -0.016367052366606641, 0.034410136017717807}},
        Reference{255, 0.23999999999999999, 0.26000000000000001, {1.0107264164727199, 0.98940766692846405, 0.98927358352728401, 0.0010726416472719436, -0.0007662460708280916, 0.001072641647274164}}
    };
    for (const auto& reference : references) {
        remez::Exchange engine({reference.taps, reference.pass*pi, reference.stop*pi, {1,10}, 16});
        const auto result = engine.run();
        check(result.status == Status::ok, "reference exchange converges");
        std::vector<double> h(reference.taps);
        check(engine.taps(h) == Status::ok, "reference tap recovery");
        if (reference.taps == 15) {
            constexpr std::array expected{-0.028458770784620859, -0.015306822035826059, 0.0062897668149520872, 0.05226038713404093, 0.11707879787873703, 0.18565661178492684, 0.23816124007428527, 0.25784217374858048};
            for (std::size_t j = 0; j < expected.size(); ++j)
                check(std::abs(h[j]-expected[j]) <= 2e-6*result.ripple, "independent SciPy tap coefficients");
        }
        const std::array frequencies{0.0, reference.pass/2, reference.pass, reference.stop, (reference.stop+1)/2, 1.0};
        for (std::size_t k = 0; k < frequencies.size(); ++k) {
            const auto value = engine.evaluate(std::cos(pi*frequencies[k]));
            const double weight = k < 3 ? 1.0 : 10.0;
            check(value.status == Status::ok && weight*std::abs(value.value-reference.values[k]) <= 3e-5*result.ripple+1e-13,
                  "independent SciPy minimax response");
            check(weight*std::abs(tap_response(h, pi*frequencies[k])-reference.values[k]) <=
                  3e-5*result.ripple+1e-13+weight*tap_tolerance(reference.taps), "taps reproduce SciPy minimax response");
        }
    }
}

void exchange() {
    for (const auto taps : {std::size_t{1}, std::size_t{3}, std::size_t{5}, std::size_t{15}, std::size_t{31}, std::size_t{63}, std::size_t{65}, std::size_t{127}, std::size_t{129}, std::size_t{255}}) {
      for (const auto edges : {std::pair{0.2, 0.3}, std::pair{0.24, 0.26}, std::pair{0.01, 0.03}, std::pair{0.97, 0.99}}) {
        const Lowpass spec{taps, edges.first*pi, edges.second*pi, {1, 10}, 32};
        remez::Exchange engine(spec);
        std::vector<double> h(taps + 1, sentinel.omega);
        check(engine.evaluate(0).status == Status::invalid_nodes && engine.extrema().empty() &&
              engine.taps(h) == Status::invalid_nodes, "unfitted exchange engine");
        const auto result = engine.run();
        if (taps == 255 && edges.first == 0.2 && result.status != Status::ok) {
            check((result.status == Status::numerical_failure || result.status == Status::nonconvergence) &&
                  engine.extrema().empty() && engine.evaluate(0).status == Status::invalid_nodes,
                  "unresolved high-order exchange exposes no polynomial");
            continue;
        }
        check(result.status == Status::ok && result.iterations <= 64 && result.flatness <= 1e-6, "exchange converges");
        const auto extrema = engine.extrema();
        check(extrema.size() == taps/2+2, "exchange preserves extremum count");
        // evaluate() rounds A to binary64, so W*(D-A) cannot resolve below a few W*eps*|A|.
        constexpr double eps = std::numeric_limits<double>::epsilon();
        const double rounding = 16 * eps * std::ranges::max(spec.weights);
        double previous = 0;
        for (const auto& point : extrema) {
            const auto value = engine.evaluate(point.x);
            const double error = point.weight*(point.desired-value.value);
            check(value.status == Status::ok && std::abs(std::abs(error)-result.ripple) <= 1e-10*result.ripple + rounding,
                  "exchange fitted extrema have equal ripple");
            check(previous == 0 || std::signbit(previous) != std::signbit(error), "exchange signs alternate");
            previous = error;
        }
        check(engine.taps(std::span{h}.first(taps - 1)) == Status::insufficient_storage, "short tap buffer rejected");
        check(engine.taps(h) == Status::ok && h.back() == sentinel.omega, "tap recovery leaves extra storage untouched");
        remez::Interpolator reference(extrema.size());
        std::vector<double> direct(taps);
        check(reference.fit(extrema).status == Status::ok && reference.taps(direct) == Status::ok, "independent tap recovery from final fit");
        for (std::size_t k = 0; k < taps; ++k)
            check(std::abs(h[k]-direct[k]) <= tap_tolerance(taps), "cached coefficients agree with direct recovery");
        const auto filter = std::span<const double>{h}.first(taps);
        check(std::ranges::all_of(filter, [](double v) { return std::isfinite(v); }), "finite taps");
        for (std::size_t k = 0; k < taps; ++k)
            check(filter[k] == filter[taps - 1 - k], "Type-I tap symmetry");
        for (const double omega : {0.0, 0.123*pi, spec.pass_edge, spec.stop_edge, 0.731*pi, pi}) {
            std::complex<double> response{};
            for (std::size_t j = 0; j < taps; ++j) response += filter[j]*std::polar(1.0, -omega*static_cast<double>(j));
            const auto amplitude = response*std::polar(1.0, omega*static_cast<double>(taps/2));
            const auto expected = engine.evaluate(std::cos(omega));
            check(expected.status == Status::ok && std::abs(amplitude.real()-expected.value) <= tap_tolerance(taps) &&
                  std::abs(amplitude.imag()) <= tap_tolerance(taps), "full complex FIR response has the expected amplitude and group delay");
        }
        // Independent band sampling is denser than the design grid; the tap response must match the
        // internal amplitude and stay within the same minimax bound.
        double dense_max = 0, tap_max = 0;
        for (std::size_t band = 0; band < 2; ++band) {
            const double lower = band == 0 ? 0 : spec.stop_edge;
            const double upper = band == 0 ? spec.pass_edge : pi;
            const double desired = band == 0 ? 1.0 : 0.0;
            for (std::size_t j = 0; j <= 8192; ++j) {
                const double omega = std::lerp(lower, upper, static_cast<double>(j)/8192);
                const auto value = engine.evaluate(std::cos(omega));
                check(value.status == Status::ok, "dense independent exchange evaluation");
                const double response = tap_response(filter, omega);
                check(std::abs(response - value.value) <= tap_tolerance(taps), "tap response matches internal amplitude");
                dense_max = std::max(dense_max, spec.weights[band]*std::abs(desired-value.value));
                tap_max = std::max(tap_max, spec.weights[band]*std::abs(desired-response));
            }
        }
        check(dense_max <= result.max_error*(1.0+2e-6) + rounding, "no large missed off-grid peak");
        check(tap_max <= result.max_error*(1.0+2e-6) + rounding + spec.weights[1]*tap_tolerance(taps),
              "tap response meets the minimax bound");
        if (taps == 1) {
            check(std::abs(engine.evaluate(0).value-1.0/11) <= tolerance && std::abs(result.ripple-10.0/11) <= tolerance,
                  "constant minimax has analytic weighted solution");
        }
        const auto cached = h;
        const auto repeated = engine.run();
        check(repeated.status == Status::ok && repeated.iterations == result.iterations && repeated.ripple == result.ripple,
              "exchange restart is deterministic");
        std::ranges::fill(h, sentinel.omega);
        check(engine.taps(h) == Status::ok && std::ranges::equal(h, cached), "restart rebuilds identical taps and preserves extra output");
        auto transferred = std::move(engine);
        check(engine.taps(h) == Status::invalid_nodes && transferred.taps(h) == Status::ok, "tap cache moves with its fitted polynomial");
        engine = std::move(transferred);
        check(transferred.taps(h) == Status::invalid_nodes && engine.taps(h) == Status::ok, "tap cache move assignment");
        check(engine.run(0).status == Status::nonconvergence && engine.extrema().empty() && engine.evaluate(0).status == Status::invalid_nodes &&
              engine.taps(h) == Status::invalid_nodes, "iteration exhaustion invalidates exchange output");
        check(engine.run(64, 0).status == Status::invalid_tolerance, "invalid exchange tolerance");
      }
    }
    const Lowpass spec{31, 0.24*pi, 0.26*pi, {1, 10}};
    remez::Exchange baseline(spec);
    const auto expected = baseline.run();
    check(expected.status == Status::ok, "exchange baseline");
    for (const auto density : {std::size_t{8}, std::size_t{16}, std::size_t{32}}) {
        for (const double scale : {1e-200, 1.0, 1e200}) {
            auto scaled = spec; scaled.density = density; scaled.weights = {scale, 10*scale};
            remez::Exchange engine(scaled);
            const auto result = engine.run();
            check(result.status == Status::ok && std::abs(result.ripple/scale-expected.ripple) <= 2e-6*expected.ripple,
                  "grid-density and common-weight scaling preserve refined minimax");
            check(std::abs(engine.evaluate(0.5).value-baseline.evaluate(0.5).value) <= 2e-6, "exchange scaling preserves amplitude");
        }
    }
    remez::Exchange limited(spec);
    check(limited.run(1).status == Status::nonconvergence && limited.evaluate(0).status == Status::invalid_nodes,
          "positive iteration budget exhaustion");
    for (const double bad : {-1.0, 0.0, 1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(limited.run(64, bad).status == Status::invalid_tolerance, "invalid tolerances");
    check(limited.run().status == Status::ok, "exchange recovers after failure");
    auto moved = std::move(limited);
    check(moved.evaluate(0).status == Status::ok && limited.evaluate(0).status == Status::invalid_nodes &&
          limited.extrema().empty() && limited.run().status == Status::invalid_nodes, "exchange move transfers ownership");
    remez::Exchange assigned(spec);
    assigned = std::move(moved);
    check(assigned.evaluate(0).status == Status::ok && moved.evaluate(0).status == Status::invalid_nodes, "exchange move assignment");
    remez::Exchange invalid({0, 0.2*pi, 0.3*pi});
    check(invalid.run().status == Status::invalid_taps, "invalid exchange setup");
    remez::Exchange capped(spec, 2);
    check(capped.run().status == Status::grid_limit, "exchange grid budget");
    auto extreme = spec; extreme.weights = {std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::max()};
    remez::Exchange unresolved(extreme);
    check(unresolved.run().status == Status::numerical_failure && unresolved.extrema().empty(), "unrepresentable exchange weights");
}


void peak_selection() {
    for (std::size_t K = 2; K <= 10; ++K) for (std::size_t n = 2; n <= K; ++n) for (std::size_t trial = 0; trial < 32; ++trial) {
        std::array<remez::detail::Peak, 10> peaks{};
        std::array<double, 10> scratch{};
        for (std::size_t k = 0; k < K; ++k) {
            const double magnitude = static_cast<double>(1+(k*k*7+trial*11+k*trial)%9);
            peaks[k] = {k, static_cast<double>(k), k%2 == 0 ? magnitude : -magnitude};
        }
        double best = 0.0;
        for (unsigned mask = 0; mask < (1u << K); ++mask) {
            std::size_t count = 0;
            bool sign = false, valid = true;
            double minimum = 10.0;
            for (std::size_t k = 0; k < K; ++k) if ((mask & (1u << k)) != 0) {
                if (count != 0 && sign == std::signbit(peaks[k].error)) valid = false;
                sign = std::signbit(peaks[k].error); ++count;
                minimum = std::min(minimum, std::abs(peaks[k].error));
            }
            if (valid && count == n) best = std::max(best, minimum);
        }
        check(remez::detail::select_peaks(std::span{peaks}.first(K), n, scratch) == n, "general peak selector returns exact count");
        double minimum = 10.0;
        for (std::size_t k = 0; k < n; ++k) {
            minimum = std::min(minimum, std::abs(peaks[k].error));
            if (k != 0) check(peaks[k-1].index < peaks[k].index && std::signbit(peaks[k-1].error) != std::signbit(peaks[k].error),
                              "selected peaks are ordered and alternating");
        }
        check(minimum == best, "peak selector matches exhaustive max-min subset oracle");
    }
}

void general_recovery() {
    for (std::size_t d : {0u,1u,7u,31u,127u}) for (auto symmetry : {remez::Symmetry::symmetric,remez::Symmetry::antisymmetric})
        for (bool even : {false,true}) {
            const auto N = symmetry == remez::Symmetry::symmetric ? 2*d+(even?2:1) : 2*d+(even?2:3);
            std::vector<remez::GridPoint> points(d+2);
            std::vector<double> h(N+1,sentinel.omega);
            remez::Interpolator polynomial(points.size());
            for (auto degree : {std::size_t{0},d/2,d}) {
                for (std::size_t j=0;j<points.size();++j) {
                    const double x=std::cos(pi*static_cast<double>(j)/static_cast<double>(points.size()-1));
                    points[j]={0,x,chebyshev(degree,x),1};
                }
                check(polynomial.fit(points).status == Status::ok,"general basis polynomial fit");
                check(polynomial.taps(std::span{h}.first(N-1),symmetry,N) == Status::insufficient_storage,"general direct recovery rejects short buffer");
                check(polynomial.taps(h,symmetry,N) == Status::ok && h.back() == sentinel.omega,"general direct basis recovery");
                for(std::size_t j=0;j<=256;++j) {
                    const double w=pi*static_cast<double>(j)/256;
                    double q=1;
                    if(symmetry==remez::Symmetry::symmetric && even) q=std::cos(w/2);
                    if(symmetry==remez::Symmetry::antisymmetric) q=even?std::sin(w/2):std::sin(w);
                    std::complex<double> value{};
                    for(std::size_t k=0;k<N;++k) value+=h[k]*std::polar(1.0,-w*static_cast<double>(k));
                    value*=std::polar(1.0,w*static_cast<double>(N-1)/2);
                    const std::complex<double> expected=symmetry==remez::Symmetry::symmetric?std::complex<double>{q*chebyshev(degree,std::cos(w)),0}:std::complex<double>{0,q*chebyshev(degree,std::cos(w))};
                    check(std::abs(value-expected)<=64*static_cast<double>(N)*std::numeric_limits<double>::epsilon(),"direct basis recovery matches analytic complex response for every form");
                }
            }
            check(polynomial.taps(h,symmetry,N+2)==Status::invalid_taps,"recovery degree mismatch rejected");
        }
}

void generalized() {
    using remez::Band; using remez::Specification; using remez::Symmetry; using remez::Response;
    const std::array lp{Band{0,.2*pi,1,1}, Band{.3*pi,pi,0,10}};
    const std::array hp{Band{0,.2*pi,0,10}, Band{.3*pi,pi,1,1}};
    const std::array bp{Band{0,.2*pi,0,10}, Band{.3*pi,.6*pi,1,1}, Band{.7*pi,pi,0,10}};
    const std::array bs{Band{0,.2*pi,1,1}, Band{.3*pi,.6*pi,0,10}, Band{.7*pi,pi,1,1}};
    const std::array multi{Band{0,.1*pi,.2,1}, Band{.2*pi,.3*pi,1,2}, Band{.4*pi,.5*pi,-.5,3}, Band{.6*pi,.7*pi,.8,4}, Band{.8*pi,pi,0,10}};
    const std::array hilbert{Band{.05*pi,.95*pi,-1,1}};
    const std::array diff3{Band{0,.8*pi,1,1}, Band{.9*pi,pi,0,10}};
    const std::array diff4{Band{0,pi,1,1}};
    // SciPy 1.18.1, density 512, maxiter 128. Differentiator slopes/weights converted from radians.
    constexpr std::array<double, 31> fixture0{-0.020934364437286704, 0.022393799499398452, 0.016594211795943312, 0.0069770856260239307, -0.0079594698470608159, -0.021778493670776066, -0.0241603603971648, -0.0086784968279752829, 0.020273503272464462, 0.046359155709411584, 0.048005945387335734, 0.010186383897047298, -0.06518463055894691, -0.15676028186370189, -0.23169478062173821, 0.73941489514358882, -0.23169478062173821, -0.15676028186370189, -0.06518463055894691, 0.010186383897047298, 0.048005945387335734, 0.046359155709411584, 0.020273503272464462, -0.0086784968279752829, -0.0241603603971648, -0.021778493670776066, -0.0079594698470608159, 0.0069770856260239307, 0.016594211795943312, 0.022393799499398452, -0.020934364437286704};
    constexpr std::array<double, 31> fixture1{-0.0085573186696835373, 0.0075097831224829516, 0.0183875770295574, 0.00072008855702256915, 0.0024710132944305596, -0.0041691981718649188, -0.053462928553118648, -0.018966368705298265, 0.065850503634368493, 0.022913613267214243, 0.008959854358972133, 0.084518681791226125, -0.087986916240881799, -0.27766409600142217, 0.054338215146355291, 0.37910873669969486, 0.054338215146355291, -0.27766409600142217, -0.087986916240881799, 0.084518681791226125, 0.008959854358972133, 0.022913613267214243, 0.065850503634368493, -0.018966368705298265, -0.053462928553118648, -0.0041691981718649188, 0.0024710132944305596, 0.00072008855702256915, 0.0183875770295574, 0.0075097831224829516, -0.0085573186696835373};
    constexpr std::array<double, 31> fixture2{-0.015735819024373918, 0.014835479339145067, -0.015260471671871252, 0.019044936255938588, 0.039794869035146015, 0.0035974877830199621, 0.019290082378912575, 0.019885688352301756, -0.067293868436692147, -0.041499868787664602, 0.014698710399737186, -0.057404169672093788, 0.084683458717752055, 0.29014327849504723, -0.060176961398609868, 0.57860163881005544, -0.060176961398609868, 0.29014327849504723, 0.084683458717752055, -0.057404169672093788, 0.014698710399737186, -0.041499868787664602, -0.067293868436692147, 0.019885688352301756, 0.019290082378912575, 0.0035974877830199621, 0.039794869035146015, 0.019044936255938588, -0.015260471671871252, 0.014835479339145067, -0.015735819024373918};
    constexpr std::array<double, 15> fixture3{0.085736008332173402, 0.13333084976016513, -0.15160227102350737, -0.16960606175362483, 0.1088321176387604, 0.055332501880287696, 0.1106009511449377, 0.31504872400767725, 0.1106009511449377, 0.055332501880287696, 0.1088321176387604, -0.16960606175362483, -0.15160227102350737, 0.13333084976016513, 0.085736008332173402};
    constexpr std::array<double, 32> fixture4{-0.0064902982704557882, -0.0049652330188863069, -0.00088958916661798166, 0.0082690213201461298, 0.019417750785666502, 0.026365055544241726, 0.022656393339840427, 0.0057963362663038627, -0.019430374432658718, -0.040917773802925161, -0.043448047379980842, -0.015798906135882074, 0.042568547742062521, 0.11891573199913763, 0.19071173716754025, 0.23418328226423002, 0.23418328226423002, 0.19071173716754025, 0.11891573199913763, 0.042568547742062521, -0.015798906135882074, -0.043448047379980842, -0.040917773802925161, -0.019430374432658718, 0.0057963362663038627, 0.022656393339840427, 0.026365055544241726, 0.019417750785666502, 0.0082690213201461298, -0.00088958916661798166, -0.0049652330188863069, -0.0064902982704557882};
    constexpr std::array<double, 31> fixture5{-0.03192367836850326, 1.6506630344403635e-15, -0.02604586500491711, -7.6878209964369998e-16, -0.037320916902723426, 7.4163962901512533e-17, -0.053092776177022238, -4.5948402104485316e-16, -0.076688664348758373, 9.6705760467325067e-17, -0.1168448104037082, -6.7313205849769761e-17, -0.20579817845397408, -2.7475916837786641e-16, -0.63446321100815162, 0, 0.63446321100815162, 2.7475916837786641e-16, 0.20579817845397408, 6.7313205849769761e-17, 0.1168448104037082, -9.6705760467325067e-17, 0.076688664348758373, 4.5948402104485316e-16, 0.053092776177022238, -7.4163962901512533e-17, 0.037320916902723426, 7.6878209964369998e-16, 0.02604586500491711, -1.6506630344403635e-15, 0.03192367836850326};
    constexpr std::array<double, 31> fixture6{-0.045934609631201745, 0.020702179097299284, 0.056590005659056572, -0.060766779698937956, 0.051189643736576695, -0.029151469258650732, -0.0059182678567156338, 0.050238177362022435, -0.094680944827040009, 0.12562421936989346, -0.12662444644123005, 0.079487824660651671, 0.038808272576328934, -0.27797780682306306, 0.88060998519685807, 0, -0.88060998519685807, 0.27797780682306306, -0.038808272576328934, -0.079487824660651671, 0.12662444644123005, -0.12562421936989346, 0.094680944827040009, -0.050238177362022435, 0.0059182678567156338, 0.029151469258650732, -0.051189643736576695, 0.060766779698937956, -0.056590005659056572, -0.020702179097299284, 0.045934609631201745};
    constexpr std::array<double, 32> fixture7{-0.026493621260073, -0.0057653007898626069, -0.016568074278603753, -0.0095335770782756027, -0.02247200528666584, -0.014961114020507524, -0.030313758532595853, -0.022995444428882961, -0.041419851777148436, -0.035735226446599294, -0.059119382195648218, -0.059012422412988874, -0.09448672865361657, -0.11737541546315011, -0.21800661329690929, -0.62866354497609622, 0.62866354497609622, 0.21800661329690929, 0.11737541546315011, 0.09448672865361657, 0.059012422412988874, 0.059119382195648218, 0.035735226446599294, 0.041419851777148436, 0.022995444428882961, 0.030313758532595853, 0.014961114020507524, 0.02247200528666584, 0.0095335770782756027, 0.016568074278603753, 0.0057653007898626069, 0.026493621260073};
    constexpr std::array<double, 32> fixture8{-0.0039359603211662582, 0.0053802639494991153, -0.0026663524159915072, 0.0025047985479988199, -0.0027306590883832815, 0.0031458377569417214, -0.0037493377129640621, 0.0046049611849922466, -0.0058422348236550972, 0.0077058498842837884, -0.010685990088053862, 0.015875971655971054, -0.026136741498605901, 0.051078527312964903, -0.14161788085897356, 1.2733853224358698, -1.2733853224358698, 0.14161788085897356, -0.051078527312964903, 0.026136741498605901, -0.015875971655971054, 0.010685990088053862, -0.0077058498842837884, 0.0058422348236550972, -0.0046049611849922466, 0.0037493377129640621, -0.0031458377569417214, 0.0027306590883832815, -0.0025047985479988199, 0.0026663524159915072, -0.0053802639494991153, 0.0039359603211662582};
    struct Reference { std::string_view name; std::size_t taps; std::span<const double> h; };
    const std::array references{
        Reference{"hp", 31, fixture0},
        Reference{"bp", 31, fixture1},
        Reference{"bs", 31, fixture2},
        Reference{"multi", 15, fixture3},
        Reference{"II-lp", 32, fixture4},
        Reference{"III-hilbert", 31, fixture5},
        Reference{"III-diff", 31, fixture6},
        Reference{"IV-hilbert", 32, fixture7},
        Reference{"IV-diff", 32, fixture8}
    };
    const auto test = [&](std::string_view name, Specification spec, double flatness = 1e-6) {
        const auto size = remez::grid_size(spec);
        check(size.status == Status::ok, "generic specification sizing");
        const auto grid = remez::make_grid(spec);
        check(grid.status == Status::ok && grid.points.size() == size.points && grid.bands.size() == spec.bands.size(), "generic grid owns exact ranges");
        for (std::size_t b = 0; b < spec.bands.size(); ++b) {
            const auto view = grid.band(b);
            check(!view.empty(), "nonempty generic band view");
            for (const auto& p : view) check(p.omega >= spec.bands[b].lower && p.omega <= spec.bands[b].upper && p.weight > 0 && std::isfinite(p.desired),
                                           "transformed points stay in owning band");
        }
        check(grid.band(spec.bands.size()).empty(), "invalid band index returns empty view");
        remez::Exchange engine(spec);
        const auto result = engine.run(128, flatness);
        if (result.status != Status::ok) std::printf("generic case %.*s N=%zu status=%d\n", static_cast<int>(name.size()), name.data(), spec.taps, static_cast<int>(result.status));
        check(result.status == Status::ok && result.flatness <= flatness, "generic exchange converges");
        std::vector<double> h(spec.taps+1, sentinel.omega);
        check(engine.taps(std::span{h}.first(spec.taps-1)) == Status::insufficient_storage && h[0] == sentinel.omega, "generic short tap buffer unchanged");
        check(engine.taps(h) == Status::ok && h.back() == sentinel.omega, "generic tap recovery preserves extra output");
        remez::Interpolator direct_model(size.degree+2);
        std::vector<double> direct(spec.taps+1,sentinel.omega);
        check(direct_model.fit(engine.extrema()).status == Status::ok && direct_model.taps(direct,spec.symmetry,spec.taps) == Status::ok && direct.back() == sentinel.omega,
              "direct recovery independently reconstructs each form from final extrema");
        for (std::size_t k=0;k<spec.taps;++k) check(std::abs(h[k]-direct[k]) <= 64*static_cast<double>(spec.taps)*std::numeric_limits<double>::epsilon(), "cached mapping agrees with direct inverse DFT");
        const bool anti = spec.symmetry == Symmetry::antisymmetric;
        for (std::size_t k = 0; k < spec.taps; ++k)
            check(std::isfinite(h[k]) && h[k] == (anti ? -h[spec.taps-1-k] : h[spec.taps-1-k]), "exact generic tap symmetry");
        if (anti && spec.taps%2 != 0) check(h[spec.taps/2] == 0.0, "Type III centre is exactly zero");
        if (anti) check(engine.amplitude(0).value == 0.0, "antisymmetric DC is exactly zero");
        if ((!anti && spec.taps%2 == 0) || (anti && spec.taps%2 != 0)) check(engine.amplitude(pi).value == 0.0, "forced Nyquist zero");
        const double alpha = static_cast<double>(spec.taps-1)/2;
        double slope = 0.0;
        for (std::size_t k = 0; k < spec.taps; ++k) slope -= static_cast<double>(k)*h[k];
        double peak = 0.0;
        for (const auto& band : spec.bands) for (std::size_t j = 0; j <= 1024; ++j) {
            const double w = std::lerp(band.lower, band.upper, static_cast<double>(j)/1024);
            std::complex<double> H{};
            for (std::size_t k = 0; k < spec.taps; ++k) H += h[k]*std::polar(1.0,-w*static_cast<double>(k));
            H *= std::polar(1.0,alpha*w);
            const double amplitude = anti ? H.imag() : H.real();
            const auto expected = engine.amplitude(w);
            check(expected.status == Status::ok && std::abs(amplitude-expected.value) <= 64*static_cast<double>(spec.taps)*std::numeric_limits<double>::epsilon(),
                  "full complex response matches physical amplitude");
            check(std::abs(anti ? H.real() : H.imag()) <= 64*static_cast<double>(spec.taps)*std::numeric_limits<double>::epsilon(), "generic phase and group delay");
            const bool relative = spec.response == Response::differentiator && band.target != 0;
            const double error = relative ? band.weight*(band.target-(w == 0 ? slope : amplitude/w))/std::abs(band.target) : band.weight*(band.target-amplitude);
            peak = std::max(peak,std::abs(error));
        }
        check(peak <= result.max_error*(1+3*flatness)+1e-10, "independent dense physical weighted error including DC slope limit");
        const auto reference = std::ranges::find_if(references, [&](const Reference& r) { return r.name == name && r.taps == spec.taps; });
        if (reference != references.end()) for (std::size_t k = 0; k < spec.taps; ++k)
            check(std::abs(h[k]-reference->h[k]) <= 3e-5*result.ripple+1e-12, "independent generalized SciPy coefficients");
        const auto saved = h;
        check(engine.run(128, flatness).status == Status::ok && engine.taps(h) == Status::ok && h == saved, "general restart deterministic");
        auto moved = std::move(engine);
        check(engine.taps(h) == Status::invalid_nodes && moved.taps(h) == Status::ok, "general move transfers tap cache");
        check(moved.run(0).status == Status::nonconvergence && moved.amplitude(.5).status == Status::invalid_nodes && moved.taps(h) == Status::invalid_nodes, "general failure invalidates all outputs");
    };
    for (auto N : {1u,3u,5u,15u,31u,63u}) {
        test("hp", {N,hp,Symmetry::symmetric,Response::constant,32});
        test("bp", {N,bp,Symmetry::symmetric,Response::constant,32});
        test("bs", {N,bs,Symmetry::symmetric,Response::constant,32});
        test("multi", {N,multi,Symmetry::symmetric,Response::constant,32});
    }
    for (auto N : {2u,4u,16u,32u,64u}) {
        test("II-lp", {N,lp,Symmetry::symmetric,Response::constant,32});
        test("II-bp", {N,bp,Symmetry::symmetric,Response::constant,32});
        test("IV-hilbert", {N,hilbert,Symmetry::antisymmetric,Response::constant,32});
        test("IV-diff", {N,diff4,Symmetry::antisymmetric,Response::differentiator,32});
    }
    for (auto N : {3u,5u,15u,31u,63u}) {
        test("III-hilbert", {N,hilbert,Symmetry::antisymmetric,Response::constant,32});
        test("III-diff", {N,diff3,Symmetry::antisymmetric,Response::differentiator,32});
    }
    // Ripples near 1e-7: scanning in the transformed domain and charging the transformed/physical
    // rounding to the reported flatness (instead of a 100x margin) lets these converge.
    test("II-lp", {192,lp,Symmetry::symmetric,Response::constant,16});
    test("III-hilbert", {191,hilbert,Symmetry::antisymmetric,Response::constant,16});
    test("IV-hilbert", {192,hilbert,Symmetry::antisymmetric,Response::constant,16});
    // Ripple 5e-9: flatness about 2e-6 is binary64-limited, so the tolerance is 1e-5.
    test("III-hilbert", {223,hilbert,Symmetry::antisymmetric,Response::constant,16}, 1e-5);
    // Equispaced seeds start these at ripples of 1e-11..1e-20 and stall; the equilibrium-measure
    // seed reaches them in 5-12 iterations. Ripples 7e-9..3e-12; the last ones are binary64-limited.
    test("lp", {223,lp,Symmetry::symmetric,Response::constant,16});
    test("bp", {255,bp,Symmetry::symmetric,Response::constant,16});
    test("lp", {319,lp,Symmetry::symmetric,Response::constant,16});
    test("II-lp", {224,lp,Symmetry::symmetric,Response::constant,16});
    test("IV-hilbert", {224,hilbert,Symmetry::antisymmetric,Response::constant,16}, 1e-5);
    test("III-diff", {223,diff3,Symmetry::antisymmetric,Response::differentiator,16});
    test("III-diff", {255,diff3,Symmetry::antisymmetric,Response::differentiator,16}, 1e-5);
    // SciPy 1.18.1 (density 16) stops at max weighted errors 9.125001e-9 (223) and 8.738776e-9 (224),
    // measured on 40001 points per band. A minimax design cannot do worse; a 10% gap would mean
    // converging to the wrong alternation set.
    for (const auto [taps, scipy] : {std::pair{223u, 9.125001e-9}, std::pair{224u, 8.738776e-9}}) {
        remez::Exchange engine(Specification{taps,lp,Symmetry::symmetric,Response::constant,16});
        const auto result = engine.run(128, 1e-6);
        check(result.status == Status::ok && result.max_error <= scipy && result.ripple >= 0.9*scipy, "long design no worse than SciPy");
    }
    auto input = lp;
    remez::Exchange owned(Specification{31,input});
    check(owned.run().status == Status::ok, "generic engine initially fitted");
    std::array<double,31> before{}, after{};
    check(owned.taps(before) == Status::ok, "owned taps before mutation");
    input[0].target = -100; input[1].weight = 0;
    check(owned.run().status == Status::ok && owned.taps(after) == Status::ok && before == after, "constructor owns band parameters independently of input");
    remez::Exchange legacy(Lowpass{31,.2*pi,.3*pi,{1,10}});
    check(legacy.run().status == Status::ok && legacy.taps(after) == Status::ok && before == after, "lowpass adapter is identical to general specification");
    const std::array same{Band{0,.2*pi,.95*std::numeric_limits<double>::max(),1},Band{.3*pi,pi,.95*std::numeric_limits<double>::max(),10}};
    remez::Exchange exact(Specification{31,same});
    check(exact.run().status == Status::ok && exact.taps(after) == Status::ok && after[15] == same[0].target && exact.extrema().empty(), "analytic exact constant response avoids overflow and fictitious alternation");
    for (std::size_t k = 0; k < after.size(); ++k) if (k != 15) check(after[k] == 0.0, "exact constant has one impulse tap");
    const double huge = .95*std::numeric_limits<double>::max();
    const std::array opposite{Band{0,.2*pi,-huge,1},Band{.3*pi,pi,huge,1}};
    remez::Exchange scalar(Specification{1,opposite});
    std::array<double,1> impulse{};
    const auto scalar_result = scalar.run();
    check(scalar_result.status == Status::ok && scalar_result.ripple == huge && scalar.taps(impulse) == Status::ok && impulse[0] == 0.0,
          "scalar minimax handles opposite large targets without overflowing their difference");
    for (auto sym : {Symmetry::symmetric,Symmetry::antisymmetric}) for (auto N : {30u,31u}) {
        const std::array zero{Band{0,pi,0,1}};
        remez::Exchange e(Specification{N,zero,sym});
        check(e.run().status == Status::ok, "exact zero supported in every form");
        std::vector<double> h(N,1);
        check(e.taps(h) == Status::ok && std::ranges::all_of(h,[](double a){return a==0.0;}), "exact zero impulse response");
    }
    for (double scale : {1e-100,1e100}) {
        auto band = hilbert; band[0].weight *= scale;
        remez::Exchange e(Specification{31,band,Symmetry::antisymmetric});
        const auto r = e.run();
        check(r.status == Status::ok && e.taps(after) == Status::ok, "general weight scaling preserves convergence");
        for (std::size_t k=0;k<31;++k) check(std::abs(after[k]-fixture5[k]) <= 3e-5*r.ripple/scale+1e-12, "general weight scaling preserves taps");
    }
    remez::Exchange tiny(Specification{2,diff4,Symmetry::antisymmetric,Response::differentiator});
    check(tiny.run().status == Status::ok && tiny.amplitude(std::numeric_limits<double>::denorm_min()).value == std::numeric_limits<double>::denorm_min(),
          "half-angle underflow does not erase representable physical amplitude");
    const std::array wrong{Band{0,pi,1,1}};
    check(remez::grid_size(Specification{32,wrong}).status == Status::incompatible_endpoint, "Type II rejects nonzero forced endpoint");
    check(remez::grid_size(Specification{31,wrong,Symmetry::antisymmetric}).status == Status::incompatible_endpoint, "Hilbert DC endpoint incompatible");
    check(remez::grid_size(Specification{31,wrong,Symmetry::antisymmetric,Response::differentiator}).status == Status::incompatible_endpoint, "Type III differentiator Nyquist incompatible");
    check(remez::grid_size(Specification{31,wrong,Symmetry::symmetric,Response::differentiator}).status == Status::invalid_specification, "differentiator requires antisymmetry");
    check(remez::grid_size(Specification{1,hilbert,Symmetry::antisymmetric}).status == Status::invalid_taps, "minimum antisymmetric length");
    check(remez::grid_size(Specification{31,{}}).status == Status::invalid_bands, "empty band specification");
    check(remez::grid_size(Specification{std::numeric_limits<std::size_t>::max(),lp,Symmetry::symmetric,Response::constant,1},std::numeric_limits<std::size_t>::max()).status == Status::grid_limit,
          "interpolation storage byte bound rejected before allocation");
    input = lp; input[1].lower = input[0].upper;
    check(remez::grid_size(Specification{31,input}).status == Status::invalid_edges, "touching bands rejected");
    input = lp; input[0].target = std::numeric_limits<double>::infinity();
    check(remez::grid_size(Specification{31,input}).status == Status::invalid_target, "nonfinite target rejected");
}

#ifdef REMEZ_CRT_CHECKS
bool counting = false;
bool reject = false;
std::size_t allocations = 0, allocation_bytes = 0, accepted_allocations = 0;

int __cdecl allocation_hook(int kind, void*, std::size_t bytes, int block, long,
                            const unsigned char*, int) {
    if (counting && (block & 0xffff) != _CRT_BLOCK && (kind == _HOOK_ALLOC || kind == _HOOK_REALLOC)) {
        ++allocations;
        allocation_bytes += bytes;
        return reject && allocations > accepted_allocations ? 0 : 1;
    }
    return 1;
}

bool outstanding(const _CrtMemState& before, const _CrtMemState& after) {
    _CrtMemState difference{};
    return _CrtMemDifference(&difference, &before, &after) != 0;
}

void memory_contract() {
    const Lowpass spec{31, pi / 4, pi / 2};
    std::array<remez::GridPoint, 194> storage{};
    const auto previous = _CrtSetAllocHook(allocation_hook);
    allocations = 0; counting = true;
    const auto first_size = remez::grid_size(spec);
    const auto first_fill = remez::fill_grid(spec, storage);
    counting = false;
    check(first_size.status == Status::ok && first_fill == Status::ok && allocations == 0,
          "first-use sizing/fill makes zero allocations");
    _CrtMemState before{}, after{};
    _CrtMemCheckpoint(&before);
    {
        allocations = 0; allocation_bytes = 0; counting = true;
        remez::Interpolator interpolator(3);
        counting = false;
        check(allocations == 1 && allocation_bytes == 120, "three-node workspace uses one 120-byte allocation");
        const std::array<remez::GridPoint, 3> extrema{{{0, 1, 1, 1}, {0, 0, 1, 1}, {0, -1, 0, 1}}};
        allocations = 0; counting = true;
        const auto first_fit = interpolator.fit(extrema);
        const auto first_value = interpolator.evaluate(0.25);
        const auto first_error = interpolator.weighted_error({0, 0.25, 1, 2});
        std::array<double, 3> taps{};
        const auto recovered = interpolator.taps(taps);
        counting = false;
        check(allocations == 0 && first_fit.status == Status::ok && first_value.status == Status::ok &&
              first_error.status == Status::ok && std::abs(first_error.value-0.25) <= tolerance &&
              recovered == Status::ok && std::abs(taps[0]-0.25) <= tolerance && std::abs(taps[1]-0.75) <= tolerance && taps[0] == taps[2] &&
              std::abs(first_fit.value + 0.25) <= tolerance && std::abs(first_value.value - 0.875) <= tolerance,
              "first-use interpolation fit/evaluation makes zero allocations");
        std::array<remez::Value, 33> values{};
        for (int repeat = 0; repeat < 100; ++repeat) {
            allocations = 0; counting = true;
            const auto fit = interpolator.fit(extrema);
            for (std::size_t j = 0; j < values.size(); ++j)
                values[j] = interpolator.evaluate(-1.0 + static_cast<double>(j) / 16);
            counting = false;
            check(allocations == 0 && fit.status == Status::ok &&
                  std::ranges::all_of(values, [](const auto& v) { return v.status == Status::ok; }),
                  "repeated interpolation fit/evaluation makes zero allocations");
        }
        for (const double gain : {std::numeric_limits<double>::denorm_min(), 0.95*std::numeric_limits<double>::max()}) {
            auto constant = extrema;
            for (auto& point : constant) point.desired = gain;
            allocations = 0; counting = true;
            const auto fit = interpolator.fit(constant);
            const auto status = interpolator.taps(taps);
            counting = false;
            check(allocations == 0 && fit.status == Status::ok && status == Status::ok,
                  "extreme-range tap normalization makes zero allocations");
        }
        auto invalid = extrema; invalid[0].weight = 0.0;
        allocations = 0; counting = true;
        const auto rejected = interpolator.fit(invalid);
        const auto stale = interpolator.evaluate(0.0);
        counting = false;
        check(allocations == 0 && rejected.status == Status::invalid_weights && stale.status == Status::invalid_nodes,
              "rejected interpolation fit makes zero allocations");
        invalid = extrema; invalid[0].weight = std::numeric_limits<double>::denorm_min();
        invalid[1].weight = std::numeric_limits<double>::max();
        allocations = 0; counting = true;
        const auto failure = interpolator.fit(invalid);
        counting = false;
        check(allocations == 0 && failure.status == Status::numerical_failure, "numerical interpolation failure makes zero allocations");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "interpolation workspace released on success and failure");
    _CrtMemCheckpoint(&before);
    {
        allocations = 0; allocation_bytes = 0; counting = true;
        remez::Interpolator interpolator(129);
        counting = false;
        check(allocations == 1 && allocation_bytes >= 5160 && allocation_bytes <= 5224,
              "129-node workspace: 5160-byte payload plus at most 64 bytes of STL alignment overhead");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "129-node workspace released");
    _CrtMemCheckpoint(&before);
    {
        allocations = 0; allocation_bytes = 0; counting = true;
        remez::Exchange engine(spec);
        counting = false;
        const std::array bands{remez::Band{0,spec.pass_edge,1,spec.weights[0]}, remez::Band{spec.stop_edge,pi,0,spec.weights[1]}};
        const auto size = remez::grid_size(remez::Specification{spec.taps,bands});
        const auto payload = 32*size.points + 48*bands.size() + 72*(size.degree+2) + 24*size.peak_capacity + 8*size.workspace;
        check(allocations == 6 && allocation_bytes >= payload && allocation_bytes <= payload + 384,
              "exchange uses six explicit setup buffers");
        std::array<double, 31> taps{};
        for (int repeat = 0; repeat < 10; ++repeat) {
            allocations = 0; counting = true;
            const auto result = engine.run();
            const auto value = engine.evaluate(0.25);
            const auto extrema = engine.extrema();
            const auto recovered = engine.taps(taps);
            counting = false;
            check(allocations == 0 && result.status == Status::ok && value.status == Status::ok && !extrema.empty() &&
                  recovered == Status::ok, "first-use and repeated exchange and tap recovery make zero allocations");
        }
        allocations = 0; counting = true;
        const auto exhausted = engine.run(1);
        const auto invalid = engine.run(64, 0);
        counting = false;
        check(allocations == 0 && exhausted.status == Status::nonconvergence && invalid.status == Status::invalid_tolerance,
              "exchange runtime failures make zero allocations");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "exchange success and failure release all buffers");
    _CrtMemCheckpoint(&before);
    {
        allocations = 0; counting = true;
        auto temporary = std::make_unique<int>(42);
        temporary.reset();
        counting = false;
        check(allocations > 0, "temporary allocation detector positive control");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "temporary allocation freed");
    // Intentionally leak across the snapshot, then clean up the fault probe.
    _CrtMemCheckpoint(&before);
    int* leaked = std::make_unique<int>(42).release();
    _CrtMemCheckpoint(&after);
    check(outstanding(before, after), "leak detector positive control");
    delete leaked;
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "leak probe cleanup");
    _CrtMemCheckpoint(&before);
    for (int repeat = 0; repeat < 100; ++repeat) {
        allocations = 0; counting = true;
        auto result = remez::make_grid(spec);
        counting = false;
        check(result.status == Status::ok && allocations == 1, "one explicit owning setup allocation");
        allocations = 0; counting = true;
        const auto size = remez::grid_size(spec);
        const auto status = remez::fill_grid(spec, result.points);
        counting = false;
        check(size.status == Status::ok && status == Status::ok && allocations == 0, "zero sizing/fill allocations");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "repeated grids release all storage");
    _CrtMemCheckpoint(&before);
    {
        const auto unresolved = remez::make_grid(Lowpass{1, std::numeric_limits<double>::denorm_min(), pi / 2});
        check(unresolved.status == Status::numerical_resolution && unresolved.points.empty(), "numerical failure releases output");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before, after), "numerical failure has no leak");
    allocations = 0; counting = true;
    const auto failed = remez::make_grid(Lowpass{0, 0.2 * pi, 0.3 * pi});
    counting = false;
    check(failed.status == Status::invalid_taps && allocations == 0, "invalid setup does not allocate");
    for (std::size_t operation = 0; operation < 8; ++operation) {
        _CrtMemCheckpoint(&before);
        bool caught = false;
        allocations = 0; accepted_allocations = operation < 2 ? 0 : operation - 2;
        reject = true; counting = true;
        try {
            if (operation == 0) {
                const auto grid = remez::make_grid(spec);
                (void)grid;
            } else if (operation == 1) {
                remez::Interpolator interpolator(3);
            } else {
                remez::Exchange engine(spec);
            }
        } catch (const std::bad_alloc&) {
            caught = true;
        }
        counting = false; reject = false; accepted_allocations = 0;
        check(caught, "setup allocation failure propagates");
        _CrtMemCheckpoint(&after);
        check(!outstanding(before, after), "setup allocation failure has no leak");
    }
    _CrtMemCheckpoint(&before);
    {
        const std::array lowpass{remez::Band{0,.2*pi,1,1},remez::Band{.3*pi,pi,0,10}};
        const std::array hilbert{remez::Band{.05*pi,.95*pi,-1,1}};
        const std::array diff{remez::Band{0,pi,1,1}};
        const std::array specs{
            remez::Specification{32,lowpass},
            remez::Specification{31,hilbert,remez::Symmetry::antisymmetric},
            remez::Specification{32,diff,remez::Symmetry::antisymmetric,remez::Response::differentiator}};
        std::array<double,64> h{};
        for (const auto& design : specs) {
            const auto size = remez::grid_size(design);
            std::vector<remez::GridPoint> points(size.points);
            std::vector<remez::BandRange> ranges(design.bands.size());
            allocations = 0; counting = true;
            const auto sized = remez::grid_size(design);
            const auto filled = remez::fill_grid(design,points,ranges);
            const auto short_points = remez::fill_grid(design,std::span{points}.first(points.size()-1),ranges);
            const auto short_ranges = remez::fill_grid(design,points,std::span{ranges}.first(ranges.size()-1));
            counting = false;
            check(allocations == 0 && sized.status == Status::ok && filled == Status::ok && short_points == Status::insufficient_storage && short_ranges == Status::insufficient_storage,
                  "general sizing/fill and rejected short buffers allocate nothing");
            const auto payload = 32*size.points + 48*design.bands.size() + 72*(size.degree+2) + 24*size.peak_capacity + 8*size.workspace;
            allocations = 0; allocation_bytes = 0; counting = true;
            remez::Exchange engine(design);
            counting = false;
            check(allocations == 6 && allocation_bytes >= payload && allocation_bytes <= payload+384, "generic setup owns six bounded buffers");
            allocations = 0; counting = true;
            const auto result = engine.run();
            const auto status = engine.taps(h);
            const auto value = engine.amplitude(.3*pi);
            const auto rerun = engine.run();
            counting = false;
            check(allocations == 0 && result.status == Status::ok && status == Status::ok && value.status == Status::ok && rerun.status == Status::ok,
                  "first-use and repeated general design/recovery/evaluation allocate nothing");
        }
        std::array<remez::detail::Peak,10> peaks{};
        std::array<double,10> scratch{};
        for (std::size_t k=0;k<peaks.size();++k) peaks[k]={k,static_cast<double>(k),k%2==0?static_cast<double>(k+1):-static_cast<double>(k+1)};
        allocations = 0; counting = true;
        const auto count = remez::detail::select_peaks(peaks,4,scratch);
        counting = false;
        check(allocations == 0 && count == 4,"threshold sorting and general peak selection allocate nothing");
    }
    _CrtMemCheckpoint(&after);
    check(!outstanding(before,after),"all general workspace released");
    _CrtSetAllocHook(previous);
}
#endif
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view{argv[1]} == "--exchange-probe") { exchange_probe(); return 0; }
    if (argc == 2 && std::string_view{argv[1]} == "--probe") { probe(); return 0; }
    if (argc == 2 && std::string_view{argv[1]} == "--benchmark") { benchmark(); return 0; }
#ifdef REMEZ_CRT_CHECKS
    memory_contract();
#endif
    invalid_inputs();
    for (auto taps : {std::size_t{1}, std::size_t{3}, std::size_t{15}, std::size_t{31}, std::size_t{63}, std::size_t{127}, std::size_t{255}}) {
        for (auto density : {std::size_t{4}, std::size_t{8}, std::size_t{16}, std::size_t{32}}) {
            for (const auto edges : {std::pair{0.01, 0.03}, std::pair{0.2, 0.3}, std::pair{0.49, 0.51}, std::pair{0.97, 0.99}}) {
                verify({taps, edges.first * pi, edges.second * pi, {1e-100, 1e100}, density});
            }
        }
    }
    verify({31, 0.5 * pi, (0.5 + 1e-8) * pi});
    verify({15, 1e-6, 0.9 * pi, {std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::max()}});
    const auto exact = remez::grid_size({31, pi / 4, pi / 2});
    check(exact.status == Status::ok && exact.band_points == std::array<std::size_t, 2>{65, 129}, "analytic grid count");
    const auto landmarks = remez::make_grid({3, pi / 3, 2 * pi / 3, {1.0, 1.0}, 1});
    constexpr std::array known_cosines{1.0, 0.5, -0.5, -1.0};
    check(landmarks.status == Status::ok && landmarks.points.size() == known_cosines.size(), "landmark grid");
    for (std::size_t j = 0; j < known_cosines.size(); ++j)
        check(std::abs(landmarks.points[j].x - known_cosines[j]) <= tolerance, "analytic cosine landmark");
    interpolation();
    invalid_interpolation();
    range_regression();
    grid_interpolation();
    precision_regression();
    exchange();
    exchange_reference();
    tap_recovery();
    peak_selection();
    general_recovery();
    generalized();
    std::puts("PASS Phases 1, 2, 3, 4 and 5");
}



