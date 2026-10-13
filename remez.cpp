#include "remez.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <numeric>

namespace remez {
namespace {
constexpr double pi = std::numbers::pi_v<double>;

unsigned form(Symmetry symmetry, std::size_t taps) noexcept {
    return symmetry == Symmetry::symmetric ? (taps % 2 == 0 ? 2u : 1u) : (taps % 2 == 0 ? 4u : 3u);
}

// A(w) = Q(w) P(cos w); Q forces A(pi) = 0 for II and III and A(0) = 0 for III and IV.
bool forced(unsigned type, double omega) noexcept {
    return (omega == 0.0 && type >= 3) || (omega == pi && (type == 2 || type == 3));
}

double factor(unsigned type, double omega) noexcept {
    if (forced(type, omega)) return 0.0;
    switch (type) {
    case 1: return 1.0;
    case 2: return std::cos(omega / 2);
    case 3: return std::sin(omega);
    default: return std::sin(omega / 2);
    }
}

// Q(w)/w, finite at w = 0 for the antisymmetric forms that differentiators use.
double factor_over_omega(unsigned type, double omega) noexcept {
    const double t = type == 3 ? omega : omega / 2;
    return (type == 3 ? 1.0 : 0.5) * (t == 0.0 ? 1.0 : std::sin(t) / t);
}

// Grid points where A is forced to zero carry no information and are left out, except DC of a
// differentiator band, whose relative error has a finite limit.
bool omitted(unsigned type, Response response, const Band& band, double omega) noexcept {
    return forced(type, omega) && !(omega == 0.0 && response == Response::differentiator && band.target != 0.0);
}

std::size_t intervals(const Band& band, double step, std::size_t limit) noexcept {
    const double count = std::max(1.0, std::ceil((band.upper - band.lower) / step));
    return count < static_cast<double>(limit) ? static_cast<std::size_t>(count) : limit;
}

struct Layout {
    Status status;
    std::size_t points{}, degree{};
};

// Validates the specification and sizes the grid without allocating. Buffer elements are at most
// 48 bytes, so capping every count at PTRDIFF_MAX/64 keeps all byte sizes representable.
Layout layout(const Specification& spec, std::size_t max_points) noexcept {
    max_points = std::min(max_points, static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / 64);
    if ((spec.symmetry != Symmetry::symmetric && spec.symmetry != Symmetry::antisymmetric) ||
        (spec.response != Response::constant && spec.response != Response::differentiator) ||
        (spec.response == Response::differentiator && spec.symmetry != Symmetry::antisymmetric))
        return {Status::invalid_specification};
    if (spec.taps == 0 || spec.taps > max_points || (spec.symmetry == Symmetry::antisymmetric && spec.taps < 2))
        return {Status::invalid_taps};
    if (spec.bands.empty() || spec.bands.size() > max_points) return {Status::invalid_bands};
    if (spec.density == 0) return {Status::invalid_density};
    const auto type = form(spec.symmetry, spec.taps);
    const auto degree = spec.symmetry == Symmetry::symmetric ? (spec.taps - 1) / 2 : spec.taps / 2 - 1;
    if (spec.density > max_points / (degree + 1)) return {Status::grid_limit};
    const double step = pi / static_cast<double>(spec.density * (degree + 1));
    std::size_t points = 0;
    double previous = -1.0;
    for (const auto& band : spec.bands) {
        if (!(std::isfinite(band.lower) && std::isfinite(band.upper) && previous < band.lower &&
              band.lower < band.upper && band.upper <= pi)) return {Status::invalid_bands};
        if (!(std::isfinite(band.weight) && band.weight > 0.0)) return {Status::invalid_weights};
        if (!std::isfinite(band.target)) return {Status::invalid_target};
        if (band.target != 0.0 && (omitted(type, spec.response, band, band.lower) || forced(type, band.upper)))
            return {Status::incompatible_endpoint};
        const auto count = intervals(band, step, max_points) + 1 - omitted(type, spec.response, band, band.lower) -
                           omitted(type, spec.response, band, band.upper);
        if (count > max_points - points) return {Status::grid_limit};
        points += count;
        previous = band.upper;
    }
    if (points < degree + 2) return {Status::invalid_density};
    return {Status::ok, points, degree};
}

std::array<Band, 2> lowpass_bands(const Lowpass& spec) noexcept {
    return {{{0.0, spec.pass_edge, 1.0, spec.weights[0]}, {spec.stop_edge, pi, 0.0, spec.weights[1]}}};
}

// Alternation points of best approximations follow the equilibrium measure of the bands in
// x = cos(omega) (Kadec). Equispaced references act like equispaced interpolation: their levelled
// error collapses toward rounding level as the degree grows and the exchange stalls. For k bands
// with endpoints e(0) > ... > e(2k-1) the density is |q(x)|/(pi sqrt|R(x)|), R = prod (x - e(j)),
// q monic with one root per gap where the gap's mass vanishes (Akhiezer). Interval i (even: band
// i/2, odd: a gap) maps by x = c + h cos(t), which absorbs its endpoint factors: dx/sqrt(...) = dt.
// Writes the reference from grid points; roots holds the k-1 roots. False means seed otherwise.
template <class Range>
bool equilibrium(std::span<const Point> grid, std::span<const Range> bands, std::span<Point> reference,
                 std::span<double> roots) noexcept {
    const auto k = bands.size(), n = reference.size();
    if (n < k || roots.size() + 1 < k ||
        std::ranges::any_of(bands, [](const Range& b) { return b.end - b.begin < 2; })) return false;
    roots = roots.first(k - 1);
    const auto endpoint = [&](std::size_t j) { return grid[j % 2 == 0 ? bands[j / 2].begin : bands[j / 2].end - 1].x; };
    // |q(x)|/sqrt|R(x)| without interval i's own endpoints, and without root `skip` if skip < k-1.
    const auto weight = [&](double x, std::size_t i, std::size_t skip) {
        double value = 1.0;
        for (std::size_t g = 0; g < roots.size(); ++g) if (g != skip) value *= x - roots[g];
        for (std::size_t j = 0; j < 2 * k; ++j) if (j != i && j != i + 1) value /= std::sqrt(std::abs(x - endpoint(j)));
        return std::abs(value);
    };
    // Gauss-Chebyshev nodes on interval i; the measure only steers the seed, so 64 suffice.
    constexpr std::size_t nodes = 64;
    const auto node = [&](std::size_t i, std::size_t l) {
        const double t = pi * (static_cast<double>(l) + 0.5) / nodes;
        return std::midpoint(endpoint(i), endpoint(i + 1)) + (endpoint(i) - endpoint(i + 1)) / 2 * std::cos(t);
    };
    // Gauss-Seidel over the gaps: with the other roots fixed each gap condition is linear in its
    // root, the weighted mean of x over the gap. One sweep is exact for two bands.
    for (std::size_t g = 0; g < roots.size(); ++g) roots[g] = std::midpoint(endpoint(2 * g + 1), endpoint(2 * g + 2));
    for (int sweep = 0; sweep < 64; ++sweep) {
        double change = 0.0;
        for (std::size_t g = 0; g < roots.size(); ++g) {
            double moment = 0.0, mass = 0.0;
            for (std::size_t l = 0; l < nodes; ++l) {
                const double x = node(2 * g + 1, l), w = weight(x, 2 * g + 1, g);
                moment += x * w;
                mass += w;
            }
            if (!(mass > 0.0 && std::isfinite(moment / mass))) return false;
            change = std::max(change, std::abs(moment / mass - roots[g]));
            roots[g] = moment / mass;
        }
        if (change <= 0x1p-40) break;
    }
    const auto band_mass = [&](std::size_t b) {
        double mass = 0.0;
        for (std::size_t l = 0; l < nodes; ++l) mass += weight(node(2 * b, l), 2 * b, k);
        return mass;
    };
    double total = 0.0;
    for (std::size_t b = 0; b < k; ++b) total += band_mass(b);
    if (!(total > 0.0 && std::isfinite(total))) return false;
    // One point per band plus a cumulative-floor share of the rest by mass, placed at equal steps of
    // the band's own measure (trapezoid CDF in t) on distinct grid points. Band endpoints are
    // included except the last band's upper end: a mirror-symmetric reference for a mirror-symmetric
    // specification (a Hilbert band centred on pi/2) makes the levelled error of an even target vanish.
    const auto spare = n - k;
    std::size_t out = 0, used = 0;
    double before = 0.0;
    for (std::size_t b = 0; b < k; ++b) {
        before += band_mass(b);
        const auto share = b + 1 == k ? spare : std::min(spare, static_cast<std::size_t>(static_cast<double>(spare) * before / total));
        const auto m = 1 + share - used;
        used = share;
        const auto band = grid.subspan(bands[b].begin, bands[b].end - bands[b].begin);
        if (m > band.size()) return false;
        const double c = std::midpoint(band.front().x, band.back().x), h = (band.front().x - band.back().x) / 2;
        const auto density = [&](std::size_t j) { return weight(band[j].x, 2 * b, k); };
        const auto angle = [&](std::size_t j) { return std::acos(std::clamp((band[j].x - c) / h, -1.0, 1.0)); };
        const auto step = [&](std::size_t j) { return (density(j - 1) + density(j)) / 2 * (angle(j) - angle(j - 1)); };
        double area = 0.0;
        for (std::size_t j = 1; j < band.size(); ++j) area += step(j);
        if (!(area > 0.0 && std::isfinite(area))) return false;
        const auto steps = b + 1 == k ? m : m - 1;
        const auto target = [&](std::size_t j) { return steps == 0 ? area / 2 : area * static_cast<double>(j) / static_cast<double>(steps); };
        std::size_t placed = 0, previous = 0;
        const auto place = [&](std::size_t j) {
            j = std::min(std::max(j, placed == 0 ? 0 : previous + 1), band.size() - (m - placed));
            reference[out++] = band[j];
            previous = j;
            ++placed;
        };
        double cdf = 0.0, last = 0.0;
        for (std::size_t j = 0; j < band.size() && placed < m; ++j) {
            if (j != 0) cdf += step(j);
            while (placed < m && cdf >= target(placed))
                place(j != 0 && target(placed) - last < cdf - target(placed) ? j - 1 : j);
            last = cdf;
        }
        while (placed < m) place(band.size() - 1);
    }
    return true;
}
}

Exchange::Exchange(const Lowpass& spec, std::size_t max_points)
    : Exchange(Specification{spec.taps % 2 == 0 ? 0 : spec.taps, lowpass_bands(spec), Symmetry::symmetric,
                             Response::constant, spec.density}, max_points) {}

Exchange::Exchange(const Specification& spec, std::size_t max_points)
    : symmetry_(spec.symmetry), response_(spec.response) {
    const auto size = layout(spec, max_points);
    status_ = size.status;
    if (status_ != Status::ok) return;
    type_ = form(spec.symmetry, spec.taps);
    grid_.resize(size.points);
    reference_.resize(size.degree + 2);
    bands_.resize(spec.bands.size());
    nodes_.resize(size.degree + 2);
    peaks_.resize(size.points);  // Every sign run holds at least one grid point.
    samples_.resize(spec.taps / 2 + 1);
    taps_.resize(spec.taps);
    const double step = pi / static_cast<double>(spec.density * (size.degree + 1));
    std::size_t count = 0;
    for (std::size_t b = 0; b < spec.bands.size(); ++b) {
        const auto& band = spec.bands[b];
        const auto L = intervals(band, step, max_points);
        bands_[b] = {band, count, count};
        for (std::size_t j = 0; j <= L; ++j) {
            const double omega = std::lerp(band.lower, band.upper, static_cast<double>(j) / static_cast<double>(L));
            if (omitted(type_, response_, band, omega)) continue;
            const auto p = point(band, omega);
            if (!(std::isfinite(p.desired) && std::isfinite(p.weight) && p.weight > 0.0) ||
                (count != 0 && !(p.x < grid_[count - 1].x))) { status_ = Status::numerical_failure; return; }
            grid_[count++] = p;
        }
        bands_[b].end = count;
    }
}

// Transformed target and weight: W(D - QP) = WQ(D/Q - P). For a differentiator band with slope s,
// the relative error W(s w - QP)/(|s| w) = (W r/|s|)(s/r - P) with r = Q/w.
Point Exchange::point(const Band& band, double omega) const noexcept {
    const double x = std::cos(omega);
    if (response_ == Response::differentiator && band.target != 0.0) {
        const double r = factor_over_omega(type_, omega);
        return {omega, x, band.target / r, band.weight * r / std::abs(band.target)};
    }
    const double q = factor(type_, omega);
    return {omega, x, band.target / q, band.weight * q};
}

const Exchange::Range& Exchange::range(std::size_t index) const noexcept {
    return *std::ranges::lower_bound(bands_, index, {}, [](const Range& r) { return r.end - 1; });
}

double Exchange::interpolate(double x, double Node::*value) const noexcept {
    double numerator = 0.0, denominator = 0.0;
    for (const auto& node : std::span{nodes_}.first(nodes_.size() - 1)) {
        if (x == node.x) return node.*value;
        const double t = node.w / (x - node.x);
        numerator += t * (node.*value);
        denominator += t;
    }
    return numerator / denominator;
}

// Barycentric weights w_k = 1/prod_{j != k} 2(x_k - x_j); doubling maps [-1, 1] to capacity 1, so
// products stay moderate for well-spread nodes, and frexp every 16 factors keeps any set in range.
// The levelled error follows from sum w_k P(x_k) = 0 for deg P < n-1 with W_k(D_k - P(x_k)) =
// (-1)^k delta, and the node values are P(x_k) = D_k - (-1)^k delta/W_k.
Status Exchange::fit(double& delta) noexcept {
    const auto n = reference_.size();
    const auto difference = [&](std::size_t k, std::size_t j) { return 2 * (reference_[k].x - reference_[j].x); };
    int largest = std::numeric_limits<int>::min();
    for (std::size_t k = 0; k < n; ++k) {
        double product = 1.0;
        int exponent = 0, e = 0;
        for (std::size_t j = 0, factors = 0; j < n; ++j) {
            if (j == k) continue;
            product *= difference(k, j);
            if (++factors % 16 == 0) { product = std::frexp(product, &e); exponent += e; }
        }
        product = std::frexp(product, &e);
        exponent += e;
        if (product == 0.0 || !std::isfinite(product)) return Status::numerical_failure;
        // y holds the weight's binary exponent until the node values are known.
        nodes_[k] = {reference_[k].omega, reference_[k].x, 1.0 / product, static_cast<double>(-exponent), 0.0};
        largest = std::max(largest, -exponent);
    }
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        auto& node = nodes_[k];
        node.w = std::ldexp(node.w, static_cast<int>(node.y) - largest);
        if (node.w == 0.0) return Status::numerical_failure;
        numerator += node.w * reference_[k].desired;
        denominator += (k % 2 == 0 ? node.w : -node.w) / reference_[k].weight;
    }
    delta = numerator / denominator;
    if (!std::isfinite(delta)) return Status::numerical_failure;
    for (std::size_t k = 0; k < n; ++k)
        nodes_[k].y = reference_[k].desired - (k % 2 == 0 ? delta : -delta) / reference_[k].weight;
    // All n nodes define a degree n-1 interpolant whose top coefficient sum w_k y_k vanishes only up
    // to rounding. Evaluating on n-1 nodes is exactly degree n-2; the weights of the reduced set are
    // w_k (x_k - x_m). Dropping the largest weight keeps the remaining terms well scaled.
    const auto m = static_cast<std::size_t>(std::ranges::max_element(nodes_, {}, [](const Node& node) { return std::abs(node.w); }) - nodes_.begin());
    for (std::size_t k = 0; k < n; ++k) if (k != m) nodes_[k].w *= difference(k, m);
    std::swap(nodes_[m], nodes_.back());
    return Status::ok;
}

// Successive parabolic interpolation on the grid bracket around a sampled peak. Band-edge peaks
// stay put: the edge is a grid point and the bracket must not cross into the transition.
void Exchange::refine(Peak& peak) const noexcept {
    const auto j = peak.index;
    const auto& r = range(j);
    if (j == r.begin || j + 1 == r.end) return;
    const double sign = peak.error > 0 ? 1.0 : -1.0;
    double a = grid_[j - 1].omega, b = peak.omega, c = grid_[j + 1].omega;
    double fa = sign * error(grid_[j - 1]), fb = sign * peak.error, fc = sign * error(grid_[j + 1]);
    const double resolution = 1e-6 * (c - a);
    for (int iteration = 0; iteration < 16; ++iteration) {
        const double p = (b - a) * (fb - fc), q = (b - c) * (fb - fa);
        if (p == q) break;
        const double u = b - ((b - a) * p - (b - c) * q) / (2 * (p - q));
        if (!(u > a && u < c) || std::abs(u - b) <= resolution) break;
        const double fu = sign * error(point(r.band, u));
        if (fu > fb) {
            (u < b ? c : a) = b;
            (u < b ? fc : fa) = fb;
            b = u;
            fb = fu;
        } else {
            (u < b ? a : c) = u;
            (u < b ? fa : fc) = fu;
        }
    }
    peak.omega = b;
    peak.error = sign * fb;
}

Result Exchange::run(std::size_t max_iterations, double tolerance) noexcept {
    ready_ = exact_ = false;
    if (status_ != Status::ok) return {status_};
    if (grid_.empty()) return {Status::invalid_specification};  // A moved-from engine.
    if (!(tolerance > 0.0 && tolerance < 1.0)) return {Status::invalid_tolerance};
    if (max_iterations == 0) return {Status::nonconvergence};
    // Zero targets, or one shared Type-I constant, are met exactly by a constant amplitude.
    constant_ = bands_.front().band.target;
    const bool zero = std::ranges::all_of(bands_, [](const Range& r) { return r.band.target == 0.0; });
    if (zero || (type_ == 1 && response_ == Response::constant &&
                 std::ranges::all_of(bands_, [&](const Range& r) { return r.band.target == constant_; }))) {
        std::ranges::fill(taps_, 0.0);
        if (type_ == 1) taps_[taps_.size() / 2] = constant_;
        ready_ = exact_ = true;
        return {Status::ok};
    }
    const auto n = reference_.size();
    // With more bands than points, uniform indices may all share one target, which would make the
    // levelled error vanish identically: then one point moves to the most different target.
    if (!equilibrium<Range>(grid_, bands_, reference_, samples_)) {
        const auto step = (grid_.size() - 1) / (n - 1), extra = (grid_.size() - 1) % (n - 1);
        for (std::size_t k = 0, index = 0, carry = 0; k < n; ++k) {
            reference_[k] = grid_[index];
            index += step;
            carry += extra;
            if (carry >= n - 1) { ++index; carry -= n - 1; }
        }
        const auto [low, high] = std::ranges::minmax(reference_, {}, &Point::desired);
        if (low.desired == high.desired) {
            const auto [lowest, highest] = std::ranges::minmax(grid_, {}, &Point::desired);
            reference_[n / 2] = std::abs(lowest.desired - low.desired) > std::abs(highest.desired - low.desired) ? lowest : highest;
            std::ranges::sort(reference_, {}, &Point::omega);
        }
    }
    const auto magnitude = [](const Peak& p) { return std::abs(p.error); };
    Result result{Status::nonconvergence};
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        result.iterations = iteration + 1;
        double delta = 0.0;
        if (const auto status = fit(delta); status != Status::ok) { result.status = status; return result; }
        result.ripple = std::abs(delta);
        // Each sign run of the sampled error contributes its largest point.
        std::size_t count = 0;
        for (std::size_t j = 0; j < grid_.size(); ++j) {
            const double e = error(grid_[j]);
            if (!std::isfinite(e)) { result.status = Status::numerical_failure; return result; }
            if (e == 0.0) continue;
            if (count != 0 && std::signbit(e) == std::signbit(peaks_[count - 1].error)) {
                if (std::abs(e) > std::abs(peaks_[count - 1].error)) peaks_[count - 1] = {j, grid_[j].omega, e};
            } else {
                peaks_[count++] = {j, grid_[j].omega, e};
            }
        }
        if (count < n) { result.status = Status::numerical_failure; return result; }
        auto peaks = std::span{peaks_}.first(count);
        for (auto& peak : peaks) refine(peak);
        result.max_error = magnitude(std::ranges::max(peaks, {}, magnitude));
        // Keep n alternating peaks: the weakest goes, an interior one with its weaker equal-sign
        // neighbour; with n+1 left, the weaker end goes.
        while (peaks.size() > n) {
            std::size_t k = peaks.size() == n + 1 ? (magnitude(peaks.front()) < magnitude(peaks.back()) ? 0 : n)
                : static_cast<std::size_t>(std::ranges::min_element(peaks, {}, magnitude) - peaks.begin());
            const std::size_t removed = k == 0 || k + 1 == peaks.size() ? 1 : 2;
            if (removed == 2 && magnitude(peaks[k + 1]) > magnitude(peaks[k - 1])) peaks[k - 1] = peaks[k + 1];
            std::shift_left(peaks.begin() + static_cast<std::ptrdiff_t>(k), peaks.end(), static_cast<std::ptrdiff_t>(removed));
            peaks = peaks.first(peaks.size() - removed);
        }
        // Refined peaks that meet mean the error is unresolved in binary64.
        if (std::ranges::adjacent_find(peaks, std::ranges::greater_equal{}, &Peak::omega) != peaks.end()) {
            result.status = Status::numerical_failure;
            return result;
        }
        result.flatness = (result.max_error - magnitude(std::ranges::min(peaks, {}, magnitude))) / result.max_error;
        if (result.flatness <= tolerance) {
            result.status = recover();
            ready_ = result.status == Status::ok;
            return result;
        }
        for (std::size_t k = 0; k < n; ++k) reference_[k] = point(range(peaks[k].index).band, peaks[k].omega);
    }
    return result;
}

// Direct real inverse DFT of A = Q P, with P the interpolant of the given node values, sampled at
// w_m = 2 pi m/N. With t = n - (N-1)/2 and c_m = 1 at m = 0 and m = N/2, else 2:
// h[n] = (1/N) sum c_m A_m cos(w_m t), or -(1/N) sum c_m A_m sin(w_m t) when antisymmetric.
// w_m t = pi q/N with integer q = m(2n+1-N), reduced mod 2N to stay exact.
void Exchange::inverse(double Node::*value, bool accumulate) noexcept {
    const auto N = taps_.size();
    for (std::size_t m = 0; m < samples_.size(); ++m) {
        const double omega = 2 * pi * static_cast<double>(m) / static_cast<double>(N);
        samples_[m] = forced(type_, omega) ? 0.0 : factor(type_, omega) * interpolate(std::cos(omega), value);
    }
    const bool symmetric = symmetry_ == Symmetry::symmetric;
    const auto period = static_cast<std::ptrdiff_t>(2 * N);
    for (std::size_t n = 0; n < (N + 1) / 2; ++n) {
        const auto t = static_cast<std::ptrdiff_t>(2 * n + 1) - static_cast<std::ptrdiff_t>(N);
        double sum = 0.0;
        for (std::size_t m = 0; m < samples_.size(); ++m) {
            const auto q = (static_cast<std::ptrdiff_t>(m) * t % period + period) % period;
            const double angle = pi * static_cast<double>(q) / static_cast<double>(N);
            sum += (m == 0 || 2 * m == N ? 1.0 : 2.0) * samples_[m] * (symmetric ? std::cos(angle) : -std::sin(angle));
        }
        // The centre of an odd antisymmetric filter is exactly zero (t = 0, sin = 0).
        const double h = t == 0 && !symmetric ? 0.0 : sum / static_cast<double>(N);
        taps_[n] = accumulate ? taps_[n] + h : h;
        taps_[N - 1 - n] = symmetric ? taps_[n] : -taps_[n];
    }
}

// P(cos w) = A(w)/Q(w) of the current taps, with A = sum h[k] cos(w t_k) or -sum h[k] sin(w t_k),
// t_k = k - (N-1)/2. At DC of an antisymmetric form both vanish; the ratio is A'(0)/Q'(0).
double Exchange::tap_polynomial(double omega) const noexcept {
    const auto N = taps_.size();
    const bool symmetric = symmetry_ == Symmetry::symmetric;
    double a = 0.0;
    for (std::size_t k = 0; k < N; ++k) {
        const double t = (static_cast<double>(2 * k + 1) - static_cast<double>(N)) / 2;
        a += omega == 0.0 && !symmetric ? -taps_[k] * t : taps_[k] * (symmetric ? std::cos(omega * t) : -std::sin(omega * t));
    }
    return omega == 0.0 && !symmetric ? a / factor_over_omega(type_, 0.0) : a / factor(type_, omega);
}

// Samples inside transition bands come from an extrapolating interpolant: their rounding is
// amplified by its Lebesgue function there and leaks into the bands through the transform. One step
// of iterative refinement removes it. The taps' polynomial and P both have degree d, so their
// difference is reproduced exactly by interpolating the residuals at the d+1 active nodes, where
// the values are exact data; the correction's own rounding scales with those tiny residuals.
Status Exchange::recover() noexcept {
    inverse(&Node::y, false);
    for (auto& node : std::span{nodes_}.first(nodes_.size() - 1)) node.residual = node.y - tap_polynomial(node.omega);
    inverse(&Node::residual, true);
    return std::ranges::all_of(taps_, [](double h) { return std::isfinite(h); }) ? Status::ok : Status::numerical_failure;
}

std::span<const double> Exchange::taps() const noexcept {
    return ready_ ? std::span<const double>{taps_} : std::span<const double>{};
}

double Exchange::amplitude(double omega) const noexcept {
    if (!ready_ || !(omega >= 0.0 && omega <= pi)) return std::numeric_limits<double>::quiet_NaN();
    if (exact_) return type_ == 1 ? constant_ : 0.0;
    return forced(type_, omega) ? 0.0 : factor(type_, omega) * interpolate(std::cos(omega));
}

std::span<const Point> Exchange::extrema() const noexcept {
    return ready_ && !exact_ ? std::span<const Point>{reference_} : std::span<const Point>{};
}

}
