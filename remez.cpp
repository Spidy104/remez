#include "remez.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>

namespace remez {
namespace {
constexpr double pi = std::numbers::pi_v<double>;
constexpr auto size_max = std::numeric_limits<std::size_t>::max();

std::array<Band, 2> bands(const Lowpass& spec) noexcept {
    return {{{0.0, spec.pass_edge, 1.0, spec.weights[0]}, {spec.stop_edge, pi, 0.0, spec.weights[1]}}};
}

std::size_t intervals(const Band& band, double step) noexcept {
    const double count = std::max(1.0, std::ceil((band.upper - band.lower) / step));
    return std::isfinite(count) && count < static_cast<double>(size_max) ? static_cast<std::size_t>(count) : 0;
}
}

GridSize grid_size(const Lowpass& spec, std::size_t max_points) noexcept {
    if (spec.taps == 0 || spec.taps%2 == 0) return {Status::invalid_taps};
    const auto active = bands(spec);
    const auto size = grid_size(Specification{spec.taps, active, Symmetry::symmetric, Response::constant, spec.density}, max_points);
    if (size.status != Status::ok) return {size.status};
    const double step = pi/static_cast<double>(spec.density*(size.degree+1));
    GridSize result{Status::ok};
    std::ranges::transform(active, result.band_points.begin(), [step](const Band& band) { return intervals(band, step) + 1; });
    return result;
}

Status fill_grid(const Lowpass& spec, std::span<GridPoint> output) noexcept {
    const auto size = grid_size(spec, size_max);
    if (size.status != Status::ok) return size.status;
    const auto active = bands(spec);
    std::array<BandRange,2> ranges{};
    return fill_grid(Specification{spec.taps, active, Symmetry::symmetric, Response::constant, spec.density}, output, ranges);
}

Grid make_grid(const Lowpass& spec, std::size_t max_points) {
    const auto size = grid_size(spec, max_points);
    if (size.status != Status::ok)
        return {size.status, {}, 0};
    Grid grid{Status::ok, std::vector<GridPoint>(size.count()), size.band_points[0]};
    grid.status = fill_grid(spec, grid.points);
    if (grid.status != Status::ok)
        return {grid.status, {}, 0};
    return grid;
}

namespace {
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53 &&
              std::numeric_limits<double>::max_exponent == 1024);

struct DoubleDouble {
    double hi{}, lo{};
};

[[nodiscard]] DoubleDouble sum(double a, double b) noexcept {
    const double hi = a + b, back = hi - a;
    return {hi, (a - (hi - back)) + (b - back)};
}

[[nodiscard]] DoubleDouble operator+(DoubleDouble a, DoubleDouble b) noexcept {
    const auto high = sum(a.hi, b.hi);
    return sum(high.hi, high.lo + (a.lo + b.lo));
}

[[nodiscard]] DoubleDouble operator-(DoubleDouble a, DoubleDouble b) noexcept {
    return a + DoubleDouble{-b.hi, -b.lo};
}

[[nodiscard]] DoubleDouble operator*(DoubleDouble a, DoubleDouble b) noexcept {
    const double hi = a.hi * b.hi;
    return sum(hi, std::fma(a.hi, b.hi, -hi) + (a.hi * b.lo + a.lo * b.hi));
}

[[nodiscard]] DoubleDouble operator/(DoubleDouble a, DoubleDouble b) noexcept {
    const double first = a.hi / b.hi;
    const auto remainder = a - b * DoubleDouble{first};
    return sum(first, (remainder.hi + remainder.lo) / b.hi);
}

[[nodiscard]] DoubleDouble scale(DoubleDouble a, int exponent) noexcept {
    return {std::scalbn(a.hi, exponent), std::scalbn(a.lo, exponent)};
}

[[nodiscard]] DoubleDouble product(DoubleDouble a, double b, double c, bool* lost = nullptr) noexcept {
    int ea = 0, eb = 0, ec = 0;
    const double ma = std::frexp(a.hi, &ea), mb = std::frexp(b, &eb), mc = std::frexp(c, &ec);
    const auto normalized = DoubleDouble{ma, std::scalbn(a.lo, -ea)} * DoubleDouble{mb} * DoubleDouble{mc};
    const int exponent = ea + eb + ec;
    const auto value = scale(normalized, exponent);
    if (lost != nullptr) {
        const auto restored = scale(value, -exponent);
        *lost = restored.hi != normalized.hi || restored.lo != normalized.lo;
    }
    return value;
}

void normalize(DoubleDouble& factor, long long& exponent) noexcept {
    int e = 0;
    factor.hi = std::frexp(factor.hi, &e);
    factor.lo = std::scalbn(factor.lo, -e);
    exponent += e;
}

[[nodiscard]] bool finite(DoubleDouble a) noexcept {
    return std::isfinite(a.hi) && std::isfinite(a.lo);
}

[[nodiscard]] Value rounded(DoubleDouble a) noexcept {
    const double value = a.hi + a.lo;
    return finite(a) && std::isfinite(value) ? Value{Status::ok, value} : Value{Status::numerical_failure};
}

[[nodiscard]] double amplitude_bound(double value, double bound = 0.0) noexcept {
    // Relative bounds alone vanish under gradual underflow; round the absolute allowance upward.
    return std::nextafter(bound + std::numeric_limits<double>::epsilon()*std::abs(value) +
                          std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::infinity());
}

// Compensated sum with the same error as Neumaier's, but Knuth's branch-free two-sum: the
// magnitude test mispredicts on alternating barycentric terms and dominated evaluation time.
struct Sum {
    double total{}, correction{};
    void add(double term) noexcept {
        const auto exact = sum(total, term);
        total = exact.hi;
        correction += exact.lo;
    }
    [[nodiscard]] double value() const noexcept { return total + correction; }
};

struct Estimate { double error, bound; };

// E = W(D-A) from a bounded binary64 amplitude, with |E - exact| <= bound when amplitude is ok.
[[nodiscard]] Estimate bounded_error(const GridPoint& point, const Bounded& amplitude) noexcept {
    const double error = point.weight * (point.desired - amplitude.value);
    return {error, point.weight * amplitude.bound + 2 * std::numeric_limits<double>::epsilon() * std::abs(error) +
                   std::numeric_limits<double>::denorm_min()};
}

// The estimate, if its bound keeps the sign exact and is within max(relative*|E|, absolute);
// otherwise the caller needs the two-component path.
[[nodiscard]] std::optional<double> certified_error(const GridPoint& point, const Bounded& amplitude,
                                      double relative, double absolute) noexcept {
    const auto [error, bound] = bounded_error(point, amplitude);
    if (amplitude.status == Status::ok && std::isfinite(error) && bound < std::abs(error) &&
        bound <= std::max(relative * std::abs(error), absolute))
        return error;
    return std::nullopt;
}
}

namespace {
unsigned form(Symmetry symmetry, std::size_t taps) noexcept {
    return symmetry == Symmetry::symmetric ? (taps % 2 == 0 ? 2u : 1u) : (taps % 2 == 0 ? 4u : 3u);
}
unsigned form(const Specification& spec) noexcept { return form(spec.symmetry, spec.taps); }

bool forced(unsigned type, double omega) noexcept {
    return (omega == 0.0 && type >= 3) || (omega == pi && (type == 2 || type == 3));
}

double factor(unsigned type, double omega) noexcept {
    if (type == 1) return 1.0;
    if (forced(type, omega)) return 0.0;
    if (type == 2) return std::cos(omega / 2);
    return type == 3 ? std::sin(omega) : std::sin(omega / 2);
}

double ratio(unsigned type, double omega) noexcept {
    const double x = type == 3 ? omega : omega / 2;
    return (type == 3 ? 1.0 : 0.5) * (x == 0.0 ? 1.0 : std::sin(x) / x);
}

// a*b/c without overflowing a reciprocal or an intermediate product.
Value scaled_ratio(double a, double b, double c) noexcept {
    if (a == 0.0 || b == 0.0) return {Status::ok};
    int ea{}, eb{}, ec{};
    const auto x = DoubleDouble{std::frexp(a, &ea)} * DoubleDouble{std::frexp(b, &eb)} /
                   DoubleDouble{std::frexp(c, &ec)};
    return rounded(scale(x, ea + eb - ec));
}

struct Transformed { Status status; GridPoint point; };
Transformed transform(unsigned type, Response response, const Band& band, double omega) noexcept {
    const double q = factor(type, omega);
    Value desired{}, weight{};
    if (response == Response::differentiator && band.target != 0.0) {
        const double a = ratio(type, omega);
        desired = scaled_ratio(band.target, 1.0, a);
        weight = scaled_ratio(band.weight, a, std::abs(band.target));
    } else {
        if (q == 0.0) return {Status::numerical_resolution, {}};
        desired = scaled_ratio(band.target, 1.0, q);
        weight = scaled_ratio(band.weight, q, 1.0);
    }
    if (desired.status != Status::ok || weight.status != Status::ok || weight.value <= 0.0)
        return {Status::numerical_failure, {}};
    return {Status::ok, {omega, std::cos(omega), desired.value, weight.value}};
}

bool omit(unsigned type, Response response, const Band& band, double omega) noexcept {
    return forced(type, omega) && !(omega == 0.0 && response == Response::differentiator && band.target != 0.0);
}
}

DesignSize grid_size(const Specification& spec, std::size_t max_points) noexcept {
    if ((spec.symmetry != Symmetry::symmetric && spec.symmetry != Symmetry::antisymmetric) ||
        (spec.response != Response::constant && spec.response != Response::differentiator) ||
        (spec.response == Response::differentiator && spec.symmetry != Symmetry::antisymmetric))
        return {Status::invalid_specification};
    if (spec.taps == 0 || (spec.symmetry == Symmetry::antisymmetric && spec.taps < 2)) return {Status::invalid_taps};
    if (spec.bands.empty()) return {Status::invalid_bands};
    const auto type = form(spec);
    const auto d = spec.symmetry == Symmetry::symmetric ? (spec.taps - 1) / 2 : spec.taps / 2 - 1;
    if (spec.density == 0) return {Status::invalid_density};
    constexpr auto limit = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max());
    constexpr auto node_bytes = sizeof(double) + 2*sizeof(DoubleDouble);
    if (d+2 > limit/node_bytes || spec.bands.size() > limit/sizeof(BandRange)) return {Status::grid_limit};
    if (spec.density > size_max / (d + 1)) return {Status::grid_limit};
    const double step = pi / static_cast<double>(spec.density * (d + 1));
    std::size_t M = 0;
    double previous = -1.0;
    for (const auto& band : spec.bands) {
        if (!(std::isfinite(band.lower) && std::isfinite(band.upper) && band.lower >= 0.0 &&
              band.lower < band.upper && band.upper <= pi && previous < band.lower)) return {Status::invalid_edges};
        if (!std::isfinite(band.target)) return {Status::invalid_target};
        if (!(std::isfinite(band.weight) && band.weight > 0.0)) return {Status::invalid_weights};
        if (band.target != 0.0 && ((forced(type, band.lower) &&
             !(band.lower == 0.0 && spec.response == Response::differentiator)) || forced(type, band.upper)))
            return {Status::incompatible_endpoint};
        const auto L = intervals(band, step);
        if (L == 0 || L == size_max) return {Status::grid_limit};
        const auto skip = static_cast<std::size_t>(omit(type, spec.response, band, band.lower)) +
                          static_cast<std::size_t>(omit(type, spec.response, band, band.upper));
        const auto count = L + 1 - skip;
        if (count == 0) return {Status::insufficient_points};
        if (count > max_points - M) return {Status::grid_limit};
        M += count;
        previous = band.upper;
    }
    if (M < d + 2) return {Status::insufficient_points};
    if (M > limit / sizeof(GridPoint)) return {Status::grid_limit};
    const auto stationary = type == 1 ? (d == 0 ? 0 : d - 1) : (type == 3 ? d + 1 : d);
    const auto K = spec.response == Response::differentiator ? M + d + 2 : std::max(d + 2, stationary + 2 * spec.bands.size());
    const auto W = 2 * (d + 1) + std::max({d + 1, K, spec.taps - (d + 1)});
    if (K > limit / sizeof(detail::Peak) || W > limit / sizeof(double)) return {Status::grid_limit};
    return {Status::ok, M, d, K, W};
}

Status fill_grid(const Specification& spec, std::span<GridPoint> output, std::span<BandRange> ranges) noexcept {
    const auto size = grid_size(spec, size_max);
    if (size.status != Status::ok) return size.status;
    if (output.size() < size.points || ranges.size() < spec.bands.size()) return Status::insufficient_storage;
    const auto type = form(spec);
    const double step = pi / static_cast<double>(spec.density * (size.degree + 1));
    std::size_t offset = 0;
    double previous = 2.0;
    for (std::size_t b = 0; b < spec.bands.size(); ++b) {
        const auto& band = spec.bands[b];
        const auto L = intervals(band, step), begin = offset;
        for (std::size_t j = 0; j <= L; ++j) {
            const double omega = std::lerp(band.lower, band.upper, static_cast<double>(j) / static_cast<double>(L));
            if (omit(type, spec.response, band, omega)) continue;
            const auto value = transform(type, spec.response, band, omega);
            if (value.status != Status::ok) return value.status;
            if (!(value.point.x < previous)) return Status::numerical_resolution;
            output[offset++] = value.point;
            previous = value.point.x;
        }
        ranges[b] = {band, begin, offset};
    }
    return Status::ok;
}

DesignGrid make_grid(const Specification& spec, std::size_t max_points) {
    const auto size = grid_size(spec, max_points);
    if (size.status != Status::ok) return {size.status, {}, {}};
    DesignGrid grid{Status::ok, std::vector<GridPoint>(size.points), std::vector<BandRange>(spec.bands.size()),
                    size.degree, size.peak_capacity, size.workspace};
    grid.status = fill_grid(spec, grid.points, grid.bands);
    return grid.status == Status::ok ? std::move(grid) : DesignGrid{grid.status, {}, {}};
}

std::size_t detail::select_peaks(std::span<Peak> peaks, std::size_t n, std::span<double> scratch) noexcept {
    if (n == 0 || peaks.size() < n || scratch.size() < peaks.size()) return 0;
    if (peaks.size() == n) return n;
    if (peaks.size() == n + 1) {
        if (std::abs(peaks.front().error) < std::abs(peaks.back().error)) std::shift_left(peaks.begin(), peaks.end(), 1);
        return n;
    }
    const auto levels = scratch.first(peaks.size());
    std::ranges::transform(peaks, levels.begin(), [](const Peak& peak) { return std::abs(peak.error); });
    std::ranges::sort(levels);
    const auto feasible = [&](double threshold) noexcept {
        std::size_t count = 0;
        bool sign = false;
        for (const auto& peak : peaks) {
            if (std::abs(peak.error) < threshold) continue;
            if (count == 0 || sign != std::signbit(peak.error)) { ++count; sign = std::signbit(peak.error); }
        }
        return count >= n;
    };
    const auto boundary = std::ranges::partition_point(levels, feasible);
    if (boundary == levels.begin()) return 0;
    const double threshold = *(boundary - 1);
    std::size_t count = 0;
    for (const auto peak : peaks) {
        if (std::abs(peak.error) < threshold) continue;
        if (count != 0 && std::signbit(peak.error) == std::signbit(peaks[count - 1].error)) {
            if (std::abs(peak.error) > std::abs(peaks[count - 1].error)) peaks[count - 1] = peak;
        } else peaks[count++] = peak;
    }
    const auto largest = static_cast<std::size_t>(std::ranges::max_element(peaks.first(count), {},
        [](const detail::Peak& peak) { return std::abs(peak.error); }) - peaks.begin());
    const auto first = static_cast<std::ptrdiff_t>(std::min(largest > n / 2 ? largest - n / 2 : 0, count - n));
    std::shift_left(peaks.begin(), peaks.begin() + first + static_cast<std::ptrdiff_t>(n), first);
    return n;
}

struct Interpolator::Node {
    double x{};
    DoubleDouble gamma{}, value{};
};

struct Interpolator::PreciseValue {
    Status status;
    DoubleDouble value{};
    double bound{};  // Value and bound are normalized by value_scale_.
};

Interpolator::Interpolator(std::size_t node_count) : nodes_(node_count) {}
Interpolator::~Interpolator() = default;
Interpolator::Interpolator(Interpolator&&) noexcept = default;
Interpolator& Interpolator::operator=(Interpolator&&) noexcept = default;

Value Interpolator::fit(std::span<const GridPoint> extrema) noexcept {
    ready_ = false;
    const auto n = nodes_.size();
    if (n < 2 || extrema.size() != n)
        return {Status::invalid_nodes};
    // Scratch exponents must be exactly representable in binary64.
    if (n > (1ULL << std::numeric_limits<double>::digits) / 1078)
        return {Status::numerical_failure};
    double min_weight = std::numeric_limits<double>::max(), desired_scale = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const auto& point = extrema[k];
        if (!(std::isfinite(point.x) && point.x >= -1.0 && point.x <= 1.0 && std::isfinite(point.desired)))
            return {Status::invalid_nodes};
        if (!(std::isfinite(point.weight) && point.weight > 0.0))
            return {Status::invalid_weights};
        if (k != 0 && point.x == extrema[k - 1].x)
            return {Status::numerical_resolution};
        if (k != 0 && point.x > extrema[k - 1].x)
            return {Status::invalid_nodes};
        min_weight = std::min(min_weight, point.weight);
        desired_scale = std::max(desired_scale, std::abs(point.desired));
        nodes_[k].x = point.x;
    }
    if (desired_scale == 0.0) desired_scale = 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        DoubleDouble factor{1.0};
        long long exponent = 0;
        int pending = 0;
        for (std::size_t j = 0; j < n; ++j) {
            if (j == k) continue;
            auto difference = sum(nodes_[k].x, -nodes_[j].x);
            if (std::abs(difference.hi) < 0x1p-100) {
                int e = 0;
                difference.hi = std::frexp(difference.hi, &e);
                difference.lo = std::scalbn(difference.lo, -e);
                exponent += e;
            }
            factor = factor * difference;
            // Eight bounded factors stay in [2^-801, 2^8], leaving room for both components.
            if (++pending == 8) {
                normalize(factor, exponent);
                pending = 0;
            }
        }
        if (pending != 0) normalize(factor, exponent);
        nodes_[k].gamma = DoubleDouble{1.0} / factor;
        // Values are unused until the ripple is solved; reuse their storage.
        nodes_[k].value.hi = static_cast<double>(-exponent);
    }
    const auto max_exponent = static_cast<long long>(std::ranges::max_element(nodes_, {},
        [](const Node& node) { return node.value.hi; })->value.hi);
    DoubleDouble numerator{}, denominator{};
    for (std::size_t k = 0; k < n; ++k) {
        const auto shift = static_cast<long long>(nodes_[k].value.hi) - max_exponent;
        if (shift < -1074) return {Status::numerical_failure};
        nodes_[k].gamma = scale(nodes_[k].gamma, static_cast<int>(shift));
        const auto q = DoubleDouble{min_weight} / DoubleDouble{extrema[k].weight};
        const auto desired = DoubleDouble{extrema[k].desired} / DoubleDouble{desired_scale};
        if (!finite(nodes_[k].gamma) || nodes_[k].gamma.hi == 0.0 || !finite(q) || q.hi <= 0.0 ||
            (extrema[k].desired != 0.0 && desired.hi == 0.0))
            return {Status::numerical_failure};
        numerator = numerator + nodes_[k].gamma * desired;
        const auto magnitude = k % 2 == 0 ? nodes_[k].gamma : DoubleDouble{-nodes_[k].gamma.hi, -nodes_[k].gamma.lo};
        denominator = denominator + magnitude * q;
    }
    const auto delta_scaled = numerator / denominator;
    const auto delta = product(delta_scaled, desired_scale, min_weight);
    if (!finite(delta_scaled) || !finite(delta) || !std::isfinite(delta.hi + delta.lo) || (delta_scaled.hi != 0.0 && delta.hi == 0.0))
        return {Status::numerical_failure};
    value_scale_ = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const auto correction = product(delta_scaled * (DoubleDouble{min_weight} / DoubleDouble{extrema[k].weight}), desired_scale, 1.0);
        nodes_[k].value = k % 2 == 0 ? DoubleDouble{extrema[k].desired} - correction : DoubleDouble{extrema[k].desired} + correction;
        if (!finite(nodes_[k].value)) return {Status::numerical_failure};
        value_scale_ = std::max(value_scale_, std::abs(nodes_[k].value.hi));
    }
    if (value_scale_ != 0.0) {
        const int exponent = std::ilogb(value_scale_);
        value_scale_ = std::scalbn(1.0, exponent);
        if (exponent != 0) {
            for (auto& node : nodes_) {
                const auto value = scale(node.value, -exponent);
                const auto restored = scale(value, exponent);
                if (!finite(value) || restored.hi != node.value.hi || restored.lo != node.value.lo)
                    return {Status::numerical_failure};
                node.value = value;
            }
        }
    }
    // Omitting the largest weight bounds its evaluation condition by n-1.
    const auto omitted = static_cast<std::size_t>(std::ranges::max_element(nodes_, {},
        [](const Node& node) { return std::abs(node.gamma.hi); }) - nodes_.begin());
    const auto tail = std::span{nodes_}.subspan(omitted);
    std::ranges::rotate(tail, tail.begin() + 1);
    // Removing an interpolation node requires removing its weight factor.
    for (std::size_t k = 0; k + 1 < n; ++k) {
        nodes_[k].gamma = nodes_[k].gamma * sum(nodes_[k].x, -nodes_.back().x);
        if (nodes_[k].gamma.hi == 0.0) return {Status::numerical_failure};
    }
    const auto active = std::span{nodes_}.first(n - 1);
    const double max_gamma = std::abs(std::ranges::max_element(active, {}, [](const Node& node) { return std::abs(node.gamma.hi); })->gamma.hi);
    for (auto& node : active) node.gamma = node.gamma / DoubleDouble{max_gamma};
    ready_ = true;
    const auto last = precise(nodes_.back().x);
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double tolerance = 128 * eps * eps * static_cast<double>(n) * std::max(1.0, std::abs(delta_scaled.hi));
    const auto value_ratio = DoubleDouble{value_scale_} / DoubleDouble{desired_scale};
    for (std::size_t k = 0; k < n; ++k) {
        const auto value = (k == omitted ? last.value : nodes_[k < omitted ? k : k - 1].value) * value_ratio;
        const auto error = (DoubleDouble{extrema[k].desired} / DoubleDouble{desired_scale} - value) /
                           (DoubleDouble{min_weight} / DoubleDouble{extrema[k].weight});
        const auto expected = k % 2 == 0 ? delta_scaled : DoubleDouble{-delta_scaled.hi, -delta_scaled.lo};
        if (last.status != Status::ok || !finite(error) || std::abs((error - expected).hi) > tolerance) {
            ready_ = false;
            return {Status::numerical_failure};
        }
    }
    return {Status::ok, delta.hi + delta.lo};
}

struct Interpolator::Location {
    const Node* node;   // Exact active node, if x is one.
    double distance{};  // Nearest active-node distance; scales barycentric terms toward unity.
};

Interpolator::Location Interpolator::locate(double x) const noexcept {
    const auto active = std::span{nodes_}.first(nodes_.size() - 1);
    const auto next = std::ranges::lower_bound(active, x, std::ranges::greater{}, &Node::x);
    if (next != active.end() && x == next->x) return {&*next};
    double distance = next == active.end() ? std::numeric_limits<double>::infinity() : std::abs(x - next->x);
    if (next != active.begin()) distance = std::min(distance, std::abs(x - (next - 1)->x));
    return {nullptr, distance};
}

// Binary64 barycentric A(x) with its error bound. Each term carries at most 3 eps relative error
// (rounded gamma, value, difference, quotient, product) and compensated sums add eps, so
// |dA| <= 8 eps (M_N + |A| M_D)/|D| with a 2x margin; M_D >= |D| covers rounding A itself.
// The bound is infinite when binary64 cannot resolve x.
Bounded Interpolator::fast(double x) const noexcept {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const auto at = locate(x);
    if (at.node != nullptr) {
        const double value = std::fma(at.node->value.hi, value_scale_, at.node->value.lo*value_scale_);
        return std::isfinite(value) ? Bounded{Status::ok, value, amplitude_bound(value)} : Bounded{Status::numerical_failure};
    }
    if (value_scale_ == 0.0) return {Status::ok};
    Sum numerator, denominator;
    double numerator_magnitude = 0.0, denominator_magnitude = 0.0;
    for (const auto& node : std::span{nodes_}.first(nodes_.size() - 1)) {
        const double term = node.gamma.hi * (at.distance / (x - node.x));
        const double value_term = term * node.value.hi;
        numerator.add(value_term);
        denominator.add(term);
        numerator_magnitude += std::abs(value_term);
        denominator_magnitude += std::abs(term);
    }
    const double scaled = numerator.value() / denominator.value();
    const double bound = 8 * eps * (numerator_magnitude + std::abs(scaled) * denominator_magnitude) /
                         std::abs(denominator.value());
    if (!std::isfinite(scaled) || !std::isfinite(bound)) return {Status::ok, 0.0, std::numeric_limits<double>::infinity()};
    const double value = scaled * value_scale_;
    return std::isfinite(value) ? Bounded{Status::ok, value, amplitude_bound(value, bound*value_scale_)} :
                                 Bounded{Status::numerical_failure};
}

// Two-component barycentric A(x)/value_scale_ with its error bound.
Interpolator::PreciseValue Interpolator::precise(double x) const noexcept {
    const auto at = locate(x);
    if (at.node != nullptr) return {Status::ok, at.node->value};
    if (value_scale_ == 0.0) return {Status::ok};
    const double distance = at.distance;
    DoubleDouble numerator{}, denominator{};
    double magnitude = 0.0, numerator_magnitude = 0.0;
    for (const auto& node : std::span{nodes_}.first(nodes_.size() - 1)) {
        const auto difference = sum(x, -node.x);
        const double first = distance / difference.hi;
        // The scalar numerator needs only an FMA remainder and the denominator's low component.
        const auto quotient = sum(first, (std::fma(-first, difference.hi, distance) - first*difference.lo) / difference.hi);
        const auto term = node.gamma * quotient;
        numerator = numerator + term * node.value;
        denominator = denominator + term;
        magnitude += std::abs(term.hi);
        numerator_magnitude += std::abs(term.hi * node.value.hi);
    }
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const auto n = static_cast<double>(nodes_.size());
    if (!finite(denominator) || std::abs(denominator.hi) <= 128 * eps * eps * n * magnitude)
        return {Status::numerical_failure};
    const auto scaled = numerator / denominator;
    if (!finite(scaled)) return {Status::numerical_failure};
    // Two-component terms and sloppy additions each err by a few eps^2 of their operands, so
    // recursive summation gives |dA| <= 16 (n+4) eps^2 (M_N + |A| M_D)/|D| with margin.
    const double bound = 16 * (n + 4) * eps * eps *
        ((numerator_magnitude + std::abs(scaled.hi) * magnitude) / std::abs(denominator.hi));
    return {Status::ok, scaled, bound};
}

bool Interpolator::usable(double x) const noexcept {
    return ready_ && nodes_.size() >= 2 && std::isfinite(x) && x >= -1.0 && x <= 1.0;
}

// The public evaluations differ only in when the binary64 result is accurate enough.
Value Interpolator::evaluate(double x) const noexcept {
    if (!usable(x)) return {Status::invalid_nodes};
    if (const auto a = fast(x); a.status == Status::ok && a.bound <= 0x1p-42 * std::abs(a.value)) return {Status::ok, a.value};
    const auto a = precise(x);
    return a.status == Status::ok ? rounded(product(a.value, value_scale_, 1.0)) : Value{a.status};
}

Bounded Interpolator::estimate(double x) const noexcept {
    if (!usable(x)) return {Status::invalid_nodes};
    // Node values are normalized to at most 2, so 2^-46 is a few eps of the amplitude scale.
    if (const auto a = fast(x); a.status == Status::ok && a.bound <= 0x1p-46 * value_scale_) return a;
    const auto a = precise(x);
    if (a.status != Status::ok) return {a.status};
    const auto value = rounded(product(a.value, value_scale_, 1.0));
    return value.status == Status::ok ? Bounded{Status::ok, value.value, amplitude_bound(value.value, a.bound*value_scale_)} :
                                       Bounded{value.status};
}

Value Interpolator::weighted_error(const GridPoint& point, double relative_accuracy,
                                  double absolute_accuracy) const noexcept {
    if (!usable(point.x) || !std::isfinite(point.desired)) return {Status::invalid_nodes};
    if (!(std::isfinite(point.weight) && point.weight > 0.0)) return {Status::invalid_weights};
    if (!(relative_accuracy > 0.0 && relative_accuracy < 1.0) || !(absolute_accuracy >= 0.0 && std::isfinite(absolute_accuracy)))
        return {Status::invalid_tolerance};
    // Zero crossings and tiny ripples take the two-component path.
    if (const auto error = certified_error(point, fast(point.x), relative_accuracy, absolute_accuracy))
        return {Status::ok, *error};
    const auto amplitude = precise(point.x);
    if (amplitude.status != Status::ok) return {amplitude.status};
    const double error_scale = std::max(std::abs(point.desired), value_scale_);
    if (error_scale == 0.0) return {Status::ok};
    // Subtract in a common scale, then combine scale and weight before rounding tiny amplitudes.
    const auto difference = DoubleDouble{point.desired}/DoubleDouble{error_scale} -
                            amplitude.value*(DoubleDouble{value_scale_}/DoubleDouble{error_scale});
    bool lost = false;
    const auto value = rounded(product(difference, error_scale, point.weight, &lost));
    if (value.status == Status::ok && lost &&
        std::numeric_limits<double>::denorm_min() > relative_accuracy*std::abs(value.value))
        return {Status::numerical_failure};
    return value;
}

namespace {
struct Minimum {
    Status status;
    double x, value;
};

// Brent's minimizer on [a, b], seeded with the sampled peak (x, fx): parabolic steps through the
// three best points, golden-section steps when the parabola leaves the bracket or stalls.
template <class Function>
[[nodiscard]] Minimum brent_minimize(const Function& f, double a, double b, double x, double fx) noexcept {
    constexpr double golden = 0.38196601125010515;  // (3 - sqrt(5))/2
    // Below sqrt(eps)|x| a smooth peak is flat to rounding; 2^-20 of the bracket bounds the
    // remaining peak error by about 2^-40 of the variation across one grid interval.
    constexpr double root_eps = 0x1p-26;
    const double tolerance = 0x1p-20 * (b - a);
    double w = x, v = x, fw = fx, fv = fx, d = 0.0, e = 0.0;
    for (int iteration = 0; iteration < 64; ++iteration) {
        const double m = std::midpoint(a, b);
        const double tol1 = root_eps * std::abs(x) + tolerance, tol2 = 2 * tol1;
        if (std::abs(x - m) <= tol2 - (b - a) / 2) break;
        bool golden_step = true;
        if (std::abs(e) > tol1) {
            const double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if (q > 0) p = -p; else q = -q;
            const double previous = e;
            e = d;
            if (std::abs(p) < std::abs(q * previous / 2) && p > q * (a - x) && p < q * (b - x)) {
                d = p / q;
                if (x + d - a < tol2 || b - (x + d) < tol2) d = x < m ? tol1 : -tol1;
                golden_step = false;
            }
        }
        if (golden_step) {
            e = (x < m ? b : a) - x;
            d = golden * e;
        }
        const double u = std::abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1);
        const auto fu = f(u);
        if (fu.status != Status::ok) return {fu.status, x, fx};
        if (fu.value <= fx) {
            (u < x ? b : a) = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu.value;
        } else {
            (u < x ? a : b) = u;
            if (fu.value <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu.value;
            } else if (fu.value <= fv || v == x || v == w) {
                v = u; fv = fu.value;
            }
        }
    }
    return {Status::ok, x, fx};
}

struct Series {
    double tail{};   // sum_{k>=1} |c_k|
    double error{};  // Bound on reconstructed amplitude error; infinite disables the series.
    double magnitude{};
};

// The fitted amplitude has degree r, so its r+1 values at y_m = cos(pi m/r) determine it exactly
// (DCT-I): c_k = (2/r) sum''_m A(y_m) cos(pi m k/r), halved for k = 0 and k = r, where sum''
// halves the m = 0 and m = r terms.
[[nodiscard]] Series expand(const Interpolator& polynomial, std::span<const double> cosines,
              std::span<double> samples, std::span<double> c) noexcept {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr Series unusable{0.0, std::numeric_limits<double>::infinity()};
    const auto r = c.size() - 1;
    if (r == 0) {
        const auto value = polynomial.estimate(0.0);
        if (value.status != Status::ok) return unusable;
        c[0] = value.value;
        return {0.0, value.bound, std::abs(value.value)};
    }
    double sample_error = 0.0, sample_size = 0.0, magnitude = 0.0;
    for (std::size_t m = 0; m <= r; ++m) {
        const auto value = polynomial.estimate(cosines[m]);
        if (value.status != Status::ok || !std::isfinite(value.bound)) return unusable;
        samples[m] = value.value;
        sample_error = std::max(sample_error, value.bound);
        sample_size += std::abs(value.value);
        magnitude = std::max(magnitude, std::abs(value.value));
    }
    // Stage even/odd pairs in the existing sample buffer; the middle sample is unpaired.
    const double first = samples[0], last = samples[r];
    samples[0] = 0.5 * (first + last);
    samples[r] = 0.5 * (first - last);
    for (std::size_t m = 1; m <= (r - 1) / 2; ++m) {
        const double a = samples[m], b = samples[r - m];
        samples[m] = a + b;
        samples[r - m] = a - b;
    }
    const double scale = 2.0 / static_cast<double>(r);
    double tail = 0.0, size = 0.0;
    for (std::size_t k = 0; k <= r; ++k) {
        Sum sum;
        // cos(pi*(r-m)*k/r) = (-1)^k cos(pi*m*k/r): pair samples before summing.
        const bool even = k % 2 == 0;
        sum.add(samples[even ? 0 : r]);
        for (std::size_t m = 1, q = k; m <= (r - 1) / 2; ++m) {
            sum.add(samples[even ? m : r - m] * cosines[q <= r ? q : 2 * r - q]);
            q += k;
            if (q >= 2 * r) q -= 2 * r;
        }
        if (r % 2 == 0 && k % 2 == 0) sum.add((k % 4 == 0 ? 1.0 : -1.0) * samples[r / 2]);
        c[k] = (k == 0 || k == r ? 0.5 : 1.0) * scale * sum.value();
        size += std::abs(c[k]);
        if (k != 0) tail += std::abs(c[k]);
    }
    // Sample errors reach A(x) as the interpolant of those errors, bounded by the Lebesgue
    // constant of the r+1 Chebyshev extreme points, (2/pi) ln(r+1) + 1. Each c_k also carries
    // (2/r) 5 eps sum|A_m| from the folded cosine table (argument and value roundings) and the
    // products and pair addition, plus 3 eps |c_k| from compensated summation and scaling. 2x margin.
    const double lebesgue = 2 / pi * std::log(static_cast<double>(r + 1)) + 1;
    const double error = 2 * (lebesgue * sample_error +
                              scale * static_cast<double>(r + 1) * 5 * eps * sample_size + 3 * eps * size);
    return {tail, error, magnitude};
}

// Clenshaw's recurrence b_k = c_k + 2x b_{k+1} - b_{k+2}, A = c_0 + x b_1 - b_2, at four points
// at once: one recurrence is latency-bound (no faster than barycentric), four interleaved are
// about 3.5x faster. A residual at step k is a perturbation of c_k, which enters A multiplied by
// |T_k(x)| <= 1, so |dA| <= 2 eps (sum_{k>=1}|c_k| + 5 sum|b_k| + |c_0| + 2|b_1| + |A|) + series
// error (2x margin).
constexpr std::size_t lanes = 4;

[[nodiscard]] std::array<Bounded, lanes> clenshaw(std::span<const double> c, const Series& series,
                                    const std::array<double, lanes>& x) noexcept {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    std::array<double, lanes> b1{}, b2{}, magnitude{};
    for (std::size_t k = c.size() - 1; k >= 1; --k) {
        for (std::size_t l = 0; l < lanes; ++l) {
            const double b = (c[k] + (2 * x[l]) * b1[l]) - b2[l];
            b2[l] = b1[l];
            b1[l] = b;
            magnitude[l] += std::abs(b);
        }
    }
    std::array<Bounded, lanes> result{};
    for (std::size_t l = 0; l < lanes; ++l) {
        const double value = (c[0] + x[l] * b1[l]) - b2[l];
        result[l] = {Status::ok, value, 2 * eps * (series.tail + 5 * magnitude[l] + std::abs(c[0]) +
                                                   2 * std::abs(b1[l]) + std::abs(value)) + series.error};
    }
    return result;
}

// Alternation points of best approximations on a set E follow its equilibrium measure (Kadec);
// equispaced references have exponentially growing Lebesgue constants, so their levelled error
// can start at rounding level and never resolve the sign pattern. For E a union of k intervals
// with endpoints e(0) > ... > e(2k-1) in x = cos(omega), the density is |q(x)|/(pi*sqrt|R(x)|),
// R = prod (x - e(j)), q monic with one root per gap where the gap's mass vanishes (Akhiezer).
// Interval i (even: band i/2, odd: a gap) maps by x = c + h*cos(t), absorbing its own endpoint
// factors: dx/sqrt((e(i)-x)(x-e(i+1))) = dt, so every integrand below is smooth in t.
// Writes reference[] from grid points and roots[0..k-2]; false leaves the caller to seed otherwise.
[[nodiscard]] bool equilibrium_reference(std::span<const GridPoint> points, std::span<const BandRange> bands,
                                         std::span<GridPoint> reference, std::span<double> roots) noexcept {
    const auto k = bands.size(), n = reference.size();
    if (k == 0 || n < k || roots.size() + 1 < k ||
        std::ranges::any_of(bands, [](const BandRange& b) { return b.end - b.begin < 2; })) return false;
    roots = roots.first(k - 1);
    const auto endpoint = [&](std::size_t j) noexcept {
        const auto& b = bands[j / 2];
        return points[j % 2 == 0 ? b.begin : b.end - 1].x;
    };
    // |q(x)|/sqrt|R(x)| without interval i's endpoints, and without root `skip` when skip < k-1.
    const auto weight = [&](double x, std::size_t i, std::size_t skip) noexcept {
        double value = 1.0;
        for (std::size_t g = 0; g < roots.size(); ++g) if (g != skip) value *= x - roots[g];
        for (std::size_t j = 0; j < 2 * k; ++j) if (j != i && j != i + 1) value /= std::sqrt(std::abs(x - endpoint(j)));
        return std::abs(value);
    };
    // Gauss-Chebyshev nodes on interval i; the measure only steers the seed, so 64 suffice.
    constexpr std::size_t nodes = 64;
    const auto node = [&](std::size_t i, std::size_t l) noexcept {
        const double t = pi * (static_cast<double>(l) + 0.5) / static_cast<double>(nodes);
        return std::midpoint(endpoint(i), endpoint(i + 1)) + (endpoint(i) - endpoint(i + 1)) / 2 * std::cos(t);
    };
    // Gauss-Seidel over the gaps: with the other roots fixed, the gap condition is linear in this
    // root, which is the weighted mean of x over its gap. One sweep is exact for two bands.
    for (std::size_t g = 0; g < roots.size(); ++g) roots[g] = std::midpoint(endpoint(2 * g + 1), endpoint(2 * g + 2));
    for (std::size_t sweep = 0; sweep < 64; ++sweep) {
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
    const auto band_mass = [&](std::size_t b) noexcept {
        double mass = 0.0;
        for (std::size_t l = 0; l < nodes; ++l) mass += weight(node(2 * b, l), 2 * b, k);
        return mass;
    };
    double total = 0.0;
    for (std::size_t b = 0; b < k; ++b) total += band_mass(b);
    if (!(total > 0.0 && std::isfinite(total))) return false;
    // One point per band plus a cumulative-floor share of the rest by mass; each band then places
    // its points at equal steps of its own measure, endpoints included, on the trapezoid CDF in t.
    const auto spare = n - k;
    std::size_t out = 0, used = 0;
    double before = 0.0;
    for (std::size_t b = 0; b < k; ++b) {
        before += band_mass(b);
        const auto share = b + 1 == k ? spare : std::min(spare, static_cast<std::size_t>(static_cast<double>(spare) * before / total));
        const auto m = 1 + share - used;
        used = share;
        const auto band = points.subspan(bands[b].begin, bands[b].end - bands[b].begin);
        if (m > band.size()) return false;
        const double c = std::midpoint(band.front().x, band.back().x), h = (band.front().x - band.back().x) / 2;
        const auto density = [&](std::size_t j) noexcept { return weight(band[j].x, 2 * b, k); };
        const auto angle = [&](std::size_t j) noexcept { return std::acos(std::clamp((band[j].x - c) / h, -1.0, 1.0)); };
        double area = 0.0;
        for (std::size_t j = 1; j < band.size(); ++j) area += (density(j - 1) + density(j)) / 2 * (angle(j) - angle(j - 1));
        if (!(area > 0.0 && std::isfinite(area))) return false;
        const auto target = [&](std::size_t j) noexcept {
            return m == 1 ? area / 2 : area * static_cast<double>(j) / static_cast<double>(m - 1);
        };
        std::size_t placed = 0, previous = 0;
        const auto place = [&](std::size_t j) noexcept {
            // Distinct increasing indices with room left for the band's remaining points.
            j = std::min(std::max(j, placed == 0 ? 0 : previous + 1), band.size() - (m - placed));
            reference[out++] = band[j];
            previous = j;
            ++placed;
        };
        double cdf = 0.0, last = 0.0;
        for (std::size_t j = 0; j < band.size() && placed < m; ++j) {
            if (j != 0) cdf += (density(j - 1) + density(j)) / 2 * (angle(j) - angle(j - 1));
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
    : Exchange(Specification{spec.taps % 2 == 0 ? 0 : spec.taps,
        bands(spec),
        Symmetry::symmetric, Response::constant, spec.density}, max_points) {}

Exchange::Exchange(const Specification& spec, std::size_t max_points)
    : grid_(make_grid(spec, max_points)),
      polynomial_(grid_.status == Status::ok ? grid_.degree + 2 : 0),
      extrema_(grid_.status == Status::ok ? grid_.degree + 2 : 0),
      peaks_(grid_.status == Status::ok ? grid_.peak_capacity : 0),
      chebyshev_(grid_.status == Status::ok ? grid_.workspace : 0),
      taps_(spec.taps), symmetry_(spec.symmetry), response_(spec.response), type_(form(spec)) {
    if (chebyshev_.empty()) return;
    min_weight_ = std::ranges::min_element(grid_.points, {}, &GridPoint::weight)->weight;
    const auto d = grid_.degree;
    if (d == 0) { chebyshev_[0] = 1.0; return; }
    for (std::size_t q = 0; q <= d; ++q)
        chebyshev_[q] = 2*q <= d ? std::cos(pi * static_cast<double>(q) / static_cast<double>(d)) :
                                -std::cos(pi * static_cast<double>(d-q) / static_cast<double>(d));
}

Status Exchange::recover(std::span<const double> c, bool certified) noexcept {
    const auto d = grid_.degree, L = d + 1;
    const auto h = std::span{chebyshev_}.subspan(L, taps_);
    if (certified) {
        const auto coefficient = [&](std::size_t j) { return j < c.size() ? c[j] : 0.0; };
        const auto m = taps_ / 2;
        if (type_ == 1) {
            h[m] = c[0];
            for (std::size_t j = 1; j <= m; ++j) h[m-j] = h[m+j] = 0.5*c[j];
        } else if (type_ == 3) {
            h[m] = 0.0;
            for (std::size_t j = 1; j <= m; ++j) {
                const double b = j == 1 ? c[0]-0.5*coefficient(2) : 0.5*(coefficient(j-1)-coefficient(j+1));
                h[m-j] = 0.5*b;
                h[m+j] = -h[m-j];
            }
        } else {
            for (std::size_t j = 0; j < m; ++j) {
                const double sign = type_ == 2 ? 1.0 : -1.0;
                const double b = j == 0 ? c[0]+sign*0.5*coefficient(1) : 0.5*(c[j]+sign*coefficient(j+1));
                h[m-1-j] = 0.5*b;
                h[m+j] = type_ == 2 ? h[m-1-j] : -h[m-1-j];
            }
        }
        return std::ranges::all_of(h, [](double a) { return std::isfinite(a); }) ? Status::ok : Status::numerical_failure;
    }
    return polynomial_.taps(h, symmetry_, taps_);
}

Status Interpolator::taps(std::span<double> output, Symmetry symmetry, std::size_t count) const noexcept {
    if (!ready_ || nodes_.size() < 2) return Status::invalid_nodes;
    if ((symmetry != Symmetry::symmetric && symmetry != Symmetry::antisymmetric) || count == 0 ||
        (symmetry == Symmetry::antisymmetric && count < 2)) return Status::invalid_taps;
    const auto degree = symmetry == Symmetry::symmetric ? (count-1)/2 : count/2-1;
    if (degree != nodes_.size()-2) return Status::invalid_taps;
    if (output.size() < count) return Status::insufficient_storage;
    const auto type = form(symmetry, count);
    const auto h = output.first(count);
    if (type == 1) return this->taps(output);
    // Real inverse DFT of H=i*exp(-i*alpha*w)*A for antisymmetry, exp(-i*alpha*w)*A otherwise.
    const auto m = count/2, start = count%2 == 0 ? m-1 : m;
    const auto samples = h.subspan(start, m+1);
    double magnitude = 0.0;
    for (std::size_t j = 0; j <= m; ++j) {
        const double omega = 2*pi*static_cast<double>(j)/static_cast<double>(count);
        const auto a = estimate(std::cos(omega));
        if (a.status != Status::ok || !std::isfinite(a.bound)) return Status::numerical_failure;
        const auto value = rounded(product(DoubleDouble{a.value}, factor(type, omega), 1.0));
        if (value.status != Status::ok) return value.status;
        samples[j] = value.value;
        magnitude = std::max(magnitude, std::abs(value.value));
    }
    const int exponent = magnitude == 0.0 ? 0 : std::ilogb(magnitude);
    for (auto& a : samples) {
        const double normalized = std::scalbn(a, -exponent);
        if (std::scalbn(normalized, exponent) != a) return Status::numerical_failure;
        a = normalized;
    }
    const auto tap = [&](std::size_t n) noexcept {
        Sum sum;
        if (symmetry == Symmetry::symmetric) sum.add(samples[0]);
        if (count%2 == 0 && type == 4) sum.add((m+1+n)%2 == 0 ? samples[m] : -samples[m]);
        const auto step = count-1-2*n;
        for (std::size_t j = 1, q = step; j <= (count-1)/2; ++j) {
            const double angle = pi * static_cast<double>(std::min(q, 2*count-q))/static_cast<double>(count);
            const double wave = symmetry == Symmetry::symmetric ? std::cos(angle) : (q <= count ? std::sin(angle) : -std::sin(angle));
            sum.add(2*samples[j]*wave);
            q += step;
            if (q >= 2*count) q -= 2*count;
        }
        return std::scalbn(sum.value()/static_cast<double>(count), exponent);
    };
    const double held = count%2 == 0 ? tap(m-1) : 0.0;
    for (std::size_t n = 0; n < start; ++n) h[n] = tap(n);
    h[start] = held;
    for (std::size_t n = 0; n < m; ++n) h[count-1-n] = symmetry == Symmetry::symmetric ? h[n] : -h[n];
    return std::ranges::all_of(h, [](double a) { return std::isfinite(a); }) ? Status::ok : Status::numerical_failure;
}

ExchangeResult Exchange::run(std::size_t max_iterations, double tolerance) noexcept {
    ready_ = false;
    exact_ = false;
    tap_status_ = Status::invalid_nodes;
    if (grid_.status != Status::ok) return {grid_.status};
    if (extrema_.empty()) return {Status::invalid_nodes};
    if (!(std::isfinite(tolerance) && tolerance > 0.0 && tolerance < 1.0))
        return {Status::invalid_tolerance};
    if (max_iterations == 0) return {Status::nonconvergence};
    const auto n = extrema_.size();
    const auto L = grid_.degree + 1;
    const auto workspace = std::span{chebyshev_};
    const auto scratch = workspace.subspan(L, grid_.peak_capacity);
    const auto coefficients = workspace.last(L);
    const bool zero = std::ranges::all_of(grid_.bands, [](const BandRange& b) { return b.band.target == 0.0; });
    constant_ = grid_.bands.front().band.target;
    if (zero || (type_ == 1 && response_ == Response::constant && std::ranges::all_of(grid_.bands,
        [&](const BandRange& b) { return b.band.target == constant_; }))) {
        if (zero) constant_ = 0.0;
        const auto h = workspace.subspan(L, taps_);
        std::ranges::fill(h, 0.0);
        if (type_ == 1) h[taps_/2] = constant_;
        exact_ = ready_ = true;
        tap_status_ = Status::ok;
        return {Status::ok};
    }
    // Equilibrium seeds give every band a point. With more bands than points, uniform grid indices
    // may all land on one target level, which would make the levelled error identically zero.
    if (!equilibrium_reference(grid_.points, grid_.bands, extrema_, scratch)) {
        const auto intervals_count = n - 1;
        const auto step = (grid_.points.size()-1)/intervals_count;
        const auto extra = (grid_.points.size()-1)%intervals_count;
        std::size_t index = 0, remainder = 0;
        for (auto& point : extrema_) {
            point = grid_.points[index];
            index += step;
            remainder += extra;
            if (remainder >= intervals_count) { ++index; remainder -= intervals_count; }
        }
        const auto [low, high] = std::ranges::minmax_element(extrema_, {}, &GridPoint::desired);
        if (high->desired-low->desired <= 32*std::numeric_limits<double>::epsilon()*std::max(std::abs(low->desired), std::abs(high->desired))) {
            const auto [minimum, maximum] = std::ranges::minmax_element(grid_.points, {}, &GridPoint::desired);
            const auto different = std::abs(minimum->desired-low->desired) > std::abs(maximum->desired-low->desired) ? minimum : maximum;
            if (different->desired != low->desired) extrema_[n/2] = *different;
            std::ranges::sort(extrema_, {}, &GridPoint::omega);
        }
    }
    if (n == 2 && type_ == 1) {
        double largest = -1.0;
        for (std::size_t i = 0; i < grid_.bands.size(); ++i) for (std::size_t j = i+1; j < grid_.bands.size(); ++j) {
            const auto& a = grid_.bands[i].band;
            const auto& b = grid_.bands[j].band;
            const double low = std::min(a.weight, b.weight), high = std::max(a.weight, b.weight);
            const double gain = std::max(std::abs(a.target),std::abs(b.target));
            const auto difference = gain == 0.0 ? DoubleDouble{} : DoubleDouble{a.target}/DoubleDouble{gain}-DoubleDouble{b.target}/DoubleDouble{gain};
            const auto candidate = rounded(product(difference/(DoubleDouble{1.0}+DoubleDouble{low}/DoubleDouble{high}), gain, low));
            if (candidate.status != Status::ok) return {candidate.status};
            if (std::abs(candidate.value) > largest) {
                largest = std::abs(candidate.value);
                extrema_[0] = grid_.points[grid_.bands[i].begin];
                extrema_[1] = grid_.points[grid_.bands[j].begin];
            }
        }
    }
    const auto band_index = [&](double omega) noexcept {
        const auto found = std::ranges::lower_bound(grid_.bands, omega, {}, [](const BandRange& b) { return b.band.upper; });
        return static_cast<std::size_t>(found-grid_.bands.begin());
    };
    // Convergence compares errors at relative resolution `tolerance`; certifying each error to
    // tolerance/256 keeps max_error, flatness and peak choice within 1% of that resolution while
    // letting well-conditioned binary64 evaluations skip the two-component path.
    const double accuracy = tolerance * 0x1p-8;
    double rounding_floor = 0.0;
    if (type_ != 1) for (const auto& range : grid_.bands) {
        const auto& band = range.band;
        const auto scale_value = scaled_ratio(band.weight,
            response_ == Response::differentiator && band.target != 0.0 ? 1.0 : std::abs(band.target), 1.0);
        if (scale_value.status != Status::ok) return {scale_value.status};
        rounding_floor = std::max(rounding_floor, 16*std::numeric_limits<double>::epsilon()*scale_value.value + std::numeric_limits<double>::denorm_min());
    }
    ExchangeResult result{Status::nonconvergence};
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        result.iterations = iteration + 1;
        const auto fit = polynomial_.fit(extrema_);
        if (fit.status != Status::ok) { result.status = fit.status; return result; }
        result.ripple = std::abs(fit.value);
        result.max_error = 0.0;
        // Chebyshev form of this fit; points its bound cannot serve use the barycentric path.
        const auto series = expand(polynomial_, workspace.first(L), scratch.first(L), coefficients);
        // Errors near peaks are about the ripple, so a series floor above accuracy*ripple even at
        // the smaller band weight would certify almost nothing in the fully certified scan.
        const bool chebyshev = min_weight_ * series.error <= accuracy * result.ripple;
        // max_error >= ripple (the fitted nodes sit at +-ripple) and convergence compares errors
        // only against it, so an error is good enough once its sign is certain and its bound is
        // within accuracy*max(|E|, ripple). Points near zero crossings then stay in binary64.
        // Errors are measured in the transformed problem the fit solves, W*q*(D/q - P), which equals
        // the physical W*(D - q*P) up to rounding q (about eps*W*|D|, checked at convergence).
        const double floor = accuracy * result.ripple;
        const auto error_at = [&](const GridPoint& point) noexcept {
            return polynomial_.weighted_error(point, accuracy, floor);
        };
        const auto series_error_at = [&](const GridPoint& point, const Bounded& amplitude) noexcept {
            const auto error = certified_error(point, amplitude, accuracy, floor);
            return error ? Value{Status::ok, *error} : error_at(point);
        };
        // Each sign run contributes its largest peak, including band endpoints. Capacity includes
        // band endpoints and type-dependent stationary points; differentiators use a grid bound.
        // Only peaks need the accuracy above: the lazy scan ranks points by binary64 estimates with
        // rigorous bounds (Clenshaw while its bound is far below the ripple, else barycentric) and
        // certifies uncertain signs and each run's two largest estimates. If any other point of a run
        // could still exceed its certified peak, the scan restarts fully certified, so the peaks and
        // max_error equal those of certifying every point.
        struct Candidate { const GridPoint* point; std::size_t index; double error{}, bound{}; };  // bound 0: certified
        std::size_t count = 0, held = 0;
        std::array<Candidate, 2> top{};
        double rest = 0.0;  // Largest |error| + bound of the run's other points.
        bool ambiguous = false;
        const auto close_run = [&]() noexcept -> Status {
            if (held == 0) return Status::ok;
            for (std::size_t i = 0; i < held; ++i) {
                auto& c = top[i];
                if (c.bound == 0.0 || (i == 1 && std::abs(c.error) + c.bound < std::abs(top[0].error))) continue;
                const auto value = error_at(*c.point);
                if (value.status != Status::ok) return value.status;
                ambiguous |= std::signbit(value.value) != std::signbit(c.error) || value.value == 0.0;
                c.error = value.value;
                c.bound = 0.0;
            }
            const bool second = held == 2 && std::abs(top[1].error) > std::abs(top[0].error);
            const auto& best = top[second ? 1 : 0];
            const auto& other = top[second ? 0 : 1];
            ambiguous |= std::max(rest, held == 2 ? std::abs(other.error) + other.bound : 0.0) > std::abs(best.error);
            if (count == peaks_.size()) return Status::numerical_failure;
            peaks_[count++] = {best.index, best.point->omega, best.error};
            result.max_error = std::max(result.max_error, std::abs(best.error));
            held = 0;
            rest = 0.0;
            return Status::ok;
        };
        const auto visit = [&](Candidate c) noexcept -> Status {
            if (c.bound != 0.0 && !(c.bound < std::abs(c.error))) {
                const auto value = error_at(*c.point);
                if (value.status != Status::ok) return value.status;
                c = {c.point, c.index, value.value, 0.0};
            }
            if (c.error == 0.0) return Status::ok;
            if (held != 0 && std::signbit(c.error) != std::signbit(top[0].error))
                if (const auto status = close_run(); status != Status::ok) return status;
            const auto larger = [](const Candidate& a, const Candidate& b) { return std::abs(a.error) > std::abs(b.error); };
            if (held < 2) top[held++] = c;
            else if (larger(c, top[1])) { rest = std::max(rest, std::abs(top[1].error) + top[1].bound); top[1] = c; }
            else rest = std::max(rest, std::abs(c.error) + c.bound);
            if (held == 2 && larger(top[1], top[0])) std::swap(top[0], top[1]);
            return Status::ok;
        };
        const auto estimate_at = [&](const GridPoint& point, std::size_t index, const Bounded& amplitude) noexcept {
            // Estimates accurate enough to certify need no second evaluation.
            if (const auto error = certified_error(point, amplitude, accuracy, floor)) return Candidate{&point, index, *error, 0.0};
            const auto [error, bound] = bounded_error(point, amplitude);
            // A failed or unbounded estimate gets an infinite bound, which forces certification.
            return Candidate{&point, index, error, amplitude.status == Status::ok && std::isfinite(error) && std::isfinite(bound) ?
                                                  bound : std::numeric_limits<double>::infinity()};
        };
        const bool estimates = min_weight_ * series.error <= result.ripple / 16;
        // Grid points go through Clenshaw `lanes` at a time; previous extrema interleave singly.
        const auto scan = [&](bool lazy) noexcept -> Status {
            count = held = 0;
            rest = result.max_error = 0.0;
            ambiguous = false;
            const auto size = grid_.points.size();
            std::size_t extremum = 0;
            for (std::size_t base = 0; base < size; base += lanes) {
                std::array<Bounded, lanes> amplitudes{};
                if (lazy ? estimates : chebyshev) {
                    std::array<double, lanes> x{};
                    for (std::size_t l = 0; l < lanes; ++l) x[l] = grid_.points[std::min(base + l, size - 1)].x;
                    amplitudes = clenshaw(coefficients, series, x);
                }
                for (std::size_t l = 0; l < lanes && base + l < size; ++l) {
                    const auto j = base + l;
                    while (extremum < n && extrema_[extremum].omega <= grid_.points[j].omega) {
                        const auto& point = extrema_[extremum++];
                        if (point.omega == grid_.points[j].omega) continue;
                        const auto nearest = j != 0 && point.omega-grid_.points[j-1].omega < grid_.points[j].omega-point.omega ? j-1 : j;
                        Candidate c{&point, nearest};
                        if (lazy) c = estimate_at(point, nearest, polynomial_.fast(point.x));
                        else if (const auto value = error_at(point); value.status == Status::ok) c.error = value.value;
                        else return value.status;
                        if (const auto status = visit(c); status != Status::ok) return status;
                    }
                    const auto& point = grid_.points[j];
                    Candidate c{&point, j};
                    if (lazy) c = estimate_at(point, j, estimates ? amplitudes[l] : polynomial_.fast(point.x));
                    else if (const auto value = chebyshev ? series_error_at(point, amplitudes[l]) : error_at(point);
                             value.status == Status::ok) c.error = value.value;
                    else return value.status;
                    if (const auto status = visit(c); status != Status::ok) return status;
                }
            }
            return close_run();
        };
        // When Clenshaw certifies the whole iteration, ranking first would only add bookkeeping.
        auto status = scan(!chebyshev);
        if (status == Status::ok && ambiguous) status = scan(false);
        if (status != Status::ok) { result.status = status; return result; }
        if (count < n) { result.status = Status::numerical_failure; return result; }
        for (auto& peak : std::span{peaks_}.first(count)) {
            const auto& point = grid_.points[peak.index];
            const auto& range = grid_.bands[band_index(point.omega)];
            const auto begin = range.begin, end = range.end;
            const double lower = peak.index > begin ? grid_.points[peak.index-1].omega : range.band.lower;
            const double upper = peak.index+1 < end ? grid_.points[peak.index+1].omega : range.band.upper;
            const double sign = std::signbit(peak.error) ? -1.0 : 1.0;
            // Minimize -sign*E inside this band; the bracket never crosses the transition.
            // Type I targets are the band constants; other forms divide by q(omega) at each probe.
            const auto negated_error = [&](double omega) noexcept -> Value {
                const auto target = type_ == 1 ? Transformed{Status::ok, {omega, std::cos(omega), point.desired, point.weight}}
                                               : transform(type_, response_, range.band, omega);
                if (target.status != Status::ok) return {target.status};
                const auto value = error_at(target.point);
                return {value.status, -sign * value.value};
            };
            // Seed with the recorded peak, which may be a previous extremum nearest this grid index.
            const auto refined = brent_minimize(negated_error, lower, upper, peak.omega, -sign * peak.error);
            if (refined.status != Status::ok) { result.status = refined.status; return result; }
            if (-refined.value > sign * peak.error) { peak.omega = refined.x; peak.error = -sign * refined.value; }
            result.max_error = std::max(result.max_error, std::abs(peak.error));
        }
        if (detail::select_peaks(std::span{peaks_}.first(count), n, scratch) != n) {
            result.status = Status::numerical_failure; return result;
        }
        const auto selected = std::span{peaks_}.first(n);
        // Opposite-sign peaks are separated by a zero; refined peaks that meet or cross mean the
        // error is unresolved in binary64, not that the caller supplied unordered nodes.
        if (std::ranges::adjacent_find(selected, std::ranges::greater_equal{}, &detail::Peak::omega) != selected.end()) {
            result.status = Status::numerical_failure;
            return result;
        }
        const double minimum = std::abs(std::ranges::min_element(selected, {},
            [](const detail::Peak& peak) { return std::abs(peak.error); })->error);
        result.flatness = (result.max_error - minimum) / result.max_error;
        if (result.flatness <= tolerance && result.ripple / result.max_error >= 1.0 / (1.0 + tolerance)) {
            // Physical and transformed errors differ by rounding q and D/q, at most rounding_floor
            // per point, so the physical spread can exceed the transformed one by 2*floor/max_error.
            // More iterations cannot shrink it: report it, and fail if it breaks the tolerance.
            result.flatness += 2 * rounding_floor / result.max_error;
            if (result.flatness > tolerance) { result.status = Status::numerical_failure; return result; }
            tap_status_ = recover(coefficients, std::isfinite(series.error) && series.error <=
                8 * std::numeric_limits<double>::epsilon() * static_cast<double>(taps_) * series.magnitude);
            ready_ = true;
            result.status = Status::ok;
            return result;
        }
        // Type I targets and weights are the band constants; the transform's exact scaling costs
        // about 5% of a 15-tap run.
        for (std::size_t k = 0; k < n; ++k) {
            const double omega = selected[k].omega;
            const auto p = type_ == 1 ? Transformed{Status::ok, {omega, std::cos(omega), grid_.points[selected[k].index].desired,
                                                                 grid_.points[selected[k].index].weight}}
                                      : transform(type_, response_, grid_.bands[band_index(omega)].band, omega);
            if (p.status != Status::ok) { result.status = p.status; return result; }
            extrema_[k] = p.point;
        }
    }
    return result;
}

Value Exchange::evaluate(double x) const noexcept {
    if (!ready_ || grid_.points.empty() || !std::isfinite(x) || x < -1.0 || x > 1.0) return {Status::invalid_nodes};
    return exact_ ? Value{Status::ok, constant_} : polynomial_.evaluate(x);
}

Value Exchange::amplitude(double omega) const noexcept {
    if (!ready_ || grid_.points.empty() || !std::isfinite(omega) || omega < 0.0 || omega > pi) return {Status::invalid_nodes};
    if (exact_) return {Status::ok, constant_};
    if (forced(type_, omega)) return {Status::ok};
    if (type_ == 1) return polynomial_.evaluate(std::cos(omega));
    const auto a = polynomial_.precise(std::cos(omega));
    if (a.status != Status::ok) return {a.status};
    // Q=w*(Q/w) keeps half-angle underflow from destroying a representable amplitude.
    return symmetry_ == Symmetry::antisymmetric ? rounded(product(a.value * DoubleDouble{ratio(type_, omega)}, polynomial_.value_scale_, omega)) :
        rounded(product(a.value, polynomial_.value_scale_, factor(type_, omega)));
}

std::span<const GridPoint> Exchange::extrema() const noexcept {
    return ready_ && !grid_.points.empty() && !exact_ ? std::span<const GridPoint>{extrema_} : std::span<const GridPoint>{};
}

Status Exchange::taps(std::span<double> output) const noexcept {
    if (!ready_ || grid_.points.empty()) return Status::invalid_nodes;
    if (output.size() < taps_) return Status::insufficient_storage;
    if (tap_status_ != Status::ok) return tap_status_;
    std::ranges::copy(std::span{chebyshev_}.subspan(grid_.degree+1, taps_), output.begin());
    return Status::ok;
}

Status Interpolator::taps(std::span<double> output) const noexcept {
    if (!ready_ || nodes_.size() < 2) return Status::invalid_nodes;
    // A(w) = h[r] + 2 sum_k h[r+k] cos(k w) has degree r, so N = 2r+1 samples w_m = 2 pi m/N invert
    // it exactly: h[r+k] = (A(0) + 2 sum_m A(w_m) cos(w_m k))/N.
    const auto r = nodes_.size() - 2, N = 2 * r + 1;
    if (output.size() < N) return Status::insufficient_storage;
    // Reducing m*k mod N keeps every cosine argument in [0, pi] and exact up to one rounding.
    const auto cosine = [N](std::size_t q) noexcept {
        q %= N;
        return std::cos(2 * pi * static_cast<double>(std::min(q, N - q)) / static_cast<double>(N));
    };
    // Samples occupy output[r..2r]; taps r-1..0 fill the free lower half, the centre is held
    // until the samples are no longer needed, then the upper half mirrors the lower.
    double magnitude = 0.0;
    for (std::size_t m = 0; m <= r; ++m) {
        const auto amplitude = estimate(cosine(m));
        if (amplitude.status != Status::ok) return amplitude.status;
        if (!std::isfinite(amplitude.bound)) return Status::numerical_failure;
        output[r + m] = amplitude.value;
        magnitude = std::max(magnitude, std::abs(amplitude.value));
    }
    // Exact binary scaling keeps the cosine sum finite even when unscaled 2*A would overflow.
    const int exponent = magnitude == 0.0 ? 0 : std::ilogb(magnitude);
    const auto samples = output.subspan(r, r + 1);
    for (auto& value : samples) {
        const double normalized = std::scalbn(value, -exponent);
        if (std::scalbn(normalized, exponent) != value) return Status::numerical_failure;
        value = normalized;
    }
    const auto tap = [&](std::size_t k) noexcept {
        Sum sum;
        sum.add(samples[0]);
        if (k == 0) {
            for (std::size_t m = 1; m <= r; ++m) sum.add(2 * samples[m]);
            return std::scalbn(sum.value() / static_cast<double>(N), exponent);
        }
        const double angle = 2 * pi * static_cast<double>(k) / static_cast<double>(N);
        const double c = std::cos(angle), s = std::sin(angle);
        double cm = c, sm = s;
        for (std::size_t m = 1, q = k; m <= r; ++m) {
            // Short rotations avoid a cosine per term; exact-index restarts limit drift.
            if (m % 32 == 0) {
                const double phase = 2 * pi * static_cast<double>(q) / static_cast<double>(N);
                cm = std::cos(phase);
                sm = std::sin(phase);
            } else if (m != 1) {
                const double next = cm * c - sm * s;
                sm = sm * c + cm * s;
                cm = next;
            }
            sum.add(2 * samples[m] * cm);
            q += k;
            if (q >= N) q -= N;
        }
        return std::scalbn(sum.value() / static_cast<double>(N), exponent);
    };
    for (std::size_t k = 1; k <= r; ++k) output[r - k] = tap(k);
    output[r] = tap(0);
    std::ranges::reverse_copy(output.first(r), output.subspan(r + 1).begin());
    return std::ranges::all_of(output.first(N), [](double h) { return std::isfinite(h); }) ? Status::ok : Status::numerical_failure;
}

}

