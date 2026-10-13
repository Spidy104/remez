#include <remez.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace remez {
namespace {

constexpr double pi = std::numbers::pi_v<double>;

[[nodiscard]] constexpr unsigned form(Symmetry symmetry, std::size_t taps) noexcept {
    if (symmetry == Symmetry::symmetric) return taps % 2 == 0 ? 2u : 1u;
    return taps % 2 == 0 ? 4u : 3u;
}

[[nodiscard]] constexpr bool forced(unsigned type, double omega) noexcept {
    return (omega == 0.0 && type >= 3) || (omega == pi && (type == 2 || type == 3));
}

[[nodiscard]] double factor(unsigned type, double omega) noexcept {
    if (forced(type, omega)) return 0.0;
    switch (type) {
    case 1: return 1.0;
    case 2: return std::cos(omega / 2);
    case 3: return std::sin(omega);
    default: return std::sin(omega / 2);
    }
}

[[nodiscard]] double factor_over_omega(unsigned type, double omega) noexcept {
    const double t = type == 3 ? omega : omega / 2;
    const double scale = type == 3 ? 1.0 : 0.5;
    return t == 0.0 ? scale : scale * std::sin(t) / t;
}

[[nodiscard]] constexpr bool omitted(unsigned type, Response response, const Band& band, double omega) noexcept {
    const bool slope_at_dc = omega == 0.0 && response == Response::differentiator && band.target != 0.0;
    return forced(type, omega) && !slope_at_dc;
}

constexpr std::size_t element_ceiling = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / 64;

[[nodiscard]] std::size_t intervals(const Band& band, double step, std::size_t limit) noexcept {
    const double count = std::max(1.0, std::ceil((band.upper - band.lower) / step));
    return count < static_cast<double>(limit) ? static_cast<std::size_t>(count) : limit;
}

[[nodiscard]] std::size_t degree_of(const Specification& spec) noexcept {
    return spec.symmetry == Symmetry::symmetric ? (spec.taps - 1) / 2 : spec.taps / 2 - 1;
}

[[nodiscard]] double grid_step(const Specification& spec) noexcept {
    return pi / static_cast<double>(spec.density * (degree_of(spec) + 1));
}

[[nodiscard]] bool valid_enums(const Specification& spec) noexcept {
    const bool symmetry = spec.symmetry == Symmetry::symmetric || spec.symmetry == Symmetry::antisymmetric;
    const bool response = spec.response == Response::constant || spec.response == Response::differentiator;
    return symmetry && response && (spec.response == Response::constant || spec.symmetry == Symmetry::antisymmetric);
}

struct Layout {
    Status status;
    std::size_t points{};
};

[[nodiscard]] Layout layout(const Specification& spec, std::size_t max_points) noexcept {
    max_points = std::min(max_points, element_ceiling);
    if (!valid_enums(spec)) return {Status::invalid_specification};
    if (spec.taps == 0 || spec.taps > max_points || (spec.symmetry == Symmetry::antisymmetric && spec.taps < 2))
        return {Status::invalid_taps};
    if (spec.bands.empty() || spec.bands.size() > max_points) return {Status::invalid_bands};
    if (spec.density == 0) return {Status::invalid_density};
    const auto type = form(spec.symmetry, spec.taps);
    const auto degree = degree_of(spec);
    if (spec.density > max_points / (degree + 1)) return {Status::grid_limit};
    const double step = grid_step(spec);
    std::size_t points = 0;
    double previous = -1.0;
    for (const auto& band : spec.bands) {
        const bool ordered = std::isfinite(band.lower) && std::isfinite(band.upper) && previous < band.lower &&
                             band.lower < band.upper && band.upper <= pi;
        if (!ordered) return {Status::invalid_bands};
        if (!(std::isfinite(band.weight) && band.weight > 0.0)) return {Status::invalid_weights};
        if (!std::isfinite(band.target)) return {Status::invalid_target};
        const bool lower_omitted = omitted(type, spec.response, band, band.lower);
        const bool upper_omitted = omitted(type, spec.response, band, band.upper);
        if (band.target != 0.0 && (lower_omitted || forced(type, band.upper))) return {Status::incompatible_endpoint};
        const auto count = intervals(band, step, max_points) + 1 - static_cast<std::size_t>(lower_omitted) -
                           static_cast<std::size_t>(upper_omitted);
        if (count > max_points - points) return {Status::grid_limit};
        points += count;
        previous = band.upper;
    }
    if (points < degree + 2) return {Status::invalid_density};
    return {Status::ok, points};
}

[[nodiscard]] std::array<Band, 2> lowpass_bands(const Lowpass& spec) noexcept {
    return {Band{0.0, spec.pass_edge, 1.0, spec.weights[0]}, Band{spec.stop_edge, pi, 0.0, spec.weights[1]}};
}

[[nodiscard]] Specification lowpass_specification(const Lowpass& spec, std::span<const Band> bands) noexcept {
    const std::size_t taps = spec.taps % 2 == 0 ? 0 : spec.taps;
    return {taps, bands, Symmetry::symmetric, Response::constant, spec.density};
}

[[nodiscard]] Specification resolve(Specification spec, std::size_t max_points) noexcept {
    if (spec.symmetry != Symmetry::automatic) return spec;
    if (spec.response == Response::differentiator) {
        spec.symmetry = Symmetry::antisymmetric;
        return spec;
    }
    spec.symmetry = Symmetry::symmetric;
    if (layout(spec, max_points).status != Status::incompatible_endpoint) return spec;
    auto alternative = spec;
    alternative.symmetry = Symmetry::antisymmetric;
    return layout(alternative, max_points).status == Status::ok ? alternative : spec;
}

}

Exchange::Exchange(const Lowpass& spec, std::size_t max_points)
    : Exchange(lowpass_specification(spec, lowpass_bands(spec)), max_points) {}

Exchange::Exchange(const Specification& requested, std::size_t max_points) {
    const auto spec = resolve(requested, max_points);
    symmetry_ = spec.symmetry;
    response_ = spec.response;
    const auto size = layout(spec, max_points);
    status_ = size.status;
    if (status_ != Status::ok) return;
    type_ = form(spec.symmetry, spec.taps);
    const auto n = degree_of(spec) + 2;
    grid_.resize(size.points);
    reference_.resize(n);
    bands_.resize(spec.bands.size());
    nodes_.resize(n);
    peaks_.resize(size.points);
    samples_.resize(spec.taps / 2 + 1);
    taps_.resize(spec.taps);
    const double step = grid_step(spec);
    std::size_t count = 0;
    for (std::size_t b = 0; b < spec.bands.size(); ++b) {
        const auto& band = spec.bands[b];
        const auto segments = intervals(band, step, max_points);
        bands_[b] = {band, count, count};
        for (std::size_t j = 0; j <= segments; ++j) {
            const double omega = std::lerp(band.lower, band.upper, static_cast<double>(j) / static_cast<double>(segments));
            if (omitted(type_, response_, band, omega)) continue;
            const auto p = point(band, omega);
            const bool finite = std::isfinite(p.desired) && std::isfinite(p.weight) && p.weight > 0.0;
            const bool descending = count == 0 || p.x < grid_[count - 1].x;
            if (!finite || !descending) {
                status_ = Status::numerical_failure;
                return;
            }
            grid_[count++] = p;
        }
        bands_[b].end = count;
    }
}

Exchange::Exchange(Exchange&& other) noexcept
    : grid_(std::move(other.grid_)),
      reference_(std::move(other.reference_)),
      bands_(std::move(other.bands_)),
      nodes_(std::move(other.nodes_)),
      peaks_(std::move(other.peaks_)),
      samples_(std::move(other.samples_)),
      taps_(std::move(other.taps_)),
      symmetry_(other.symmetry_),
      response_(other.response_),
      type_(other.type_),
      status_(std::exchange(other.status_, Status::invalid_specification)),
      ready_(std::exchange(other.ready_, false)),
      exact_(std::exchange(other.exact_, false)),
      constant_(other.constant_) {}

Exchange& Exchange::operator=(Exchange&& other) noexcept {
    grid_ = std::move(other.grid_);
    reference_ = std::move(other.reference_);
    bands_ = std::move(other.bands_);
    nodes_ = std::move(other.nodes_);
    peaks_ = std::move(other.peaks_);
    samples_ = std::move(other.samples_);
    taps_ = std::move(other.taps_);
    symmetry_ = other.symmetry_;
    response_ = other.response_;
    type_ = other.type_;
    constant_ = other.constant_;
    status_ = std::exchange(other.status_, Status::invalid_specification);
    ready_ = std::exchange(other.ready_, false);
    exact_ = std::exchange(other.exact_, false);
    return *this;
}

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

double Exchange::interpolate(double x, Values values) const noexcept {
    double numerator = 0.0;
    double denominator = 0.0;
    for (const auto& node : std::span{nodes_}.first(nodes_.size() - 1)) {
        const double value = values == Values::fitted ? node.y : node.residual;
        if (x == node.x) return value;
        const double t = node.w / (x - node.x);
        numerator += t * value;
        denominator += t;
    }
    return numerator / denominator;
}

double Exchange::error(const Point& p) const noexcept {
    double numerator = 0.0;
    double denominator = 0.0;
    for (const auto& node : std::span{nodes_}.first(nodes_.size() - 1)) {
        const double e = (p.desired - node.desired) - node.shift;
        if (p.x == node.x) return p.weight * e;
        const double t = node.w / (p.x - node.x);
        numerator += t * e;
        denominator += t;
    }
    return p.weight * numerator / denominator;
}

Status Exchange::fit(double& delta) noexcept {
    const auto n = reference_.size();
    const auto difference = [&](std::size_t k, std::size_t j) { return 2 * (reference_[k].x - reference_[j].x); };
    int largest = std::numeric_limits<int>::min();
    for (std::size_t k = 0; k < n; ++k) {
        double product = 1.0;
        int exponent = 0;
        int e = 0;
        for (std::size_t j = 0, factors = 0; j < n; ++j) {
            if (j == k) continue;
            product *= difference(k, j);
            if (++factors % 16 == 0) {
                product = std::frexp(product, &e);
                exponent += e;
            }
        }
        product = std::frexp(product, &e);
        exponent += e;
        if (product == 0.0 || !std::isfinite(product)) return Status::numerical_failure;
        nodes_[k] = {reference_[k].omega, reference_[k].x, 1.0 / product, 0.0, 0.0, reference_[k].desired, 0.0};
        nodes_[k].y = -exponent;
        largest = std::max(largest, -exponent);
    }
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        auto& node = nodes_[k];
        node.w = std::ldexp(node.w, static_cast<int>(node.y) - largest);
        if (node.w == 0.0) return Status::numerical_failure;
        numerator += node.w * reference_[k].desired;
        denominator += (k % 2 == 0 ? node.w : -node.w) / reference_[k].weight;
    }
    delta = numerator / denominator;
    if (!std::isfinite(delta)) return Status::numerical_failure;
    for (std::size_t k = 0; k < n; ++k) {
        nodes_[k].shift = -(k % 2 == 0 ? delta : -delta) / reference_[k].weight;
        nodes_[k].y = nodes_[k].desired + nodes_[k].shift;
    }
    const auto dropped = std::ranges::max_element(nodes_, {}, [](const Node& node) { return std::abs(node.w); });
    const auto m = static_cast<std::size_t>(dropped - nodes_.begin());
    for (std::size_t k = 0; k < n; ++k)
        if (k != m) nodes_[k].w *= difference(k, m);
    std::swap(nodes_[m], nodes_.back());
    return Status::ok;
}

void Exchange::refine(Peak& peak) const noexcept {
    const auto j = peak.index;
    const auto& r = range(j);
    if (j == r.begin || j + 1 == r.end) return;
    const double sign = peak.error > 0 ? 1.0 : -1.0;
    double a = grid_[j - 1].omega;
    double b = peak.omega;
    double c = grid_[j + 1].omega;
    double fa = sign * error(grid_[j - 1]);
    double fb = sign * peak.error;
    double fc = sign * error(grid_[j + 1]);
    const double resolution = 1e-6 * (c - a);
    for (int iteration = 0; iteration < 16; ++iteration) {
        const double p = (b - a) * (fb - fc);
        const double q = (b - c) * (fb - fa);
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

void Exchange::select(std::span<Peak>& peaks) const noexcept {
    const auto n = reference_.size();
    const auto magnitude = [](const Peak& p) { return std::abs(p.error); };
    while (peaks.size() > n) {
        std::size_t k = 0;
        if (peaks.size() == n + 1)
            k = magnitude(peaks.front()) < magnitude(peaks.back()) ? 0 : n;
        else
            k = static_cast<std::size_t>(std::ranges::min_element(peaks, {}, magnitude) - peaks.begin());
        const bool interior = k != 0 && k + 1 != peaks.size();
        const std::size_t removed = interior ? 2 : 1;
        if (interior && magnitude(peaks[k + 1]) > magnitude(peaks[k - 1])) peaks[k - 1] = peaks[k + 1];
        std::shift_left(peaks.begin() + static_cast<std::ptrdiff_t>(k), peaks.end(), static_cast<std::ptrdiff_t>(removed));
        peaks = peaks.first(peaks.size() - removed);
    }
}

Result Exchange::run(std::size_t max_iterations, double tolerance) noexcept {
    ready_ = false;
    exact_ = false;
    if (status_ != Status::ok) return {status_};
    if (!(tolerance > 0.0 && tolerance < 1.0)) return {Status::invalid_tolerance};
    if (max_iterations == 0) return {Status::nonconvergence};
    constant_ = bands_.front().band.target;
    const bool zero = std::ranges::all_of(bands_, [](const Range& r) { return r.band.target == 0.0; });
    const bool flat = std::ranges::all_of(bands_, [&](const Range& r) { return r.band.target == constant_; });
    if (zero || (type_ == 1 && response_ == Response::constant && flat)) {
        std::ranges::fill(taps_, 0.0);
        if (type_ == 1) taps_[taps_.size() / 2] = constant_;
        ready_ = true;
        exact_ = true;
        return {Status::ok};
    }
    seed();
    const auto n = reference_.size();
    const auto magnitude = [](const Peak& p) { return std::abs(p.error); };
    Result result{Status::nonconvergence};
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        result.iterations = iteration + 1;
        double delta = 0.0;
        if (const auto status = fit(delta); status != Status::ok) {
            result.status = status;
            return result;
        }
        result.ripple = std::abs(delta);
        std::size_t count = 0;
        for (std::size_t j = 0; j < grid_.size(); ++j) {
            const double e = error(grid_[j]);
            if (!std::isfinite(e)) {
                result.status = Status::numerical_failure;
                return result;
            }
            if (e == 0.0) continue;
            if (count != 0 && std::signbit(e) == std::signbit(peaks_[count - 1].error)) {
                if (std::abs(e) > std::abs(peaks_[count - 1].error)) peaks_[count - 1] = {j, grid_[j].omega, e};
            } else {
                peaks_[count++] = {j, grid_[j].omega, e};
            }
        }
        if (count < n) {
            result.status = Status::numerical_failure;
            return result;
        }
        auto peaks = std::span{peaks_}.first(count);
        for (auto& peak : peaks) refine(peak);
        result.max_error = magnitude(std::ranges::max(peaks, {}, magnitude));
        select(peaks);
        const auto crossed = [](const Peak& left, const Peak& right) { return left.omega >= right.omega; };
        if (std::ranges::adjacent_find(peaks, crossed) != peaks.end()) {
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

void Exchange::inverse(Values values, bool accumulate) noexcept {
    const auto N = taps_.size();
    const bool symmetric = symmetry_ == Symmetry::symmetric;
    for (std::size_t m = 0; m < samples_.size(); ++m) {
        const double omega = 2 * pi * static_cast<double>(m) / static_cast<double>(N);
        samples_[m] = forced(type_, omega) ? 0.0 : factor(type_, omega) * interpolate(std::cos(omega), values);
    }
    const auto period = static_cast<std::ptrdiff_t>(2 * N);
    for (std::size_t n = 0; n < (N + 1) / 2; ++n) {
        const auto t = static_cast<std::ptrdiff_t>(2 * n + 1) - static_cast<std::ptrdiff_t>(N);
        double sum = 0.0;
        for (std::size_t m = 0; m < samples_.size(); ++m) {
            const auto q = (static_cast<std::ptrdiff_t>(m) * t % period + period) % period;
            const double angle = pi * static_cast<double>(q) / static_cast<double>(N);
            const double multiplicity = m == 0 || 2 * m == N ? 1.0 : 2.0;
            sum += multiplicity * samples_[m] * (symmetric ? std::cos(angle) : -std::sin(angle));
        }
        const double h = t == 0 && !symmetric ? 0.0 : sum / static_cast<double>(N);
        taps_[n] = accumulate ? taps_[n] + h : h;
        taps_[N - 1 - n] = symmetric ? taps_[n] : -taps_[n];
    }
}

double Exchange::tap_polynomial(double omega) const noexcept {
    const auto N = taps_.size();
    const bool symmetric = symmetry_ == Symmetry::symmetric;
    const bool slope = omega == 0.0 && !symmetric;
    double a = 0.0;
    for (std::size_t k = 0; k < N; ++k) {
        const double t = (static_cast<double>(2 * k + 1) - static_cast<double>(N)) / 2;
        if (slope)
            a -= taps_[k] * t;
        else
            a += taps_[k] * (symmetric ? std::cos(omega * t) : -std::sin(omega * t));
    }
    return slope ? a / factor_over_omega(type_, 0.0) : a / factor(type_, omega);
}

Status Exchange::recover() noexcept {
    inverse(Values::fitted, false);
    for (auto& node : std::span{nodes_}.first(nodes_.size() - 1)) node.residual = node.y - tap_polynomial(node.omega);
    inverse(Values::residual, true);
    const bool finite = std::ranges::all_of(taps_, [](double h) { return std::isfinite(h); });
    return finite ? Status::ok : Status::numerical_failure;
}

std::span<const double> Exchange::taps() const noexcept {
    if (!ready_) return {};
    return taps_;
}

std::optional<double> Exchange::amplitude(double omega) const noexcept {
    if (!ready_ || !(omega >= 0.0 && omega <= pi)) return std::nullopt;
    if (exact_) return type_ == 1 ? constant_ : 0.0;
    if (forced(type_, omega)) return 0.0;
    return factor(type_, omega) * interpolate(std::cos(omega));
}

std::span<const Point> Exchange::extrema() const noexcept {
    if (!ready_ || exact_) return {};
    return reference_;
}

std::span<const Point> Exchange::grid() const noexcept {
    if (status_ != Status::ok) return {};
    return grid_;
}

}
