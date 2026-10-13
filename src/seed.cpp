#include <remez.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>

namespace remez {

bool Exchange::equilibrium() noexcept {
    constexpr std::size_t nodes = 64;
    constexpr double settled = 0x1p-40;
    const auto k = bands_.size();
    const auto n = reference_.size();
    const bool too_narrow = std::ranges::any_of(bands_, [](const Range& b) { return b.end - b.begin < 2; });
    if (n < k || samples_.size() + 1 < k || too_narrow) return false;
    const auto roots = std::span{samples_}.first(k - 1);
    const auto endpoint = [&](std::size_t j) { return grid_[j % 2 == 0 ? bands_[j / 2].begin : bands_[j / 2].end - 1].x; };
    const auto weight = [&](double x, std::size_t interval, std::size_t skipped_root) {
        double value = 1.0;
        for (std::size_t g = 0; g < roots.size(); ++g)
            if (g != skipped_root) value *= x - roots[g];
        for (std::size_t j = 0; j < 2 * k; ++j)
            if (j != interval && j != interval + 1) value /= std::sqrt(std::abs(x - endpoint(j)));
        return std::abs(value);
    };
    const auto node = [&](std::size_t interval, std::size_t l) {
        const double t = std::numbers::pi_v<double> * (static_cast<double>(l) + 0.5) / static_cast<double>(nodes);
        const double half = (endpoint(interval) - endpoint(interval + 1)) / 2;
        return std::midpoint(endpoint(interval), endpoint(interval + 1)) + half * std::cos(t);
    };
    for (std::size_t g = 0; g < roots.size(); ++g) roots[g] = std::midpoint(endpoint(2 * g + 1), endpoint(2 * g + 2));
    for (int sweep = 0; sweep < 64; ++sweep) {
        double change = 0.0;
        for (std::size_t g = 0; g < roots.size(); ++g) {
            double moment = 0.0;
            double mass = 0.0;
            for (std::size_t l = 0; l < nodes; ++l) {
                const double x = node(2 * g + 1, l);
                const double w = weight(x, 2 * g + 1, g);
                moment += x * w;
                mass += w;
            }
            const double root = moment / mass;
            if (!(mass > 0.0 && std::isfinite(root))) return false;
            change = std::max(change, std::abs(root - roots[g]));
            roots[g] = root;
        }
        if (change <= settled) break;
    }
    const auto band_mass = [&](std::size_t b) {
        double mass = 0.0;
        for (std::size_t l = 0; l < nodes; ++l) mass += weight(node(2 * b, l), 2 * b, k);
        return mass;
    };
    double total = 0.0;
    for (std::size_t b = 0; b < k; ++b) total += band_mass(b);
    if (!(total > 0.0 && std::isfinite(total))) return false;
    const auto spare = n - k;
    std::size_t out = 0;
    std::size_t used = 0;
    double before = 0.0;
    for (std::size_t b = 0; b < k; ++b) {
        before += band_mass(b);
        const bool last_band = b + 1 == k;
        const auto proportional = static_cast<std::size_t>(static_cast<double>(spare) * before / total);
        const auto share = last_band ? spare : std::min(spare, proportional);
        const auto m = 1 + share - used;
        used = share;
        const auto band = std::span{grid_}.subspan(bands_[b].begin, bands_[b].end - bands_[b].begin);
        if (m > band.size()) return false;
        const double centre = std::midpoint(band.front().x, band.back().x);
        const double half = (band.front().x - band.back().x) / 2;
        const auto density = [&](std::size_t j) { return weight(band[j].x, 2 * b, k); };
        const auto angle = [&](std::size_t j) { return std::acos(std::clamp((band[j].x - centre) / half, -1.0, 1.0)); };
        const auto slice = [&](std::size_t j) { return (density(j - 1) + density(j)) / 2 * (angle(j) - angle(j - 1)); };
        double area = 0.0;
        for (std::size_t j = 1; j < band.size(); ++j) area += slice(j);
        if (!(area > 0.0 && std::isfinite(area))) return false;
        const auto steps = last_band ? m : m - 1;
        const auto target = [&](std::size_t j) {
            return steps == 0 ? area / 2 : area * static_cast<double>(j) / static_cast<double>(steps);
        };
        std::size_t placed = 0;
        std::size_t previous = 0;
        const auto place = [&](std::size_t j) {
            j = std::min(std::max(j, placed == 0 ? 0 : previous + 1), band.size() - (m - placed));
            reference_[out++] = band[j];
            previous = j;
            ++placed;
        };
        double cdf = 0.0;
        double last = 0.0;
        for (std::size_t j = 0; j < band.size() && placed < m; ++j) {
            if (j != 0) cdf += slice(j);
            while (placed < m && cdf >= target(placed))
                place(j != 0 && target(placed) - last < cdf - target(placed) ? j - 1 : j);
            last = cdf;
        }
        while (placed < m) place(band.size() - 1);
    }
    return true;
}

void Exchange::seed() noexcept {
    if (equilibrium()) return;
    const auto n = reference_.size();
    const auto step = (grid_.size() - 1) / (n - 1);
    const auto extra = (grid_.size() - 1) % (n - 1);
    for (std::size_t k = 0, index = 0, carry = 0; k < n; ++k) {
        reference_[k] = grid_[index];
        index += step;
        carry += extra;
        if (carry >= n - 1) {
            ++index;
            carry -= n - 1;
        }
    }
    const auto desired = [](const Point& p) { return p.desired; };
    const auto [low, high] = std::ranges::minmax(reference_, {}, desired);
    if (low.desired != high.desired) return;
    const auto [lowest, highest] = std::ranges::minmax(grid_, {}, desired);
    const bool lower_differs_more = std::abs(lowest.desired - low.desired) > std::abs(highest.desired - low.desired);
    reference_[n / 2] = lower_differs_more ? lowest : highest;
    std::ranges::sort(reference_, {}, [](const Point& p) { return p.omega; });
}

}
