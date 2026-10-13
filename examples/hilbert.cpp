#include <remez.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

int main() {
    constexpr double pi = std::numbers::pi_v<double>;
    const std::array band{remez::Band{.lower = 0.05 * pi, .upper = 0.95 * pi, .target = -1.0}};
    remez::Exchange designer(remez::Specification{.taps = 31, .bands = band, .symmetry = remez::Symmetry::antisymmetric});
    const auto result = designer.run();
    if (result.status != remez::Status::ok) return 1;
    const auto h = designer.taps();
    const std::size_t delay = (h.size() - 1) / 2;

    const double omega = 0.3 * pi;
    std::vector<double> x(512);
    for (std::size_t n = 0; n < x.size(); ++n) x[n] = std::cos(omega * static_cast<double>(n));
    double worst = 0.0;
    for (std::size_t n = h.size(); n < x.size(); ++n) {
        double y = 0.0;
        for (std::size_t k = 0; k < h.size(); ++k) y += h[k] * x[n - k];
        const double expected = std::sin(omega * static_cast<double>(n - delay));
        worst = std::fmax(worst, std::fabs(y - expected));
    }
    std::printf("31-tap Hilbert transformer, ripple %.2e\n", result.ripple);
    std::printf("cos(0.3 pi n) in, sin(0.3 pi (n - %zu)) out, worst deviation %.2e\n", delay, worst);
}
