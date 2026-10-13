#include <remez.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

int main() {
    constexpr double pi = std::numbers::pi_v<double>;
    const std::array band{remez::Band{.lower = 0.0, .upper = 0.8 * pi, .target = 1.0}, remez::Band{.lower = 0.9 * pi, .upper = pi, .target = 0.0}};
    remez::Exchange designer(remez::Specification{
        .taps = 31, .bands = band, .symmetry = remez::Symmetry::automatic, .response = remez::Response::differentiator});
    const auto result = designer.run();
    if (result.status != remez::Status::ok) return 1;
    const auto h = designer.taps();
    const auto delay = static_cast<double>(h.size() - 1) / 2;

    std::printf("31-tap differentiator, relative error %.2e\n", result.ripple);
    for (const double omega : {0.1 * pi, 0.4 * pi, 0.7 * pi}) {
        std::vector<double> x(512);
        for (std::size_t n = 0; n < x.size(); ++n) x[n] = std::sin(omega * static_cast<double>(n));
        double worst = 0.0;
        for (std::size_t n = h.size(); n < x.size(); ++n) {
            double y = 0.0;
            for (std::size_t k = 0; k < h.size(); ++k) y += h[k] * x[n - k];
            const double derivative = omega * std::cos(omega * (static_cast<double>(n) - delay));
            worst = std::fmax(worst, std::fabs(y - derivative) / omega);
        }
        std::printf("sin(%.1f pi n): output matches the derivative to %.2e relative\n", omega / pi, worst);
    }
}
