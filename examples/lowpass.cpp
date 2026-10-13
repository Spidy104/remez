#include <remez.hpp>

#include <cstdio>
#include <numbers>

int main() {
    constexpr double pi = std::numbers::pi_v<double>;
    remez::Exchange designer(remez::Lowpass{.taps = 31, .pass_edge = 0.2 * pi, .stop_edge = 0.3 * pi, .weights = {1.0, 10.0}});
    const auto result = designer.run();
    if (result.status != remez::Status::ok) return 1;
    std::printf("31-tap lowpass, %zu iterations, passband ripple %.4f, stopband ripple %.5f\n", result.iterations, result.ripple,
                result.ripple / 10.0);
    for (const double h : designer.taps()) std::printf("% .12f\n", h);
}
