#include <remez.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <span>
#include <string_view>
#include <vector>

namespace {

constexpr double pi = std::numbers::pi_v<double>;

int fail(std::string_view stage, remez::Status status) {
    const auto name = remez::to_string(status);
    std::fprintf(stderr, "%.*s failed: %.*s\n", static_cast<int>(stage.size()), stage.data(),
                 static_cast<int>(name.size()), name.data());
    return 1;
}

double tone_gain(std::span<const double> taps, double omega) {
    std::vector<double> input(4096);
    for (std::size_t n = 0; n < input.size(); ++n) input[n] = std::sin(omega * static_cast<double>(n));
    double energy = 0.0;
    for (std::size_t n = taps.size(); n < input.size(); ++n) {
        double output = 0.0;
        for (std::size_t k = 0; k < taps.size(); ++k) output += taps[k] * input[n - k];
        energy += output * output;
    }
    return std::sqrt(2 * energy / static_cast<double>(input.size() - taps.size()));
}

}

int main() {
    const std::array bands{
        remez::Band{.lower = 0.0, .upper = 0.2 * pi, .target = 0.0, .weight = 10.0},
        remez::Band{.lower = 0.3 * pi, .upper = 0.5 * pi, .target = 1.0, .weight = 1.0},
        remez::Band{.lower = 0.6 * pi, .upper = pi, .target = 0.0, .weight = 10.0},
    };
    remez::Exchange designer(remez::Specification{.taps = 63, .bands = bands});
    if (designer.status() != remez::Status::ok) return fail("setup", designer.status());

    const remez::Result result = designer.run();
    if (result.status != remez::Status::ok) return fail("design", result.status);

    const std::vector<double> taps(designer.taps().begin(), designer.taps().end());
    const double stopband_db = 20 * std::log10(result.ripple / bands[0].weight);
    std::printf("taps %zu, iterations %zu, passband ripple %.2e, stopband %.1f dB, flatness %.1e\n", taps.size(),
                result.iterations, result.ripple, stopband_db, result.flatness);
    for (const double omega : {0.1 * pi, 0.4 * pi, 0.8 * pi})
        std::printf("tone at %.1f pi: gain %.6f\n", omega / pi, tone_gain(taps, omega));
}
