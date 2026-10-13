#include <remez.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace {

constexpr double sample_rate = 48000.0;

void report(std::string_view stage, remez::Status status) {
    const auto name = remez::to_string(status);
    std::fprintf(stderr, "%.*s: %.*s\n", static_cast<int>(stage.size()), stage.data(), static_cast<int>(name.size()), name.data());
}

double tone_gain(std::span<const double> taps, double hz) {
    const double omega = 2 * std::numbers::pi_v<double> * hz / sample_rate;
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
    using remez::Region;
    const std::array spec{
        remez::DecibelBand{.region = Region::stop, .lower_hz = 0.0, .upper_hz = 6000.0, .decibels = 60.0},
        remez::DecibelBand{.region = Region::pass, .lower_hz = 8000.0, .upper_hz = 12000.0, .decibels = 0.5},
        remez::DecibelBand{.region = Region::stop, .lower_hz = 14000.0, .upper_hz = 24000.0, .decibels = 60.0},
    };
    std::array<remez::Band, spec.size()> bands{};
    if (const auto status = remez::convert(spec, sample_rate, bands); status != remez::Status::ok) {
        report("convert", status);
        return 1;
    }

    std::optional<std::vector<double>> taps;
    remez::Result result{remez::Status::nonconvergence};
    for (std::size_t length = 9; length <= 255 && !taps; length += 2) {
        remez::Exchange designer(remez::Specification{.taps = length, .bands = bands});
        result = designer.run();
        if (result.status != remez::Status::ok) {
            report("design", result.status);
            return 1;
        }
        if (result.max_error <= 1.0) taps.emplace(designer.taps().begin(), designer.taps().end());
    }
    if (!taps) {
        std::fprintf(stderr, "no length up to 255 taps meets the specification\n");
        return 1;
    }

    std::printf("%zu taps meet the specification after %zu iterations\n", taps->size(), result.iterations);
    std::printf("passband ripple %.3f dB (asked 0.5), stopband %.1f dB (asked 60)\n",
                remez::achieved_decibels(spec[1], result.max_error), remez::achieved_decibels(spec[0], result.max_error));
    for (const double hz : {3000.0, 10000.0, 18000.0})
        std::printf("%6.0f Hz tone: gain %.6f (%.1f dB)\n", hz, tone_gain(*taps, hz), 20 * std::log10(tone_gain(*taps, hz)));
}
