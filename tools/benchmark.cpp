#include <remez.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <string_view>

namespace {

constexpr double pi = std::numbers::pi_v<double>;

struct Case {
    std::string_view name;
    remez::Specification spec;
};

double median_ms(remez::Exchange& designer) {
    std::array<double, 5> times{};
    for (auto& time : times) {
        const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < 10; ++repeat)
            if (designer.run().status != remez::Status::ok) std::exit(1);
        time = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 10;
    }
    std::ranges::sort(times);
    return times[2];
}

}

int main() {
    const std::array lowpass{remez::Band{0, .2 * pi, 1, 1}, remez::Band{.3 * pi, pi, 0, 10}};
    const std::array narrow{remez::Band{0, .24 * pi, 1, 1}, remez::Band{.26 * pi, pi, 0, 10}};
    const std::array bandpass{remez::Band{0, .2 * pi, 0, 10}, remez::Band{.3 * pi, .6 * pi, 1, 1}, remez::Band{.7 * pi, pi, 0, 10}};
    const std::array hilbert{remez::Band{.05 * pi, .95 * pi, -1, 1}};
    const std::array slope{remez::Band{0, .8 * pi, 1, 1}, remez::Band{.9 * pi, pi, 0, 10}};
    const auto anti = remez::Symmetry::antisymmetric;
    const std::array cases{
        Case{"lowpass", {15, lowpass}},
        Case{"lowpass", {63, lowpass}},
        Case{"lowpass", {127, lowpass}},
        Case{"lowpass narrow", {255, narrow}},
        Case{"lowpass", {223, lowpass}},
        Case{"lowpass", {255, lowpass}},
        Case{"bandpass", {223, bandpass}},
        Case{"hilbert III", {223, hilbert, anti}},
        Case{"differentiator III", {223, slope, anti, remez::Response::differentiator}},
    };
    for (const auto& c : cases) {
        remez::Exchange designer(c.spec);
        const auto result = designer.run();
        if (result.status != remez::Status::ok) return 1;
        std::printf("%-20.*s taps=%-4zu iterations=%-3zu ripple=%.3e flatness=%.1e ms=%.3f\n", static_cast<int>(c.name.size()),
                    c.name.data(), c.spec.taps, result.iterations, result.ripple, result.flatness, median_ms(designer));
    }
}
