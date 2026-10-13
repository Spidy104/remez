#include <remez.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <span>
#include <string_view>

namespace {

constexpr double pi = std::numbers::pi_v<double>;

void design(std::string_view name, std::size_t taps, std::span<const remez::Band> bands) {
    remez::Exchange designer(remez::Specification{.taps = taps, .bands = bands, .symmetry = remez::Symmetry::automatic});
    const auto result = designer.run();
    const auto status = remez::to_string(result.status);
    const auto shape = designer.symmetry() == remez::Symmetry::antisymmetric ? "antisymmetric" : "symmetric";
    std::printf("%-9.*s %2zu taps  %-22.*s", static_cast<int>(name.size()), name.data(), taps, static_cast<int>(status.size()),
                status.data());
    if (result.status == remez::Status::ok)
        std::printf("%-13s ripple %.1f dB", shape, 20 * std::log10(result.ripple));
    std::printf("\n");
}

}

int main() {
    const std::array highpass{remez::Band{0, 0.4 * pi, 0, 1}, remez::Band{0.55 * pi, pi, 1, 1}};
    const std::array bandstop{remez::Band{0, 0.2 * pi, 1, 1}, remez::Band{0.35 * pi, 0.55 * pi, 0, 1}, remez::Band{0.7 * pi, pi, 1, 1}};
    design("highpass", 8, highpass);
    design("highpass", 9, highpass);
    design("bandstop", 8, bandstop);
    design("bandstop", 9, bandstop);
}
