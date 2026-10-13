#include <remez.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <optional>
#include <source_location>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef REMEZ_CRT_CHECKS
#include <crtdbg.h>
#include <memory>
#include <new>
#endif

namespace {

using remez::Band;
using remez::Exchange;
using remez::Lowpass;
using remez::Response;
using remez::Specification;
using remez::Status;
using remez::Symmetry;

constexpr double pi = std::numbers::pi_v<double>;
constexpr double eps = std::numeric_limits<double>::epsilon();

static_assert(!std::is_copy_constructible_v<Exchange> && !std::is_copy_assignable_v<Exchange>);
static_assert(std::is_nothrow_move_constructible_v<Exchange> && std::is_nothrow_move_assignable_v<Exchange>);
static_assert(remez::to_string(Status::nonconvergence) == "nonconvergence");
static_assert(remez::to_string(Status::ok) == "ok");

void check(bool condition, std::string_view message, std::source_location at = std::source_location::current()) {
    if (condition) return;
    std::fprintf(stderr, "%s:%u: %.*s\n", at.file_name(), static_cast<unsigned>(at.line()), static_cast<int>(message.size()),
                 message.data());
    std::exit(1);
}

double amplitude(const Exchange& engine, double omega, std::source_location at = std::source_location::current()) {
    const auto value = engine.amplitude(omega);
    check(value.has_value(), "amplitude available", at);
    return value.value_or(0.0);
}

std::complex<double> shifted_response(std::span<const double> h, double omega) {
    std::complex<double> sum{};
    for (std::size_t k = 0; k < h.size(); ++k) sum += h[k] * std::polar(1.0, -omega * static_cast<double>(k));
    return sum * std::polar(1.0, omega * static_cast<double>(h.size() - 1) / 2);
}

double weighted_error(const Band& band, Response response, double omega, double a, double slope) {
    if (response == Response::differentiator && band.target != 0.0)
        return band.weight * (band.target - (omega == 0.0 ? slope : a / omega)) / std::abs(band.target);
    return band.weight * (band.target - a);
}

double factor(const Specification& spec, double omega) {
    if (spec.symmetry == Symmetry::symmetric) return spec.taps % 2 == 0 ? std::cos(omega / 2) : 1.0;
    return spec.taps % 2 == 0 ? std::sin(omega / 2) : std::sin(omega);
}

class EvaluationBound {
public:
    EvaluationBound(std::span<const remez::Point> reference, double ripple) {
        const auto n = reference.size();
        std::vector<double> log_weight(n);
        for (std::size_t k = 0; k < n; ++k) {
            x_.push_back(reference[k].x);
            y_.push_back(std::abs(reference[k].desired) + ripple / reference[k].weight);
            double sign = 1.0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j == k) continue;
                const double a = reference[k].omega;
                const double b = reference[j].omega;
                const double d = -4 * std::sin((a + b) / 2) * std::sin((a - b) / 2);
                log_weight[k] -= std::log(std::abs(d));
                if (d < 0) sign = -sign;
            }
            w_.push_back(sign);
        }
        const double top = std::ranges::max(log_weight);
        for (std::size_t k = 0; k < n; ++k) w_[k] *= std::exp(log_weight[k] - top);
    }

    [[nodiscard]] double operator()(double omega, double q, double value) const {
        const double at = std::cos(omega);
        double sum = 0.0;
        for (std::size_t k = 0; k < x_.size(); ++k) {
            if (at == x_[k]) return 0.0;
            sum += w_[k] / (at - x_[k]);
        }
        double lebesgue = 0.0;
        double weighted = 0.0;
        for (std::size_t k = 0; k < x_.size(); ++k) {
            const double l = std::abs(w_[k] / (at - x_[k]) / sum);
            lebesgue += l;
            weighted += l * y_[k];
        }
        const auto n = static_cast<double>(x_.size());
        return eps * ((3 * n + 4) * weighted * std::abs(q) + (3 * n + 2) * lebesgue * std::abs(value));
    }

private:
    std::vector<double> x_;
    std::vector<double> w_;
    std::vector<double> y_;
};

remez::Result verify(const Specification& spec, double tolerance = 1e-6) {
    Exchange engine(spec);
    check(engine.status() == Status::ok && engine.taps().empty() && !engine.amplitude(0.5), "nothing before run");
    const auto result = engine.run(128, tolerance);
    if (result.status != Status::ok) {
        const auto name = remez::to_string(result.status);
        std::fprintf(stderr, "taps=%zu status=%.*s\n", spec.taps, static_cast<int>(name.size()), name.data());
    }
    check(result.status == Status::ok && result.flatness <= tolerance, "exchange converges");
    const auto h = engine.taps();
    const bool anti = spec.symmetry == Symmetry::antisymmetric;
    check(h.size() == spec.taps, "tap count");
    for (std::size_t k = 0; k < h.size(); ++k)
        check(std::isfinite(h[k]) && h[k] == (anti ? -h[h.size() - 1 - k] : h[h.size() - 1 - k]), "finite taps with exact symmetry");
    double scale = 1.0;
    for (const auto& band : spec.bands)
        scale = std::max(scale, std::abs(band.target) * (spec.response == Response::differentiator ? pi : 1.0));
    const double tap_tolerance = 64 * static_cast<double>(spec.taps) * eps * scale;
    double slope = 0.0;
    for (std::size_t k = 0; k < h.size(); ++k) slope -= static_cast<double>(k) * h[k];
    const EvaluationBound evaluation(engine.extrema(), result.ripple);
    double excess = 0.0;
    for (const auto& band : spec.bands) {
        for (std::size_t j = 0; j <= 2048; ++j) {
            const double omega = std::lerp(band.lower, band.upper, static_cast<double>(j) / 2048);
            const auto H = shifted_response(h, omega);
            const double a = anti ? H.imag() : H.real();
            const double internal = amplitude(engine, omega);
            check(std::abs(a - internal) <= tap_tolerance + evaluation(omega, factor(spec, omega), internal) &&
                      std::abs(anti ? H.real() : H.imag()) <= tap_tolerance,
                  "complex tap response matches amplitude, linear phase");
            const bool relative = spec.response == Response::differentiator && band.target != 0.0;
            const double per_omega = 1.0 / std::max(omega, 1.0 / static_cast<double>(spec.taps)) / std::abs(band.target);
            const double allowance = band.weight * tap_tolerance * (relative ? per_omega : 1.0);
            excess = std::max(excess, std::abs(weighted_error(band, spec.response, omega, a, slope)) - allowance);
        }
    }
    check(excess <= result.max_error * (1 + 3 * tolerance) + 1e-12, "no missed peak in dense physical search");
    const auto extrema = engine.extrema();
    check(extrema.size() == (anti ? spec.taps / 2 + 1 : (spec.taps - 1) / 2 + 2), "reference has degree+2 points");
    double previous = 0.0;
    for (const auto& p : extrema) {
        const auto band = std::ranges::find_if(spec.bands, [&](const Band& b) { return p.omega >= b.lower && p.omega <= b.upper; });
        check(band != spec.bands.end(), "extremum inside an active band");
        const double e = weighted_error(*band, spec.response, p.omega, amplitude(engine, p.omega), slope);
        check(std::abs(std::abs(e) - result.ripple) <= 1e-6 * result.ripple + 1e-12, "equal ripple at the reference");
        check(previous == 0.0 || std::signbit(previous) != std::signbit(e), "alternating signs");
        previous = e;
    }
    const std::vector<double> saved(h.begin(), h.end());
    const auto again = engine.run(128, tolerance);
    check(again.iterations == result.iterations && std::ranges::equal(engine.taps(), saved), "deterministic restart");
    check(engine.run(0).status == Status::nonconvergence && engine.taps().empty() && engine.extrema().empty() &&
              !engine.amplitude(0.5),
          "failed run withdraws every output");
    return result;
}

void lowpass() {
    const std::array edges{std::pair{0.2, 0.3}, std::pair{0.24, 0.26}, std::pair{0.01, 0.03}, std::pair{0.97, 0.99}};
    for (const std::size_t taps : {1u, 3u, 5u, 15u, 31u, 63u, 65u, 127u, 129u, 255u}) {
        for (const auto& [pass, stop] : edges) {
            const Lowpass spec{taps, pass * pi, stop * pi, {1, 10}, 32};
            const std::array bands{Band{0, spec.pass_edge, 1, 1}, Band{spec.stop_edge, pi, 0, 10}};
            const auto result = verify({taps, bands, Symmetry::symmetric, Response::constant, 32});
            Exchange adapter(spec);
            check(adapter.run().ripple == result.ripple, "lowpass adapter equals the general specification");
            if (taps == 1) check(std::abs(result.ripple - 10.0 / 11) <= 1e-15, "one tap: analytic weighted constant");
        }
    }
}

void scipy_references() {
    struct Reference {
        std::size_t taps;
        double pass;
        double stop;
        std::array<double, 6> values;
    };
    const std::array references{
        Reference{15, 0.2, 0.3, {1.369204595481571, 1.2240989852775044, 0.63079540451842875, 0.036920459548157314, -0.0067866262183710335, 0.03692045954815687}},
        Reference{63, 0.24, 0.26, {1.3441013601771825, 0.85340865793979082, 0.65589863982281948, 0.034410136017717363, -0.016367052366606641, 0.034410136017717807}},
        Reference{255, 0.24, 0.26, {1.0107264164727199, 0.98940766692846405, 0.98927358352728401, 0.0010726416472719436, -0.0007662460708280916, 0.001072641647274164}}};
    constexpr std::array taps15{-0.028458770784620859, -0.015306822035826059, 0.0062897668149520872, 0.05226038713404093,
                                0.11707879787873703,   0.18565661178492684,   0.23816124007428527,    0.25784217374858048};
    for (const auto& ref : references) {
        Exchange engine(Lowpass{ref.taps, ref.pass * pi, ref.stop * pi, {1, 10}, 16});
        const auto result = engine.run();
        check(result.status == Status::ok, "reference design converges");
        const std::array at{0.0, ref.pass / 2, ref.pass, ref.stop, (ref.stop + 1) / 2, 1.0};
        for (std::size_t k = 0; k < at.size(); ++k) {
            const double weight = k < 3 ? 1.0 : 10.0;
            check(weight * std::abs(amplitude(engine, at[k] * pi) - ref.values[k]) <= 3e-5 * result.ripple + 1e-12,
                  "SciPy minimax amplitude");
        }
        if (ref.taps == 15)
            for (std::size_t k = 0; k < taps15.size(); ++k)
                check(std::abs(engine.taps()[k] - taps15[k]) <= 2e-6 * result.ripple, "SciPy coefficients");
    }
}

void general() {
    const std::array lp{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, pi, 0, 10}};
    const std::array hp{Band{0, .2 * pi, 0, 10}, Band{.3 * pi, pi, 1, 1}};
    const std::array bp{Band{0, .2 * pi, 0, 10}, Band{.3 * pi, .6 * pi, 1, 1}, Band{.7 * pi, pi, 0, 10}};
    const std::array bs{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, .6 * pi, 0, 10}, Band{.7 * pi, pi, 1, 1}};
    const std::array multi{Band{0, .1 * pi, .2, 1}, Band{.2 * pi, .3 * pi, 1, 2}, Band{.4 * pi, .5 * pi, -.5, 3},
                           Band{.6 * pi, .7 * pi, .8, 4}, Band{.8 * pi, pi, 0, 10}};
    const std::array hilbert{Band{.05 * pi, .95 * pi, -1, 1}};
    const std::array diff3{Band{0, .8 * pi, 1, 1}, Band{.9 * pi, pi, 0, 10}};
    const std::array diff4{Band{0, pi, 1, 1}};
    constexpr std::array<double, 31> f_hp{-0.020934364437286704, 0.022393799499398452, 0.016594211795943312, 0.0069770856260239307, -0.0079594698470608159, -0.021778493670776066, -0.0241603603971648, -0.0086784968279752829, 0.020273503272464462, 0.046359155709411584, 0.048005945387335734, 0.010186383897047298, -0.06518463055894691, -0.15676028186370189, -0.23169478062173821, 0.73941489514358882, -0.23169478062173821, -0.15676028186370189, -0.06518463055894691, 0.010186383897047298, 0.048005945387335734, 0.046359155709411584, 0.020273503272464462, -0.0086784968279752829, -0.0241603603971648, -0.021778493670776066, -0.0079594698470608159, 0.0069770856260239307, 0.016594211795943312, 0.022393799499398452, -0.020934364437286704};
    constexpr std::array<double, 31> f_bp{-0.0085573186696835373, 0.0075097831224829516, 0.0183875770295574, 0.00072008855702256915, 0.0024710132944305596, -0.0041691981718649188, -0.053462928553118648, -0.018966368705298265, 0.065850503634368493, 0.022913613267214243, 0.008959854358972133, 0.084518681791226125, -0.087986916240881799, -0.27766409600142217, 0.054338215146355291, 0.37910873669969486, 0.054338215146355291, -0.27766409600142217, -0.087986916240881799, 0.084518681791226125, 0.008959854358972133, 0.022913613267214243, 0.065850503634368493, -0.018966368705298265, -0.053462928553118648, -0.0041691981718649188, 0.0024710132944305596, 0.00072008855702256915, 0.0183875770295574, 0.0075097831224829516, -0.0085573186696835373};
    constexpr std::array<double, 31> f_bs{-0.015735819024373918, 0.014835479339145067, -0.015260471671871252, 0.019044936255938588, 0.039794869035146015, 0.0035974877830199621, 0.019290082378912575, 0.019885688352301756, -0.067293868436692147, -0.041499868787664602, 0.014698710399737186, -0.057404169672093788, 0.084683458717752055, 0.29014327849504723, -0.060176961398609868, 0.57860163881005544, -0.060176961398609868, 0.29014327849504723, 0.084683458717752055, -0.057404169672093788, 0.014698710399737186, -0.041499868787664602, -0.067293868436692147, 0.019885688352301756, 0.019290082378912575, 0.0035974877830199621, 0.039794869035146015, 0.019044936255938588, -0.015260471671871252, 0.014835479339145067, -0.015735819024373918};
    constexpr std::array<double, 15> f_multi{0.085736008332173402, 0.13333084976016513, -0.15160227102350737, -0.16960606175362483, 0.1088321176387604, 0.055332501880287696, 0.1106009511449377, 0.31504872400767725, 0.1106009511449377, 0.055332501880287696, 0.1088321176387604, -0.16960606175362483, -0.15160227102350737, 0.13333084976016513, 0.085736008332173402};
    constexpr std::array<double, 32> f_lp2{-0.0064902982704557882, -0.0049652330188863069, -0.00088958916661798166, 0.0082690213201461298, 0.019417750785666502, 0.026365055544241726, 0.022656393339840427, 0.0057963362663038627, -0.019430374432658718, -0.040917773802925161, -0.043448047379980842, -0.015798906135882074, 0.042568547742062521, 0.11891573199913763, 0.19071173716754025, 0.23418328226423002, 0.23418328226423002, 0.19071173716754025, 0.11891573199913763, 0.042568547742062521, -0.015798906135882074, -0.043448047379980842, -0.040917773802925161, -0.019430374432658718, 0.0057963362663038627, 0.022656393339840427, 0.026365055544241726, 0.019417750785666502, 0.0082690213201461298, -0.00088958916661798166, -0.0049652330188863069, -0.0064902982704557882};
    constexpr std::array<double, 31> f_hil3{-0.03192367836850326, 1.6506630344403635e-15, -0.02604586500491711, -7.6878209964369998e-16, -0.037320916902723426, 7.4163962901512533e-17, -0.053092776177022238, -4.5948402104485316e-16, -0.076688664348758373, 9.6705760467325067e-17, -0.1168448104037082, -6.7313205849769761e-17, -0.20579817845397408, -2.7475916837786641e-16, -0.63446321100815162, 0, 0.63446321100815162, 2.7475916837786641e-16, 0.20579817845397408, 6.7313205849769761e-17, 0.1168448104037082, -9.6705760467325067e-17, 0.076688664348758373, 4.5948402104485316e-16, 0.053092776177022238, -7.4163962901512533e-17, 0.037320916902723426, 7.6878209964369998e-16, 0.02604586500491711, -1.6506630344403635e-15, 0.03192367836850326};
    constexpr std::array<double, 31> f_diff3{-0.045934609631201745, 0.020702179097299284, 0.056590005659056572, -0.060766779698937956, 0.051189643736576695, -0.029151469258650732, -0.0059182678567156338, 0.050238177362022435, -0.094680944827040009, 0.12562421936989346, -0.12662444644123005, 0.079487824660651671, 0.038808272576328934, -0.27797780682306306, 0.88060998519685807, 0, -0.88060998519685807, 0.27797780682306306, -0.038808272576328934, -0.079487824660651671, 0.12662444644123005, -0.12562421936989346, 0.094680944827040009, -0.050238177362022435, 0.0059182678567156338, 0.029151469258650732, -0.051189643736576695, 0.060766779698937956, -0.056590005659056572, -0.020702179097299284, 0.045934609631201745};
    constexpr std::array<double, 32> f_hil4{-0.026493621260073, -0.0057653007898626069, -0.016568074278603753, -0.0095335770782756027, -0.02247200528666584, -0.014961114020507524, -0.030313758532595853, -0.022995444428882961, -0.041419851777148436, -0.035735226446599294, -0.059119382195648218, -0.059012422412988874, -0.09448672865361657, -0.11737541546315011, -0.21800661329690929, -0.62866354497609622, 0.62866354497609622, 0.21800661329690929, 0.11737541546315011, 0.09448672865361657, 0.059012422412988874, 0.059119382195648218, 0.035735226446599294, 0.041419851777148436, 0.022995444428882961, 0.030313758532595853, 0.014961114020507524, 0.02247200528666584, 0.0095335770782756027, 0.016568074278603753, 0.0057653007898626069, 0.026493621260073};
    constexpr std::array<double, 32> f_diff4{-0.0039359603211662582, 0.0053802639494991153, -0.0026663524159915072, 0.0025047985479988199, -0.0027306590883832815, 0.0031458377569417214, -0.0037493377129640621, 0.0046049611849922466, -0.0058422348236550972, 0.0077058498842837884, -0.010685990088053862, 0.015875971655971054, -0.026136741498605901, 0.051078527312964903, -0.14161788085897356, 1.2733853224358698, -1.2733853224358698, 0.14161788085897356, -0.051078527312964903, 0.026136741498605901, -0.015875971655971054, 0.010685990088053862, -0.0077058498842837884, 0.0058422348236550972, -0.0046049611849922466, 0.0037493377129640621, -0.0031458377569417214, 0.0027306590883832815, -0.0025047985479988199, 0.0026663524159915072, -0.0053802639494991153, 0.0039359603211662582};
    const auto sym = Symmetry::symmetric;
    const auto anti = Symmetry::antisymmetric;
    const auto constant = Response::constant;
    const auto slope = Response::differentiator;
    const auto against = [](const Specification& spec, std::span<const double> fixture) {
        Exchange engine(spec);
        const auto r = engine.run();
        check(r.status == Status::ok && engine.taps().size() == fixture.size(), "fixture design converges");
        for (std::size_t k = 0; k < fixture.size(); ++k)
            check(std::abs(engine.taps()[k] - fixture[k]) <= 3e-5 * r.ripple + 1e-12, "independent SciPy coefficients");
    };
    against({31, hp, sym, constant, 32}, f_hp);
    against({31, bp, sym, constant, 32}, f_bp);
    against({31, bs, sym, constant, 32}, f_bs);
    against({15, multi, sym, constant, 32}, f_multi);
    against({32, lp, sym, constant, 32}, f_lp2);
    against({31, hilbert, anti, constant, 32}, f_hil3);
    against({31, diff3, anti, slope, 32}, f_diff3);
    against({32, hilbert, anti, constant, 32}, f_hil4);
    against({32, diff4, anti, slope, 32}, f_diff4);
    for (const std::size_t N : {1u, 3u, 5u, 15u, 31u, 63u}) {
        verify({N, hp, sym, constant, 32});
        verify({N, bp, sym, constant, 32});
        verify({N, bs, sym, constant, 32});
        verify({N, multi, sym, constant, 32});
    }
    for (const std::size_t N : {2u, 4u, 16u, 32u, 64u}) {
        verify({N, lp, sym, constant, 32});
        verify({N, bp, sym, constant, 32});
        verify({N, hilbert, anti, constant, 32});
        verify({N, diff4, anti, slope, 32});
    }
    for (const std::size_t N : {3u, 5u, 15u, 31u, 63u}) {
        verify({N, hilbert, anti, constant, 32});
        verify({N, diff3, anti, slope, 32});
    }
    verify({223, lp, sym, constant, 16});
    verify({224, lp, sym, constant, 16});
    verify({255, bp, sym, constant, 16});
    verify({223, hilbert, anti, constant, 16});
    verify({224, hilbert, anti, constant, 16});
    verify({223, diff3, anti, slope, 16});
    verify({255, diff3, anti, slope, 16});
    for (const auto& spec : {Specification{63, hilbert, anti}, Specification{223, hilbert, anti}, Specification{64, hilbert, anti},
                             Specification{223, lp}, Specification{255, bp}, Specification{223, diff3, anti, slope},
                             Specification{31, lp}}) {
        Exchange engine(spec);
        const double first = engine.run(1).ripple;
        const double final = engine.run(128, 1e-5).ripple;
        check(first >= 0.1 * final && first <= final * (1 + 1e-5), "first levelled error is a near-minimax lower bound");
    }
    for (const auto& [N, scipy] : {std::pair{223u, 9.125001e-9}, std::pair{224u, 8.738776e-9}}) {
        Exchange engine(Specification{N, lp, sym, constant, 16});
        const auto r = engine.run();
        check(r.status == Status::ok && r.max_error <= scipy && r.ripple >= 0.9 * scipy, "no worse than SciPy");
    }
    verify({401, std::array{Band{0, .24 * pi, 1, 1}, Band{.26 * pi, pi, 0, 10}}, sym, constant, 16});
    verify({402, std::array{Band{0, .24 * pi, 1, 1}, Band{.26 * pi, pi, 0, 10}}, sym, constant, 16});
    verify({801, std::array{Band{0, .245 * pi, 1, 1}, Band{.25 * pi, pi, 0, 10}}, sym, constant, 16});
    verify({511, std::array{Band{.01 * pi, .99 * pi, -1, 1}}, anti, constant, 16});
    verify({512, std::array{Band{.01 * pi, .99 * pi, -1, 1}}, anti, constant, 16});
    verify({319, lp, sym, constant, 16}, 1e-5);
    Exchange deep(Specification{320, lp, sym, constant, 16});
    check(deep.run(128, 1e-5).status != Status::ok && deep.taps().empty(), "binary64 limit reported");
    Exchange base(Specification{31, hilbert, anti});
    const auto unscaled = base.run();
    for (const double scale : {1e-100, 1e100}) {
        auto band = hilbert;
        band[0].weight *= scale;
        Exchange engine(Specification{31, band, anti});
        const auto r = engine.run();
        check(r.status == Status::ok && std::abs(r.ripple / scale - unscaled.ripple) <= 1e-10 * unscaled.ripple, "weight scaling");
        for (std::size_t k = 0; k < 31; ++k)
            check(std::abs(engine.taps()[k] - base.taps()[k]) <= 1e-12, "weight scaling keeps taps");
    }
    verify({63, std::array{Band{0, .2 * pi, 1, 1e-3}, Band{.3 * pi, pi, 0, 1e3}}, sym, constant, 32});
    const std::array same{Band{0, .2 * pi, 0.5, 1}, Band{.3 * pi, pi, 0.5, 10}};
    Exchange exact(Specification{31, same});
    check(exact.run().status == Status::ok && exact.taps()[15] == 0.5 && exact.extrema().empty() && amplitude(exact, 1.0) == 0.5,
          "exact constant response");
    for (const auto symmetry : {sym, anti}) {
        for (const std::size_t N : {30u, 31u}) {
            Exchange zero(Specification{N, std::array{Band{0, pi, 0, 1}}, symmetry});
            check(zero.run().status == Status::ok && std::ranges::all_of(zero.taps(), [](double v) { return v == 0.0; }),
                  "exact zero in every form");
        }
    }
}

void ownership() {
    const std::array lp{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, pi, 0, 10}};
    auto input = lp;
    Exchange owned(Specification{31, input});
    check(owned.run().status == Status::ok, "owned design");
    const std::vector<double> first(owned.taps().begin(), owned.taps().end());
    input[0].target = -100;
    check(owned.run().status == Status::ok && std::ranges::equal(owned.taps(), first), "input bands copied at setup");
    Exchange moved(std::move(owned));
    check(std::ranges::equal(moved.taps(), first) && amplitude(moved, .1) == amplitude(moved, .1), "move keeps the design");
    check(owned.status() == Status::invalid_specification && owned.taps().empty() && owned.extrema().empty() &&
              owned.grid().empty() && !owned.amplitude(.1),
          "moved-from engine exposes nothing");
    check(owned.run().status == Status::invalid_specification, "moved-from engine refuses to run");
    Exchange source(Specification{31, lp});
    check(source.run().status == Status::ok, "assignment source designed");
    Exchange sink(Lowpass{15, .2 * pi, .3 * pi});
    sink = std::move(source);
    check(!source.amplitude(.1) && source.taps().empty() && source.extrema().empty(), "move-assigned-from engine exposes nothing");
    Exchange target(Lowpass{15, .2 * pi, .3 * pi});
    target = std::move(moved);
    check(std::ranges::equal(target.taps(), first) && moved.status() == Status::invalid_specification && moved.taps().empty(),
          "move assignment transfers the design");
    owned = Exchange(Specification{31, lp});
    check(owned.run().status == Status::ok && std::ranges::equal(owned.taps(), first), "moved-from engine accepts a new design");
}

void invalid_inputs() {
    const std::array lp{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, pi, 0, 10}};
    const std::array full{Band{0, pi, 1, 1}};
    const auto expect = [](const Specification& spec, Status status, std::size_t max_points = Exchange::default_max_points) {
        Exchange engine(spec, max_points);
        check(engine.status() == status && engine.run().status == status && engine.grid().empty(), "invalid setup rejected without a grid");
    };
    const auto with = [&](auto edit) {
        auto bands = lp;
        edit(bands);
        return bands;
    };
    expect({0, lp}, Status::invalid_taps);
    expect({1, std::array{Band{.05 * pi, .95 * pi, 1, 1}}, Symmetry::antisymmetric}, Status::invalid_taps);
    expect({31, {}}, Status::invalid_bands);
    expect({31, lp, Symmetry::symmetric, Response::constant, 0}, Status::invalid_density);
    expect({31, lp, Symmetry::symmetric, Response::differentiator}, Status::invalid_specification);
    expect({31, lp, static_cast<Symmetry>(7)}, Status::invalid_specification);
    for (const auto& bands : {with([](auto& b) { b[1].lower = b[0].upper; }), with([](auto& b) { b[0].upper = b[0].lower; }),
                              with([](auto& b) { b[1].upper = 4.0; }), with([](auto& b) { b[0].lower = -1.0; }),
                              with([](auto& b) { b[0].upper = std::numeric_limits<double>::quiet_NaN(); }),
                              with([](auto& b) { std::swap(b[0], b[1]); })})
        expect({31, bands}, Status::invalid_bands);
    for (const double w : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        expect({31, with([&](auto& b) { b[1].weight = w; })}, Status::invalid_weights);
    expect({31, with([](auto& b) { b[0].target = std::numeric_limits<double>::infinity(); })}, Status::invalid_target);
    expect({32, full}, Status::incompatible_endpoint);
    expect({31, full, Symmetry::antisymmetric}, Status::incompatible_endpoint);
    expect({31, full, Symmetry::antisymmetric, Response::differentiator}, Status::incompatible_endpoint);
    expect({31, std::array{Band{.5 * pi, .5 * pi + 1e-9, 1, 1}}, Symmetry::symmetric, Response::constant, 1}, Status::invalid_density);
    expect({31, lp}, Status::grid_limit, 100);
    expect({std::numeric_limits<std::size_t>::max(), lp}, Status::invalid_taps, std::numeric_limits<std::size_t>::max());
    Exchange engine(Lowpass{31, .2 * pi, .3 * pi});
    for (const double bad : {-1.0, 0.0, 1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        check(engine.run(64, bad).status == Status::invalid_tolerance, "invalid tolerance");
    check(engine.run(1).status == Status::nonconvergence && engine.taps().empty(), "iteration budget exhaustion");
    check(engine.run().status == Status::ok, "recovers after failure");
    check(!engine.amplitude(-0.1) && !engine.amplitude(4.0) && !engine.amplitude(std::numeric_limits<double>::quiet_NaN()),
          "amplitude outside [0, pi] is empty");
    check(Exchange(Lowpass{32, .2 * pi, .3 * pi}).status() == Status::invalid_taps, "lowpass is odd-length Type I");
    check(Exchange(Lowpass{31, .3 * pi, .2 * pi}).status() == Status::invalid_bands, "lowpass edges must be ordered");
    Exchange type3(Specification{31, std::array{Band{0, .8 * pi, 1, 1}, Band{.9 * pi, pi, 0, 1}}, Symmetry::antisymmetric,
                                 Response::differentiator});
    const auto points = type3.grid();
    const auto at = [](double omega) { return [omega](const remez::Point& p) { return p.omega == omega; }; };
    check(points.front().omega == 0.0 && points.back().omega < pi && points.back().omega > .9 * pi, "DC kept, forced Nyquist left out");
    check(std::ranges::count_if(points, at(.8 * pi)) == 1 && std::ranges::count_if(points, at(.9 * pi)) == 1, "exact interior band endpoints");
    for (std::size_t j = 1; j < points.size(); ++j)
        check(points[j].x < points[j - 1].x && points[j].x == std::cos(points[j].omega), "decreasing x = cos(omega)");
}

#ifdef REMEZ_CRT_CHECKS
bool counting = false;
bool reject = false;
std::size_t allocations = 0;
std::size_t allocation_bytes = 0;
std::size_t accepted = 0;

int __cdecl allocation_hook(int kind, void*, std::size_t bytes, int block, long, const unsigned char*, int) {
    if (!counting || (block & 0xffff) == _CRT_BLOCK || (kind != _HOOK_ALLOC && kind != _HOOK_REALLOC)) return 1;
    ++allocations;
    allocation_bytes += bytes;
    return reject && allocations > accepted ? 0 : 1;
}

bool outstanding(const _CrtMemState& before) {
    _CrtMemState after{};
    _CrtMemState difference{};
    _CrtMemCheckpoint(&after);
    return _CrtMemDifference(&difference, &before, &after) != 0;
}

void memory_contract() {
    const auto previous = _CrtSetAllocHook(allocation_hook);
    _CrtMemState before{};
    _CrtMemCheckpoint(&before);
    {
        const std::array lp{Band{0, .2 * pi, 1, 1}, Band{.3 * pi, pi, 0, 10}};
        const std::array hilbert{Band{.05 * pi, .95 * pi, -1, 1}};
        const std::array diff{Band{0, pi, 1, 1}};
        for (const auto& spec : {Specification{31, lp}, Specification{32, lp}, Specification{31, hilbert, Symmetry::antisymmetric},
                                 Specification{32, diff, Symmetry::antisymmetric, Response::differentiator}}) {
            allocations = 0;
            allocation_bytes = 0;
            counting = true;
            Exchange engine(spec);
            counting = false;
            const auto M = engine.grid().size();
            const auto n = spec.symmetry == Symmetry::symmetric ? (spec.taps - 1) / 2 + 2 : spec.taps / 2 + 1;
            const auto payload = 32 * M + 32 * n + 48 * spec.bands.size() + 56 * n + 24 * M + 8 * (spec.taps / 2 + 1) + 8 * spec.taps;
            check(allocations == 7 && allocation_bytes >= payload && allocation_bytes <= payload + 7 * 64, "seven bounded setup buffers");
            for (int repeat = 0; repeat < 3; ++repeat) {
                allocations = 0;
                counting = true;
                const auto result = engine.run();
                const auto value = engine.amplitude(0.3 * pi);
                const auto taps = engine.taps();
                const auto extrema = engine.extrema();
                const auto exhausted = engine.run(1);
                const auto invalid = engine.run(64, 0.0);
                Exchange moved(std::move(engine));
                engine = std::move(moved);
                counting = false;
                check(allocations == 0 && result.status == Status::ok && value.has_value() && taps.size() == spec.taps &&
                          !extrema.empty() && exhausted.status == Status::nonconvergence && invalid.status == Status::invalid_tolerance,
                      "runs, queries and moves allocate nothing, on success or failure");
            }
        }
        allocations = 0;
        counting = true;
        Exchange invalid(Specification{0, lp});
        counting = false;
        check(allocations == 0 && invalid.status() == Status::invalid_taps, "invalid setup allocates nothing");
        for (std::size_t operation = 0; operation < 7; ++operation) {
            allocations = 0;
            accepted = operation;
            reject = true;
            counting = true;
            bool caught = false;
            try {
                Exchange engine(Specification{31, lp});
            } catch (const std::bad_alloc&) {
                caught = true;
            }
            reject = false;
            counting = false;
            check(caught, "allocation failure at setup propagates");
        }
    }
    check(!outstanding(before), "no leak after success, failure or allocation failure");
    allocations = 0;
    counting = true;
    auto temporary = std::make_unique<int>(42);
    temporary.reset();
    counting = false;
    check(allocations > 0, "allocation detector positive control");
    _CrtMemCheckpoint(&before);
    {
        const auto held = std::make_unique<int>(42);
        check(outstanding(before), "leak detector positive control");
    }
    check(!outstanding(before), "leak probe released");
    _CrtSetAllocHook(previous);
}
#endif

}

int main() {
#ifdef REMEZ_CRT_CHECKS
    memory_contract();
#endif
    invalid_inputs();
    ownership();
    scipy_references();
    lowpass();
    general();
    std::puts("PASS");
}
