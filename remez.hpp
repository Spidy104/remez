#ifndef REMEZ_HPP
#define REMEZ_HPP

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace remez {

enum class Status {
    ok, invalid_specification, invalid_taps, invalid_bands, invalid_weights, invalid_target,
    incompatible_endpoint, invalid_density, grid_limit, invalid_tolerance, numerical_failure, nonconvergence
};

enum class Symmetry { symmetric, antisymmetric };
enum class Response { constant, differentiator };

// Frequencies are radians/sample in [0, pi]. A constant band asks for A(w) = target with error
// weight*(target - A). A differentiator band with target s != 0 asks for A(w) = s*w with relative
// error weight*(s*w - A)/(|s|*w); zero-target bands stay ordinary stopbands.
struct Band {
    double lower, upper, target, weight{1.0};
};

// Odd/even length and symmetry select Types I-IV. Grid spacing is pi/(density*(degree+1)).
struct Specification {
    std::size_t taps;
    std::span<const Band> bands;
    Symmetry symmetry{Symmetry::symmetric};
    Response response{Response::constant};
    std::size_t density{16};
};

// Odd-length Type-I lowpass: gain 1 on [0, pass_edge], 0 on [stop_edge, pi], weights {pass, stop}.
struct Lowpass {
    std::size_t taps;
    double pass_edge, stop_edge;
    std::array<double, 2> weights{1.0, 1.0};
    std::size_t density{16};
};

// A grid or reference point of the transformed problem P(x) ~ desired at x = cos(omega), where
// A(omega) = Q(omega) P(x) and Q is 1, cos(w/2), sin(w), sin(w/2) for Types I-IV.
struct Point {
    double omega, x, desired, weight;
};

struct Result {
    Status status;
    std::size_t iterations{};
    double ripple{}, max_error{}, flatness{};
};

// Parks-McClellan exchange. The constructor validates the specification and makes every
// allocation; run(), taps(), amplitude() and extrema() never allocate.
class Exchange {
public:
    explicit Exchange(const Specification& spec, std::size_t max_points = 1'000'000);
    explicit Exchange(const Lowpass& spec, std::size_t max_points = 1'000'000);

    [[nodiscard]] Status status() const noexcept { return status_; }
    [[nodiscard]] Result run(std::size_t max_iterations = 64, double tolerance = 1e-6) noexcept;
    // Impulse response h[0..taps-1]; empty unless the last run converged.
    [[nodiscard]] std::span<const double> taps() const noexcept;
    // Real amplitude A(omega) with H = exp(-i*alpha*w)*A (times i if antisymmetric); NaN unless the
    // last run converged or omega is outside [0, pi].
    [[nodiscard]] double amplitude(double omega) const noexcept;
    // Final reference in transformed units; empty unless converged to a nonzero error.
    [[nodiscard]] std::span<const Point> extrema() const noexcept;
    [[nodiscard]] std::span<const Point> grid() const noexcept { return grid_; }

private:
    struct Range { Band band; std::size_t begin, end; };
    struct Node { double omega, x, w, y, residual, desired, shift; };
    struct Peak { std::size_t index; double omega, error; };

    [[nodiscard]] double interpolate(double x, double Node::*value = &Node::y) const noexcept;
    [[nodiscard]] double error(const Point& p) const noexcept;
    [[nodiscard]] double tap_polynomial(double omega) const noexcept;
    void inverse(double Node::*value, bool accumulate) noexcept;
    [[nodiscard]] Point point(const Band& band, double omega) const noexcept;
    [[nodiscard]] const Range& range(std::size_t index) const noexcept;
    [[nodiscard]] Status fit(double& delta) noexcept;
    void refine(Peak& peak) const noexcept;
    [[nodiscard]] Status recover() noexcept;

    std::vector<Point> grid_, reference_;
    std::vector<Range> bands_;
    std::vector<Node> nodes_;
    std::vector<Peak> peaks_;
    std::vector<double> samples_, taps_;
    Symmetry symmetry_{};
    Response response_{};
    unsigned type_{};
    Status status_{Status::invalid_specification};
    bool ready_{}, exact_{};
    double constant_{};
};

}

#endif
