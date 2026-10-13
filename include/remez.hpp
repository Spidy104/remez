#ifndef REMEZ_HPP
#define REMEZ_HPP

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace remez {

enum class Status {
    ok,
    invalid_specification,
    invalid_taps,
    invalid_bands,
    invalid_weights,
    invalid_target,
    incompatible_endpoint,
    invalid_density,
    grid_limit,
    invalid_tolerance,
    numerical_failure,
    nonconvergence
};

[[nodiscard]] constexpr std::string_view to_string(Status status) noexcept {
    switch (status) {
    case Status::ok: return "ok";
    case Status::invalid_specification: return "invalid_specification";
    case Status::invalid_taps: return "invalid_taps";
    case Status::invalid_bands: return "invalid_bands";
    case Status::invalid_weights: return "invalid_weights";
    case Status::invalid_target: return "invalid_target";
    case Status::incompatible_endpoint: return "incompatible_endpoint";
    case Status::invalid_density: return "invalid_density";
    case Status::grid_limit: return "grid_limit";
    case Status::invalid_tolerance: return "invalid_tolerance";
    case Status::numerical_failure: return "numerical_failure";
    case Status::nonconvergence: return "nonconvergence";
    }
    return "unknown";
}

enum class Symmetry { symmetric, antisymmetric, automatic };
enum class Response { constant, differentiator };

struct Band {
    double lower;
    double upper;
    double target;
    double weight{1.0};
};

struct Specification {
    std::size_t taps;
    std::span<const Band> bands;
    Symmetry symmetry{Symmetry::symmetric};
    Response response{Response::constant};
    std::size_t density{16};
};

struct Lowpass {
    std::size_t taps;
    double pass_edge;
    double stop_edge;
    std::array<double, 2> weights{1.0, 1.0};
    std::size_t density{16};
};

struct Point {
    double omega;
    double x;
    double desired;
    double weight;
};

struct Result {
    Status status;
    std::size_t iterations{};
    double ripple{};
    double max_error{};
    double flatness{};
};

class Exchange {
public:
    static constexpr std::size_t default_max_points = 1'000'000;

    explicit Exchange(const Specification& spec, std::size_t max_points = default_max_points);
    explicit Exchange(const Lowpass& spec, std::size_t max_points = default_max_points);
    Exchange(const Exchange&) = delete;
    Exchange& operator=(const Exchange&) = delete;
    Exchange(Exchange&& other) noexcept;
    Exchange& operator=(Exchange&& other) noexcept;
    ~Exchange() = default;

    [[nodiscard]] Status status() const noexcept { return status_; }
    [[nodiscard]] Symmetry symmetry() const noexcept { return symmetry_; }
    [[nodiscard]] Result run(std::size_t max_iterations = 64, double tolerance = 1e-6) noexcept;
    [[nodiscard]] std::span<const double> taps() const noexcept;
    [[nodiscard]] std::optional<double> amplitude(double omega) const noexcept;
    [[nodiscard]] std::span<const Point> extrema() const noexcept;
    [[nodiscard]] std::span<const Point> grid() const noexcept;

private:
    struct Range {
        Band band;
        std::size_t begin;
        std::size_t end;
    };
    struct Node {
        double omega;
        double x;
        double w;
        double y;
        double residual;
        double desired;
        double shift;
    };
    struct Peak {
        std::size_t index;
        double omega;
        double error;
    };
    enum class Values { fitted, residual };

    [[nodiscard]] Point point(const Band& band, double omega) const noexcept;
    [[nodiscard]] const Range& range(std::size_t index) const noexcept;
    [[nodiscard]] double interpolate(double x, Values values = Values::fitted) const noexcept;
    [[nodiscard]] double error(const Point& p) const noexcept;
    [[nodiscard]] bool equilibrium() noexcept;
    void seed() noexcept;
    [[nodiscard]] Status fit(double& delta) noexcept;
    void refine(Peak& peak) const noexcept;
    void select(std::span<Peak>& peaks) const noexcept;
    void inverse(Values values, bool accumulate) noexcept;
    [[nodiscard]] double tap_polynomial(double omega) const noexcept;
    [[nodiscard]] Status recover() noexcept;

    std::vector<Point> grid_;
    std::vector<Point> reference_;
    std::vector<Range> bands_;
    std::vector<Node> nodes_;
    std::vector<Peak> peaks_;
    std::vector<double> samples_;
    std::vector<double> taps_;
    Symmetry symmetry_{};
    Response response_{};
    unsigned type_{};
    Status status_{Status::invalid_specification};
    bool ready_{};
    bool exact_{};
    double constant_{};
};

}

#endif
