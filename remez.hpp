#ifndef REMEZ_HPP
#define REMEZ_HPP

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace remez {

enum class Status {
    ok, invalid_taps, invalid_edges, invalid_weights, invalid_density,
    grid_limit, insufficient_points, insufficient_storage, numerical_resolution,
    invalid_nodes, numerical_failure, invalid_tolerance, nonconvergence,
    invalid_bands, invalid_target, incompatible_endpoint, invalid_specification
};

// Odd-length (taps = 2r+1) symmetric Type-I lowpass. Edges are angular frequencies in radians per
// sample, 0 < pass_edge < stop_edge < pi (pi is Nyquist). Desired gain is 1 on [0, pass_edge] and
// 0 on [stop_edge, pi]; weights are {passband, stopband}. Grid spacing is pi/(density*(r+1)).
struct Lowpass {
    std::size_t taps;
    double pass_edge;
    double stop_edge;
    std::array<double, 2> weights{1.0, 1.0};
    std::size_t density{16};
};

enum class Symmetry { symmetric, antisymmetric };
enum class Response { constant, differentiator };

struct Band {
    double lower, upper, target, weight{1.0};
};

struct Specification {
    std::size_t taps;
    std::span<const Band> bands;
    Symmetry symmetry{Symmetry::symmetric};
    Response response{Response::constant};
    std::size_t density{16};
};

struct BandRange {
    Band band;
    std::size_t begin{}, end{};
};

struct GridPoint {
    double omega, x, desired, weight;
};

struct GridSize {
    Status status;
    std::array<std::size_t, 2> band_points{};
    [[nodiscard]] std::size_t count() const noexcept {
        return band_points[0] + band_points[1];
    }
};

struct Grid {
    Status status{Status::invalid_taps};
    std::vector<GridPoint> points;
    std::size_t passband_points{};
    [[nodiscard]] std::span<const GridPoint> passband() const noexcept {
        return std::span{points}.first(passband_points);
    }
    [[nodiscard]] std::span<const GridPoint> stopband() const noexcept {
        return std::span{points}.subspan(passband_points);
    }
};

[[nodiscard]] GridSize grid_size(const Lowpass& spec,
                                std::size_t max_points = 1'000'000) noexcept;
[[nodiscard]] Status fill_grid(const Lowpass& spec,
                               std::span<GridPoint> output) noexcept;
[[nodiscard]] Grid make_grid(const Lowpass& spec,
                              std::size_t max_points = 1'000'000);

struct DesignSize {
    Status status;
    std::size_t points{}, degree{}, peak_capacity{}, workspace{};
};

struct DesignGrid {
    Status status{Status::invalid_specification};
    std::vector<GridPoint> points;
    std::vector<BandRange> bands;
    std::size_t degree{}, peak_capacity{}, workspace{};
    [[nodiscard]] std::span<const GridPoint> band(std::size_t i) const noexcept {
        return i < bands.size() ? std::span{points}.subspan(bands[i].begin, bands[i].end-bands[i].begin) : std::span<const GridPoint>{};
    }
};

[[nodiscard]] DesignSize grid_size(const Specification& spec, std::size_t max_points = 1'000'000) noexcept;
[[nodiscard]] Status fill_grid(const Specification& spec, std::span<GridPoint> output, std::span<BandRange> ranges) noexcept;
[[nodiscard]] DesignGrid make_grid(const Specification& spec, std::size_t max_points = 1'000'000);

namespace detail {
struct Peak { std::size_t index; double omega, error; };
// Internal kernel: finite nonzero peaks, ordered and alternating. Scratch has peaks.size() entries.
[[nodiscard]] std::size_t select_peaks(std::span<Peak> peaks, std::size_t n, std::span<double> scratch) noexcept;
}

struct Value {
    Status status;
    double value{};
};

struct Bounded {
    Status status;
    double value{}, bound{};  // |value - exact| <= bound; bound may be infinite.
};

class Interpolator {
    friend class Exchange;
    struct Node;
    struct PreciseValue;
    struct Location;
    std::vector<Node> nodes_;
    double value_scale_{};
    bool ready_{};
    [[nodiscard]] bool usable(double x) const noexcept;
    [[nodiscard]] Location locate(double x) const noexcept;
    [[nodiscard]] Bounded fast(double x) const noexcept;
    [[nodiscard]] PreciseValue precise(double x) const noexcept;

public:
    explicit Interpolator(std::size_t node_count);
    ~Interpolator();
    Interpolator(const Interpolator&) = delete;
    Interpolator& operator=(const Interpolator&) = delete;
    Interpolator(Interpolator&&) noexcept;
    Interpolator& operator=(Interpolator&&) noexcept;

    [[nodiscard]] Value fit(std::span<const GridPoint> extrema) noexcept;
    [[nodiscard]] Value evaluate(double x) const noexcept;
    // A(x) rounded to binary64 with a bound on its error from the exact interpolant. Uses the
    // two-component path when the binary64 bound exceeds 2^-46 of the node-value scale.
    [[nodiscard]] Bounded estimate(double x) const noexcept;
    // W*(D-A). A binary64 result is returned only when its rounding bound is below |E| and within
    // max(relative_accuracy*|E|, absolute_accuracy), with relative_accuracy in (0, 1) and
    // absolute_accuracy >= 0; otherwise it is computed in two components.
    [[nodiscard]] Value weighted_error(const GridPoint& point, double relative_accuracy = 0x1p-32,
                                       double absolute_accuracy = 0.0) const noexcept;
    // Recover 2*node_count-3 Type-I taps from the fitted degree-(node_count-2) amplitude.
    // No allocations; extra output elements are untouched. Discard partial output on failure.
    [[nodiscard]] Status taps(std::span<double> output) const noexcept;
    // Direct real inverse DFT for the form inferred from symmetry and count. Degree must match.
    [[nodiscard]] Status taps(std::span<double> output, Symmetry symmetry, std::size_t count) const noexcept;
};

struct ExchangeResult {
    Status status;
    std::size_t iterations{};
    double ripple{}, max_error{}, flatness{};
};

class Exchange {
    DesignGrid grid_;
    Interpolator polynomial_;
    std::vector<GridPoint> extrema_;
    std::vector<detail::Peak> peaks_;
    // Cosine table, sample/selection scratch, coefficients; after convergence scratch holds taps.
    std::vector<double> chebyshev_;
    std::size_t taps_{};
    Symmetry symmetry_{};
    Response response_{};
    unsigned type_{};
    Status tap_status_{Status::invalid_nodes};
    bool ready_{}, exact_{};
    double constant_{}, min_weight_{};
    [[nodiscard]] Status recover(std::span<const double> coefficients, bool certified) noexcept;

public:
    explicit Exchange(const Lowpass& spec, std::size_t max_points = 1'000'000);
    explicit Exchange(const Specification& spec, std::size_t max_points = 1'000'000);
    Exchange(const Exchange&) = delete;
    Exchange& operator=(const Exchange&) = delete;
    Exchange(Exchange&&) noexcept = default;
    Exchange& operator=(Exchange&&) noexcept = default;

    [[nodiscard]] ExchangeResult run(std::size_t max_iterations = 64,
                                     double tolerance = 1e-6) noexcept;
    // Transformed polynomial P(x); for Type I it is the physical amplitude.
    [[nodiscard]] Value evaluate(double x) const noexcept;
    // Signed amplitude in radians/sample: H=exp(-i*alpha*w)*A, or i*exp(-i*alpha*w)*A for antisymmetry.
    [[nodiscard]] Value amplitude(double omega) const noexcept;
    // Borrowed fitted points contain transformed desired values and weights.
    [[nodiscard]] std::span<const GridPoint> extrema() const noexcept;
    // Writes the spec.taps impulse-response taps h[0..taps-1] into output[0..taps-1]
    // without allocating; extra elements are untouched. run() prepares taps in existing workspace.
    // Requires a successful run(); a recovery failure does not invalidate the fitted amplitude.
    // On failure, discard any partial output.
    [[nodiscard]] Status taps(std::span<double> output) const noexcept;
};

}

#endif


