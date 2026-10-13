#include <climits>
#include <cstdio>
#include <memory>
#include <span>
#include <string_view>

int main(int argc, char** argv) {
    const std::span arguments(argv, static_cast<std::size_t>(argc));
    const std::string_view mode = arguments.size() > 1 ? arguments[1] : "clean";
    if (mode == "overflow") {
        volatile int value = INT_MAX;
        std::printf("%d\n", value + 1);
    } else if (mode == "bounds") {
        auto data = std::make_unique<int[]>(4);
        volatile int index = argc + 2;
        data[index] = 42;
        std::printf("%d\n", data[index]);
    } else if (mode == "use-after-free") {
        auto data = std::make_unique<int[]>(4);
        int& borrowed = data[0];
        data.reset();
        borrowed = 42;
        std::printf("%d\n", borrowed);
    } else if (mode == "clean") {
        auto data = std::make_unique<int[]>(4);
        data[0] = 42;
        std::printf("%d\n", data[0]);
    } else {
        return 2;
    }
}
