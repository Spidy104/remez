set_project("remez")
set_version("1.0.0")
set_languages("c++20")
add_rules("mode.debug", "mode.release")
set_warnings("allextra", "error")
add_cxxflags("/permissive-", "/fp:precise", {tools = "cl"})
add_cxxflags("-Wpedantic", "-fno-fast-math", "-ffp-contract=off", {tools = "clang"})

option("memory_checks")
    set_default(false)
    set_showmenu(true)
option_end()

if has_config("memory_checks") then
    set_runtimes("MDd")
    add_defines("_DEBUG", "_ITERATOR_DEBUG_LEVEL=0")
else
    set_runtimes("MD")
end

target("remez")
    set_kind("static")
    add_files("src/*.cpp")
    add_includedirs("include", {public = true})
    add_headerfiles("include/remez.hpp")

target("tests")
    set_kind("binary")
    set_default(not has_config("memory_checks"))
    add_deps("remez")
    add_files("tools/tests.cpp")

target("example")
    set_kind("binary")
    set_default(not has_config("memory_checks"))
    add_deps("remez")
    add_files("examples/design.cpp")

target("benchmark")
    set_kind("binary")
    set_default(false)
    add_deps("remez")
    add_files("tools/benchmark.cpp")

if has_config("memory_checks") then
    target("memory_tests")
        set_kind("binary")
        add_deps("remez")
        add_defines("REMEZ_CRT_CHECKS")
        add_files("tools/tests.cpp")
end

target("sanitizer_probe")
    set_kind("binary")
    set_default(false)
    add_files("tools/sanitizer_probe.cpp")
