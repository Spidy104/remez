set_project("remez")
set_languages("c++20")
add_rules("mode.debug", "mode.release")
set_warnings("allextra", "error")
set_runtimes("MD")
add_cxxflags("/permissive-", "/fp:precise", {tools = "cl"})
add_cxxflags("-Wpedantic", "-fno-fast-math", "-ffp-contract=off", {tools = "clang"})

option("memory_checks")
    set_default(false)
    set_showmenu(true)
option_end()
target("tests")
    set_kind("binary")
    add_files("remez.cpp", "tools/tests.cpp")

if has_config("memory_checks") then
    target("memory_tests")
        set_kind("binary")
        set_default(false)
        set_runtimes("MDd")
        add_defines("_DEBUG", "REMEZ_CRT_CHECKS", "_ITERATOR_DEBUG_LEVEL=0")
        add_files("remez.cpp", "tools/tests.cpp")
end

target("sanitizer_probe")
    set_kind("binary")
    set_default(false)
    add_files("tools/sanitizer_probe.cpp")

