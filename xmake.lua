add_rules("mode.debug", "mode.release")
set_languages("c++23","c17")

set_encodings("utf-8")

set_version("0.0.1")-- 设置工程版本

add_requires("raylib") 

target("泊松分布")
    set_kind("binary")
    add_includedirs("include/")
    add_files("src/*.cpp")
    add_packages("raylib")
