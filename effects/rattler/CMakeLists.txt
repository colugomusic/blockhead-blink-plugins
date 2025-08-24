cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-rattler VERSION 0.3.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
    ./src/dsp.hpp
    ./src/model.h
    ./src/plugin.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("rattler" "${SRC}" "${RES}")
target_include_directories(effect_rattler PRIVATE ../..)
