cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-zap VERSION 0.5.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
    ./src/dsp.hpp
    ./src/model.h
    ./src/plugin.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("zap" "${SRC}" "${RES}")
target_include_directories(blockhead-blink-plugins-effect-zap PRIVATE ../..)
