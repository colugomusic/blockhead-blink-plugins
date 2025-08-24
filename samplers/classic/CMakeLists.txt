cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-sampler-classic VERSION 0.5.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
    ./src/draw.hpp
    ./src/dsp.hpp
    ./src/model.h
    ./src/plugin.cpp
	../../shared/noise_generator.h
	../../shared/noise_generator.cpp
)

blink_plugin_add_sampler("classic" "${SRC}" "")
target_include_directories(sampler_classic PRIVATE ../..)
