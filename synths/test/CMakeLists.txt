cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-synth-test VERSION 0.5.1)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
	./src/dsp.hpp
	./src/model.h
	./src/plugin.cpp
	../../shared/noise_generator.h
	../../shared/noise_generator.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_synth("test" "${SRC}" "${RES}")
target_include_directories(synth_test PRIVATE ../..)
