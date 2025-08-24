cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-filter_lowpass VERSION 0.5.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
	./src/dsp.hpp
	./src/model.h
	./src/plugin.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("filter.lowpass" "${SRC}" "${RES}")
target_include_directories(effect_filter.lowpass PRIVATE ../..)