cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-template VERSION 0.1.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
	./src/model.h
	./src/dsp.hpp
	./src/plugin.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("template" "${SRC}" "${RES}")
target_include_directories(blockhead-blink-plugins-effect-template PRIVATE ../..)
