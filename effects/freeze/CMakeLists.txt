cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-freeze VERSION 0.6.0)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
	./src/dsp.hpp
	./src/model.h
	./src/buffer.h
	./src/plugin.cpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("freeze" "${SRC}" "${RES}")
target_include_directories(blockhead-blink-plugins-effect-freeze PRIVATE ../..)