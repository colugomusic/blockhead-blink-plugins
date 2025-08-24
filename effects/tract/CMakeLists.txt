cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-effect-tract VERSION 0.5.1)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
  ./src/dsp.hpp
  ./src/model.h
  ./src/plugin.cpp
	../../shared/tract_params.hpp
)

list(APPEND RES ./res/icon.png)
blink_plugin_add_effect("tract" "${SRC}" "${RES}")
target_include_directories(blockhead-blink-plugins-effect-tract PRIVATE ../..)
