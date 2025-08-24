cmake_minimum_required(VERSION 3.30)
project(blockhead-blink-plugins-sampler-fudge VERSION 0.5.1)

find_package(blink-plugin REQUIRED CONFIG)

list(APPEND SRC
    ./src/audio_data.h
    ./src/convert.h
    ./src/draw.hpp
    ./src/dsp.hpp
    ./src/grain.h
    ./src/model.h
    ./src/plugin.cpp
	./src/sample_analysis.h
	./src/sample_analysis.cpp
	../../shared/noise_generator.h
	../../shared/noise_generator.cpp
)

blink_plugin_add_sampler("fudge" "${SRC}" "")
target_include_directories(blockhead-blink-plugins-sampler-fudge PRIVATE ..  ../..)
