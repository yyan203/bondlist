cmake_minimum_required(VERSION 3.9)
project(bondlist)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(bondlist
        atom.h
        bond.h
        main.cpp
        readdata.h
        simulationbox.h
        structureanalysis.h)
