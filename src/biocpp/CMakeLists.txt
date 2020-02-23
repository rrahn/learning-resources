project (biostec CXX)
cmake_minimum_required(VERSION 3.4)

find_package(SeqAn3 3.0.1 REQUIRED HINTS ../../submodules/seqan3/build_system/)

add_executable(hello_world hello_world.cpp)
target_link_libraries(hello_world seqan3::seqan3)
