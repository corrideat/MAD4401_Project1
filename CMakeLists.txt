cmake_minimum_required(VERSION 2.6)
project(PROJECT1)

set(CMAKE_C_FLAGS "-std=c99 -Wall -Wstrict-aliasing=2 -DDEBUG -g -ggdb -fstrict-aliasing -O1")

add_executable(project1 src/matrix.c src/utilities.c src/project1.c src/main.c)

target_link_libraries(project1 m)
