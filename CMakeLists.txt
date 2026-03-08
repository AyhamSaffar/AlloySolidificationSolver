cmake_minimum_required(VERSION 3.28)
project(AlloySolidificationSolver)
set(CMAKE_CXX_STANDARD 20)


add_subdirectory("libs/Enzyme/enzyme")
# find_package(Enzyme REQUIRED CONFIG)


add_executable(AlloySolidificationSolver src/main.cpp)
target_include_directories(AlloySolidificationSolver PUBLIC include)
target_link_libraries(AlloySolidificationSolver PUBLIC Enzyme)

# target_link_libraries(AlloySolidificationSolver PUBLIC Enzyme)

