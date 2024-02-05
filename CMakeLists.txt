cmake_minimum_required(VERSION 3.27)
project(CHEM_179_HW1)


set(CMAKE_CXX_STANDARD 23)
#set(ARMADILLO_LIBRARY "/Users/vittor/Downloads/armadillo-12.6.7/lib")
#set(ARMADILLO_INCLUDE_DIR "/Users/vittor/Downloads/armadillo-12.6.7/include")

include_directories(/Users/vittor/Downloads/armadillo-12.6.7/include)



add_executable(CHEM_179_HW1
        HW1/hw1.cpp
        HW1/hw1.h
)

#Find LAPACK
find_package(LAPACK REQUIRED)
find_library(ARMA_LIBRARY armadillo)
target_link_libraries(CHEM_179_HW1 ${ARMA_LIBRARY} ${LAPACK_LIBRARIES})
