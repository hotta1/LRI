cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(scheduler NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C++ compilers
enable_language(CXX)

# LRSW
add_executable(LRSW LRSW.C LRSW_worker.C)
#add_executable(histogram histogram.C)
target_link_libraries(LRSW ${ALPS_LIBRARIES})
#target_link_libraries(histogram ${ALPS_LIBRARIES})
add_alps_test(LRSW)
