cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# file(STRINGS "VERSION" VERSION_STRING)
project(EHT-C LANGUAGES C)

# help macro
macro(log var)
  message(STATUS "${var}: ${${var}}")
endmacro()

# set building options
set(CMAKE_VERBOSE_MAKEFILE "FALSE" CACHE STRING "Cmake verbose output")
option(BUILD_TESTING "Build tests" "ON")

if("${CMAKE_BUILD_TYPE} " STREQUAL " ") # workaround to keep Release as default
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Compiler flags mode" FORCE) # the altenratives are Test and Coverage
endif()

if("${FULL_STACK}" STREQUAL "ON") # workaround to keep Release as default
  set(MEMORY_MODE "FULL_STACK" CACHE STRING "Memory mode" FORCE)
else()
  set(MEMORY_MODE "STACK_HEAP" CACHE STRING "Memory mode" FORCE)
endif()

log(PROTOCOL_VERSION)
log(CMAKE_COMPILER_IS_GNUCC)
log(CMAKE_BUILD_TYPE)
log(CMAKE_VERBOSE_MAKEFILE)
log(BUILD_TESTING)
log(BUILD_BENCHMARK)

# set directories names
set(EHT_LIGHT_A "${CMAKE_BINARY_DIR}/EHT_light_A")
set(EHT_LIGHT_B "${CMAKE_BINARY_DIR}/EHT_light_B")
set(EHT_MEDIUM_A "${CMAKE_BINARY_DIR}/EHT_medium_A")
set(EHT_MEDIUM_B "${CMAKE_BINARY_DIR}/EHT_medium_B")
set(EHT_HIGH_A "${CMAKE_BINARY_DIR}/EHT_high_A")
set(EHT_HIGH_B "${CMAKE_BINARY_DIR}/EHT_high_B")


set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")
set(BINARY_DIR "lib")
set(OUTPUT_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include")
set(TEST_BINARY_DIR "${CMAKE_BINARY_DIR}/test")
set(COVERAGE_DIR "${CMAKE_BINARY_DIR}/coverage")
set(TEST_PREFIX "${CMAKE_BINARY_DIR}/Testing")

log(TEST_VECTORS_DIR)

file(GLOB SOURCES ${SOURCE_DIR}/*.c)

# Other not used flags: -fsanitize=address -Wformat -Wformat-security -Wall -Wextra
if(CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -O3")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Test")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -std=c11 -g -O0 -m64 -D_FILE_OFFSET_BITS=64 -static-libgcc -Wformat=0 -Werror -Wno-unused-function -Wno-unused-result -Wno-strict-prototypes -Wcast-align -Wunused-value -Wunused-variable -Wundef -pedantic")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -std=c11 -g --coverage -O0 -m64 -D_FILE_OFFSET_BITS=64 -static-libgcc -Wformat=0 -Werror -Wno-unused-function -Wno-unused-result -Wno-strict-prototypes -Wcast-align -Wunused-value -Wunused-variable -Wundef -pedantic")
    else()
    	message( FATAL_ERROR "CMAKE_BUILD_TYPE not valid")
    endif()
endif(CMAKE_COMPILER_IS_GNUCC)
log(CMAKE_C_COMPILER)
log(CMAKE_C_FLAGS)


file(MAKE_DIRECTORY "${EHT_LIGHT_A}")
file(MAKE_DIRECTORY "${EHT_LIGHT_B}")
file(MAKE_DIRECTORY "${EHT_MEDIUM_A}")
file(MAKE_DIRECTORY "${EHT_MEDIUM_B}")
file(MAKE_DIRECTORY "${EHT_HIGH_A}")
file(MAKE_DIRECTORY "${EHT_HIGH_B}")

# copy configure file
set(PROTOCOL_VERSION "EHT_light_A")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_LIGHT_A}/config.h")
set(PROTOCOL_VERSION "EHT_light_B")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_LIGHT_B}/config.h")
set(PROTOCOL_VERSION "EHT_medium_A")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_MEDIUM_A}/config.h")
set(PROTOCOL_VERSION "EHT_medium_B")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_MEDIUM_B}/config.h")
set(PROTOCOL_VERSION "EHT_high_A")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_HIGH_A}/config.h")
set(PROTOCOL_VERSION "EHT_high_B")
configure_file("${INCLUDE_DIR}/config.h.in" "${EHT_HIGH_B}/config.h")

# make include directory
file(MAKE_DIRECTORY "${OUTPUT_INCLUDE_DIR}")

# copy other source files to directory
file(GLOB PUBLIC_HEADERS
  "${INCLUDE_DIR}/*.h"
)

# copy headers
file(COPY ${PUBLIC_HEADERS}/
	DESTINATION ${OUTPUT_INCLUDE_DIR}/
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	)

# add library
add_library(ehtla ${SOURCES})
add_library(ehtlb ${SOURCES})
add_library(ehtma ${SOURCES})
add_library(ehtmb ${SOURCES})
add_library(ehtha ${SOURCES})
add_library(ehthb ${SOURCES})

# include directory
target_include_directories(ehtla PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_LIGHT_A}")
target_include_directories(ehtlb PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_LIGHT_B}")
target_include_directories(ehtma PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_MEDIUM_A}")
target_include_directories(ehtmb PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_MEDIUM_B}")
target_include_directories(ehtha PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_HIGH_A}")
target_include_directories(ehthb PUBLIC ${OUTPUT_INCLUDE_DIR} "${EHT_HIGH_B}")

# target math library
target_link_libraries(ehtla m)
target_link_libraries(ehtlb m)
target_link_libraries(ehtma m)
target_link_libraries(ehtmb m)
target_link_libraries(ehtha m)
target_link_libraries(ehthb m)

#set output library location
set_target_properties(ehtla
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_LIGHT_A}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_LIGHT_A}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_LIGHT_A}/${BINARY_DIR}"
)
set_target_properties(ehtlb
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_LIGHT_B}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_LIGHT_B}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_LIGHT_B}/${BINARY_DIR}"
)
set_target_properties(ehtma
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_MEDIUM_A}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_MEDIUM_A}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_MEDIUM_A}/${BINARY_DIR}"
)
set_target_properties(ehtmb
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_MEDIUM_B}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_MEDIUM_B}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_MEDIUM_B}/${BINARY_DIR}"
)
set_target_properties(ehtha
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_HIGH_A}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_HIGH_A}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_HIGH_A}/${BINARY_DIR}"
)
set_target_properties(ehthb
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${EHT_HIGH_B}/${BINARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${EHT_HIGH_B}/${BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${EHT_HIGH_B}/${BINARY_DIR}"
)


if(BUILD_TESTING STREQUAL "ON")
  message(STATUS "Build tests")
  enable_testing()
  add_subdirectory(test)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  # target coverage is only build when `make coverage` is issued
  add_custom_target(coverage
      COMMAND ${CMAKE_COMMAND} -E make_directory ${COVERAGE_DIR}
      # COMMAND lcov --capture -i --directory . --output-file ${COVERAGE_DIR}/coverage_base.info
      COMMAND lcov --capture --directory . --output-file ${COVERAGE_DIR}/coverage.info
      # COMMAND lcov -a ${COVERAGE_DIR}/coverage_base.info -a ${COVERAGE_DIR}/coverage.info -o ${COVERAGE_DIR}/coverage.info
      COMMAND lcov --remove ${COVERAGE_DIR}/coverage.info -o ${COVERAGE_DIR}/coverage.info "${TEST_DIR}/*" "${OUTPUT_INCLUDE_DIR}/*"
      COMMAND genhtml ${COVERAGE_DIR}/coverage.info --output-directory ${COVERAGE_DIR}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

