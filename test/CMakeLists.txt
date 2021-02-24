cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# custom function for testing
function(my_add_test name source dependency1 dependency2 expected_response )

  add_executable(${name} 
  	"${source}")

  target_link_libraries(${name} 
  	PRIVATE ${dependency1} ${dependency2}
  	)

  add_test(test_${name} ${name} 
  	WORKING_DIRECTORY ${TEST_BINARY_DIR}
  	)

  set_tests_properties(test_${name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${expected_response}"
  	)

endfunction()

link_directories(${BINARY_DIR})

# test_correctness
my_add_test(EHT_light_A_correctness "${TEST_DIR}/test_correctness.c" m ehtla "Test passed")
my_add_test(EHT_light_B_correctness "${TEST_DIR}/test_correctness.c" m ehtlb "Test passed")
my_add_test(EHT_medium_A_correctness "${TEST_DIR}/test_correctness.c" m ehtma "Test passed")
my_add_test(EHT_medium_B_correctness "${TEST_DIR}/test_correctness.c" m ehtmb "Test passed")
my_add_test(EHT_high_A_correctness "${TEST_DIR}/test_correctness.c" m ehtha "Test passed")
my_add_test(EHT_high_B_correctness "${TEST_DIR}/test_correctness.c" m ehthb "Test passed")

# bench_cycles
add_executable (EHT_light_A_bench_cycles "${TEST_DIR}/bench_cycles.c")
add_executable (EHT_light_B_bench_cycles "${TEST_DIR}/bench_cycles.c")
add_executable (EHT_medium_A_bench_cycles "${TEST_DIR}/bench_cycles.c")
add_executable (EHT_medium_B_bench_cycles "${TEST_DIR}/bench_cycles.c")
add_executable (EHT_high_A_bench_cycles "${TEST_DIR}/bench_cycles.c")
add_executable (EHT_high_B_bench_cycles "${TEST_DIR}/bench_cycles.c")

target_link_libraries (EHT_light_A_bench_cycles LINK_PUBLIC m ehtla)
target_link_libraries (EHT_light_B_bench_cycles LINK_PUBLIC m ehtlb)
target_link_libraries (EHT_medium_A_bench_cycles LINK_PUBLIC m ehtma)
target_link_libraries (EHT_medium_B_bench_cycles LINK_PUBLIC m ehtmb)
target_link_libraries (EHT_high_A_bench_cycles LINK_PUBLIC m ehtha)
target_link_libraries (EHT_high_B_bench_cycles LINK_PUBLIC m ehthb)
