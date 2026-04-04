set(PRIVATE_COMPILE_OPTIONS -Wall -Wextra -Wno-sign-compare -Wno-implicit-fallthrough -fPIC)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND PRIVATE_COMPILE_OPTIONS -fext-numeric-literals -Wno-deprecated-declarations)
endif()

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Core library -- all production source except main.cpp and *Test.*
add_library(reseq_lib STATIC
  AdapterStats.cpp
  archive_format.cpp
  BamIngestionEngine.cpp
  CoverageStats.cpp
  DataStats.cpp
  DataStatsInterface.cpp
  ErrorStats.cpp
  FragmentDistributionStats.cpp
  FragmentDuplicationStats.cpp
  ProbabilityEstimates.cpp
  QualityStats.cpp
  ReadSequenceStats.cpp
  Reference.cpp
  Simulator.cpp
  Surrounding.cpp
  TileStats.cpp
)
target_include_directories(reseq_lib PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/reseq"
)
target_link_libraries(reseq_lib PUBLIC
  skewer_matrix nlopt reseq_seqan GTest::gtest
  ${Boost_LIBRARIES} ZLIB::ZLIB BZip2::BZip2 rt
)
target_compile_options(reseq_lib PRIVATE ${PRIVATE_COMPILE_OPTIONS})

# Production executable
add_executable(reseq main.cpp)
target_sources(reseq PRIVATE cli/cli_common.cpp cli/convert_profile.cpp cli/illumina_pe.cpp cli/query_profile.cpp cli/replace_n.cpp cli/seq_to_illumina.cpp)
target_link_libraries(reseq PRIVATE reseq_lib)
target_compile_options(reseq PRIVATE ${PRIVATE_COMPILE_OPTIONS})

install(TARGETS reseq
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

option(RESEQ_BUILD_TESTS "Build test binary" ON)
if(RESEQ_BUILD_TESTS)
  add_executable(reseq_test
    test_main.cpp
    AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp
    ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp
    FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp
    QualityStatsTest.cpp ReadSequenceStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp
    SimulatorTest.cpp SurroundingTest.cpp TileStatsTest.cpp
    utilitiesTest.cpp VectTest.cpp
    RegressionTest.cpp
  )
  target_link_libraries(reseq_test PRIVATE reseq_lib GTest::gtest)
  target_compile_options(reseq_test PRIVATE ${PRIVATE_COMPILE_OPTIONS})
  # Run all tests as a single invocation — tests have inter-test dependencies
  # (shared Register() state) that require running together, not individually.
  add_test(NAME reseq_tests COMMAND reseq_test)

  if(CODE_COVERAGE)
    target_compile_options(reseq_test PRIVATE --coverage -fprofile-update=atomic)
    target_link_options(reseq_test PRIVATE --coverage)
  endif()
endif()

if(CODE_COVERAGE)
  target_compile_options(reseq_lib PRIVATE --coverage -fprofile-update=atomic)
  target_link_options(reseq_lib PRIVATE --coverage)
  target_compile_options(reseq PRIVATE --coverage -fprofile-update=atomic)
  target_link_options(reseq PRIVATE --coverage)
endif()
