# Own libraries
set(PRIVATE_COMPILE_OPTIONS -Wall -Wextra -Wno-sign-compare -Wno-implicit-fallthrough)

add_library(DataStats STATIC DataStats.cpp AdapterStats.cpp CoverageStats.cpp ErrorStats.cpp FragmentDistributionStats.cpp FragmentDuplicationStats.cpp QualityStats.cpp TileStats.cpp Reference.cpp utilities.cpp)
target_link_libraries(DataStats skewer_matrix)
target_link_libraries(DataStats nlopt)
target_link_libraries(DataStats ${SEQAN_LIBRARIES})
target_link_libraries(DataStats ${Boost_LIBRARIES})
target_compile_options(DataStats PRIVATE ${PRIVATE_COMPILE_OPTIONS})

add_library(DataStatsInterface STATIC DataStatsInterface.cpp)
target_link_libraries(DataStatsInterface DataStats)
target_compile_options(DataStatsInterface PRIVATE ${PRIVATE_COMPILE_OPTIONS})

# Executables
add_executable(reseq main.cpp ProbabilityEstimates.cpp Simulator.cpp AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp QualityStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp SimulatorTest.cpp TileStatsTest.cpp utilitiesTest.cpp VectTest.cpp)
target_link_libraries(reseq DataStats)
target_link_libraries(reseq gtest)
#target_link_libraries(reseq nlopt)
target_link_libraries(reseq ${SEQAN_LIBRARIES})
target_link_libraries(reseq ${Boost_LIBRARIES})
target_compile_options(reseq PRIVATE ${PRIVATE_COMPILE_OPTIONS})
