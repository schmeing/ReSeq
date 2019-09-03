# Own libraries
add_library(skewer_fastq src/fastq.cpp)

add_library(skewer_matrix src/matrix.cpp)
target_link_libraries(skewer_matrix skewer_fastq)