cmake_minimum_required(VERSION 2.6)

set ( snappy_src ../lib/snappy/snappy.cc ../lib/snappy/snappy-stubs-internal.cc ../lib/snappy/snappy-sinksource.cc )
add_executable(mtrap ${snappy_src} mtrap.cpp estimatefamily.cpp distance.cpp scorematrix.cpp tqmatrix.cpp weightedmatrix.cpp pairwise.cpp tree.cpp utility.cpp entropy.cpp primarylibrary.cpp search.cpp sequences.cpp globaloption.cpp probalign.cpp)
add_executable(test_tree ${snappy_src} test_tree.cpp tree.cpp utility.cpp scorematrix.cpp tqmatrix.cpp weightedmatrix.cpp)

# LD_FLAG
# set_target_properties(mtrap PROPERTIES LINK_FLAGS -lc)

# compiler option
# set_target_properties(mtrap PROPERTIES COMPILE_FLAGS -Wall)

# compress library
#target_link_libraries(mtrap snappy)

