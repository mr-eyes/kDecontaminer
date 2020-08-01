cmake_minimum_required (VERSION 3.4)
project (kDecontaminer C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")


# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(libs/kProcessor)
set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/libs/kProcessor/include")

# Add include directories.
include_directories("${PROJECT_SOURCE_DIR}/libs/kProcessor/ThirdParty/MQF/include")
include_directories (${kProcessor_INCLUDE_PATH})


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------


include_directories(libs/kseq/)
include_directories(includes)
include_directories(libs/sqlite3pp/)


add_executable (peReadsStats reads_stats.cpp src/kseqReader.cpp)
target_link_libraries (peReadsStats kProcessor z)
target_include_directories(peReadsStats INTERFACE ${kProcessor_INCLUDE_PATH})

#add_executable (analyzeIdx analyze_index.cpp)
#target_link_libraries (analyzeIdx kProcessor z)
#target_include_directories(analyzeIdx INTERFACE ${kProcessor_INCLUDE_PATH})

#add_executable (testKseq test_kseq.cpp)
#target_link_libraries (testKseq z)



# genesIndexing executable (Have an issue after some time in the execution)
#add_executable (genesIndexing genes_indexing.cpp)
#target_link_libraries (genesIndexing kProcessor pthread z)
#target_include_directories(genesIndexing INTERFACE ${kProcessor_INCLUDE_PATH})