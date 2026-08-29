cmake_minimum_required(VERSION 3.4)
project(coco_oracle)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -w")
# KSeqWrapper is not linked in: the oracle never reads a sequence file, and
# its no-zlib fallback branch references a symbol that does not exist.

# The CoCo checkout this directory lives in *is* the reference; there is no
# copy of the C++ sources under Coco-rust.
set(REF ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(${REF}/lib/kseq ${REF}/lib/kmer-cnt ${REF}/src)

SET(GATB_CORE_EXCLUDE_TOOLS 1)
SET(GATB_CORE_EXCLUDE_TESTS 1)
SET(GATB_CORE_EXCLUDE_EXAMPLES 1)
add_subdirectory(${REF}/lib/gatb-core/gatb-core ${CMAKE_CURRENT_BINARY_DIR}/gatb)
add_definitions(${gatb-core-flags})
include_directories(${gatb-core-includes})
link_directories(${gatb-core-extra-libraries-path})

# CountProfile.cpp is #included by oracle_main.cpp, so it is not listed here.
add_executable(oracle
  oracle_main.cpp
  ${REF}/src/types.cpp
  ${REF}/src/kmer.cpp
  ${REF}/src/uint128.cpp
  ${REF}/src/KmerTranslator.cpp
  ${REF}/src/Info.cpp
  ${REF}/src/filehandling.cpp
)
target_link_libraries(oracle ${gatb-core-libraries})
