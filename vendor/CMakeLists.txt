# Eigen 3.3.3
add_library(distort-eigen-vendored INTERFACE)
target_include_directories(distort-eigen-vendored INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
target_compile_definitions(distort-eigen-vendored INTERFACE EIGEN_MPL2_ONLY)
