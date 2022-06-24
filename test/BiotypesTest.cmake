include(FetchContent)

find_package(Catch2 3 QUIET)
if (NOT Catch2_FOUND)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.0.1)

  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

add_executable(biotypes_test ${CMAKE_CURRENT_LIST_DIR}/bioseq_test.cc)
target_link_libraries(biotypes_test
        PRIVATE
        biotypes
        Catch2::Catch2WithMain)
