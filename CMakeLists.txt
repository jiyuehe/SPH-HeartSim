cmake_minimum_required(VERSION 3.16)
project(SPHinXsys VERSION 1.0.0 LANGUAGES CXX)

set(SPHINXSYS_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake)
include(Common) # brings macro to current namespace

# ------ Options
option(SPHINXSYS_3D "Build sphinxsys_3d library" ON)
option(SPHINXSYS_BUILD_TESTS "Build sim" ON)

# ------ Global properties (Some cannot be set on INTERFACE targets)
set(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Enable verbose compilation commands for Makefile and Ninja" FORCE) # Extra fluff needed for Ninja: https://github.com/ninja-build/ninja/issues/900
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate `compile_commands.json` file which is necessary for code completioner, static analyzer, etc

# set(CMAKE_CXX_STANDARD_REQUIRED ON) # C++ standard requirement for all targets is NOT optional
set(CMAKE_CXX_EXTENSIONS OFF) # Disable any non-standard compiler specific C++ extensions for all targets.

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE) # Enable for all targets, needed for python modules

# ------ Declare core library
add_library(sphinxsys_core INTERFACE)

# ------ Local properties
target_compile_features(sphinxsys_core INTERFACE cxx_std_17) # Requires a compiler supporting at least C++17

# ## Simbody
find_package(Simbody CONFIG REQUIRED)
set(Simbody_LIBS
    $<$<BOOL:${Simbody_LIBRARIES}>:${Simbody_LIBRARIES}> # contains "SimTKcommon SimTKmath SimTKsimbody" in triplet with dynamic linking, empty otherwise
    $<$<BOOL:${Simbody_STATIC_LIBRARIES}>:${Simbody_STATIC_LIBRARIES}> # contains "SimTKcommon_static SimTKmath_static SimTKsimbody_static" in triplet with static linking, empty otherwise
)
target_link_libraries(sphinxsys_core INTERFACE ${Simbody_LIBS})
target_compile_definitions(sphinxsys_core INTERFACE _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # Because of Simbody headers

# ## Eigen
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE Eigen3::Eigen)

# ## TBB
find_package(TBB CONFIG REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE TBB::tbb TBB::tbbmalloc $<$<PLATFORM_ID:Windows>:TBB::tbbmalloc_proxy>)

# ## Threads
find_package(Threads REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE Threads::Threads)

# ## Boost
set(Boost_NO_WARN_NEW_VERSIONS TRUE) # In case your CMake version is older than the release of Boost found

# Boost geometry and program-options are old enough to afford skipping this warning
find_package(Boost REQUIRED)
find_path(BOOST_INCLUDE_DIR boost/geometry) # Header-only Boost libraries are not components

if(NOT BOOST_INCLUDE_DIR)
    message(FATAL_ERROR "Please install Boost.Geometry library")
endif()

target_link_libraries(sphinxsys_core INTERFACE Boost::boost)
find_package(Boost QUIET COMPONENTS program_options)

if(TARGET Boost::program_options)
    target_compile_definitions(sphinxsys_core INTERFACE BOOST_AVAILABLE)
    target_link_libraries(sphinxsys_core INTERFACE Boost::program_options)
endif()

# ------ Setup the concrete libraries
add_subdirectory(src)
#add_subdirectory(modules)

# ------ Setup the tests
if(SPHINXSYS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(sim)
endif()
