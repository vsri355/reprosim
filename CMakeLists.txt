

cmake_minimum_required(VERSION 3.20)
project(Reprosim VERSION 1.0.0 LANGUAGES Fortran C CXX)

# Policy settings
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif ()
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(IS_MULTI_CONFIG)
  set(CONFIG_DIR_SUFFIX "/$<CONFIG>")
endif()

if (MSVC)
  cmake_path(GET CMAKE_Fortran_COMPILER PARENT_PATH FORTRAN_RUNTIME_PATH)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Required for creating documentation from restructured text files
find_package(Sphinx)

# Helpful functions
include(miscfunctions)

# Create additional build types
include(buildtypes)

# Mark some variables as advanced hiding them from basic users
tidy_gui()

# Build shared or static library
if(NOT BUILD_SHARED)
  set(BUILD_SHARED ON)
endif()
set(REPROSIM_BUILD_SHARED ${BUILD_SHARED} CACHE BOOL "Build shared library (true) or static library (false)")
set(BUILD_SHARED_LIBS ${REPROSIM_BUILD_SHARED} CACHE INTERNAL "Internalise BUILD_SHARED_LIBS, manipulate via REPROSIM_BUILD_SHARED" FORCE)

if(IS_MULTI_CONFIG)
  mark_as_advanced(CMAKE_CONFIGURATION_TYPES)
else()
  # Set the build type Release, Debug plus more if you define them
  if(NOT BUILD_TYPE)
    set(BUILD_TYPE Release)
  endif()
  set(REPROSIM_BUILD_TYPE ${BUILD_TYPE} CACHE STRING "Set the build type; [Release], Debug, Pedantic (GNU Fortran only)")
  set(CMAKE_BUILD_TYPE ${REPROSIM_BUILD_TYPE} CACHE INTERNAL "Internalise CMAKE_BUILD_TYPE, manipulate via REPROSIM_BUILD_TYPE" FORCE)
endif()

# Set whether to build bindgins or not, may require SWIG
option(REPROSIM_BUILD_BINDINGS "Build bindings for ${PROJECT}, may require SWIG" YES)

# Set base install location
if(NOT INSTALL_PREFIX)
  set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(REPROSIM_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE STRING  "Set the base install prefix")
set(CMAKE_INSTALL_PREFIX ${REPROSIM_INSTALL_PREFIX} CACHE INTERNAL "Internlaise CMAKE_INSTALL_PREFIX, manipulate via REPROSIM_INSTALL_PREFIX" FORCE)
mark_as_advanced(REPROSIM_INSTALL_PREFIX)

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# All the source files for the library and bindings are located under this directory.
add_subdirectory(src)

# Must enable testing at the lowest level in CMake.
enable_testing()

# All the tests for the library are located under this directory.
add_subdirectory(tests)

set(ADDITIONAL_GENERATED_FILES)
if(SPHINX_FOUND)
  if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
  endif()

  if(NOT DEFINED SPHINX_THEME_DIR)
    set(SPHINX_THEME_DIR)
  endif()

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  # Sphinx source directory
  set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/documentation")

  add_custom_target(docs
    ${SPHINX_EXECUTABLE}
      -q -b html
      -d "${SPHINX_CACHE_DIR}"
      "${SPHINX_SOURCE_DIR}"
      "${SPHINX_HTML_DIR}"
    WORKING_DIRECTORY "${SPHINX_SOURCE_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

  list(APPEND ADDITIONAL_GENERATED_FILES html _doctrees)

endif()

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_GENERATED_FILES}")

