cmake_minimum_required(VERSION 3.14)

project(qemdilib)

# download/ build MDI library
# always build static library with -fpic
# support cross-compilation and ninja-build
include(ExternalProject)
ExternalProject_Add(mdi_build
  URL     "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.4.1.tar.gz"
  URL_MD5 "f9505fccd4c79301a619f6452dad4ad9"
  CMAKE_ARGS ${CMAKE_REQUEST_PIC}
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
  -Dlanguage=Fortran
  -Dlibtype=SHARED
  -Dpython_plugins=OFF
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS "<BINARY_DIR>/MDI_Library/libmdi.a"
  )

# where is the compiled library?
ExternalProject_get_property(mdi_build BINARY_DIR)
set(MDI_BINARY_DIR "${BINARY_DIR}/MDI_Library")
# workaround for older CMake versions
file(MAKE_DIRECTORY ${MDI_BINARY_DIR})

# create imported target for the MDI library
add_library(mdi-lib UNKNOWN IMPORTED)
add_dependencies(mdi-lib mdi_build)
set_target_properties(mdi-lib PROPERTIES
  IMPORTED_LOCATION "${MDI_BINARY_DIR}/libmdi.a"
  INTERFACE_INCLUDE_DIRECTORIES ${MDI_BINARY_DIR}
  )
