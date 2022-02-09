#link_directories( ${CMAKE_BINARY_DIR} )
#include_directories( ${CMAKE_BINARY_DIR} )



# Compile the driver

add_executable(driver_plug_c
               driver_plug_c.c)
target_link_libraries(driver_plug_c PRIVATE mdi)
set_target_properties(driver_plug_c PROPERTIES
    LINKER_LANGUAGE C
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

link_against_mpi_c(driver_plug_c)
