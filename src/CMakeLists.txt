#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (fpgpu_sources "")
list (APPEND fpgpu_sources
    main.f90
)

add_executable (xfpgpu ${fpgpu_sources})
set_property(TARGET xfpgpu PROPERTY LINKER_LANGUAGE Fortran)