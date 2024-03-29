file(GLOB files_decomp decomp_2d.f90
	               log.f90
                       glassman.f90)
include_directories(${CMAKE_SOURCE_DIR}/decomp2d)

#if(${FFT_Choice} MATCHES "generic")
#  file(GLOB files_fft fft_generic.f90)
#endif(${FFT_Choice} MATCHES "generic")

if(FFTW_FOUND)
  message (STATUS "Compiling using FFTW3")
  file(GLOB files_fft fft_fftw3.f90)
else(FFTW_FOUND)
  message (STATUS "Compiling using Generic FFT")
  file(GLOB files_fft fft_generic.f90)
endif(FFTW_FOUND)


set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
if (MPI_FOUND)
  target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(FFTW_FOUND)
  target_include_directories(decomp2d PRIVATE ${FFTW_INCLUDE_DIRS})
  #target_include_directories(decomp2d PRIVATE ${PKG_FFTW_LIBRARY_DIRS})
endif(FFTW_FOUND)

install(TARGETS decomp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
)
