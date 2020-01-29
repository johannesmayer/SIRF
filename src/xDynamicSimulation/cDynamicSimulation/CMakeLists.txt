#========================================================================
# Author: Johannes Mayer
# Date: 15.03.2018
#=========================================================================



add_definitions(-DBOOST_ALL_NO_LIB)

if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
    set(Boost_DLL_DIR ${Boost_LIBRARY_DIR_RELEASE})
    message(STATUS "Install boost shared libraries from ${Boost_DLL_DIR} ")
	foreach (__boost_lib system filesystem thread date_time chrono)
	   file(GLOB Boost_DLL "${Boost_DLL_DIR}/boost_${__boost_lib}*.dll")
       install( FILES ${Boost_DLL} DESTINATION bin )
	endforeach()
  endif()

# Should be using target_link_libraries
# include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
# include_directories(${PROJECT_SOURCE_DIR}/src/xGadgetron/cGadgetron)
# include_directories(${PROJECT_SOURCE_DIR}/src/xSTIR/cSTIR)
# include_directories(${PROJECT_SOURCE_DIR}/src/Registration/cReg)
# include_directories(${PROJECT_SOURCE_DIR}/src/iUtilities)
# include_directories("${GADGETRON_INCLUDE_DIR}/../")

include_directories(./Testing)

set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynamicsimulation_x.cpp  Testing/tests_dynamicsimulation.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_mr_acquisition_model.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} auxiliary_input_output.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_auxiliary_testing_functions.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_auxiliary_input_output.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/auxiliary_testing_functions.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/all_simulation_tests.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissueparameters.cpp Testing/tests_tissueparameters.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} tissuelabelmapper.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} contrastgenerator.cpp  Testing/tests_contrastgenerator.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} phantom_input.cpp Testing/tests_phantom_input.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_encoding.cpp )
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynamics.cpp Testing/tests_dynamics.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynsim_noisegenerator.cpp Testing/tests_noisegenerator.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} dynsim_deformer.cpp Testing/tests_dynsim_deformer.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} volume_orientator.cpp Testing/tests_volume_orientator.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/tests_memory_usage.cpp)
set( SRC_FILES_TO_COMPILE ${SRC_FILES_TO_COMPILE} Testing/app_mracquisitiondata.cpp)


add_library( cdynamicsimulation ${SRC_FILES_TO_COMPILE})

set (cdynamicsimulation "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope

TARGET_INCLUDE_DIRECTORIES(cdynamicsimulation PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)
target_include_directories(cdynamicsimulation PRIVATE "${FFTW3_INCLUDE_DIR}")
target_include_directories(cdynamicsimulation PUBLIC "${HDF5_INCLUDE_DIRS}")

target_link_libraries(cdynamicsimulation PUBLIC iutilities)
# Add boost library dependencies
if((CMAKE_VERSION VERSION_LESS 3.5.0) OR (NOT _Boost_IMPORTED_TARGETS))
  # This is harder than it should be on older CMake versions to be able to cope with
  # spaces in filenames.
  foreach(C SYSTEM FILESYSTEM THREAD DATE_TIME CHRONO)
    target_link_libraries(cdynamicsimulation PUBLIC optimized "${Boost_${C}_LIBRARY_RELEASE}")
    target_link_libraries(cdynamicsimulation PUBLIC debug  "${Boost_${C}_LIBRARY_DEBUG}")
  endforeach()
else()
  # Nice and simple for recent CMake (which knows about your Boost version)
  target_link_libraries(cdynamicsimulation PUBLIC Boost::system Boost::filesystem Boost::thread Boost::date_time Boost::chrono)
endif()

target_link_libraries(cdynamicsimulation PUBLIC cgadgetron)
target_link_libraries(cdynamicsimulation PUBLIC cstir)
target_link_libraries(cdynamicsimulation PUBLIC Reg)

message(STATUS "XXXXXXXXXXXXXXXXXXSTIR_LIB ${STIR_LIBRARIES}")
add_executable( TestJobs main.cpp ${STIR_REGISTRIES})
target_link_libraries( TestJobs cdynamicsimulation ${STIR_LIBRARIES} ${HDF5_LIBRARIES} ${ismrmrd_full_path} ${Gadge_libs_full_path})
ADD_TEST(NAME SIMU_TEST_CPLUSPLUS COMMAND TestJobs WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
