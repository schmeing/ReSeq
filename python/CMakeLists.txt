find_program(PYTHON "python")

if (PYTHON)
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	# Python packages
	FIND_PACKAGE(SWIG 3 REQUIRED)
	include(${SWIG_USE_FILE})
	FIND_PACKAGE(PythonLibs 3 REQUIRED)
	include_directories(${PYTHON_INCLUDE_PATH})
	
	include_directories("${PROJECT_SOURCE_DIR}/python")
	
	SET(CMAKE_SWIG_FLAGS "")
	
	# Set a variable where to put the python modules
	SET( PYTHON_MODULE_DIRECTORY  "${PROJECT_BINARY_DIR}/pyMods" )
	
	# Define the swig outdir to be the one from the previously set variable (.cxx and .py file, don't know how to define the outdir only for the .py file)
	set(CMAKE_SWIG_OUTDIR ${PYTHON_MODULE_DIRECTORY})
	
	SET_SOURCE_FILES_PROPERTIES(DataStats.i PROPERTIES CPLUSPLUS ON)
	SET_SOURCE_FILES_PROPERTIES(DataStats.i PROPERTIES SWIG_FLAGS "-Wall")
	SET_SOURCE_FILES_PROPERTIES(DataStats.i PROPERTIES SWIG_FLAGS "-modern")
	SET_SOURCE_FILES_PROPERTIES(DataStats.i PROPERTIES SWIG_FLAGS "-DSWIGWORDSIZE64")
	SWIG_ADD_MODULE(DataStats python DataStats.i)
	set_target_properties(_DataStats PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYTHON_MODULE_DIRECTORY}) # Define the outdir for the .so python module file
	SWIG_LINK_LIBRARIES(DataStats ${PYTHON_LIBRARIES})
	SWIG_LINK_LIBRARIES(DataStats ${Boost_LIBRARIES})
	SWIG_LINK_LIBRARIES(DataStats ${SEQAN_LIBRARIES})
	SWIG_LINK_LIBRARIES(DataStats DataStatsInterface)

	INSTALL(TARGETS _DataStats DESTINATION lib) 
	install(FILES ${PROJECT_BINARY_DIR}/pyMods/DataStats.py DESTINATION lib)

	INSTALL(FILES ${PROJECT_SOURCE_DIR}/python/plotDataStats.py ${PROJECT_SOURCE_DIR}/python/reseq-prepare-names.py
		DESTINATION bin
		PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
		
	# Solve issue that the wrapper file created by swig is not recreated if the included headers change
	# https://stackoverflow.com/questions/39358420/cmake-swig-dependicies
	
	# Removes generated file (if any of the dependent files are changed)
	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp
		COMMAND ${CMAKE_COMMAND} -E remove ${swig_generated_file_fullname}
		COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp
		DEPENDS ../reseq/DataStatsInterface.h ../reseq/Reference.h ../reseq/Surrounding.h ../reseq/utilities.hpp ../reseq/Vect.hpp  # The dependent files
		COMMENT "Removing old SWIG generated file" VERBATIM)

	# Custom target for establishing dependency
	add_custom_target(
		swigtrick
		DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp)

	add_dependencies(_DataStats swigtrick)
endif()