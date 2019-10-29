find_program(PYTHON "python")

if (PYTHON)
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	# Python packages
	FIND_PACKAGE(SWIG 3 REQUIRED)
	include(${SWIG_USE_FILE})
	message(${SWIG_USE_FILE})
	FIND_PACKAGE(PythonLibs 2 EXACT REQUIRED)
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
endif()
