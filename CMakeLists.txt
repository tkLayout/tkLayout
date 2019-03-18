#############################################################################
#
# CMAKE file for building tkLayout software
#
#  Author: Z. Drasal, CERN
#
#############################################################################

#----------------------------------------------------------------------------
# CMAKE compatibility issues
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

#----------------------------------------------------------------------------
# Project name and version
PROJECT(tklayout)

#----------------------------------------------------------------------------
# Project options
OPTION(INSTALL_DOC "Set to ON to build/install Documentation" OFF )


# using cache lets users override with " cmake -DCMAKE_CXX_STANDARD=14 .. "
SET(CMAKE_CXX_STANDARD 17 CACHE STRING "")

# Output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )


#----------------------------------------------------------------------------
# Project dependencies

FIND_PACKAGE( Boost COMPONENTS program_options filesystem system REQUIRED )

find_package(ROOT COMPONENTS Geom RIO HistPainter)
include_directories(${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Include, source, dirs 
#




INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR}
                     ${ROOT_INCLUDE_DIRS} 
                     ${PROJECT_SOURCE_DIR}/include/
                     ${PROJECT_SOURCE_DIR}/include/AnalyzerVisitors
                     ${PROJECT_SOURCE_DIR}/include/InnerCabling
                     ${PROJECT_SOURCE_DIR}/include/OuterCabling
                     )

FILE( GLOB all_sources ${PROJECT_SOURCE_DIR}/src/*.cc 
                       ${PROJECT_SOURCE_DIR}/src/AnalyzerVisitors/*.cc 
                       ${PROJECT_SOURCE_DIR}/src/InnerCabling/*.cc 
                       ${PROJECT_SOURCE_DIR}/src/OuterCabling/*.cc)
FILE( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh 
                   ${PROJECT_SOURCE_DIR}/include/AnalyzerVisitors*.hh 
                   ${PROJECT_SOURCE_DIR}/include/InnerCabling*.hh 
                   ${PROJECT_SOURCE_DIR}/include/OuterCabling*.hh)

# Exlude files with main function defined
SET ( sources "" )
FOREACH( file ${all_sources} )
 IF ( ${file} MATCHES "MaterialSection.cc" OR ${file} MATCHES "HoughTrack.cc" OR
      ${file} MATCHES "TrackShooter.cc" OR ${file} MATCHES "tunePtParam.cc" OR ${file} MATCHES "diskPlace.cc")
   SET ( APPEND source_other ${file} )
   MESSAGE( STATUS "Omitting the following ?buggy? file: ${file} !!!" ) 
 ELSEIF( ${file} MATCHES "tklayout.cc" OR ${file} MATCHES "setup.cc" OR ${file} MATCHES "delphize.cc" )
   IF ( ${file} MATCHES "tklayout.cc" ) 
     SET( source_tklayout ${file} )
   ENDIF()
   IF ( ${file} MATCHES "setup.cc" ) 
     SET( source_setup ${file} )
   ENDIF()
   IF ( ${file} MATCHES "delphize.cc" )
     SET( source_delphize ${file} )
   ENDIF()
 ELSE()
   IF( ${file} MATCHES "MainConfigHandler.cc" )
     SET( source_mainhandler ${file} )
   ENDIF()
   IF( ${file} MATCHES "global_funcs.cc" )
     SET( source_globalfunctions ${file} )
   ENDIF()
   IF( ${file} MATCHES "GraphVizCreator.cc" )
     SET( source_graphvizcreator ${file} )
   ENDIF()
   LIST( APPEND sources ${file} )
 ENDIF()
ENDFOREACH()

#----------------------------------------------------------------------------
# Set revision
# a custom target that is always built
ADD_CUSTOM_TARGET(revisiontag ALL)

# creates svnversion.hh using cmake script
ADD_CUSTOM_COMMAND(TARGET revisiontag COMMAND ${CMAKE_COMMAND}
   -DSOURCE_DIR=${CMAKE_SOURCE_DIR} 
   -P ${CMAKE_SOURCE_DIR}/GetGitVersion.cmake)

#--------------------------------------------------------
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)            # Set automatically rpath for dynamic linking of external libraries
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") # Set automatically rpath for dynamic linking of internally created & installed libraries

ADD_EXECUTABLE(tklayout ${source_tklayout} ${sources} ${headers} )
ADD_EXECUTABLE(setup.bin ${source_setup} ${source_graphvizcreator} ${source_mainhandler} ${source_globalfunctions} ${headers} )
ADD_EXECUTABLE(delphize ${source_delphize} )

# explicitly say that the executable depends on custom target
ADD_DEPENDENCIES(tklayout revisiontag)

TARGET_LINK_LIBRARIES(tklayout ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(setup.bin ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(delphize ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
INSTALL(TARGETS tklayout  RUNTIME DESTINATION bin)
INSTALL(TARGETS setup.bin RUNTIME DESTINATION bin)
INSTALL(TARGETS delphize  RUNTIME DESTINATION bin)

INSTALL(DIRECTORY style DESTINATION share/tkLayout)
INSTALL(DIRECTORY xml DESTINATION share/tkLayout)
INSTALL(DIRECTORY geometries DESTINATION share/tkLayout)
INSTALL(DIRECTORY config DESTINATION share/tkLayout)

IF(CMAKE_HOST_UNIX)
    
    # Set the default version of SvnRevision.cc file
    INSTALL( CODE "MESSAGE(STATUS \"Setting SvnRevision.cc to default\" )
                   EXECUTE_PROCESS(COMMAND bash -c \"if [ -e SvnRevision.orig.cc ]; then rm -f ${PROJECT_SOURCE_DIR}/src/SvnRevision.cc; mv SvnRevision.orig.cc ${PROJECT_SOURCE_DIR}/src/SvnRevision.cc; fi\")" )
    

    # Configure tkLayout
    #INSTALL( CODE "MESSAGE(STATUS \"Configuring tkLayout software...\" )
    #               MESSAGE(STATUS \" Target directories setup:\" )
    #		   EXECUTE_PROCESS(COMMAND bash -c \"if [ -e ${PROJECT_BINARY_DIR}/bin/setup.bin ]; then ${PROJECT_BINARY_DIR}/bin/setup.bin --dirNames; fi\")" ) 
    #INSTALL( CODE "MESSAGE(STATUS \" Installation status:\" )
    #               EXECUTE_PROCESS(COMMAND bash -c \"if ! ${PROJECT_SOURCE_DIR}/bin/setup.bin --checkDir; then echo \\\"Problem during installation\\\"; else echo \\\"Installation successful\\\"; fi\")" )

ENDIF()

#----------------------------------------------------------------------------
# Documentation
#
# Find Doxygen
IF( INSTALL_DOC )
FIND_PACKAGE( Doxygen )
IF( DOXYGEN_FOUND )

    # Add doc target
    ADD_CUSTOM_TARGET( doc ) 
    ADD_CUSTOM_COMMAND( TARGET doc
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${PROJECT_SOURCE_DIR}/README.md" "${PROJECT_SOURCE_DIR}/doc/README.md"    
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS "Doxyfile"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation in ${PROJECT_SOURCE_DIR}/doc... "
        VERBATIM )

ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# Install documentation
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${PROJECT_BINARY_DIR} ARGS doc)" )
    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE )
ENDIF()


#----------------------------------------------------------------------------
# Uninstall 
#
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    
