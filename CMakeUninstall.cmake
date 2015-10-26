#############################################################################
#
# CMAKE uninstall module
#
#############################################################################

# Check that software was installed
IF( NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )
  MESSAGE( FATAL_ERROR "Cannot find install manifest: @CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )
ENDIF( NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )

# Read install_manifest.txt and uninstall all what was automatically installed
FILE( READ "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" files )

STRING( REGEX REPLACE "\n" ";" files "${files}" )
FOREACH( file ${files} )
  MESSAGE( STATUS "Uninstalling $ENV{DESTDIR}${file}" )
  IF( IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}" )
    EXEC_PROGRAM(
      "@CMAKE_COMMAND@" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF( NOT "${rm_retval}" STREQUAL 0)
      MESSAGE( FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}" )
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ELSE(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  ENDIF(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
ENDFOREACH( file )

# Remove symlink
MESSAGE( STATUS "Uninstalling symlink $ENV{HOME}/bin/tklayout" )
IF( IS_SYMLINK "$ENV{HOME}/bin/tklayout")
  EXEC_PROGRAM(
         "@CMAKE_COMMAND@" ARGS "-E remove \"$ENV{HOME}/bin/tklayout\""
         OUTPUT_VARIABLE rm_out
         RETURN_VALUE rm_retval
         )
  IF( NOT "${rm_retval}" STREQUAL 0)
      MESSAGE( FATAL_ERROR "Problem when removing symlink $ENV{HOME}/bin/tklayout" )
  ENDIF(NOT "${rm_retval}" STREQUAL 0)
ELSE()
MESSAGE(STATUS "Symlink $ENV{HOME}/bin/tklayout does not exist.")
ENDIF()


MESSAGE( STATUS "Uninstalling symlink $ENV{HOME}/bin/delphize" )
IF( IS_SYMLINK "$ENV{HOME}/bin/delphize")
  EXEC_PROGRAM(
         "@CMAKE_COMMAND@" ARGS "-E remove \"$ENV{HOME}/bin/delphize\""
         OUTPUT_VARIABLE rm_out
         RETURN_VALUE rm_retval
         )
  IF( NOT "${rm_retval}" STREQUAL 0)
      MESSAGE( FATAL_ERROR "Problem when removing symlink $ENV{HOME}/bin/delphize" )
  ENDIF(NOT "${rm_retval}" STREQUAL 0)
ELSE()
MESSAGE(STATUS "Symlink $ENV{HOME}/bin/delphize does not exist.")
ENDIF()

# Remove documentation
MESSAGE( STATUS "Removing generated doc html directory @CMAKE_CURRENT_SOURCE_DIR@/doc/html" )
IF( EXISTS "@CMAKE_CURRENT_SOURCE_DIR@/doc/html" )
  FILE( REMOVE_RECURSE "@CMAKE_CURRENT_SOURCE_DIR@/doc/html" )
ELSE ()
  MESSAGE(STATUS "Documentation html directory does not exist.")
ENDIF()
IF( EXISTS "@CMAKE_CURRENT_SOURCE_DIR@/doc/README.md")
  FILE( REMOVE "@CMAKE_CURRENT_SOURCE_DIR@/doc/README.md" )
ELSE ()
  MESSAGE(STATUS "Documentation file README.md was not copied here.")
ENDIF()
