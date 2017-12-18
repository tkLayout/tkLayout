#############################################################################
#
# CMAKE module to find Git version
#
#  Author: Z. Drasal, CERN
#
#############################################################################
FIND_PACKAGE(Git)
IF(GIT_FOUND)

  # Get Git branch
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD RESULT_VARIABLE res_var_branch OUTPUT_VARIABLE GIT_BRANCH )
  IF( NOT ${res_var_branch} EQUAL 0 )
    SET( GIT_BRANCH "Git_Command_Failed")
    MESSAGE( WARNING "Git branch failed (not a repo, or no tags). Build will not contain git revision info." )
  ENDIF()
 
  # Get Git version
  IF(CMAKE_HOST_UNIX) 
    # Unix like system
#   EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-list HEAD | wc -l RESULT_VARIABLE res_var OUTPUT_VARIABLE GIT_VERSION )
    EXECUTE_PROCESS(COMMAND bash -c "git rev-list HEAD | wc -l" RESULT_VARIABLE res_var_version OUTPUT_VARIABLE GIT_VERSION )
    
  ELSE()
    # Non-unix system
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD RESULT_VARIABLE res_var_version OUTPUT_VARIABLE GIT_VERSION )
  ENDIF()
  IF( NOT ${res_var_version} EQUAL 0 )
    SET( GIT_VERSION "Git_Command_Failed")
    MESSAGE( WARNING "Git version failed (not shell, not a repo, or no tags). Build will not contain git revision info." )
  ENDIF()

  # Get Git author
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} config --get remote.origin.url RESULT_VARIABLE res_var_author OUTPUT_VARIABLE GIT_AUTHOR )
  IF( NOT ${res_var_author} EQUAL 0 )
    SET( GIT_AUTHOR "Git_Command_Failed")
    MESSAGE( WARNING "Git config url failed (not a repo, or no tags). Build will not contain git revision info." )
  ENDIF()

  STRING( REPLACE "\n" "" GIT_BRANCH  ${GIT_BRANCH}  )
  STRING( REPLACE "\n" "" GIT_VERSION ${GIT_VERSION} )
  STRING( REPLACE "\n" "" GIT_AUTHOR  ${GIT_AUTHOR}  )
 
  # All information found
  IF ( (${res_var_branch} EQUAL 0) AND (${res_var_version} EQUAL 0) AND (${res_var_author} EQUAL 0) ) 
    #STRING( CONCAT GIT_REVISION ${GIT_BRANCH} "-" ${GIT_VERSION} " on " ${GIT_AUTHOR} )
    SET( GIT_REVISION "${GIT_BRANCH}-${GIT_VERSION} on ${GIT_AUTHOR}" )

    MESSAGE( STATUS "")
    MESSAGE( STATUS "Git Revision: ${GIT_REVISION} " )
    MESSAGE( STATUS "" )

    # Unix like system
    IF(CMAKE_HOST_UNIX)
  
      # Create new SvnRevision.cc
      SET( vstring "#include <SvnRevision.hh>\n\n"
                   "//SvnRevision.cc - written by cmake. changes will be lost!\n"
                   "const std::string SvnRevision::revisionNumber = \"${GIT_REVISION}\"\;\n")

      FILE(WRITE SvnRevision.cc.txt ${vstring} )
      EXECUTE_PROCESS(COMMAND bash -c "if [ ! -e SvnRevision.orig.cc ]; then cp -p ${CMAKE_CURRENT_SOURCE_DIR}/../src/SvnRevision.cc SvnRevision.orig.cc; fi" )
      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
                              SvnRevision.cc.txt ${CMAKE_CURRENT_SOURCE_DIR}/../src/SvnRevision.cc)
    ENDIF()
  ENDIF()
ENDIF()


