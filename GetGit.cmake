FIND_PACKAGE(Git)
IF(GIT_FOUND)
  #EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags RESULT_VARIABLE res_var OUTPUT_VARIABLE GIT_COM_ID )
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD RESULT_VARIABLE res_var OUTPUT_VARIABLE GIT_COM_ID )
  IF( NOT ${res_var} EQUAL 0 )
    SET( GIT_COMMIT_ID "git commit id unknown")
    MESSAGE( WARNING "Git failed (not a repo, or no tags). Build will not contain git revision info." )
  ENDIF()
  STRING( REPLACE "\n" "" GIT_COMMIT_ID ${GIT_COM_ID} )
ELSE()
  SET( GIT_COMMIT_ID "unknown (git not found!)")
  MESSAGE( WARNING "Git not found. Build will not contain git revision info." ) 
ENDIF()

SET( vstring "#include <SvnRevision.h>\n\n"
             "//SvnRevision.cpp - written by cmake. changes will be lost!\n"
             "const std::string SvnRevision::revisionNumber = \"${GIT_COMMIT_ID}\"\;\n")

FILE(WRITE SvnRevision.cpp.txt ${vstring} )
# copy the file to the final header only if the version changes
# reduces needless rebuilds
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        SvnRevision.cpp.txt ${CMAKE_CURRENT_SOURCE_DIR}/../src/SvnRevision.cpp)
