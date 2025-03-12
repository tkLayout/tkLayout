cmake_minimum_required(VERSION 3.10)
project(TkLayout)

# -------------------------------------------------------
# 1) C++ standard and compiler options
# -------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    -Wall
    -Werror
    -fmax-errors=2
    -Wno-error=deprecated-declarations
    -Wno-error=infinite-recursion
    -Wno-error=unused-result
    -Wno-error=uninitialized
    -Wno-error=maybe-uninitialized
    -O3
)

# If you absolutely need -O5, note that typical compilers treat -O3 as maximum
# or have -Ofast. Modify as needed:
# add_compile_options(-Ofast)

# Other flags from Makefile can be appended here if you like:
#   -g, -pg, etc.
# -------------------------------------------------------
# 2) Find required packages (ROOT, Boost)
# -------------------------------------------------------
# Root
find_package(ROOT REQUIRED COMPONENTS Hist RIO Core HistPainter Geom)
# => 'ROOT_INCLUDE_DIRS' and 'ROOT_LIBRARIES' become available

# Boost
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
# => 'Boost_INCLUDE_DIRS' and 'Boost_LIBRARIES' become available

# -------------------------------------------------------
# 3) Include directories
# -------------------------------------------------------
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# -------------------------------------------------------
# 4) (Optional) Capture revision macro from getRevisionDefine script
# -------------------------------------------------------
execute_process(
    COMMAND ./getRevisionDefine
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE SVN_REV_DEFINES
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# -------------------------------------------------------
# 5) Gather source files
# -------------------------------------------------------
# Adjust paths/filenames if they differ in your project
set(SRC_FILES
    src/Analyzer.cc
    src/AnalyzerTools.cc
    src/AnalyzerVisitor.cc
    src/Bag.cc
    src/Barrel.cc
    src/capabilities.cc
    src/ConversionStation.cc
    src/CoordinateOperations.cc
    src/DetectorModule.cc
    src/DetIdBuilder.cc
    src/Disk.cc
    src/Endcap.cc
    src/Extractor.cc
    src/GeometricModule.cc
    src/global_funcs.cc
    src/GraphVizCreator.cc
    src/Histo.cc
    src/Hit.cc
    src/InactiveElement.cc
    src/InactiveRing.cc
    src/InactiveSurfaces.cc
    src/InactiveTube.cc
    src/IrradiationMap.cc
    src/IrradiationMapsManager.cc
    src/Layer.cc
    src/MainConfigHandler.cc
    src/MatCalc.cc
    src/MatCalcDummy.cc
    src/MaterialBudget.cc
    src/MaterialObject.cc
    src/MaterialProperties.cc
    src/MaterialTab.cc
    src/MaterialTable.cc
    src/Materialway.cc
    src/MatParser.cc
    src/MessageLogger.cc
    src/ModuleCap.cc
    src/Module.cc
    src/Palette.cc
    src/PlotDrawer.cc
    src/Polygon3d.cc
    src/Property.cc
    src/PtError.cc
    src/PtErrorAdapter.cc
    src/ReportIrradiation.cc
    src/ReportModuleCount.cc
    src/Ring.cc
    src/RodPair.cc
    src/RootWeb.cc
    src/Sensor.cc
    src/SimParms.cc
    src/Squid.cc
    src/StopWatch.cc
    src/SummaryTable.cc
    src/SupportStructure.cc
    src/tk2CMSSW.cc
    src/Tracker.cc
    src/Track.cc
    src/Usher.cc
    src/Vizard.cc
    src/WeightDistributionGrid.cc
    src/XMLWriter.cc
    src/SvnRevision.cc
)

# Analyzer Visitors
set(SRC_ANALYZER_VISITORS
    src/AnalyzerVisitors/Bandwidth.cc
    src/AnalyzerVisitors/GeometricInfo.cc
    src/AnalyzerVisitors/IrradiationPower.cc
    src/AnalyzerVisitors/MaterialBillAnalyzer.cc
    src/AnalyzerVisitors/TriggerDistanceTuningPlots.cc
    src/AnalyzerVisitors/TriggerFrequency.cc
    src/AnalyzerVisitors/TriggerProcessorBandwidth.cc
    src/AnalyzerVisitors/ModuleCount.cc
)

# OuterCabling
set(SRC_OUTER_CABLING
    src/OuterCabling/OuterBundle.cc
    src/OuterCabling/OuterCable.cc
    src/OuterCabling/outer_cabling_constants.cc
    src/OuterCabling/outer_cabling_functions.cc
    src/OuterCabling/OuterCablingMap.cc
    src/OuterCabling/OuterDTC.cc
    src/OuterCabling/OuterGBT.cc
    src/OuterCabling/ModulesToBundlesConnector.cc
    src/OuterCabling/PhiPosition.cc
    src/OuterCabling/ServicesChannel.cc
)

# InnerCabling
set(SRC_INNER_CABLING
    src/InnerCabling/GBT.cc
    src/InnerCabling/HvLine.cc
    src/InnerCabling/inner_cabling_constants.cc
    src/InnerCabling/inner_cabling_functions.cc
    src/InnerCabling/InnerBundle.cc
    src/InnerCabling/InnerDTC.cc
    src/InnerCabling/InnerCablingMap.cc
    src/InnerCabling/ModulesToPowerChainsConnector.cc
    src/InnerCabling/PowerChain.cc
)

# -------------------------------------------------------
# 6) Define the library that holds most of your code
# -------------------------------------------------------
add_library(tklayout_core
    ${SRC_FILES}
    ${SRC_ANALYZER_VISITORS}
    ${SRC_OUTER_CABLING}
    ${SRC_INNER_CABLING}
)

# Apply the revision macro (from getRevisionDefine) to this library build
target_compile_options(tklayout_core PRIVATE ${SVN_REV_DEFINES})

# Link it with ROOT & Boost
target_link_libraries(tklayout_core
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

# -------------------------------------------------------
# 7) Define your main executables
# -------------------------------------------------------
# A) tklayout
add_executable(tklayout src/tklayout.cc)
target_link_libraries(tklayout
    tklayout_core
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)
set_target_properties(tklayout PROPERTIES
    OUTPUT_NAME "tklayout"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# B) setup
add_executable(setup_exe src/setup.cc)
target_link_libraries(setup_exe
    tklayout_core
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)
set_target_properties(setup_exe PROPERTIES
    OUTPUT_NAME "setup"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# C) diskPlace
add_executable(diskPlace src/diskPlace.cc)
target_link_libraries(diskPlace
)
set_target_properties(diskPlace PROPERTIES
    OUTPUT_NAME "diskPlace"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# -------------------------------------------------------
# 8) Doxygen docs
# -------------------------------------------------------
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/tkdoc.doxy
                   ${CMAKE_CURRENT_BINARY_DIR}/tkdoc.doxy
                   COPYONLY)
    add_custom_target(doxydoc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/tkdoc.doxy
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen..."
        VERBATIM
    )
endif()

# -------------------------------------------------------
# 9) "clean_all" target to mimic 'make clean'
# -------------------------------------------------------
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMENT "Remove the entire build directory."
)
