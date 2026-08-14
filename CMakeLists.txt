cmake_minimum_required(VERSION 3.16)
project(xring_proxy_master)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT DEFINED ENV{GUROBI_HOME} AND NOT DEFINED GUROBI_HOME)
  message(FATAL_ERROR "Set GUROBI_HOME (CMake var or env) to your Gurobi installation.")
endif()
if(DEFINED ENV{GUROBI_HOME})
  set(GUROBI_HOME $ENV{GUROBI_HOME})
endif()

include_directories(${GUROBI_HOME}/include)
link_directories(${GUROBI_HOME}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src)

set(COMMON_SOURCES
        src/DemandMatrices.cpp
        src/InstanceSetup.cpp
        src/MILPSolver.cpp
        src/RingLayout.cpp
        src/RingDemandRouting.cpp
        src/RingGeometry.cpp
        src/LShapeCrossing.cpp
        src/ProxyMasterMilp.cpp
        src/StageE.cpp
        src/ProxyMasterLoop.cpp
        src/ShortcutGrid.cpp
        src/ShortcutRouter.cpp
        src/ShortcutOrchestrator.cpp
        src/ShortcutMethodD.cpp
        src/ShortcutMethods.cpp
        src/WavelengthLoadBalance.cpp
        src/PlotExport.cpp
        src/ShortcutExport.cpp
)

add_executable(proxy_master
        src/main.cpp
        ${COMMON_SOURCES}
)
target_link_libraries(proxy_master PRIVATE gurobi_c++ gurobi130)

if(EXISTS "${CMAKE_SOURCE_DIR}/tools/snr_crossing_benchmark.cpp")
  add_executable(snr_crossing_benchmark
          tools/snr_crossing_benchmark.cpp
          ${COMMON_SOURCES}
  )
  target_link_libraries(snr_crossing_benchmark PRIVATE gurobi_c++ gurobi130)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/tools/benchmark_batch.cpp")
  add_executable(benchmark_batch
          tools/benchmark_batch.cpp
          ${COMMON_SOURCES}
  )
  target_link_libraries(benchmark_batch PRIVATE gurobi_c++ gurobi130)
endif()
