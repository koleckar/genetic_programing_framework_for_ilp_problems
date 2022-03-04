cmake_minimum_required(VERSION 3.20)

project(semestral_work)

set(CMAKE_CXX_STANDARD 17)

include_directories(src/problemSpecificBackend)
include_directories(src/genericSolverEA)
include_directories(src/ILP)

add_executable(main
        src/main.cpp
        src/problemSpecificBackend/VRPTWInstance.cpp
        src/problemSpecificBackend/VRPTWInstance.h
        src/problemSpecificBackend/CVRPInstance.cpp
        src/problemSpecificBackend/CVRPInstance.h
        src/problemSpecificBackend/VrpRepXmlReader.cpp
        src/problemSpecificBackend/VrpRepXmlReader.h
        src/genericSolverEA/GenericEASolver.cpp
        src/genericSolverEA/GenericEASolver.h
        src/genericSolverEA/Instance.h
        src/genericSolverEA/Population.h
        src/genericSolverEA/Population.cpp
        src/genericSolverEA/Solution.h
        lib/tinyxml/tinyxml.cpp
        lib/tinyxml/tinyxmlerror.cpp
        lib/tinyxml/tinyxmlparser.cpp
        lib/tinyxml/tinystr.cpp)

add_executable(test_code
        src/testCode.cpp
        src/problemSpecificBackend/VrpRepXmlReader.cpp
        src/problemSpecificBackend/VrpRepXmlReader.h
        lib/tinyxml/tinyxml.cpp
        lib/tinyxml/tinyxmlerror.cpp
        lib/tinyxml/tinyxmlparser.cpp
        lib/tinyxml/tinystr.cpp
        src/genericSolverEA/Population.cpp
        src/genericSolverEA/Population.h)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -march=native -pthread")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIR})

add_executable(ILP
        src/ILP/ILPModel.h
        src/ILP/ILPModel.cpp
        )



target_link_libraries(ILP
        ${GUROBI_LIBRARIES}
        )


