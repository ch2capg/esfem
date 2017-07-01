cmake_minimum_required(VERSION 3.0)
project(diss_surface VERSION 1.0)

set(HOME $ENV{HOME})
set(DUNE_DIR ${HOME}/cpp/dune)
set(ALU_DIR ${HOME}/cpp/ALUGrid-1.52/build)
set(ALBERTA_DIR ${HOME}/cpp/alberta-3.0.1/build)

# Get local folder
foreach(local_dir include src test src/config)
  set(LOCAL_INCLUDE ${LOCAL_INCLUDE} ${CMAKE_SOURCE_DIR}/${local_dir})
endforeach(local_dir)

# Get dynamic assert and ODE_Solver
set(LOCAL_INCLUDE ${LOCAL_INCLUDE} ${HOME}/cpp/dassert/include 
 ${HOME}/cpp/ODE_Solver)

# Get dune includes
foreach(dune_dir build/include dune-common dune-geometry 
    dune-istl dune-localfunctions dune-grid dune-fem dune-alugrid)
  set(DUNE_INCLUDE ${DUNE_INCLUDE} ${DUNE_DIR}/${dune_dir})
endforeach(dune_dir)

include_directories(${LOCAL_INCLUDE} ${DUNE_INCLUDE})

set(LINK_DIRS "${DUNE_DIR}/build/lib" "${ALU_DIR}/lib"
  "${ALBERTA_DIR}/lib" "${HOME}/cpp/alberta-3.0.1/build/lib"
  "${HOME}/cpp/dune/build/lib" "${HOME}/dune/install-prefix/lib/")
foreach(lib dunecommon dunegeometry dunegrid 
    dunealugrid dunefem alugrid alugrid_serial dunealbertagrid_3d 
    lapack blas z m gmp gmpxx alberta_utilities)
  find_library("${lib}LIB" "${lib}" ${LINK_DIRS})
  set(LIBS "${LIBS}" "${${lib}LIB}")
endforeach(lib)

set(PRIVATE_DEFI HAVE_CONFIG_H ENABLE_GMP=1 ALUGRID_SIMPLEX 
  ALBERTA_DIM=3 ENABLE_ALBERTA ENABLE_ZLIB=1 
  GRIDDIM=2 WORLDDIM=3 POLORDER=1)
set(PUBLIC_DEFI DEBUG=1 CURRENT_MODE=Mode::throw_
  PFILE="/Users/christianpower/cpp/DISS_surface/data/tumor_parameter.txt"
  FEF_PATH="/Users/christianpower/cpp/DISS_surface/output/"
  SPHEREPATH="/Users/christianpower/cpp/DISS_surface/data/sphere")

add_executable(max_h src/maxH_main.cpp)
target_compile_definitions(max_h PUBLIC ${PUBLIC_DEFI})

foreach(src brusselator_algo.cpp 
  io_parameter_impl.cpp secOrd_op_identity.cpp secOrd_op_rhs.cpp
  brusselator_algo_impl.cpp grid_fef.cpp io_paraview.cpp secOrd_op_initData.cpp
  secOrd_op_rhs_impl.cpp brusselator_main.cpp io_dgf.cpp 
  secOrd_op_initData_impl.cpp secOrd_op_rhs_u.cpp grid_GridAndTime.cpp
  io_errorStream.cpp secOrd_op_initData_u.cpp secOrd_op_rhs_w.cpp
  grid_GridAndTime_impl.cpp io_l2h1Calculator.cpp secOrd_op_brusselator.cpp
  secOrd_op_initData_w.cpp secOrd_op_solutionDriven.cpp grid_deformation.cpp
  io_parameter.cpp secOrd_op_brusselator_impl.cpp secOrd_op_linearHeat.cpp
  secOrd_op_solutionDriven_impl.cpp)
  set(BR_SRC "${BR_SRC}" "src/${src}")
endforeach(src)
add_executable(brusselator ${BR_SRC})
target_link_libraries(brusselator ${LIBS})
target_compile_definitions(brusselator
  PUBLIC ${PUBLIC_DEFI} PRIVATE ${PRIVATE_DEFI})
# to delete
# linHeat_main.cpp grid_evolution.cpp
