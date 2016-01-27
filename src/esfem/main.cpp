#include <config.h>
#include "esfem.h"
#include <iostream>
#include <memory>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/common/exceptions.hh>
// #include <dune/fem/solver/oemsolver.hh>
// #include <dune/fem/quadrature/quadrature.hh>
// #include <dune/common/fmatrix.hh>

class Error_logger{
public:
  explicit Error_logger(const Discrete_function::FE_function& exact_solution,
			const Discrete_function::FE_function& numerical_solution,
			std::ostream& err_stream);

  void write_error();
private:
  const Discrete_function::FE_function& u;
  const Discrete_function::FE_function& uN;
  std::ostream& log;
};

int main(int argc, char** argv) try{
  Dune::Fem::MPIManager::initialize(argc, argv);

  const auto parameter_file =
    "/Users/christianpower/cpp/DISS_surfaces/data/tumor_parameter.txt";  
  Parameter::PDE_data data {argc, argv, parameter_file};
#ifdef DEBUG
  std::clog << data << std::endl;
#endif

  FE_grid::init_grid_and_timeProvider(data);
  auto& time_provider = FE_grid::time_provider();
  
  Discrete_function::init_fef_and_norms(data);
  auto& exact_solution = Discrete_function::exact_solution();  
  auto& uN = Discrete_function::numerical_solution();
  auto& les_rhs = Discrete_function::rhs();
  auto& f_rhs = Discrete_function::load_vector();
  
  const Operator::Linear_heat fem_op {data, time_provider};
  const Dune::Fem::CGInverseOperator<Discrete_function::FE_function>
    solver {fem_op, data.eps(), data.eps()};

  auto& err_stream = Parameter::err_stream(data);
  Error_logger logger {exact_solution, uN, err_stream};
  
  Operator::interpolate_initData(time_provider, exact_solution);
  Operator::interpolate_initData(time_provider, uN);
  fem_op.mass_matrix(uN, les_rhs);
  
  time_provider.next(data.global_timeStep());
  for(long it = 0; it < data.max_timeSteps(); ++it){
    // Mu^n+1 + tau Au^n+1 = M^n + tau f^n+1
    Operator::assemble_RHS(time_provider, f_rhs);
    les_rhs.axpy(data.global_timeStep(), f_rhs);
    solver(les_rhs, uN);

    Operator::interpolate_initData(time_provider, exact_solution);
    logger.write_error();

    fem_op.mass_matrix(uN, les_rhs);
    
    time_provider.next(data.global_timeStep());
  }
  err_stream.close();
} 
catch(const Dune::Exception& e){
  std::cerr << "Dune error: " << e << std::endl;
  return 1;
}
catch(const std::exception& e){
  std::cerr << e.what() << std::endl;
  return 2;
}

// ----------------------------------------------------------------------
// Internal Implementation

Error_logger::Error_logger(const Discrete_function::FE_function& exact_solution,
			   const Discrete_function::FE_function& numerical_solution,
			   std::ostream& err_stream)
  : u {exact_solution}, uN {numerical_solution}, log {err_stream}
    {
      log << std::scientific;      
    }
void Error_logger::write_error(){
  static const auto& l2 = Discrete_function::l2norm();
  static const auto& h1 = Discrete_function::h1norm();
  log << l2.distance(u, uN) << ' '
      << h1.distance(u, uN) << std::endl;
}
