#include "discrete_function.h"
#include "check_pseudo_externals.h"
#include <stdexcept>

constexpr int no_externalVars = 8;
static Check_pseudo_externals::Control_table<no_externalVars> control_vars {};

void Discrete_function::init_fef_and_norms(const Parameter::PDE_data& d) try{
  auto& grid_part = FE_grid::grid_part();

  l2norm(&grid_part);
  h1norm(&grid_part);

  auto& fe_sp = fe_space(&grid_part);

  const std::string Un_name {"U_n"};
  numerical_solution(&Un_name, &fe_sp);
  const std::string U_name {"exact_solution"};
  exact_solution(&U_name, &fe_sp);
  const std::string rhs_name {"rhs"};
  rhs(&rhs_name, &fe_sp);
  const std::string xi_name {"xi"};
  fef_feOperator(&xi_name, &fe_sp);
  const std::string loadVector_name {"load_vector"};
  load_vector(&loadVector_name, &fe_sp);
}
catch(const std::exception& e){
  std::string err_msg {"Error in init_fef_and_norms().\n"};
  err_msg += e.what();
  throw std::runtime_error {err_msg};
}
Discrete_function::L2_norm&
Discrete_function::l2norm(FE_grid::Grid_part* const gp_ptr){
  constexpr int var_no = 0;
  if(!gp_ptr && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in l2norm().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static L2_norm l2norm {*gp_ptr};
  control_vars.just_init<var_no>();
  return l2norm;
}
Discrete_function::H1_norm&
Discrete_function::h1norm(FE_grid::Grid_part* const gp_ptr){
  constexpr int var_no = 1;
  if(!gp_ptr && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in h1norm().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static H1_norm h1norm {*gp_ptr};
  control_vars.just_init<var_no>();
  return h1norm;  
}
Discrete_function::FE_space& Discrete_function::fe_space(FE_grid::Grid_part* gp_ptr){
  constexpr int var_no = 2;
  if(!gp_ptr && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in fe_space().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static Discrete_function::FE_space fe_space {*gp_ptr};
  control_vars.just_init<var_no>();
  return fe_space;
}
Discrete_function::FE_function& Discrete_function::
numerical_solution(const std::string* const s_ptr,
		   Discrete_function::FE_space* const fes_ptr){
  constexpr int var_no = 3;
  if((!s_ptr || !fes_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in numerical_solution().  "
	"This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_function U_n(*s_ptr, *fes_ptr);
  control_vars.just_init<var_no>();
  return U_n;
}
Discrete_function::FE_function& Discrete_function::
exact_solution(const std::string* const s_ptr,
		   Discrete_function::FE_space* const fes_ptr){
  constexpr int var_no = 4;
  if((!s_ptr || !fes_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in exact_solution().  "
	"This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_function U_n(*s_ptr, *fes_ptr);
  control_vars.just_init<var_no>();
  return U_n;
}
Discrete_function::FE_function& Discrete_function::
rhs(const std::string* const s_ptr,
    Discrete_function::FE_space* const fes_ptr){
  constexpr int var_no = 5;
  if((!s_ptr || !fes_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in rhs().  "
	"This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_function U_n(*s_ptr, *fes_ptr);
  control_vars.just_init<var_no>();
  return U_n;
}
Discrete_function::FE_function& Discrete_function::
fef_feOperator(const std::string* const s_ptr,
	       Discrete_function::FE_space* const fes_ptr){
  constexpr int var_no = 6;
  if((!s_ptr || !fes_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in fef_feOperator().  "
	"This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_function U_n(*s_ptr, *fes_ptr);
  control_vars.just_init<var_no>();
  return U_n;
}
Discrete_function::FE_function& Discrete_function::
load_vector(const std::string* const s_ptr,
	    Discrete_function::FE_space* const fes_ptr){
  constexpr int var_no = 7;
  if((!s_ptr || !fes_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in load_vector().  "
	"This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_function U_n(*s_ptr, *fes_ptr);
  control_vars.just_init<var_no>();
  return U_n;
}
