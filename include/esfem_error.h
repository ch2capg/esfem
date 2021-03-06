/*! \file esfem_error.h
    \author Christian Power
    \date 20. March 2016

    \brief All error classes for namespace `Esfem`

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Use the facilities provided in `dassert.h` to give very easy
     assertion checking.

         Created by Christian Power on 18.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.

*/

#ifndef ESFEM_ERROR_H
#define ESFEM_ERROR_H 

#include <stdexcept>
#include "dassert.h"

namespace Esfem{
  struct Parameter_error : std::runtime_error{
    explicit Parameter_error(const std::string& msg) 
      : std::runtime_error {"Esfem::Io::Parameter.\n" + msg}
    {}
  };
  /*!< \brief Used in io_parameter{,_impl}.{h,cpp} */
  struct SolutionDriven_error : std::runtime_error{
    explicit SolutionDriven_error(const std::string& msg)
      : std::runtime_error
    {"Esfem::SecOrd_op::Solution_driven.\n" + msg}
    {}
  };
  /*!< \brief Used in secOrd_op_solutionDriven{,_impl}.{h,cpp} */
  struct BrusselatorScheme_error : std::runtime_error{
    explicit BrusselatorScheme_error(const std::string& msg)
      : std::runtime_error
    {"Esfem::Brusselator_scheme.\n" + msg}
    {}
  };
  /*!< \brief Used in brusselator_algo{,_impl}.{h,cpp} */
  struct Grid_error : std::runtime_error{
    explicit Grid_error(const std::string& msg)
      : std::runtime_error
      {"Esfem::Grid::Grid_and_time.\n" + msg}
    {}
  };
  /*!< \brief Used in grid.h */
  //! Used for secOrd_op_rhs.h
  struct Rhs_error : std::runtime_error{
    //! Inherit constructor and prefix one line to `msg`
    explicit Rhs_error(const std::string& msg)
      :std::runtime_error {"Esfem::Rhs.\n" + msg}
    {}
  };
  //! Used for secOrd_op_initData.h
  struct InitData_error : std::runtime_error{
    //! Inherit constructor and prefix one line to `msg`
    explicit InitData_error(const std::string& msg)
      :std::runtime_error {"Esfem::InitData.\n" + msg}
    {}
  };
}

// ----------------------------------------------------------------------
// Assert::dynamic specialization

template<>
inline void Assert::dynamic<false, Esfem::Parameter_error>
(const bool, const std::string&){};
template<>
inline void Assert::dynamic<false, Esfem::SolutionDriven_error>
(const bool, const std::string&){};
template<>
inline void Assert::dynamic<false, Esfem::BrusselatorScheme_error>
(const bool, const std::string&){};
template<>
inline void Assert::dynamic<false, Esfem::Grid_error>
(const bool, const std::string&){};
template<>
inline void Assert::dynamic<false, Esfem::Rhs_error>
(const bool, const std::string&){};
template<>
inline void Assert::dynamic<false, Esfem::InitData_error>
(const bool, const std::string&){};

#endif // ESFEM_ERROR_H
