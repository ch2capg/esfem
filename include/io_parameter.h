/*! \file io_parameter.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_PARAMETER_H
#define IO_PARAMETER_H 

#include <string>
#include <vector>
#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    class Parameter{
    public:
      explicit Parameter(int argc, char** argv, const std::string& parameter_file_name);
      ~Parameter();
      Parameter(const Parameter&) = delete;
      Parameter& operator=(const Parameter&) = delete;

      /*! \name File names*/
      //@{
      const std::string& grid() const noexcept;
      const std::string& error_log() const noexcept;
      const std::string& paraview() const noexcept;
      //@}

      /*! \name Time control*/
      //@{
      double start_time() const noexcept;
      double global_timeStep() const noexcept;
      //@}
      
      /*! \name Flow control parameter */
      //@{
      long max_timeSteps() const; 
      long prePattern_timeSteps() const;
      long pattern_timeSteps() const; 
      //@}
      
      double eps() const noexcept;
      /*!< \brief Generic precision */

      /*! \name BDF coefficients */
      //@{
      const std::vector<double>& bdf_alphas() const noexcept;
      const std::vector<double>& bdf_gammas() const noexcept;
      //@}

      /*! \name Tumor growth PDE data */
      //@{
      double tg_a() const noexcept;
      double tg_b() const noexcept;
      double tg_Dc() const noexcept;
      double tg_gamma() const noexcept;
      double velocity_regularization() const noexcept;
      double surface_growthFactor() const noexcept;
      double mcf_regularization() const noexcept;
      //@}

      /*! \name Tumor growth initial data */
      //@{
      double u_hom_value() const noexcept;
      double w_hom_value() const noexcept;
      double u_pertubation() const noexcept;
      double w_pertubation() const noexcept;
      const std::string& u_init_dof() const noexcept;
      /*!< \brief File name for backup of initial nodal values */
      const std::string& w_init_dof() const noexcept;
      /*!< \brief File name for backup of initial nodal values */
      //@}
      
      friend std::ostream& operator<<(std::ostream&, const Parameter&);
      /*! \brief Prints out all data members. */
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // IO_PARAMETER_H

/*! Log:
 */
