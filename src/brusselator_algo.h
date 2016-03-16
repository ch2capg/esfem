/*! \file brusselator_algo.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Pseudocode:
       (1+τγ)(Mu)ⁿ⁺¹ + τ (Au)ⁿ⁺¹ - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ = (Mu)ⁿ + τγa Mⁿ⁺¹1
       (Mw)ⁿ⁺¹ + τDᶜ (Aw)ⁿ⁺¹ + τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹ = (Mw)ⁿ + τγb Mⁿ⁺¹1

     Created by Christian Power on 11.03.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_H
#define BRUSSELATOR_ALGO_H 

#include <string>
#include <memory>
#include "esfem.h"


namespace Esfem{
  void brusselator_algo(int argc, char** argv);
  // ESFEM algorithm.  Only this should be invoked by main

  class Brusselator_scheme{
  public:
    explicit Brusselator_scheme(int argc, char** argv, const std::string& parameter_fname);
    ~Brusselator_scheme();

    void next_timeStep(); 
    long prePattern_timeSteps() const; 
    long pattern_timeSteps() const; 

    void pre_pattern_action();
    // void intermediate_action();
    void pattern_action();
    void final_action();
  
  private:
    struct Data;
    std::unique_ptr<Data> d_ptr;
    void pre_loop_action(); // to be invoked only in the constructor
  };

}

#endif // BRUSSELATOR_ALGO_H

/*! Log:
 */
