// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Fri 9 Jan 2015 15:02:25

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of msI2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_msI2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_msI2;

   beta_msI2 = 0.5*gN*oneOver16PiSqr*(3.1622776601683795*Tr14 - 10*gN*
      AbsSqr(MassBp))*UNITMATRIX(2);


   return beta_msI2;
}

/**
 * Calculates the two-loop beta function of msI2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_msI2_two_loop(const Soft_traces& soft_traces) const
{
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,2,2> beta_msI2;

   beta_msI2 = 0.25*gN*twoLoop*(20*gN*Tr2U144 + 25.298221281347036*Tr34 +
      639*Power(gN,3)*AbsSqr(MassBp))*UNITMATRIX(2);


   return beta_msI2;
}

} // namespace flexiblesusy
