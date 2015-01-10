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

// File generated at Fri 9 Jan 2015 15:02:28

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassBp.
 *
 * @return one-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassBp_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassBp;

   beta_MassBp = 18.8*MassBp*oneOver16PiSqr*Sqr(gN);


   return beta_MassBp;
}

/**
 * Calculates the two-loop beta function of MassBp.
 *
 * @return two-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassBp_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_MassBp;

   beta_MassBp = 0.04*twoLoop*Sqr(gN)*(285*traceAdjKappaTKappa + 190*
      traceAdjLambda12TLambda12 + 210*traceAdjYdTYd + 70*traceAdjYeTYe + 90*
      traceAdjYuTYu - 285*MassBp*traceKappaAdjKappa - 190*MassBp*
      traceLambda12AdjLambda12 - 210*MassBp*traceYdAdjYd - 70*MassBp*
      traceYeAdjYe - 90*MassBp*traceYuAdjYu + 162*MassB*Sqr(g1) + 162*MassBp*
      Sqr(g1) + 510*MassBp*Sqr(g2) + 510*MassWB*Sqr(g2) + 1200*MassBp*Sqr(g3) +
      1200*MassG*Sqr(g3) + 916*MassBp*Sqr(gN) - 190*Conj(Lambdax)*(MassBp*
      Lambdax - TLambdax));


   return beta_MassBp;
}

} // namespace flexiblesusy
