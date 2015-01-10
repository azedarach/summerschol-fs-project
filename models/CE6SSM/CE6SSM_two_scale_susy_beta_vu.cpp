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

// File generated at Fri 9 Jan 2015 15:01:33

#include "CE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = 0.1*oneOver16PiSqr*vu*(-30*traceYuAdjYu - 10*AbsSqr(Lambdax)
      + 3*Sqr(g1) + 15*Sqr(g2) + 2*Sqr(gN));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = -0.005*twoLoop*vu*(297*Power(g1,4) + 725*Power(g2,4) + 192*
      Power(gN,4) - 600*traceYdAdjYuYuAdjYd - 1800*traceYuAdjYuYuAdjYu + 90*Sqr
      (g1)*Sqr(g2) + 36*Sqr(g1)*Sqr(gN) + 60*Sqr(g2)*Sqr(gN) + 20*traceYuAdjYu*
      (17*Sqr(g1) + 45*Sqr(g2) + 160*Sqr(g3) + 3*Sqr(gN)) + 20*AbsSqr(Lambdax)*
      (-30*traceKappaAdjKappa - 20*traceLambda12AdjLambda12 - 30*traceYdAdjYd -
      10*traceYeAdjYe + 3*Sqr(g1) + 15*Sqr(g2) + 17*Sqr(gN)) - 600*Sqr(Conj(
      Lambdax))*Sqr(Lambdax));


   return beta_vu;
}

} // namespace flexiblesusy
