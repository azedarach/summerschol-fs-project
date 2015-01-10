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

// File generated at Fri 9 Jan 2015 15:01:32

#include "CE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = 0.1*oneOver16PiSqr*Lambdax*(30*traceKappaAdjKappa + 20*
      traceLambda12AdjLambda12 + 30*traceYdAdjYd + 10*traceYeAdjYe + 30*
      traceYuAdjYu + 40*AbsSqr(Lambdax) - 6*Sqr(g1) - 30*Sqr(g2) - 19*Sqr(gN));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = -0.005*twoLoop*Lambdax*(-1188*Power(g1,4) - 3300*Power(
      g2,4) - 3933*Power(gN,4) + 1200*traceKappaAdjKappaKappaAdjKappa + 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1800*traceYdAdjYdYdAdjYd +
      1200*traceYdAdjYuYuAdjYd + 600*traceYeAdjYeYeAdjYe + 1800*
      traceYuAdjYuYuAdjYu - 160*traceKappaAdjKappa*Sqr(g1) - 240*
      traceLambda12AdjLambda12*Sqr(g1) - 240*traceYeAdjYe*Sqr(g1) - 160*
      traceYuAdjYu*Sqr(g1) - 1200*traceLambda12AdjLambda12*Sqr(g2) - 360*Sqr(g1
      )*Sqr(g2) - 3200*traceKappaAdjKappa*Sqr(g3) - 3200*traceYuAdjYu*Sqr(g3) +
      360*traceKappaAdjKappa*Sqr(gN) + 240*traceLambda12AdjLambda12*Sqr(gN) +
      40*traceYeAdjYe*Sqr(gN) + 60*traceYuAdjYu*Sqr(gN) - 54*Sqr(g1)*Sqr(gN) -
      390*Sqr(g2)*Sqr(gN) + 40*traceYdAdjYd*(2*Sqr(g1) - 80*Sqr(g3) + 3*Sqr(gN)
      ) - 20*AbsSqr(Lambdax)*(-60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 90*traceYdAdjYd - 30*traceYeAdjYe - 90*
      traceYuAdjYu + 12*Sqr(g1) + 60*Sqr(g2) + 13*Sqr(gN)) + 2000*Sqr(Conj(
      Lambdax))*Sqr(Lambdax));


   return beta_Lambdax;
}

} // namespace flexiblesusy
