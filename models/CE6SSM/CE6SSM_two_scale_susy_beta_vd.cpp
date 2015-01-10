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
 * Calculates the one-loop beta function of vd.
 *
 * @return one-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = 0.05*oneOver16PiSqr*vd*(-60*traceYdAdjYd - 20*traceYeAdjYe -
      20*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN));


   return beta_vd;
}

/**
 * Calculates the two-loop beta function of vd.
 *
 * @return two-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vd_two_loop(const Susy_traces& susy_traces) const
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


   double beta_vd;

   beta_vd = -0.00125*twoLoop*vd*(1188*Power(g1,4) + 2900*Power(g2,4) +
      1773*Power(gN,4) - 7200*traceYdAdjYdYdAdjYd - 2400*traceYdAdjYuYuAdjYd -
      2400*traceYeAdjYeYeAdjYe + 1200*traceYeAdjYe*Sqr(g1) + 1200*traceYeAdjYe*
      Sqr(g2) + 360*Sqr(g1)*Sqr(g2) + 200*traceYeAdjYe*Sqr(gN) - 36*Sqr(g1)*Sqr
      (gN) + 540*Sqr(g2)*Sqr(gN) + 200*traceYdAdjYd*(2*Sqr(g1) + 18*Sqr(g2) +
      64*Sqr(g3) + 3*Sqr(gN)) + 40*AbsSqr(Lambdax)*(-60*traceKappaAdjKappa - 40
      *traceLambda12AdjLambda12 - 60*traceYuAdjYu + 6*Sqr(g1) + 30*Sqr(g2) + 29
      *Sqr(gN)) - 2400*Sqr(Conj(Lambdax))*Sqr(Lambdax));


   return beta_vd;
}

} // namespace flexiblesusy
