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
 * Calculates the one-loop beta function of Lambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_susy_parameters::calc_beta_Lambda12_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = oneOver16PiSqr*(Lambda12*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2) -
      1.9*Sqr(gN)) + 2*(Lambda12*(Lambda12).adjoint()*Lambda12));


   return beta_Lambda12;
}

/**
 * Calculates the two-loop beta function of Lambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_susy_parameters::calc_beta_Lambda12_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = twoLoop*(0.005*Lambda12*(1188*Power(g1,4) + 3300*Power
      (g2,4) + 3933*Power(gN,4) - 1200*traceKappaAdjKappaKappaAdjKappa - 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 240*
      traceLambda12AdjLambda12*Sqr(g1) + 1200*traceLambda12AdjLambda12*Sqr(g2)
      + 360*Sqr(g1)*Sqr(g2) - 80*AbsSqr(Lambdax)*(15*traceYdAdjYd + 5*
      traceYeAdjYe - 3*(-5*traceYuAdjYu + Sqr(g1) + 5*Sqr(g2) - Sqr(gN))) + 40*
      traceKappaAdjKappa*(4*Sqr(g1) + 80*Sqr(g3) - 9*Sqr(gN)) - 240*
      traceLambda12AdjLambda12*Sqr(gN) + 54*Sqr(g1)*Sqr(gN) + 390*Sqr(g2)*Sqr(
      gN) - 800*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-6*traceKappaAdjKappa - 4*
      traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) + 2.5*Sqr(gN))*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 2*(Lambda12*(Lambda12).adjoint()*Lambda12
      *(Lambda12).adjoint()*Lambda12));


   return beta_Lambda12;
}

} // namespace flexiblesusy
