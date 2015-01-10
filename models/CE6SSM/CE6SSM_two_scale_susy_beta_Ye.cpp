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

// File generated at Fri 9 Jan 2015 15:01:31

#include "CE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 0.7*Sqr(gN)) + 3*(Ye*Ye.adjoint()*Ye
      ));


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = twoLoop*(0.025*Ye*(756*Power(g1,4) + 660*Power(g2,4) + 273*
      Power(gN,4) - 360*traceYdAdjYdYdAdjYd - 120*traceYdAdjYuYuAdjYd - 120*
      traceYeAdjYeYeAdjYe + 48*traceYeAdjYe*Sqr(g1) + 72*Sqr(g1)*Sqr(g2) - 8*
      traceYeAdjYe*Sqr(gN) + 6*Sqr(g1)*Sqr(gN) + 78*Sqr(g2)*Sqr(gN) + 40*AbsSqr
      (Lambdax)*(-3*traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 3*
      traceYuAdjYu + Sqr(gN)) - 8*traceYdAdjYd*(2*Sqr(g1) - 80*Sqr(g3) + 3*Sqr(
      gN)) - 120*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 1.5*(-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 4*Sqr(g2) + Sqr(gN))*(Ye*Ye.adjoint()*
      Ye) - 4*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye));


   return beta_Ye;
}

} // namespace flexiblesusy
