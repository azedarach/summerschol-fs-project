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
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = oneOver16PiSqr*(-0.03333333333333333*Yu*(-90*traceYuAdjYu -
      30*AbsSqr(Lambdax) + 26*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3) + 9*Sqr(gN)) +
      Yu*Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu));


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = twoLoop*(Yu*(8.695555555555556*Power(g1,4) + 16.5*Power(g2,4
      ) + 14.222222222222221*Power(g3,4) + 2.865*Power(gN,4) - 3*
      traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu + Sqr(g1)*Sqr(g2) +
      3.022222222222222*Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) + 0.1*traceYuAdjYu*
      (8*Sqr(g1) + 160*Sqr(g3) - 3*Sqr(gN)) + 0.5366666666666666*Sqr(g1)*Sqr(gN
      ) + 0.75*Sqr(g2)*Sqr(gN) + 0.5333333333333333*Sqr(g3)*Sqr(gN) + 0.5*
      AbsSqr(Lambdax)*(-6*traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 6*
      traceYdAdjYd - 2*traceYeAdjYe + 3*Sqr(gN)) - 3*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)) + 0.2*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) +
      2*Sqr(g1) + 3*Sqr(gN))*(Yu*Yd.adjoint()*Yd) - 9*traceYuAdjYu*(Yu*
      Yu.adjoint()*Yu) - 3*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*Yu) + 0.4*Sqr(g1)*(
      Yu*Yu.adjoint()*Yu) + 6*Sqr(g2)*(Yu*Yu.adjoint()*Yu) + 0.6*Sqr(gN)*(Yu*
      Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      );


   return beta_Yu;
}

} // namespace flexiblesusy
