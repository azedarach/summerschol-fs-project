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

// File generated at Fri 9 Jan 2015 15:02:11

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd -
      5.333333333333333*Sqr(g3)*TYd - 0.7*Sqr(gN)*TYd + Yd*(6*traceAdjYdTYd + 2
      *traceAdjYeTYe + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 1.4*MassBp*Sqr(gN) + 2*Conj(Lambdax)*
      TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*
      Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu);


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = twoLoop*(4.588888888888889*Power(g1,4)*TYd + 16.5*Power(g2,
      4)*TYd + 14.222222222222221*Power(g3,4)*TYd + 6.825*Power(gN,4)*TYd - 9*
      traceYdAdjYdYdAdjYd*TYd - 3*traceYdAdjYuYuAdjYd*TYd - 3*
      traceYeAdjYeYeAdjYe*TYd - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd - 2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYd - 3*traceYuAdjYu*AbsSqr(
      Lambdax)*TYd - 0.4*traceYdAdjYd*Sqr(g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*
      TYd + Sqr(g1)*Sqr(g2)*TYd + 16*traceYdAdjYd*Sqr(g3)*TYd +
      0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 8*Sqr(g2)*Sqr(g3)*TYd - 0.6*
      traceYdAdjYd*Sqr(gN)*TYd - 0.2*traceYeAdjYe*Sqr(gN)*TYd + AbsSqr(Lambdax)
      *Sqr(gN)*TYd + 0.32666666666666666*Sqr(g1)*Sqr(gN)*TYd + 1.5*Sqr(g2)*Sqr(
      gN)*TYd + 1.3333333333333333*Sqr(g3)*Sqr(gN)*TYd - 3*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)*TYd - 0.0022222222222222222*Yd*(8260*Power(g1,4)*MassB +
      12285*Power(gN,4)*MassBp + 25600*Power(g3,4)*MassG + 29700*Power(g2,4)*
      MassWB + 16200*traceYdAdjYdTYdAdjYd + 2700*traceYdAdjYuTYuAdjYd + 5400*
      traceYeAdjYeTYeAdjYe + 2700*traceYuAdjYdTYdAdjYu + 360*traceAdjYdTYd*Sqr(
      g1) - 1080*traceAdjYeTYe*Sqr(g1) + 1080*MassB*traceYeAdjYe*Sqr(g1) + 900*
      MassB*Sqr(g1)*Sqr(g2) + 900*MassWB*Sqr(g1)*Sqr(g2) - 14400*traceAdjYdTYd*
      Sqr(g3) + 800*MassB*Sqr(g1)*Sqr(g3) + 800*MassG*Sqr(g1)*Sqr(g3) + 7200*
      MassG*Sqr(g2)*Sqr(g3) + 7200*MassWB*Sqr(g2)*Sqr(g3) + 540*traceAdjYdTYd*
      Sqr(gN) + 180*traceAdjYeTYe*Sqr(gN) - 180*MassBp*traceYeAdjYe*Sqr(gN) +
      294*MassB*Sqr(g1)*Sqr(gN) + 294*MassBp*Sqr(g1)*Sqr(gN) + 1350*MassBp*Sqr(
      g2)*Sqr(gN) + 1350*MassWB*Sqr(g2)*Sqr(gN) + 1200*MassBp*Sqr(g3)*Sqr(gN) +
      1200*MassG*Sqr(g3)*Sqr(gN) - 180*traceYdAdjYd*(2*MassB*Sqr(g1) - 80*
      MassG*Sqr(g3) + 3*MassBp*Sqr(gN)) + 5400*Lambdax*Sqr(Conj(Lambdax))*
      TLambdax + 900*Conj(Lambdax)*(Lambdax*(3*traceAdjKappaTKappa + 2*
      traceAdjLambda12TLambda12 + 3*traceAdjYuTYu + MassBp*Sqr(gN)) + (3*
      traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYuAdjYu - Sqr(gN
      ))*TLambdax)) - 0.4*(45*traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr(g1
      ) + 30*MassWB*Sqr(g2) + 6*MassBp*Sqr(gN) + 15*Conj(Lambdax)*TLambdax)*(Yd
      *Yd.adjoint()*Yd) - 12*traceYdAdjYd*(Yd*Yd.adjoint()*TYd) - 4*
      traceYeAdjYe*(Yd*Yd.adjoint()*TYd) - 4*AbsSqr(Lambdax)*(Yd*Yd.adjoint()*
      TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()*TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd
      ) + 1.8*Sqr(gN)*(Yd*Yd.adjoint()*TYd) - 6*traceAdjYuTYu*(Yd*Yu.adjoint()*
      Yu) - 1.6*MassB*Sqr(g1)*(Yd*Yu.adjoint()*Yu) - 0.4*MassBp*Sqr(gN)*(Yd*
      Yu.adjoint()*Yu) - 2*Conj(Lambdax)*TLambdax*(Yd*Yu.adjoint()*Yu) - 6*
      traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*AbsSqr(Lambdax)*(Yd*Yu.adjoint()*
      TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) + 0.4*Sqr(gN)*(Yd*Yu.adjoint()*
      TYu) - 15*traceYdAdjYd*(TYd*Yd.adjoint()*Yd) - 5*traceYeAdjYe*(TYd*
      Yd.adjoint()*Yd) - 5*AbsSqr(Lambdax)*(TYd*Yd.adjoint()*Yd) + 1.2*Sqr(g1)*
      (TYd*Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*Yd) + 1.8*Sqr(gN)*(
      TYd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(TYd*Yu.adjoint()*Yu) - AbsSqr(
      Lambdax)*(TYd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(TYd*Yu.adjoint()*Yu) + 0.2*
      Sqr(gN)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) -
      8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*
      Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu));


   return beta_TYd;
}

} // namespace flexiblesusy
