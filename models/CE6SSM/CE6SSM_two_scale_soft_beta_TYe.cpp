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

// File generated at Fri 9 Jan 2015 15:02:12

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = oneOver16PiSqr*(3*traceYdAdjYd*TYe + traceYeAdjYe*TYe +
      AbsSqr(Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe - 0.7*Sqr(gN)*TYe +
      Ye*(6*traceAdjYdTYd + 2*traceAdjYeTYe + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr
      (g2) + 1.4*MassBp*Sqr(gN) + 2*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint(
      )*TYe) + 5*(TYe*Ye.adjoint()*Ye));


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = twoLoop*(18.9*Power(g1,4)*TYe + 16.5*Power(g2,4)*TYe +
      6.825*Power(gN,4)*TYe - 9*traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd
      *TYe - 3*traceYeAdjYeYeAdjYe*TYe - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*
      TYe - 2*traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYe - 3*traceYuAdjYu*
      AbsSqr(Lambdax)*TYe - 0.4*traceYdAdjYd*Sqr(g1)*TYe + 1.2*traceYeAdjYe*Sqr
      (g1)*TYe + 1.8*Sqr(g1)*Sqr(g2)*TYe + 16*traceYdAdjYd*Sqr(g3)*TYe - 0.6*
      traceYdAdjYd*Sqr(gN)*TYe - 0.2*traceYeAdjYe*Sqr(gN)*TYe + AbsSqr(Lambdax)
      *Sqr(gN)*TYe + 0.51*Sqr(g1)*Sqr(gN)*TYe + 1.95*Sqr(g2)*Sqr(gN)*TYe - 3*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe - 0.02*Ye*(3780*Power(g1,4)*MassB +
      1365*Power(gN,4)*MassBp + 3300*Power(g2,4)*MassWB + 1800*
      traceYdAdjYdTYdAdjYd + 300*traceYdAdjYuTYuAdjYd + 600*
      traceYeAdjYeTYeAdjYe + 300*traceYuAdjYdTYdAdjYu + 40*traceAdjYdTYd*Sqr(g1
      ) - 120*traceAdjYeTYe*Sqr(g1) + 120*MassB*traceYeAdjYe*Sqr(g1) + 180*
      MassB*Sqr(g1)*Sqr(g2) + 180*MassWB*Sqr(g1)*Sqr(g2) - 1600*traceAdjYdTYd*
      Sqr(g3) + 60*traceAdjYdTYd*Sqr(gN) + 20*traceAdjYeTYe*Sqr(gN) - 20*MassBp
      *traceYeAdjYe*Sqr(gN) + 51*MassB*Sqr(g1)*Sqr(gN) + 51*MassBp*Sqr(g1)*Sqr(
      gN) + 195*MassBp*Sqr(g2)*Sqr(gN) + 195*MassWB*Sqr(g2)*Sqr(gN) - 20*
      traceYdAdjYd*(2*MassB*Sqr(g1) - 80*MassG*Sqr(g3) + 3*MassBp*Sqr(gN)) +
      600*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 100*Conj(Lambdax)*(Lambdax*(3*
      traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 + 3*traceAdjYuTYu +
      MassBp*Sqr(gN)) + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*
      traceYuAdjYu - Sqr(gN))*TLambdax)) - 3*(6*traceAdjYdTYd + 2*traceAdjYeTYe
      + 4*MassWB*Sqr(g2) + MassBp*Sqr(gN) + 2*Conj(Lambdax)*TLambdax)*(Ye*
      Ye.adjoint()*Ye) - 12*traceYdAdjYd*(Ye*Ye.adjoint()*TYe) - 4*traceYeAdjYe
      *(Ye*Ye.adjoint()*TYe) - 4*AbsSqr(Lambdax)*(Ye*Ye.adjoint()*TYe) + 1.2*
      Sqr(g1)*(Ye*Ye.adjoint()*TYe) + 6*Sqr(g2)*(Ye*Ye.adjoint()*TYe) + 1.8*Sqr
      (gN)*(Ye*Ye.adjoint()*TYe) - 15*traceYdAdjYd*(TYe*Ye.adjoint()*Ye) - 5*
      traceYeAdjYe*(TYe*Ye.adjoint()*Ye) - 5*AbsSqr(Lambdax)*(TYe*Ye.adjoint()*
      Ye) - 1.2*Sqr(g1)*(TYe*Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()*Ye
      ) + 2.7*Sqr(gN)*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint(
      )*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*Ye
      *Ye.adjoint()*Ye));


   return beta_TYe;
}

} // namespace flexiblesusy
