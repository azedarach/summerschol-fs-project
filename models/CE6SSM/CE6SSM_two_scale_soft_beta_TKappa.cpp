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

// File generated at Fri 9 Jan 2015 15:02:13

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TKappa_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = oneOver16PiSqr*(3*traceKappaAdjKappa*TKappa + 2*
      traceLambda12AdjLambda12*TKappa + 2*AbsSqr(Lambdax)*TKappa -
      0.26666666666666666*Sqr(g1)*TKappa - 5.333333333333333*Sqr(g3)*TKappa -
      1.9*Sqr(gN)*TKappa + Kappa*(6*traceAdjKappaTKappa + 4*
      traceAdjLambda12TLambda12 + 0.5333333333333333*MassB*Sqr(g1) +
      10.666666666666666*MassG*Sqr(g3) + 3.8*MassBp*Sqr(gN) + 4*Conj(Lambdax)*
      TLambdax) + 3*(Kappa*(Kappa).adjoint()*TKappa) + 3*(TKappa*(Kappa)
      .adjoint()*Kappa));


   return beta_TKappa;
}

/**
 * Calculates the two-loop beta function of TKappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TKappa_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = twoLoop*(2.5955555555555554*Power(g1,4)*TKappa +
      14.222222222222221*Power(g3,4)*TKappa + 19.665*Power(gN,4)*TKappa - 6*
      traceKappaAdjKappaKappaAdjKappa*TKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TKappa - 6*traceYdAdjYd*
      AbsSqr(Lambdax)*TKappa - 2*traceYeAdjYe*AbsSqr(Lambdax)*TKappa - 6*
      traceYuAdjYu*AbsSqr(Lambdax)*TKappa + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TKappa + 1.2*traceLambda12AdjLambda12*Sqr(g1)*TKappa + 1.2*AbsSqr(Lambdax
      )*Sqr(g1)*TKappa + 6*traceLambda12AdjLambda12*Sqr(g2)*TKappa + 6*AbsSqr(
      Lambdax)*Sqr(g2)*TKappa + 16*traceKappaAdjKappa*Sqr(g3)*TKappa +
      1.4222222222222223*Sqr(g1)*Sqr(g3)*TKappa - 1.8*traceKappaAdjKappa*Sqr(gN
      )*TKappa - 1.2*traceLambda12AdjLambda12*Sqr(gN)*TKappa - 1.2*AbsSqr(
      Lambdax)*Sqr(gN)*TKappa + 0.17333333333333334*Sqr(g1)*Sqr(gN)*TKappa +
      3.466666666666667*Sqr(g3)*Sqr(gN)*TKappa - 4*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TKappa - 0.0022222222222222222*Kappa*(4672*Power(g1,4)*MassB +
      35397*Power(gN,4)*MassBp + 25600*Power(g3,4)*MassG + 10800*
      traceKappaAdjKappaTKappaAdjKappa + 7200*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 720*traceAdjKappaTKappa*
      Sqr(g1) - 1080*traceAdjLambda12TLambda12*Sqr(g1) + 1080*MassB*
      traceLambda12AdjLambda12*Sqr(g1) - 5400*traceAdjLambda12TLambda12*Sqr(g2)
      + 5400*MassWB*traceLambda12AdjLambda12*Sqr(g2) - 14400*
      traceAdjKappaTKappa*Sqr(g3) + 1280*MassB*Sqr(g1)*Sqr(g3) + 1280*MassG*Sqr
      (g1)*Sqr(g3) + 1620*traceAdjKappaTKappa*Sqr(gN) + 1080*
      traceAdjLambda12TLambda12*Sqr(gN) - 1080*MassBp*traceLambda12AdjLambda12*
      Sqr(gN) + 156*MassB*Sqr(g1)*Sqr(gN) + 156*MassBp*Sqr(g1)*Sqr(gN) + 3120*
      MassBp*Sqr(g3)*Sqr(gN) + 3120*MassG*Sqr(g3)*Sqr(gN) + 180*
      traceKappaAdjKappa*(4*MassB*Sqr(g1) + 80*MassG*Sqr(g3) - 9*MassBp*Sqr(gN)
      ) + 7200*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 360*Conj(Lambdax)*(Lambdax
      *(15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*Sqr(g1)
      + 15*MassWB*Sqr(g2) - 3*MassBp*Sqr(gN)) + (15*traceYdAdjYd + 5*
      traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) + 3*Sqr(gN))*
      TLambdax)) - (12*traceAdjKappaTKappa + 8*traceAdjLambda12TLambda12 + 5*
      MassBp*Sqr(gN) + 8*Conj(Lambdax)*TLambdax)*(Kappa*(Kappa).adjoint()*Kappa
      ) - 9*traceKappaAdjKappa*(Kappa*(Kappa).adjoint()*TKappa) - 6*
      traceLambda12AdjLambda12*(Kappa*(Kappa).adjoint()*TKappa) - 6*AbsSqr(
      Lambdax)*(Kappa*(Kappa).adjoint()*TKappa) + 3.5*Sqr(gN)*(Kappa*(Kappa)
      .adjoint()*TKappa) - 9*traceKappaAdjKappa*(TKappa*(Kappa).adjoint()*Kappa
      ) - 6*traceLambda12AdjLambda12*(TKappa*(Kappa).adjoint()*Kappa) - 6*
      AbsSqr(Lambdax)*(TKappa*(Kappa).adjoint()*Kappa) + 4*Sqr(gN)*(TKappa*(
      Kappa).adjoint()*Kappa) - 3*(Kappa*(Kappa).adjoint()*Kappa*(Kappa)
      .adjoint()*TKappa) - 4*(Kappa*(Kappa).adjoint()*TKappa*(Kappa).adjoint()*
      Kappa) - 3*(TKappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*Kappa));


   return beta_TKappa;
}

} // namespace flexiblesusy
