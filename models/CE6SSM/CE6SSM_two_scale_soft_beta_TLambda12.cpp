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
 * Calculates the one-loop beta function of TLambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_TLambda12_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = oneOver16PiSqr*(Lambda12*(6*traceAdjKappaTKappa + 4*
      traceAdjLambda12TLambda12 + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 3.8*
      MassBp*Sqr(gN) + 4*Conj(Lambdax)*TLambdax) + 3*traceKappaAdjKappa*
      TLambda12 + 2*traceLambda12AdjLambda12*TLambda12 + 2*AbsSqr(Lambdax)*
      TLambda12 - 0.6*Sqr(g1)*TLambda12 - 3*Sqr(g2)*TLambda12 - 1.9*Sqr(gN)*
      TLambda12 + 3*(Lambda12*(Lambda12).adjoint()*TLambda12) + 3*(TLambda12*(
      Lambda12).adjoint()*Lambda12));


   return beta_TLambda12;
}

/**
 * Calculates the two-loop beta function of TLambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_TLambda12_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = twoLoop*(-0.02*Lambda12*(1188*Power(g1,4)*MassB +
      3933*Power(gN,4)*MassBp + 3300*Power(g2,4)*MassWB + 1200*
      traceKappaAdjKappaTKappaAdjKappa + 800*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 80*traceAdjKappaTKappa*Sqr
      (g1) - 120*traceAdjLambda12TLambda12*Sqr(g1) + 120*MassB*
      traceLambda12AdjLambda12*Sqr(g1) - 600*traceAdjLambda12TLambda12*Sqr(g2)
      + 600*MassWB*traceLambda12AdjLambda12*Sqr(g2) + 180*MassB*Sqr(g1)*Sqr(g2)
      + 180*MassWB*Sqr(g1)*Sqr(g2) - 1600*traceAdjKappaTKappa*Sqr(g3) + 180*
      traceAdjKappaTKappa*Sqr(gN) + 120*traceAdjLambda12TLambda12*Sqr(gN) - 120
      *MassBp*traceLambda12AdjLambda12*Sqr(gN) + 39*MassB*Sqr(g1)*Sqr(gN) + 39*
      MassBp*Sqr(g1)*Sqr(gN) + 195*MassBp*Sqr(g2)*Sqr(gN) + 195*MassWB*Sqr(g2)*
      Sqr(gN) + 20*traceKappaAdjKappa*(4*MassB*Sqr(g1) + 80*MassG*Sqr(g3) - 9*
      MassBp*Sqr(gN)) + 800*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 40*Conj(
      Lambdax)*(Lambdax*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu
      + 3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2) - 3*MassBp*Sqr(gN)) + (15*
      traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2)
      + 3*Sqr(gN))*TLambdax)) + 5.94*Power(g1,4)*TLambda12 + 16.5*Power(g2,4)*
      TLambda12 + 19.665*Power(gN,4)*TLambda12 - 6*
      traceKappaAdjKappaKappaAdjKappa*TLambda12 - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TLambda12 - 6*traceYdAdjYd*
      AbsSqr(Lambdax)*TLambda12 - 2*traceYeAdjYe*AbsSqr(Lambdax)*TLambda12 - 6*
      traceYuAdjYu*AbsSqr(Lambdax)*TLambda12 + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TLambda12 + 1.2*traceLambda12AdjLambda12*Sqr(g1)*TLambda12 + 1.2*AbsSqr(
      Lambdax)*Sqr(g1)*TLambda12 + 6*traceLambda12AdjLambda12*Sqr(g2)*TLambda12
      + 6*AbsSqr(Lambdax)*Sqr(g2)*TLambda12 + 1.8*Sqr(g1)*Sqr(g2)*TLambda12 +
      16*traceKappaAdjKappa*Sqr(g3)*TLambda12 - 1.8*traceKappaAdjKappa*Sqr(gN)*
      TLambda12 - 1.2*traceLambda12AdjLambda12*Sqr(gN)*TLambda12 - 1.2*AbsSqr(
      Lambdax)*Sqr(gN)*TLambda12 + 0.39*Sqr(g1)*Sqr(gN)*TLambda12 + 1.95*Sqr(g2
      )*Sqr(gN)*TLambda12 - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambda12 - (12*
      traceAdjKappaTKappa + 8*traceAdjLambda12TLambda12 + 5*MassBp*Sqr(gN) + 8*
      Conj(Lambdax)*TLambdax)*(Lambda12*(Lambda12).adjoint()*Lambda12) - 9*
      traceKappaAdjKappa*(Lambda12*(Lambda12).adjoint()*TLambda12) - 6*
      traceLambda12AdjLambda12*(Lambda12*(Lambda12).adjoint()*TLambda12) - 6*
      AbsSqr(Lambdax)*(Lambda12*(Lambda12).adjoint()*TLambda12) + 3.5*Sqr(gN)*(
      Lambda12*(Lambda12).adjoint()*TLambda12) - 9*traceKappaAdjKappa*(
      TLambda12*(Lambda12).adjoint()*Lambda12) - 6*traceLambda12AdjLambda12*(
      TLambda12*(Lambda12).adjoint()*Lambda12) - 6*AbsSqr(Lambdax)*(TLambda12*(
      Lambda12).adjoint()*Lambda12) + 4*Sqr(gN)*(TLambda12*(Lambda12).adjoint()
      *Lambda12) - 3*(Lambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint
      ()*TLambda12) - 4*(Lambda12*(Lambda12).adjoint()*TLambda12*(Lambda12)
      .adjoint()*Lambda12) - 3*(TLambda12*(Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12));


   return beta_TLambda12;
}

} // namespace flexiblesusy
