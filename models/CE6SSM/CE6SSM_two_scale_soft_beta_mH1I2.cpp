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

// File generated at Fri 9 Jan 2015 15:02:25

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mH1I2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_mH1I2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = oneOver16PiSqr*(2*ms2*((Lambda12).adjoint()*Lambda12) + 2
      *((TLambda12).adjoint()*TLambda12) + mH1I2*(Lambda12).adjoint()*Lambda12
      + 2*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) + (Lambda12)
      .adjoint()*Lambda12*mH1I2 - 0.7745966692414834*g1*Tr11*UNITMATRIX(2) -
      0.9486832980505138*gN*Tr14*UNITMATRIX(2) - 1.2*AbsSqr(MassB)*Sqr(g1)*
      UNITMATRIX(2) - 6*AbsSqr(MassWB)*Sqr(g2)*UNITMATRIX(2) - 1.8*AbsSqr(
      MassBp)*Sqr(gN)*UNITMATRIX(2));


   return beta_mH1I2;
}

/**
 * Calculates the two-loop beta function of mH1I2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> CE6SSM_soft_parameters::calc_beta_mH1I2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = twoLoop*(-6*traceconjTKappaTpTKappa*((Lambda12).adjoint()
      *Lambda12) - 4*traceconjTLambda12TpTLambda12*((Lambda12).adjoint()*
      Lambda12) - 12*ms2*traceKappaAdjKappa*((Lambda12).adjoint()*Lambda12) - 6
      *traceKappaAdjKappaconjmDx2*((Lambda12).adjoint()*Lambda12) - 6*
      traceKappaconjmDxbar2AdjKappa*((Lambda12).adjoint()*Lambda12) - 8*ms2*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*Lambda12) - 4*
      traceLambda12AdjLambda12conjmH2I2*((Lambda12).adjoint()*Lambda12) - 4*
      tracemH1I2AdjLambda12Lambda12*((Lambda12).adjoint()*Lambda12) - 4*mHd2*
      AbsSqr(Lambdax)*((Lambda12).adjoint()*Lambda12) - 4*mHu2*AbsSqr(Lambdax)*
      ((Lambda12).adjoint()*Lambda12) - 8*ms2*AbsSqr(Lambdax)*((Lambda12)
      .adjoint()*Lambda12) - 4*AbsSqr(TLambdax)*((Lambda12).adjoint()*Lambda12)
      + 2*ms2*Sqr(gN)*((Lambda12).adjoint()*Lambda12) - 6*
      traceconjTKappaTpKappa*((Lambda12).adjoint()*TLambda12) - 4*
      traceconjTLambda12TpLambda12*((Lambda12).adjoint()*TLambda12) - 4*Conj(
      TLambdax)*Lambdax*((Lambda12).adjoint()*TLambda12) - 6*
      traceAdjKappaTKappa*((TLambda12).adjoint()*Lambda12) - 4*
      traceAdjLambda12TLambda12*((TLambda12).adjoint()*Lambda12) - 2*MassBp*Sqr
      (gN)*((TLambda12).adjoint()*Lambda12) - 4*Conj(Lambdax)*TLambdax*((
      TLambda12).adjoint()*Lambda12) - 6*traceKappaAdjKappa*((TLambda12)
      .adjoint()*TLambda12) - 4*traceLambda12AdjLambda12*((TLambda12).adjoint()
      *TLambda12) - 4*AbsSqr(Lambdax)*((TLambda12).adjoint()*TLambda12) + 2*Sqr
      (gN)*((TLambda12).adjoint()*TLambda12) - 3*traceKappaAdjKappa*(mH1I2*(
      Lambda12).adjoint()*Lambda12) - 2*traceLambda12AdjLambda12*(mH1I2*(
      Lambda12).adjoint()*Lambda12) - 2*AbsSqr(Lambdax)*(mH1I2*(Lambda12)
      .adjoint()*Lambda12) + Sqr(gN)*(mH1I2*(Lambda12).adjoint()*Lambda12) - 6*
      traceKappaAdjKappa*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - 4*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12
      ) - 4*AbsSqr(Lambdax)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) +
      2*Sqr(gN)*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - 3*
      traceKappaAdjKappa*((Lambda12).adjoint()*Lambda12*mH1I2) - 2*
      traceLambda12AdjLambda12*((Lambda12).adjoint()*Lambda12*mH1I2) - 2*AbsSqr
      (Lambdax)*((Lambda12).adjoint()*Lambda12*mH1I2) + Sqr(gN)*((Lambda12)
      .adjoint()*Lambda12*mH1I2) - 4*ms2*((Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*Lambda12*(
      TLambda12).adjoint()*TLambda12) - 2*((Lambda12).adjoint()*TLambda12*(
      TLambda12).adjoint()*Lambda12) - 2*((TLambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*TLambda12) - 2*((TLambda12).adjoint()*TLambda12*(
      Lambda12).adjoint()*Lambda12) - mH1I2*(Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12 - 2*((Lambda12).adjoint()*mH2I2.conjugate()*
      Lambda12*(Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*mH1I2*(Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*(Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - (Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12*mH1I2 + 6*Power(g2,4)*
      Tr22*UNITMATRIX(2) + 1.4696938456699067*g1*gN*Tr2U114*UNITMATRIX(2) +
      1.4696938456699067*g1*gN*Tr2U141*UNITMATRIX(2) - 3.0983866769659336*g1*
      Tr31*UNITMATRIX(2) - 3.794733192202055*gN*Tr34*UNITMATRIX(2) + 87*Power(
      g2,4)*AbsSqr(MassWB)*UNITMATRIX(2) + 1.2*Tr2U111*Sqr(g1)*UNITMATRIX(2) +
      3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(2) + 1.8*MassB*Conj(MassWB)
      *Sqr(g1)*Sqr(g2)*UNITMATRIX(2) + 1.8*Tr2U144*Sqr(gN)*UNITMATRIX(2) + 5.4*
      AbsSqr(MassWB)*Sqr(g2)*Sqr(gN)*UNITMATRIX(2) + 2.7*MassBp*Conj(MassWB)*
      Sqr(g2)*Sqr(gN)*UNITMATRIX(2) + 0.18*Conj(MassB)*Sqr(g1)*(198*MassB*Sqr(
      g1) + 10*(2*MassB + MassWB)*Sqr(g2) - (2*MassB + MassBp)*Sqr(gN))*
      UNITMATRIX(2) + 0.01*Conj(MassBp)*Sqr(gN)*(200*(2*MassBp*((Lambda12)
      .adjoint()*Lambda12) - (Lambda12).adjoint()*TLambda12) - 9*(2*(MassB + 2*
      MassBp)*Sqr(g1) - 3*(10*(2*MassBp + MassWB)*Sqr(g2) + 197*MassBp*Sqr(gN))
      )*UNITMATRIX(2)));


   return beta_mH1I2;
}

} // namespace flexiblesusy
