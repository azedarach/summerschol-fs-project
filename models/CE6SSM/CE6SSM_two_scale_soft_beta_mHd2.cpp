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

// File generated at Fri 9 Jan 2015 15:02:19

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 -
      0.9486832980505138*gN*Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe +
      6*tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*
      tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*
      AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*
      AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) -
      1.8*AbsSqr(MassBp)*Sqr(gN));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   beta_mHd2 = twoLoop*(6*Power(g2,4)*Tr22 + 1.4696938456699067*g1*gN*
      Tr2U114 + 1.4696938456699067*g1*gN*Tr2U141 - 3.0983866769659336*g1*Tr31 -
      3.794733192202055*gN*Tr34 - 36*tracemd2YdAdjYdYdAdjYd - 6*
      tracemd2YdAdjYuYuAdjYd - 12*traceme2YeAdjYeYeAdjYe - 12*
      traceml2AdjYeYeAdjYeYe - 36*tracemq2AdjYdYdAdjYdYd - 6*
      tracemq2AdjYdYdAdjYuYu - 6*tracemq2AdjYuYuAdjYdYd - 6*
      tracemu2YuAdjYdYdAdjYu - 36*traceYdAdjTYdTYdAdjYd - 6*
      traceYdAdjTYuTYuAdjYd - 36*traceYdAdjYdTYdAdjTYd - 36*mHd2*
      traceYdAdjYdYdAdjYd - 6*traceYdAdjYuTYuAdjTYd - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 12*
      traceYeAdjTYeTYeAdjYe - 12*traceYeAdjYeTYeAdjTYe - 12*mHd2*
      traceYeAdjYeYeAdjYe - 6*traceYuAdjTYdTYdAdjYu - 6*traceYuAdjYdTYdAdjTYu +
      87*Power(g2,4)*AbsSqr(MassWB) - 6*traceconjTKappaTpTKappa*AbsSqr(Lambdax
      ) - 4*traceconjTLambda12TpTLambda12*AbsSqr(Lambdax) - 6*traceconjTYuTpTYu
      *AbsSqr(Lambdax) - 6*mHd2*traceKappaAdjKappa*AbsSqr(Lambdax) - 6*mHu2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 12*ms2*traceKappaAdjKappa*AbsSqr(
      Lambdax) - 6*traceKappaAdjKappaconjmDx2*AbsSqr(Lambdax) - 6*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 4*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 4*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 8*ms2*traceLambda12AdjLambda12
      *AbsSqr(Lambdax) - 4*traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) -
      4*tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 6*tracemq2AdjYuYu*
      AbsSqr(Lambdax) - 6*tracemu2YuAdjYu*AbsSqr(Lambdax) - 6*mHd2*traceYuAdjYu
      *AbsSqr(Lambdax) - 12*mHu2*traceYuAdjYu*AbsSqr(Lambdax) - 6*ms2*
      traceYuAdjYu*AbsSqr(Lambdax) - 6*traceKappaAdjKappa*AbsSqr(TLambdax) - 4*
      traceLambda12AdjLambda12*AbsSqr(TLambdax) - 6*traceYuAdjYu*AbsSqr(
      TLambdax) - 24*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 6*traceAdjKappaTKappa*
      Conj(TLambdax)*Lambdax - 4*traceAdjLambda12TLambda12*Conj(TLambdax)*
      Lambdax - 6*traceAdjYuTYu*Conj(TLambdax)*Lambdax + 1.2*Tr2U111*Sqr(g1) -
      0.8*traceconjTYdTpTYd*Sqr(g1) + 0.8*MassB*traceconjTYdTpYd*Sqr(g1) + 2.4*
      traceconjTYeTpTYe*Sqr(g1) - 2.4*MassB*traceconjTYeTpYe*Sqr(g1) - 0.8*
      tracemd2YdAdjYd*Sqr(g1) + 2.4*traceme2YeAdjYe*Sqr(g1) + 2.4*
      traceml2AdjYeYe*Sqr(g1) - 0.8*tracemq2AdjYdYd*Sqr(g1) - 0.8*mHd2*
      traceYdAdjYd*Sqr(g1) + 2.4*mHd2*traceYeAdjYe*Sqr(g1) + 3.6*AbsSqr(MassWB)
      *Sqr(g1)*Sqr(g2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 32*
      traceconjTYdTpTYd*Sqr(g3) - 32*MassG*traceconjTYdTpYd*Sqr(g3) + 32*
      tracemd2YdAdjYd*Sqr(g3) + 32*tracemq2AdjYdYd*Sqr(g3) + 32*mHd2*
      traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) + 1.8*Tr2U144*Sqr(gN) - 1.2*
      traceconjTYdTpTYd*Sqr(gN) + 1.2*MassBp*traceconjTYdTpYd*Sqr(gN) - 0.4*
      traceconjTYeTpTYe*Sqr(gN) + 0.4*MassBp*traceconjTYeTpYe*Sqr(gN) - 1.2*
      tracemd2YdAdjYd*Sqr(gN) - 0.4*traceme2YeAdjYe*Sqr(gN) - 0.4*
      traceml2AdjYeYe*Sqr(gN) - 1.2*tracemq2AdjYdYd*Sqr(gN) - 1.2*mHd2*
      traceYdAdjYd*Sqr(gN) - 0.4*mHd2*traceYeAdjYe*Sqr(gN) + 2*mHd2*AbsSqr(
      Lambdax)*Sqr(gN) + 2*mHu2*AbsSqr(Lambdax)*Sqr(gN) + 2*ms2*AbsSqr(Lambdax)
      *Sqr(gN) + 2*AbsSqr(TLambdax)*Sqr(gN) - 2*MassBp*Conj(TLambdax)*Lambdax*
      Sqr(gN) + 5.4*AbsSqr(MassWB)*Sqr(g2)*Sqr(gN) + 2.7*MassBp*Conj(MassWB)*
      Sqr(g2)*Sqr(gN) + 0.02*Conj(MassB)*Sqr(g1)*(40*traceAdjYdTYd - 120*
      traceAdjYeTYe - 80*MassB*traceYdAdjYd + 240*MassB*traceYeAdjYe + 1782*
      MassB*Sqr(g1) + 180*MassB*Sqr(g2) + 90*MassWB*Sqr(g2) - 18*MassB*Sqr(gN)
      - 9*MassBp*Sqr(gN)) - 12*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*mHu2*
      Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      + 0.01*Conj(MassBp)*Sqr(gN)*(120*traceAdjYdTYd + 40*traceAdjYeTYe - 240*
      MassBp*traceYdAdjYd - 80*MassBp*traceYeAdjYe - 18*MassB*Sqr(g1) - 36*
      MassBp*Sqr(g1) + 540*MassBp*Sqr(g2) + 270*MassWB*Sqr(g2) + 5319*MassBp*
      Sqr(gN) + 200*Conj(Lambdax)*(2*MassBp*Lambdax - TLambdax)) - 6*
      traceconjTKappaTpKappa*Conj(Lambdax)*TLambdax - 4*
      traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax - 6*traceconjTYuTpYu*
      Conj(Lambdax)*TLambdax);


   return beta_mHd2;
}

} // namespace flexiblesusy
