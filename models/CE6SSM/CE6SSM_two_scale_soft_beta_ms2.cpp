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

// File generated at Fri 9 Jan 2015 15:02:24

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
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
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = oneOver16PiSqr*(1.5811388300841898*gN*Tr14 + 6*
      traceconjTKappaTpTKappa + 4*traceconjTLambda12TpTLambda12 + 6*ms2*
      traceKappaAdjKappa + 6*traceKappaAdjKappaconjmDx2 + 6*
      traceKappaconjmDxbar2AdjKappa + 4*ms2*traceLambda12AdjLambda12 + 4*
      traceLambda12AdjLambda12conjmH2I2 + 4*tracemH1I2AdjLambda12Lambda12 + 4*(
      mHd2 + mHu2 + ms2)*AbsSqr(Lambdax) + 4*AbsSqr(TLambdax) - 5*AbsSqr(MassBp
      )*Sqr(gN));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
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
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12;
   const double traceLambda12AdjTLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12;
   const double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
   const double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = 0.05*twoLoop*(3*Conj(MassBp)*Sqr(gN)*(24*
      traceAdjKappaTKappa + 16*traceAdjLambda12TLambda12 - 48*MassBp*
      traceKappaAdjKappa - 32*MassBp*traceLambda12AdjLambda12 + 1065*MassBp*Sqr
      (gN) + 16*Conj(Lambdax)*(-2*MassBp*Lambdax + TLambdax)) + 4*(
      31.622776601683796*gN*Tr34 - 60*traceKappaAdjKappaconjmDx2KappaAdjKappa -
      120*ms2*traceKappaAdjKappaKappaAdjKappa - 60*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 - 60*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa - 120*
      traceKappaAdjKappaTKappaAdjTKappa - 120*traceKappaAdjTKappaTKappaAdjKappa
      - 60*traceKappaconjmDxbar2AdjKappaKappaAdjKappa - 40*
      traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 - 80*ms2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 40*
      traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 - 80*
      traceLambda12AdjLambda12TLambda12AdjTLambda12 - 80*
      traceLambda12AdjTLambda12TLambda12AdjLambda12 - 80*
      tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 - 8*MassB*
      traceconjTKappaTpKappa*Sqr(g1) + 8*traceconjTKappaTpTKappa*Sqr(g1) - 12*
      MassB*traceconjTLambda12TpLambda12*Sqr(g1) + 12*
      traceconjTLambda12TpTLambda12*Sqr(g1) + 8*ms2*traceKappaAdjKappa*Sqr(g1)
      + 8*traceKappaAdjKappaconjmDx2*Sqr(g1) + 8*traceKappaconjmDxbar2AdjKappa*
      Sqr(g1) + 12*ms2*traceLambda12AdjLambda12*Sqr(g1) + 12*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g1) + 12*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1) + 16*traceKappaAdjKappa*AbsSqr(
      MassB)*Sqr(g1) + 24*traceLambda12AdjLambda12*AbsSqr(MassB)*Sqr(g1) - 8*
      traceAdjKappaTKappa*Conj(MassB)*Sqr(g1) - 12*traceAdjLambda12TLambda12*
      Conj(MassB)*Sqr(g1) - 60*MassWB*traceconjTLambda12TpLambda12*Sqr(g2) + 60
      *traceconjTLambda12TpTLambda12*Sqr(g2) + 60*ms2*traceLambda12AdjLambda12*
      Sqr(g2) + 60*traceLambda12AdjLambda12conjmH2I2*Sqr(g2) + 60*
      tracemH1I2AdjLambda12Lambda12*Sqr(g2) + 120*traceLambda12AdjLambda12*
      AbsSqr(MassWB)*Sqr(g2) - 60*traceAdjLambda12TLambda12*Conj(MassWB)*Sqr(g2
      ) - 160*MassG*traceconjTKappaTpKappa*Sqr(g3) + 160*
      traceconjTKappaTpTKappa*Sqr(g3) + 160*ms2*traceKappaAdjKappa*Sqr(g3) +
      160*traceKappaAdjKappaconjmDx2*Sqr(g3) + 160*
      traceKappaconjmDxbar2AdjKappa*Sqr(g3) + 320*traceKappaAdjKappa*AbsSqr(
      MassG)*Sqr(g3) - 160*traceAdjKappaTKappa*Conj(MassG)*Sqr(g3) + 25*Tr2U144
      *Sqr(gN) + 18*MassBp*traceconjTKappaTpKappa*Sqr(gN) - 18*
      traceconjTKappaTpTKappa*Sqr(gN) + 12*MassBp*traceconjTLambda12TpLambda12*
      Sqr(gN) - 12*traceconjTLambda12TpTLambda12*Sqr(gN) - 18*ms2*
      traceKappaAdjKappa*Sqr(gN) - 18*traceKappaAdjKappaconjmDx2*Sqr(gN) - 18*
      traceKappaconjmDxbar2AdjKappa*Sqr(gN) - 12*ms2*traceLambda12AdjLambda12*
      Sqr(gN) - 12*traceLambda12AdjLambda12conjmH2I2*Sqr(gN) - 12*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN) - 80*(mHd2 + mHu2 + ms2)*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 4*Conj(Lambdax)*(-15*traceconjTYdTpTYd*Lambdax -
      5*traceconjTYeTpTYe*Lambdax - 15*traceconjTYuTpTYu*Lambdax - 15*
      tracemd2YdAdjYd*Lambdax - 5*traceme2YeAdjYe*Lambdax - 5*traceml2AdjYeYe*
      Lambdax - 15*tracemq2AdjYdYd*Lambdax - 15*tracemq2AdjYuYu*Lambdax - 15*
      tracemu2YuAdjYu*Lambdax - 30*mHd2*traceYdAdjYd*Lambdax - 15*mHu2*
      traceYdAdjYd*Lambdax - 15*ms2*traceYdAdjYd*Lambdax - 10*mHd2*traceYeAdjYe
      *Lambdax - 5*mHu2*traceYeAdjYe*Lambdax - 5*ms2*traceYeAdjYe*Lambdax - 15*
      mHd2*traceYuAdjYu*Lambdax - 30*mHu2*traceYuAdjYu*Lambdax - 15*ms2*
      traceYuAdjYu*Lambdax - 40*AbsSqr(TLambdax)*Lambdax + 3*mHd2*Lambdax*Sqr(
      g1) + 3*mHu2*Lambdax*Sqr(g1) + 3*ms2*Lambdax*Sqr(g1) + 15*mHd2*Lambdax*
      Sqr(g2) + 15*mHu2*Lambdax*Sqr(g2) + 15*ms2*Lambdax*Sqr(g2) - 3*mHd2*
      Lambdax*Sqr(gN) - 3*mHu2*Lambdax*Sqr(gN) - 3*ms2*Lambdax*Sqr(gN) + 3*Conj
      (MassB)*Sqr(g1)*(2*MassB*Lambdax - TLambdax) + 15*Conj(MassWB)*Sqr(g2)*(2
      *MassWB*Lambdax - TLambdax) - 15*traceconjTYdTpYd*TLambdax - 5*
      traceconjTYeTpYe*TLambdax - 15*traceconjTYuTpYu*TLambdax) - 4*Conj(
      TLambdax)*(Lambdax*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu
      + 3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2) - 3*MassBp*Sqr(gN)) + (15*
      traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2)
      + 3*Sqr(gN))*TLambdax)));


   return beta_ms2;
}

} // namespace flexiblesusy
