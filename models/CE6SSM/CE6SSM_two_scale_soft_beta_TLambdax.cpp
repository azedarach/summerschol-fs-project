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

// File generated at Fri 9 Jan 2015 15:02:14

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
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


   double beta_TLambdax;

   beta_TLambdax = oneOver16PiSqr*(0.2*Lambdax*(30*traceAdjKappaTKappa +
      20*traceAdjLambda12TLambda12 + 30*traceAdjYdTYd + 10*traceAdjYeTYe + 30*
      traceAdjYuTYu + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 19*MassBp*Sqr(gN))
      + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(
      g2) - 1.9*Sqr(gN))*TLambdax);


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;


   double beta_TLambdax;

   beta_TLambdax = twoLoop*(-0.02*Lambdax*(1188*Power(g1,4)*MassB + 3933*
      Power(gN,4)*MassBp + 3300*Power(g2,4)*MassWB + 1200*
      traceKappaAdjKappaTKappaAdjKappa + 800*
      traceLambda12AdjLambda12TLambda12AdjLambda12 + 1800*traceYdAdjYdTYdAdjYd
      + 600*traceYdAdjYuTYuAdjYd + 600*traceYeAdjYeTYeAdjYe + 600*
      traceYuAdjYdTYdAdjYu + 1800*traceYuAdjYuTYuAdjYu - 80*traceAdjKappaTKappa
      *Sqr(g1) - 120*traceAdjLambda12TLambda12*Sqr(g1) + 40*traceAdjYdTYd*Sqr(
      g1) - 120*traceAdjYeTYe*Sqr(g1) - 80*traceAdjYuTYu*Sqr(g1) + 80*MassB*
      traceKappaAdjKappa*Sqr(g1) + 120*MassB*traceLambda12AdjLambda12*Sqr(g1) +
      80*MassB*traceYuAdjYu*Sqr(g1) - 600*traceAdjLambda12TLambda12*Sqr(g2) +
      600*MassWB*traceLambda12AdjLambda12*Sqr(g2) + 180*MassB*Sqr(g1)*Sqr(g2) +
      180*MassWB*Sqr(g1)*Sqr(g2) - 1600*traceAdjKappaTKappa*Sqr(g3) - 1600*
      traceAdjYdTYd*Sqr(g3) - 1600*traceAdjYuTYu*Sqr(g3) + 1600*MassG*
      traceKappaAdjKappa*Sqr(g3) + 1600*MassG*traceYuAdjYu*Sqr(g3) + 180*
      traceAdjKappaTKappa*Sqr(gN) + 120*traceAdjLambda12TLambda12*Sqr(gN) + 60*
      traceAdjYdTYd*Sqr(gN) + 20*traceAdjYeTYe*Sqr(gN) + 30*traceAdjYuTYu*Sqr(
      gN) - 180*MassBp*traceKappaAdjKappa*Sqr(gN) - 120*MassBp*
      traceLambda12AdjLambda12*Sqr(gN) - 30*MassBp*traceYuAdjYu*Sqr(gN) + 39*
      MassB*Sqr(g1)*Sqr(gN) + 39*MassBp*Sqr(g1)*Sqr(gN) + 195*MassBp*Sqr(g2)*
      Sqr(gN) + 195*MassWB*Sqr(g2)*Sqr(gN) + 20*traceYeAdjYe*(6*MassB*Sqr(g1) -
      MassBp*Sqr(gN)) - 20*traceYdAdjYd*(2*MassB*Sqr(g1) - 80*MassG*Sqr(g3) +
      3*MassBp*Sqr(gN))) + (5.94*Power(g1,4) + 16.5*Power(g2,4) + 19.665*Power(
      gN,4) - 6*traceKappaAdjKappaKappaAdjKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 9*traceYdAdjYdYdAdjYd - 6*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - 9*traceYuAdjYuYuAdjYu + 0.8
      *traceKappaAdjKappa*Sqr(g1) + 1.2*traceLambda12AdjLambda12*Sqr(g1) + 0.8*
      traceYuAdjYu*Sqr(g1) + 6*traceLambda12AdjLambda12*Sqr(g2) + 1.8*Sqr(g1)*
      Sqr(g2) + 16*traceKappaAdjKappa*Sqr(g3) + 16*traceYuAdjYu*Sqr(g3) + 0.2*
      traceYeAdjYe*(6*Sqr(g1) - Sqr(gN)) - 1.8*traceKappaAdjKappa*Sqr(gN) - 1.2
      *traceLambda12AdjLambda12*Sqr(gN) - 0.3*traceYuAdjYu*Sqr(gN) + 0.39*Sqr(
      g1)*Sqr(gN) + 1.95*Sqr(g2)*Sqr(gN) - 0.2*traceYdAdjYd*(2*Sqr(g1) - 80*Sqr
      (g3) + 3*Sqr(gN)))*TLambdax - 50*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambdax
      - 0.1*AbsSqr(Lambdax)*(2*Lambdax*(60*traceAdjKappaTKappa + 40*
      traceAdjLambda12TLambda12 + 90*traceAdjYdTYd + 30*traceAdjYeTYe + 90*
      traceAdjYuTYu + 12*MassB*Sqr(g1) + 60*MassWB*Sqr(g2) + 13*MassBp*Sqr(gN))
      - 3*(-60*traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 90*
      traceYdAdjYd - 30*traceYeAdjYe - 90*traceYuAdjYu + 12*Sqr(g1) + 60*Sqr(g2
      ) + 13*Sqr(gN))*TLambdax));


   return beta_TLambdax;
}

} // namespace flexiblesusy
