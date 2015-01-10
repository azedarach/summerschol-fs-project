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

// File generated at Fri 9 Jan 2015 15:06:53

#include "CNMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappa.
 *
 * @return one-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_TKappa_one_loop(const Soft_traces& soft_traces) const
{


   double beta_TKappa;

   beta_TKappa = 6*oneOver16PiSqr*(3*AbsSqr(Kappa)*TKappa + Conj(Lambdax)
      *(Lambdax*TKappa + 2*Kappa*TLambdax));


   return beta_TKappa;
}

/**
 * Calculates the two-loop beta function of TKappa.
 *
 * @return two-loop beta function
 */
double CNMSSM_soft_parameters::calc_beta_TKappa_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TKappa;

   beta_TKappa = -1.2*twoLoop*(100*Sqr(Conj(Kappa))*Sqr(Kappa)*TKappa +
      10*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TKappa + 4*Kappa*TLambdax) + Conj(
      Lambdax)*(Lambdax*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu +
      60*AbsSqr(Kappa) - 3*Sqr(g1) - 15*Sqr(g2))*TKappa + 2*Kappa*(Lambdax*(15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*Sqr(g1) + 15
      *MassWB*Sqr(g2)) + (15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu +
      20*AbsSqr(Kappa) - 3*Sqr(g1) - 15*Sqr(g2))*TLambdax)));


   return beta_TKappa;
}

} // namespace flexiblesusy
