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

// File generated at Fri 9 Jan 2015 15:02:15

#include "CE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of BMuPr.
 *
 * @return one-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = 0.2*oneOver16PiSqr*(-(BMuPr*(3*Sqr(g1) + 15*Sqr(g2) + 2*
      Sqr(gN))) + 2*MuPr*(3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2) + 2*MassBp*Sqr(gN
      )));


   return beta_BMuPr;
}

/**
 * Calculates the two-loop beta function of BMuPr.
 *
 * @return two-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_two_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = -0.06*twoLoop*(-(BMuPr*(99*Power(g1,4) + 275*Power(g2,4)
      + 64*Power(gN,4) + 20*Sqr(g2)*Sqr(gN) + Sqr(g1)*(30*Sqr(g2) + 4*Sqr(gN)))
      ) + 4*MuPr*(99*Power(g1,4)*MassB + 64*Power(gN,4)*MassBp + 275*Power(g2,4
      )*MassWB + 10*(MassBp + MassWB)*Sqr(g2)*Sqr(gN) + Sqr(g1)*(15*(MassB +
      MassWB)*Sqr(g2) + 2*(MassB + MassBp)*Sqr(gN))));


   return beta_BMuPr;
}

} // namespace flexiblesusy
