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

// File generated at Fri 9 Jan 2015 15:02:28

#ifndef CE6SSM_INPUT_PARAMETERS_H
#define CE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct CE6SSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   double Azero;
   double LambdaInput;
   double KappaInput;
   double muPrimeInput;
   double BmuPrimeInput;
   double vSInput;
   double Lambda12Input;

   CE6SSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), Azero(0), LambdaInput(0), KappaInput(0),
   muPrimeInput(0), BmuPrimeInput(0), vSInput(0), Lambda12Input(0)

   {}
};

} // namespace flexiblesusy

#endif
