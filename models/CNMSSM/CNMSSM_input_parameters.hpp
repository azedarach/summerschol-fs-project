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

// File generated at Fri 9 Jan 2015 15:07:01

#ifndef CNMSSM_INPUT_PARAMETERS_H
#define CNMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct CNMSSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   int SignvS;
   double Azero;
   double LambdaInput;

   CNMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), SignvS(1), Azero(0), LambdaInput(0)

   {}
};

} // namespace flexiblesusy

#endif
