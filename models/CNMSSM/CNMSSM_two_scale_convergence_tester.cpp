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

#include "CNMSSM_two_scale_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

namespace flexiblesusy {

#define OLD1(p) ol.get_##p()
#define NEW1(p) ne.get_##p()

#define OLD(p,i) ol.get_##p()(i)
#define NEW(p,i) ne.get_##p()(i)

CNMSSM_convergence_tester<Two_scale>::CNMSSM_convergence_tester(CNMSSM<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CNMSSM<Two_scale> >(model, accuracy_goal)
{
}

CNMSSM_convergence_tester<Two_scale>::~CNMSSM_convergence_tester()
{
}

double CNMSSM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CNMSSM<Two_scale>& ol = get_last_iteration_model();
   const CNMSSM<Two_scale>& ne = get_model();

   double diff[37] = { 0 };

   diff[0] = MaxRelDiff(OLD1(MGlu),NEW1(MGlu));
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 1] = MaxRelDiff(OLD(MSd,i),NEW(MSd,i));
   }
   for (unsigned i = 0; i < 3; ++i) {
      diff[i + 7] = MaxRelDiff(OLD(MSv,i),NEW(MSv,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 10] = MaxRelDiff(OLD(MSu,i),NEW(MSu,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 16] = MaxRelDiff(OLD(MSe,i),NEW(MSe,i));
   }
   for (unsigned i = 0; i < 3; ++i) {
      diff[i + 22] = MaxRelDiff(OLD(Mhh,i),NEW(Mhh,i));
   }
   for (unsigned i = 1; i < 3; ++i) {
      diff[i + 25] = MaxRelDiff(OLD(MAh,i),NEW(MAh,i));
   }
   for (unsigned i = 1; i < 2; ++i) {
      diff[i + 28] = MaxRelDiff(OLD(MHpm,i),NEW(MHpm,i));
   }
   for (unsigned i = 0; i < 5; ++i) {
      diff[i + 30] = MaxRelDiff(OLD(MChi,i),NEW(MChi,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 35] = MaxRelDiff(OLD(MCha,i),NEW(MCha,i));
   }

   return *std::max_element(diff, diff + 37);

}

} // namespace flexiblesusy
