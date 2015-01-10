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

// File generated at Fri 9 Jan 2015 15:02:29

#include "CE6SSM_two_scale_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

namespace flexiblesusy {

#define OLD1(p) ol.get_##p()
#define NEW1(p) ne.get_##p()

#define OLD(p,i) ol.get_##p()(i)
#define NEW(p,i) ne.get_##p()(i)

CE6SSM_convergence_tester<Two_scale>::CE6SSM_convergence_tester(CE6SSM<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CE6SSM<Two_scale> >(model, accuracy_goal)
{
}

CE6SSM_convergence_tester<Two_scale>::~CE6SSM_convergence_tester()
{
}

double CE6SSM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CE6SSM<Two_scale>& ol = get_last_iteration_model();
   const CE6SSM<Two_scale>& ne = get_model();

   double diff[73] = { 0 };

   diff[0] = MaxRelDiff(OLD1(MGlu),NEW1(MGlu));
   diff[1] = MaxRelDiff(OLD1(MChaP),NEW1(MChaP));
   diff[2] = MaxRelDiff(OLD1(MVZp),NEW1(MVZp));
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 3] = MaxRelDiff(OLD(MSd,i),NEW(MSd,i));
   }
   for (unsigned i = 0; i < 3; ++i) {
      diff[i + 9] = MaxRelDiff(OLD(MSv,i),NEW(MSv,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 12] = MaxRelDiff(OLD(MSu,i),NEW(MSu,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 18] = MaxRelDiff(OLD(MSe,i),NEW(MSe,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 24] = MaxRelDiff(OLD(MSDX,i),NEW(MSDX,i));
   }
   for (unsigned i = 0; i < 3; ++i) {
      diff[i + 30] = MaxRelDiff(OLD(Mhh,i),NEW(Mhh,i));
   }
   for (unsigned i = 2; i < 3; ++i) {
      diff[i + 33] = MaxRelDiff(OLD(MAh,i),NEW(MAh,i));
   }
   for (unsigned i = 1; i < 2; ++i) {
      diff[i + 36] = MaxRelDiff(OLD(MHpm,i),NEW(MHpm,i));
   }
   for (unsigned i = 0; i < 6; ++i) {
      diff[i + 38] = MaxRelDiff(OLD(MChi,i),NEW(MChi,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 44] = MaxRelDiff(OLD(MCha,i),NEW(MCha,i));
   }
   for (unsigned i = 0; i < 3; ++i) {
      diff[i + 46] = MaxRelDiff(OLD(MFDX,i),NEW(MFDX,i));
   }
   for (unsigned i = 0; i < 4; ++i) {
      diff[i + 49] = MaxRelDiff(OLD(MSHI0,i),NEW(MSHI0,i));
   }
   for (unsigned i = 0; i < 4; ++i) {
      diff[i + 53] = MaxRelDiff(OLD(MSHIp,i),NEW(MSHIp,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 57] = MaxRelDiff(OLD(MChaI,i),NEW(MChaI,i));
   }
   for (unsigned i = 0; i < 4; ++i) {
      diff[i + 59] = MaxRelDiff(OLD(MChiI,i),NEW(MChiI,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 63] = MaxRelDiff(OLD(MSSI0,i),NEW(MSSI0,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 65] = MaxRelDiff(OLD(MFSI,i),NEW(MFSI,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 67] = MaxRelDiff(OLD(MSHp0,i),NEW(MSHp0,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 69] = MaxRelDiff(OLD(MSHpp,i),NEW(MSHpp,i));
   }
   for (unsigned i = 0; i < 2; ++i) {
      diff[i + 71] = MaxRelDiff(OLD(MChiP,i),NEW(MChiP,i));
   }

   return *std::max_element(diff, diff + 73);

}

} // namespace flexiblesusy
