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

// File generated at Fri 9 Jan 2015 15:02:31

#ifndef CE6SSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define CE6SSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "CE6SSM_high_scale_constraint.hpp"
#include "CE6SSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CE6SSM;

class Two_scale;

template<>
class CE6SSM_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CE6SSM_high_scale_constraint();
   CE6SSM_high_scale_constraint(const CE6SSM_input_parameters&);
   virtual ~CE6SSM_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   void initialize();
   void set_input_parameters(const CE6SSM_input_parameters&);
   void set_scale(double); ///< fix unification scale (0 = unfixed)

private:
   double scale;
   double initial_scale_guess;
   double fixed_scale; ///< fixed unification scale (0 = unfixed)
   CE6SSM<Two_scale>* model;
   CE6SSM_input_parameters inputPars;

   void update_scale();
};

} // namespace flexiblesusy

#endif
