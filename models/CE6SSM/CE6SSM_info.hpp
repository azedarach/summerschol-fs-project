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

// File generated at Fri 9 Jan 2015 15:02:30

#ifndef CE6SSM_INFO_H
#define CE6SSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace CE6SSM_info {
   enum Particles : unsigned {VG, Glu, Fv, ChaP, VP, VZ, VZp, Sd, Sv, Su, Se,
      SDX, hh, Ah, Hpm, Chi, Cha, Fe, Fd, Fu, FDX, SHI0, SHIp, ChaI, ChiI, SSI0,
      FSI, SHp0, SHpp, ChiP, VWm, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, Kappa00, Kappa01
      , Kappa02, Kappa10, Kappa11, Kappa12, Kappa20, Kappa21, Kappa22, Lambda1200,
      Lambda1201, Lambda1210, Lambda1211, Lambdax, Yu00, Yu01, Yu02, Yu10, Yu11,
      Yu12, Yu20, Yu21, Yu22, MuPr, g1, g2, g3, gN, vd, vu, vs, TYd00, TYd01,
      TYd02, TYd10, TYd11, TYd12, TYd20, TYd21, TYd22, TYe00, TYe01, TYe02, TYe10,
      TYe11, TYe12, TYe20, TYe21, TYe22, TKappa00, TKappa01, TKappa02, TKappa10,
      TKappa11, TKappa12, TKappa20, TKappa21, TKappa22, TLambda1200, TLambda1201,
      TLambda1210, TLambda1211, TLambdax, TYu00, TYu01, TYu02, TYu10, TYu11, TYu12
      , TYu20, TYu21, TYu22, BMuPr, mq200, mq201, mq202, mq210, mq211, mq212,
      mq220, mq221, mq222, ml200, ml201, ml202, ml210, ml211, ml212, ml220, ml221,
      ml222, mHd2, mHu2, md200, md201, md202, md210, md211, md212, md220, md221,
      md222, mu200, mu201, mu202, mu210, mu211, mu212, mu220, mu221, mu222, me200,
      me201, me202, me210, me211, me212, me220, me221, me222, ms2, mH1I200,
      mH1I201, mH1I210, mH1I211, mH2I200, mH2I201, mH2I210, mH2I211, msI200,
      msI201, msI210, msI211, mDx200, mDx201, mDx202, mDx210, mDx211, mDx212,
      mDx220, mDx221, mDx222, mDxbar200, mDxbar201, mDxbar202, mDxbar210,
      mDxbar211, mDxbar212, mDxbar220, mDxbar221, mDxbar222, mHp2, mHpbar2, MassB,
      MassWB, MassG, MassBp, NUMBER_OF_PARAMETERS};

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const char* model_name;
   extern const bool is_low_energy_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
