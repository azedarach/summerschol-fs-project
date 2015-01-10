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

// File generated at Fri 9 Jan 2015 15:01:31

#include "CE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME CE6SSM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

CE6SSM_susy_parameters::CE6SSM_susy_parameters(const CE6SSM_input_parameters& input_)
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), Kappa(Eigen::Matrix<double,3,3>::Zero()), Lambda12(Eigen::Matrix<double,
   2,2>::Zero()), Lambdax(0), Yu(Eigen::Matrix<double,3,3>::Zero()), MuPr(0),
   g1(0), g2(0), g3(0), gN(0), vd(0), vu(0), vs(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

CE6SSM_susy_parameters::CE6SSM_susy_parameters(
   double scale_, double loops_, double thresholds_,
   const CE6SSM_input_parameters& input_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , const Eigen::Matrix<double,3,3>& Kappa_, const Eigen::Matrix<double,2,2>&
   Lambda12_, double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, double
   MuPr_, double g1_, double g2_, double g3_, double gN_, double vd_, double
   vu_, double vs_

)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), Kappa(Kappa_), Lambda12(Lambda12_), Lambdax(Lambdax_),
   Yu(Yu_), MuPr(MuPr_), g1(g1_), g2(g2_), g3(g3_), gN(gN_), vd(vd_), vu(vu_),
   vs(vs_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd CE6SSM_susy_parameters::beta() const
{
   return calc_beta().get();
}

CE6SSM_susy_parameters CE6SSM_susy_parameters::calc_beta() const
{
   Susy_traces susy_traces;
   calc_susy_traces(susy_traces);

   Eigen::Matrix<double,3,3> beta_Yd(calc_beta_Yd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Ye(calc_beta_Ye_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Kappa(calc_beta_Kappa_one_loop(TRACE_STRUCT))
      ;
   Eigen::Matrix<double,2,2> beta_Lambda12(calc_beta_Lambda12_one_loop(
      TRACE_STRUCT));
   double beta_Lambdax(calc_beta_Lambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yu(calc_beta_Yu_one_loop(TRACE_STRUCT));
   double beta_MuPr(calc_beta_MuPr_one_loop(TRACE_STRUCT));
   double beta_g1(calc_beta_g1_one_loop(TRACE_STRUCT));
   double beta_g2(calc_beta_g2_one_loop(TRACE_STRUCT));
   double beta_g3(calc_beta_g3_one_loop(TRACE_STRUCT));
   double beta_gN(calc_beta_gN_one_loop(TRACE_STRUCT));
   double beta_vd(calc_beta_vd_one_loop(TRACE_STRUCT));
   double beta_vu(calc_beta_vu_one_loop(TRACE_STRUCT));
   double beta_vs(calc_beta_vs_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
      beta_Kappa += calc_beta_Kappa_two_loop(TRACE_STRUCT);
      beta_Lambda12 += calc_beta_Lambda12_two_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
      beta_MuPr += calc_beta_MuPr_two_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
      beta_gN += calc_beta_gN_two_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_two_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_two_loop(TRACE_STRUCT);
      beta_vs += calc_beta_vs_two_loop(TRACE_STRUCT);

   }


   return CE6SSM_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_Yd, beta_Ye, beta_Kappa, beta_Lambda12, beta_Lambdax, beta_Yu, beta_MuPr, beta_g1, beta_g2, beta_g3, beta_gN, beta_vd, beta_vu, beta_vs);
}

void CE6SSM_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Kappa = Eigen::Matrix<double,3,3>::Zero();
   Lambda12 = Eigen::Matrix<double,2,2>::Zero();
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   MuPr = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   gN = 0.;
   vd = 0.;
   vu = 0.;
   vs = 0.;

}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu -
      0.016666666666666666*(2*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3) + 3*Sqr(gN))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(0.2*(-10*(Yd.adjoint()*Yd*Yd.adjoint()*Yd +
         Yu.adjoint()*Yu*Yu.adjoint()*Yu) + Yd.adjoint()*Yd*(-5*AbsSqr(Lambdax)
         + 2*Sqr(g1) + 3*Sqr(gN) - 15*(Yd*Yd.adjoint()).trace() - 5*(Ye*
         Ye.adjoint()).trace()) + Yu.adjoint()*Yu*(-5*AbsSqr(Lambdax) + 4*Sqr(
         g1) + Sqr(gN) - 15*(Yu*Yu.adjoint()).trace())) + 0.0002777777777777778
         *(1156*Power(g1,4) + 29700*Power(g2,4) + 25600*Power(g3,4) + 1701*
         Power(gN,4) + 4*Sqr(g1)*(90*Sqr(g2) + 160*Sqr(g3) - 33*Sqr(gN)) + 960*
         Sqr(g3)*Sqr(gN) + 180*Sqr(g2)*(160*Sqr(g3) + 3*Sqr(gN)))*UNITMATRIX(3)
         );
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*(Ye.adjoint()*Ye - 0.1*(3*Sqr(g1) + 15*Sqr(g2
      ) + 2*Sqr(gN))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += 0.01*twoLoop*(10*(-20*(Ye.adjoint()*Ye*Ye.adjoint()*
         Ye) + Ye.adjoint()*Ye*(-10*AbsSqr(Lambdax) + 12*Sqr(g1) + 3*Sqr(gN) -
         30*(Yd*Yd.adjoint()).trace() - 10*(Ye*Ye.adjoint()).trace())) + 3*(99*
         Power(g1,4) + 275*Power(g2,4) + 64*Power(gN,4) + 20*Sqr(g2)*Sqr(gN) +
         6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN)))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2)
      - 0.45*Sqr(gN) + 3*(Yd*Yd.adjoint()).trace() + (Ye*Ye.adjoint()).trace())
      ;

   if (get_loops() > 1) {
      anomDim += twoLoop*(2.97*Power(g1,4) + 8.25*Power(g2,4) + 4.4325
         *Power(gN,4) + 0.9*Sqr(g1)*Sqr(g2) - 0.09*Sqr(g1)*Sqr(gN) + 1.35*Sqr(
         g2)*Sqr(gN) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 0.2*(2*Sqr(g1) - 80*
         Sqr(g3) + 3*Sqr(gN))*(Yd*Yd.adjoint()).trace() + 1.2*Sqr(g1)*(Ye*
         Ye.adjoint()).trace() - 0.2*Sqr(gN)*(Ye*Ye.adjoint()).trace() + AbsSqr
         (Lambdax)*(Sqr(gN) - 3*(Yu*Yu.adjoint()).trace() - 3*(Kappa*(Kappa)
         .adjoint()).trace() - 2*(Lambda12*(Lambda12).adjoint()).trace()) - 9*(
         Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*
         Yd.adjoint()).trace() - 3*(Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace());
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2)
      - 0.2*Sqr(gN) + 3*(Yu*Yu.adjoint()).trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(-3*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.5*
         AbsSqr(Lambdax)*(3*Sqr(gN) - 6*(Yd*Yd.adjoint()).trace() - 2*(Ye*
         Ye.adjoint()).trace() - 6*(Kappa*(Kappa).adjoint()).trace() - 4*(
         Lambda12*(Lambda12).adjoint()).trace()) + 0.01*(10*(8*Sqr(g1) + 160*
         Sqr(g3) - 3*Sqr(gN))*(Yu*Yu.adjoint()).trace() + 3*(99*Power(g1,4) +
         275*Power(g2,4) + 64*Power(gN,4) + 30*Sqr(g1)*Sqr(g2) + 12*Sqr(g1)*Sqr
         (gN) + 20*Sqr(g2)*Sqr(gN) - 100*(Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace() - 300*(Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace())));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = 0.06666666666666667*oneOver16PiSqr*(30*(Yd.conjugate()*
      Yd.transpose()) - (2*Sqr(g1) + 40*Sqr(g3) + 3*Sqr(gN))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Yd.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yd.transpose() + Yd.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yd.transpose()) + Yd.conjugate()*Yd.transpose()*(-2*
         AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6*Sqr(g2) + 0.6*Sqr(gN) - 6*(Yd*
         Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) +
         0.017777777777777778*(73*Power(g1,4) + 400*Power(g3,4) + 108*Power(gN,
         4) + Sqr(g1)*(40*Sqr(g3) - 6*Sqr(gN)) + 60*Sqr(g3)*Sqr(gN))*UNITMATRIX
         (3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.016666666666666666*(32*Sqr(g1) + 160*Sqr(g3) + 3*Sqr(gN))*UNITMATRIX(3)
      );

   if (get_loops() > 1) {
      anomDim += twoLoop*(-0.4*(5*(Yu.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yu.transpose() + Yu.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yu.transpose()) + Yu.conjugate()*Yu.transpose()*(5*
         AbsSqr(Lambdax) + Sqr(g1) - 15*Sqr(g2) - Sqr(gN) + 15*(Yu*Yu.adjoint()
         ).trace())) + 0.0002777777777777778*(19456*Power(g1,4) + 25600*Power(
         g3,4) + 1701*Power(gN,4) + 960*Sqr(g3)*Sqr(gN) + 256*Sqr(g1)*(40*Sqr(
         g3) + 3*Sqr(gN)))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*(2*(Ye.conjugate()*Ye.transpose()) - 0.05*(24
      *Sqr(g1) + Sqr(gN))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-0.4*(5*(Ye.conjugate()*Ye.transpose()*
         Ye.conjugate()*Ye.transpose()) + Ye.conjugate()*Ye.transpose()*(5*
         AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr(g2) - 3*Sqr(gN) + 15*(Yd*
         Yd.adjoint()).trace() + 5*(Ye*Ye.adjoint()).trace())) + 0.0075*(1728*
         Power(g1,4) + 63*Power(gN,4) - 16*Sqr(g1)*Sqr(gN))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;

   anomDim = oneOver16PiSqr*(2*AbsSqr(Lambdax) - 1.25*Sqr(gN) + 3*(Kappa*
      (Kappa).adjoint()).trace() + 2*(Lambda12*(Lambda12).adjoint()).trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(13.3125*Power(gN,4) - 4*Sqr(Conj(Lambdax))*
         Sqr(Lambdax) - 0.4*AbsSqr(Lambdax)*(15*(Yd*Yd.adjoint()).trace() + 5*(
         Ye*Ye.adjoint()).trace() - 3*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN) - 5*(Yu*
         Yu.adjoint()).trace())) + 0.2*(4*Sqr(g1) + 80*Sqr(g3) - 9*Sqr(gN))*(
         Kappa*(Kappa).adjoint()).trace() + 1.2*Sqr(g1)*(Lambda12*(Lambda12)
         .adjoint()).trace() + 6*Sqr(g2)*(Lambda12*(Lambda12).adjoint()).trace(
         ) - 1.2*Sqr(gN)*(Lambda12*(Lambda12).adjoint()).trace() - 6*(Kappa*(
         Kappa).adjoint()*Kappa*(Kappa).adjoint()).trace() - 4*(Lambda12*(
         Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace());
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH1ISH1I() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;

   anomDim = oneOver16PiSqr*((Lambda12).adjoint()*Lambda12 - 0.15*(2*Sqr(
      g1) + 10*Sqr(g2) + 3*Sqr(gN))*UNITMATRIX(2));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-((Lambda12).adjoint()*Lambda12*(Lambda12)
         .adjoint()*Lambda12) + (Lambda12).adjoint()*Lambda12*(-2*AbsSqr(
         Lambdax) + Sqr(gN) - 3*(Kappa*(Kappa).adjoint()).trace() - 2*(Lambda12
         *(Lambda12).adjoint()).trace()) + 0.0075*(396*Power(g1,4) + 1100*Power
         (g2,4) + 591*Power(gN,4) + 12*Sqr(g1)*(10*Sqr(g2) - Sqr(gN)) + 180*Sqr
         (g2)*Sqr(gN))*UNITMATRIX(2));
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH2ISH2I() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;

   anomDim = oneOver16PiSqr*(Lambda12.conjugate()*(Lambda12).transpose()
      - 0.1*(3*Sqr(g1) + 15*Sqr(g2) + 2*Sqr(gN))*UNITMATRIX(2));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(Lambda12.conjugate()*(Lambda12).transpose(
         )*Lambda12.conjugate()*(Lambda12).transpose()) + Lambda12.conjugate()*
         (Lambda12).transpose()*(-2*AbsSqr(Lambdax) + 1.5*Sqr(gN) - 3*(Kappa*(
         Kappa).adjoint()).trace() - 2*(Lambda12*(Lambda12).adjoint()).trace())
         + 0.03*(99*Power(g1,4) + 275*Power(g2,4) + 64*Power(gN,4) + 20*Sqr(g2
         )*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN)))*UNITMATRIX(2));
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SsIRSsIR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;

   anomDim = -1.25*oneOver16PiSqr*Sqr(gN)*UNITMATRIX(2);

   if (get_loops() > 1) {
      anomDim += 13.3125*Power(gN,4)*twoLoop*UNITMATRIX(2);
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxLSDxL() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*(Kappa.conjugate()*(Kappa).transpose() -
      0.06666666666666667*(2*Sqr(g1) + 40*Sqr(g3) + 3*Sqr(gN))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(Kappa.conjugate()*(Kappa).transpose()*
         Kappa.conjugate()*(Kappa).transpose()) + Kappa.conjugate()*(Kappa)
         .transpose()*(-2*AbsSqr(Lambdax) + 1.5*Sqr(gN) - 3*(Kappa*(Kappa)
         .adjoint()).trace() - 2*(Lambda12*(Lambda12).adjoint()).trace()) +
         0.017777777777777778*(73*Power(g1,4) + 400*Power(g3,4) + 108*Power(gN,
         4) + Sqr(g1)*(40*Sqr(g3) - 6*Sqr(gN)) + 60*Sqr(g3)*Sqr(gN))*UNITMATRIX
         (3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxbarRSDxbarR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = oneOver16PiSqr*((Kappa).adjoint()*Kappa -
      0.016666666666666666*(8*Sqr(g1) + 160*Sqr(g3) + 27*Sqr(gN))*UNITMATRIX(3)
      );

   if (get_loops() > 1) {
      anomDim += twoLoop*(-((Kappa).adjoint()*Kappa*(Kappa).adjoint()*
         Kappa) + (Kappa).adjoint()*Kappa*(-2*AbsSqr(Lambdax) + Sqr(gN) - 3*(
         Kappa*(Kappa).adjoint()).trace() - 2*(Lambda12*(Lambda12).adjoint())
         .trace()) + 0.0002777777777777778*(4672*Power(g1,4) + 25600*Power(g3,4
         ) + 15957*Power(gN,4) + 8640*Sqr(g3)*Sqr(gN) + 16*Sqr(g1)*(160*Sqr(g3)
         + 81*Sqr(gN)))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SHpSHp() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;

   anomDim = 0.1*oneOver16PiSqr*(-3*Sqr(g1) - 15*Sqr(g2) - 2*Sqr(gN));

   if (get_loops() > 1) {
      anomDim += 0.03*twoLoop*(99*Power(g1,4) + 275*Power(g2,4) + 64*
         Power(gN,4) + 20*Sqr(g2)*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN)));
   }

   return anomDim;
}

double CLASSNAME::get_SHpbarSHpbar() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;

   anomDim = 0.1*oneOver16PiSqr*(-3*Sqr(g1) - 15*Sqr(g2) - 2*Sqr(gN));

   if (get_loops() > 1) {
      anomDim += 0.03*twoLoop*(99*Power(g1,4) + 275*Power(g2,4) + 64*
         Power(gN,4) + 20*Sqr(g2)*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN)));
   }

   return anomDim;
}


const Eigen::ArrayXd CE6SSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = Ye(0,0);
   pars(10) = Ye(0,1);
   pars(11) = Ye(0,2);
   pars(12) = Ye(1,0);
   pars(13) = Ye(1,1);
   pars(14) = Ye(1,2);
   pars(15) = Ye(2,0);
   pars(16) = Ye(2,1);
   pars(17) = Ye(2,2);
   pars(18) = Kappa(0,0);
   pars(19) = Kappa(0,1);
   pars(20) = Kappa(0,2);
   pars(21) = Kappa(1,0);
   pars(22) = Kappa(1,1);
   pars(23) = Kappa(1,2);
   pars(24) = Kappa(2,0);
   pars(25) = Kappa(2,1);
   pars(26) = Kappa(2,2);
   pars(27) = Lambda12(0,0);
   pars(28) = Lambda12(0,1);
   pars(29) = Lambda12(1,0);
   pars(30) = Lambda12(1,1);
   pars(31) = Lambdax;
   pars(32) = Yu(0,0);
   pars(33) = Yu(0,1);
   pars(34) = Yu(0,2);
   pars(35) = Yu(1,0);
   pars(36) = Yu(1,1);
   pars(37) = Yu(1,2);
   pars(38) = Yu(2,0);
   pars(39) = Yu(2,1);
   pars(40) = Yu(2,2);
   pars(41) = MuPr;
   pars(42) = g1;
   pars(43) = g2;
   pars(44) = g3;
   pars(45) = gN;
   pars(46) = vd;
   pars(47) = vu;
   pars(48) = vs;


   return pars;
}

void CE6SSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Kappa = " << Kappa << '\n';
   ostr << "Lambda12 = " << Lambda12 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "MuPr = " << MuPr << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "gN = " << gN << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vs = " << vs << '\n';

}

void CE6SSM_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   Yd(0,0) = pars(0);
   Yd(0,1) = pars(1);
   Yd(0,2) = pars(2);
   Yd(1,0) = pars(3);
   Yd(1,1) = pars(4);
   Yd(1,2) = pars(5);
   Yd(2,0) = pars(6);
   Yd(2,1) = pars(7);
   Yd(2,2) = pars(8);
   Ye(0,0) = pars(9);
   Ye(0,1) = pars(10);
   Ye(0,2) = pars(11);
   Ye(1,0) = pars(12);
   Ye(1,1) = pars(13);
   Ye(1,2) = pars(14);
   Ye(2,0) = pars(15);
   Ye(2,1) = pars(16);
   Ye(2,2) = pars(17);
   Kappa(0,0) = pars(18);
   Kappa(0,1) = pars(19);
   Kappa(0,2) = pars(20);
   Kappa(1,0) = pars(21);
   Kappa(1,1) = pars(22);
   Kappa(1,2) = pars(23);
   Kappa(2,0) = pars(24);
   Kappa(2,1) = pars(25);
   Kappa(2,2) = pars(26);
   Lambda12(0,0) = pars(27);
   Lambda12(0,1) = pars(28);
   Lambda12(1,0) = pars(29);
   Lambda12(1,1) = pars(30);
   Lambdax = pars(31);
   Yu(0,0) = pars(32);
   Yu(0,1) = pars(33);
   Yu(0,2) = pars(34);
   Yu(1,0) = pars(35);
   Yu(1,1) = pars(36);
   Yu(1,2) = pars(37);
   Yu(2,0) = pars(38);
   Yu(2,1) = pars(39);
   Yu(2,2) = pars(40);
   MuPr = pars(41);
   g1 = pars(42);
   g2 = pars(43);
   g3 = pars(44);
   gN = pars(45);
   vd = pars(46);
   vu = pars(47);
   vs = pars(48);

}

const CE6SSM_input_parameters& CE6SSM_susy_parameters::get_input() const
{
   return input;
}

void CE6SSM_susy_parameters::set_input_parameters(const CE6SSM_input_parameters& input_)
{
   input = input_;
}

void CE6SSM_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappa = (Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceLambda12AdjLambda12 = (Lambda12*(Lambda12).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = (Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12 = (Lambda12*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;

}

std::ostream& operator<<(std::ostream& ostr, const CE6SSM_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
