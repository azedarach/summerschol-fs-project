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

// File generated at Fri 9 Jan 2015 15:06:51

#include "CNMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

CNMSSM_soft_parameters::CNMSSM_soft_parameters(const CNMSSM_input_parameters& input_)
   : CNMSSM_susy_parameters(input_)
   , TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TLambdax(0), TKappa(0), TYu(Eigen::Matrix<double,3,3>::Zero()),
   mq2(Eigen::Matrix<double,3,3>::Zero()), ml2(Eigen::Matrix<double,3,3>::Zero(
   )), mHd2(0), mHu2(0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(
   Eigen::Matrix<double,3,3>::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()),
   ms2(0), MassB(0), MassWB(0), MassG(0)

{
   set_number_of_parameters(numberOfParameters);
}

CNMSSM_soft_parameters::CNMSSM_soft_parameters(
   const CNMSSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, double TLambdax_, double TKappa_, const Eigen::Matrix<double,3,3>&
   TYu_, const Eigen::Matrix<double,3,3>& mq2_, const Eigen::Matrix<double,3,3>
   & ml2_, double mHd2_, double mHu2_, const Eigen::Matrix<double,3,3>& md2_,
   const Eigen::Matrix<double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_
   , double ms2_, double MassB_, double MassWB_, double MassG_

)
   : CNMSSM_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TLambdax(TLambdax_), TKappa(TKappa_), TYu(TYu_), mq2
   (mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_)
   , ms2(ms2_), MassB(MassB_), MassWB(MassWB_), MassG(MassG_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd CNMSSM_soft_parameters::beta() const
{
   return calc_beta().get();
}

CNMSSM_soft_parameters CNMSSM_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   Eigen::Matrix<double,3,3> beta_TYd(calc_beta_TYd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYe(calc_beta_TYe_one_loop(TRACE_STRUCT));
   double beta_TLambdax(calc_beta_TLambdax_one_loop(TRACE_STRUCT));
   double beta_TKappa(calc_beta_TKappa_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYu(calc_beta_TYu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mq2(calc_beta_mq2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_ml2(calc_beta_ml2_one_loop(TRACE_STRUCT));
   double beta_mHd2(calc_beta_mHd2_one_loop(TRACE_STRUCT));
   double beta_mHu2(calc_beta_mHu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_md2(calc_beta_md2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_me2(calc_beta_me2_one_loop(TRACE_STRUCT));
   double beta_ms2(calc_beta_ms2_one_loop(TRACE_STRUCT));
   double beta_MassB(calc_beta_MassB_one_loop(TRACE_STRUCT));
   double beta_MassWB(calc_beta_MassWB_one_loop(TRACE_STRUCT));
   double beta_MassG(calc_beta_MassG_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_TYd += calc_beta_TYd_two_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_two_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_two_loop(TRACE_STRUCT);
      beta_TKappa += calc_beta_TKappa_two_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_two_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_two_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_two_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_two_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_two_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_two_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_two_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_two_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);

   }


   const CNMSSM_susy_parameters susy_betas(CNMSSM_susy_parameters::calc_beta());

   return CNMSSM_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TLambdax, beta_TKappa, beta_TYu, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_MassB, beta_MassWB, beta_MassG);
}

void CNMSSM_soft_parameters::clear()
{
   CNMSSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TLambdax = 0.;
   TKappa = 0.;
   TYu = Eigen::Matrix<double,3,3>::Zero();
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;

}

const Eigen::ArrayXd CNMSSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(CNMSSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(35) = TYd(0,0);
   pars(36) = TYd(0,1);
   pars(37) = TYd(0,2);
   pars(38) = TYd(1,0);
   pars(39) = TYd(1,1);
   pars(40) = TYd(1,2);
   pars(41) = TYd(2,0);
   pars(42) = TYd(2,1);
   pars(43) = TYd(2,2);
   pars(44) = TYe(0,0);
   pars(45) = TYe(0,1);
   pars(46) = TYe(0,2);
   pars(47) = TYe(1,0);
   pars(48) = TYe(1,1);
   pars(49) = TYe(1,2);
   pars(50) = TYe(2,0);
   pars(51) = TYe(2,1);
   pars(52) = TYe(2,2);
   pars(53) = TLambdax;
   pars(54) = TKappa;
   pars(55) = TYu(0,0);
   pars(56) = TYu(0,1);
   pars(57) = TYu(0,2);
   pars(58) = TYu(1,0);
   pars(59) = TYu(1,1);
   pars(60) = TYu(1,2);
   pars(61) = TYu(2,0);
   pars(62) = TYu(2,1);
   pars(63) = TYu(2,2);
   pars(64) = mq2(0,0);
   pars(65) = mq2(0,1);
   pars(66) = mq2(0,2);
   pars(67) = mq2(1,0);
   pars(68) = mq2(1,1);
   pars(69) = mq2(1,2);
   pars(70) = mq2(2,0);
   pars(71) = mq2(2,1);
   pars(72) = mq2(2,2);
   pars(73) = ml2(0,0);
   pars(74) = ml2(0,1);
   pars(75) = ml2(0,2);
   pars(76) = ml2(1,0);
   pars(77) = ml2(1,1);
   pars(78) = ml2(1,2);
   pars(79) = ml2(2,0);
   pars(80) = ml2(2,1);
   pars(81) = ml2(2,2);
   pars(82) = mHd2;
   pars(83) = mHu2;
   pars(84) = md2(0,0);
   pars(85) = md2(0,1);
   pars(86) = md2(0,2);
   pars(87) = md2(1,0);
   pars(88) = md2(1,1);
   pars(89) = md2(1,2);
   pars(90) = md2(2,0);
   pars(91) = md2(2,1);
   pars(92) = md2(2,2);
   pars(93) = mu2(0,0);
   pars(94) = mu2(0,1);
   pars(95) = mu2(0,2);
   pars(96) = mu2(1,0);
   pars(97) = mu2(1,1);
   pars(98) = mu2(1,2);
   pars(99) = mu2(2,0);
   pars(100) = mu2(2,1);
   pars(101) = mu2(2,2);
   pars(102) = me2(0,0);
   pars(103) = me2(0,1);
   pars(104) = me2(0,2);
   pars(105) = me2(1,0);
   pars(106) = me2(1,1);
   pars(107) = me2(1,2);
   pars(108) = me2(2,0);
   pars(109) = me2(2,1);
   pars(110) = me2(2,2);
   pars(111) = ms2;
   pars(112) = MassB;
   pars(113) = MassWB;
   pars(114) = MassG;


   return pars;
}

void CNMSSM_soft_parameters::print(std::ostream& ostr) const
{
   CNMSSM_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "TKappa = " << TKappa << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';

}

void CNMSSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   CNMSSM_susy_parameters::set(pars);

   TYd(0,0) = pars(35);
   TYd(0,1) = pars(36);
   TYd(0,2) = pars(37);
   TYd(1,0) = pars(38);
   TYd(1,1) = pars(39);
   TYd(1,2) = pars(40);
   TYd(2,0) = pars(41);
   TYd(2,1) = pars(42);
   TYd(2,2) = pars(43);
   TYe(0,0) = pars(44);
   TYe(0,1) = pars(45);
   TYe(0,2) = pars(46);
   TYe(1,0) = pars(47);
   TYe(1,1) = pars(48);
   TYe(1,2) = pars(49);
   TYe(2,0) = pars(50);
   TYe(2,1) = pars(51);
   TYe(2,2) = pars(52);
   TLambdax = pars(53);
   TKappa = pars(54);
   TYu(0,0) = pars(55);
   TYu(0,1) = pars(56);
   TYu(0,2) = pars(57);
   TYu(1,0) = pars(58);
   TYu(1,1) = pars(59);
   TYu(1,2) = pars(60);
   TYu(2,0) = pars(61);
   TYu(2,1) = pars(62);
   TYu(2,2) = pars(63);
   mq2(0,0) = pars(64);
   mq2(0,1) = pars(65);
   mq2(0,2) = pars(66);
   mq2(1,0) = pars(67);
   mq2(1,1) = pars(68);
   mq2(1,2) = pars(69);
   mq2(2,0) = pars(70);
   mq2(2,1) = pars(71);
   mq2(2,2) = pars(72);
   ml2(0,0) = pars(73);
   ml2(0,1) = pars(74);
   ml2(0,2) = pars(75);
   ml2(1,0) = pars(76);
   ml2(1,1) = pars(77);
   ml2(1,2) = pars(78);
   ml2(2,0) = pars(79);
   ml2(2,1) = pars(80);
   ml2(2,2) = pars(81);
   mHd2 = pars(82);
   mHu2 = pars(83);
   md2(0,0) = pars(84);
   md2(0,1) = pars(85);
   md2(0,2) = pars(86);
   md2(1,0) = pars(87);
   md2(1,1) = pars(88);
   md2(1,2) = pars(89);
   md2(2,0) = pars(90);
   md2(2,1) = pars(91);
   md2(2,2) = pars(92);
   mu2(0,0) = pars(93);
   mu2(0,1) = pars(94);
   mu2(0,2) = pars(95);
   mu2(1,0) = pars(96);
   mu2(1,1) = pars(97);
   mu2(1,2) = pars(98);
   mu2(2,0) = pars(99);
   mu2(2,1) = pars(100);
   mu2(2,2) = pars(101);
   me2(0,0) = pars(102);
   me2(0,1) = pars(103);
   me2(0,2) = pars(104);
   me2(1,0) = pars(105);
   me2(1,1) = pars(106);
   me2(1,2) = pars(107);
   me2(2,0) = pars(108);
   me2(2,1) = pars(109);
   me2(2,2) = pars(110);
   ms2 = pars(111);
   MassB = pars(112);
   MassWB = pars(113);
   MassG = pars(114);

}

void CNMSSM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   TRACE_STRUCT.traceAdjYdTYd = (Yd.adjoint()*TYd).trace();
   TRACE_STRUCT.traceAdjYeTYe = (Ye.adjoint()*TYe).trace();
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceAdjYuTYu = (Yu.adjoint()*TYu).trace();
   TRACE_STRUCT.traceYdAdjYdTYdAdjYd = (Yd*Yd.adjoint()*TYd*Yd.adjoint()).trace
      ();
   TRACE_STRUCT.traceYdAdjYuTYuAdjYd = (Yd*Yu.adjoint()*TYu*Yd.adjoint()).trace
      ();
   TRACE_STRUCT.traceYeAdjYeTYeAdjYe = (Ye*Ye.adjoint()*TYe*Ye.adjoint()).trace
      ();
   TRACE_STRUCT.traceYuAdjYdTYdAdjYu = (Yu*Yd.adjoint()*TYd*Yu.adjoint()).trace
      ();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYuAdjYuTYuAdjYu = (Yu*Yu.adjoint()*TYu*Yu.adjoint()).trace
      ();
   TRACE_STRUCT.traceconjTYdTpTYd = (TYd.conjugate()*(TYd).transpose()).trace()
      ;
   TRACE_STRUCT.traceconjTYeTpTYe = (TYe.conjugate()*(TYe).transpose()).trace()
      ;
   TRACE_STRUCT.tracemd2YdAdjYd = (md2*Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceme2YeAdjYe = (me2*Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceml2AdjYeYe = (ml2*Ye.adjoint()*Ye).trace();
   TRACE_STRUCT.tracemq2AdjYdYd = (mq2*Yd.adjoint()*Yd).trace();
   TRACE_STRUCT.traceconjTYdTpYd = (TYd.conjugate()*Yd.transpose()).trace();
   TRACE_STRUCT.traceconjTYeTpYe = (TYe.conjugate()*Ye.transpose()).trace();
   TRACE_STRUCT.traceconjTYuTpTYu = (TYu.conjugate()*(TYu).transpose()).trace()
      ;
   TRACE_STRUCT.tracemq2AdjYuYu = (mq2*Yu.adjoint()*Yu).trace();
   TRACE_STRUCT.tracemu2YuAdjYu = (mu2*Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceconjTYuTpYu = (TYu.conjugate()*Yu.transpose()).trace();
   TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = (Yd*Yd.adjoint()*TYd*(TYd).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = (Yd*Yu.adjoint()*TYu*(TYd).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = (Yd*(TYd).adjoint()*TYd*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = (Yd*(TYu).adjoint()*TYu*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = (Ye*Ye.adjoint()*TYe*(TYe).adjoint())
      .trace();
   TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = (Ye*(TYe).adjoint()*TYe*Ye.adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = (Yu*Yd.adjoint()*TYd*(TYu).adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = (Yu*(TYd).adjoint()*TYd*Yu.adjoint())
      .trace();
   TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = (md2*Yd*Yd.adjoint()*Yd*Yd.adjoint())
      .trace();
   TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = (md2*Yd*Yu.adjoint()*Yu*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = (me2*Ye*Ye.adjoint()*Ye*Ye.adjoint())
      .trace();
   TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = (ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye)
      .trace();
   TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = (mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      .trace();
   TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = (mq2*Yd.adjoint()*Yd*Yu.adjoint()*Yu)
      .trace();
   TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = (mq2*Yu.adjoint()*Yu*Yd.adjoint()*Yd)
      .trace();
   TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = (mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = (Yu*Yu.adjoint()*TYu*(TYu).adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = (Yu*(TYu).adjoint()*TYu*Yu.adjoint())
      .trace();
   TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = (mq2*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      .trace();
   TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = (mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint())
      .trace();


   TRACE_STRUCT.Tr11 = 0.7745966692414834*g1*(-mHd2 + mHu2 + (md2).trace() + (
      me2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace());
   TRACE_STRUCT.Tr2U111 = 0.1*Sqr(g1)*(3*mHd2 + 3*mHu2 + 2*(md2).trace() + 6*(
      me2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*(mu2).trace());
   TRACE_STRUCT.Tr31 = 0.012909944487358056*g1*(30*(mHd2 - mHu2)*AbsSqr(Lambdax
      ) - 9*mHd2*Sqr(g1) + 9*mHu2*Sqr(g1) - 45*mHd2*Sqr(g2) + 45*mHu2*Sqr(g2) + 4*
      (Sqr(g1) + 20*Sqr(g3))*(md2).trace() + 36*Sqr(g1)*(me2).trace() - 9*Sqr(g1)*
      (ml2).trace() - 45*Sqr(g2)*(ml2).trace() + Sqr(g1)*(mq2).trace() + 45*Sqr(g2
      )*(mq2).trace() + 80*Sqr(g3)*(mq2).trace() - 32*Sqr(g1)*(mu2).trace() - 160*
      Sqr(g3)*(mu2).trace() + 90*mHd2*(Yd*Yd.adjoint()).trace() + 30*mHd2*(Ye*
      Ye.adjoint()).trace() - 90*mHu2*(Yu*Yu.adjoint()).trace() - 60*(Yd*
      Yd.adjoint()*md2.conjugate()).trace() - 30*(Yd*mq2.conjugate()*Yd.adjoint())
      .trace() - 60*(Ye*Ye.adjoint()*me2.conjugate()).trace() + 30*(Ye*
      ml2.conjugate()*Ye.adjoint()).trace() + 120*(Yu*Yu.adjoint()*mu2.conjugate()
      ).trace() - 30*(Yu*mq2.conjugate()*Yu.adjoint()).trace());
   TRACE_STRUCT.Tr22 = 0.5*(mHd2 + mHu2 + (ml2).trace() + 3*(mq2).trace());
   TRACE_STRUCT.Tr23 = 0.5*((md2).trace() + 2*(mq2).trace() + (mu2).trace());

}

std::ostream& operator<<(std::ostream& ostr, const CNMSSM_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
