// ====================================================================
// Test suite for implementation of third generation sfermion 
// calculation in the case of a U(1) extended model with one
// SM singlet
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "CE6SSM_two_scale_model.hpp"

#include "linalg2.hpp"
#include "sfermions.hpp"
#include "wrappers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CE6SSM_3rd_gen_sfermions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace flexiblesusy;

CE6SSM_input_parameters get_default_input()
{
   CE6SSM_input_parameters input;

   input.m0 = 500.;
   input.m12 = 1000.;
   input.TanBeta = 10.;
   input.Azero = 500.;
   input.LambdaInput = 0.4;
   input.KappaInput = 0.2;
   input.muPrimeInput = 500.;
   input.BmuPrimeInput = 5000.;
   input.vSInput = 6700.;
   input.Lambda12Input = 0.1;

   return input;
}

void set_default_susy_parameters(CE6SSM<Two_scale>& model)
{
   Eigen::Matrix<double,3,3> YdInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> YeInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> KappaInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,2,2> Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,3,3> YuInput = Eigen::Matrix<double,3,3>::Zero();

   YdInput(2,2) = 0.1;
   YeInput(2,2) = 0.05;

   KappaInput(0,0) = 0.2;
   KappaInput(1,1) = 0.2;
   KappaInput(2,2) = 0.2;

   Lambda12Input(0,0) = 0.1;
   Lambda12Input(1,1) = 0.1;

   const double LambdaxInput = 0.4;

   YuInput(2,2) = 0.8;

   const double MuPrInput = 500.;
   const double g1Input = 0.46;
   const double g2Input = 0.64;
   const double g3Input = 1.03;
   const double gNInput = 0.47;
   const double vdInput = 24.5;
   const double vuInput = 244.8;
   const double vsInput = 6700.;

   model.set_Yd(YdInput);
   model.set_Ye(YeInput);
   model.set_Kappa(KappaInput);
   model.set_Lambda12(Lambda12Input);
   model.set_Lambdax(LambdaxInput);
   model.set_Yu(YuInput);
   model.set_MuPr(MuPrInput);
   model.set_g1(g1Input);
   model.set_g2(g2Input);
   model.set_g3(g3Input);
   model.set_gN(gNInput);
   model.set_vd(vdInput);
   model.set_vu(vuInput);
   model.set_vs(vsInput);
}

void set_default_soft_parameters(CE6SSM<Two_scale>& model)
{
   Eigen::Matrix<double,3,3> TYdInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> TYeInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> TKappaInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,2,2> TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,3,3> TYuInput = Eigen::Matrix<double,3,3>::Zero();

   TYdInput(2,2) = 100.;
   TYeInput(2,2) = 50.;

   TKappaInput(0,0) = 200.;
   TKappaInput(1,1) = 200.;
   TKappaInput(2,2) = 200.;

   TLambda12Input(0,0) = 50.;
   TLambda12Input(1,1) = 50.;

   const double TLambdaxInput = 100.;

   TYuInput(2,2) = 500.;

   const double BMuPrInput = 5000.;

   Eigen::Matrix<double,3,3> mq2Input = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> ml2Input = Eigen::Matrix<double,3,3>::Zero();

   mq2Input(0,0) = 2.5e7;
   mq2Input(1,1) = 1.0e6;
   mq2Input(2,2) = 5.0e5;

   ml2Input(0,0) = 3.0e7;
   ml2Input(1,1) = 4.0e5;
   ml2Input(2,2) = 1.0e5;

   const double mHd2Input = 1.0e6;
   const double mHu2Input = -1.0e5;

   Eigen::Matrix<double,3,3> md2Input = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> mu2Input = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> me2Input = Eigen::Matrix<double,3,3>::Zero();

   md2Input(0,0) = 2.5e7;
   md2Input(1,1) = 2.5e7;
   md2Input(2,2) = 9.0e4;

   mu2Input(0,0) = 4.0e7;
   mu2Input(1,1) = 8.0e6;
   mu2Input(2,2) = 7.0e5;

   me2Input(0,0) = 2.5e7;
   me2Input(1,1) = 2.5e7;
   me2Input(2,2) = 5.0e4;

   const double ms2Input = -1.0e6;

   Eigen::Matrix<double,2,2> mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,2,2> mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,2,2> msI2Input = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,3,3> mDx2Input = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();

   mH1I2Input(0,0) = 2.5e7;
   mH1I2Input(1,1) = 2.5e7;

   mH2I2Input(0,0) = 2.5e7;
   mH2I2Input(1,1) = 2.5e7;

   msI2Input(0,0) = 2.5e7;
   msI2Input(1,1) = 2.5e7;

   mDx2Input(0,0) = 2.5e7;
   mDx2Input(1,1) = 2.5e7;
   mDx2Input(2,2) = 2.5e7;

   mDxbar2Input(0,0) = 2.5e7;
   mDxbar2Input(1,1) = 2.5e7;
   mDxbar2Input(2,2) = 2.5e7;

   const double mHp2Input = 2.5e5;
   const double mHpbar2Input = 2.5e7;

   const double MassBInput = 100.;
   const double MassWBInput = 1100.;
   const double MassGInput = 2000.;
   const double MassBpInput = 150.;

   model.set_TYd(TYdInput);
   model.set_TYe(TYeInput);
   model.set_TKappa(TKappaInput);
   model.set_TLambda12(TLambda12Input);
   model.set_TLambdax(TLambdaxInput);
   model.set_TYu(TYuInput);
   model.set_BMuPr(BMuPrInput);
   model.set_mq2(mq2Input);
   model.set_ml2(ml2Input);
   model.set_mHd2(mHd2Input);
   model.set_mHu2(mHu2Input);
   model.set_md2(md2Input);
   model.set_mu2(mu2Input);
   model.set_me2(me2Input);
   model.set_ms2(ms2Input);
   model.set_mH1I2(mH1I2Input);
   model.set_mH2I2(mH2I2Input);
   model.set_msI2(msI2Input);
   model.set_mDx2(mDx2Input);
   model.set_mDxbar2(mDxbar2Input);
   model.set_mHp2(mHp2Input);
   model.set_mHpbar2(mHpbar2Input);
   model.set_MassB(MassBInput);
   model.set_MassWB(MassWBInput);
   model.set_MassG(MassGInput);
   model.set_MassBp(MassBpInput);
}

CE6SSM<Two_scale> get_model()
{
   CE6SSM<Two_scale> model;

   CE6SSM_input_parameters input = get_default_input();

   model.set_input_parameters(input);

   set_default_susy_parameters(model);
   set_default_soft_parameters(model);

   return model;
}
Eigen::Array<double,2,1> calculate_sfermion_3rd_generation(const CE6SSM<Two_scale>& model, CE6SSM_info::Particles p)
{
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> full_mass_matrix;
   switch(p) {
   case CE6SSM_info::Su:
      full_mass_matrix = model.get_mass_matrix_Su();
      break;
   case CE6SSM_info::Sd:
      full_mass_matrix = model.get_mass_matrix_Sd();
      break;
   case CE6SSM_info::Se:
      full_mass_matrix = model.get_mass_matrix_Se();
      break;
   case CE6SSM_info::Sv:
      full_mass_matrix = model.get_mass_matrix_Sv();
      break;
   default:
      full_mass_matrix = Eigen::Matrix<double,6,6>::Zero();
   }

   Eigen::Matrix<double,2,2> mass_matrix_3rd_gen;
   if (p != CE6SSM_info::Sv) {
      mass_matrix_3rd_gen(0,0) = full_mass_matrix(2,2);
      mass_matrix_3rd_gen(0,1) = full_mass_matrix(2,5);
      mass_matrix_3rd_gen(1,0) = full_mass_matrix(5,2);
      mass_matrix_3rd_gen(1,1) = full_mass_matrix(5,5);
   } else {
      mass_matrix_3rd_gen(0,0) = full_mass_matrix(2,2);
      mass_matrix_3rd_gen(0,1) = 0.;
      mass_matrix_3rd_gen(1,0) = 0.;
      mass_matrix_3rd_gen(1,1) = 0.;
   }
   Eigen::Array<double,2,1> MSf;
   Eigen::Matrix<double,2,2> Zf;

   fs_diagonalize_hermitian(mass_matrix_3rd_gen, MSf, Zf);
   
   return AbsSqrt(MSf);
}
// Test that the 3rd generation calculation is correct
// for stop masses
BOOST_AUTO_TEST_CASE( test_CE6SSM_stop_masses )
{
   const double max_err = 1.0e-6;

   CE6SSM<Two_scale> model = get_model();

   double mst1 = 0.;
   double mst2 = 0.;
   double thetat = 0.;

   model.calculate_MSu_3rd_generation(mst1, mst2, thetat);

   if (mst1 > mst2) {
      double temp = mst1;
      mst1 = mst2;
      mst2 = temp;
   }

   Eigen::Array<double,2,1> MSu = calculate_sfermion_3rd_generation(model, CE6SSM_info::Su);

   BOOST_CHECK_LE(Abs(MSu(0) - mst1), max_err);
   BOOST_CHECK_LE(Abs(MSu(1) - mst2), max_err);
}

// Test that the 3rd generation calculation is correct
// for sbottom masses
BOOST_AUTO_TEST_CASE( test_CE6SSM_sbottom_masses )
{
   const double max_err = 1.0e-6;

   CE6SSM<Two_scale> model = get_model();

   double msb1 = 0.;
   double msb2 = 0.;
   double thetab = 0.;

   model.calculate_MSd_3rd_generation(msb1, msb2, thetab);

   if (msb1 > msb2) {
      double temp = msb1;
      msb1 = msb2;
      msb2 = temp;
   }

   Eigen::Array<double,2,1> MSd = calculate_sfermion_3rd_generation(model, CE6SSM_info::Sd);

   BOOST_CHECK_LE(Abs(MSd(0) - msb1), max_err);
   BOOST_CHECK_LE(Abs(MSd(1) - msb2), max_err);
}


// Test that the 3rd generation calculation is correct
// for stau masses
BOOST_AUTO_TEST_CASE( test_CE6SSM_stau_masses )
{
   const double max_err = 1.0e-6;

   CE6SSM<Two_scale> model = get_model();

   double mstau1 = 0.;
   double mstau2 = 0.;
   double thetatau = 0.;

   model.calculate_MSe_3rd_generation(mstau1, mstau2, thetatau);

   if (mstau1 > mstau2) {
      double temp = mstau1;
      mstau1 = mstau2;
      mstau2 = temp;
   }

   Eigen::Array<double,2,1> MSe = calculate_sfermion_3rd_generation(model, CE6SSM_info::Se);

   BOOST_CHECK_LE(Abs(MSe(0) - mstau1), max_err);
   BOOST_CHECK_LE(Abs(MSe(1) - mstau2), max_err);
}


// Test that the 3rd generation calculation is correct
// for tau sneutrino masses
BOOST_AUTO_TEST_CASE( test_CE6SSM_tau_sneutrino_masses )
{
   const double max_err = 1.0e-6;

   CE6SSM<Two_scale> model = get_model();

   double msv1 = 0.;
   double msv2 = 0.;
   double thetav = 0.;

   model.calculate_MSv_3rd_generation(msv1, msv2, thetav);

   if (msv1 > msv2) {
      double temp = msv1;
      msv1 = msv2;
      msv2 = temp;
   }

   Eigen::Array<double,2,1> MSv = calculate_sfermion_3rd_generation(model, CE6SSM_info::Sv);

   BOOST_CHECK_LE(Abs(MSv(0) - msv1), max_err);
   BOOST_CHECK_LE(Abs(MSv(1) - msv2), max_err);
}
