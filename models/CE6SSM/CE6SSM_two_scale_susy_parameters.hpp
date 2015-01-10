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

#ifndef CE6SSM_TWO_SCALE_susy_parameters_H
#define CE6SSM_TWO_SCALE_susy_parameters_H

#include "betafunction.hpp"
#include "CE6SSM_input_parameters.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Susy_traces

class CE6SSM_susy_parameters : public Beta_function {
public:
   explicit CE6SSM_susy_parameters(const CE6SSM_input_parameters& input_ = CE6SSM_input_parameters());
   CE6SSM_susy_parameters(double scale_, double loops_, double thresholds_, const CE6SSM_input_parameters& input_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , const Eigen::Matrix<double,3,3>& Kappa_, const Eigen::Matrix<double,2,2>&
   Lambda12_, double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, double
   MuPr_, double g1_, double g2_, double g3_, double gN_, double vd_, double
   vu_, double vs_
);
   virtual ~CE6SSM_susy_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual const Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);
   const CE6SSM_input_parameters& get_input() const;
   void set_input_parameters(const CE6SSM_input_parameters&);

   CE6SSM_susy_parameters calc_beta() const;
   virtual void clear();

   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, double value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, double value) { Ye(i,k) = value; }
   void set_Kappa(const Eigen::Matrix<double,3,3>& Kappa_) { Kappa = Kappa_; }
   void set_Kappa(int i, int k, double value) { Kappa(i,k) = value; }
   void set_Lambda12(const Eigen::Matrix<double,2,2>& Lambda12_) { Lambda12 = Lambda12_; }
   void set_Lambda12(int i, int k, double value) { Lambda12(i,k) = value; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, double value) { Yu(i,k) = value; }
   void set_MuPr(double MuPr_) { MuPr = MuPr_; }
   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_gN(double gN_) { gN = gN_; }
   void set_vd(double vd_) { vd = vd_; }
   void set_vu(double vu_) { vu = vu_; }
   void set_vs(double vs_) { vs = vs_; }

   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   const Eigen::Matrix<double,3,3>& get_Kappa() const { return Kappa; }
   double get_Kappa(int i, int k) const { return Kappa(i,k); }
   const Eigen::Matrix<double,2,2>& get_Lambda12() const { return Lambda12; }
   double get_Lambda12(int i, int k) const { return Lambda12(i,k); }
   double get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   double get_MuPr() const { return MuPr; }
   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_gN() const { return gN; }
   double get_vd() const { return vd; }
   double get_vu() const { return vu; }
   double get_vs() const { return vs; }

   Eigen::Matrix<double,3,3> get_SqSq() const;
   Eigen::Matrix<double,3,3> get_SlSl() const;
   double get_SHdSHd() const;
   double get_SHuSHu() const;
   Eigen::Matrix<double,3,3> get_SdRSdR() const;
   Eigen::Matrix<double,3,3> get_SuRSuR() const;
   Eigen::Matrix<double,3,3> get_SeRSeR() const;
   double get_SsRSsR() const;
   Eigen::Matrix<double,2,2> get_SH1ISH1I() const;
   Eigen::Matrix<double,2,2> get_SH2ISH2I() const;
   Eigen::Matrix<double,2,2> get_SsIRSsIR() const;
   Eigen::Matrix<double,3,3> get_SDxLSDxL() const;
   Eigen::Matrix<double,3,3> get_SDxbarRSDxbarR() const;
   double get_SHpSHp() const;
   double get_SHpbarSHpbar() const;


protected:
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   Eigen::Matrix<double,3,3> Kappa;
   Eigen::Matrix<double,2,2> Lambda12;
   double Lambdax;
   Eigen::Matrix<double,3,3> Yu;
   double MuPr;
   double g1;
   double g2;
   double g3;
   double gN;
   double vd;
   double vu;
   double vs;

   CE6SSM_input_parameters input;

private:
   static const int numberOfParameters = 49;

   struct Susy_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceKappaAdjKappa;
      double traceLambda12AdjLambda12;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceKappaAdjKappaKappaAdjKappa;
      double traceLambda12AdjLambda12Lambda12AdjLambda12;
      double traceYuAdjYuYuAdjYu;

   };
   void calc_susy_traces(Susy_traces&) const;

   Eigen::Matrix<double,3,3> calc_beta_Yd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Kappa_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Kappa_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_Lambda12_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,2,2> calc_beta_Lambda12_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MuPr_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MuPr_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gN_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gN_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vd_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vd_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vu_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vs_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vs_two_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const CE6SSM_susy_parameters&);

} // namespace flexiblesusy

#endif
