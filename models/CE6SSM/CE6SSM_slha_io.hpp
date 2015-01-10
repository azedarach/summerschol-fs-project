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

#ifndef CE6SSM_SLHA_IO_H
#define CE6SSM_SLHA_IO_H

#include "CE6SSM_two_scale_model.hpp"
#include "CE6SSM_info.hpp"
#include "CE6SSM_physical.hpp"
#include "slha_io.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_POLE_MASS(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(PHYSICAL(p))>::type>::type p;
#define SM(p) Electroweak_constants::p

namespace flexiblesusy {

struct CE6SSM_input_parameters;
class Spectrum_generator_settings;

class CE6SSM_slha_io {
public:
   CE6SSM_slha_io();
   ~CE6SSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(CE6SSM_input_parameters&) const;
   template <class T> void fill(CE6SSM<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_input_scale() const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const CE6SSM_input_parameters&);
   template <class T> void set_extra(const CE6SSM<T>&);
   void set_minpar(const CE6SSM_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const CE6SSM<T>&);
   void set_spinfo(const Problems<CE6SSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(CE6SSM_input_parameters&, int, double);
   static void fill_extpar_tuple(CE6SSM_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 24;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   static void convert_to_slha_convention(CE6SSM_physical&);
   void set_mass(const CE6SSM_physical&, bool);
   void set_mixing_matrices(const CE6SSM_physical&, bool);
   template <class T> void set_model_parameters(const CE6SSM<T>&);
   double read_scale() const;
   template <class T> void fill_drbar_parameters(CE6SSM<T>&) const;
   template <class T> void fill_physical(CE6SSM<T>&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void CE6SSM_slha_io::fill(CE6SSM<T>& model) const
{
   fill_drbar_parameters(model);
   fill_physical(model);
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
template <class T>
void CE6SSM_slha_io::fill_drbar_parameters(CE6SSM<T>& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   model.set_gN(slha_io.read_entry("gauge", 4));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   {
      DEFINE_PARAMETER(TYe);
      slha_io.read_block("Te", TYe);
      model.set_TYe(TYe);
   }
   {
      DEFINE_PARAMETER(TYd);
      slha_io.read_block("Td", TYd);
      model.set_TYd(TYd);
   }
   {
      DEFINE_PARAMETER(TYu);
      slha_io.read_block("Tu", TYu);
      model.set_TYu(TYu);
   }
   {
      DEFINE_PARAMETER(mq2);
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      DEFINE_PARAMETER(me2);
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      DEFINE_PARAMETER(ml2);
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      DEFINE_PARAMETER(mu2);
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      DEFINE_PARAMETER(md2);
      slha_io.read_block("MSD2", md2);
      model.set_md2(md2);
   }
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   {
      DEFINE_PARAMETER(mH1I2);
      slha_io.read_block("mHdInert2", mH1I2);
      model.set_mH1I2(mH1I2);
   }
   {
      DEFINE_PARAMETER(mH2I2);
      slha_io.read_block("mHuInert2", mH2I2);
      model.set_mH2I2(mH2I2);
   }
   {
      DEFINE_PARAMETER(mDx2);
      slha_io.read_block("mX2", mDx2);
      model.set_mDx2(mDx2);
   }
   {
      DEFINE_PARAMETER(mDxbar2);
      slha_io.read_block("mXBar2", mDxbar2);
      model.set_mDxbar2(mDxbar2);
   }
   model.set_ms2(slha_io.read_entry("MSOFT", 23));
   {
      DEFINE_PARAMETER(msI2);
      slha_io.read_block("msInert2", msI2);
      model.set_msI2(msI2);
   }
   model.set_mHp2(slha_io.read_entry("MSOFT", 24));
   model.set_mHpbar2(slha_io.read_entry("MSOFT", 25));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_MassBp(slha_io.read_entry("MSOFT", 4));
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));
   model.set_vs(slha_io.read_entry("ESIXRUN", 11));
   {
      DEFINE_PARAMETER(Kappa);
      slha_io.read_block("ESIXKAPPA", Kappa);
      model.set_Kappa(Kappa);
   }
   {
      DEFINE_PARAMETER(TKappa);
      slha_io.read_block("ESIXTKAPPA", TKappa);
      model.set_TKappa(TKappa);
   }
   model.set_Lambdax(slha_io.read_entry("ESIXRUN", 1));
   model.set_TLambdax(slha_io.read_entry("ESIXRUN", 2));
   {
      DEFINE_PARAMETER(Lambda12);
      slha_io.read_block("ESIXLAMBDA", Lambda12);
      model.set_Lambda12(Lambda12);
   }
   {
      DEFINE_PARAMETER(TLambda12);
      slha_io.read_block("ESIXTLAMBDA", TLambda12);
      model.set_TLambda12(TLambda12);
   }
   model.set_MuPr(slha_io.read_entry("ESIXRUN", 0));
   model.set_BMuPr(slha_io.read_entry("ESIXRUN", 101));


   model.set_scale(read_scale());
}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
template <class T>
void CE6SSM_slha_io::fill_physical(CE6SSM<T>& model) const
{
   {
      DEFINE_PARAMETER(ZD);
      slha_io.read_block("DSQMIX", ZD);
      PHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PARAMETER(ZU);
      slha_io.read_block("USQMIX", ZU);
      PHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PARAMETER(ZE);
      slha_io.read_block("SELMIX", ZE);
      PHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PARAMETER(ZV);
      slha_io.read_block("SNUMIX", ZV);
      PHYSICAL(ZV) = ZV;
   }
   {
      DEFINE_PARAMETER(ZDX);
      slha_io.read_block("ESIXZDX", ZDX);
      PHYSICAL(ZDX) = ZDX;
   }
   {
      DEFINE_PARAMETER(ZH);
      slha_io.read_block("NMHMIX", ZH);
      PHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PARAMETER(ZA);
      slha_io.read_block("NMAMIX", ZA);
      PHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      PHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PARAMETER(ZN);
      slha_io.read_block("NMNMIX", ZN);
      PHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PARAMETER(ZNI);
      slha_io.read_block("ESIXZNI", ZNI);
      PHYSICAL(ZNI) = ZNI;
   }
   {
      DEFINE_PARAMETER(ZNp);
      slha_io.read_block("ZNPMIX", ZNp);
      PHYSICAL(ZNp) = ZNp;
   }
   {
      DEFINE_PARAMETER(ZMI);
      slha_io.read_block("ESIXZMI", ZMI);
      PHYSICAL(ZMI) = ZMI;
   }
   {
      DEFINE_PARAMETER(ZPI);
      slha_io.read_block("ESIXZPI", ZPI);
      PHYSICAL(ZPI) = ZPI;
   }
   {
      DEFINE_PARAMETER(ZFSI);
      slha_io.read_block("ESIXZSI", ZFSI);
      PHYSICAL(ZFSI) = ZFSI;
   }
   {
      DEFINE_PARAMETER(ZSSI);
      slha_io.read_block("ZSSI", ZSSI);
      PHYSICAL(ZSSI) = ZSSI;
   }
   {
      DEFINE_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      PHYSICAL(UP) = UP;
   }
   {
      DEFINE_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      PHYSICAL(UM) = UM;
   }
   {
      DEFINE_PARAMETER(UHI0);
      slha_io.read_block("INHMIX", UHI0);
      PHYSICAL(UHI0) = UHI0;
   }
   {
      DEFINE_PARAMETER(UHIp);
      slha_io.read_block("ICHMIX", UHIp);
      PHYSICAL(UHIp) = UHIp;
   }
   {
      DEFINE_PARAMETER(ZEL);
      slha_io.read_block("UELMIX", ZEL);
      PHYSICAL(ZEL) = ZEL;
   }
   {
      DEFINE_PARAMETER(ZER);
      slha_io.read_block("UERMIX", ZER);
      PHYSICAL(ZER) = ZER;
   }
   {
      DEFINE_PARAMETER(ZDL);
      slha_io.read_block("UDLMIX", ZDL);
      PHYSICAL(ZDL) = ZDL;
   }
   {
      DEFINE_PARAMETER(ZDR);
      slha_io.read_block("UDRMIX", ZDR);
      PHYSICAL(ZDR) = ZDR;
   }
   {
      DEFINE_PARAMETER(ZUL);
      slha_io.read_block("UULMIX", ZUL);
      PHYSICAL(ZUL) = ZUL;
   }
   {
      DEFINE_PARAMETER(ZUR);
      slha_io.read_block("UURMIX", ZUR);
      PHYSICAL(ZUR) = ZUR;
   }
   {
      DEFINE_PARAMETER(ZDXL);
      slha_io.read_block("ESIXZXL", ZDXL);
      PHYSICAL(ZDXL) = ZDXL;
   }
   {
      DEFINE_PARAMETER(ZDXR);
      slha_io.read_block("ESIXZXR", ZDXR);
      PHYSICAL(ZDXR) = ZDXR;
   }
   {
      DEFINE_PARAMETER(UHp0);
      slha_io.read_block("UHNPMIX", UHp0);
      PHYSICAL(UHp0) = UHp0;
   }
   {
      DEFINE_PARAMETER(UHpp);
      slha_io.read_block("UHPPMIX", UHpp);
      PHYSICAL(UHpp) = UHpp;
   }

   PHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   PHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   PHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   PHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   PHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   PHYSICAL(MChaP) = slha_io.read_entry("MASS", 1000091);
   PHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   PHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   PHYSICAL(MVZp) = slha_io.read_entry("MASS", 31);
   PHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   PHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 1000003);
   PHYSICAL(MSd)(2) = slha_io.read_entry("MASS", 1000005);
   PHYSICAL(MSd)(3) = slha_io.read_entry("MASS", 2000001);
   PHYSICAL(MSd)(4) = slha_io.read_entry("MASS", 2000003);
   PHYSICAL(MSd)(5) = slha_io.read_entry("MASS", 2000005);
   PHYSICAL(MSv)(0) = slha_io.read_entry("MASS", 1000012);
   PHYSICAL(MSv)(1) = slha_io.read_entry("MASS", 1000014);
   PHYSICAL(MSv)(2) = slha_io.read_entry("MASS", 1000016);
   PHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   PHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 1000004);
   PHYSICAL(MSu)(2) = slha_io.read_entry("MASS", 1000006);
   PHYSICAL(MSu)(3) = slha_io.read_entry("MASS", 2000002);
   PHYSICAL(MSu)(4) = slha_io.read_entry("MASS", 2000004);
   PHYSICAL(MSu)(5) = slha_io.read_entry("MASS", 2000006);
   PHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   PHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 1000013);
   PHYSICAL(MSe)(2) = slha_io.read_entry("MASS", 1000015);
   PHYSICAL(MSe)(3) = slha_io.read_entry("MASS", 2000011);
   PHYSICAL(MSe)(4) = slha_io.read_entry("MASS", 2000013);
   PHYSICAL(MSe)(5) = slha_io.read_entry("MASS", 2000015);
   PHYSICAL(MSDX)(0) = slha_io.read_entry("MASS", 1000051);
   PHYSICAL(MSDX)(1) = slha_io.read_entry("MASS", 2000051);
   PHYSICAL(MSDX)(2) = slha_io.read_entry("MASS", 1000052);
   PHYSICAL(MSDX)(3) = slha_io.read_entry("MASS", 2000052);
   PHYSICAL(MSDX)(4) = slha_io.read_entry("MASS", 1000053);
   PHYSICAL(MSDX)(5) = slha_io.read_entry("MASS", 2000053);
   PHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   PHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   PHYSICAL(Mhh)(2) = slha_io.read_entry("MASS", 45);
   PHYSICAL(MAh)(2) = slha_io.read_entry("MASS", 36);
   PHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   PHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   PHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   PHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   PHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   PHYSICAL(MChi)(4) = slha_io.read_entry("MASS", 1000045);
   PHYSICAL(MChi)(5) = slha_io.read_entry("MASS", 1000055);
   PHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   PHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   PHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   PHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   PHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   PHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   PHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   PHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   PHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   PHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   PHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   PHYSICAL(MFDX)(0) = slha_io.read_entry("MASS", 51);
   PHYSICAL(MFDX)(1) = slha_io.read_entry("MASS", 52);
   PHYSICAL(MFDX)(2) = slha_io.read_entry("MASS", 53);
   PHYSICAL(MSHI0)(0) = slha_io.read_entry("MASS", 82);
   PHYSICAL(MSHI0)(1) = slha_io.read_entry("MASS", 86);
   PHYSICAL(MSHI0)(2) = slha_io.read_entry("MASS", 84);
   PHYSICAL(MSHI0)(3) = slha_io.read_entry("MASS", 88);
   PHYSICAL(MSHIp)(0) = slha_io.read_entry("MASS", 81);
   PHYSICAL(MSHIp)(1) = slha_io.read_entry("MASS", 85);
   PHYSICAL(MSHIp)(2) = slha_io.read_entry("MASS", 83);
   PHYSICAL(MSHIp)(3) = slha_io.read_entry("MASS", 87);
   PHYSICAL(MChaI)(0) = slha_io.read_entry("MASS", 1000085);
   PHYSICAL(MChaI)(1) = slha_io.read_entry("MASS", 1000086);
   PHYSICAL(MChiI)(0) = slha_io.read_entry("MASS", 1000081);
   PHYSICAL(MChiI)(1) = slha_io.read_entry("MASS", 1000082);
   PHYSICAL(MChiI)(2) = slha_io.read_entry("MASS", 1000083);
   PHYSICAL(MChiI)(3) = slha_io.read_entry("MASS", 1000084);
   PHYSICAL(MSSI0)(0) = slha_io.read_entry("MASS", 89);
   PHYSICAL(MSSI0)(1) = slha_io.read_entry("MASS", 90);
   PHYSICAL(MFSI)(0) = slha_io.read_entry("MASS", 1000089);
   PHYSICAL(MFSI)(1) = slha_io.read_entry("MASS", 1000090);
   PHYSICAL(MSHp0)(0) = slha_io.read_entry("MASS", 92);
   PHYSICAL(MSHp0)(1) = slha_io.read_entry("MASS", 94);
   PHYSICAL(MSHpp)(0) = slha_io.read_entry("MASS", 91);
   PHYSICAL(MSHpp)(1) = slha_io.read_entry("MASS", 93);
   PHYSICAL(MChiP)(0) = slha_io.read_entry("MASS", 1000092);
   PHYSICAL(MChiP)(1) = slha_io.read_entry("MASS", 1000094);
   PHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CE6SSM_slha_io::set_model_parameters(const CE6SSM<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(gN)), "gN")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", MODELPARAMETER(Yu), "Yu", model.get_scale());
   slha_io.set_block("Yd", MODELPARAMETER(Yd), "Yd", model.get_scale());
   slha_io.set_block("Ye", MODELPARAMETER(Ye), "Ye", model.get_scale());
   slha_io.set_block("Te", MODELPARAMETER(TYe), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu), "TYu", model.get_scale());
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (MODELPARAMETER(mHd2)), "mHd2")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(mHu2)), "mHu2")
            << FORMAT_ELEMENT(23, (MODELPARAMETER(ms2)), "ms2")
            << FORMAT_ELEMENT(24, (MODELPARAMETER(mHp2)), "mHp2")
            << FORMAT_ELEMENT(25, (MODELPARAMETER(mHpbar2)), "mHpbar2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(MassBp)), "MassBp")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("mHdInert2", MODELPARAMETER(mH1I2), "mH1I2", model.get_scale());
   slha_io.set_block("mHuInert2", MODELPARAMETER(mH2I2), "mH2I2", model.get_scale());
   slha_io.set_block("mX2", MODELPARAMETER(mDx2), "mDx2", model.get_scale());
   slha_io.set_block("mXBar2", MODELPARAMETER(mDxbar2), "mDxbar2", model.get_scale());
   slha_io.set_block("msInert2", MODELPARAMETER(msI2), "msI2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(102, (MODELPARAMETER(vd)), "vd")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(vu)), "vu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block ESIXRUN Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(11, (MODELPARAMETER(vs)), "vs")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Lambdax)), "Lambdax")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(TLambdax)), "TLambdax")
            << FORMAT_ELEMENT(0, (MODELPARAMETER(MuPr)), "MuPr")
            << FORMAT_ELEMENT(101, (MODELPARAMETER(BMuPr)), "BMuPr")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXTKAPPA", MODELPARAMETER(TKappa), "TKappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());
   slha_io.set_block("ESIXTLAMBDA", MODELPARAMETER(TLambda12), "TLambda12", model.get_scale());

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void CE6SSM_slha_io::set_extra(const CE6SSM<T>& model)
{
   CE6SSM_physical physical(model.get_physical());
   convert_to_slha_convention(physical);


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class
 */
template <class T>
void CE6SSM_slha_io::set_spectrum(const CE6SSM<T>& model)
{
   CE6SSM_physical physical(model.get_physical());
   convert_to_slha_convention(physical);

   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);
}

} // namespace flexiblesusy

#endif
