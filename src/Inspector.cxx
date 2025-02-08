#include "Cuts/Inspector.hxx"

#include "Utilities/Logger.hxx"

#include "Cuts/Default.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

void Inspector::InitDefaultCuts_V0() {
    //
    AddCut(V0::Species::Lambda, "EtaCut", Default::Lambda::AbsMaxEta, Cut::Limit::AbsoluteMax);
    AddCut(V0::Species::Lambda, "PtCut", Default::Lambda::MinPt, Cut::Limit::Minimum);
    AddCut(V0::Species::Lambda, "MassRange", Default::Lambda::MinMass, Default::Lambda::MaxMass);
    AddCut(V0::Species::Lambda, "RadiusCut", Default::Lambda::MinRadius, Cut::Limit::Minimum);
    AddCut(V0::Species::Lambda, "DistFromPVCut", Default::Lambda::MinDistFromPV, Cut::Limit::Minimum);
    AddCut(V0::Species::Lambda, "CPAwrtPVCut", Default::Lambda::MinCPAwrtPV, Default::Lambda::MaxCPAwrtPV);
    AddCut(V0::Species::Lambda, "DCAwrtPVCut", Default::Lambda::MinDCAwrtPV, Cut::Limit::Minimum);
    AddCut(V0::Species::Lambda, "DCAbtwDauCut", Default::Lambda::MaxDCAbtwDau, Cut::Limit::Maximum);
    AddCut(V0::Species::Lambda, "DCAnegV0Cut", Default::Lambda::MaxDCAnegV0, Cut::Limit::Maximum);
    AddCut(V0::Species::Lambda, "DCAposV0Cut", Default::Lambda::MaxDCAposV0, Cut::Limit::Maximum);
    AddCut(V0::Species::Lambda, "ArmQtOverAlphaCut", Default::Lambda::MaxArmQtOverAlpha, Cut::Limit::Maximum);
    //
    AddCut(V0::Species::KaonZeroShort, "EtaCut", Default::KaonZeroShort::AbsMaxEta, Cut::Limit::AbsoluteMax);
    AddCut(V0::Species::KaonZeroShort, "PtCut", Default::KaonZeroShort::MinPt, Cut::Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "MassRange", Default::KaonZeroShort::MinMass, Default::KaonZeroShort::MaxMass);
    AddCut(V0::Species::KaonZeroShort, "RadiusCut", Default::KaonZeroShort::MinRadius, Cut::Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "DistFromPVCut", Default::KaonZeroShort::MinDistFromPV, Default::KaonZeroShort::MaxDistFromPV);
    AddCut(V0::Species::KaonZeroShort, "CPAwrtPVCut", Default::KaonZeroShort::MinCPAwrtPV, Default::KaonZeroShort::MaxCPAwrtPV);
    AddCut(V0::Species::KaonZeroShort, "DCAwrtPVCut", Default::KaonZeroShort::MinDCAwrtPV, Cut::Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "DCAbtwDauCut", Default::KaonZeroShort::MaxDCAbtwDau, Cut::Limit::Maximum);
    AddCut(V0::Species::KaonZeroShort, "DCAnegV0Cut", Default::KaonZeroShort::MaxDCAnegV0, Cut::Limit::Maximum);
    AddCut(V0::Species::KaonZeroShort, "DCAposV0Cut", Default::KaonZeroShort::MaxDCAposV0, Cut::Limit::Maximum);
    //
    AddCut(V0::Species::PionPair, "EtaCut", Default::PionPair::AbsMaxEta, Cut::Limit::AbsoluteMax);
    AddCut(V0::Species::PionPair, "PtCut", Default::PionPair::MinPt, Cut::Limit::Minimum);
    AddCut(V0::Species::PionPair, "RadiusCut", Default::PionPair::MinRadius, Cut::Limit::Minimum);
    AddCut(V0::Species::PionPair, "CPAwrtPVCut", Default::PionPair::MinCPAwrtPV, Default::PionPair::MaxCPAwrtPV);
    AddCut(V0::Species::PionPair, "DCAwrtPVCut", Default::PionPair::MinDCAwrtPV, Cut::Limit::Minimum);
    AddCut(V0::Species::PionPair, "DCAbtwDauCut", Default::PionPair::MaxDCAbtwDau, Cut::Limit::Maximum);
    AddCut(V0::Species::PionPair, "DCAnegV0Cut", Default::PionPair::MaxDCAnegV0, Cut::Limit::Maximum);
    AddCut(V0::Species::PionPair, "DCAposV0Cut", Default::PionPair::MaxDCAposV0, Cut::Limit::Maximum);
}

void Inspector::InitDefaultCuts_Sexaquark() {
    //
    AddCut(Sexaquark::Channel::A, "Radius", Default::ChannelA::MinRadius, Cut::Limit::Minimum);
    AddCut(Sexaquark::Channel::A, "Rapidity", Default::ChannelA::MaxRapidity, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "CPAwrtPV", Default::ChannelA::MinCPAwrtPV, Default::ChannelA::MaxCPAwrtPV);
    AddCut(Sexaquark::Channel::A, "DCALaSV", Default::ChannelA::MaxDCALaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCALaNegSV", Default::ChannelA::MaxDCALaNegSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCALaPosSV", Default::ChannelA::MaxDCALaPosSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCAK0SV", Default::ChannelA::MaxDCAK0SV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCAK0NegSV", Default::ChannelA::MaxDCAK0NegSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCAK0PosSV", Default::ChannelA::MaxDCAK0PosSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::A, "DCAbtwV0s", Default::ChannelA::MaxDCAbtwV0s, Cut::Limit::Maximum);
    //
    AddCut(Sexaquark::Channel::D, "Radius", Default::ChannelD::MinRadius, Cut::Limit::Minimum);
    AddCut(Sexaquark::Channel::D, "Rapidity", Default::ChannelD::MaxRapidity, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::D, "CPAwrtPV", Default::ChannelD::MinCPAwrtPV, Default::ChannelD::MaxCPAwrtPV);
    AddCut(Sexaquark::Channel::D, "DCALaSV", Default::ChannelD::MaxDCALaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::D, "DCALaNegSV", Default::ChannelD::MaxDCALaNegSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::D, "DCALaPosSV", Default::ChannelD::MaxDCALaPosSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::D, "DCAKaSV", Default::ChannelD::MaxDCAKaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::D, "DCAKaLa", Default::ChannelD::MaxDCAKaLa, Cut::Limit::Maximum);
    //
    AddCut(Sexaquark::Channel::E, "Radius", Default::ChannelE::MinRadius, Cut::Limit::Minimum);
    AddCut(Sexaquark::Channel::E, "Rapidity", Default::ChannelE::MaxRapidity, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "CPAwrtPV", Default::ChannelE::MinCPAwrtPV, Default::ChannelE::MaxCPAwrtPV);
    AddCut(Sexaquark::Channel::E, "DCALaSV", Default::ChannelE::MaxDCALaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCALaNegSV", Default::ChannelE::MaxDCALaNegSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCALaPosSV", Default::ChannelE::MaxDCALaPosSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCAKaSV", Default::ChannelE::MaxDCAKaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCAKaLa", Default::ChannelE::MaxDCAKaLa, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCApmSV", Default::ChannelE::MaxDCApmSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCAppSV", Default::ChannelE::MaxDCAppSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCApmLa", Default::ChannelE::MaxDCApmLa, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCAppLa", Default::ChannelE::MaxDCAppLa, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCApmKa", Default::ChannelE::MaxDCApmKa, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::E, "DCAppKa", Default::ChannelE::MaxDCAppKa, Cut::Limit::Maximum);
    //
    AddCut(Sexaquark::Channel::H, "Radius", Default::KaonPair::MinRadius, Cut::Limit::Minimum);
    AddCut(Sexaquark::Channel::H, "DCAbtwKK", Default::KaonPair::MaxDCAbtwKK, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::H, "DCAkaSV", Default::KaonPair::MaxDCAkaSV, Cut::Limit::Maximum);
    AddCut(Sexaquark::Channel::H, "DCAkbSV", Default::KaonPair::MaxDCAkbSV, Cut::Limit::Maximum);
}

Bool_t Inspector::Approve(Candidate::V0 v0) {
    //
    if (!Check(v0, "EtaCut", v0.Eta())) return kFALSE;
    if (!Check(v0, "PtCut", v0.Pt())) return kFALSE;
    if (!Check(v0, "MassRange", v0.Mass())) return kFALSE;
    if (!Check(v0, "RadiusCut", v0.Radius())) return kFALSE;
    if (!Check(v0, "DistFromPVCut", v0.DistFromPV())) return kFALSE;
    if (!Check(v0, "CPAwrtPVCut", v0.CPAwrtPV())) return kFALSE;
    if (!Check(v0, "DCAwrtPVCut", v0.DCAwrtPV())) return kFALSE;
    if (!Check(v0, "DCAbtwDauCut", v0.DCAbtwDau())) return kFALSE;
    if (!Check(v0, "DCAnegV0Cut", v0.DCAnegV0())) return kFALSE;
    if (!Check(v0, "DCAposV0Cut", v0.DCAposV0())) return kFALSE;
    if (!Check(v0, "ArmQtOverAlphaCut", v0.ArmQtOverAlpha())) return kFALSE;
    return kTRUE;
}

Bool_t Inspector::Approve(Candidate::ChannelA sexaquark) {
    //
    if (!Check(sexaquark, "RadiusCut", sexaquark.Radius())) return kFALSE;
    if (!Check(sexaquark, "RapidityCut", sexaquark.Rapidity())) return kFALSE;
    if (!Check(sexaquark, "CPAwrtPVCut", sexaquark.CPAwrtPV())) return kFALSE;
    if (!Check(sexaquark, "DCALaSVCut", sexaquark.DCALaSV())) return kFALSE;
    if (!Check(sexaquark, "DCALaNegSVCut", sexaquark.DCALaNegSV())) return kFALSE;
    if (!Check(sexaquark, "DCALaPosSVCut", sexaquark.DCALaPosSV())) return kFALSE;
    if (!Check(sexaquark, "DCAK0SVCut", sexaquark.DCAK0SV())) return kFALSE;
    if (!Check(sexaquark, "DCAK0NegSVCut", sexaquark.DCAK0NegSV())) return kFALSE;
    if (!Check(sexaquark, "DCAK0PosSVCut", sexaquark.DCAK0PosSV())) return kFALSE;
    if (!Check(sexaquark, "DCAbtwV0sCut", sexaquark.DCAbtwV0s())) return kFALSE;
    return kTRUE;
}

Bool_t Inspector::Approve(Candidate::ChannelD sexaquark) {
    //
    if (!Check(sexaquark, "RadiusCut", sexaquark.Radius())) return kFALSE;
    if (!Check(sexaquark, "RapidityCut", sexaquark.Rapidity())) return kFALSE;
    if (!Check(sexaquark, "CPAwrtPVCut", sexaquark.CPAwrtPV())) return kFALSE;
    if (!Check(sexaquark, "DCALaSVCut", sexaquark.DCALaSV())) return kFALSE;
    if (!Check(sexaquark, "DCALaNegSVCut", sexaquark.DCALaNegSV())) return kFALSE;
    if (!Check(sexaquark, "DCALaPosSVCut", sexaquark.DCALaPosSV())) return kFALSE;
    if (!Check(sexaquark, "DCAKaSVCut", sexaquark.DCAKaSV())) return kFALSE;
    if (!Check(sexaquark, "DCAKaLaCut", sexaquark.DCAKaLa())) return kFALSE;
    return kTRUE;
}

Bool_t Inspector::Approve(Candidate::ChannelE sexaquark) {
    //
    if (!Check(sexaquark, "RadiusCut", sexaquark.Radius())) return kFALSE;
    if (!Check(sexaquark, "RapidityCut", sexaquark.Rapidity())) return kFALSE;
    if (!Check(sexaquark, "CPAwrtPVCut", sexaquark.CPAwrtPV())) return kFALSE;
    if (!Check(sexaquark, "DCAKaSVCut", sexaquark.DCAKaSV())) return kFALSE;
    if (!Check(sexaquark, "DCAKaLaCut", sexaquark.DCAKaLa())) return kFALSE;
    if (!Check(sexaquark, "DCApmSVCut", sexaquark.DCApmSV())) return kFALSE;
    if (!Check(sexaquark, "DCAppSVCut", sexaquark.DCAppSV())) return kFALSE;
    if (!Check(sexaquark, "DCApmLaCut", sexaquark.DCApmLa())) return kFALSE;
    if (!Check(sexaquark, "DCAppLaCut", sexaquark.DCAppLa())) return kFALSE;
    if (!Check(sexaquark, "DCApmKaCut", sexaquark.DCApmKa())) return kFALSE;
    if (!Check(sexaquark, "DCAppKaCut", sexaquark.DCAppKa())) return kFALSE;
    return kTRUE;
}

Bool_t Inspector::Approve(Candidate::KaonPair sexaquark) {
    //
    if (!Check(sexaquark, "RadiusCut", sexaquark.Radius())) return kFALSE;
    if (!Check(sexaquark, "RapidityCut", sexaquark.Rapidity())) return kFALSE;
    if (!Check(sexaquark, "CPAwrtPVCut", sexaquark.CPAwrtPV())) return kFALSE;
    if (!Check(sexaquark, "DCAbtwKKCut", sexaquark.DCAbtwKK())) return kFALSE;
    if (!Check(sexaquark, "DCAkaSVCut", sexaquark.DCAkaSV())) return kFALSE;
    if (!Check(sexaquark, "DCAkbSVCut", sexaquark.DCAkbSV())) return kFALSE;
    return kTRUE;
}

void Inspector::PrintAllCuts() {
    //
    InfoF("LAMBDA CUTS %s", "");
    for (auto& [cut_name, cut] : Lambda_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("KAON ZERO SHORT CUTS %s", "");
    for (auto& [cut_name, cut] : KaonZeroShort_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("PION PAIR CUTS %s", "");
    for (auto& [cut_name, cut] : PionPair_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("SEXAQUARK CHANNEL A CUTS %s", "");
    for (auto& [cut_name, cut] : ChannelA_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("SEXAQUARK CHANNEL D CUTS %s", "");
    for (auto& [cut_name, cut] : ChannelD_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("SEXAQUARK CHANNEL E CUTS %s", "");
    for (auto& [cut_name, cut] : ChannelE_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
    InfoF("KAON PAIR CUTS %s", "");
    for (auto& [cut_name, cut] : KaonPair_) InfoF(">> %s [%s]: %s", cut_name.c_str(), cut.IsOn() ? "ON" : "OFF", cut.ToString().c_str());
}

}  // namespace Cuts
}  // namespace Tree2Sexaquark
