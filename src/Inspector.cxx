#include "Cuts/Inspector.hxx"

#include "Cuts/Default.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*  */
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelA::MemFn expression, Double_t value, Limit limit_type) {
    Cut<Sexaquark::ChannelA> cut(name, expression, value, limit_type);
    fCuts_ChannelA.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelD::MemFn expression, Double_t value, Limit limit_type) {
    Cut<Sexaquark::ChannelD> cut(name, expression, value, limit_type);
    fCuts_ChannelD.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelE::MemFn expression, Double_t value, Limit limit_type) {
    Cut<Sexaquark::ChannelE> cut(name, expression, value, limit_type);
    fCuts_ChannelE.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::KaonPair::MemFn expression, Double_t value, Limit limit_type) {
    Cut<Sexaquark::KaonPair> cut(name, expression, value, limit_type);
    fCuts_KaonPair.push_back(cut);
}
/*  */
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelA::MemFn expression, Double_t min, Double_t max) {
    Cut<Sexaquark::ChannelA> cut(name, expression, min, max);
    fCuts_ChannelA.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelD::MemFn expression, Double_t min, Double_t max) {
    Cut<Sexaquark::ChannelD> cut(name, expression, min, max);
    fCuts_ChannelD.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::ChannelE::MemFn expression, Double_t min, Double_t max) {
    Cut<Sexaquark::ChannelE> cut(name, expression, min, max);
    fCuts_ChannelE.push_back(cut);
}
template <>
void Inspector::AddCut(TString name, Sexaquark::KaonPair::MemFn expression, Double_t min, Double_t max) {
    Cut<Sexaquark::KaonPair> cut(name, expression, min, max);
    fCuts_KaonPair.push_back(cut);
}

/*
 *
 */
void Inspector::SetLambdaDefaultCuts() {
    //
    AddCut(V0::Species::Lambda, "Eta", &Particle::V0::Eta, Default::Lambda::AbsMaxEta, Limit::AbsoluteMax);
    AddCut(V0::Species::Lambda, "Pt", &Particle::V0::Pt, Default::Lambda::MinPt, Limit::Minimum);
    AddCut(V0::Species::Lambda, "Mass", &Particle::V0::Mass, Default::Lambda::MinMass, Default::Lambda::MaxMass);
    AddCut(V0::Species::Lambda, "Radius", &Particle::V0::Radius, Default::Lambda::MinRadius, Limit::Minimum);
    AddCut(V0::Species::Lambda, "DistFromPV", &Particle::V0::DistFromPV, Default::Lambda::MinDistFromPV, Limit::Minimum);
    AddCut(V0::Species::Lambda, "CPAwrtPV", &Particle::V0::CPAwrtPV, Default::Lambda::MinCPAwrtPV, Default::Lambda::MaxCPAwrtPV);
    AddCut(V0::Species::Lambda, "DCAwrtPV", &Particle::V0::DCAwrtPV, Default::Lambda::MinDCAwrtPV, Limit::Minimum);
    AddCut(V0::Species::Lambda, "DCAbtwDau", &Particle::V0::DCAbtwDau, Default::Lambda::MaxDCAbtwDau, Limit::Maximum);
    AddCut(V0::Species::Lambda, "DCAnegV0", &Particle::V0::DCAnegV0, Default::Lambda::MaxDCAnegV0, Limit::Maximum);
    AddCut(V0::Species::Lambda, "DCAposV0", &Particle::V0::DCAposV0, Default::Lambda::MaxDCAposV0, Limit::Maximum);
    AddCut(V0::Species::Lambda, "ArmQtOverAlpha", &Particle::V0::ArmQtOverAlpha, Default::Lambda::MaxArmQtOverAlpha, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetKaonZeroDefaultCuts() {
    //
    AddCut(V0::Species::KaonZeroShort, "Eta", &Particle::V0::Eta, Default::KaonZero::AbsMaxEta, Limit::AbsoluteMax);
    AddCut(V0::Species::KaonZeroShort, "Pt", &Particle::V0::Pt, Default::KaonZero::MinPt, Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "Mass", &Particle::V0::Mass, Default::KaonZero::MinMass, Default::KaonZero::MaxMass);
    AddCut(V0::Species::KaonZeroShort, "Radius", &Particle::V0::Radius, Default::KaonZero::MinRadius, Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "DistFromPV", &Particle::V0::DistFromPV, Default::KaonZero::MinDistFromPV, Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "CPAwrtPV", &Particle::V0::CPAwrtPV, Default::KaonZero::MinCPAwrtPV, Default::KaonZero::MaxCPAwrtPV);
    AddCut(V0::Species::KaonZeroShort, "DCAwrtPV", &Particle::V0::DCAwrtPV, Default::KaonZero::MinDCAwrtPV, Limit::Minimum);
    AddCut(V0::Species::KaonZeroShort, "DCAbtwDau", &Particle::V0::DCAbtwDau, Default::KaonZero::MaxDCAbtwDau, Limit::Maximum);
    AddCut(V0::Species::KaonZeroShort, "DCAnegV0", &Particle::V0::DCAnegV0, Default::KaonZero::MaxDCAnegV0, Limit::Maximum);
    AddCut(V0::Species::KaonZeroShort, "DCAposV0", &Particle::V0::DCAposV0, Default::KaonZero::MaxDCAposV0, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetPionPairDefaultCuts() {
    //
    AddCut(V0::Species::PionPair, "Eta", &Particle::V0::Eta, Default::PionPair::AbsMaxEta, Limit::AbsoluteMax);
    AddCut(V0::Species::PionPair, "Pt", &Particle::V0::Pt, Default::PionPair::MinPt, Limit::Minimum);
    AddCut(V0::Species::PionPair, "Radius", &Particle::V0::Radius, Default::PionPair::MinRadius, Limit::Minimum);
    AddCut(V0::Species::PionPair, "CPAwrtPV", &Particle::V0::CPAwrtPV, Default::PionPair::MinCPAwrtPV, Default::PionPair::MaxCPAwrtPV);
    AddCut(V0::Species::PionPair, "DCAwrtPV", &Particle::V0::DCAwrtPV, Default::PionPair::MinDCAwrtPV, Limit::Minimum);
    AddCut(V0::Species::PionPair, "DCAbtwDau", &Particle::V0::DCAbtwDau, Default::PionPair::MaxDCAbtwDau, Limit::Maximum);
    AddCut(V0::Species::PionPair, "DCAnegV0", &Particle::V0::DCAnegV0, Default::PionPair::MaxDCAnegV0, Limit::Maximum);
    AddCut(V0::Species::PionPair, "DCAposV0", &Particle::V0::DCAposV0, Default::PionPair::MaxDCAposV0, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetSexaquarkDefaultCuts_ChannelA() {
    //
    /* -- Common */
    AddCut("Radius", &Sexaquark::ChannelA::Radius, Default::ChannelA::MinRadius, Limit::Minimum);
    AddCut("Rapidity", &Sexaquark::ChannelA::Rapidity, Default::ChannelA::MaxRapidity, Limit::Maximum);
    AddCut("CPAwrtPV", &Sexaquark::ChannelA::CPAwrtPV, Default::ChannelA::MinCPAwrtPV, Default::ChannelA::MaxCPAwrtPV);
    AddCut("DCALaSV", &Sexaquark::ChannelA::DCALaSV, Default::ChannelA::MaxDCALaSV, Limit::Maximum);
    AddCut("DCALaNegSV", &Sexaquark::ChannelA::DCALaNegSV, Default::ChannelA::MaxDCALaNegSV, Limit::Maximum);
    AddCut("DCALaPosSV", &Sexaquark::ChannelA::DCALaPosSV, Default::ChannelA::MaxDCALaPosSV, Limit::Maximum);
    AddCut("DCAK0SV", &Sexaquark::ChannelA::DCAK0SV, Default::ChannelA::MaxDCAK0SV, Limit::Maximum);
    AddCut("DCAK0NegSV", &Sexaquark::ChannelA::DCAK0NegSV, Default::ChannelA::MaxDCAK0NegSV, Limit::Maximum);
    AddCut("DCAK0PosSV", &Sexaquark::ChannelA::DCAK0PosSV, Default::ChannelA::MaxDCAK0PosSV, Limit::Maximum);
    AddCut("DCAbtwV0s", &Sexaquark::ChannelA::DCAbtwV0s, Default::ChannelA::MaxDCAbtwV0s, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetSexaquarkDefaultCuts_ChannelD() {
    //
    AddCut("Radius", &Sexaquark::ChannelD::Radius, Default::ChannelD::MinRadius, Limit::Minimum);
    AddCut("Rapidity", &Sexaquark::ChannelD::Rapidity, Default::ChannelD::MaxRapidity, Limit::Maximum);
    AddCut("CPAwrtPV", &Sexaquark::ChannelD::CPAwrtPV, Default::ChannelD::MinCPAwrtPV, Default::ChannelD::MaxCPAwrtPV);
    AddCut("DCALaSV", &Sexaquark::ChannelD::DCALaSV, Default::ChannelD::MaxDCALaSV, Limit::Maximum);
    AddCut("DCALaNegSV", &Sexaquark::ChannelD::DCALaNegSV, Default::ChannelD::MaxDCALaNegSV, Limit::Maximum);
    AddCut("DCALaPosSV", &Sexaquark::ChannelD::DCALaPosSV, Default::ChannelD::MaxDCALaPosSV, Limit::Maximum);
    AddCut("DCAKaSV", &Sexaquark::ChannelD::DCAKaSV, Default::ChannelD::MaxDCAKaSV, Limit::Maximum);
    AddCut("DCAKaLa", &Sexaquark::ChannelD::DCAKaLa, Default::ChannelD::MaxDCAKaLa, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetSexaquarkDefaultCuts_ChannelE() {
    //
    AddCut("Radius", &Sexaquark::ChannelE::Radius, Default::ChannelE::MinRadius, Limit::Minimum);
    AddCut("Rapidity", &Sexaquark::ChannelE::Rapidity, Default::ChannelE::MaxRapidity, Limit::Maximum);
    AddCut("CPAwrtPV", &Sexaquark::ChannelE::CPAwrtPV, Default::ChannelE::MinCPAwrtPV, Default::ChannelE::MaxCPAwrtPV);
    AddCut("DCAKaSV", &Sexaquark::ChannelE::DCAKaSV, Default::ChannelE::MaxDCAKaSV, Limit::Maximum);
    AddCut("DCAKaLa", &Sexaquark::ChannelE::DCAKaLa, Default::ChannelE::MaxDCAKaLa, Limit::Maximum);
    AddCut("DCApmSV", &Sexaquark::ChannelE::DCApmSV, Default::ChannelE::MaxDCApmSV, Limit::Maximum);
    AddCut("DCAppSV", &Sexaquark::ChannelE::DCAppSV, Default::ChannelE::MaxDCAppSV, Limit::Maximum);
    AddCut("DCApmLa", &Sexaquark::ChannelE::DCApmLa, Default::ChannelE::MaxDCApmLa, Limit::Maximum);
    AddCut("DCAppLa", &Sexaquark::ChannelE::DCAppLa, Default::ChannelE::MaxDCAppLa, Limit::Maximum);
    AddCut("DCApmKa", &Sexaquark::ChannelE::DCApmKa, Default::ChannelE::MaxDCApmKa, Limit::Maximum);
    AddCut("DCAppKa", &Sexaquark::ChannelE::DCAppKa, Default::ChannelE::MaxDCAppKa, Limit::Maximum);
}

/*
 *
 */
void Inspector::SetKaonPairDefaultCuts() {
    //
    AddCut("Radius", &Sexaquark::KaonPair::Radius, Default::KaonPair::MinRadius, Limit::Minimum);
    AddCut("Rapidity", &Sexaquark::KaonPair::Rapidity, Default::KaonPair::MaxRapidity, Limit::Maximum);
    AddCut("CPAwrtPV", &Sexaquark::KaonPair::CPAwrtPV, Default::KaonPair::MinCPAwrtPV, Default::KaonPair::MaxCPAwrtPV);
    AddCut("DCAbtwKK", &Sexaquark::KaonPair::DCAbtwKK, Default::KaonPair::MaxDCAbtwKK, Limit::Maximum);
    AddCut("DCAkaSV", &Sexaquark::KaonPair::DCAkaSV, Default::KaonPair::MaxDCAkaSV, Limit::Maximum);
    AddCut("DCAkbSV", &Sexaquark::KaonPair::DCAkbSV, Default::KaonPair::MaxDCAkbSV, Limit::Maximum);
}

/*
 *
 */
void Inspector::PrintAllCuts() {
    //
    InfoF("Lambda cuts: %s", "");
    for (auto cut : fCuts_Lambda) cut.Print();
    InfoF("KaonZero cuts: %s", "");
    for (auto cut : fCuts_KaonZero) cut.Print();
    InfoF("PionPair cuts: %s", "");
    for (auto cut : fCuts_PionPair) cut.Print();
    InfoF("ChannelA cuts: %s", "");
    for (auto cut : fCuts_ChannelA) cut.Print();
    InfoF("ChannelD cuts: %s", "");
    for (auto cut : fCuts_ChannelD) cut.Print();
    InfoF("ChannelE cuts: %s", "");
    for (auto cut : fCuts_ChannelE) cut.Print();
    InfoF("KaonPair cuts: %s", "");
    for (auto cut : fCuts_KaonPair) cut.Print();
}

}  // namespace Cuts
}  // namespace Tree2Sexaquark
