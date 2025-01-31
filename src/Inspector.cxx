#include "Cuts/Inspector.hxx"

#include "Cuts/Default.hxx"
#include "Particles/Sexaquark.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*
 *
 */
void Inspector::SetLambdaDefaultCuts() {
    //
    AddCut(
        kLambda, "eta", [](Particle::V0 thisV0) { return (Float_t)thisV0.Eta(); }, Default::Lambda_AbsMaxEta, Cuts::kAbsoluteMax);
    AddCut(
        kLambda, "pt", [](Particle::V0 thisV0) { return (Float_t)thisV0.Pt(); }, Default::Lambda_MinPt, Cuts::kMinimum);
    AddCut(
        kLambda, "mass", [](Particle::V0 thisV0) { return (Float_t)thisV0.Mass(); }, Default::Lambda_MinMass, Default::Lambda_MaxMass);
    AddCut(
        kLambda, "radius", [](Particle::V0 thisV0) { return (Float_t)thisV0.Radius(); }, Default::Lambda_MinRadius, Cuts::kMinimum);
    AddCut(
        kLambda, "dist_from_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.DistFromPV(); }, Default::Lambda_MinDistFromPV, Cuts::kMinimum);
    AddCut(
        kLambda, "cpa_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.CPAwrtPV(); }, Default::Lambda_MinCPAwrtPV,
        Default::Lambda_MaxCPAwrtPV);
    AddCut(
        kLambda, "dca_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAwrtPV(); }, Default::Lambda_MinDCAwrtPV, Cuts::kMinimum);
    AddCut(
        kLambda, "dca_btw_dau", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAbtwDau(); }, Default::Lambda_MaxDCAbtwDau, Cuts::kMaximum);
    AddCut(
        kLambda, "dca_neg_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAnegV0(); }, Default::Lambda_MaxDCAnegV0, Cuts::kMaximum);
    AddCut(
        kLambda, "dca_pos_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAposV0(); }, Default::Lambda_MaxDCAposV0, Cuts::kMaximum);
    AddCut(
        kLambda, "arm_qt_over_alpha", [](Particle::V0 thisV0) { return (Float_t)thisV0.ArmQtOverAlpha(); }, Default::Lambda_MaxArmPtOverAlpha,
        Cuts::kMaximum);
}

/*
 *
 */
void Inspector::SetKaonZeroDefaultCuts() {
    //
    AddCut(
        kKaonZero, "eta", [](Particle::V0 thisV0) { return (Float_t)thisV0.Eta(); }, Default::KaonZero_AbsMaxEta, Cuts::kAbsoluteMax);
    AddCut(
        kKaonZero, "pt", [](Particle::V0 thisV0) { return (Float_t)thisV0.Pt(); }, Default::KaonZero_MinPt, Cuts::kMinimum);
    AddCut(
        kKaonZero, "mass", [](Particle::V0 thisV0) { return (Float_t)thisV0.Mass(); }, Default::KaonZero_MinMass, Default::KaonZero_MaxMass);
    AddCut(
        kKaonZero, "radius", [](Particle::V0 thisV0) { return (Float_t)thisV0.Radius(); }, Default::KaonZero_MinRadius, Cuts::kMinimum);
    AddCut(
        kKaonZero, "dist_from_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.DistFromPV(); }, Default::KaonZero_MinDistFromPV, Cuts::kMinimum);
    AddCut(
        kKaonZero, "cpa_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.CPAwrtPV(); }, Default::KaonZero_MinCPAwrtPV,
        Default::KaonZero_MaxCPAwrtPV);
    AddCut(
        kKaonZero, "dca_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAwrtPV(); }, Default::KaonZero_MinDCAwrtPV, Cuts::kMinimum);
    AddCut(
        kKaonZero, "dca_btw_dau", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAbtwDau(); }, Default::KaonZero_MaxDCAbtwDau, Cuts::kMaximum);
    AddCut(
        kKaonZero, "dca_neg_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAnegV0(); }, Default::KaonZero_MaxDCAnegV0, Cuts::kMaximum);
    AddCut(
        kKaonZero, "dca_pos_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAposV0(); }, Default::KaonZero_MaxDCAposV0, Cuts::kMaximum);
}

/*
 *
 */
void Inspector::SetPionPairDefaultCuts() {
    //
    AddCut(
        kPionPair, "eta", [](Particle::V0 thisV0) { return (Float_t)thisV0.Eta(); }, Default::PionPair_AbsMaxEta, Cuts::kAbsoluteMax);
    AddCut(
        kPionPair, "pt", [](Particle::V0 thisV0) { return (Float_t)thisV0.Pt(); }, Default::PionPair_MinPt, Cuts::kMinimum);
    AddCut(
        kPionPair, "radius", [](Particle::V0 thisV0) { return (Float_t)thisV0.Radius(); }, Default::PionPair_MinRadius, Cuts::kMinimum);
    AddCut(
        kPionPair, "cpa_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.CPAwrtPV(); }, Default::PionPair_MinCPAwrtPV,
        Default::PionPair_MaxCPAwrtPV);
    AddCut(
        kPionPair, "dca_wrt_pv", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAwrtPV(); }, Default::PionPair_MinDCAwrtPV, Cuts::kMinimum);
    AddCut(
        kPionPair, "dca_btw_dau", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAbtwDau(); }, Default::PionPair_MaxDCAbtwDau, Cuts::kMaximum);
    AddCut(
        kPionPair, "dca_neg_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAnegV0(); }, Default::PionPair_MaxDCAnegV0, Cuts::kMaximum);
    AddCut(
        kPionPair, "dca_pos_v0", [](Particle::V0 thisV0) { return (Float_t)thisV0.DCAposV0(); }, Default::PionPair_MaxDCAposV0, Cuts::kMaximum);
}

/*
 *
 */

void Inspector::SetSexaquarkDefaultCuts_ChannelA() {
    //
}

/*
 *
 */
void Inspector::SetSexaquarkDefaultCuts_ChannelD() {
    //
}

/*
 *
 */
void Inspector::SetSexaquarkDefaultCuts_ChannelE() {
    //
}

/*
 *
 */
void Inspector::SetKaonPairDefaultCuts() {
    //
}

/*
 *
 */
Bool_t Inspector::Approve(Particle::V0& thisV0) {
    //
    std::vector<Cut<Particle::V0>> CutsCollection;
    /*  */
    if (TMath::Abs(thisV0.PdgCode) == 3122) CutsCollection = fCuts_Lambda;
    if (TMath::Abs(thisV0.PdgCode) == 310) CutsCollection = fCuts_KaonZero;
    if (TMath::Abs(thisV0.PdgCode) == 422) CutsCollection = fCuts_PionPair;
    /*  */
    for (Cut<Particle::V0>& cut : CutsCollection) {
        if (!cut.Check(thisV0)) return kFALSE;
    }
    InfoF("Particle approved. Glory to Arstotzka. %s", "");
    return kTRUE;
}

/*
 *
 */
Bool_t Inspector::Approve(Particle::Sexaquark& candidate) {
    //
    std::vector<Cut<Particle::Sexaquark>> CutsCollection;
    /*  */
    if (TMath::Abs(candidate.Channel) == kChannelA) CutsCollection = fCuts_ChannelA;
    if (TMath::Abs(candidate.Channel) == kChannelD) CutsCollection = fCuts_ChannelD;
    if (TMath::Abs(candidate.Channel) == kChannelE) CutsCollection = fCuts_ChannelE;
    if (TMath::Abs(candidate.Channel) == kChannelH) CutsCollection = fCuts_KaonPair;
    /*  */
    for (Cut<Particle::Sexaquark>& cut : CutsCollection) {
        if (!cut.Check(candidate)) return kFALSE;
    }
    InfoF("Particle approved. Glory to Arstotzka. %s", "");
    return kTRUE;
}

/*
 *
 */
void Inspector::PrintAllCuts() {
    //
    InfoF("Lambda cuts: %s", "");
    for (Cut<Particle::V0>& cut : fCuts_Lambda) cut.Print();
    InfoF("KaonZero cuts: %s", "");
    for (Cut<Particle::V0>& cut : fCuts_KaonZero) cut.Print();
    InfoF("PionPair cuts: %s", "");
    for (Cut<Particle::V0>& cut : fCuts_PionPair) cut.Print();
    InfoF("ChannelA cuts: %s", "");
    for (Cut<Particle::Sexaquark>& cut : fCuts_ChannelA) cut.Print();
    InfoF("ChannelD cuts: %s", "");
    for (Cut<Particle::Sexaquark>& cut : fCuts_ChannelD) cut.Print();
    InfoF("ChannelE cuts: %s", "");
    for (Cut<Particle::Sexaquark>& cut : fCuts_ChannelE) cut.Print();
    InfoF("KaonPair cuts: %s", "");
    for (Cut<Particle::Sexaquark>& cut : fCuts_KaonPair) cut.Print();
}

}  // namespace Cuts
}  // namespace Tree2Sexaquark
