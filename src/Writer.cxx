#include "Trees/Writer.hxx"

#include "Analysis/Settings.hxx"

namespace Tree2Sexaquark {

/*         */
/**  V0s  **/
/*** === ***/

void Writer::InitV0sBranches() {
    //
    fTree_[TreeName::V0s]->Branch("Idx", &V0.Idx);
    fTree_[TreeName::V0s]->Branch("Neg_EsdIdx", &V0.Neg_EsdIdx);
    fTree_[TreeName::V0s]->Branch("Pos_EsdIdx", &V0.Pos_EsdIdx);
    fTree_[TreeName::V0s]->Branch("PID", &V0.PID);
    fTree_[TreeName::V0s]->Branch("Px", &V0.Px);
    fTree_[TreeName::V0s]->Branch("Py", &V0.Py);
    fTree_[TreeName::V0s]->Branch("Pz", &V0.Pz);
    fTree_[TreeName::V0s]->Branch("E", &V0.E);
    fTree_[TreeName::V0s]->Branch("Xv", &V0.Xv);
    fTree_[TreeName::V0s]->Branch("Yv", &V0.Yv);
    fTree_[TreeName::V0s]->Branch("Zv", &V0.Zv);
    fTree_[TreeName::V0s]->Branch("Neg_Px", &V0.Neg_Px);
    fTree_[TreeName::V0s]->Branch("Neg_Py", &V0.Neg_Py);
    fTree_[TreeName::V0s]->Branch("Neg_Pz", &V0.Neg_Pz);
    fTree_[TreeName::V0s]->Branch("Pos_Px", &V0.Pos_Px);
    fTree_[TreeName::V0s]->Branch("Pos_Py", &V0.Pos_Py);
    fTree_[TreeName::V0s]->Branch("Pos_Pz", &V0.Pos_Pz);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_[TreeName::V0s]->Branch("McIdx", &V0.McIdx);
        fTree_[TreeName::V0s]->Branch("True_PdgCode", &V0.True_PdgCode);
        fTree_[TreeName::V0s]->Branch("IsSecondary", &V0.IsSecondary);
        fTree_[TreeName::V0s]->Branch("IsSignal", &V0.IsSignal);
        fTree_[TreeName::V0s]->Branch("ReactionID", &V0.ReactionID);
        fTree_[TreeName::V0s]->Branch("IsHybrid", &V0.IsHybrid);
    }
}

void Writer::FillV0(UInt_t idx_v0, Candidate::V0 new_v0) {
    //
    V0.Idx = idx_v0;
    V0.Neg_EsdIdx = new_v0.EsdIdxNeg;
    V0.Pos_EsdIdx = new_v0.EsdIdxNeg;
    V0.PID = new_v0.PdgCode;
    V0.Px = (Float_t)new_v0.Px();
    V0.Py = (Float_t)new_v0.Py();
    V0.Pz = (Float_t)new_v0.Pz();
    V0.E = (Float_t)new_v0.E();
    V0.Xv = new_v0.Xv();
    V0.Yv = new_v0.Yv();
    V0.Zv = new_v0.Zv();
    V0.Neg_Px = (Float_t)new_v0.NegPx();
    V0.Neg_Py = (Float_t)new_v0.NegPy();
    V0.Neg_Pz = (Float_t)new_v0.NegPz();
    V0.Pos_Px = (Float_t)new_v0.PosPx();
    V0.Pos_Py = (Float_t)new_v0.PosPy();
    V0.Pos_Pz = (Float_t)new_v0.PosPz();
    /* True Information */
    if (Analysis::Settings::IsMC) {
        V0.McIdx = new_v0.McIdxV0;
        V0.True_PdgCode = new_v0.McPdgCode;
        V0.IsSecondary = new_v0.IsSecondary;
        V0.IsSignal = new_v0.IsSignal;
        V0.ReactionID = new_v0.ReactionID;
        V0.IsHybrid = new_v0.IsHybrid;
    }
    fTree_[TreeName::V0s]->Fill();
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

void Writer::InitSexaquarkBranches_TypeA(TreeName tree_name) {
    //
    SexaquarkA_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALK0)
        Sexaquark = &Sexaquark_ALK0;
    else if (tree_name == TreeName::Sexaquarks_LK0)
        Sexaquark = &Sexaquark_LK0;
    else
        return;
    /* Common Properties */
    fTree_[tree_name]->Branch("Px", &Sexaquark->Px);
    fTree_[tree_name]->Branch("Py", &Sexaquark->Py);
    fTree_[tree_name]->Branch("Pz", &Sexaquark->Pz);
    fTree_[tree_name]->Branch("E", &Sexaquark->E);
    fTree_[tree_name]->Branch("Xv", &Sexaquark->Xv);
    fTree_[tree_name]->Branch("Yv", &Sexaquark->Yv);
    fTree_[tree_name]->Branch("Zv", &Sexaquark->Zv);
    fTree_[tree_name]->Branch("DistFromPV", &Sexaquark->DistFromPV);
    fTree_[tree_name]->Branch("CPAwrtPV", &Sexaquark->CPAwrtPV);
    fTree_[tree_name]->Branch("DCAwrtPV", &Sexaquark->DCAwrtPV);
    fTree_[tree_name]->Branch("Chi2ndf", &Sexaquark->Chi2ndf);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_[tree_name]->Branch("IsSignal", &Sexaquark->IsSignal);
        fTree_[tree_name]->Branch("ReactionID", &Sexaquark->ReactionID);
        fTree_[tree_name]->Branch("IsHybrid", &Sexaquark->IsHybrid);
        fTree_[tree_name]->Branch("NonCombBkg_PdgCode", &Sexaquark->NonCombBkg_PdgCode);
    }
    /* Shared with Channels "A"+"D"+"E" */
    fTree_[tree_name]->Branch("E_asDecay", &Sexaquark->E_asDecay);
    fTree_[tree_name]->Branch("Lambda_Idx", &Sexaquark->Lambda_Idx);
    fTree_[tree_name]->Branch("Lambda_Neg_EsdIdx", &Sexaquark->Lambda_Neg_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_Pos_EsdIdx", &Sexaquark->Lambda_Pos_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_DecayLength", &Sexaquark->Lambda_DecayLength);
    fTree_[tree_name]->Branch("DCALaSV", &Sexaquark->DCALaSV);
    fTree_[tree_name]->Branch("DCALaNegSV", &Sexaquark->DCALaNegSV);
    fTree_[tree_name]->Branch("DCALaPosSV", &Sexaquark->DCALaPosSV);
    /* Shared with Channels "A"+"D"+"H" */
    fTree_[tree_name]->Branch("OpeningAngle", &Sexaquark->OpeningAngle);
    /* Specific to Channel "A" */
    fTree_[tree_name]->Branch("K0S_Idx", &Sexaquark->K0S_Idx);
    fTree_[tree_name]->Branch("K0S_Neg_EsdIdx", &Sexaquark->K0S_Neg_EsdIdx);
    fTree_[tree_name]->Branch("K0S_Pos_EsdIdx", &Sexaquark->K0S_Pos_EsdIdx);
    fTree_[tree_name]->Branch("K0S_DecayLength", &Sexaquark->K0S_DecayLength);
    fTree_[tree_name]->Branch("DCAK0SV", &Sexaquark->DCAK0SV);
    fTree_[tree_name]->Branch("DCAK0NegSV", &Sexaquark->DCAK0NegSV);
    fTree_[tree_name]->Branch("DCAK0PosSV", &Sexaquark->DCAK0PosSV);
    fTree_[tree_name]->Branch("DCAbtwV0s", &Sexaquark->DCAbtwV0s);
}

void Writer::InitSexaquarkBranches_TypeD(TreeName tree_name) {
    //
    SexaquarkD_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALPK)
        Sexaquark = &Sexaquark_ALPK;
    else if (tree_name == TreeName::Sexaquarks_LNK)
        Sexaquark = &Sexaquark_LNK;
    else
        return;
    /* Common Properties */
    fTree_[tree_name]->Branch("Px", &Sexaquark->Px);
    fTree_[tree_name]->Branch("Py", &Sexaquark->Py);
    fTree_[tree_name]->Branch("Pz", &Sexaquark->Pz);
    fTree_[tree_name]->Branch("E", &Sexaquark->E);
    fTree_[tree_name]->Branch("Xv", &Sexaquark->Xv);
    fTree_[tree_name]->Branch("Yv", &Sexaquark->Yv);
    fTree_[tree_name]->Branch("Zv", &Sexaquark->Zv);
    fTree_[tree_name]->Branch("DistFromPV", &Sexaquark->DistFromPV);
    fTree_[tree_name]->Branch("CPAwrtPV", &Sexaquark->CPAwrtPV);
    fTree_[tree_name]->Branch("DCAwrtPV", &Sexaquark->DCAwrtPV);
    fTree_[tree_name]->Branch("Chi2ndf", &Sexaquark->Chi2ndf);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_[tree_name]->Branch("IsSignal", &Sexaquark->IsSignal);
        fTree_[tree_name]->Branch("ReactionID", &Sexaquark->ReactionID);
        fTree_[tree_name]->Branch("IsHybrid", &Sexaquark->IsHybrid);
        fTree_[tree_name]->Branch("NonCombBkg_PdgCode", &Sexaquark->NonCombBkg_PdgCode);
    }
    /* Shared with Channels "A"+"D"+"E" */
    fTree_[tree_name]->Branch("E_asDecay", &Sexaquark->E_asDecay);
    fTree_[tree_name]->Branch("Lambda_Idx", &Sexaquark->Lambda_Idx);
    fTree_[tree_name]->Branch("Lambda_Neg_EsdIdx", &Sexaquark->Lambda_Neg_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_Pos_EsdIdx", &Sexaquark->Lambda_Pos_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_DecayLength", &Sexaquark->Lambda_DecayLength);
    fTree_[tree_name]->Branch("DCALaSV", &Sexaquark->DCALaSV);
    fTree_[tree_name]->Branch("DCALaNegSV", &Sexaquark->DCALaNegSV);
    fTree_[tree_name]->Branch("DCALaPosSV", &Sexaquark->DCALaPosSV);
    /* Shared with Channels "A"+"D"+"H" */
    fTree_[tree_name]->Branch("OpeningAngle", &Sexaquark->OpeningAngle);
    /* Specific to Channel "D" */
    fTree_[tree_name]->Branch("Kaon_EsdIdx", &Sexaquark->Kaon_EsdIdx);
    fTree_[tree_name]->Branch("DCAKaSV", &Sexaquark->DCAKaSV);
    fTree_[tree_name]->Branch("DCAKaLa", &Sexaquark->DCAKaLa);
    fTree_[tree_name]->Branch("DCALaNegKa", &Sexaquark->DCALaNegKa);
    fTree_[tree_name]->Branch("DCALaPosKa", &Sexaquark->DCALaPosKa);
}

void Writer::InitSexaquarkBranches_TypeE(TreeName tree_name) {
    //
    SexaquarkE_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALPKPP)
        Sexaquark = &Sexaquark_ALPKPP;
    else if (tree_name == TreeName::Sexaquarks_LNKPP)
        Sexaquark = &Sexaquark_LNKPP;
    else
        return;
    /* Common Properties */
    fTree_[tree_name]->Branch("Px", &Sexaquark->Px);
    fTree_[tree_name]->Branch("Py", &Sexaquark->Py);
    fTree_[tree_name]->Branch("Pz", &Sexaquark->Pz);
    fTree_[tree_name]->Branch("E", &Sexaquark->E);
    fTree_[tree_name]->Branch("Xv", &Sexaquark->Xv);
    fTree_[tree_name]->Branch("Yv", &Sexaquark->Yv);
    fTree_[tree_name]->Branch("Zv", &Sexaquark->Zv);
    fTree_[tree_name]->Branch("DistFromPV", &Sexaquark->DistFromPV);
    fTree_[tree_name]->Branch("CPAwrtPV", &Sexaquark->CPAwrtPV);
    fTree_[tree_name]->Branch("DCAwrtPV", &Sexaquark->DCAwrtPV);
    fTree_[tree_name]->Branch("Chi2ndf", &Sexaquark->Chi2ndf);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_[tree_name]->Branch("IsSignal", &Sexaquark->IsSignal);
        fTree_[tree_name]->Branch("ReactionID", &Sexaquark->ReactionID);
        fTree_[tree_name]->Branch("IsHybrid", &Sexaquark->IsHybrid);
        fTree_[tree_name]->Branch("NonCombBkg_PdgCode", &Sexaquark->NonCombBkg_PdgCode);
    }
    /* Shared with Channels "A"+"D"+"E" */
    fTree_[tree_name]->Branch("E_asDecay", &Sexaquark->E_asDecay);
    fTree_[tree_name]->Branch("Lambda_Idx", &Sexaquark->Lambda_Idx);
    fTree_[tree_name]->Branch("Lambda_Neg_EsdIdx", &Sexaquark->Lambda_Neg_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_Pos_EsdIdx", &Sexaquark->Lambda_Pos_EsdIdx);
    fTree_[tree_name]->Branch("Lambda_DecayLength", &Sexaquark->Lambda_DecayLength);
    fTree_[tree_name]->Branch("DCALaSV", &Sexaquark->DCALaSV);
    fTree_[tree_name]->Branch("DCALaNegSV", &Sexaquark->DCALaNegSV);
    fTree_[tree_name]->Branch("DCALaPosSV", &Sexaquark->DCALaPosSV);
    /* Specific to Channel "E" */
    fTree_[tree_name]->Branch("Kaon_EsdIdx", &Sexaquark->Kaon_EsdIdx);
    fTree_[tree_name]->Branch("PionPair_Idx", &Sexaquark->PionPair_Idx);
    fTree_[tree_name]->Branch("PiMinus_EsdIdx", &Sexaquark->PiMinus_EsdIdx);
    fTree_[tree_name]->Branch("PiPlus_EsdIdx", &Sexaquark->PiPlus_EsdIdx);
    fTree_[tree_name]->Branch("DCAKaSV", &Sexaquark->DCAKaSV);
    fTree_[tree_name]->Branch("DCAKaLa", &Sexaquark->DCAKaLa);
    fTree_[tree_name]->Branch("DCApmSV", &Sexaquark->DCApmSV);
    fTree_[tree_name]->Branch("DCAppSV", &Sexaquark->DCAppSV);
    fTree_[tree_name]->Branch("DCApmLa", &Sexaquark->DCApmLa);
    fTree_[tree_name]->Branch("DCAppLa", &Sexaquark->DCAppLa);
    fTree_[tree_name]->Branch("DCApmKa", &Sexaquark->DCApmKa);
    fTree_[tree_name]->Branch("DCAppKa", &Sexaquark->DCAppKa);
}

void Writer::InitKaonPairBranches(TreeName tree_name) {
    //
    KaonPair_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_PKPKX)
        Sexaquark = &Sexaquark_PKPKX;
    else if (tree_name == TreeName::Sexaquarks_NKNKX)
        Sexaquark = &Sexaquark_NKNKX;
    else
        return;
    /* Common Properties */
    fTree_[tree_name]->Branch("Px", &Sexaquark->Px);
    fTree_[tree_name]->Branch("Py", &Sexaquark->Py);
    fTree_[tree_name]->Branch("Pz", &Sexaquark->Pz);
    fTree_[tree_name]->Branch("E", &Sexaquark->E);
    fTree_[tree_name]->Branch("Xv", &Sexaquark->Xv);
    fTree_[tree_name]->Branch("Yv", &Sexaquark->Yv);
    fTree_[tree_name]->Branch("Zv", &Sexaquark->Zv);
    fTree_[tree_name]->Branch("DistFromPV", &Sexaquark->DistFromPV);
    fTree_[tree_name]->Branch("CPAwrtPV", &Sexaquark->CPAwrtPV);
    fTree_[tree_name]->Branch("DCAwrtPV", &Sexaquark->DCAwrtPV);
    fTree_[tree_name]->Branch("Chi2ndf", &Sexaquark->Chi2ndf);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_[tree_name]->Branch("IsSignal", &Sexaquark->IsSignal);
        fTree_[tree_name]->Branch("ReactionID", &Sexaquark->ReactionID);
        fTree_[tree_name]->Branch("IsHybrid", &Sexaquark->IsHybrid);
        fTree_[tree_name]->Branch("NonCombBkg_PdgCode", &Sexaquark->NonCombBkg_PdgCode);
    }
    /* Shared with Channels "A"+"D"+"H" */
    fTree_[tree_name]->Branch("OpeningAngle", &Sexaquark->OpeningAngle);
    /* Specific to Channel "H" */
    fTree_[tree_name]->Branch("KaonA_EsdIdx", &Sexaquark->KaonA_EsdIdx);
    fTree_[tree_name]->Branch("KaonB_EsdIdx", &Sexaquark->KaonB_EsdIdx);
    fTree_[tree_name]->Branch("DCAbtwKK", &Sexaquark->DCAbtwKK);
    fTree_[tree_name]->Branch("DCAkaSV", &Sexaquark->DCAkaSV);
    fTree_[tree_name]->Branch("DCAkbSV", &Sexaquark->DCAkbSV);
}

void Writer::FillSexaquark(TreeName tree_name, Candidate::ChannelA new_sexaquark) {
    //
    SexaquarkA_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALK0)
        Sexaquark = &Sexaquark_ALK0;
    else if (tree_name == TreeName::Sexaquarks_LK0)
        Sexaquark = &Sexaquark_LK0;
    else
        return;
    /* Common Properties */
    Sexaquark->Px = new_sexaquark.Px();
    Sexaquark->Py = new_sexaquark.Py();
    Sexaquark->Pz = new_sexaquark.Pz();
    Sexaquark->E = new_sexaquark.E();
    Sexaquark->Xv = new_sexaquark.Xv();
    Sexaquark->Yv = new_sexaquark.Yv();
    Sexaquark->Zv = new_sexaquark.Zv();
    Sexaquark->DistFromPV = new_sexaquark.DistFromPV();
    Sexaquark->CPAwrtPV = new_sexaquark.CPAwrtPV();
    Sexaquark->DCAwrtPV = new_sexaquark.DCAwrtPV();
    Sexaquark->Chi2ndf = new_sexaquark.Chi2ndf();
    /* True Information */
    Sexaquark->IsSignal = new_sexaquark.IsSignal;
    Sexaquark->ReactionID = new_sexaquark.ReactionID;
    Sexaquark->IsHybrid = new_sexaquark.IsHybrid;
    Sexaquark->NonCombBkg_PdgCode = new_sexaquark.NonCombBkg_PdgCode;
    /* Shared with Channels "A"+"D"+"E" */
    Sexaquark->E_asDecay = new_sexaquark.E_asDecay();
    Sexaquark->Lambda_Idx = new_sexaquark.Lambda_Idx;
    Sexaquark->Lambda_Neg_EsdIdx = new_sexaquark.Lambda_Neg_EsdIdx;
    Sexaquark->Lambda_Pos_EsdIdx = new_sexaquark.Lambda_Pos_EsdIdx;
    Sexaquark->Lambda_DecayLength = new_sexaquark.Lambda_DecayLength();
    Sexaquark->DCALaSV = new_sexaquark.DCALaSV();
    Sexaquark->DCALaNegSV = new_sexaquark.DCALaNegSV();
    Sexaquark->DCALaPosSV = new_sexaquark.DCALaPosSV();
    /* Shared with Channels "A"+"D"+"H" */
    Sexaquark->OpeningAngle = new_sexaquark.OpeningAngle();
    /* Specific to Channel "A" */
    Sexaquark->K0S_Idx = new_sexaquark.K0S_Idx;
    Sexaquark->K0S_Neg_EsdIdx = new_sexaquark.K0S_Neg_EsdIdx;
    Sexaquark->K0S_Pos_EsdIdx = new_sexaquark.K0S_Pos_EsdIdx;
    Sexaquark->K0S_DecayLength = new_sexaquark.K0S_DecayLength();
    Sexaquark->DCAK0SV = new_sexaquark.DCAK0SV();
    Sexaquark->DCAK0NegSV = new_sexaquark.DCAK0NegSV();
    Sexaquark->DCAK0PosSV = new_sexaquark.DCAK0PosSV();
    Sexaquark->DCAbtwV0s = new_sexaquark.DCAbtwV0s();
    /* Fill */
    fTree_[tree_name]->Fill();
}

void Writer::FillSexaquark(TreeName tree_name, Candidate::ChannelD new_sexaquark) {
    //
    SexaquarkD_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALPK)
        Sexaquark = &Sexaquark_ALPK;
    else if (tree_name == TreeName::Sexaquarks_LNK)
        Sexaquark = &Sexaquark_LNK;
    else
        return;
    /* Common Properties */
    Sexaquark->Px = new_sexaquark.Px();
    Sexaquark->Py = new_sexaquark.Py();
    Sexaquark->Pz = new_sexaquark.Pz();
    Sexaquark->E = new_sexaquark.E();
    Sexaquark->Xv = new_sexaquark.Xv();
    Sexaquark->Yv = new_sexaquark.Yv();
    Sexaquark->Zv = new_sexaquark.Zv();
    Sexaquark->DistFromPV = new_sexaquark.DistFromPV();
    Sexaquark->CPAwrtPV = new_sexaquark.CPAwrtPV();
    Sexaquark->DCAwrtPV = new_sexaquark.DCAwrtPV();
    Sexaquark->Chi2ndf = new_sexaquark.Chi2ndf();
    /* True Information */
    if (Analysis::Settings::IsMC) {
        Sexaquark->IsSignal = new_sexaquark.IsSignal;
        Sexaquark->ReactionID = new_sexaquark.ReactionID;
        Sexaquark->IsHybrid = new_sexaquark.IsHybrid;
        Sexaquark->NonCombBkg_PdgCode = new_sexaquark.NonCombBkg_PdgCode;
    }
    /* Shared with Channels "A"+"D"+"E" */
    Sexaquark->E_asDecay = new_sexaquark.E_asDecay();
    Sexaquark->Lambda_Idx = new_sexaquark.Lambda_Idx;
    Sexaquark->Lambda_Neg_EsdIdx = new_sexaquark.Lambda_Neg_EsdIdx;
    Sexaquark->Lambda_Pos_EsdIdx = new_sexaquark.Lambda_Pos_EsdIdx;
    Sexaquark->Lambda_DecayLength = new_sexaquark.Lambda_DecayLength();
    Sexaquark->DCALaSV = new_sexaquark.DCALaSV();
    Sexaquark->DCALaNegSV = new_sexaquark.DCALaNegSV();
    Sexaquark->DCALaPosSV = new_sexaquark.DCALaPosSV();
    /* Shared with Channels "A"+"D"+"H" */
    Sexaquark->OpeningAngle = new_sexaquark.OpeningAngle();
    /* Specific to Channel "D" */
    Sexaquark->Kaon_EsdIdx = new_sexaquark.Kaon_EsdIdx;
    Sexaquark->DCAKaSV = new_sexaquark.DCAKaSV();
    Sexaquark->DCAKaLa = new_sexaquark.DCAKaLa();
    Sexaquark->DCALaNegKa = new_sexaquark.DCALaNegKa();
    Sexaquark->DCALaPosKa = new_sexaquark.DCALaPosKa();
    /* Fill */
    fTree_[tree_name]->Fill();
}

void Writer::FillSexaquark(TreeName tree_name, Candidate::ChannelE new_sexaquark) {
    //
    SexaquarkE_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_ALPKPP)
        Sexaquark = &Sexaquark_ALPKPP;
    else if (tree_name == TreeName::Sexaquarks_LNKPP)
        Sexaquark = &Sexaquark_LNKPP;
    else
        return;
    /* Common Properties */
    Sexaquark->Px = new_sexaquark.Px();
    Sexaquark->Py = new_sexaquark.Py();
    Sexaquark->Pz = new_sexaquark.Pz();
    Sexaquark->E = new_sexaquark.E();
    Sexaquark->Xv = new_sexaquark.Xv();
    Sexaquark->Yv = new_sexaquark.Yv();
    Sexaquark->Zv = new_sexaquark.Zv();
    Sexaquark->DistFromPV = new_sexaquark.DistFromPV();
    Sexaquark->CPAwrtPV = new_sexaquark.CPAwrtPV();
    Sexaquark->DCAwrtPV = new_sexaquark.DCAwrtPV();
    Sexaquark->Chi2ndf = new_sexaquark.Chi2ndf();
    /* True Information */
    if (Analysis::Settings::IsMC) {
        Sexaquark->IsSignal = new_sexaquark.IsSignal;
        Sexaquark->ReactionID = new_sexaquark.ReactionID;
        Sexaquark->IsHybrid = new_sexaquark.IsHybrid;
        Sexaquark->NonCombBkg_PdgCode = new_sexaquark.NonCombBkg_PdgCode;
    }
    /* Shared with Channels "A"+"D"+"E" */
    Sexaquark->E_asDecay = new_sexaquark.E_asDecay();
    Sexaquark->Lambda_Idx = new_sexaquark.Lambda_Idx;
    Sexaquark->Lambda_Neg_EsdIdx = new_sexaquark.Lambda_Neg_EsdIdx;
    Sexaquark->Lambda_Pos_EsdIdx = new_sexaquark.Lambda_Pos_EsdIdx;
    Sexaquark->Lambda_DecayLength = new_sexaquark.Lambda_DecayLength();
    Sexaquark->DCALaSV = new_sexaquark.DCALaSV();
    Sexaquark->DCALaNegSV = new_sexaquark.DCALaNegSV();
    Sexaquark->DCALaPosSV = new_sexaquark.DCALaPosSV();
    /* Specific to Channel "E" */
    Sexaquark->Kaon_EsdIdx = new_sexaquark.Kaon_EsdIdx;
    Sexaquark->PionPair_Idx = new_sexaquark.PionPair_Idx;
    Sexaquark->PiMinus_EsdIdx = new_sexaquark.PiMinus_EsdIdx;
    Sexaquark->PiPlus_EsdIdx = new_sexaquark.PiPlus_EsdIdx;
    Sexaquark->DCAKaSV = new_sexaquark.DCAKaSV();
    Sexaquark->DCAKaLa = new_sexaquark.DCAKaLa();
    Sexaquark->DCApmSV = new_sexaquark.DCApmSV();
    Sexaquark->DCAppSV = new_sexaquark.DCAppSV();
    Sexaquark->DCApmLa = new_sexaquark.DCApmLa();
    Sexaquark->DCAppLa = new_sexaquark.DCAppLa();
    Sexaquark->DCApmKa = new_sexaquark.DCApmKa();
    Sexaquark->DCAppKa = new_sexaquark.DCAppKa();
    /* Fill */
    fTree_[tree_name]->Fill();
}

void Writer::FillKaonPair(TreeName tree_name, Candidate::KaonPair new_kk) {
    //
    KaonPair_tt *Sexaquark;
    if (tree_name == TreeName::Sexaquarks_PKPKX)
        Sexaquark = &Sexaquark_PKPKX;
    else if (tree_name == TreeName::Sexaquarks_NKNKX)
        Sexaquark = &Sexaquark_NKNKX;
    else
        return;
    /* Common Properties */
    Sexaquark->Px = new_kk.Px();
    Sexaquark->Py = new_kk.Py();
    Sexaquark->Pz = new_kk.Pz();
    Sexaquark->E = new_kk.E();
    Sexaquark->Xv = new_kk.Xv();
    Sexaquark->Yv = new_kk.Yv();
    Sexaquark->Zv = new_kk.Zv();
    Sexaquark->DistFromPV = new_kk.DistFromPV();
    Sexaquark->CPAwrtPV = new_kk.CPAwrtPV();
    Sexaquark->DCAwrtPV = new_kk.DCAwrtPV();
    Sexaquark->Chi2ndf = new_kk.Chi2ndf();
    /* True Information */
    if (Analysis::Settings::IsMC) {
        Sexaquark->IsSignal = new_kk.IsSignal;
        Sexaquark->ReactionID = new_kk.ReactionID;
        Sexaquark->IsHybrid = new_kk.IsHybrid;
        Sexaquark->NonCombBkg_PdgCode = new_kk.NonCombBkg_PdgCode;
    }
    /* Shared with Channels "A"+"D"+"H" */
    Sexaquark->OpeningAngle = new_kk.OpeningAngle();
    /* Specific to Channel "H" */
    Sexaquark->KaonA_EsdIdx = new_kk.KaonA_EsdIdx;
    Sexaquark->KaonB_EsdIdx = new_kk.KaonB_EsdIdx;
    Sexaquark->DCAbtwKK = new_kk.DCAbtwKK();
    Sexaquark->DCAkaSV = new_kk.DCAkaSV();
    Sexaquark->DCAkbSV = new_kk.DCAkbSV();
    /* Fill */
    fTree_[tree_name]->Fill();
}

}  // namespace Tree2Sexaquark
