#include "Trees/Writer.hxx"

#include "Analysis/Settings.hxx"

namespace Tree2Sexaquark {

/*         */
/**  V0s  **/
/*** === ***/

/*
 *
 */
void Writer::InitV0sBranches() {
    //
    fTree_V0s->Branch("Idx", &V0.Idx);
    fTree_V0s->Branch("Idx_Neg", &V0.Idx_Neg);
    fTree_V0s->Branch("Idx_Pos", &V0.Idx_Pos);
    fTree_V0s->Branch("PID", &V0.PID);
    fTree_V0s->Branch("Px", &V0.Px);
    fTree_V0s->Branch("Py", &V0.Py);
    fTree_V0s->Branch("Pz", &V0.Pz);
    fTree_V0s->Branch("E", &V0.E);
    fTree_V0s->Branch("Xv", &V0.Xv);
    fTree_V0s->Branch("Yv", &V0.Yv);
    fTree_V0s->Branch("Zv", &V0.Zv);
    fTree_V0s->Branch("Neg_Px", &V0.Neg_Px);
    fTree_V0s->Branch("Neg_Py", &V0.Neg_Py);
    fTree_V0s->Branch("Neg_Pz", &V0.Neg_Pz);
    fTree_V0s->Branch("Pos_Px", &V0.Pos_Px);
    fTree_V0s->Branch("Pos_Py", &V0.Pos_Py);
    fTree_V0s->Branch("Pos_Pz", &V0.Pos_Pz);
    /* True Information */
    if (Analysis::Settings::IsMC) {
        fTree_V0s->Branch("Idx_True", &V0.Idx_True);
        fTree_V0s->Branch("True_PdgCode", &V0.True_PdgCode);
        fTree_V0s->Branch("IsSecondary", &V0.IsSecondary);
        fTree_V0s->Branch("IsSignal", &V0.IsSignal);
        fTree_V0s->Branch("ReactionID", &V0.ReactionID);
        fTree_V0s->Branch("IsHybrid", &V0.IsHybrid);
    }
}

/*
 *
 */
void Writer::FillV0(UInt_t idxV0, Particle::V0 thisV0) {
    //
    V0.Idx = idxV0;
    V0.Idx_Neg = thisV0.EsdIdxNeg;
    V0.Idx_Pos = thisV0.EsdIdxNeg;
    V0.PID = thisV0.PdgCode;
    V0.Px = (Float_t)thisV0.Px();
    V0.Py = (Float_t)thisV0.Py();
    V0.Pz = (Float_t)thisV0.Pz();
    V0.E = (Float_t)thisV0.E();
    V0.Xv = thisV0.Xv();
    V0.Yv = thisV0.Yv();
    V0.Zv = thisV0.Zv();
    V0.Neg_Px = (Float_t)thisV0.NegPx();
    V0.Neg_Py = (Float_t)thisV0.NegPy();
    V0.Neg_Pz = (Float_t)thisV0.NegPz();
    V0.Pos_Px = (Float_t)thisV0.PosPx();
    V0.Pos_Py = (Float_t)thisV0.PosPy();
    V0.Pos_Pz = (Float_t)thisV0.PosPz();
    if (Analysis::Settings::IsMC) {
        V0.Idx_True = thisV0.McIdxV0;
        V0.True_PdgCode = thisV0.McPdgCode;
        V0.IsSecondary = thisV0.IsSecondary;
        V0.IsSignal = thisV0.IsSignal;
        V0.ReactionID = thisV0.ReactionID;
        V0.IsHybrid = thisV0.IsHybrid;
    }
    fTree_V0s->Fill();
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/*
 *
 */
void Writer::InitSexaquarkBranches_TypeA() {
    //
    std::vector<SexaquarkA_tt> Sexaquark = {Sexaquark_ALK0, Sexaquark_LK0};
    std::vector<TTree*> Tree_Sexaquarks = {fTree_Sexaquarks_ALK0, fTree_Sexaquarks_LK0};
    for (Int_t i = 0; i < 2; i++) {
        /* Common Properties */
        Tree_Sexaquarks[i]->Branch("Px", &Sexaquark[i].Px);
        Tree_Sexaquarks[i]->Branch("Py", &Sexaquark[i].Py);
        Tree_Sexaquarks[i]->Branch("Pz", &Sexaquark[i].Pz);
        Tree_Sexaquarks[i]->Branch("E", &Sexaquark[i].E);
        Tree_Sexaquarks[i]->Branch("E_asDecay", &Sexaquark[i].E_asDecay);
        Tree_Sexaquarks[i]->Branch("Xv", &Sexaquark[i].Xv);
        Tree_Sexaquarks[i]->Branch("Yv", &Sexaquark[i].Yv);
        Tree_Sexaquarks[i]->Branch("Zv", &Sexaquark[i].Zv);
        Tree_Sexaquarks[i]->Branch("DistFromPV", &Sexaquark[i].DistFromPV);
        Tree_Sexaquarks[i]->Branch("CPAwrtPV", &Sexaquark[i].CPAwrtPV);
        Tree_Sexaquarks[i]->Branch("DCAwrtPV", &Sexaquark[i].DCAwrtPV);
        Tree_Sexaquarks[i]->Branch("Chi2ndf", &Sexaquark[i].Chi2ndf);
        /* -- True Information */
        if (Analysis::Settings::IsMC) {
            Tree_Sexaquarks[i]->Branch("IsSignal", &Sexaquark[i].IsSignal);
            Tree_Sexaquarks[i]->Branch("ReactionID", &Sexaquark[i].ReactionID);
            Tree_Sexaquarks[i]->Branch("IsHybrid", &Sexaquark[i].IsHybrid);
            Tree_Sexaquarks[i]->Branch("IsNonCombBkg", &Sexaquark[i].IsNonCombBkg);
            Tree_Sexaquarks[i]->Branch("AncestorIdx", &Sexaquark[i].AncestorIdx);
        }
        /* Shared with Channels "A"+"D"+"E" */
        Tree_Sexaquarks[i]->Branch("Idx_Lambda", &Sexaquark[i].Idx_Lambda);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Neg", &Sexaquark[i].Idx_Lambda_Neg);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Pos", &Sexaquark[i].Idx_Lambda_Pos);
        Tree_Sexaquarks[i]->Branch("Lambda_DecayLength", &Sexaquark[i].Lambda_DecayLength);
        Tree_Sexaquarks[i]->Branch("DCALaSV", &Sexaquark[i].DCALaSV);
        Tree_Sexaquarks[i]->Branch("DCALaNegSV", &Sexaquark[i].DCALaNegSV);
        Tree_Sexaquarks[i]->Branch("DCALaPosSV", &Sexaquark[i].DCALaPosSV);
        /* Shared with Channels "A"+"D"+"H" */
        Tree_Sexaquarks[i]->Branch("OpeningAngle", &Sexaquark[i].OpeningAngle);
        /* Specific to Channel "A" */
        Tree_Sexaquarks[i]->Branch("Idx_K0S", &Sexaquark[i].Idx_K0S);
        Tree_Sexaquarks[i]->Branch("Idx_K0S_Neg", &Sexaquark[i].Idx_K0S_Neg);
        Tree_Sexaquarks[i]->Branch("Idx_K0S_Pos", &Sexaquark[i].Idx_K0S_Pos);
        Tree_Sexaquarks[i]->Branch("K0S_DecayLength", &Sexaquark[i].K0S_DecayLength);
        Tree_Sexaquarks[i]->Branch("DCAK0SV", &Sexaquark[i].DCAK0SV);
        Tree_Sexaquarks[i]->Branch("DCAK0NegSV", &Sexaquark[i].DCAK0NegSV);
        Tree_Sexaquarks[i]->Branch("DCAK0PosSV", &Sexaquark[i].DCAK0PosSV);
        Tree_Sexaquarks[i]->Branch("DCAbtwV0s", &Sexaquark[i].DCAbtwV0s);
    }
}

/*
 *
 */
void Writer::InitSexaquarkBranches_TypeD() {
    //
    std::vector<SexaquarkD_tt> Sexaquark = {Sexaquark_ALPK, Sexaquark_LNK};
    std::vector<TTree*> Tree_Sexaquarks = {fTree_Sexaquarks_ALPK, fTree_Sexaquarks_LNK};
    for (Int_t i = 0; i < 2; i++) {
        /* Common Properties */
        Tree_Sexaquarks[i]->Branch("Px", &Sexaquark[i].Px);
        Tree_Sexaquarks[i]->Branch("Py", &Sexaquark[i].Py);
        Tree_Sexaquarks[i]->Branch("Pz", &Sexaquark[i].Pz);
        Tree_Sexaquarks[i]->Branch("E", &Sexaquark[i].E);
        Tree_Sexaquarks[i]->Branch("E_asDecay", &Sexaquark[i].E_asDecay);
        Tree_Sexaquarks[i]->Branch("Xv", &Sexaquark[i].Xv);
        Tree_Sexaquarks[i]->Branch("Yv", &Sexaquark[i].Yv);
        Tree_Sexaquarks[i]->Branch("Zv", &Sexaquark[i].Zv);
        Tree_Sexaquarks[i]->Branch("DistFromPV", &Sexaquark[i].DistFromPV);
        Tree_Sexaquarks[i]->Branch("CPAwrtPV", &Sexaquark[i].CPAwrtPV);
        Tree_Sexaquarks[i]->Branch("DCAwrtPV", &Sexaquark[i].DCAwrtPV);
        Tree_Sexaquarks[i]->Branch("Chi2ndf", &Sexaquark[i].Chi2ndf);
        /* -- True Information */
        if (Analysis::Settings::IsMC) {
            Tree_Sexaquarks[i]->Branch("IsSignal", &Sexaquark[i].IsSignal);
            Tree_Sexaquarks[i]->Branch("ReactionID", &Sexaquark[i].ReactionID);
            Tree_Sexaquarks[i]->Branch("IsHybrid", &Sexaquark[i].IsHybrid);
            Tree_Sexaquarks[i]->Branch("IsNonCombBkg", &Sexaquark[i].IsNonCombBkg);
            Tree_Sexaquarks[i]->Branch("AncestorIdx", &Sexaquark[i].AncestorIdx);
        }
        /* Shared with Channels "A"+"D"+"E" */
        Tree_Sexaquarks[i]->Branch("Idx_Lambda", &Sexaquark[i].Idx_Lambda);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Neg", &Sexaquark[i].Idx_Lambda_Neg);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Pos", &Sexaquark[i].Idx_Lambda_Pos);
        Tree_Sexaquarks[i]->Branch("Lambda_DecayLength", &Sexaquark[i].Lambda_DecayLength);
        Tree_Sexaquarks[i]->Branch("DCALaSV", &Sexaquark[i].DCALaSV);
        Tree_Sexaquarks[i]->Branch("DCALaNegSV", &Sexaquark[i].DCALaNegSV);
        Tree_Sexaquarks[i]->Branch("DCALaPosSV", &Sexaquark[i].DCALaPosSV);
        /* Shared with Channels "A"+"D"+"H" */
        Tree_Sexaquarks[i]->Branch("OpeningAngle", &Sexaquark[i].OpeningAngle);
        /* Specific to Channel "D" */
        Tree_Sexaquarks[i]->Branch("Idx_Kaon", &Sexaquark[i].Idx_Kaon);
        Tree_Sexaquarks[i]->Branch("DCAKaSV", &Sexaquark[i].DCAKaSV);
        Tree_Sexaquarks[i]->Branch("DCAKaLa", &Sexaquark[i].DCAKaLa);
        Tree_Sexaquarks[i]->Branch("DCALaNegKa", &Sexaquark[i].DCALaNegKa);
        Tree_Sexaquarks[i]->Branch("DCALaPosKa", &Sexaquark[i].DCALaPosKa);
    }
}

/*
 *
 */
void Writer::InitSexaquarkBranches_TypeE() {
    //
    std::vector<SexaquarkE_tt> Sexaquark = {Sexaquark_ALPKPP, Sexaquark_LNKPP};
    std::vector<TTree*> Tree_Sexaquarks = {fTree_Sexaquarks_ALPKPP, fTree_Sexaquarks_LNKPP};
    for (Int_t i = 0; i < 2; i++) {
        /* Common Properties */
        Tree_Sexaquarks[i]->Branch("Px", &Sexaquark[i].Px);
        Tree_Sexaquarks[i]->Branch("Py", &Sexaquark[i].Py);
        Tree_Sexaquarks[i]->Branch("Pz", &Sexaquark[i].Pz);
        Tree_Sexaquarks[i]->Branch("E", &Sexaquark[i].E);
        Tree_Sexaquarks[i]->Branch("E_asDecay", &Sexaquark[i].E_asDecay);
        Tree_Sexaquarks[i]->Branch("Xv", &Sexaquark[i].Xv);
        Tree_Sexaquarks[i]->Branch("Yv", &Sexaquark[i].Yv);
        Tree_Sexaquarks[i]->Branch("Zv", &Sexaquark[i].Zv);
        Tree_Sexaquarks[i]->Branch("DistFromPV", &Sexaquark[i].DistFromPV);
        Tree_Sexaquarks[i]->Branch("CPAwrtPV", &Sexaquark[i].CPAwrtPV);
        Tree_Sexaquarks[i]->Branch("DCAwrtPV", &Sexaquark[i].DCAwrtPV);
        Tree_Sexaquarks[i]->Branch("Chi2ndf", &Sexaquark[i].Chi2ndf);
        /* -- True Information */
        if (Analysis::Settings::IsMC) {
            Tree_Sexaquarks[i]->Branch("IsSignal", &Sexaquark[i].IsSignal);
            Tree_Sexaquarks[i]->Branch("ReactionID", &Sexaquark[i].ReactionID);
            Tree_Sexaquarks[i]->Branch("IsHybrid", &Sexaquark[i].IsHybrid);
            Tree_Sexaquarks[i]->Branch("IsNonCombBkg", &Sexaquark[i].IsNonCombBkg);
            Tree_Sexaquarks[i]->Branch("AncestorIdx", &Sexaquark[i].AncestorIdx);
        }
        /* Shared with Channels "A"+"D"+"E" */
        Tree_Sexaquarks[i]->Branch("Idx_Lambda", &Sexaquark[i].Idx_Lambda);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Neg", &Sexaquark[i].Idx_Lambda_Neg);
        Tree_Sexaquarks[i]->Branch("Idx_Lambda_Pos", &Sexaquark[i].Idx_Lambda_Pos);
        Tree_Sexaquarks[i]->Branch("Lambda_DecayLength", &Sexaquark[i].Lambda_DecayLength);
        Tree_Sexaquarks[i]->Branch("DCALaSV", &Sexaquark[i].DCALaSV);
        Tree_Sexaquarks[i]->Branch("DCALaNegSV", &Sexaquark[i].DCALaNegSV);
        Tree_Sexaquarks[i]->Branch("DCALaPosSV", &Sexaquark[i].DCALaPosSV);
        /* Specific to Channel "E" */
        Tree_Sexaquarks[i]->Branch("Idx_Kaon", &Sexaquark[i].Idx_Kaon);
        Tree_Sexaquarks[i]->Branch("Idx_PP", &Sexaquark[i].Idx_PP);
        Tree_Sexaquarks[i]->Branch("Idx_PiMinus", &Sexaquark[i].Idx_PiMinus);
        Tree_Sexaquarks[i]->Branch("Idx_PiPlus", &Sexaquark[i].Idx_PiPlus);
        Tree_Sexaquarks[i]->Branch("DCAKaSV", &Sexaquark[i].DCAKaSV);
        Tree_Sexaquarks[i]->Branch("DCAKaLa", &Sexaquark[i].DCAKaLa);
        Tree_Sexaquarks[i]->Branch("DCApmSV", &Sexaquark[i].DCApmSV);
        Tree_Sexaquarks[i]->Branch("DCAppSV", &Sexaquark[i].DCAppSV);
        Tree_Sexaquarks[i]->Branch("DCApmLa", &Sexaquark[i].DCApmLa);
        Tree_Sexaquarks[i]->Branch("DCAppLa", &Sexaquark[i].DCAppLa);
        Tree_Sexaquarks[i]->Branch("DCApmKa", &Sexaquark[i].DCApmKa);
        Tree_Sexaquarks[i]->Branch("DCAppKa", &Sexaquark[i].DCAppKa);
    }
}

/*
 *
 */
void Writer::InitKaonPairBranches() {
    //
    std::vector<KaonPair_tt> Sexaquark = {Sexaquark_PKPKX, Sexaquark_NKNKX};
    std::vector<TTree*> Tree_Sexaquarks = {fTree_Sexaquarks_PKPKX, fTree_Sexaquarks_NKNKX};
    for (Int_t i = 0; i < 2; i++) {
        /* Common Properties */
        Tree_Sexaquarks[i]->Branch("Px", &Sexaquark[i].Px);
        Tree_Sexaquarks[i]->Branch("Py", &Sexaquark[i].Py);
        Tree_Sexaquarks[i]->Branch("Pz", &Sexaquark[i].Pz);
        Tree_Sexaquarks[i]->Branch("E", &Sexaquark[i].E);
        Tree_Sexaquarks[i]->Branch("E_asDecay", &Sexaquark[i].E_asDecay);
        Tree_Sexaquarks[i]->Branch("Xv", &Sexaquark[i].Xv);
        Tree_Sexaquarks[i]->Branch("Yv", &Sexaquark[i].Yv);
        Tree_Sexaquarks[i]->Branch("Zv", &Sexaquark[i].Zv);
        Tree_Sexaquarks[i]->Branch("DistFromPV", &Sexaquark[i].DistFromPV);
        Tree_Sexaquarks[i]->Branch("CPAwrtPV", &Sexaquark[i].CPAwrtPV);
        Tree_Sexaquarks[i]->Branch("DCAwrtPV", &Sexaquark[i].DCAwrtPV);
        Tree_Sexaquarks[i]->Branch("Chi2ndf", &Sexaquark[i].Chi2ndf);
        /* -- True Information */
        if (Analysis::Settings::IsMC) {
            Tree_Sexaquarks[i]->Branch("IsSignal", &Sexaquark[i].IsSignal);
            Tree_Sexaquarks[i]->Branch("ReactionID", &Sexaquark[i].ReactionID);
            Tree_Sexaquarks[i]->Branch("IsHybrid", &Sexaquark[i].IsHybrid);
            Tree_Sexaquarks[i]->Branch("IsNonCombBkg", &Sexaquark[i].IsNonCombBkg);
            Tree_Sexaquarks[i]->Branch("AncestorIdx", &Sexaquark[i].AncestorIdx);
        }
        /* Shared with Channels "A"+"D"+"H" */
        Tree_Sexaquarks[i]->Branch("OpeningAngle", &Sexaquark[i].OpeningAngle);
        /* Specific to Channel "H" */
        Tree_Sexaquarks[i]->Branch("Idx_KaonA", &Sexaquark[i].Idx_KaonA);
        Tree_Sexaquarks[i]->Branch("Idx_KaonB", &Sexaquark[i].Idx_KaonB);
        Tree_Sexaquarks[i]->Branch("DCAbtwKK", &Sexaquark[i].DCAbtwKK);
        Tree_Sexaquarks[i]->Branch("DCAkaSV", &Sexaquark[i].DCAkaSV);
        Tree_Sexaquarks[i]->Branch("DCAkbSV", &Sexaquark[i].DCAkbSV);
    }
}

}  // namespace Tree2Sexaquark
