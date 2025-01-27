#include "Analysis/Manager.hxx"

#include <unordered_set>

// ROOT::Math libraries
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Math/KalmanFilter.hxx"

using namespace ROOT;

/*
 * Default constructor (with arguments)
 */
AnalysisManager::AnalysisManager(Settings_tt Settings)  //
    : Reader(),                                         //
      /*  */
      Settings(Settings),
      /*  */
      fPDG(),
      InputFile(nullptr),
      Event_UID(""),
      Event_Dir(nullptr),
      /*  */
      getMcEntry_fromMcIdx(),
      getPdgCode_fromMcIdx(),
      isMcIdxSignal(),
      isMcIdxSecondary(),
      getReactionID_fromMcIdx(),
      getMcIndices_fromReactionID(),
      getMotherMcIdx_fromMcIdx(),
      getAncestorMcIdx_fromMcIdx(),
      getNegDauMcIdx_fromMcIdx(),
      getPosDauMcIdx_fromMcIdx(),
      /*  */
      getTrackEntry_fromEsdIdx(),
      getMcIdx_fromEsdIdx(),
      /*  */
      esdIndicesOfAntiProtonTracks(),
      esdIndicesOfProtonTracks(),
      esdIndicesOfNegKaonTracks(),
      esdIndicesOfPosKaonTracks(),
      esdIndicesOfPiMinusTracks(),
      esdIndicesOfPiPlusTracks(),
      /*  */
      mcIndicesOfTrueV0s(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx(),
      /*  */
      getEsdIndices_fromReactionID() {}

/*
 *
 */
void AnalysisManager::Print() {
    //
    InfoF("IsMC           = %i", (Int_t)Settings.IsMC);
    InfoF("IsSignalMC     = %i", (Int_t)Settings.IsSignalMC);
    InfoF("InputFile      = %s", Settings.PathInputFile.c_str());
    InfoF("OutputFile     = %s", Settings.PathOutputFile.c_str());
    InfoF("LimitToNEvents = %lld", Settings.LimitToNEvents);
}

/*
 *
 */
Bool_t AnalysisManager::OpenInputFile() {
    //
    InputFile = std::unique_ptr<TFile>(TFile::Open((TString)Settings.PathInputFile, "READ"));
    if (!InputFile || InputFile->IsZombie()) {
        ErrorF("TFile %s couldn't be opened", Settings.PathInputFile.c_str());
        return kFALSE;
    }
    DebugF("TFile %s opened successfully", Settings.PathInputFile.c_str());

    TTree* EventsTree = FindTreeInFile("Events");
    if (!EventsTree) return kFALSE;

    SetEventsTree(EventsTree);
    ConnectEventBranches(Settings.IsMC);

    return kTRUE;
}

/*
 *
 */
Bool_t AnalysisManager::PrepareOutputFile() {
    //
    OutputFile = std::unique_ptr<TFile>(TFile::Open((TString)Settings.PathOutputFile, "RECREATE"));
    if (!OutputFile) {
        ErrorF("TFile %s couldn't be created", Settings.PathOutputFile.c_str());
        return kFALSE;
    }
    DebugF("TFile %s (re)created successfully", Settings.PathOutputFile.c_str());

    InitV0sTree();
    InitV0sBranches(Settings.IsMC);

    return kTRUE;
}

/*
 *
 */
void AnalysisManager::ProcessInjected() {
    //
    TTree* InjectedTree = FindTreeInEventDir("Injected");
    if (!InjectedTree) return;

    SetInjectedTree(InjectedTree);
    ConnectInjectedBranches();

    for (Long64_t sexa_entry = 0; sexa_entry < GetN_Injected(); sexa_entry++) {
        if (!ReadInjected(sexa_entry)) continue;
        // DebugF("%i, %f, %f, %f", Injected.ReactionID, Injected.Px, Injected.Py, Injected.Pz);
    }
}

/*
 *
 */
void AnalysisManager::ProcessMCParticles() {
    //
    TTree* MCTree = FindTreeInEventDir("MC");
    if (!MCTree) return;

    SetMCTree(MCTree);
    ConnectMCBranches();

    for (Long64_t mc_entry = 0; mc_entry < GetN_MCParticles(); mc_entry++) {
        if (!ReadMCParticle(mc_entry)) continue;
        /* Fill containers */
        getMcEntry_fromMcIdx[MC.Idx] = mc_entry;
        getPdgCode_fromMcIdx[MC.Idx] = MC.PdgCode;
        isMcIdxSignal[MC.Idx] = MC.Generator == 2 && MC.Idx_Ancestor == -1;
        isMcIdxSecondary[MC.Idx] = MC.IsSecFromMat || MC.IsSecFromWeak || isMcIdxSignal[MC.Idx];
        getReactionID_fromMcIdx[MC.Idx] = MC.ReactionID;
        getMcIndices_fromReactionID[MC.ReactionID].push_back(MC.Idx);
        getMotherMcIdx_fromMcIdx[MC.Idx] = MC.Idx_Mother;
        getAncestorMcIdx_fromMcIdx[MC.Idx] = MC.Idx_Ancestor;
        if (MC.Idx_Mother >= 0) {
            if (MC.PdgCode == -211 || MC.PdgCode == -2212) getNegDauMcIdx_fromMcIdx[MC.Idx_Mother] = MC.Idx;
            if (MC.PdgCode == 211 || MC.PdgCode == 2212) getPosDauMcIdx_fromMcIdx[MC.Idx_Mother] = MC.Idx;
        }
        if (MC.PdgCode == 310 || TMath::Abs(MC.PdgCode) == 3122) mcIndicesOfTrueV0s.push_back(MC.Idx);
        /*  */
        if (MC.Generator != 2 || MC.Idx_Ancestor != -1) continue;
        // InfoF("%i, %i, %i, %i, %i", MC.Idx, MC.Idx_Mother, MC.Idx_Ancestor, MC.PdgCode, MC.ReactionID);
    }
}

/*
 *
 */
void AnalysisManager::ProcessTracks() {
    //
    TTree* TracksTree = FindTreeInEventDir("Tracks");
    if (!TracksTree) return;

    SetTracksTree(TracksTree);
    ConnectTracksBranches();

    Int_t True_PdgCode;
    UInt_t MotherMcIdx;
    UInt_t Aux_McIdx;
    UInt_t ReactionID;

    for (Long64_t track_entry = 0; track_entry < GetN_Tracks(); track_entry++) {
        if (!ReadTrack(track_entry)) continue;
        /* Fill containers */
        getTrackEntry_fromEsdIdx[Track.Idx] = track_entry;
        getMcIdx_fromEsdIdx[Track.Idx] = Track.Idx_True;
        if (GetMotherMcIdx(Track.Idx_True, MotherMcIdx)) {
            if (GetNegDauMcIdx(MotherMcIdx, Aux_McIdx)) {
                if (Aux_McIdx == Track.Idx_True) getNegDauEsdIdx_fromMcIdx[MotherMcIdx] = Track.Idx;
            }
            if (GetPosDauMcIdx(MotherMcIdx, Aux_McIdx)) {
                if (Aux_McIdx == Track.Idx_True) getPosDauEsdIdx_fromMcIdx[MotherMcIdx] = Track.Idx;
            }
        }
        if (GetReactionID(Track.Idx_True, ReactionID)) {
            getEsdIndices_fromReactionID[ReactionID].push_back(Track.Idx);
        }
        /* PID */
        if (!GetPdgCode(Track.Idx_True, True_PdgCode)) continue;
        if (True_PdgCode == -2212) esdIndicesOfAntiProtonTracks.push_back(Track.Idx);
        if (True_PdgCode == 2212) esdIndicesOfProtonTracks.push_back(Track.Idx);
        if (True_PdgCode == -321) esdIndicesOfNegKaonTracks.push_back(Track.Idx);
        if (True_PdgCode == 321) esdIndicesOfPosKaonTracks.push_back(Track.Idx);
        if (True_PdgCode == -211) esdIndicesOfPiMinusTracks.push_back(Track.Idx);
        if (True_PdgCode == 211) esdIndicesOfPiPlusTracks.push_back(Track.Idx);
    }
}

/*         */
/**  V0s  **/
/*** === ***/

/*
 * Find all true (primary, secondary, signal) V0s for which both of their daughters were reconstructed and passed track selection
 */
void AnalysisManager::ProcessFindableV0s() {
    //
    Int_t V0_PdgCode;
    UInt_t Neg_McIdx, Pos_McIdx;
    Int_t Neg_PdgCode, Pos_PdgCode;
    UInt_t Neg_EsdIdx, Pos_EsdIdx;

    for (UInt_t& V0_McIdx : mcIndicesOfTrueV0s) {
        //
        if (!GetPdgCode(V0_McIdx, V0_PdgCode)) continue;
        if (!GetNegDauMcIdx(V0_McIdx, Neg_McIdx)) continue;
        if (!GetPosDauMcIdx(V0_McIdx, Pos_McIdx)) continue;
        if (!GetPdgCode(Neg_McIdx, Neg_PdgCode)) continue;
        if (!GetPdgCode(Pos_McIdx, Pos_PdgCode)) continue;
        if (!GetNegDauEsdIdx(V0_McIdx, Neg_EsdIdx)) continue;
        if (!GetPosDauEsdIdx(V0_McIdx, Pos_EsdIdx)) continue;
        //
        InfoF("mom: %u, %i, %i, neg: %u, %i, %u, pos: %u, %i, %u",  //
              V0_McIdx, V0_PdgCode, (Int_t)IsSignal(V0_McIdx),      //
              Neg_McIdx, Neg_PdgCode, Neg_EsdIdx,                   //
              Pos_McIdx, Pos_PdgCode, Pos_EsdIdx);
    }
}

/*
 * Find all V0s via Kalman Filter.
 * Note: if `pdgV0 == -1` (default value), it will store K0S candidates and secondary pion pairs.
 */
void AnalysisManager::KalmanV0Finder(Int_t pdgNegDaughter, Int_t pdgPosDaughter, Int_t pdgV0) {

    KFParticle::SetField(Event.MagneticField);

    Track_tt TrackNeg, TrackPos;
    Int_t mcIdxNeg, mcIdxPos;
    Int_t mcIdxV0;

    /* Declare KFParticle objects */

    KFParticle kfDaughterNeg, kfDaughterPos;
    KFParticle kfTransportedNeg, kfTransportedPos;

    /* Declare 4-momentum vectors */

    Math::PxPyPzEVector lvTrackNeg, lvTrackPos;
    Math::PxPyPzEVector lvV0;

    /* Information from MC */

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;
    Int_t reaction_id;
    Bool_t is_hybrid;

    /* Choose tracks species to loop over */

    std::vector<UInt_t> esdIndicesNegTracks;
    std::vector<UInt_t> esdIndicesPosTracks;
    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    } else if (pdgV0 == 3122) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfProtonTracks;
    } else {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    }

    /* Loop over all possible pairs of tracks */

    for (UInt_t& esdIdxNeg : esdIndicesNegTracks) {
        for (UInt_t& esdIdxPos : esdIndicesPosTracks) {

            /* Sanity check: prevent tracks from being repeated */

            if (esdIdxNeg == esdIdxPos) continue;

            /* Get tracks */

            GetTrack(esdIdxNeg, TrackNeg);
            GetTrack(esdIdxPos, TrackPos);

            /* Kalman Filter */

            kfDaughterNeg = CreateKFParticle(TrackNeg, fPDG.GetParticle(pdgNegDaughter)->Mass());
            kfDaughterPos = CreateKFParticle(TrackPos, fPDG.GetParticle(pdgPosDaughter)->Mass());

            KFParticle kfV0;
            kfV0.AddDaughter(kfDaughterNeg);
            kfV0.AddDaughter(kfDaughterPos);

            /* Transport V0 and daughters */

            kfV0.TransportToDecayVertex();
            kfTransportedNeg = TransportKFParticle(kfDaughterNeg, kfDaughterPos, fPDG.GetParticle(pdgNegDaughter)->Mass(),  //
                                                   (Int_t)TrackNeg.Charge);
            kfTransportedPos = TransportKFParticle(kfDaughterPos, kfDaughterNeg, fPDG.GetParticle(pdgPosDaughter)->Mass(),  //
                                                   (Int_t)TrackPos.Charge);

            /* Reconstruct V0 */

            lvTrackNeg = Math::PxPyPzMVector(kfDaughterNeg.Px(), kfDaughterNeg.Py(), kfDaughterNeg.Pz(),  //
                                             fPDG.GetParticle(pdgNegDaughter)->Mass());
            lvTrackPos = Math::PxPyPzMVector(kfDaughterPos.Px(), kfDaughterPos.Py(), kfDaughterPos.Pz(),  //
                                             fPDG.GetParticle(pdgPosDaughter)->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            /* Optimization: if `pdgV0 == -1`, use the same pi+pi- loop to process both K0S and pi+pi- coming from a signal reaction */

            std::vector<Int_t> pdgV0s = {pdgV0};
            if (pdgV0 == -1) pdgV0s = {422, 310};

            for (Int_t& auxPdgV0 : pdgV0s) {

                /* Collect true information */

                UInt_t mc_neg, mc_pos;
                UInt_t mc_neg_mc_mother, mc_pos_mc_mother;
                Int_t mc_pdg_mother;
                if (!GetMcIdx(esdIdxNeg, mc_neg) || !GetMcIdx(esdIdxPos, mc_pos)) continue;
                if (!GetMotherMcIdx(mc_neg, mc_neg_mc_mother) || !GetMotherMcIdx(mc_pos, mc_pos_mc_mother)) continue;

                /* DEBUG: select only true V0s */

                if (mc_neg_mc_mother != mc_pos_mc_mother) continue;
                if (!GetPdgCode(mc_neg_mc_mother, mc_pdg_mother)) continue;
                if (mc_pdg_mother != pdgV0) continue;

                InfoF("%u, %u, %u, %u, %u, %i, %f", esdIdxNeg, esdIdxPos, mc_neg, mc_pos, mc_neg_mc_mother, mc_pdg_mother, lvV0.M());
                FillV0(1, esdIdxNeg, esdIdxPos, auxPdgV0, lvV0, kfV0, lvTrackNeg, lvTrackPos, Settings.IsMC);
                /*
                mcIdxPos = getMcIdx_fromEsdIdx[esdIdxPos];
                mcIdxV0 = doesMcIdxHaveMother[mcIdxNeg] && doesMcIdxHaveMother[mcIdxPos] &&
                                    getMotherMcIdx_fromMcIdx[mcIdxNeg] == getMotherMcIdx_fromMcIdx[mcIdxPos]
                                ? getMotherMcIdx_fromMcIdx[mcIdxNeg]
                                : -1;
                mc_pdg_code = mcIdxV0 >= 0 ? getPdgCode_fromMcIdx[mcIdxV0] : -1;
                */
                /** Note: different treatment for pion pairs, since they don't have a true MC V0 **/
                /*
                if (auxPdgV0 == 422) {
                    is_true = !doesMcIdxHaveMother[mcIdxNeg] && !doesMcIdxHaveMother[mcIdxPos] &&  //
                                getPdgCode_fromMcIdx[mcIdxNeg] == pdgNegDaughter && getPdgCode_fromMcIdx[mcIdxPos] == pdgPosDaughter;
                    is_secondary = isMcIdxSecondary[mcIdxNeg] && isMcIdxSecondary[mcIdxPos] && is_true;
                    is_signal = isMcIdxSignal[mcIdxNeg] && isMcIdxSignal[mcIdxPos] &&                      //
                                getReactionID_fromMcIdx[mcIdxNeg] == getReactionID_fromMcIdx[mcIdxPos] &&  //
                                is_secondary;
                    reaction_id = is_signal ? getReactionID_fromMcIdx[mcIdxNeg] : -1;
                } else {
                    is_true = doesMcIdxHaveMother[mcIdxNeg] && doesMcIdxHaveMother[mcIdxPos] &&
                                getMotherMcIdx_fromMcIdx[mcIdxNeg] == getMotherMcIdx_fromMcIdx[mcIdxPos] && //
                                getPdgCode_fromMcIdx[mcIdxNeg] == pdgNegDaughter && getPdgCode_fromMcIdx[mcIdxPos] == pdgPosDaughter &&
                // mc_pdg_code == auxPdgV0; is_secondary = isMcIdxSecondary[mcIdxV0] && is_true; is_signal = isMcIdxSignal[mcIdxV0] &&
                is_secondary; reaction_id = getReactionID_fromMcIdx[mcIdxV0];
                            }
                            is_hybrid = (isMcIdxSignal[mcIdxNeg] || isMcIdxSignal[mcIdxPos]) && !is_signal;
                */
                /* Apply cuts and store V0 */
                /*
                if (PassesV0CutsAs(auxPdgV0, kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos, is_true, is_secondary, is_signal)) {
                    StoreV0As(auxPdgV0, kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos, esdIdxNeg, esdIdxPos, mcIdxV0, mc_pdg_code,
                              is_true, is_secondary, is_signal, reaction_id, is_hybrid);
                }
                */
            }
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/*
 *
 */
void AnalysisManager::ProcessFindableSexaquarks() {
    //
    Int_t PdgCode_StruckNucleon = 2112;
    std::vector<Int_t> PdgCodes_FinalStateProducts = {-2212, 211, -211, 211};

    Long64_t track_entry;
    Int_t pdg_code;

    for (Long64_t sexa_entry = 0; sexa_entry < GetN_Injected(); sexa_entry++) {
        if (!ReadInjected(sexa_entry)) continue;

        /* All final state products should have been reconstructed and passed track selection */

        std::unordered_multiset<Int_t> pdg_reconstructed_particles;
        for (UInt_t& EsdIdx : getEsdIndices_fromReactionID[Injected.ReactionID]) {
            pdg_reconstructed_particles.insert(getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[EsdIdx]]);
        }

        std::unordered_multiset<Int_t> expected_pdg_fs_particles(PdgCodes_FinalStateProducts.begin(), PdgCodes_FinalStateProducts.end());

        if (expected_pdg_fs_particles != pdg_reconstructed_particles) continue;

        /** Loop, to get the reconstructed final-state particles info **/

        for (UInt_t& EsdIdx : getEsdIndices_fromReactionID[Injected.ReactionID]) {
            if (!GetTrackEntry(EsdIdx, track_entry)) continue;
            ReadTrack(track_entry);
            if (!GetPdgCode(Track.Idx_True, pdg_code)) continue;
            InfoF("%i, %u, %u, %i, %f, %f, %f", Injected.ReactionID, EsdIdx, Track.Idx_True, pdg_code, Track.Px, Track.Py, Track.Pz);
        }
    }
}

/*
 *
 */
void AnalysisManager::KalmanSexaquarkFinder(Int_t pdgStruckNucleon, std::vector<Int_t> pdgReactionProducts) {
    //
    if (pdgReactionProducts.size() > 2) {
        KalmanSexaquarkFinder_TypeDE(pdgReactionProducts);
    } else {
        if (TMath::Abs(pdgStruckNucleon) == 2112) {
            KalmanSexaquarkFinder_TypeA(pdgReactionProducts);
        } else {
            KalmanSexaquarkFinder_TypeH(pdgReactionProducts);
        }
    }
}

/*
 *
 */
void AnalysisManager::KalmanSexaquarkFinder_TypeA(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*
 *
 */
void AnalysisManager::KalmanSexaquarkFinder_TypeDE(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*
 *
 */
void AnalysisManager::KalmanSexaquarkFinder_TypeH(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*               */
/**  Utilities  **/
/*** ========= ***/

/*
 *
 */
void AnalysisManager::CleanContainers() {
    //
    getMcEntry_fromMcIdx.clear();
    getPdgCode_fromMcIdx.clear();
    isMcIdxSignal.clear();
    isMcIdxSecondary.clear();
    getReactionID_fromMcIdx.clear();
    getMcIndices_fromReactionID.clear();
    getMotherMcIdx_fromMcIdx.clear();
    getAncestorMcIdx_fromMcIdx.clear();
    getNegDauMcIdx_fromMcIdx.clear();
    getPosDauMcIdx_fromMcIdx.clear();
    //
    getTrackEntry_fromEsdIdx.clear();
    getMcIdx_fromEsdIdx.clear();
    //
    esdIndicesOfAntiProtonTracks.clear();
    esdIndicesOfProtonTracks.clear();
    esdIndicesOfNegKaonTracks.clear();
    esdIndicesOfPosKaonTracks.clear();
    esdIndicesOfPiMinusTracks.clear();
    esdIndicesOfPiPlusTracks.clear();
    //
    mcIndicesOfTrueV0s.clear();
    getNegDauEsdIdx_fromMcIdx.clear();
    getPosDauEsdIdx_fromMcIdx.clear();
    //
    getEsdIndices_fromReactionID.clear();
}

/*
 *
 */
void AnalysisManager::EndOfEvent() {
    //
    DisconnectInjectedBranches();
    DisconnectMCBranches();
    DisconnectTracksBranches();
    CleanContainers();
}

/*
 *
 */
void AnalysisManager::EndOfAnalysis() {
    //
    DisconnectEventBranches();
    //
    if (OutputFile) {
        OutputFile->cd();
        WriteV0sTree();
        OutputFile->Write();
    }
}
