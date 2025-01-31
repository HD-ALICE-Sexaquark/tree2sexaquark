#include "Analysis/Manager.hxx"

#include <unordered_set>

#include "Math/Vector4D.h"

#include "KFParticle.h"

#include "Math/KalmanFilter.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

/*
 *
 */
Bool_t Manager::OpenInputFile() {
    //
    InputFile = std::unique_ptr<TFile>(TFile::Open((TString)Settings::PathInputFile, "READ"));
    if (!InputFile || InputFile->IsZombie()) {
        ErrorF("TFile %s couldn't be opened", Settings::PathInputFile.c_str());
        return kFALSE;
    }
    DebugF("TFile %s opened successfully", Settings::PathInputFile.c_str());

    TTree* EventsTree = FindTreeInFile("Events");
    if (!EventsTree) return kFALSE;

    SetEventsTree(EventsTree);
    ConnectEventBranches();

    Inspector.SetLambdaDefaultCuts();
    Inspector.SetKaonZeroDefaultCuts();
    Inspector.SetPionPairDefaultCuts();

    Inspector.SetSexaquarkDefaultCuts_ChannelA();
    Inspector.SetSexaquarkDefaultCuts_ChannelD();
    Inspector.SetSexaquarkDefaultCuts_ChannelE();
    Inspector.SetKaonPairDefaultCuts();

    Inspector.PrintAllCuts();

    return kTRUE;
}

/*
 *
 */
Bool_t Manager::PrepareOutputFile() {
    //
    OutputFile = std::unique_ptr<TFile>(TFile::Open((TString)Settings::PathOutputFile, "RECREATE"));
    if (!OutputFile) {
        ErrorF("TFile %s couldn't be created", Settings::PathOutputFile.c_str());
        return kFALSE;
    }
    DebugF("TFile %s (re)created successfully", Settings::PathOutputFile.c_str());

    InitV0sTree();
    InitV0sBranches();

    return kTRUE;
}

/*
 *
 */
Bool_t Manager::GetEvent(Long64_t evt_idx) {
    //
    if (!ReadEvent(evt_idx)) {
        DebugF("Event # %lld couldn't be read, moving on...", evt_idx);
        return kFALSE;
    }
    if (Settings::IsMC) {
        if (Settings::IsSignalMC)
            Event_UID = TString::Format("A18_%u_%04u_%03u", Event.RunNumber, Event.DirNumber, Event.EventNumber);
        else
            Event_UID = TString::Format("MC_%6u_%04u_%03u", Event.RunNumber, Event.DirNumber, Event.EventNumber);
    } else {
        Event_UID = TString::Format("DATA_%6u_%03u_%u_%03u", Event.RunNumber, Event.DirNumber, Event.DirNumberB, Event.EventNumber);
    }
    InfoF("Processing Event # %lld (UID = %s)", evt_idx, Event_UID.Data());
    // InfoF(">> Centrality = %f, PV = (%f, %f, %f), B = %f", Event.Centrality, Event.PV_Xv, Event.PV_Yv, Event.PV_Zv, Event.MagneticField);

    Event_Dir = std::unique_ptr<TDirectoryFile>(InputFile->Get<TDirectoryFile>(Event_UID));
    if (!Event_Dir) {
        DebugF("TDirectoryFile %s couldn't be found, moving on...", Event_UID.Data());
        return kFALSE;
    }

    /* Set Event Properties */
    /* -- Magnetic Field */
    KFParticle::SetField(Event.MagneticField);
    /* -- Primary Vertex */
    Float_t XYZ[3] = {Event.PV_Xv, Event.PV_Yv, Event.PV_Zv};
    Float_t CovMatrix[6] = {Event.PV_CovMatrix[0], Event.PV_CovMatrix[1], Event.PV_CovMatrix[2],
                            Event.PV_CovMatrix[3], Event.PV_CovMatrix[4], Event.PV_CovMatrix[5]};
    kfPrimaryVertex = Math::CreateKFVertex(XYZ, CovMatrix);

    return kTRUE;
}

/*
 *
 */
void Manager::ProcessInjected() {
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
void Manager::ProcessMCParticles() {
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
        if (isMcIdxSignal[MC.Idx]) {
            getReactionID_fromMcIdx[MC.Idx] = MC.ReactionID;
            getMcIndices_fromReactionID[MC.ReactionID].push_back(MC.Idx);
        } else {
            getReactionID_fromMcIdx[MC.Idx] = -1;
        }
        getMotherMcIdx_fromMcIdx[MC.Idx] = MC.Idx_Mother;
        getAncestorMcIdx_fromMcIdx[MC.Idx] = MC.Idx_Ancestor;
        if (MC.Idx_Mother >= 0) {
            if (MC.PdgCode == -211 || MC.PdgCode == -2212) getNegDauMcIdx_fromMcIdx[MC.Idx_Mother] = MC.Idx;
            if (MC.PdgCode == 211 || MC.PdgCode == 2212) getPosDauMcIdx_fromMcIdx[MC.Idx_Mother] = MC.Idx;
        }
        if (MC.PdgCode == 310 || TMath::Abs(MC.PdgCode) == 3122) mcIndicesOfTrueV0s.push_back(MC.Idx);
        /*  */
        // if (isMcIdxSignal[MC.Idx]) continue;
        // InfoF("entry=%lld, mcIdx=%u, pdg=%i, primary=%i, secfrommat=%i, secfromweak=%i, mother=%i, ancestor=%i",  //
        //   mc_entry, MC.Idx, MC.PdgCode, MC.IsPrimary, MC.IsSecFromMat, MC.IsSecFromWeak, MC.Idx_Mother, MC.Idx_Ancestor);
    }
}

/*
 *
 */
void Manager::ProcessTracks() {
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
void Manager::ProcessFindableV0s() {
    //
    Int_t V0_PdgCode, Neg_PdgCode, Pos_PdgCode;
    UInt_t Neg_McIdx, Pos_McIdx;
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
void Manager::KalmanV0Finder(Int_t pdgNeg, Int_t pdgPos, Int_t pdgV0) {

    Track_tt TrackNeg, TrackPos;
    UInt_t mcIdxNeg, mcIdxPos;
    Int_t mcIdxV0;

    /* Declare KFParticle objects */

    KFParticle kfNeg, kfPos;
    KFParticle kfTransportedNeg, kfTransportedPos;

    /* Declare 4-momentum vectors */

    ROOT::Math::PxPyPzEVector lvNeg, lvPos;
    ROOT::Math::PxPyPzEVector lvV0;

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

            kfNeg = Math::CreateKFParticle(TrackNeg, fPDG.GetParticle(pdgNeg)->Mass());
            kfPos = Math::CreateKFParticle(TrackPos, fPDG.GetParticle(pdgPos)->Mass());

            KFParticle kfV0;
            kfV0.AddDaughter(kfNeg);
            kfV0.AddDaughter(kfPos);

            /* Transport V0 and daughters */

            kfV0.TransportToDecayVertex();
            kfTransportedNeg = Math::TransportKFParticle(kfNeg, kfPos, fPDG.GetParticle(pdgNeg)->Mass(),  //
                                                         (Int_t)TrackNeg.Charge);
            kfTransportedPos = Math::TransportKFParticle(kfPos, kfNeg, fPDG.GetParticle(pdgPos)->Mass(),  //
                                                         (Int_t)TrackPos.Charge);

            /* Reconstruct V0 */

            lvNeg = ROOT::Math::PxPyPzMVector(kfNeg.Px(), kfNeg.Py(), kfNeg.Pz(), fPDG.GetParticle(pdgNeg)->Mass());
            lvPos = ROOT::Math::PxPyPzMVector(kfPos.Px(), kfPos.Py(), kfPos.Pz(), fPDG.GetParticle(pdgPos)->Mass());
            lvV0 = lvNeg + lvPos;

            /* Optimization: if `pdgV0 == -1`, use the same pi+pi- loop to process both K0S and pi+pi- coming from a signal reaction */

            std::vector<Int_t> pdgV0s = {pdgV0};
            if (pdgV0 == -1) pdgV0s = {310, 422};

            for (Int_t& auxPdgV0 : pdgV0s) {

                /* Collect true information */

                UInt_t mc_neg_mc_mother, mc_pos_mc_mother;
                Int_t mc_pdg_mother;

                if (!GetMcIdx(esdIdxNeg, mcIdxNeg) || !GetMcIdx(esdIdxPos, mcIdxPos)) continue;
                if (!GetMotherMcIdx(mcIdxNeg, mc_neg_mc_mother) || !GetMotherMcIdx(mcIdxPos, mc_pos_mc_mother)) continue;

                /* DEBUG: select only true V0s */

                if (mc_neg_mc_mother != mc_pos_mc_mother) continue;
                if (!GetPdgCode(mc_neg_mc_mother, mc_pdg_mother)) continue;
                if (mc_pdg_mother != auxPdgV0) continue;

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

                Particle::V0 ThisV0(auxPdgV0, esdIdxNeg, esdIdxPos,  //
                                    lvV0, lvNeg, lvPos,              //
                                    kfV0, kfTransportedNeg, kfTransportedPos, kfPrimaryVertex);

                if (!Inspector.Approve(ThisV0)) continue;

                FillV0(1, ThisV0);
                InfoF("%u, %u, %u, %u, %u, %i, %f", esdIdxNeg, esdIdxPos, mcIdxNeg, mcIdxPos, mc_neg_mc_mother, mc_pdg_mother, lvV0.M());
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
void Manager::ProcessFindableSexaquarks() {
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
void Manager::KalmanSexaquarkFinder(Int_t pdgStruckNucleon, std::vector<Int_t> pdgReactionProducts) {
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
void Manager::KalmanSexaquarkFinder_TypeA(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*
 *
 */
void Manager::KalmanSexaquarkFinder_TypeDE(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*
 *
 */
void Manager::KalmanSexaquarkFinder_TypeH(std::vector<Int_t> pdgReactionProducts) {
    //
}

/*               */
/**  Utilities  **/
/*** ========= ***/

/*
 *
 */
void Manager::CleanContainers() {
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
void Manager::EndOfEvent() {
    //
    DisconnectInjectedBranches();
    DisconnectMCBranches();
    DisconnectTracksBranches();
    CleanContainers();
}

/*
 *
 */
void Manager::EndOfAnalysis() {
    //
    DisconnectEventBranches();
    //
    if (OutputFile) {
        OutputFile->cd();
        WriteV0sTree();
        OutputFile->Write();
    }
}

}  // namespace Analysis
}  // namespace Tree2Sexaquark
