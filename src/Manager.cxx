#include "Analysis/Manager.hxx"

#include <unordered_set>

#include "Math/Vector4D.h"

#include "KFParticle.h"

#include "Cuts/Default.hxx"
#include "Math/KalmanFilter.hxx"

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
            Event_UID = TString::Format("A18_%u_%04u_%03u", Event.RunNumber, Event.DirNumber, Event.EventNumber);  // PENDING
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
        InfoF("%i, %f, %f, %f", Injected.ReactionID, Injected.Px, Injected.Py, Injected.Pz);
    }  // end of loop over injected reactions
}

/*
 *
 */
void Manager::ProcessMCParticles() {
    //
    TTree* InputTree_MC = FindTreeInEventDir("MC");
    if (!InputTree_MC) return;

    SetMCTree(InputTree_MC);
    ConnectMCBranches();

    fTree_MC->SetBranchStatus("*", 0);
    fTree_MC->SetBranchStatus("Idx", 1);
    fTree_MC->SetBranchStatus("PdgCode", 1);
    fTree_MC->SetBranchStatus("Idx_Mother", 1);
    fTree_MC->SetBranchStatus("Idx_Ancestor", 1);
    fTree_MC->SetBranchStatus("Generator", 1);
    fTree_MC->SetBranchStatus("ReactionID", 1);

    for (Long64_t MC_Entry = 0; MC_Entry < GetN_MCParticles(); MC_Entry++) {
        if (!ReadMCParticle(MC_Entry)) continue;
        /* Fill containers */
        getMcIdx_FromMcEntry[MC_Entry] = MC.Idx;
        getMcEntry_fromMcIdx[MC.Idx] = MC_Entry;
        if (MC.PdgCode == 310 || TMath::Abs(MC.PdgCode) == 3122) mcEntriesOfTrueV0s.push_back(MC_Entry);
        if (MC.Idx_Mother >= 0) {
            if (MC.PdgCode == -211 || MC.PdgCode == -2212) getMcIdxNeg_fromMcIdx[MC.Idx_Mother] = MC.Idx;
            if (MC.PdgCode == 211 || MC.PdgCode == 2212) getMcIdxPos_fromMcIdx[MC.Idx_Mother] = MC.Idx;
        }
        if (MC_IsFinalStateSignal()) getMcEntries_fromReactionID[(UInt_t)MC.ReactionID].push_back(MC_Entry);
    }  // end of loop over MC particles
    InfoF("N Found True V0s: %lu", mcEntriesOfTrueV0s.size());
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

    fTree_Tracks->SetBranchStatus("*", 0);
    fTree_Tracks->SetBranchStatus("Idx", 1);
    fTree_Tracks->SetBranchStatus("Charge", 1);
    fTree_Tracks->SetBranchStatus("NSigmaProton", 1);
    fTree_Tracks->SetBranchStatus("NSigmaKaon", 1);
    fTree_Tracks->SetBranchStatus("NSigmaPion", 1);
    fTree_Tracks->SetBranchStatus("Idx_True", 1);

    for (Long64_t Track_Entry = 0; Track_Entry < GetN_Tracks(); Track_Entry++) {
        if (!ReadTrack(Track_Entry)) continue;
        /* Fill containers */
        getTrackEntry_fromEsdIdx[Track.Idx] = Track_Entry;
        getEsdIdx_fromTrackEntry[Track_Entry] = Track.Idx;
        if (Settings::IsMC) {
            getMcIdx_fromEsdIdx[Track.Idx] = Track.Idx_True;
            getEsdIdx_fromMcIdx[Track.Idx_True] = Track.Idx;
        }
        /* PID */
        if (TMath::Abs(Track.NSigmaProton) < Default::Proton::MaxNSigma) {
            if (Track.Charge < 0)
                TrackEntries_AntiProton.push_back(Track_Entry);
            else
                TrackEntries_Proton.push_back(Track_Entry);
        }
        if (TMath::Abs(Track.NSigmaKaon) < Default::Kaon::MaxNSigma) {
            if (Track.Charge < 0)
                TrackEntries_NegKaon.push_back(Track_Entry);
            else
                TrackEntries_PosKaon.push_back(Track_Entry);
        }
        if (TMath::Abs(Track.NSigmaPion) < Default::Pion::MaxNSigma) {
            if (Track.Charge < 0)
                TrackEntries_PiMinus.push_back(Track_Entry);
            else
                TrackEntries_PiPlus.push_back(Track_Entry);
        }
    }  // end of loop over tracks
    InfoF("N Found AntiProtons: %lu", TrackEntries_AntiProton.size());
    InfoF("N Found Protons: %lu", TrackEntries_Proton.size());
    InfoF("N Found Negative Kaons: %lu", TrackEntries_NegKaon.size());
    InfoF("N Found Positive Kaons: %lu", TrackEntries_PosKaon.size());
    InfoF("N Found Negative Pions: %lu", TrackEntries_PiMinus.size());
    InfoF("N Found Positive Pions: %lu", TrackEntries_PiPlus.size());
}

/*         */
/**  V0s  **/
/*** === ***/

/*
 * Find all true (primary, secondary, signal) V0s for which both of their daughters were reconstructed and passed track selection
 */
void Manager::ProcessFindableV0s() {
    //
    fTree_MC->SetBranchStatus("*", 0);
    fTree_MC->SetBranchStatus("Idx", 1);
    fTree_MC->SetBranchStatus("PdgCode", 1);
    fTree_MC->SetBranchStatus("Generator", 1);
    fTree_MC->SetBranchStatus("Idx_Ancestor", 1);

    fTree_Tracks->SetBranchStatus("*", 0);
    fTree_Tracks->SetBranchStatus("Idx", 1);

    Bool_t V0_IsSignal;
    UInt_t V0_McIdx;
    Int_t V0_PdgCode;

    Long64_t Neg_McEntry, Pos_McEntry;
    UInt_t Neg_McIdx, Pos_McIdx;
    Int_t Neg_PdgCode, Pos_PdgCode;

    Long64_t Neg_TrackEntry, Pos_TrackEntry;
    UInt_t Neg_EsdIdx, Pos_EsdIdx;

    for (Long64_t& V0_McEntry : mcEntriesOfTrueV0s) {
        if (!ReadMCParticle(V0_McEntry)) continue;
        V0_McIdx = MC.Idx;
        V0_PdgCode = MC.PdgCode;
        V0_IsSignal = MC_IsSignal();
        //
        Neg_McEntry = GetMcEntryNeg(V0_McEntry);
        if (!ReadMCParticle(Neg_McEntry)) continue;
        Neg_McIdx = MC.Idx;
        Neg_PdgCode = MC.PdgCode;
        Neg_TrackEntry = GetTrackEntryFromMcEntry(Neg_McEntry);
        if (!ReadTrack(Neg_TrackEntry)) continue;
        Neg_EsdIdx = Track.Idx;
        //
        Pos_McEntry = GetMcEntryPos(V0_McEntry);
        if (!ReadMCParticle(Pos_McEntry)) continue;
        Pos_McIdx = MC.Idx;
        Pos_PdgCode = MC.PdgCode;
        Pos_TrackEntry = GetTrackEntryFromMcEntry(Pos_McEntry);
        if (!ReadTrack(Pos_TrackEntry)) continue;
        Pos_EsdIdx = Track.Idx;
        //
        InfoF("v0: mcidx=%u, pdg=%i, signal=%i, neg: mcidx=%u, pdg=%i, esdidx=%u, pos: mcidx=%u, pdg=%i, esdidx=%u",  //
              V0_McIdx, V0_PdgCode, V0_IsSignal,                                                                      //
              Neg_McIdx, Neg_PdgCode, Neg_EsdIdx,                                                                     //
              Pos_McIdx, Pos_PdgCode, Pos_EsdIdx);
    }  // end of loop over true V0s
}

/*
 * Find all V0s via Kalman Filter.
 */
void Manager::KalmanV0Finder(Int_t pdgNeg, Int_t pdgPos, Int_t pdgV0) {
    //
    Track_tt TrackNeg, TrackPos;
    fTree_Tracks->SetBranchStatus("*", 0);
    fTree_Tracks->SetBranchStatus("Idx", 1);
    fTree_Tracks->SetBranchStatus("Px", 1);
    fTree_Tracks->SetBranchStatus("Py", 1);
    fTree_Tracks->SetBranchStatus("Pz", 1);
    fTree_Tracks->SetBranchStatus("X", 1);
    fTree_Tracks->SetBranchStatus("Y", 1);
    fTree_Tracks->SetBranchStatus("Z", 1);
    fTree_Tracks->SetBranchStatus("Charge", 1);
    fTree_Tracks->SetBranchStatus("Alpha", 1);
    fTree_Tracks->SetBranchStatus("Snp", 1);
    fTree_Tracks->SetBranchStatus("Tgl", 1);
    fTree_Tracks->SetBranchStatus("Signed1Pt", 1);
    fTree_Tracks->SetBranchStatus("CovMatrix", 1);
    /* Declare KFParticle objects */
    KFParticle kfNeg, kfPos;
    KFParticle kfTransportedNeg, kfTransportedPos;
    /* Declare 4-momentum vectors */
    ROOT::Math::PxPyPzEVector lvNeg, lvPos;
    ROOT::Math::PxPyPzEVector lvV0;
    /* Declare the V0 object */
    Particle::V0 newV0;
    newV0.SetPrimaryVertex(kfPrimaryVertex);
    /* Information from MC */
    Bool_t is_true;
    Int_t mc_idx_v0;
    Int_t mc_pdg_code;
    Bool_t is_secondary;
    Bool_t is_signal;
    Int_t reaction_id;
    Bool_t is_hybrid;
    /* Choose tracks species to loop over */
    std::vector<Long64_t> TrackEntries_Neg;
    std::vector<Long64_t> TrackEntries_Pos;
    if (pdgV0 == -3122) {
        TrackEntries_Neg = TrackEntries_AntiProton;
        TrackEntries_Pos = TrackEntries_PiPlus;
    } else if (pdgV0 == 3122) {
        TrackEntries_Neg = TrackEntries_PiMinus;
        TrackEntries_Pos = TrackEntries_Proton;
    } else {
        /* For KaonsZeroShort and PionPairs */
        TrackEntries_Neg = TrackEntries_PiMinus;
        TrackEntries_Pos = TrackEntries_PiPlus;
    }
    /* Loop over all possible pairs of tracks */
    for (Long64_t& Neg_Entry : TrackEntries_Neg) {
        for (Long64_t& Pos_Entry : TrackEntries_Pos) {
            /* Sanity check: prevent tracks from being repeated */
            if (Neg_Entry == Pos_Entry) continue;
            /* Get tracks */
            if (!CopyTrack(Neg_Entry, TrackNeg)) continue;
            if (!CopyTrack(Pos_Entry, TrackPos)) continue;
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
            lvNeg = ROOT::Math::PxPyPzMVector(kfTransportedNeg.Px(), kfTransportedNeg.Py(), kfTransportedNeg.Pz(), fPDG.GetParticle(pdgNeg)->Mass());
            lvPos = ROOT::Math::PxPyPzMVector(kfTransportedPos.Px(), kfTransportedPos.Py(), kfTransportedPos.Pz(), fPDG.GetParticle(pdgPos)->Mass());
            lvV0 = lvNeg + lvPos;
            /* Prepare V0 object */
            newV0.SetKinematics(lvV0, lvNeg, lvPos);
            newV0.SetGeometry(kfV0, kfTransportedNeg, kfTransportedPos);
            newV0.SetV0Info(pdgV0, TrackNeg.Idx, TrackPos.Idx);
            // if (Settings::IsMC) {
            // CollectTrueV0(pdgV0, TrackNeg.Idx, TrackPos.Idx, is_true, mc_idx_v0, mc_pdg_code, is_secondary, is_signal, reaction_id,
            //   is_hybrid);
            // newV0.SetTrueInfo(is_true, mc_idx_v0, mc_pdg_code, is_secondary, is_signal, reaction_id, is_hybrid);
            // }
            /* Apply cuts and store V0 */
            if (!Inspector.Approve(newV0)) continue;
            if (pdgV0 == -3122)
                AntiLambdas.push_back(newV0);
            else if (pdgV0 == 3122)
                Lambdas.push_back(newV0);
            else if (pdgV0 == 310)
                KaonsZeroShort.push_back(newV0);
            else if (pdgV0 == 422)
                PionPairs.push_back(newV0);
            FillV0(1, newV0);
            // InfoF("(pdg=%i) mass=%f, radius=%f", pdgV0, newV0.Mass(), newV0.Radius());
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
    if (pdgV0 == -3122)
        InfoF("N Found AntiLambdas: %lu", AntiLambdas.size());
    else if (pdgV0 == 3122)
        InfoF("N Found Lambdas: %lu", Lambdas.size());
    else if (pdgV0 == 310)
        InfoF("N Found KaonsZeroShort: %lu", KaonsZeroShort.size());
    else if (pdgV0 == 422)
        InfoF("N Found PionPairs: %lu", PionPairs.size());
}

/*
 *
 */
void Manager::CollectTrueV0(Int_t pdgHypothesis, UInt_t esdIdxNeg, UInt_t esdIdxPos, Bool_t& isTrue, Int_t& mcIdxV0, Int_t& mcPdgCode,
                            Bool_t& isSecondary, Bool_t& isSignal, Int_t& reactionID, Bool_t& isHybrid) {
    // /* -- By default */
    // isTrue = kFALSE;
    // mcIdxV0 = -1;
    // mcPdgCode = 0;
    // isSecondary = kFALSE;
    // isSignal = kFALSE;
    // reactionID = -1;
    // isHybrid = kFALSE;
    // /* -- Get MC particles linked to tracks */
    // UInt_t mcIdxNeg, mcIdxPos;
    // Bool_t hasMcNeg = GetMcIdx(esdIdxNeg, mcIdxNeg);
    // Bool_t hasMcPos = GetMcIdx(esdIdxPos, mcIdxPos);
    // /* -- Get mothers of MC particles */
    // UInt_t mcIdxMotherNeg, mcIdxMotherPos;
    // Bool_t hasMotherNeg = hasMcNeg && GetMotherMcIdx(mcIdxNeg, mcIdxMotherNeg);
    // Bool_t hasMotherPos = hasMcPos && GetMotherMcIdx(mcIdxPos, mcIdxMotherPos);
    // /* Note: different treatment for Pion Pairs */
    // if (pdgHypothesis == 422) {
    // /* -- Determine if it's a true V0 */
    // isTrue = !hasMotherNeg && !hasMotherPos && GetPdgCode(mcIdxNeg) == -211 && GetPdgCode(mcIdxPos) == 211;
    // if (isTrue) {
    // isSecondary = IsSecondary(mcIdxNeg) && IsSecondary(mcIdxPos);
    // isSignal = IsSignal(mcIdxNeg) && IsSignal(mcIdxPos) && GetReactionID(mcIdxNeg) == GetReactionID(mcIdxPos);
    // if (isSignal)
    // reactionID = GetReactionID(mcIdxNeg);
    // else
    // isHybrid = (IsSignal(mcIdxNeg) && !IsSignal(mcIdxPos)) || (!IsSignal(mcIdxNeg) && IsSignal(mcIdxPos));
    // }
    // } else {
    // /* -- Determine if it's a true V0 */
    // isTrue = hasMotherNeg && hasMotherPos && mcIdxMotherNeg == mcIdxMotherPos && getPdgCode_fromMcIdx[mcIdxMotherNeg] == pdgHypothesis;
    // if (isTrue) {
    // mcIdxV0 = mcIdxMotherNeg;
    // mcPdgCode = getPdgCode_fromMcIdx[mcIdxV0];
    // isSecondary = isMcIdxSecondary[mcIdxV0];
    // isSignal = IsSignal(mcIdxV0);
    // if (isSignal)
    // reactionID = GetReactionID(mcIdxV0);
    // else
    // isHybrid = (IsSignal(mcIdxNeg) && !IsSignal(mcIdxPos)) || (!IsSignal(mcIdxNeg) && IsSignal(mcIdxPos));
    // }
    // }
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/*
 *
 */
void Manager::ProcessFindableSexaquarks() {
    //
    std::unordered_multiset<Int_t> PdgCodes_Expected = {-2212, 211, -211, 211};
    /*
    // PENDING
    Int_t PdgCode_StruckNucleon = 2112;
    if (fReactionChannel == "ALK0") {
        PdgCode_StruckNucleon = 2112;
        PdgCodes_FinalStateProducts = {-2212, 211, -211, 211};
    } else if (fReactionChannel == "ALPK") {
        PdgCodes_FinalStateProducts = {-2212, 211, 321};
    } else if (fReactionChannel == "ALPKPP") {
        PdgCodes_FinalStateProducts = {-2212, 211, 321, -211, 211};
    } else if (fReactionChannel == "PKPKX") {
        PdgCodes_FinalStateProducts = {321, 321};
    }
    */

    fTree_Injected->SetBranchStatus("*", 0);
    fTree_Injected->SetBranchStatus("ReactionID", 1);

    fTree_MC->SetBranchStatus("*", 0);
    fTree_MC->SetBranchStatus("PdgCode", 1);

    fTree_Tracks->SetBranchStatus("*", 0);
    fTree_Tracks->SetBranchStatus("Px", 1);
    fTree_Tracks->SetBranchStatus("Py", 1);
    fTree_Tracks->SetBranchStatus("Pz", 1);

    for (Long64_t Reaction_Entry = 0; Reaction_Entry < GetN_Injected(); Reaction_Entry++) {
        if (!ReadInjected(Reaction_Entry)) continue;
        std::vector<Long64_t> Track_Entries;
        std::unordered_multiset<Int_t> PdgCodes_Reconstructed;
        /* Loop over MC */
        for (Long64_t& MC_Entry : getMcEntries_fromReactionID[Injected.ReactionID]) {
            if (!ReadMCParticle(MC_Entry)) continue;
            /* -- Condition: should've been reconstructed */
            Long64_t Track_Entry = GetTrackEntryFromMcEntry(MC_Entry);
            if (!ReadTrack(Track_Entry)) continue;
            Track_Entries.push_back(Track_Entry);
            PdgCodes_Reconstructed.insert(MC.PdgCode);
        }
        /* It should match the expected final state products */
        if (PdgCodes_Reconstructed != PdgCodes_Expected) continue;
        /* Loop over Tracks */
        for (Long64_t& Track_Entry : Track_Entries) {
            if (!ReadTrack(Track_Entry)) continue;
            if (!ReadMCParticle(GetMcEntryFromTrackEntry(Track_Entry))) continue;
            InfoF("ReactionID=%u, PdgCode=%i, Px=%f, Py=%f, Pz=%f", Injected.ReactionID, MC.PdgCode, Track.Px, Track.Py, Track.Pz);
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
    getMcIdx_FromMcEntry.clear();
    getMcEntry_fromMcIdx.clear();
    mcEntriesOfTrueV0s.clear();
    getMcIdxNeg_fromMcIdx.clear();
    getMcIdxPos_fromMcIdx.clear();
    getMcEntries_fromReactionID.clear();
    //
    getTrackEntry_fromEsdIdx.clear();
    getEsdIdx_fromTrackEntry.clear();
    getMcIdx_fromEsdIdx.clear();
    getEsdIdx_fromMcIdx.clear();
    //
    TrackEntries_AntiProton.clear();
    TrackEntries_Proton.clear();
    TrackEntries_NegKaon.clear();
    TrackEntries_PosKaon.clear();
    TrackEntries_PiMinus.clear();
    TrackEntries_PiPlus.clear();
    //
    AntiLambdas.clear();
    Lambdas.clear();
    KaonsZeroShort.clear();
    PionPairs.clear();
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
