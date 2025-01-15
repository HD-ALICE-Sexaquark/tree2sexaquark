#include "Analysis/Manager.hxx"

/*
 * Default constructor (with arguments)
 */
AnalysisManager::AnalysisManager(TString InputFile_Path, Bool_t IsMC, Long64_t LimitToNEvents)
    : Reader(),  //
      /*  */
      fInputFile_Path(InputFile_Path),
      fIsMC(IsMC),
      fLimitToNEvents(LimitToNEvents),
      /*  */
      fPDG(),
      fInputFile(nullptr),
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
      getMcIdx_fromEsdIdx(),
      getEsdIndices_fromReactionID(),
      mcIndicesOfTrueV0s(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx() {}

/*
 *
 */
void AnalysisManager::Print() {
    //
    InfoF("IsMC           = %i", (Int_t)fIsMC);
    InfoF("InputFile_Path = %s", fInputFile_Path.Data());
    InfoF("LimitToNEvents = %lld", fLimitToNEvents);
}

/*
 *
 */
Bool_t AnalysisManager::OpenInputFile() {
    //
    fInputFile = TFile::Open(fInputFile_Path, "READ");
    if (!fInputFile || fInputFile->IsZombie()) {
        ErrorF("TFile %s couldn't be opened", fInputFile_Path.Data());
        return kFALSE;
    }
    DebugF("TFile %s opened successfully", fInputFile_Path.Data());

    TTree* EventsTree = FindTreeIn(fInputFile, "Events");
    if (!EventsTree) return kFALSE;

    SetEventsTree(EventsTree);
    ConnectEventBranches(IsMC());

    return kTRUE;
}

/*
 *
 */
Bool_t AnalysisManager::GetEvent(Long64_t evt_idx) {
    //
    if (!ReadEvent(evt_idx)) {
        DebugF("Event # %lld couldn't be read, moving on...", evt_idx);
        return kFALSE;
    }
    Event_UID = TString::Format("A18_%i_%04u_%03u", Event.RunNumber, Event.EventNumber, Event.DirNumber);
    InfoF("Processing Event # %lld (UID = %s)", evt_idx, Event_UID.Data());

    Event_Dir = fInputFile->Get<TDirectoryFile>(Event_UID);
    if (!Event_Dir) {
        DebugF("TDirectoryFile %s couldn't be found, moving on...", Event_UID.Data());
        return kFALSE;
    }

    return kTRUE;
}

/*
 *
 */
void AnalysisManager::ProcessInjected() {
    //
    TTree* InjectedTree = FindTreeIn(Event_Dir, "Injected");
    if (!InjectedTree) return;

    SetInjectedTree(InjectedTree);
    ConnectInjectedBranches();

    for (Long64_t sexa_entry = 0; sexa_entry < GetN_Injected(); sexa_entry++) {
        if (!ReadInjected(sexa_entry)) continue;
        InfoF("%i, %f, %f, %f", Injected.ReactionID, Injected.Px, Injected.Py, Injected.Pz);
    }
}

/*
 *
 */
void AnalysisManager::ProcessMCParticles() {
    //
    TTree* MCTree = FindTreeIn(Event_Dir, "MC");
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
        InfoF("%i, %i, %i, %i, %i", MC.Idx, MC.Idx_Mother, MC.Idx_Ancestor, MC.PdgCode, MC.ReactionID);
    }
}

/*
 *
 */
void AnalysisManager::ProcessTracks() {
    //
    TTree* TracksTree = FindTreeIn(Event_Dir, "Tracks");
    if (!TracksTree) return;

    SetTracksTree(TracksTree);
    ConnectTracksBranches();

    Int_t True_PdgCode;

    for (Long64_t track_entry = 0; track_entry < GetN_Tracks(); track_entry++) {
        if (!ReadTrack(track_entry)) continue;
        /* Fill containers */
        getMcIdx_fromEsdIdx[Track.Idx] = Track.Idx_True;
        /*  */
        if (!GetPdgCode(Track.Idx_True, True_PdgCode)) continue;
        if (True_PdgCode != 2212) continue;
        InfoF("%i, %i, %i, %f, %f, %f, %f", Track.Idx, Track.Idx_True, True_PdgCode, Track.Px, Track.Py, Track.Pz, Track.NSigmaProton);
    }
}

/*
 *
 */
void AnalysisManager::ProcessFindableV0s() {
    //
}

/*
 *
 */
void AnalysisManager::ProcessFindableSexaquarks() {
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
    getMcIdx_fromEsdIdx.clear();
    //
    getEsdIndices_fromReactionID.clear();
    //
    mcIndicesOfTrueV0s.clear();
    getNegDauEsdIdx_fromMcIdx.clear();
    getPosDauEsdIdx_fromMcIdx.clear();
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
    fInputFile->Close();
    delete fInputFile;
}
