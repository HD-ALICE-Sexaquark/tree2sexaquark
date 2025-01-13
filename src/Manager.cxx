
#include "Analysis/Manager.hxx"

AnalysisManager* AnalysisManager::Instance = nullptr;

/*
 * Get AnalysisManager singleton instance
 */
AnalysisManager* AnalysisManager::GetInstance() {
    if (Instance == nullptr) Instance = new AnalysisManager();
    return Instance;
}

/*
 * Delete AnalysisManager singleton instance
 */
void AnalysisManager::DeleteInstance() {
    if (Instance != nullptr) {
        delete Instance;
        Instance = nullptr;
    }
}

/*
 * Default private constructor
 */
AnalysisManager::AnalysisManager()
    : Reader(),  //
      fInputFile(nullptr),
      Event_UID(""),
      Event_Dir(nullptr) {
    // replace previous instance by this one
    if (Instance) delete Instance;
    Instance = this;
}

/*
 * Private destructor
 */
AnalysisManager::~AnalysisManager() { Instance = nullptr; }

/*
 *
 */
Bool_t AnalysisManager::Configure(Int_t argc, char* argv[]) {
    // Parser* ThisParser = new Parser();
    // IsMC = ThisParser->IsMC();
    InfoF("IsMC          = %i", (Int_t)fIsMC);
    InfoF("InputFile_Path = %s", fInputFile_Path.Data());
    InfoF("LimitToNEvents = %u", fLimitToNEvents);
    return kTRUE;
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
Bool_t AnalysisManager::GetEvent(UInt_t evt_idx) {
    //
    if (!ReadEvent(evt_idx)) {
        DebugF("Event # %u couldn't be read, moving on...", evt_idx);
        return kFALSE;
    }
    Event_UID = TString::Format("A18_%i_%04u_%03u", Event.RunNumber, Event.EventNumber, Event.DirNumber);
    InfoF("Processing Event # %u (UID = %s)", evt_idx, Event_UID.Data());

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

    for (UInt_t sexa_idx = 0; sexa_idx < GetN_Injected(); sexa_idx++) {
        if (!ReadInjected(sexa_idx)) continue;
        InfoF("%i, %f, %f, %f", Injected.ReactionID, Injected.Px, Injected.Py, Injected.Pz);
    }

    DisconnectInjectedBranches();
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

    for (UInt_t mc_idx = 0; mc_idx < GetN_MCParticles(); mc_idx++) {
        if (!ReadMCParticle(mc_idx)) continue;
        if (MC.Generator != 2 || MC.Idx_Ancestor != -1) continue;
        InfoF("%i, %i, %i, %i, %i", MC.Idx, MC.Idx_Mother, MC.Idx_Ancestor, MC.PdgCode, MC.ReactionID);
    }

    DisconnectMCBranches();
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

    for (UInt_t track_idx = 0; track_idx < GetN_Tracks(); track_idx++) {
        if (!ReadTrack(track_idx)) continue;
        InfoF("%i, %f, %f, %f", Track.Idx, Track.Px, Track.Py, Track.Pz);
    }

    DisconnectTracksBranches();
}
