#ifndef T2S_TREES_READER_HXX
#define T2S_TREES_READER_HXX

#include "TTree.h"

#include "Utilities/Logger.hxx"

#include "Trees/Structure.hxx"

class Reader {
   public:
    Reader();
    virtual ~Reader();

    TTree* GetEventsTree() { return fTree_Events; }
    void SetEventsTree(TTree* tree) { fTree_Events = tree; }
    void ConnectEventBranches(Bool_t IsMC);
    virtual Long64_t GetN_Events() { return fTree_Events->GetEntries(); }
    Int_t ReadEvent(Long64_t event_idx) { return fTree_Events->GetEntry(event_idx); }
    Event_tt Event;

    void ConnectInjectedBranches();

    void SetInjectedTree(TTree* tree) { fTree_Injected = tree; }
    void DisconnectInjectedBranches() { fTree_Injected->ResetBranchAddresses(); }
    Long64_t GetN_Injected() { return fTree_Injected->GetEntries(); }
    Int_t ReadInjected(Long64_t sexa_idx) { return fTree_Injected->GetEntry(sexa_idx); }
    Injected_tt Injected;

    void ConnectMCBranches();

    void SetMCTree(TTree* tree) { fTree_MC = tree; }
    void DisconnectMCBranches() { fTree_MC->ResetBranchAddresses(); }
    Long64_t GetN_MCParticles() { return fTree_MC->GetEntries(); }
    Int_t ReadMCParticle(Long64_t mc_idx) { return fTree_MC->GetEntry(mc_idx); }
    MC_tt MC;

    void ConnectTracksBranches();

    void SetTracksTree(TTree* tree) { fTree_Tracks = tree; }
    void DisconnectTracksBranches() { fTree_Tracks->ResetBranchAddresses(); }
    Long64_t GetN_Tracks() { return fTree_Tracks->GetEntries(); }
    Int_t ReadTrack(Long64_t track_idx) { return fTree_Tracks->GetEntry(track_idx); }
    Track_tt Track;

   private:
    TTree* fTree_Events;
    TTree* fTree_Injected;
    TTree* fTree_MC;
    TTree* fTree_Tracks;
};

#endif  // T2S_TREES_READER_HXX
