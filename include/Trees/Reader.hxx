#ifndef T2S_TREES_READER_HXX
#define T2S_TREES_READER_HXX

#include "TTree.h"

#include "Trees/Structure.hxx"

class Reader {
   public:
    Reader() = default;
    ~Reader() = default;

    Event_tt Event;
    Injected_tt Injected;
    MC_tt MC;
    Track_tt Track;

    void SetEventsTree(TTree* tree) { fTree_Events = tree; }
    TTree* GetEventsTree() { return fTree_Events; }
    void ConnectEventBranches(Bool_t IsMC);
    void DisconnectEventBranches() {
        if (fTree_Events) fTree_Events->ResetBranchAddresses();
    }
    virtual Long64_t GetN_Events() { return fTree_Events->GetEntries(); }
    Int_t ReadEvent(Long64_t event_entry) { return fTree_Events->GetEntry(event_entry); }

    void SetInjectedTree(TTree* tree) { fTree_Injected = tree; }
    void ConnectInjectedBranches();
    void DisconnectInjectedBranches() {
        if (fTree_Injected) fTree_Injected->ResetBranchAddresses();
    }
    Long64_t GetN_Injected() { return fTree_Injected->GetEntries(); }
    Int_t ReadInjected(Long64_t sexa_entry) { return fTree_Injected->GetEntry(sexa_entry); }

    void SetMCTree(TTree* tree) { fTree_MC = tree; }
    void ConnectMCBranches();
    void DisconnectMCBranches() {
        if (fTree_MC) fTree_MC->ResetBranchAddresses();
    }
    Long64_t GetN_MCParticles() { return fTree_MC->GetEntries(); }
    Int_t ReadMCParticle(Long64_t mc_entry) { return fTree_MC->GetEntry(mc_entry); }

    void SetTracksTree(TTree* tree) { fTree_Tracks = tree; }
    void ConnectTracksBranches();
    void DisconnectTracksBranches() {
        if (fTree_Tracks) fTree_Tracks->ResetBranchAddresses();
    }
    Long64_t GetN_Tracks() { return fTree_Tracks->GetEntries(); }
    Int_t ReadTrack(Long64_t track_entry) { return fTree_Tracks->GetEntry(track_entry); }

   private:
    TTree* fTree_Events;
    TTree* fTree_Injected;
    TTree* fTree_MC;
    TTree* fTree_Tracks;
};

#endif  // T2S_TREES_READER_HXX
