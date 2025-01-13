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
    virtual UInt_t GetN_Events() { return fTree_Events->GetEntries(); }
    Int_t ReadEvent(UInt_t event_idx) { return fTree_Events->GetEntry(event_idx); }
    Event_tt Event;

    void ConnectInjectedBranches();

    void SetInjectedTree(TTree* tree) { fTree_Injected = tree; }
    void DisconnectInjectedBranches() { fTree_Injected->ResetBranchAddresses(); }
    UInt_t GetN_Injected() { return fTree_Injected->GetEntries(); }
    Int_t ReadInjected(UInt_t sexa_idx) { return fTree_Injected->GetEntry(sexa_idx); }
    Injected_tt Injected;

    void ConnectMCBranches();

    void SetMCTree(TTree* tree) { fTree_MC = tree; }
    void DisconnectMCBranches() { fTree_MC->ResetBranchAddresses(); }
    UInt_t GetN_MCParticles() { return fTree_MC->GetEntries(); }
    Int_t ReadMCParticle(UInt_t mc_idx) { return fTree_MC->GetEntry(mc_idx); }
    MC_tt MC;

    void ConnectTracksBranches();

    void SetTracksTree(TTree* tree) { fTree_Tracks = tree; }
    void DisconnectTracksBranches() { fTree_Tracks->ResetBranchAddresses(); }
    UInt_t GetN_Tracks() { return fTree_Tracks->GetEntries(); }
    Int_t ReadTrack(UInt_t track_idx) { return fTree_Tracks->GetEntry(track_idx); }
    Track_tt Track;

   private:
    TTree* fTree_Events;
    TTree* fTree_Injected;
    TTree* fTree_MC;
    TTree* fTree_Tracks;
};

#endif  // T2S_TREES_READER_HXX
