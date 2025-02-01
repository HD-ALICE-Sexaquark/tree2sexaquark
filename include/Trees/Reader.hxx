#ifndef T2S_TREES_READER_HXX
#define T2S_TREES_READER_HXX

#include "TTree.h"

#include "Trees/Input.hxx"

namespace Tree2Sexaquark {

class Reader {
   public:
    Reader() = default;
    ~Reader() = default;

    void SetEventsTree(TTree* tree) { fTree_Events = tree; }
    TTree* GetEventsTree() { return fTree_Events; }
    void ConnectEventBranches();
    void DisconnectEventBranches() {
        if (fTree_Events) fTree_Events->ResetBranchAddresses();
    }
    inline virtual Long64_t GetN_Events() { return fTree_Events->GetEntries(); }
    inline Int_t ReadEvent(Long64_t eventEntry) { return fTree_Events->GetEntry(eventEntry); }

    void SetInjectedTree(TTree* tree) { fTree_Injected = tree; }
    void ConnectInjectedBranches();
    void DisconnectInjectedBranches() {
        if (fTree_Injected) fTree_Injected->ResetBranchAddresses();
    }
    inline Long64_t GetN_Injected() { return fTree_Injected->GetEntries(); }
    inline Int_t ReadInjected(Long64_t reactionEntry) { return fTree_Injected->GetEntry(reactionEntry); }

    void SetMCTree(TTree* tree) { fTree_MC = tree; }
    void ConnectMCBranches();
    void DisconnectMCBranches() {
        if (fTree_MC) fTree_MC->ResetBranchAddresses();
    }
    inline Long64_t GetN_MCParticles() { return fTree_MC->GetEntries(); }
    inline Int_t ReadMCParticle(Long64_t mcEntry) { return fTree_MC->GetEntry(mcEntry); }

    void SetTracksTree(TTree* tree) { fTree_Tracks = tree; }
    void ConnectTracksBranches();
    void DisconnectTracksBranches() {
        if (fTree_Tracks) fTree_Tracks->ResetBranchAddresses();
    }
    inline Long64_t GetN_Tracks() { return fTree_Tracks->GetEntries(); }
    inline Int_t ReadTrack(Long64_t trackEntry) { return fTree_Tracks->GetEntry(trackEntry); }

   protected:
    Event_tt Event;
    Injected_tt Injected;
    MC_tt MC;
    Track_tt Track;

    TTree* fTree_Events;
    TTree* fTree_Injected;
    TTree* fTree_MC;
    TTree* fTree_Tracks;
};

}  // namespace Tree2Sexaquark

#endif  // T2S_TREES_READER_HXX
