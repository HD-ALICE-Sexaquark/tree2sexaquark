#ifndef T2S_TREES_MANAGER_HXX
#define T2S_TREES_MANAGER_HXX

#include "TTree.h"

#include "Utilities/Logger.hxx"

#include "Trees/Structure.hxx"

class Manager {
   public:
    static Manager* GetManager();
    static void DeleteManager();
    void Speak(Int_t a);

    TTree* GetEventsTree() { return fTree_Events; }
    void SetEventsTree(TTree* tree) { fTree_Events = tree; }
    void ConnectEventBranches();
    UInt_t GetNEvents() { return fTree_Events->GetEntries(); }
    void GetEvent(UInt_t event_index) { fTree_Events->GetEntry(event_index); }
    Event_tt Event;

   private:
    Manager();
    virtual ~Manager();
    static Manager* Instance;

    TTree* fTree_Events;
};

#endif  // T2S_TREES_MANAGER_HXX
