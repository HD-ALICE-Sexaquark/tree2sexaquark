#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include "TDirectoryFile.h"
#include "TFile.h"
#include "TString.h"

#include "Trees/Reader.hxx"

class AnalysisManager : public Reader {
   public:
    AnalysisManager();
    AnalysisManager(TString InputFile_Path, Bool_t IsMC, Long64_t LimitToNEvents);
    ~AnalysisManager() = default;

    void Print();
    Bool_t IsMC() { return fIsMC; }

    Bool_t OpenInputFile();

    Long64_t GetN_Events() {
        if (!fLimitToNEvents) return GetEventsTree()->GetEntries();
        return fLimitToNEvents;
    }
    Bool_t GetEvent(Long64_t evt_idx);

    TTree* FindTreeIn(TFile* InputFile, TString tree_name) {
        TTree* AuxTree = InputFile->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TFile %s", tree_name.Data(), InputFile->GetName());
            return nullptr;
        }
        DebugF("TTree %s found in TFile %s", tree_name.Data(), InputFile->GetName());
        return AuxTree;
    }

    TTree* FindTreeIn(TDirectoryFile* InputDir, TString tree_name) {
        TTree* AuxTree = InputDir->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TDirectoryFile %s", tree_name.Data(), InputDir->GetName());
            return nullptr;
        }
        DebugF("TTree %s found in TDirectoryFile %s", tree_name.Data(), InputDir->GetName());
        return AuxTree;
    }

    void ProcessInjected();
    void ProcessMCParticles();
    void ProcessTracks();

   private:
    /* Analysis Properties */
    TString fInputFile_Path;
    Bool_t fIsMC;
    Long64_t fLimitToNEvents;

    /* -- File */
    TFile* fInputFile;
    /* -- Event */
    TString Event_UID;
    TDirectoryFile* Event_Dir;
};

#endif  // T2S_ANALYSIS_MANAGER_HXX
