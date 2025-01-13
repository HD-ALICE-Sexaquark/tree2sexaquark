#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include "TDirectoryFile.h"
#include "TFile.h"
#include "TString.h"

#include "Trees/Reader.hxx"

class AnalysisManager : public Reader {
   public:
    static AnalysisManager* GetInstance();
    static void DeleteInstance();
    static AnalysisManager* Init() { return GetInstance(); }

    Bool_t Configure(Int_t argc, char* argv[]);
    Bool_t IsMC() { return fIsMC; }

    Bool_t OpenInputFile();

    UInt_t GetN_Events() {
        if (!fLimitToNEvents) return GetEventsTree()->GetEntries();
        return fLimitToNEvents;
    }
    Bool_t GetEvent(UInt_t evt_idx);

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
    AnalysisManager();
    virtual ~AnalysisManager();
    static AnalysisManager* Instance;

    /* Analysis Properties */
    Bool_t fIsMC = kTRUE;
    UInt_t fLimitToNEvents = 5;

    /* -- File */
    TString fInputFile_Path = "/home/ceres/borquez/some/esd2tree/task/attempts/local_signalMC_A1.8_18qr_test/SimpleTrees.root";
    TFile* fInputFile;

    /* Event */
    TString Event_UID;
    TDirectoryFile* Event_Dir;
};

#endif  // T2S_ANALYSIS_MANAGER_HXX
