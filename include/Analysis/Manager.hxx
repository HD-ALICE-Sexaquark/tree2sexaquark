#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include <unordered_map>

#include "TDatabasePDG.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TString.h"

#include "Trees/Reader.hxx"

class AnalysisManager : public Reader {
   public:
    AnalysisManager() = default;
    AnalysisManager(TString InputFile_Path, Bool_t IsMC, Long64_t LimitToNEvents);
    ~AnalysisManager() = default;

    void Print();
    Bool_t OpenInputFile();
    Bool_t IsMC() { return fIsMC; }

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

    /* Injected AntiSexaquark-Nucleon Interactions */
    void ProcessInjected();

    /* MC Particles */
    void ProcessMCParticles();

    /* Tracks */
    void ProcessTracks();

    /* V0s */
    void ProcessFindableV0s();
    Int_t GetN_TrueV0s() { return (Int_t)mcIndicesOfTrueV0s.size(); }

    /* Sexaquarks */
    void ProcessFindableSexaquarks();

    /* Related to Containers */
    /* -- filled at `ProcessMCParticles()` */
    Bool_t GetMcEntry(UInt_t mcIdx, Long64_t& mcEntry) {
        if (getMcEntry_fromMcIdx.find(mcIdx) == getMcEntry_fromMcIdx.end()) return kFALSE;
        mcEntry = getMcEntry_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t GetPdgCode(UInt_t mcIdx, Int_t& pdgCode) {
        if (getPdgCode_fromMcIdx.find(mcIdx) == getPdgCode_fromMcIdx.end()) return kFALSE;
        pdgCode = getPdgCode_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t IsSignal(UInt_t mcIdx) {
        if (isMcIdxSignal.find(mcIdx) == isMcIdxSignal.end()) return kFALSE;
        return isMcIdxSignal[mcIdx];
    }
    Bool_t IsSecondary(UInt_t mcIdx) {
        if (isMcIdxSecondary.find(mcIdx) == isMcIdxSecondary.end()) return kFALSE;
        return isMcIdxSecondary[mcIdx];
    }
    Bool_t GetReactionID(UInt_t mcIdx, UInt_t& reactionID) {
        if (getReactionID_fromMcIdx.find(mcIdx) == getReactionID_fromMcIdx.end()) return kFALSE;
        reactionID = getReactionID_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t GetMcIndices(UInt_t reactionID, std::vector<UInt_t>& mcIndices) {
        if (getMcIndices_fromReactionID.find(reactionID) == getMcIndices_fromReactionID.end()) return kFALSE;
        mcIndices = getMcIndices_fromReactionID[reactionID];
        return kTRUE;
    }
    Bool_t GetMotherMcIdx(UInt_t mcIdx, UInt_t& motherMcIdx) {
        if (getMotherMcIdx_fromMcIdx.find(mcIdx) == getMotherMcIdx_fromMcIdx.end()) return kFALSE;
        if (getMotherMcIdx_fromMcIdx[mcIdx] == -1) return kFALSE;
        motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t GetAncestorMcIdx(UInt_t mcIdx, UInt_t& ancestorMcIdx) {
        if (getAncestorMcIdx_fromMcIdx.find(mcIdx) == getAncestorMcIdx_fromMcIdx.end()) return kFALSE;
        if (getAncestorMcIdx_fromMcIdx[mcIdx] == -1) return kFALSE;
        ancestorMcIdx = getAncestorMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t GetNegDauMcIdx(UInt_t mcIdx, UInt_t& negDauMcIdx) {
        if (getNegDauMcIdx_fromMcIdx.find(mcIdx) == getNegDauMcIdx_fromMcIdx.end()) return kFALSE;
        negDauMcIdx = getNegDauMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    Bool_t GetPosDauMcIdx(UInt_t mcIdx, UInt_t& posDauMcIdx) {
        if (getPosDauMcIdx_fromMcIdx.find(mcIdx) == getPosDauMcIdx_fromMcIdx.end()) return kFALSE;
        posDauMcIdx = getPosDauMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    /* -- filled at `ProcessTracks()` */
    Bool_t GetMcIdx(UInt_t esdIdx, Int_t& mcIdx) {
        if (getMcIdx_fromEsdIdx.find(esdIdx) == getMcIdx_fromEsdIdx.end()) return kFALSE;
        mcIdx = getMcIdx_fromEsdIdx[esdIdx];
        return kTRUE;
    }

    /* Utilities */
    void CleanContainers();
    void EndOfEvent();
    void EndOfAnalysis();

   private:
    /* Analysis Properties */
    TString fInputFile_Path;
    Bool_t fIsMC;
    Long64_t fLimitToNEvents;

    /* ROOT Objects */
    TDatabasePDG fPDG;
    /* -- File */
    TFile* fInputFile;
    /* -- Event */
    TString Event_UID;
    TDirectoryFile* Event_Dir;

    /* Containers */
    /* -- filled at `ProcessMCParticles()` */
    std::unordered_map<UInt_t, Long64_t> getMcEntry_fromMcIdx;                    // key: `mcIdx`, value: get MC entry (position within the tree) from
    std::unordered_map<UInt_t, Int_t> getPdgCode_fromMcIdx;                       // key: `mcIdx`
    std::unordered_map<UInt_t, Bool_t> isMcIdxSignal;                             // key: `mcIdx`
    std::unordered_map<UInt_t, Bool_t> isMcIdxSecondary;                          // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getReactionID_fromMcIdx;                   // key: `mcIdx`
    std::unordered_map<UInt_t, std::vector<UInt_t>> getMcIndices_fromReactionID;  // key: `ReactionID`
    std::unordered_map<UInt_t, UInt_t> getMotherMcIdx_fromMcIdx;                  // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getAncestorMcIdx_fromMcIdx;                // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getNegDauMcIdx_fromMcIdx;                  // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getPosDauMcIdx_fromMcIdx;                  // key: `mcIdx`
    /* -- filled at `ProcessTracks()` */
    std::unordered_map<UInt_t, UInt_t> getMcIdx_fromEsdIdx;  // key: `esdIdx`
    /* -- used in `ProcessFindableV0s()` */
    std::vector<UInt_t> mcIndicesOfTrueV0s;
    std::unordered_map<UInt_t, UInt_t> getNegDauEsdIdx_fromMcIdx;  // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getPosDauEsdIdx_fromMcIdx;  // key: `mcIdx`
    /* -- used in `ProcessFindableSexaquarks()` */
    std::unordered_map<UInt_t, std::vector<UInt_t>> getEsdIndices_fromReactionID;  // key: `ReactionID`
};

#endif  // T2S_ANALYSIS_MANAGER_HXX
