#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include <memory>
#include <unordered_map>
#include <vector>

#include "RtypesCore.h"
#include "TDatabasePDG.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TString.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Utilities/Logger.hxx"

#include "Analysis/Settings.hxx"
#include "Trees/Reader.hxx"

using namespace ROOT;

KFParticle CreateKFParticle(Track_tt track, Double_t mass);
// KFVertex CreateKFVertex(Double_t* XYZ, Double_t* CovarianceMatrix);
KFVertex CreateKFVertex(Double_t* XYZ);
KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Double_t massThis, Int_t chargeThis);

class AnalysisManager : public Reader {
   public:
    AnalysisManager(Settings_tt Settings);

    void Print();
    Bool_t OpenInputFile();
    Bool_t IsMC() { return Settings.IsMC; }
    Bool_t IsSignalMC() { return Settings.IsSignalMC; }

    inline Long64_t GetN_Events() {
        if (!Settings.LimitToNEvents) return GetEventsTree()->GetEntries();
        return Settings.LimitToNEvents;
    }

    inline Bool_t GetEvent(Long64_t evt_idx) {
        if (!ReadEvent(evt_idx)) {
            DebugF("Event # %lld couldn't be read, moving on...", evt_idx);
            return kFALSE;
        }
        if (Settings.IsMC) {
            if (Settings.IsSignalMC)
                Event_UID = TString::Format("A18_%u_%04u_%03u", Event.RunNumber, Event.DirNumber, Event.EventNumber);
            else
                Event_UID = TString::Format("BKG_%6u_%04u_%03u", Event.RunNumber, Event.DirNumber, Event.EventNumber);
        } else {
            Event_UID = TString::Format("DATA_%6u_%03u_%u_%03u", Event.RunNumber, Event.DirNumber, Event.DirNumberB, Event.EventNumber);
        }
        InfoF("Processing Event # %lld (UID = %s)", evt_idx, Event_UID.Data());
        InfoF(">> Centrality = %f, PV = (%f, %f, %f), B = %f", Event.Centrality, Event.PV_Xv, Event.PV_Yv, Event.PV_Zv, Event.MagneticField);

        Event_Dir = std::unique_ptr<TDirectoryFile>(InputFile->Get<TDirectoryFile>(Event_UID));
        if (!Event_Dir) {
            DebugF("TDirectoryFile %s couldn't be found, moving on...", Event_UID.Data());
            return kFALSE;
        }

        return kTRUE;
    }

    inline TTree* FindTreeInFile(TString tree_name) {
        TTree* AuxTree = InputFile->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TFile %s", tree_name.Data(), InputFile->GetName());
            return nullptr;
        }
        DebugF("TTree %s found in TFile %s", tree_name.Data(), InputFile->GetName());
        return AuxTree;
    }

    inline TTree* FindTreeInEventDir(TString tree_name) {
        TTree* AuxTree = Event_Dir->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TDirectoryFile %s", tree_name.Data(), Event_Dir->GetName());
            return nullptr;
        }
        DebugF("TTree %s found in TDirectoryFile %s", tree_name.Data(), Event_Dir->GetName());
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
    void KalmanV0Finder(Int_t pdgNegDaughter, Int_t pdgPosDaughter, Int_t pdgV0 = -1);

    /* Sexaquarks */
    void ProcessFindableSexaquarks();
    void KalmanSexaquarkFinder(Int_t pdgStruckNucleon, std::vector<Int_t> pdgReactionProducts);
    void KalmanSexaquarkFinder_TypeA(std::vector<Int_t> pdgReactionProducts);
    void KalmanSexaquarkFinder_TypeDE(std::vector<Int_t> pdgReactionProducts);
    void KalmanSexaquarkFinder_TypeH(std::vector<Int_t> pdgReactionProducts);

    /* Related to Containers */
    /* -- filled at `ProcessMCParticles()` */
    inline Bool_t GetMcEntry(UInt_t mcIdx, Long64_t& mcEntry) {
        if (getMcEntry_fromMcIdx.find(mcIdx) == getMcEntry_fromMcIdx.end()) return kFALSE;
        mcEntry = getMcEntry_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetPdgCode(UInt_t mcIdx, Int_t& pdgCode) {
        if (getPdgCode_fromMcIdx.find(mcIdx) == getPdgCode_fromMcIdx.end()) return kFALSE;
        pdgCode = getPdgCode_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t IsSignal(UInt_t mcIdx) {
        if (isMcIdxSignal.find(mcIdx) == isMcIdxSignal.end()) return kFALSE;
        return isMcIdxSignal[mcIdx];
    }
    inline Bool_t IsSecondary(UInt_t mcIdx) {
        if (isMcIdxSecondary.find(mcIdx) == isMcIdxSecondary.end()) return kFALSE;
        return isMcIdxSecondary[mcIdx];
    }
    inline Bool_t GetReactionID(UInt_t mcIdx, UInt_t& reactionID) {
        if (getReactionID_fromMcIdx.find(mcIdx) == getReactionID_fromMcIdx.end()) return kFALSE;
        if (getReactionID_fromMcIdx[mcIdx] == -1) return kFALSE;
        reactionID = getReactionID_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetMcIndices(UInt_t reactionID, std::vector<UInt_t>& mcIndices) {
        if (getMcIndices_fromReactionID.find(reactionID) == getMcIndices_fromReactionID.end()) return kFALSE;
        mcIndices = getMcIndices_fromReactionID[reactionID];
        return kTRUE;
    }
    inline Bool_t GetMotherMcIdx(UInt_t mcIdx, UInt_t& motherMcIdx) {
        if (getMotherMcIdx_fromMcIdx.find(mcIdx) == getMotherMcIdx_fromMcIdx.end()) return kFALSE;
        if (getMotherMcIdx_fromMcIdx[mcIdx] == -1) return kFALSE;
        motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetAncestorMcIdx(UInt_t mcIdx, UInt_t& ancestorMcIdx) {
        if (getAncestorMcIdx_fromMcIdx.find(mcIdx) == getAncestorMcIdx_fromMcIdx.end()) return kFALSE;
        if (getAncestorMcIdx_fromMcIdx[mcIdx] == -1) return kFALSE;
        ancestorMcIdx = getAncestorMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetNegDauMcIdx(UInt_t mcIdx, UInt_t& negDauMcIdx) {
        if (getNegDauMcIdx_fromMcIdx.find(mcIdx) == getNegDauMcIdx_fromMcIdx.end()) return kFALSE;
        negDauMcIdx = getNegDauMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetPosDauMcIdx(UInt_t mcIdx, UInt_t& posDauMcIdx) {
        if (getPosDauMcIdx_fromMcIdx.find(mcIdx) == getPosDauMcIdx_fromMcIdx.end()) return kFALSE;
        posDauMcIdx = getPosDauMcIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    /* -- filled at `ProcessTracks()` */
    inline Bool_t GetTrackEntry(UInt_t esdIdx, Long64_t& trackEntry) {
        if (getTrackEntry_fromEsdIdx.find(esdIdx) == getTrackEntry_fromEsdIdx.end()) return kFALSE;
        trackEntry = getTrackEntry_fromEsdIdx[esdIdx];
        return kTRUE;
    }
    inline Bool_t GetMcIdx(UInt_t esdIdx, UInt_t& mcIdx) {
        if (getMcIdx_fromEsdIdx.find(esdIdx) == getMcIdx_fromEsdIdx.end()) return kFALSE;
        mcIdx = getMcIdx_fromEsdIdx[esdIdx];
        return kTRUE;
    }
    inline Bool_t GetTrack(UInt_t esdIdx, Track_tt& track) {
        Long64_t trackEntry;
        if (!GetTrackEntry(esdIdx, trackEntry)) return kFALSE;
        if (!ReadTrack(trackEntry)) return kFALSE;
        track = Track;
        return kTRUE;
    }
    /* -- used in `ProcessFindableV0s()` */
    inline Bool_t GetNegDauEsdIdx(UInt_t mcIdx, UInt_t& negDauEsdIdx) {
        if (getNegDauEsdIdx_fromMcIdx.find(mcIdx) == getNegDauEsdIdx_fromMcIdx.end()) return kFALSE;
        negDauEsdIdx = getNegDauEsdIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }
    inline Bool_t GetPosDauEsdIdx(UInt_t mcIdx, UInt_t& posDauEsdIdx) {
        if (getPosDauEsdIdx_fromMcIdx.find(mcIdx) == getPosDauEsdIdx_fromMcIdx.end()) return kFALSE;
        posDauEsdIdx = getPosDauEsdIdx_fromMcIdx[mcIdx];
        return kTRUE;
    }

    /* Utilities */
    void CleanContainers();
    void EndOfEvent();
    void EndOfAnalysis();

   private:
    /* Analysis Properties */
    Settings_tt Settings;

    /* ROOT Objects */
    TDatabasePDG fPDG;
    /* -- File */
    std::unique_ptr<TFile> InputFile;
    /* -- Event */
    TString Event_UID;
    std::unique_ptr<TDirectoryFile> Event_Dir;

    /* Containers */
    /* -- filled at `ProcessMCParticles()` */
    std::unordered_map<UInt_t, Long64_t> getMcEntry_fromMcIdx;                    // key: `mcIdx`, value: get MC position within the tree
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
    std::unordered_map<UInt_t, Long64_t> getTrackEntry_fromEsdIdx;  // key: `esdIdx`, value: get track position within the tree
    std::unordered_map<UInt_t, UInt_t> getMcIdx_fromEsdIdx;         // key: `esdIdx`
    /* -- and looped over in `KalmanV0Finder()` and Sexaquark Finders */
    std::vector<UInt_t> esdIndicesOfAntiProtonTracks;
    std::vector<UInt_t> esdIndicesOfProtonTracks;
    std::vector<UInt_t> esdIndicesOfNegKaonTracks;
    std::vector<UInt_t> esdIndicesOfPosKaonTracks;
    std::vector<UInt_t> esdIndicesOfPiMinusTracks;
    std::vector<UInt_t> esdIndicesOfPiPlusTracks;
    /* -- used in `ProcessFindableV0s()` */
    std::vector<UInt_t> mcIndicesOfTrueV0s;
    std::unordered_map<UInt_t, UInt_t> getNegDauEsdIdx_fromMcIdx;  // key: `mcIdx`
    std::unordered_map<UInt_t, UInt_t> getPosDauEsdIdx_fromMcIdx;  // key: `mcIdx`
    /* -- used in `ProcessFindableSexaquarks()` */
    std::unordered_map<UInt_t, std::vector<UInt_t>> getEsdIndices_fromReactionID;  // key: `ReactionID`
};

#endif  // T2S_ANALYSIS_MANAGER_HXX
