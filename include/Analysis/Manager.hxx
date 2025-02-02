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
#include "KFVertex.h"

#include "Utilities/Logger.hxx"

#include "Analysis/Settings.hxx"
#include "Cuts/Inspector.hxx"
#include "Particles/V0.hxx"
#include "Trees/Reader.hxx"
#include "Trees/Writer.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

class Manager : public Reader, public Writer {
   public:
    Manager() = default;
    ~Manager() = default;

    Bool_t OpenInputFile();
    Bool_t PrepareOutputFile();

    /* Events */
    Bool_t GetEvent(Long64_t evt_idx);
    inline Long64_t GetN_Events() {
        if (!Settings::LimitToNEvents) return GetEventsTree()->GetEntries();
        return Settings::LimitToNEvents;
    }

    inline TTree* FindTreeInFile(TString tree_name) {
        TTree* AuxTree = InputFile->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TFile %s", tree_name.Data(), InputFile->GetName());
            return nullptr;
        }
        InfoF("TTree %s found in TFile %s", tree_name.Data(), InputFile->GetName());
        return AuxTree;
    }

    inline TTree* FindTreeInEventDir(TString tree_name) {
        TTree* AuxTree = Event_Dir->Get<TTree>(tree_name);
        if (!AuxTree) {
            ErrorF("TTree %s couldn't be found in TDirectoryFile %s", tree_name.Data(), Event_Dir->GetName());
            return nullptr;
        }
        InfoF("TTree %s found in TDirectoryFile %s", tree_name.Data(), Event_Dir->GetName());
        return AuxTree;
    }

    /* Injected AntiSexaquark-Nucleon Interactions */
    void ProcessInjected();

    /* MC Particles */
    void ProcessMCParticles();
    inline Bool_t MC_IsFinalStateSignal() { return MC.Generator == 2 && MC.Idx_Ancestor != -1; }
    inline Bool_t MC_IsSignal() { return MC.Generator == 2; }
    inline Bool_t MC_IsSecondary() { return MC.IsSecFromMat || MC.IsSecFromWeak || MC_IsSignal(); }
    /* -- translate between MC entries and indices */
    inline Long64_t GetMcEntry(UInt_t mc_idx) { return getMcEntry_fromMcIdx[mc_idx]; }
    inline UInt_t GetMcIdx(Long64_t mc_entry) { return getMcIdx_FromMcEntry[mc_entry]; }
    /* -- necessary for V0s */
    inline Long64_t GetMcEntryNeg(Long64_t mc_entry_mother) {
        UInt_t mcIdxMother = GetMcIdx(mc_entry_mother);
        if (getMcIdxNeg_fromMcIdx.find(mcIdxMother) == getMcIdxNeg_fromMcIdx.end()) return -1;
        UInt_t mcIdxNeg = getMcIdxNeg_fromMcIdx[mcIdxMother];
        return GetMcEntry(mcIdxNeg);
    }
    inline Long64_t GetMcEntryPos(Long64_t mc_entry_mother) {
        UInt_t mcIdxMother = GetMcIdx(mc_entry_mother);
        if (getMcIdxPos_fromMcIdx.find(mcIdxMother) == getMcIdxPos_fromMcIdx.end()) return -1;
        UInt_t mcIdxPos = getMcIdxPos_fromMcIdx[mcIdxMother];
        return GetMcEntry(mcIdxPos);
    }

    /* Tracks */
    void ProcessTracks();
    /* -- translate between Track entries and ESD indices */
    inline Long64_t GetTrackEntry(UInt_t esd_idx) { return getTrackEntry_fromEsdIdx[esd_idx]; }
    inline UInt_t GetEsdIdx(Long64_t track_entry) { return getEsdIdx_fromTrackEntry[track_entry]; }
    /* -- link between MC and Tracks */
    inline Long64_t GetTrackEntryFromMcEntry(Long64_t mc_entry) {
        UInt_t mcIdx = GetMcIdx(mc_entry);
        if (getEsdIdx_fromMcIdx.find(mcIdx) == getEsdIdx_fromMcIdx.end()) return -1;
        return GetTrackEntry(getEsdIdx_fromMcIdx[mcIdx]);
    }
    inline Long64_t GetMcEntryFromTrackEntry(Long64_t track_entry) {
        UInt_t esdIdx = GetEsdIdx(track_entry);
        if (getMcIdx_fromEsdIdx.find(esdIdx) == getMcIdx_fromEsdIdx.end()) return -1;
        return GetMcEntry(getMcIdx_fromEsdIdx[esdIdx]);
    }
    /*  */
    inline Bool_t CopyTrack(Long64_t track_entry, Track_tt& track) {
        if (!ReadTrack(track_entry)) return kFALSE;
        track = Track;
        return kTRUE;
    }

    /* V0s */
    void ProcessFindableV0s();
    void KalmanV0Finder(Int_t pdg_neg, Int_t pdg_pos, Int_t pdg_v0);
    void CollectTrueV0(Int_t pdg_hypothesis, UInt_t esd_idx_neg, UInt_t esd_idx_pos, Bool_t& is_true, Int_t& mc_idx_v0, Int_t& mc_pdg_code,
                       Bool_t& is_secondary, Bool_t& is_signal, Int_t& reaction_id, Bool_t& is_hybrid);

    /* Sexaquarks */
    void ProcessFindableSexaquarks();
    void KalmanSexaquarkFinder(Int_t pdg_struck_nucleon, std::vector<Int_t> pdg_reaction_products);
    void KalmanSexaquarkFinder_TypeA(std::vector<Int_t> pdg_reaction_products);
    void KalmanSexaquarkFinder_TypeDE(std::vector<Int_t> pdg_reaction_products);
    void KalmanSexaquarkFinder_TypeH(std::vector<Int_t> pdg_reaction_products);

    /* Cuts */
    Cuts::Inspector Inspector;

    /* Utilities */
    void CleanContainers();
    void EndOfEvent();
    void EndOfAnalysis();

   private:
    /* -- Files */
    std::unique_ptr<TFile> InputFile;
    std::unique_ptr<TFile> OutputFile;
    /* -- Event */
    TString Event_UID;
    std::unique_ptr<TDirectoryFile> Event_Dir;
    KFVertex kfPrimaryVertex;  // primary vertex

    /* ROOT Objects */
    TDatabasePDG fPDG;

    /* Containers */
    /* -- filled in `ProcessMCParticles()` */
    std::unordered_map<Long64_t, UInt_t> getMcIdx_FromMcEntry;
    std::unordered_map<UInt_t, Long64_t> getMcEntry_fromMcIdx;
    /*** -- and looped over in `ProcessFindableV0s()` */
    std::vector<Long64_t> mcEntriesOfTrueV0s;
    std::unordered_map<UInt_t, UInt_t> getMcIdxNeg_fromMcIdx;
    std::unordered_map<UInt_t, UInt_t> getMcIdxPos_fromMcIdx;
    /*** -- and looped over in `ProcessFindableSexaquarks()` */
    std::map<UInt_t, std::vector<Long64_t>> getMcEntries_fromReactionID;  // NOTE: final state particles
    /* -- filled in `ProcessTracks()` */
    std::unordered_map<UInt_t, Long64_t> getTrackEntry_fromEsdIdx;
    std::unordered_map<Long64_t, UInt_t> getEsdIdx_fromTrackEntry;
    /*** -- and looped over in `ProcessFindableV0s()` */
    std::unordered_map<UInt_t, UInt_t> getMcIdx_fromEsdIdx;
    std::unordered_map<UInt_t, UInt_t> getEsdIdx_fromMcIdx;
    /*** -- and looped over in `KalmanV0Finder()` and `KalmanSexaquarkFinder()` */
    std::vector<Long64_t> TrackEntries_AntiProton;
    std::vector<Long64_t> TrackEntries_Proton;
    std::vector<Long64_t> TrackEntries_NegKaon;
    std::vector<Long64_t> TrackEntries_PosKaon;
    std::vector<Long64_t> TrackEntries_PiMinus;
    std::vector<Long64_t> TrackEntries_PiPlus;
    /* -- filled in `KalmanV0Finder()` and looped over in `KalmanSexaquarkFinder()` */
    std::vector<Candidate::V0> AntiLambdas;
    std::vector<Candidate::V0> Lambdas;
    std::vector<Candidate::V0> KaonsZeroShort;
    std::vector<Candidate::V0> PionPairs;
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_MANAGER_HXX
