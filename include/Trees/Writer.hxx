#ifndef T2S_TREES_WRITER_HXX
#define T2S_TREES_WRITER_HXX

#include "TTree.h"

#include "Utilities/Logger.hxx"

#include "Particles/ChannelA.hxx"
#include "Particles/ChannelD.hxx"
#include "Particles/ChannelE.hxx"
#include "Particles/KaonPair.hxx"
#include "Particles/V0.hxx"
#include "Trees/Output.hxx"
#include <unordered_map>

namespace Tree2Sexaquark {

class Writer {
   public:
    enum class TreeName : UInt_t {
        V0s,
        Sexaquarks_ALK0,
        Sexaquarks_ALPK,
        Sexaquarks_ALPKPP,
        Sexaquarks_PKPKX,
        Sexaquarks_LK0,
        Sexaquarks_LNK,
        Sexaquarks_LNKPP,
        Sexaquarks_NKNKX
    };

    Writer() = default;
    ~Writer() = default;

    void InitTree(TreeName tree_name) {
        fTree_[tree_name] = new TTree(to_string[tree_name], to_string[tree_name]);
        InitBranches(tree_name);
    }

    void FillV0(UInt_t idx_v0, Candidate::V0 new_v0);
    void FillSexaquark(TreeName tree_name, Candidate::ChannelA new_sexaquark);
    void FillSexaquark(TreeName tree_name, Candidate::ChannelD new_sexaquark);
    void FillSexaquark(TreeName tree_name, Candidate::ChannelE new_sexaquark);
    void FillKaonPair(TreeName tree_name, Candidate::KaonPair new_kk);

    void WriteTree(TreeName tree_name) {
        if (!fTree_[tree_name]->GetDirectory()) {
            ErrorF("TTree %s has no directory", fTree_[tree_name]->GetName());
            return;
        }
        fTree_[tree_name]->Write();
        InfoF("TTree %s written to TFile %s", fTree_[tree_name]->GetName(), fTree_[tree_name]->GetDirectory()->GetName());
    }

   protected:
    V0_tt V0;

    SexaquarkA_tt Sexaquark_ALK0;
    SexaquarkD_tt Sexaquark_ALPK;
    SexaquarkE_tt Sexaquark_ALPKPP;
    KaonPair_tt Sexaquark_PKPKX;

    SexaquarkA_tt Sexaquark_LK0;
    SexaquarkD_tt Sexaquark_LNK;
    SexaquarkE_tt Sexaquark_LNKPP;
    KaonPair_tt Sexaquark_NKNKX;

   private:
    void InitBranches(TreeName tree_name) {
        if (tree_name == TreeName::V0s) {
            InitV0sBranches();
        } else if (tree_name == TreeName::Sexaquarks_ALK0 || tree_name == TreeName::Sexaquarks_LK0) {
            InitSexaquarkBranches_TypeA(tree_name);
        } else if (tree_name == TreeName::Sexaquarks_ALPK || tree_name == TreeName::Sexaquarks_LNK) {
            InitSexaquarkBranches_TypeD(tree_name);
        } else if (tree_name == TreeName::Sexaquarks_ALPKPP || tree_name == TreeName::Sexaquarks_LNKPP) {
            InitSexaquarkBranches_TypeE(tree_name);
        } else if (tree_name == TreeName::Sexaquarks_PKPKX || tree_name == TreeName::Sexaquarks_NKNKX) {
            InitKaonPairBranches(tree_name);
        } else {
            ErrorF("TreeName %s not recognized", to_string[tree_name].Data());
        }
    }
    void InitV0sBranches();
    void InitSexaquarkBranches_TypeA(TreeName tree_name);
    void InitSexaquarkBranches_TypeD(TreeName tree_name);
    void InitSexaquarkBranches_TypeE(TreeName tree_name);
    void InitKaonPairBranches(TreeName tree_name);

    std::unordered_map<TreeName, TTree*> fTree_;

    std::unordered_map<TreeName, TString> to_string = {{TreeName::V0s, "V0s"},
                                                       {TreeName::Sexaquarks_ALK0, "Sexaquarks_ALK0"},
                                                       {TreeName::Sexaquarks_ALPK, "Sexaquarks_ALPK"},
                                                       {TreeName::Sexaquarks_ALPKPP, "Sexaquarks_ALPKPP"},
                                                       {TreeName::Sexaquarks_PKPKX, "Sexaquarks_PKPKX"},
                                                       {TreeName::Sexaquarks_LK0, "Sexaquarks_LK0"},
                                                       {TreeName::Sexaquarks_LNK, "Sexaquarks_LNK"},
                                                       {TreeName::Sexaquarks_LNKPP, "Sexaquarks_LNKPP"},
                                                       {TreeName::Sexaquarks_NKNKX, "Sexaquarks_NKNKX"}};
};

}  // namespace Tree2Sexaquark

#endif  // T2S_TREES_READER_HXX
