#ifndef T2S_TREES_WRITER_HXX
#define T2S_TREES_WRITER_HXX

#include "TTree.h"

#include "Utilities/Logger.hxx"

#include "Particles/V0.hxx"
#include "Trees/Output.hxx"

namespace Tree2Sexaquark {

class Writer {
   public:
    Writer() = default;
    ~Writer() = default;

    /* V0s */

    V0_tt V0;

    void InitV0sTree() { fTree_V0s = new TTree("V0s", "V0s"); }
    void InitV0sBranches();
    void FillV0(UInt_t idxV0, Particle::V0 thisV0);
    void WriteV0sTree() {
        if (!fTree_V0s->GetDirectory()) {
            ErrorF("TTree %s has no directory", fTree_V0s->GetName());
            return;
        }
        fTree_V0s->Write();
        InfoF("TTree %s written to TFile %s", fTree_V0s->GetName(), fTree_V0s->GetDirectory()->GetName());
    }

    /* Sexaquarks */

    SexaquarkA_tt Sexaquark_ALK0;
    SexaquarkD_tt Sexaquark_ALPK;
    SexaquarkE_tt Sexaquark_ALPKPP;
    KaonPair_tt Sexaquark_PKPKX;

    SexaquarkA_tt Sexaquark_LK0;
    SexaquarkD_tt Sexaquark_LNK;
    SexaquarkE_tt Sexaquark_LNKPP;
    KaonPair_tt Sexaquark_NKNKX;

    void InitSexaquarkTrees() {
        //
        fTree_Sexaquarks_ALK0 = new TTree("Sexaquarks_ALK0", "Sexaquarks_ALK0");
        fTree_Sexaquarks_ALPK = new TTree("Sexaquarks_ALPK", "Sexaquarks_ALPK");
        fTree_Sexaquarks_ALPKPP = new TTree("Sexaquarks_ALPKPP", "Sexaquarks_ALPKPP");
        fTree_Sexaquarks_PKPKX = new TTree("Sexaquarks_PKPKX", "Sexaquarks_PKPKX");
        //
        fTree_Sexaquarks_LK0 = new TTree("Sexaquarks_LK0", "Sexaquarks_LK0");
        fTree_Sexaquarks_LNK = new TTree("Sexaquarks_LNK", "Sexaquarks_LNK");
        fTree_Sexaquarks_LNKPP = new TTree("Sexaquarks_LNKPP", "Sexaquarks_LNKPP");
        fTree_Sexaquarks_NKNKX = new TTree("Sexaquarks_NKNKX", "Sexaquarks_NKNKX");
    }
    void InitSexaquarkBranches_TypeA();
    void InitSexaquarkBranches_TypeD();
    void InitSexaquarkBranches_TypeE();
    void InitKaonPairBranches();

   private:
    /* -- V0s */
    TTree* fTree_V0s;
    /* -- AntiSexaquark Candidates */
    TTree* fTree_Sexaquarks_ALK0;
    TTree* fTree_Sexaquarks_ALPK;
    TTree* fTree_Sexaquarks_ALPKPP;
    TTree* fTree_Sexaquarks_PKPKX;
    /* -- Fake Sexaquark Candidates */
    TTree* fTree_Sexaquarks_LK0;
    TTree* fTree_Sexaquarks_LNK;
    TTree* fTree_Sexaquarks_LNKPP;
    TTree* fTree_Sexaquarks_NKNKX;
};

}  // namespace Tree2Sexaquark

#endif  // T2S_TREES_READER_HXX
