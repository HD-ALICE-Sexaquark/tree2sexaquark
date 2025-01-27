#ifndef T2S_TREES_WRITER_HXX
#define T2S_TREES_WRITER_HXX

#include "TTree.h"

// ROOT::Math libraries
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Utilities/Logger.hxx"

#include "Trees/Output.hxx"

using namespace ROOT;

class Writer {
   public:
    Writer() = default;
    ~Writer() = default;

    /* V0s */

    V0_tt V0;

    void InitV0sTree() { fTree_V0s = new TTree("V0s", "V0s"); }
    void InitV0sBranches(Bool_t IsMC);
    void FillV0(UInt_t idxV0, UInt_t esdIdxNegDau, UInt_t esdIdxPosDau, Int_t pdgV0, Math::PxPyPzEVector lvV0, KFParticle kfV0,
                Math::PxPyPzEVector lvNegDaughter, Math::PxPyPzEVector lvPosDaughter,  //
                Bool_t AnalysisMC = kFALSE, Bool_t isTrue = kFALSE, Int_t mcIdxV0 = -1, Int_t mcPdgCode = -1, Bool_t isSecondary = kFALSE,
                Bool_t isSignal = kFALSE, Int_t ReactionID = -1, Bool_t isHybrid = kFALSE);
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
    void InitSexaquarkBranches_TypeA(Bool_t IsMC);
    void InitSexaquarkBranches_TypeD(Bool_t IsMC);
    void InitSexaquarkBranches_TypeE(Bool_t IsMC);
    void InitKaonPairBranches(Bool_t IsMC);

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

#endif  // T2S_TREES_READER_HXX
