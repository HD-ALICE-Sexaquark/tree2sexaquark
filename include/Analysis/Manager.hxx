#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include "ROOT/RResultPtr.hxx"
#include "RtypesCore.h"
#include "TDatabasePDG.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TString.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFVertex.h"

#include "Utilities/Logger.hxx"

#include "Analysis/Settings.hxx"
// #include "Cuts/Inspector.hxx"
// #include "Particles/V0.hxx"
// #include "Trees/Writer.hxx"

using namespace ROOT::RDF;
using RDataFrame = ROOT::RDataFrame;

using namespace ROOT::VecOps;
using RVecI = ROOT::RVecI;  // reserved for masks and indices that could be -1
using cRVecI = const ROOT::RVecI &;
using RVecL = ROOT::RVecL;
using cRVecL = const ROOT::RVecL &;
using RVecU = ROOT::RVecU;
using cRVecU = const ROOT::RVecU &;
using RVecUL = ROOT::RVecUL;  // reserved for protected indices
using cRVecUL = const ROOT::RVecUL &;
using RVecF = ROOT::RVecF;
using cRVecF = const ROOT::RVecF &;

using PxPyPzMVector = ROOT::Math::PxPyPzMVector;
using XYZPoint = ROOT::Math::XYZPoint;

namespace Tree2Sexaquark {
namespace Analysis {

struct FoundV0_TrueInfo {
    ULong64_t neg, pos, entry;
    Bool_t same_mother, is_true, is_signal, is_secondary, is_hybrid;
    Int_t neg_pdg_code, pos_pdg_code, pdg_code;
    UInt_t reaction_id;
};

struct FoundV0 {
    ULong64_t idx, neg, pos;
    KFParticle kf, kf_neg, kf_pos;
    PxPyPzMVector lv, lv_neg, lv_pos;
};

struct TypeA_TrueInfo {
    FoundV0_TrueInfo v0a_mc, v0b_mc;
    Bool_t is_signal, is_hybrid;
    UInt_t reaction_id;
    // ULong64_t ancestor; // PENDING
    // Bool_t same_ancestor, is_noncomb_bkg; // PENDING
};

struct TypeA {
    FoundV0 v0a, v0b;
    KFParticle kf;
    PxPyPzMVector lv, lv_asdecay;
};

class Manager {
   public:
    Manager() = default;
    ~Manager() = default;

    void Init();
    RNode ProcessEvent(RNode df);

    /* Injected AntiSexaquark-Nucleon Interactions */
    RNode ProcessInjected(RNode df);

    /* MC Particles */
    RNode ProcessMCParticles(RNode df);

    /* Tracks */
    RNode ProcessTracks(RNode df);

    /* V0s */
    RNode FindV0s(RNode df, Int_t pdg_v0, Int_t pdg_neg, Int_t pdg_pos);

    /* Sexaquarks */
    RNode FindSexaquarks(RNode df, Int_t pdg_struck_nucleon, std::vector<Int_t> pdg_reaction_products) {
        if (TMath::Abs(pdg_struck_nucleon) == 2112) {
            return FindSexaquarks_TypeA(df, pdg_struck_nucleon < 0);
        }
        if (pdg_reaction_products.size() > 2) {
            return FindSexaquarks_TypeDE(df, pdg_struck_nucleon < 0);
        }
        return FindSexaquarks_TypeH(df, pdg_struck_nucleon < 0);
    }

    void CollectTrueInfo_ChannelA();

    /* Cuts */

    /* Vector Gymnastics */
    static RVecI Mask(cRVecUL entries, size_t reference_size);
    template <typename T>
    static RVec<T> Extract(const RVec<T> &property, cRVecUL link_);
    template <typename T>
    static RVec<T> ExtractIf(const RVec<T> &property, cRVecL link_, cRVecI link_protection, T default_value);
    template <typename T>
    static RVec<RVec<T>> ExtractVector(const RVec<T> &property, const RVec<RVecUL> &entries_);
    template <typename T>
    static RVec<T> ExtractVector_First(const RVec<T> &property, const RVec<RVecUL> &entries_);
    template <typename T>
    static RVec<T> ExtractVector_Sum(const RVec<T> &property, const RVec<RVecUL> &entries_);

    /* Utilities */
    inline Float_t GetMass(Int_t pdg_code) { return TDatabasePDG::Instance()->GetParticle(pdg_code)->Mass(); }
    void PrintAll(RNode df);
    void EndOfAnalysis(RNode df);

   private:
    /* Functions */
    RNode FindSexaquarks_TypeA(RNode df, Bool_t anti_channel);
    RNode FindSexaquarks_TypeDE(RNode df, Bool_t anti_channel);
    RNode FindSexaquarks_TypeH(RNode df, Bool_t anti_channel);

    std::vector<std::string> fAnalyzed_V0sNames;
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_MANAGER_HXX
