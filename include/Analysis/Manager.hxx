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

struct V0_tt {
    ULong64_t neg, pos;
    KFParticle kf, kf_neg, kf_pos;
    PxPyPzMVector lv, lv_neg, lv_pos;
};

class Manager /* : public Writer */ {
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
    void KalmanSexaquarkFinder(Int_t pdg_struck_nucleon, std::vector<Int_t> pdg_reaction_products) {
        if (TMath::Abs(pdg_struck_nucleon) == 2112) {
            KalmanSexaquarkFinder_TypeA(pdg_struck_nucleon > 0);
        } else if (TMath::Abs(pdg_struck_nucleon) == 2212) {
            if (pdg_reaction_products.size() > 2)
                KalmanSexaquarkFinder_TypeDE(pdg_struck_nucleon > 0);
            else
                KalmanSexaquarkFinder_TypeH(pdg_struck_nucleon > 0);
        } else {
            ErrorF("Struck nucleon %s not recognized", TDatabasePDG::Instance()->GetParticle(pdg_struck_nucleon)->GetName());
        }
    }

    void CollectTrueInfo_ChannelA();
    void CollectTrueInfo_ChannelD();
    void CollectTrueInfo_ChannelE();
    void CollectTrueInfo_ChannelH();

    /* Cuts */
    // Cuts::Inspector Inspector;

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
    void KalmanSexaquarkFinder_TypeA(Bool_t anti_channel);
    void KalmanSexaquarkFinder_TypeDE(Bool_t anti_channel);
    void KalmanSexaquarkFinder_TypeH(Bool_t anti_channel);

    /* -- Event */
    KFVertex kfPrimaryVertex;  // primary vertex
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_MANAGER_HXX
