#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include "RtypesCore.h"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"

using namespace ROOT::RDF;
using RDataFrame = ROOT::RDataFrame;

#include "ROOT/RVec.hxx"

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
using PxPyPzEVector = ROOT::Math::PxPyPzEVector;
using XYZPoint = ROOT::Math::XYZPoint;

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

namespace Tree2Sexaquark {

namespace PdgCode {
const Int_t PiMinus = -211;
const Int_t PiPlus = 211;
const Int_t NegKaon = -321;
const Int_t PosKaon = 321;
const Int_t AntiProton = -2212;
const Int_t Proton = 2212;
const Int_t AntiNeutron = -2112;
const Int_t Neutron = 2112;
const Int_t AntiLambda = -3122;
const Int_t Lambda = 3122;
const Int_t KaonZeroShort = 310;
}  // namespace PdgCode

namespace Analysis {

struct KF_Track {
    ULong64_t entry;
    KFParticle kf;
    PxPyPzMVector lv;
};

struct FoundV0 {
    KF_Track neg, pos;
    ULong64_t idx;
    KFParticle kf;
    PxPyPzMVector lv;
};

struct FoundV0_TrueInfo {
    ULong64_t neg, pos, entry;
    Bool_t same_mother, is_true, is_signal, is_secondary, is_hybrid;
    Int_t neg_pdg_code, pos_pdg_code, pdg_code;
    UInt_t reaction_id;
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

struct TypeD_TrueInfo {
    FoundV0_TrueInfo v0a_mc, v0b_mc;
    Bool_t is_signal, is_hybrid;
    UInt_t reaction_id;
    // ULong64_t ancestor; // PENDING
    // Bool_t same_ancestor, is_noncomb_bkg; // PENDING
};

struct TypeD {
    FoundV0 v0a, v0b;
    KFParticle kf;
    PxPyPzMVector lv, lv_asdecay;
};

class Manager {
   public:
    Manager() = default;
    ~Manager() = default;

    static void Init();
    static RNode ProcessEvent(RNode df);

    /* Injected AntiSexaquark-Nucleon Interactions */
    static RNode ProcessInjected(RNode df);

    /* MC Particles */
    static RNode ProcessMCParticles(RNode df);

    /* Tracks */
    static RVec<KF_Track> Tracks_KF_Creator(Double_t pdg_mass, cRVecUL entries_,                                    //
                                            cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,          //
                                            cRVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,  //
                                            const RVec<RVecF> &cov_matrix, const Float_t &magnetic_field);
    static RNode ProcessTracks(RNode df);

    /* V0s */
    static RVec<FoundV0> V0s_KF_Finder(Int_t pdg_code_v0, Double_t neg_mass, Double_t pos_mass,             //
                                       const RVec<KF_Track> &neg_tracks, const RVec<KF_Track> &pos_tracks,  //
                                       const Float_t &magnetic_field, const XYZPoint &v3_pv);
    static RVec<FoundV0_TrueInfo> V0s_TrueInfoCollector(Int_t pdg_code_v0, Double_t pdg_code_neg, Double_t pdg_code_pos,  //
                                                        const RVec<FoundV0> &found_v0s,                                   //
                                                        cRVecUL track_mc_entry_, cRVecL mother_mc_entry_,                 //
                                                        cRVecI pdg_code, cRVecI mc_is_secondary, cRVecI mc_is_signal, cRVecU mc_reaction_id);
    static Bool_t V0s_PassesCuts(Int_t pdg_code_v0, const FoundV0 &v0, const XYZPoint &v3_pv);
    RNode FindV0s(RNode df, Int_t pdg_code_v0, Int_t pdg_code_neg, Int_t pdg_code_pos);

    /* Sexaquarks */
    RNode FindSexaquarks(const RNode &data_frame, Int_t pdg_struck_nucleon, const std::vector<Int_t> &pdg_reaction_products) {
        if (TMath::Abs(pdg_struck_nucleon) == PdgCode::Neutron) {
            return FindSexaquarks_TypeA(data_frame, pdg_struck_nucleon < 0);
        }
        if (pdg_reaction_products.size() > 2) {
            return FindSexaquarks_TypeDE(data_frame, pdg_struck_nucleon < 0);
        }
        return FindSexaquarks_TypeH(data_frame, pdg_struck_nucleon < 0);
    }
    /* -- Type A */
    static RVec<TypeA> TypeA_KF_Finder(const RVec<FoundV0> &found_v0a, const RVec<FoundV0> &found_v0b,  //
                                       const Float_t &magnetic_field, const KFVertex &kf_pv);
    static RVec<TypeA_TrueInfo> TypeA_TrueInfoCollector(const RVec<TypeA> &found,  //
                                                        const RVec<FoundV0_TrueInfo> &v0a_mc, const RVec<FoundV0_TrueInfo> &v0b_mc);
    static Bool_t TypeA_PassesCuts(const TypeA &sexa, const KFVertex &kf_pv);
    /* -- Type D */
    static RVec<TypeD> TypeD_KF_Finder();
    static RVec<TypeD_TrueInfo> TypeD_TrueInfoCollector();
    static Bool_t TypeD_PassesCuts();

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
    static void PrintAll(RNode df);
    void EndOfAnalysis(RNode df);

   private:
    /* Functions */
    RNode FindSexaquarks_TypeA(RNode df, Bool_t anti_channel);
    RNode FindSexaquarks_TypeDE(RNode df, Bool_t anti_channel);
    RNode FindSexaquarks_TypeH(RNode df, Bool_t anti_channel);

    std::vector<std::string> fAnalyzed_V0sNames;
    std::vector<std::string> fAnalyzed_Channels;
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_MANAGER_HXX
