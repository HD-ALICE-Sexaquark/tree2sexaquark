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

class KFParticle;
class KFVertex;

namespace Tree2Sexaquark {

enum PdgCode : Int_t {
    PiMinus = -211,
    PiPlus = 211,
    NegKaon = -321,
    PosKaon = 321,
    AntiProton = -2212,
    Proton = 2212,
    AntiNeutron = -2112,
    Neutron = 2112,
    AntiLambda = -3122,
    Lambda = 3122,
    KaonZeroShort = 310,
    None = 0
};

namespace Analysis {

struct KF_Track {
    ULong64_t entry;
    KFParticle kf;
    PxPyPzMVector lv;
};

struct MC_Track {
    ULong64_t entry, mother_entry;
    UInt_t reaction_id;
    Int_t pdg_code;
    Bool_t is_secondary, is_signal;
};

struct KF_V0 {
    KF_Track neg, pos;
    ULong64_t idx;
    KFParticle kf;
    PxPyPzMVector lv;
};

struct MC_V0 {
    MC_Track neg, pos;
    Long64_t mc_entry;
    Int_t pdg_code;
    UInt_t reaction_id;
    Bool_t has_mc, is_true, is_signal, is_hybrid;
};

struct TypeA {
    KF_V0 v0a, v0b;
    KFParticle kf;
    PxPyPzMVector lv, lv_asdecay;
};

struct MC_TypeA {
    MC_V0 mc_v0a, mc_v0b;
    Bool_t is_signal, is_hybrid;
    UInt_t reaction_id;
    // ULong64_t ancestor; // PENDING
    // Bool_t same_ancestor, is_noncomb_bkg; // PENDING
};

struct TypeD {
    KF_V0 v0;
    KF_Track ba;
    KFParticle kf;
    PxPyPzMVector lv, lv_asdecay;
};

struct MC_TypeD {
    MC_V0 mc_v0;
    MC_Track mc_ba;
    Bool_t is_signal, is_hybrid;
    UInt_t reaction_id;
    // ULong64_t ancestor; // PENDING
    // Bool_t same_ancestor, is_noncomb_bkg; // PENDING
};

class Manager {
   public:
    Manager() = default;
    ~Manager() = default;

    void Init();
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
    static RVec<MC_Track> Linked_MC_Creator(cRVecUL linked_mc_entries_, cRVecL mother_mc_entries_,  //
                                            cRVecI pdg_code, cRVecI is_secondary, cRVecI is_signal, cRVecU reaction_id);
    static RNode ProcessTracks(RNode df);

    /* V0s */
    static RVec<KF_V0> V0s_KF_Finder(PdgCode pdg_code_v0, Double_t neg_mass, Double_t pos_mass,           //
                                     const RVec<KF_Track> &neg_tracks, const RVec<KF_Track> &pos_tracks,  //
                                     const Float_t &magnetic_field, const XYZPoint &v3_pv);
    static RVec<MC_V0> V0s_TrueInfoCollector(PdgCode pdg_code_v0, PdgCode pdg_code_neg, PdgCode pdg_code_pos,  //
                                             const RVec<KF_V0> &found_v0s, const RVec<MC_Track> &linked_mc,    //
                                             cRVecI pdg_code, cRVecI is_signal, cRVecU reaction_id);
    static Bool_t V0s_PassesCuts(PdgCode pdg_code_v0, const KF_V0 &v0, const XYZPoint &v3_pv);
    RNode FindV0s(RNode df, PdgCode pdg_code_v0, PdgCode pdg_code_neg, PdgCode pdg_code_pos);

    /* Sexaquarks */
    RNode FindSexaquarks(const RNode &data_frame, PdgCode pdg_struck_nucleon, const std::vector<PdgCode> &pdg_reaction_products) {
        if (TMath::Abs(pdg_struck_nucleon) == PdgCode::Neutron) {
            return FindSexaquarks_TypeA(data_frame, pdg_struck_nucleon, pdg_reaction_products);
        } else {  // if (TMath::Abs(struck_nucleon) == PdgCode::Proton)
            return FindSexaquarks_TypeD(data_frame, pdg_struck_nucleon, pdg_reaction_products);
        }
    }
    /* -- Type A */
    static RVec<TypeA> TypeA_KF_Finder(const RVec<KF_V0> &found_v0a, const RVec<KF_V0> &found_v0b,  //
                                       const Float_t &magnetic_field, const KFVertex &kf_pv);
    static RVec<MC_TypeA> TypeA_TrueInfoCollector(const RVec<TypeA> &found,  //
                                                  const RVec<MC_V0> &mc_v0a, const RVec<MC_V0> &mc_v0b);
    static Bool_t TypeA_PassesCuts(const TypeA &sexa, const KFVertex &kf_pv);
    /* -- Type D */
    static RVec<TypeD> TypeD_KF_Finder(const RVec<KF_V0> &found_v0s, const RVec<KF_Track> &bach_tracks,  //
                                       const Float_t &magnetic_field, const KFVertex &kf_pv);
    static RVec<MC_TypeD> TypeD_TrueInfoCollector(PdgCode pdg_code_bach, const RVec<TypeD> &found,  //
                                                  const RVec<MC_V0> &mc_v0, const RVec<MC_Track> &mc_ba);
    static Bool_t TypeD_PassesCuts(const TypeD &sexa, const KFVertex &kf_pv);

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
    RNode FindSexaquarks_TypeA(RNode df, PdgCode pdg_struck_nucleon, const std::vector<PdgCode> &pdg_reaction_products);
    RNode FindSexaquarks_TypeD(RNode df, PdgCode pdg_struck_nucleon, const std::vector<PdgCode> &pdg_reaction_products);

    /* Containers */
    std::vector<std::string> fAnalyzed_V0sNames;
    std::vector<std::string> fAnalyzed_Channels;
    std::map<PdgCode, std::string> fParticleName_;  // key: pdg_code
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_MANAGER_HXX
