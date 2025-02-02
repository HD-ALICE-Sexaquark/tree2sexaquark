#ifndef T2S_V0S_HXX
#define T2S_V0S_HXX

#include "TMath.h"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Math/Common.hxx"
#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class V0 : public Base {
   public:
    enum class Species : UInt_t {
        Lambda = 3122,
        KaonZeroShort = 310,
        PionPair = 422,
    };

    void SetV0Info(Int_t pdgHypothesis, UInt_t esdIdxNeg, UInt_t esdIdxPos) {
        PdgCode = pdgHypothesis;
        EsdIdxNeg = esdIdxNeg;
        EsdIdxPos = esdIdxPos;
    }
    void SetKinematics(ROOT::Math::PxPyPzEVector lv_v0, ROOT::Math::PxPyPzEVector lv_neg, ROOT::Math::PxPyPzEVector lv_pos) {
        lvThis = lv_v0;
        lvNeg = lv_neg;
        lvPos = lv_pos;
    }
    void SetGeometry(KFParticle kf_v0, KFParticle kf_neg, KFParticle kf_pos) {
        kfThis = kf_v0;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());
        kfNeg = kf_neg;
        kfPos = kf_pos;
    }
    void SetTrueInfo(Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t reactionID, Bool_t isHybrid) {
        IsTrue = isTrue;
        McIdxV0 = mcIdxV0;
        McPdgCode = mcPdgCode;
        IsSecondary = isSecondary;
        IsSignal = isSignal;
        ReactionID = reactionID;
        IsHybrid = isHybrid;
    }

    inline Double_t NegPx() { return lvNeg.Px(); }
    inline Double_t NegPy() { return lvNeg.Py(); }
    inline Double_t NegPz() { return lvNeg.Pz(); }

    inline Double_t NegXv() { return (Double_t)kfNeg.GetX(); }
    inline Double_t NegYv() { return (Double_t)kfNeg.GetY(); }
    inline Double_t NegZv() { return (Double_t)kfNeg.GetZ(); }

    inline Double_t PosPx() { return lvPos.Px(); }
    inline Double_t PosPy() { return lvPos.Py(); }
    inline Double_t PosPz() { return lvPos.Pz(); }

    inline Double_t PosXv() { return (Double_t)kfPos.GetX(); }
    inline Double_t PosYv() { return (Double_t)kfPos.GetY(); }
    inline Double_t PosZv() { return (Double_t)kfPos.GetZ(); }

    inline Double_t DCAbtwDau() { return TMath::Abs((Double_t)kfNeg.GetDistanceFromParticle(kfPos)); }
    inline Double_t DCAnegV0() { return TMath::Abs((Double_t)kfNeg.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCAposV0() { return TMath::Abs((Double_t)kfPos.GetDistanceFromVertex(kfThis)); }
    inline Double_t ArmQt() { return Math::ArmenterosQt(lvThis.Vect(), lvNeg.Vect()); }
    inline Double_t ArmAlpha() { return Math::ArmenterosAlpha(lvThis.Vect(), lvNeg.Vect(), lvPos.Vect()); }
    inline Double_t ArmQtOverAlpha() { return ArmQt() / TMath::Abs(ArmAlpha()); }

    V0::Species GetSpecies() { return static_cast<V0::Species>(TMath::Abs(PdgCode)); }

    Int_t PdgCode;  // hypothetical V0 PDG code
    UInt_t EsdIdxNeg;
    UInt_t EsdIdxPos;

    /* -- True Information */
    Bool_t IsTrue;
    Int_t McIdxV0;
    Int_t McPdgCode;
    Bool_t IsSecondary;
    Bool_t IsSignal;
    Int_t ReactionID;
    Bool_t IsHybrid;

   private:
    ROOT::Math::PxPyPzEVector lvNeg;
    ROOT::Math::PxPyPzEVector lvPos;

    KFParticle kfNeg;
    KFParticle kfPos;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_V0S_HXX
