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

namespace Tree2Sexaquark {
namespace Particle {

/*
 *
 */
class V0 {
   public:
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV);
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV,                                //
       Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t reactionID, Bool_t isHybrid);
    ~V0() = default;

    inline Double_t Px() { return lvV0.Px(); }
    inline Double_t Py() { return lvV0.Py(); }
    inline Double_t Pz() { return lvV0.Pz(); }
    inline Double_t E() { return lvV0.E(); }
    inline Double_t Eta() { return lvV0.Eta(); }
    inline Double_t Pt() { return lvV0.Pt(); }
    inline Double_t Mass() { return lvV0.M(); }

    inline Float_t Xv() { return kfV0.GetX(); }
    inline Float_t Yv() { return kfV0.GetY(); }
    inline Float_t Zv() { return kfV0.GetZ(); }

    inline Double_t NegPx() { return lvNeg.Px(); }
    inline Double_t NegPy() { return lvNeg.Py(); }
    inline Double_t NegPz() { return lvNeg.Pz(); }

    inline Double_t PosPx() { return lvPos.Px(); }
    inline Double_t PosPy() { return lvPos.Py(); }
    inline Double_t PosPz() { return lvPos.Pz(); }

    inline Double_t Radius() { return v3V0.Rho(); }
    inline Double_t DistFromPV() { return (v3V0 - v3PV).R(); }
    inline Double_t CPAwrtPV() { return Math::CosinePointingAngle(lvV0.Vect(), v3V0, v3PV); }
    inline Double_t DCAwrtPV() { return TMath::Abs(kfV0.GetDistanceFromVertex(kfPV)); }
    inline Double_t DCAbtwDau() { return TMath::Abs(kfNeg.GetDistanceFromParticle(kfPos)); }
    inline Double_t DCAnegV0() { return TMath::Abs(kfNeg.GetDistanceFromVertex(kfV0)); }
    inline Double_t DCAposV0() { return TMath::Abs(kfPos.GetDistanceFromVertex(kfV0)); }
    inline Double_t ArmQt() { return Math::ArmenterosQt(lvV0.Vect(), lvNeg.Vect()); }
    inline Double_t ArmAlpha() { return Math::ArmenterosAlpha(lvV0.Vect(), lvNeg.Vect(), lvPos.Vect()); }
    inline Double_t ArmQtOverAlpha() { return ArmQt() / TMath::Abs(ArmAlpha()); }

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
    ROOT::Math::PxPyPzEVector lvV0;
    ROOT::Math::PxPyPzEVector lvNeg;
    ROOT::Math::PxPyPzEVector lvPos;

    KFParticle kfV0;
    KFParticle kfNeg;
    KFParticle kfPos;
    KFVertex kfPV;

    ROOT::Math::XYZPoint v3V0;
    ROOT::Math::XYZPoint v3PV;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_V0S_HXX
