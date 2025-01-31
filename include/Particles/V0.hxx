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
namespace Particle {

/*
 *
 */
class V0 : public Base {
   public:
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV);
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV,                                //
       Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t reactionID, Bool_t isHybrid);
    ~V0() = default;

    inline Double_t NegPx() { return lvNeg.Px(); }
    inline Double_t NegPy() { return lvNeg.Py(); }
    inline Double_t NegPz() { return lvNeg.Pz(); }

    inline Double_t PosPx() { return lvPos.Px(); }
    inline Double_t PosPy() { return lvPos.Py(); }
    inline Double_t PosPz() { return lvPos.Pz(); }

    inline Double_t DCAbtwDau() { return TMath::Abs(kfNeg.GetDistanceFromParticle(kfPos)); }
    inline Double_t DCAnegV0() { return TMath::Abs(kfNeg.GetDistanceFromVertex(GetKFParticle())); }
    inline Double_t DCAposV0() { return TMath::Abs(kfPos.GetDistanceFromVertex(GetKFParticle())); }
    inline Double_t ArmQt() { return Math::ArmenterosQt(GetMomentumVector(), lvNeg.Vect()); }
    inline Double_t ArmAlpha() { return Math::ArmenterosAlpha(GetMomentumVector(), lvNeg.Vect(), lvPos.Vect()); }
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
    ROOT::Math::PxPyPzEVector lvNeg;
    ROOT::Math::PxPyPzEVector lvPos;

    KFParticle kfNeg;
    KFParticle kfPos;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_V0S_HXX
