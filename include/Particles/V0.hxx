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

namespace V0 {
enum class Species {
    Lambda = 3122,
    KaonZero = 310,
    PionPair = 422,
};
}  // namespace V0

namespace Particle {

/*
 *
 */
class V0 : public Base {
   public:
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV)
        : Base(lvV0, kfV0, kfPV),  //
          PdgCode(PdgCode),
          EsdIdxNeg(EsdIdxNeg),
          EsdIdxPos(EsdIdxPos),
          /*  */
          lvNeg(lvNeg),
          lvPos(lvPos),
          kfNeg(kfNeg),
          kfPos(kfPos),
          /*  */
          IsTrue(0),
          McIdxV0(0),
          McPdgCode(0),
          IsSecondary(0),
          IsSignal(0),
          ReactionID(0),
          IsHybrid(0) {}
    V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV,                                //
       Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t reactionID, Bool_t isHybrid)
        : Base(lvV0, kfV0, kfPV),  //
          PdgCode(PdgCode),
          EsdIdxNeg(EsdIdxNeg),
          EsdIdxPos(EsdIdxPos),
          /*  */
          lvNeg(lvNeg),
          lvPos(lvPos),
          kfNeg(kfNeg),
          kfPos(kfPos),
          /*  */
          IsTrue(isTrue),
          McIdxV0(mcIdxV0),
          McPdgCode(mcPdgCode),
          IsSecondary(isSecondary),
          IsSignal(isSignal),
          ReactionID(reactionID),
          IsHybrid(isHybrid) {}
    ~V0() = default;

    inline Double_t NegPx() { return lvNeg.Px(); }
    inline Double_t NegPy() { return lvNeg.Py(); }
    inline Double_t NegPz() { return lvNeg.Pz(); }

    inline Double_t PosPx() { return lvPos.Px(); }
    inline Double_t PosPy() { return lvPos.Py(); }
    inline Double_t PosPz() { return lvPos.Pz(); }

    inline Double_t DCAbtwDau() { return TMath::Abs((Double_t)kfNeg.GetDistanceFromParticle(kfPos)); }
    inline Double_t DCAnegV0() { return TMath::Abs((Double_t)kfNeg.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCAposV0() { return TMath::Abs((Double_t)kfPos.GetDistanceFromVertex(kfThis)); }
    inline Double_t ArmQt() { return Math::ArmenterosQt(GetMomentumVector(), lvNeg.Vect()); }
    inline Double_t ArmAlpha() { return Math::ArmenterosAlpha(GetMomentumVector(), lvNeg.Vect(), lvPos.Vect()); }
    inline Double_t ArmQtOverAlpha() { return ArmQt() / TMath::Abs(ArmAlpha()); }

    typedef Double_t (V0::*MemFn)();

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
