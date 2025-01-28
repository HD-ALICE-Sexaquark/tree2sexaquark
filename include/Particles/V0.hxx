#ifndef T2S_V0S_HXX
#define T2S_V0S_HXX

// ROOT::Math libraries
#include "Math/Vector4D.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"

namespace Tree2Sexaquark {
namespace Particle {

class V0 : public KFParticle {
   public:
    V0() = default;
    V0(KFParticle kfV0, ROOT::Math::PxPyPzEVector lvV0, Int_t pdgV0, UInt_t esdIdxNeg, ROOT::Math::PxPyPzEVector lvNeg, UInt_t esdIdxPos,
       ROOT::Math::PxPyPzEVector lvPos);
    V0(KFParticle kfV0, ROOT::Math::PxPyPzEVector lvV0, Int_t pdgV0, UInt_t esdIdxNeg, ROOT::Math::PxPyPzEVector lvNeg, UInt_t esdIdxPos,
       ROOT::Math::PxPyPzEVector lvPos, Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t ReactionID,
       Bool_t isHybrid);
    ~V0() = default;

    inline Double_t Px() { return lvV0.Px(); }
    inline Double_t Py() { return lvV0.Py(); }
    inline Double_t Pz() { return lvV0.Pz(); }
    inline Double_t E() { return lvV0.E(); }

    inline Double_t NegPx() { return lvNeg.Px(); }
    inline Double_t NegPy() { return lvNeg.Py(); }
    inline Double_t NegPz() { return lvNeg.Pz(); }

    inline Double_t PosPx() { return lvPos.Px(); }
    inline Double_t PosPy() { return lvPos.Py(); }
    inline Double_t PosPz() { return lvPos.Pz(); }

    Int_t PdgCode;  // hypothetical V0 PDG code

    UInt_t EsdIdxNeg;
    UInt_t EsdIdxPos;

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
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_V0S_HXX
