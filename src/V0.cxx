#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Particle {

/*
 *
 */
V0::V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV)
    : Base(lvV0, kfV0, kfPV),
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

/*
 *
 */
V0::V0(Int_t PdgCode, UInt_t EsdIdxNeg, UInt_t EsdIdxPos,                                                 //
       ROOT::Math::PxPyPzEVector lvV0, ROOT::Math::PxPyPzEVector lvNeg, ROOT::Math::PxPyPzEVector lvPos,  //
       KFParticle kfV0, KFParticle kfNeg, KFParticle kfPos, KFVertex kfPV,                                //
       Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t reactionID, Bool_t isHybrid)
    : Base(lvV0, kfV0, kfPV),
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

}  // namespace Particle
}  // namespace Tree2Sexaquark
