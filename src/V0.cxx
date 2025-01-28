#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Particle {

/*
 *
 */
V0::V0(KFParticle kfV0, ROOT::Math::PxPyPzEVector lvV0, Int_t pdgV0, UInt_t esdIdxNeg, ROOT::Math::PxPyPzEVector lvNeg, UInt_t esdIdxPos,
       ROOT::Math::PxPyPzEVector lvPos)
    : KFParticle(kfV0),
      lvV0(lvV0),
      PdgCode(pdgV0),
      EsdIdxNeg(esdIdxNeg),
      lvNeg(lvNeg),
      EsdIdxPos(esdIdxPos),
      lvPos(lvPos),
      IsTrue(kFALSE),
      McIdxV0(-1),
      McPdgCode(0),
      IsSecondary(kFALSE),
      IsSignal(kFALSE),
      ReactionID(-1),
      IsHybrid(kFALSE) {}

/*
 *
 */
V0::V0(KFParticle kfV0, ROOT::Math::PxPyPzEVector lvV0, Int_t pdgV0, UInt_t esdIdxNeg, ROOT::Math::PxPyPzEVector lvNeg, UInt_t esdIdxPos,
       ROOT::Math::PxPyPzEVector lvPos, Bool_t isTrue, Int_t mcIdxV0, Int_t mcPdgCode, Bool_t isSecondary, Bool_t isSignal, Int_t ReactionID,
       Bool_t isHybrid)
    : KFParticle(kfV0),
      lvV0(lvV0),
      PdgCode(pdgV0),
      EsdIdxNeg(esdIdxNeg),
      lvNeg(lvNeg),
      EsdIdxPos(esdIdxPos),
      lvPos(lvPos),
      IsTrue(isTrue),
      McIdxV0(mcIdxV0),
      McPdgCode(mcPdgCode),
      IsSecondary(isSecondary),
      IsSignal(isSignal),
      ReactionID(ReactionID),
      IsHybrid(isHybrid) {}

}  // namespace Particle
}  // namespace Tree2Sexaquark
