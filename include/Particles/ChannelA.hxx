#ifndef T2S_SEXAQUARK_CHANNEL_A_HXX
#define T2S_SEXAQUARK_CHANNEL_A_HXX

#include "Math/Vector4D.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Sexaquark {

/*
 *
 */
class ChannelA : public Particle::Sexaquark {
   public:
    inline Double_t DCAbtwV0s() { return TMath::Abs((Double_t)kfLambda.GetDistanceFromParticle(kfKaonZero)); };

    inline Double_t DCAK0SV() { return TMath::Abs((Double_t)kfKaonZero.GetDistanceFromVertex(kfThis)); };

    inline Double_t DCAK0NegSV() { return TMath::Abs((Double_t)kfKaonZeroNeg.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAK0PosSV() { return TMath::Abs((Double_t)kfKaonZeroPos.GetDistanceFromVertex(kfThis)); };

    typedef Double_t (ChannelA::*MemFn)();

   private:
    ROOT::Math::PxPyPzEVector lvKaonZero;
    KFParticle kfKaonZero;
    KFParticle kfKaonZeroNeg;
    KFParticle kfKaonZeroPos;
};

}  // namespace Sexaquark
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_A_HXX
