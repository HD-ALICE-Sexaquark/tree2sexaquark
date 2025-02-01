#ifndef T2S_SEXAQUARK_CHANNEL_D_HXX
#define T2S_SEXAQUARK_CHANNEL_D_HXX

#include "Math/Vector4D.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Sexaquark {

/*
 *
 */
class ChannelD : public Particle::Sexaquark {
   public:
    inline Double_t DCAKaSV() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAKaLa() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfLambda)); };

    typedef Double_t (ChannelD::*MemFn)();

   private:
    ROOT::Math::PxPyPzEVector lvKaon;
    KFParticle kfKaon;
};

}  // namespace Sexaquark
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_D_HXX
