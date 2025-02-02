#ifndef T2S_SEXAQUARK_CHANNEL_D_HXX
#define T2S_SEXAQUARK_CHANNEL_D_HXX

#include "Math/Vector4D.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

/*
 *
 */
class ChannelD : public Sexaquark {
   public:
    inline Double_t DCAKaSV() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAKaLa() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfLambda)); };

   private:
    ROOT::Math::PxPyPzEVector lvKaon;
    KFParticle kfKaon;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_D_HXX
