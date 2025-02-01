#ifndef T2S_KAON_PAIR_HXX
#define T2S_KAON_PAIR_HXX

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Sexaquark {

/*
 *
 */
class KaonPair : public Particle::Base {
   public:
    inline Double_t DCAbtwKK() { return TMath::Abs((Double_t)kfKaonA.GetDistanceFromParticle(kfKaonB)); };
    inline Double_t DCAkaSV() { return TMath::Abs((Double_t)kfKaonA.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAkbSV() { return TMath::Abs((Double_t)kfKaonB.GetDistanceFromVertex(kfThis)); };

    typedef Double_t (KaonPair::*MemFn)();

   private:
    ROOT::Math::PxPyPzEVector lvKaonA;
    ROOT::Math::PxPyPzEVector lvKaonB;

    KFParticle kfKaonA;
    KFParticle kfKaonB;
};

}  // namespace Sexaquark
}  // namespace Tree2Sexaquark

#endif  // T2S_KAON_PAIR_HXX
