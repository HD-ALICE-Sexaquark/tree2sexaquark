#ifndef T2S_SEXAQUARK_CHANNEL_E_HXX
#define T2S_SEXAQUARK_CHANNEL_E_HXX

#include "Math/Vector4D.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

/*
 *
 */
class ChannelE : public Sexaquark {
   public:
    inline Double_t DCAKaSV() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAKaLa() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfLambda)); };
    inline Double_t DCApmSV() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCAppSV() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCApmLa() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfLambda)); }
    inline Double_t DCAppLa() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfLambda)); }
    inline Double_t DCApmKa() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfKaon)); }
    inline Double_t DCAppKa() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfKaon)); }

   private:
    ROOT::Math::PxPyPzEVector lvKaon;
    KFParticle kfKaon;
    KFParticle kfPiMinus;
    KFParticle kfPiPlus;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_E_HXX
