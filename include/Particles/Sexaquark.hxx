#ifndef T2S_SEXAQUARK_HXX
#define T2S_SEXAQUARK_HXX

#include "Math/Vector4D.h"

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class Sexaquark : public Base {
   public:
    enum class Channel : Char_t { A = 'A', D = 'D', E = 'E', H = 'H' };

    inline Double_t E_asDecay() { return lvThis_asDecay.E(); }
    inline Double_t DCALaSV() { return TMath::Abs((Double_t)kfLambda.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCALaNegSV() { return TMath::Abs((Double_t)kfLambda_Neg.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCALaPosSV() { return TMath::Abs((Double_t)kfLambda_Pos.GetDistanceFromVertex(kfThis)); };
    inline Double_t Lambda_DecayLength() { return (Double_t)kfLambda.GetDecayLength(); }

    UInt_t Lambda_Idx;
    UInt_t Lambda_Neg_EsdIdx;
    UInt_t Lambda_Pos_EsdIdx;

    /* True Information */
    Int_t NonCombBkg_PdgCode;

   protected:
    ROOT::Math::PxPyPzEVector lvThis_asDecay;

    ROOT::Math::PxPyPzEVector lvLambda;
    KFParticle kfLambda;
    KFParticle kfLambda_Neg;
    KFParticle kfLambda_Pos;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_HXX
