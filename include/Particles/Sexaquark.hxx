#ifndef T2S_SEXAQUARK_HXX
#define T2S_SEXAQUARK_HXX

#include "Math/Vector4D.h"

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Particle {

/*
 *
 */
class Sexaquark : public Base {
   public:
    inline Double_t DCALaSV() { return TMath::Abs((Double_t)kfLambda.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCALaNegSV() { return TMath::Abs((Double_t)kfLambdaNeg.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCALaPosSV() { return TMath::Abs((Double_t)kfLambdaPos.GetDistanceFromVertex(kfThis)); };

   protected:
    KFParticle kfLambda;

   private:
    ROOT::Math::PxPyPzEVector lvLambda;
    KFParticle kfLambdaNeg;
    KFParticle kfLambdaPos;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_HXX
