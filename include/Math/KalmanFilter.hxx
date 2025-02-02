#ifndef T2S_KALMAN_FILTER_HXX
#define T2S_KALMAN_FILTER_HXX

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Trees/Input.hxx"

namespace Tree2Sexaquark {
namespace Math {

KFParticle CreateKFParticle(Track_tt track, Double_t mass);
KFVertex CreateKFVertex(Float_t* xyz, Float_t* cov_matrix);
KFParticle TransportKFParticle(KFParticle kf_this, KFParticle kf_other, Double_t mass_this, Int_t charge_this);

/*
class KFParticleMother : public KFParticle {
   public:
    Bool_t CheckDaughter(KFParticle daughter) {
        Float_t m[8], mV[36], D[3][3];
        return KFParticleBase::GetMeasurement(daughter, m, mV, D);
    }
};
*/

}  // namespace Math
}  // namespace Tree2Sexaquark

#endif  // T2S_KALMAN_FILTER_HXX
