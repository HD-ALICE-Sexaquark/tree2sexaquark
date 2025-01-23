#ifndef T2S_KALMAN_FILTER_HXX
#define T2S_KALMAN_FILTER_HXX

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Trees/Structure.hxx"

// namespace Tree2Sexaquark::Math {

KFParticle CreateKFParticle(Track_tt track, Double_t mass);
// KFVertex CreateKFVertex(Double_t* XYZ, Double_t* CovarianceMatrix);
KFVertex CreateKFVertex(Double_t* XYZ);
KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Double_t massThis, Int_t chargeThis);

// }  // namespace Tree2Sexaquark::Math
/*
class KFParticleMother : public KFParticle {
   public:
    Bool_t CheckDaughter(KFParticle daughter) {
        Float_t m[8], mV[36], D[3][3];
        return KFParticleBase::GetMeasurement(daughter, m, mV, D);
    }
};
*/
#endif  // T2S_KALMAN_FILTER_HXX
