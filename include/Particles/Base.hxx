#ifndef T2S_BASE_PARTICLE_HXX
#define T2S_BASE_PARTICLE_HXX

#include "Math/Vector3Dfwd.h"
#include "Rtypes.h"
#include "TMath.h"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "Math/Common.hxx"

namespace Tree2Sexaquark {
namespace Particle {

/*
 * Base class for all particles
 */
class Base {
   public:
    Base(ROOT::Math::PxPyPzEVector lvParticle, KFParticle kfParticle, KFVertex kfPV)
        : lvParticle(lvParticle),  //
          kfParticle(kfParticle),
          kfPV(kfPV) {
        //
        v3Particle = ROOT::Math::XYZPoint(kfParticle.GetX(), kfParticle.GetY(), kfParticle.GetZ());
        v3PV = ROOT::Math::XYZPoint(kfPV.GetX(), kfPV.GetY(), kfPV.GetZ());
    }
    ~Base() = default;

    inline Double_t Px() { return lvParticle.Px(); }
    inline Double_t Py() { return lvParticle.Py(); }
    inline Double_t Pz() { return lvParticle.Pz(); }
    inline Double_t E() { return lvParticle.E(); }
    inline Double_t Eta() { return lvParticle.Eta(); }
    inline Double_t Pt() { return lvParticle.Pt(); }
    inline Double_t Mass() { return lvParticle.M(); }

    inline Float_t Xv() { return kfParticle.GetX(); }
    inline Float_t Yv() { return kfParticle.GetY(); }
    inline Float_t Zv() { return kfParticle.GetZ(); }

    inline Double_t Radius() { return v3Particle.Rho(); }
    inline Double_t DistFromPV() { return (v3Particle - v3PV).R(); }
    inline Double_t CPAwrtPV() { return Math::CosinePointingAngle(GetMomentumVector(), v3Particle, v3PV); }
    inline Double_t DCAwrtPV() { return TMath::Abs(kfParticle.GetDistanceFromVertex(kfPV)); }

    ROOT::Math::XYZVector GetMomentumVector() { return lvParticle.Vect(); }
    KFParticle GetKFParticle() { return kfParticle; }

   private:
    ROOT::Math::PxPyPzEVector lvParticle;
    KFParticle kfParticle;
    KFVertex kfPV;

    ROOT::Math::XYZPoint v3Particle;
    ROOT::Math::XYZPoint v3PV;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_BASE_PARTICLE_HXX
