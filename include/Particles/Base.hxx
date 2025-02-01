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
        : lvThis(lvParticle),  //
          kfThis(kfParticle),
          kfPV(kfPV) {
        //
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());
        v3PV = ROOT::Math::XYZPoint(kfPV.GetX(), kfPV.GetY(), kfPV.GetZ());
    }
    ~Base() = default;

    inline Double_t Px() { return lvThis.Px(); }
    inline Double_t Py() { return lvThis.Py(); }
    inline Double_t Pz() { return lvThis.Pz(); }
    inline Double_t E() { return lvThis.E(); }
    inline Double_t Eta() { return lvThis.Eta(); }
    inline Double_t Pt() { return lvThis.Pt(); }
    inline Double_t Mass() { return lvThis.M(); }
    inline Double_t Rapidity() { return lvThis.Rapidity(); }

    inline Float_t Xv() { return kfThis.GetX(); }
    inline Float_t Yv() { return kfThis.GetY(); }
    inline Float_t Zv() { return kfThis.GetZ(); }

    inline Double_t Radius() { return v3This.Rho(); }
    inline Double_t DistFromPV() { return (v3This - v3PV).R(); }
    inline Double_t CPAwrtPV() { return Math::CosinePointingAngle(GetMomentumVector(), v3This, v3PV); }
    inline Double_t DCAwrtPV() { return TMath::Abs(kfThis.GetDistanceFromVertex(kfPV)); }

    ROOT::Math::XYZVector GetMomentumVector() { return lvThis.Vect(); }

   protected:
    KFParticle kfThis;

   private:
    ROOT::Math::PxPyPzEVector lvThis;
    KFVertex kfPV;

    ROOT::Math::XYZPoint v3This;
    ROOT::Math::XYZPoint v3PV;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_BASE_PARTICLE_HXX
