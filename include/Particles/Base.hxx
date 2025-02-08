#ifndef T2S_BASE_PARTICLE_HXX
#define T2S_BASE_PARTICLE_HXX

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
namespace Candidate {

/*
 * Base class for all particles
 */
class Base {
   public:
    Base() = default;
    ~Base() = default;

    virtual void SetKinematics(ROOT::Math::PxPyPzEVector lorentz_vector) { lvThis = lorentz_vector; }
    virtual void SetGeometry(KFParticle kf_particle) {
        kfThis = kf_particle;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());
    }

    /*
     * PENDING: I should move this one
     */
    void SetPrimaryVertex(KFVertex kf_vertex) {
        kfPV = kf_vertex;
        v3PV = ROOT::Math::XYZPoint(kfPV.GetX(), kfPV.GetY(), kfPV.GetZ());
    }

    KFParticle GetKf() { return kfThis; }

    inline Double_t Px() { return lvThis.Px(); }
    inline Double_t Py() { return lvThis.Py(); }
    inline Double_t Pz() { return lvThis.Pz(); }
    inline Double_t E() { return lvThis.E(); }
    inline Double_t Eta() { return lvThis.Eta(); }
    inline Double_t Pt() { return lvThis.Pt(); }
    inline Double_t Mass() { return lvThis.M(); }
    inline Double_t Rapidity() { return lvThis.Rapidity(); }

    inline Double_t Xv() { return (Double_t)kfThis.GetX(); }
    inline Double_t Yv() { return (Double_t)kfThis.GetY(); }
    inline Double_t Zv() { return (Double_t)kfThis.GetZ(); }
    inline Double_t Chi2() { return (Double_t)kfThis.GetChi2(); }
    inline Double_t Chi2ndf() { return Chi2() / (Double_t)kfThis.GetNDF(); }

    inline Double_t Radius() { return v3This.Rho(); }
    inline Double_t DistFromPV() { return (v3This - v3PV).R(); }
    inline Double_t CPAwrtPV() { return Math::CosinePointingAngle(lvThis.Vect(), v3This, v3PV); }
    inline Double_t DCAwrtPV() { return TMath::Abs(kfThis.GetDistanceFromVertex(kfPV)); }

    /* True Information */
    Bool_t IsSignal;
    Int_t ReactionID;
    Bool_t IsHybrid;

   protected:
    ROOT::Math::PxPyPzEVector lvThis;

    KFParticle kfThis;
    ROOT::Math::XYZPoint v3This;

    KFVertex kfPV;
    ROOT::Math::XYZPoint v3PV;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_BASE_PARTICLE_HXX
