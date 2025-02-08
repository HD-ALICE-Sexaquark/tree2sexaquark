#ifndef T2S_KAON_PAIR_HXX
#define T2S_KAON_PAIR_HXX

#include "Math/VectorUtil.h"
#include "TMath.h"

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class KaonPair : public Base {
   public:
    /* Setters */
    void SetKaonPairInfo();
    void SetKinematics(ROOT::Math::PxPyPzEVector lv_kk, ROOT::Math::PxPyPzEVector lv_ka, ROOT::Math::PxPyPzEVector lv_kb) {
        lvThis = lv_kk;
        lvKaonA = lv_ka;
        lvKaonB = lv_kb;
    }
    void SetGeometry(KFParticle kf_kk, KFParticle kf_ka, KFParticle kf_kb) {
        kfThis = kf_kk;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());

        kfKaonA = kf_ka;
        kfKaonB = kf_kb;
    }
    void SetTrueInfo(Bool_t is_signal, Int_t reaction_id, Bool_t is_hybrid, Int_t pdg_noncomb_bkg);

    inline Double_t OpeningAngle() { return ROOT::Math::VectorUtil::Angle(lvKaonA.Vect(), lvKaonB.Vect()); }
    inline Double_t DCAbtwKK() { return TMath::Abs((Double_t)kfKaonA.GetDistanceFromParticle(kfKaonB)); };
    inline Double_t DCAkaSV() { return TMath::Abs((Double_t)kfKaonA.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAkbSV() { return TMath::Abs((Double_t)kfKaonB.GetDistanceFromVertex(kfThis)); };

    UInt_t KaonA_EsdIdx;
    UInt_t KaonB_EsdIdx;

    /* True Information */
    Int_t NonCombBkg_PdgCode;

   private:
    ROOT::Math::PxPyPzEVector lvKaonA;
    ROOT::Math::PxPyPzEVector lvKaonB;

    KFParticle kfKaonA;
    KFParticle kfKaonB;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_KAON_PAIR_HXX
