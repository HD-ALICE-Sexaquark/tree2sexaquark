#ifndef T2S_SEXAQUARK_CHANNEL_D_HXX
#define T2S_SEXAQUARK_CHANNEL_D_HXX

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class ChannelD : public Sexaquark {
   public:
    /* Setters */
    void SetSexaquarkInfo();
    void SetKinematics(ROOT::Math::PxPyPzEVector lv_sexaquark, ROOT::Math::PxPyPzEVector lv_lambda, ROOT::Math::PxPyPzEVector lv_kaon) {
        lvThis = lv_sexaquark;
        lvLambda = lv_lambda;
        lvKaon = lv_kaon;
        lvThis_asDecay = lvLambda + lvKaon;
    }
    void SetGeometry(KFParticle kf_sexaquark, KFParticle kf_lambda, KFParticle kf_lambda_neg, KFParticle kf_lambda_pos, KFParticle kf_kaon) {
        kfThis = kf_sexaquark;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());

        kfLambda = kf_lambda;
        kfLambda_Neg = kf_lambda_neg;
        kfLambda_Pos = kf_lambda_pos;

        kfKaon = kf_kaon;
    }
    void SetTrueInfo(Bool_t is_signal, Int_t reaction_id, Bool_t is_hybrid, Int_t pdg_noncomb_bkg);

    inline Double_t OpeningAngle() { return ROOT::Math::VectorUtil::Angle(lvLambda.Vect(), lvKaon.Vect()); }
    inline Double_t DCAKaSV() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAKaLa() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfLambda)); };
    inline Double_t DCALaNegKa() { return TMath::Abs((Double_t)kfLambda_Neg.GetDistanceFromParticle(kfKaon)); };
    inline Double_t DCALaPosKa() { return TMath::Abs((Double_t)kfLambda_Pos.GetDistanceFromParticle(kfKaon)); };

    UInt_t Kaon_EsdIdx;

   private:
    ROOT::Math::PxPyPzEVector lvKaon;
    KFParticle kfKaon;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_D_HXX
