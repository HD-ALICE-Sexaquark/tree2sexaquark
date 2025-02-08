#ifndef T2S_SEXAQUARK_CHANNEL_A_HXX
#define T2S_SEXAQUARK_CHANNEL_A_HXX

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class ChannelA : public Sexaquark {
   public:
    /* Setters */
    void SetSexaquarkInfo();
    void SetKinematics(ROOT::Math::PxPyPzEVector lv_sexaquark, ROOT::Math::PxPyPzEVector lv_lambda, ROOT::Math::PxPyPzEVector lv_kaonzero) {
        lvThis = lv_sexaquark;
        lvLambda = lv_lambda;
        lvKaonZeroShort = lv_kaonzero;
        lvThis_asDecay = lvLambda + lvKaonZeroShort;
    }
    void SetGeometry(KFParticle kf_sexaquark, KFParticle kf_lambda, KFParticle kf_lambda_neg, KFParticle kf_lambda_pos, KFParticle kf_kaonzero,
                     KFParticle kf_kaonzero_neg, KFParticle kf_kaonzero_pos) {
        kfThis = kf_sexaquark;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());

        kfLambda = kf_lambda;
        kfLambda_Neg = kf_lambda_neg;
        kfLambda_Pos = kf_lambda_pos;

        kfKaonZeroShort = kf_kaonzero;
        kfKaonZeroShort_Neg = kf_kaonzero_neg;
        kfKaonZeroShort_Pos = kf_kaonzero_pos;
    }
    void SetTrueInfo(Bool_t is_signal, Int_t reaction_id, Bool_t is_hybrid, Int_t pdg_noncomb_bkg);

    inline Double_t OpeningAngle() { return ROOT::Math::VectorUtil::Angle(lvLambda.Vect(), lvKaonZeroShort.Vect()); }
    inline Double_t DCAbtwV0s() { return TMath::Abs((Double_t)kfLambda.GetDistanceFromParticle(kfKaonZeroShort)); };
    inline Double_t DCAK0SV() { return TMath::Abs((Double_t)kfKaonZeroShort.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAK0NegSV() { return TMath::Abs((Double_t)kfKaonZeroShort_Neg.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAK0PosSV() { return TMath::Abs((Double_t)kfKaonZeroShort_Pos.GetDistanceFromVertex(kfThis)); };
    inline Double_t K0S_DecayLength() { return (Double_t)kfKaonZeroShort.GetDecayLength(); }

    UInt_t K0S_Idx;
    UInt_t K0S_Neg_EsdIdx;
    UInt_t K0S_Pos_EsdIdx;

   private:
    ROOT::Math::PxPyPzEVector lvKaonZeroShort;
    KFParticle kfKaonZeroShort;
    KFParticle kfKaonZeroShort_Neg;
    KFParticle kfKaonZeroShort_Pos;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_A_HXX
