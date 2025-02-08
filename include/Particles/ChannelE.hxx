#ifndef T2S_SEXAQUARK_CHANNEL_E_HXX
#define T2S_SEXAQUARK_CHANNEL_E_HXX

#include "Math/Vector4D.h"

#include "Particles/Sexaquark.hxx"

namespace Tree2Sexaquark {
namespace Candidate {

class ChannelE : public Sexaquark {
   public:
    /* Setters */
    void SetSexaquarkInfo();
    void SetKinematics(ROOT::Math::PxPyPzEVector lv_sexaquark, ROOT::Math::PxPyPzEVector lv_lambda, ROOT::Math::PxPyPzEVector lv_kaon,
                       ROOT::Math::PxPyPzEVector lv_piminus, ROOT::Math::PxPyPzEVector lv_piplus) {
        lvThis = lv_sexaquark;
        lvLambda = lv_lambda;
        lvKaon = lv_kaon;
        lvPiMinus = lv_piminus;
        lvPiPlus = lv_piplus;
        lvThis_asDecay = lvLambda + lvKaon + lvPiMinus + lvPiPlus;
    }
    void SetGeometry(KFParticle kf_sexaquark, KFParticle kf_lambda, KFParticle kf_lambda_neg, KFParticle kf_lambda_pos, KFParticle kf_kaon,
                     KFParticle kf_piminus, KFParticle kf_piplus) {
        kfThis = kf_sexaquark;
        v3This = ROOT::Math::XYZPoint(kfThis.GetX(), kfThis.GetY(), kfThis.GetZ());

        kfLambda = kf_lambda;
        kfLambda_Neg = kf_lambda_neg;
        kfLambda_Pos = kf_lambda_pos;

        kfKaon = kf_kaon;
        kfPiMinus = kf_piminus;
        kfPiPlus = kf_piplus;
    }
    void SetTrueInfo(Bool_t is_signal, Int_t reaction_id, Bool_t is_hybrid, Int_t pdg_noncomb_bkg);

    inline Double_t DCAKaSV() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfThis)); };
    inline Double_t DCAKaLa() { return TMath::Abs((Double_t)kfKaon.GetDistanceFromVertex(kfLambda)); };
    inline Double_t DCApmSV() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCAppSV() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfThis)); }
    inline Double_t DCApmLa() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfLambda)); }
    inline Double_t DCAppLa() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfLambda)); }
    inline Double_t DCApmKa() { return TMath::Abs((Double_t)kfPiMinus.GetDistanceFromVertex(kfKaon)); }
    inline Double_t DCAppKa() { return TMath::Abs((Double_t)kfPiPlus.GetDistanceFromVertex(kfKaon)); }

    UInt_t Kaon_EsdIdx;
    UInt_t PionPair_Idx;
    UInt_t PiMinus_EsdIdx;
    UInt_t PiPlus_EsdIdx;

   private:
    ROOT::Math::PxPyPzEVector lvKaon;
    ROOT::Math::PxPyPzEVector lvPiMinus;
    ROOT::Math::PxPyPzEVector lvPiPlus;
    KFParticle kfKaon;
    KFParticle kfPiMinus;
    KFParticle kfPiPlus;
};

}  // namespace Candidate
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_CHANNEL_E_HXX
