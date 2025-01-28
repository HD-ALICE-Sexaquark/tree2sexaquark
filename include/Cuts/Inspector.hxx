#ifndef T2S_CUTS_INSPECTOR_HXX
#define T2S_CUTS_INSPECTOR_HXX

#include <vector>

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"

#include "Utilities/Logger.hxx"

#include "Cuts/Cut.hxx"
// #include "Cuts/Default.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

class Inspector {
   public:
    Inspector() = default;
    ~Inspector() = default;

    void Init();
    inline void AddCut(Cut cut) { fCuts.push_back(cut); }
    Bool_t Approve(Particle::V0& thisV0);

    /*
    Bool_t V0::PassesSelection() {
        //
        ROOT::Math::XYZPoint V0Vertex(GetX(), GetY(), GetZ());

        if (kMin_V0_Radius.count(PdgCode) && V0Vertex.Rho() < kMin_V0_Radius[PdgCode]) return kFALSE;

        if (kMin_V0_Mass.count(PdgCode) && lvV0.M() < kMin_V0_Mass[PdgCode]) return kFALSE;
        if (kMax_V0_Mass.count(PdgCode) && lvV0.M() > kMax_V0_Mass[PdgCode]) return kFALSE;

        if (kMin_V0_Pt.count(PdgCode) && lvV0.Pt() < kMin_V0_Pt[PdgCode]) return kFALSE;

        if (kMax_V0_Eta.count(PdgCode) && TMath::Abs(lvV0.Eta()) > kMax_V0_Eta[PdgCode]) return kFALSE;

        if (kMin_V0_DistFromPV.count(PdgCode) && (V0Vertex - v3PrimaryVertex).R() < kMin_V0_DistFromPV[PdgCode]) return kFALSE;
        if (kMax_V0_DistFromPV.count(PdgCode) && (V0Vertex - v3PrimaryVertex).R() > kMax_V0_DistFromPV[PdgCode]) return kFALSE;

        Double_t CPAwrtPV =
            Math::CosinePointingAngle(lvV0, GetX(), GetY(), GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        if (kMin_V0_CPAwrtPV.count(PdgCode) && CPAwrtPV < kMin_V0_CPAwrtPV[PdgCode]) return kFALSE;
        if (kMax_V0_CPAwrtPV.count(PdgCode) && CPAwrtPV > kMax_V0_CPAwrtPV[PdgCode]) return kFALSE;

        Double_t DCAwrtPV = Math::LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), GetX(), GetY(), GetZ(), fPrimaryVertex->GetX(),
    fPrimaryVertex->GetY(), fPrimaryVertex->GetZ()); if (kMin_V0_DCAwrtPV.count(PdgCode) && DCAwrtPV < kMin_V0_DCAwrtPV[PdgCode]) return kFALSE;

        Double_t DCAbtwDau = TMath::Abs(kfNegDaughter.GetDistanceFromParticle(kfPosDaughter));
        if (kMax_V0_DCAbtwDau.count(PdgCode) && DCAbtwDau > kMax_V0_DCAbtwDau[PdgCode]) return kFALSE;

        Double_t DCAnegV0 = TMath::Abs(kfNegDaughter.GetDistanceFromVertex(kfV0));
        if (kMax_V0_DCAnegV0.count(PdgCode) && DCAnegV0 > kMax_V0_DCAnegV0[PdgCode]) return kFALSE;

        Double_t DCAposV0 = TMath::Abs(kfPosDaughter.GetDistanceFromVertex(kfV0));
        if (kMax_V0_DCAposV0.count(PdgCode) && DCAposV0 > kMax_V0_DCAposV0[PdgCode]) return kFALSE;

        Double_t ArmPt = Math::ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvNegDaughter.Px(), lvNegDaughter.Py(), lvNegDaughter.Pz());
        Double_t ArmAlpha = Math::ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvNegDaughter.Px(), lvNegDaughter.Py(), lvNegDaughter.Pz(),
                                                  lvPosDaughter.Px(), lvPosDaughter.Py(), lvPosDaughter.Pz());
        Double_t ArmPtOverAlpha = TMath::Abs(ArmPt / ArmAlpha);
        if (kMax_V0_ArmPtOverAlpha.count(PdgCode) && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[PdgCode]) return kFALSE;

        Double_t Chi2ndf = (Double_t)GetChi2() / (Double_t)GetNDF();
        if (kMax_V0_Chi2ndf.count(PdgCode) && Chi2ndf > kMax_V0_Chi2ndf[PdgCode]) return kFALSE;

        return kTRUE;
    }
    */

   private:
    std::vector<Cut> fCuts;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_INSPECTOR_HXX
