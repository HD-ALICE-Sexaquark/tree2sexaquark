#include "Cuts/Inspector.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*
 *
 */
void Inspector::Init() {
    //
    Cut eta([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.Eta(); },  //
            0.9, Cut::kAbsoluteMax);
    Cut mass([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.Mass(); },  //
             1.08, 1.16);
    Cut radius([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.Radius(); },  //
               20., Cut::kMinimum);
    Cut dist_from_pv([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.DistFromPV(); },  //
                     40., Cut::kMinimum);
    Cut cpa_wrt_pv([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.CPAwrtPV(); },  //
                   0.1, 0.99);
    Cut dca_wrt_pv([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.DCAwrtPV(); },  //
                   4., Cut::kMinimum);
    Cut dca_btw_dau([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.DCAbtwDau(); },  //
                    2., Cut::kMaximum);
    Cut dca_neg_v0([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.DCAnegV0(); },  //
                   2., Cut::kMaximum);
    Cut dca_pos_v0([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.DCAposV0(); },  //
                   2., Cut::kMaximum);
    Cut arm_qt_over_alpha([](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.ArmQtOverAlpha(); },  //
                          0.2, Cut::kMaximum);
    //
    AddCut(eta);
    AddCut(mass);
    AddCut(radius);
    AddCut(dist_from_pv);
    AddCut(cpa_wrt_pv);
    AddCut(dca_wrt_pv);
    AddCut(dca_btw_dau);
    AddCut(dca_neg_v0);
    AddCut(dca_pos_v0);
    AddCut(arm_qt_over_alpha);
}

/*
 *
 */
Bool_t Inspector::Approve(Particle::V0& thisV0) {
    //
    for (Cut& cut : fCutsCollection) {
        if (!cut.Check(thisV0)) return kFALSE;
    }
    InfoF("V0 approved. Glory to Arstotzka. %s", "");
    return kTRUE;
}

}  // namespace Cuts
}  // namespace Tree2Sexaquark
