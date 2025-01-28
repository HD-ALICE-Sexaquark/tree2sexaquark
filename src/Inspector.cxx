#include "Cuts/Inspector.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

void Inspector::Init() {
    Cut eta_cut;
    auto eta_exp = [](Particle::V0 thisV0) -> Float_t { return (Float_t)thisV0.Eta(); };
    eta_cut.SetExpression(eta_exp);
    eta_cut.SetMaximum(0.8);
    eta_cut.SetMinimum(-0.8);
    AddCut(eta_cut);
}

Bool_t Inspector::Approve(Particle::V0& thisV0) {
    for (Cut& cut : fCuts) {
        if (!cut.Check(thisV0)) return kFALSE;
    }
    InfoF("V0 approved. Glory to Arstotzka. %s", "");
    return kTRUE;
}

}  // namespace Cuts
}  // namespace Tree2Sexaquark
