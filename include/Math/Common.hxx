#include "Rtypes.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

namespace Tree2Sexaquark {
namespace Math {

Double_t CosinePointingAngle(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZPoint pos_v0, ROOT::Math::XYZPoint pos_ref);
Double_t CosinePointingAngle(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t v0_x, Double_t v0_y, Double_t v0_z, Double_t ref_x,
                             Double_t ref_y, Double_t ref_z);

Double_t ArmenterosQt(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZVector mom_dau);
Double_t ArmenterosQt(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t dau_px, Double_t dau_py, Double_t dau_pz);

Double_t ArmenterosAlpha(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZVector mom_neg, ROOT::Math::XYZVector mom_pos);
Double_t ArmenterosAlpha(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t neg_px, Double_t neg_py, Double_t neg_pz, Double_t pos_px,
                         Double_t pos_py, Double_t pos_pz);

Double_t LinePointDCA(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZPoint pos_v0, ROOT::Math::XYZPoint pos_ref);
Double_t LinePointDCA(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t v0_x, Double_t v0_y, Double_t v0_z, Double_t ref_x, Double_t ref_y,
                      Double_t ref_z);

}  // namespace Math
}  // namespace Tree2Sexaquark
