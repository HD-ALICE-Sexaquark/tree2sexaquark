#include "Rtypes.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

namespace Tree2Sexaquark {
namespace Math {

Double_t CosinePointingAngle(ROOT::Math::XYZVector momV0, ROOT::Math::XYZPoint posV0, ROOT::Math::XYZPoint posRef);
Double_t CosinePointingAngle(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_Xv, Double_t V0_Yv, Double_t V0_Zv, Double_t Ref_Xv,
                             Double_t Ref_Yv, Double_t Ref_Zv);

Double_t ArmenterosQt(ROOT::Math::XYZVector momV0, ROOT::Math::XYZVector momDau);
Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Dau_Px, Double_t Dau_Py, Double_t Dau_Pz);

Double_t ArmenterosAlpha(ROOT::Math::XYZVector momV0, ROOT::Math::XYZVector momNeg, ROOT::Math::XYZVector momPos);
Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                         Double_t Pos_Py, Double_t Pos_Pz);

Double_t LinePointDCA(ROOT::Math::XYZVector momV0, ROOT::Math::XYZPoint posV0, ROOT::Math::XYZPoint posRef);
Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_Xv, Double_t V0_Yv, Double_t V0_Zv, Double_t Ref_Xv,
                      Double_t Ref_Yv, Double_t Ref_Zv);

}  // namespace Math
}  // namespace Tree2Sexaquark
