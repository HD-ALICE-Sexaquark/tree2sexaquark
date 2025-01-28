#include "Rtypes.h"

#include "Math/Vector4D.h"

namespace Tree2Sexaquark {
namespace Math {

Double_t CosinePointingAngle(ROOT::Math::PxPyPzEVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                             Double_t refPointZ);
Double_t CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                             Double_t refPointZ);
Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                         Double_t Pos_Py, Double_t Pos_Pz);
Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz);
Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z, Double_t refPointX,
                      Double_t refPointY, Double_t refPointZ);

}  // namespace Math
}  // namespace Tree2Sexaquark
