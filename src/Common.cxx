#include "Math/Common.hxx"

#include "TMath.h"

#include "Math/VectorUtil.h"

namespace Tree2Sexaquark {
namespace Math {

/*
 * Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 */
Double_t CosinePointingAngle(ROOT::Math::XYZVector momV0, ROOT::Math::XYZPoint posV0, ROOT::Math::XYZPoint posRef) {
    //
    return TMath::Cos(ROOT::Math::VectorUtil::Angle(momV0, posV0 - posRef));
}

/*
 * Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particles' 4-momentum.
 */
Double_t CosinePointingAngle(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_Xv, Double_t V0_Yv, Double_t V0_Zv, Double_t Ref_Xv,
                             Double_t Ref_Yv, Double_t Ref_Zv) {
    //
    ROOT::Math::XYZVector momV0(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZPoint posV0(V0_Xv, V0_Yv, V0_Zv);
    ROOT::Math::XYZPoint posRef(Ref_Xv, Ref_Yv, Ref_Zv);
    return CosinePointingAngle(momV0, posV0, posRef);
}

/*
 * Calculate Armenteros-Podolanski qT
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::PtArmV0())
 */
Double_t ArmenterosQt(ROOT::Math::XYZVector momV0, ROOT::Math::XYZVector momDau) {
    //
    return ROOT::Math::VectorUtil::Perp(momV0, momDau);
}

/*
 * Overload of `ArmenterosQt(...)`
 */
Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Dau_Px, Double_t Dau_Py, Double_t Dau_Pz) {
    //
    ROOT::Math::XYZVector momV0(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZVector momDau(Dau_Px, Dau_Py, Dau_Pz);
    return ArmenterosQt(momV0, momDau);
}

/*
 * Calculate Armenteros-Podolanski alpha
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::AlphaV0())
 */
Double_t ArmenterosAlpha(ROOT::Math::XYZVector momV0, ROOT::Math::XYZVector momNeg, ROOT::Math::XYZVector momPos) {
    //
    Double_t lQlNeg = momNeg.Dot(momV0) / momV0.R();
    Double_t lQlPos = momPos.Dot(momV0) / momV0.R();
    if (TMath::Abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 * Overload of `ArmenterosAlpha()`
 */
Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                         Double_t Pos_Py, Double_t Pos_Pz) {
    //
    ROOT::Math::XYZVector momV0(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZVector momNeg(Neg_Px, Neg_Py, Neg_Pz);
    ROOT::Math::XYZVector momPos(Pos_Px, Pos_Py, Pos_Pz);
    return ArmenterosAlpha(momV0, momNeg, momPos);
}

/*
 * Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 * This function stores the point of closest approach, and returns its distance to the PV.
 */
Double_t LinePointDCA(ROOT::Math::XYZVector momV0, ROOT::Math::XYZPoint posV0, ROOT::Math::XYZPoint posRef) {
    //
    ROOT::Math::XYZVector CrossProduct = (posRef - posV0).Cross(momV0);
    return CrossProduct.R() / momV0.R();
}

/*
 * Overload of `LinePointDCA()`
 */
Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_Xv, Double_t V0_Yv, Double_t V0_Zv, Double_t Ref_Xv,
                      Double_t Ref_Yv, Double_t Ref_Zv) {
    //
    ROOT::Math::XYZVector momV0(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZPoint posV0(V0_Xv, V0_Yv, V0_Zv);
    ROOT::Math::XYZPoint posRef(Ref_Xv, Ref_Yv, Ref_Zv);
    return LinePointDCA(momV0, posV0, posRef);
}

}  // namespace Math
}  // namespace Tree2Sexaquark
