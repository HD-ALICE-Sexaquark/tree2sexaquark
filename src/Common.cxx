#include "Math/Common.hxx"

#include "TMath.h"

#include "Math/VectorUtil.h"

namespace Tree2Sexaquark {
namespace Math {

/*
 * Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 */
Double_t CosinePointingAngle(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZPoint pos_v0, ROOT::Math::XYZPoint pos_ref) {
    //
    return TMath::Cos(ROOT::Math::VectorUtil::Angle(mom_v0, pos_v0 - pos_ref));
}

/*
 * Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particles' 4-momentum.
 */
Double_t CosinePointingAngle(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t v0_x, Double_t v0_y, Double_t v0_z, Double_t ref_x,
                             Double_t ref_y, Double_t ref_z) {
    //
    ROOT::Math::XYZVector MomV0(v0_px, v0_py, v0_pz);
    ROOT::Math::XYZPoint PosV0(v0_x, v0_y, v0_z);
    ROOT::Math::XYZPoint PosRef(ref_x, ref_y, ref_z);
    return CosinePointingAngle(MomV0, PosV0, PosRef);
}

/*
 * Calculate Armenteros-Podolanski qT
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::PtArmV0())
 */
Double_t ArmenterosQt(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZVector mom_dau) {
    //
    return ROOT::Math::VectorUtil::Perp(mom_v0, mom_dau);
}

/*
 * Overload of `ArmenterosQt(...)`
 */
Double_t ArmenterosQt(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t dau_px, Double_t dau_py, Double_t dau_pz) {
    //
    ROOT::Math::XYZVector MomV0(v0_px, v0_py, v0_pz);
    ROOT::Math::XYZVector MomDau(dau_px, dau_py, dau_pz);
    return ArmenterosQt(MomV0, MomDau);
}

/*
 * Calculate Armenteros-Podolanski alpha
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::AlphaV0())
 */
Double_t ArmenterosAlpha(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZVector mom_neg, ROOT::Math::XYZVector mom_pos) {
    //
    Double_t lQlNeg = mom_neg.Dot(mom_v0) / mom_v0.R();
    Double_t lQlPos = mom_pos.Dot(mom_v0) / mom_v0.R();
    if (TMath::Abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 * Overload of `ArmenterosAlpha()`
 */
Double_t ArmenterosAlpha(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t neg_px, Double_t neg_py, Double_t neg_pz, Double_t pos_px,
                         Double_t pos_py, Double_t pos_pz) {
    //
    ROOT::Math::XYZVector MomV0(v0_px, v0_py, v0_pz);
    ROOT::Math::XYZVector MomNeg(neg_px, neg_py, neg_pz);
    ROOT::Math::XYZVector MomPos(pos_px, pos_py, pos_pz);
    return ArmenterosAlpha(MomV0, MomNeg, MomPos);
}

/*
 * Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 * This function stores the point of closest approach, and returns its distance to the PV.
 */
Double_t LinePointDCA(ROOT::Math::XYZVector mom_v0, ROOT::Math::XYZPoint pos_v0, ROOT::Math::XYZPoint pos_ref) {
    //
    ROOT::Math::XYZVector CrossProduct = (pos_ref - pos_v0).Cross(mom_v0);
    return CrossProduct.R() / mom_v0.R();
}

/*
 * Overload of `LinePointDCA()`
 */
Double_t LinePointDCA(Double_t v0_px, Double_t v0_py, Double_t v0_pz, Double_t v0_x, Double_t v0_y, Double_t v0_z, Double_t ref_x, Double_t ref_y,
                      Double_t ref_z) {
    //
    ROOT::Math::XYZVector MomV0(v0_px, v0_py, v0_pz);
    ROOT::Math::XYZPoint PosV0(v0_x, v0_y, v0_z);
    ROOT::Math::XYZPoint PosRef(ref_x, ref_y, ref_z);
    return LinePointDCA(MomV0, PosV0, PosRef);
}

}  // namespace Math
}  // namespace Tree2Sexaquark
