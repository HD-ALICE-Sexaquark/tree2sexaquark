#include "TMath.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"

#include "Math/Common.hxx"

namespace Tree2Sexaquark {
namespace Math {

/*
 * Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle())
 */
Double_t CosinePointingAngle(ROOT::Math::PxPyPzEVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                             Double_t refPointZ) {
    ROOT::Math::XYZVector posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(ROOT::Math::VectorUtil::Angle(lvParticle, posRelativeToRef));
}

/*
 * Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's 4-momentum.
 */
Double_t CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                             Double_t refPointZ) {
    ROOT::Math::XYZVector posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    ROOT::Math::XYZVector momParticle(Px, Py, Pz);
    return TMath::Cos(ROOT::Math::VectorUtil::Angle(momParticle, posRelativeToRef));
}

/*
 * This gives the Armenteros-Podolanski alpha.
 * (Based on `AliRoot/STEER/ESD/AliESDv0::AlphaV0()`)
 */
Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                         Double_t Pos_Py, Double_t Pos_Pz) {
    ROOT::Math::XYZVector momTot(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZVector momNeg(Neg_Px, Neg_Py, Neg_Pz);
    ROOT::Math::XYZVector momPos(Pos_Px, Pos_Py, Pos_Pz);
    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.R();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.R();
    if (TMath::Abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 * This gives the Armenteros-Podolanski qT
 * Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::PtArmV0())
 */
Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
    ROOT::Math::XYZVector momTot(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZVector momNeg(Neg_Px, Neg_Py, Neg_Pz);
    return ROOT::Math::VectorUtil::Perp(momTot, momNeg);
}

/*
 * Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 * This function stores the point of closest approach, and returns its distance to the PV.
 */
Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z, Double_t refPointX,
                      Double_t refPointY, Double_t refPointZ) {
    ROOT::Math::XYZVector V0Momentum(V0_Px, V0_Py, V0_Pz);
    ROOT::Math::XYZPoint V0Vertex(V0_X, V0_Y, V0_Z);
    ROOT::Math::XYZPoint RefVertex(refPointX, refPointY, refPointZ);
    ROOT::Math::XYZVector CrossProduct = (RefVertex - V0Vertex).Cross(V0Momentum);
    return CrossProduct.R() / V0Momentum.R();
}

}  // namespace Math
}  // namespace Tree2Sexaquark
