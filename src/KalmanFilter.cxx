#include "Math/KalmanFilter.hxx"

#include "TMath.h"

namespace Tree2Sexaquark {
namespace Math {

/*
 * Correct initialization of a KFParticle.
 * (Copied from (https://github.com/alisw/AliPhysics)/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFParticle CreateKFParticle(Float_t px, Float_t py, Float_t pz, Float_t x, Float_t y, Float_t z, Int_t charge,  //
                            Float_t alpha, Float_t snp, Float_t tgl, Float_t signed1pt, Float_t* cov_matrix, Double_t mass) {

    Double_t fP[6];
    fP[0] = x;
    fP[1] = y;
    fP[2] = z;

    Int_t fQ = charge;
    fP[3] = TMath::Abs(charge) * px;
    fP[4] = TMath::Abs(charge) * py;
    fP[5] = TMath::Abs(charge) * pz;

    Double_t pt = 1. / TMath::Abs(signed1pt) * TMath::Abs(charge);
    Double_t cs = TMath::Cos(alpha);
    Double_t sn = TMath::Sin(alpha);
    Double_t r = TMath::Sqrt((1. - snp) * (1. + snp));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + snp * cs / r);
    Double_t m43 = -pt * pt * (r * cs - snp * sn);
    Double_t m24 = pt * (cs - snp * sn / r);
    Double_t m44 = -pt * pt * (r * sn + snp * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * tgl;

    m43 *= ((Double_t)charge / (Double_t)TMath::Abs(charge));
    m44 *= ((Double_t)charge / (Double_t)TMath::Abs(charge));
    m45 *= ((Double_t)charge / (Double_t)TMath::Abs(charge));

    Double_t fC[21];
    fC[0] = cov_matrix[0] * m00 * m00;
    fC[1] = cov_matrix[0] * m00 * m10;
    fC[2] = cov_matrix[0] * m10 * m10;
    fC[3] = cov_matrix[1] * m00;
    fC[4] = cov_matrix[1] * m10;
    fC[5] = cov_matrix[2];
    fC[6] = m00 * (cov_matrix[3] * m23 + cov_matrix[10] * m43);
    fC[7] = m10 * (cov_matrix[3] * m23 + cov_matrix[10] * m43);
    fC[8] = cov_matrix[4] * m23 + cov_matrix[11] * m43;
    fC[9] = m23 * (cov_matrix[5] * m23 + cov_matrix[12] * m43) + m43 * (cov_matrix[12] * m23 + cov_matrix[14] * m43);
    fC[10] = m00 * (cov_matrix[3] * m24 + cov_matrix[10] * m44);
    fC[11] = m10 * (cov_matrix[3] * m24 + cov_matrix[10] * m44);
    fC[12] = cov_matrix[4] * m24 + cov_matrix[11] * m44;
    fC[13] = m23 * (cov_matrix[5] * m24 + cov_matrix[12] * m44) + m43 * (cov_matrix[12] * m24 + cov_matrix[14] * m44);
    fC[14] = m24 * (cov_matrix[5] * m24 + cov_matrix[12] * m44) + m44 * (cov_matrix[12] * m24 + cov_matrix[14] * m44);
    fC[15] = m00 * (cov_matrix[6] * m35 + cov_matrix[10] * m45);
    fC[16] = m10 * (cov_matrix[6] * m35 + cov_matrix[10] * m45);
    fC[17] = cov_matrix[7] * m35 + cov_matrix[11] * m45;
    fC[18] = m23 * (cov_matrix[8] * m35 + cov_matrix[12] * m45) + m43 * (cov_matrix[13] * m35 + cov_matrix[14] * m45);
    fC[19] = m24 * (cov_matrix[8] * m35 + cov_matrix[12] * m45) + m44 * (cov_matrix[13] * m35 + cov_matrix[14] * m45);
    fC[20] = m35 * (cov_matrix[9] * m35 + cov_matrix[13] * m45) + m45 * (cov_matrix[13] * m35 + cov_matrix[14] * m45);

    KFParticle part;
    part.Create(fP, fC, fQ, mass);

    return part;
}

/*
 * Correct initialization of a KFVertex.
 * - Arguments: XYZ[3], CovarianceMatrix[6]
 * (Copied from (https://github.com/alisw/AliPhysics)/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFVertex CreateKFVertex(Float_t* xyz, Float_t* cov_matrix) {
    //
    KFPVertex kfpVtx;
    kfpVtx.SetXYZ(xyz);
    kfpVtx.SetCovarianceMatrix(cov_matrix);

    KFVertex KFVtx(kfpVtx);
    return KFVtx;
}

/*
 * Transport a KFParticle to the point of closest approach w.r.t. another KFParticle.
 */
KFParticle TransportKFParticle(KFParticle kf_this, KFParticle kf_other, Double_t mass_this, Int_t charge_this) {
    //
    Float_t dS[2];
    Float_t dsdr[4][6];
    kf_this.GetDStoParticle(kf_other, dS, dsdr);

    Float_t mP[8], mC[36];
    kf_this.Transport(dS[0], dsdr[0], mP, mC);

    Float_t mM = mass_this;
    Float_t mQ = charge_this;  // only valid for charged particles with Q = +/- 1
    KFParticle kfTransported;
    kfTransported.Create(mP, mC, mQ, mM);

    return kfTransported;
}

}  // namespace Math
}  // namespace Tree2Sexaquark
