#include "TMath.h"

#include "Math/KalmanFilter.hxx"

namespace Tree2Sexaquark {
namespace Math {
/*
 * Correct initialization of a KFParticle.
 * (Copied from (https://github.com/alisw/AliPhysics)/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFParticle CreateKFParticle(Track_tt track, Double_t mass) {

    Double_t fP[6];
    fP[0] = track.X;
    fP[1] = track.Y;
    fP[2] = track.Z;

    Int_t fQ = track.Charge;
    fP[3] = TMath::Abs(track.Charge) * track.Px;
    fP[4] = TMath::Abs(track.Charge) * track.Py;
    fP[5] = TMath::Abs(track.Charge) * track.Pz;

    Double_t pt = 1. / TMath::Abs(track.Signed1Pt) * TMath::Abs(track.Charge);
    Double_t cs = TMath::Cos(track.Alpha);
    Double_t sn = TMath::Sin(track.Alpha);
    Double_t r = TMath::Sqrt((1. - track.Snp) * (1. + track.Snp));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + track.Snp * cs / r);
    Double_t m43 = -pt * pt * (r * cs - track.Snp * sn);
    Double_t m24 = pt * (cs - track.Snp * sn / r);
    Double_t m44 = -pt * pt * (r * sn + track.Snp * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * track.Tgl;

    m43 *= ((Double_t)track.Charge / (Double_t)TMath::Abs(track.Charge));
    m44 *= ((Double_t)track.Charge / (Double_t)TMath::Abs(track.Charge));
    m45 *= ((Double_t)track.Charge / (Double_t)TMath::Abs(track.Charge));

    Double_t fC[21];
    fC[0] = track.CovMatrix[0] * m00 * m00;
    fC[1] = track.CovMatrix[0] * m00 * m10;
    fC[2] = track.CovMatrix[0] * m10 * m10;
    fC[3] = track.CovMatrix[1] * m00;
    fC[4] = track.CovMatrix[1] * m10;
    fC[5] = track.CovMatrix[2];
    fC[6] = m00 * (track.CovMatrix[3] * m23 + track.CovMatrix[10] * m43);
    fC[7] = m10 * (track.CovMatrix[3] * m23 + track.CovMatrix[10] * m43);
    fC[8] = track.CovMatrix[4] * m23 + track.CovMatrix[11] * m43;
    fC[9] = m23 * (track.CovMatrix[5] * m23 + track.CovMatrix[12] * m43) + m43 * (track.CovMatrix[12] * m23 + track.CovMatrix[14] * m43);
    fC[10] = m00 * (track.CovMatrix[3] * m24 + track.CovMatrix[10] * m44);
    fC[11] = m10 * (track.CovMatrix[3] * m24 + track.CovMatrix[10] * m44);
    fC[12] = track.CovMatrix[4] * m24 + track.CovMatrix[11] * m44;
    fC[13] = m23 * (track.CovMatrix[5] * m24 + track.CovMatrix[12] * m44) + m43 * (track.CovMatrix[12] * m24 + track.CovMatrix[14] * m44);
    fC[14] = m24 * (track.CovMatrix[5] * m24 + track.CovMatrix[12] * m44) + m44 * (track.CovMatrix[12] * m24 + track.CovMatrix[14] * m44);
    fC[15] = m00 * (track.CovMatrix[6] * m35 + track.CovMatrix[10] * m45);
    fC[16] = m10 * (track.CovMatrix[6] * m35 + track.CovMatrix[10] * m45);
    fC[17] = track.CovMatrix[7] * m35 + track.CovMatrix[11] * m45;
    fC[18] = m23 * (track.CovMatrix[8] * m35 + track.CovMatrix[12] * m45) + m43 * (track.CovMatrix[13] * m35 + track.CovMatrix[14] * m45);
    fC[19] = m24 * (track.CovMatrix[8] * m35 + track.CovMatrix[12] * m45) + m44 * (track.CovMatrix[13] * m35 + track.CovMatrix[14] * m45);
    fC[20] = m35 * (track.CovMatrix[9] * m35 + track.CovMatrix[13] * m45) + m45 * (track.CovMatrix[13] * m35 + track.CovMatrix[14] * m45);

    KFParticle part;
    part.Create(fP, fC, fQ, mass);

    return part;
}

/*
 * Correct initialization of a KFVertex.
 * - Arguments: XYZ[3], CovarianceMatrix[6]
 * (Copied from (https://github.com/alisw/AliPhysics)/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFVertex CreateKFVertex(Float_t* XYZ, Float_t* CovarianceMatrix) {

    KFPVertex kfpVtx;
    kfpVtx.SetXYZ(XYZ);
    kfpVtx.SetCovarianceMatrix(CovarianceMatrix);

    KFVertex KFVtx(kfpVtx);

    return KFVtx;
}

/*
 * Transport a KFParticle to the point of closest approach w.r.t. another KFParticle.
 */
KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Double_t massThis, Int_t chargeThis) {

    Float_t dS[2];
    Float_t dsdr[4][6];
    kfThis.GetDStoParticle(kfOther, dS, dsdr);
    // GetDStoParticleBz(fMagneticField, kfThis, kfOther, dS, dsdr); // TEST

    Float_t mP[8], mC[36];
    kfThis.Transport(dS[0], dsdr[0], mP, mC);

    Float_t mM = massThis;
    Float_t mQ = chargeThis;  // only valid for charged particles with Q = +/- 1

    KFParticle kfTransported;
    kfTransported.Create(mP, mC, mQ, mM);

    return kfTransported;
}

}  // namespace Math
}  // namespace Tree2Sexaquark
