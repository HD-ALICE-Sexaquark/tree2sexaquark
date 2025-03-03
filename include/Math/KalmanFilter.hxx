#ifndef T2S_KALMAN_FILTER_HXX
#define T2S_KALMAN_FILTER_HXX

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"
#include "KFVertex.h"

#include "ROOT/RVec.hxx"
using namespace ROOT::VecOps;
using RVecF = ROOT::RVecF;
using cRVecF = const ROOT::RVecF&;

namespace Tree2Sexaquark::Math {

KFParticle CreateKFParticle(Float_t px, Float_t py, Float_t pz, Float_t x, Float_t y, Float_t z,       //
                            Int_t charge, Float_t alpha, Float_t snp, Float_t tgl, Float_t signed1pt,  //
                            cRVecF cov_matrix, Double_t mass);
KFVertex CreateKFVertex(Float_t* xyz, Float_t* cov_matrix);

}  // namespace Tree2Sexaquark::Math

#endif  // T2S_KALMAN_FILTER_HXX
