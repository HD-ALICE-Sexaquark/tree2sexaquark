#ifndef T2S_DEFAULT_CUTS_HXX
#define T2S_DEFAULT_CUTS_HXX

#include <map>
#include <tuple>

#include "RtypesCore.h"
#include "TString.h"

namespace Tree2Sexaquark {
namespace Default {

const Float_t Lambda_AbsMaxEta = 0.9;
const Float_t Lambda_MinPt = 1.0;
const Float_t Lambda_MinMass = 1.08;
const Float_t Lambda_MaxMass = 1.16;
const Float_t Lambda_MinRadius = 20.;
const Float_t Lambda_MinDistFromPV = 40.;
const Float_t Lambda_MinCPAwrtPV = 0.1;
const Float_t Lambda_MaxCPAwrtPV = 0.99;
const Float_t Lambda_MinDCAwrtPV = 4.;
const Float_t Lambda_MaxDCAbtwDau = 2.;
const Float_t Lambda_MaxDCAnegV0 = 2.;
const Float_t Lambda_MaxDCAposV0 = 2.;
const Float_t Lambda_MaxArmPtOverAlpha = 0.2;

const Float_t KaonZero_AbsMaxEta = 0.8;
const Float_t KaonZero_MinPt = 1.0;
const Float_t KaonZero_MinMass = 0.475;
const Float_t KaonZero_MaxMass = 0.525;
const Float_t KaonZero_MinRadius = 15.;
const Float_t KaonZero_MinDistFromPV = 30.;
const Float_t KaonZero_MaxDistFromPV = 175.;
const Float_t KaonZero_MinCPAwrtPV = 0.2;
const Float_t KaonZero_MaxCPAwrtPV = 0.99;
const Float_t KaonZero_MinDCAwrtPV = 4.;
const Float_t KaonZero_MaxDCAbtwDau = 0.25;
const Float_t KaonZero_MaxDCAnegV0 = 0.25;
const Float_t KaonZero_MaxDCAposV0 = 0.25;

const Float_t PionPair_AbsMaxEta = 0.9;
const Float_t PionPair_MinPt = 0.1;
const Float_t PionPair_MinRadius = 15.;
const Float_t PionPair_MinCPAwrtPV = 0.8;
const Float_t PionPair_MaxCPAwrtPV = 0.98;
const Float_t PionPair_MinDCAwrtPV = 4.;
const Float_t PionPair_MaxDCAbtwDau = 0.3;
const Float_t PionPair_MaxDCAnegV0 = 0.2;
const Float_t PionPair_MaxDCAposV0 = 0.2;

}  // namespace Default
}  // namespace Tree2Sexaquark

#endif  // T2S_DEFAULT_CUTS_HXX
