#ifndef T2S_DEFAULT_CUTS_HXX
#define T2S_DEFAULT_CUTS_HXX

#include <map>
#include <tuple>

#include "RtypesCore.h"
#include "TString.h"

namespace Tree2Sexaquark {
namespace Default {

namespace Proton {
const Double_t MaxNSigma = 3.;
}  // namespace Proton

namespace Kaon {
const Double_t MaxNSigma = 3.;
}  // namespace Kaon

namespace Pion {
const Double_t MaxNSigma = 3.;
}  // namespace Pion

namespace Lambda {
const Double_t AbsMaxEta = 0.9;
const Double_t MinPt = 1.0;
const Double_t MinMass = 1.08;
const Double_t MaxMass = 1.16;
const Double_t MinRadius = 20.;
const Double_t MinDistFromPV = 40.;
const Double_t MinCPAwrtPV = 0.1;
const Double_t MaxCPAwrtPV = 0.99;
const Double_t MinDCAwrtPV = 4.;
const Double_t MaxDCAbtwDau = 2.;
const Double_t MaxDCAnegV0 = 2.;
const Double_t MaxDCAposV0 = 2.;
const Double_t MaxArmQtOverAlpha = 0.2;
}  // namespace Lambda

namespace KaonZero {
const Double_t AbsMaxEta = 0.8;
const Double_t MinPt = 1.0;
const Double_t MinMass = 0.475;
const Double_t MaxMass = 0.525;
const Double_t MinRadius = 15.;
const Double_t MinDistFromPV = 30.;
const Double_t MaxDistFromPV = 175.;
const Double_t MinCPAwrtPV = 0.2;
const Double_t MaxCPAwrtPV = 0.99;
const Double_t MinDCAwrtPV = 4.;
const Double_t MaxDCAbtwDau = 0.25;
const Double_t MaxDCAnegV0 = 0.25;
const Double_t MaxDCAposV0 = 0.25;
}  // namespace KaonZero

namespace PionPair {
const Double_t AbsMaxEta = 0.9;
const Double_t MinPt = 0.1;
const Double_t MinRadius = 15.;
const Double_t MinCPAwrtPV = 0.8;
const Double_t MaxCPAwrtPV = 0.98;
const Double_t MinDCAwrtPV = 4.;
const Double_t MaxDCAbtwDau = 0.3;
const Double_t MaxDCAnegV0 = 0.2;
const Double_t MaxDCAposV0 = 0.2;
}  // namespace PionPair

// Min_KaonPair_Pt = 0.2; // PENDING
// Max_SexaquarkD_DCAv0SV = 0.2; // PENDING
// Max_SexaquarkE_DCAv0SV = 0.1; // PENDING

namespace ChannelA {
const Double_t MinRadius = 30.;
const Double_t MaxRapidity = 0.8;
const Double_t MinCPAwrtPV = 0.99;
const Double_t MaxCPAwrtPV = 1.;
const Double_t MaxDCALaSV = 2.;
const Double_t MaxDCALaNegSV = 10.;
const Double_t MaxDCALaPosSV = 10.;
const Double_t MaxDCAK0SV = 2.;
const Double_t MaxDCAK0NegSV = 10.;
const Double_t MaxDCAK0PosSV = 10.;
const Double_t MaxDCAbtwV0s = 2.;
}  // namespace ChannelA

namespace ChannelD {
const Double_t MinRadius = 30.;
const Double_t MaxRapidity = 0.8;
const Double_t MinCPAwrtPV = 0.99;
const Double_t MaxCPAwrtPV = 1.;
const Double_t MaxDCALaSV = 2.;
const Double_t MaxDCALaNegSV = 10.;
const Double_t MaxDCALaPosSV = 10.;
const Double_t MaxDCAKaSV = 0.2;
const Double_t MaxDCAKaLa = 0.3;
}  // namespace ChannelD

namespace ChannelE {
const Double_t MinRadius = 30.;
const Double_t MaxRapidity = 0.8;
const Double_t MinCPAwrtPV = 0.99;
const Double_t MaxCPAwrtPV = 1.;
const Double_t MaxDCALaSV = 2.;
const Double_t MaxDCALaNegSV = 10.;
const Double_t MaxDCALaPosSV = 10.;
const Double_t MaxDCAKaSV = 0.1;
const Double_t MaxDCAKaLa = 0.1;
const Double_t MaxDCApmSV = 0.1;
const Double_t MaxDCAppSV = 0.1;
const Double_t MaxDCApmLa = 0.1;
const Double_t MaxDCAppLa = 0.1;
const Double_t MaxDCApmKa = 0.1;
const Double_t MaxDCAppKa = 0.1;
}  // namespace ChannelE

namespace KaonPair {
const Double_t MinRadius = 30.;
const Double_t MaxRapidity = 0.8;
const Double_t MinCPAwrtPV = 0.99;
const Double_t MaxCPAwrtPV = 1.;
const Double_t MaxDCAbtwKK = 0.2;
const Double_t MaxDCAkaSV = 0.2;
const Double_t MaxDCAkbSV = 0.2;
}  // namespace KaonPair

}  // namespace Default
}  // namespace Tree2Sexaquark

#endif  // T2S_DEFAULT_CUTS_HXX
