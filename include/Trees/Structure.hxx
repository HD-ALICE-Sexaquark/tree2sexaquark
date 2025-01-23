/***********************************************************************************************************/
/* !!                                                                                                   !! */
/* !! The variables of this header file should be equivalent to the ones in `AliAnalysisTaskEsd2Tree.h` !! */
/* !!                                                                                                   !! */
/***********************************************************************************************************/

#ifndef T2S_TREES_STRUCTURE_HXX
#define T2S_TREES_STRUCTURE_HXX

#include "Rtypes.h"

/*
 * Simple Tree name: "Events"
 */
struct Event_tt {
    UInt_t RunNumber;             //
    UInt_t DirNumber;             //
    UInt_t DirNumberB;            //
    UInt_t EventNumber;           //
    Float_t Centrality;           //
    Float_t PV_TrueXv;            //
    Float_t PV_TrueYv;            //
    Float_t PV_TrueZv;            //
    Bool_t IsGenPileup;           //
    Bool_t IsSBCPileup;           //
    Float_t PV_RecXv;             //
    Float_t PV_RecYv;             //
    Float_t PV_RecZv;             //
    Int_t PV_NContributors;       //
    Float_t PV_ZvErr_FromSPD;     //
    Float_t PV_ZvErr_FromTracks;  //
    Float_t PV_Zv_FromSPD;        //
    Float_t PV_Zv_FromTracks;     //
    Float_t PV_Dispersion;        //
    UInt_t NTracks;               //
    Int_t NTPCClusters;           //
    Bool_t IsMB;                  //
    Bool_t IsHighMultV0;          //
    Bool_t IsHighMultSPD;         //
    Bool_t IsCentral;             //
    Bool_t IsSemiCentral;         //
};

/*
 * Simple Tree name: "Injected"
 */
struct Injected_tt {
    UInt_t ReactionID;   //
    Float_t Px;          //
    Float_t Py;          //
    Float_t Pz;          //
    Float_t Nucleon_Px;  //
    Float_t Nucleon_Py;  //
    Float_t Nucleon_Pz;  //
};

/*
 * Simple Tree name: "MC"
 */
struct MC_tt {
    UInt_t Idx;            //
    Int_t PdgCode;         //
    Int_t Idx_Mother;      //
    Int_t Idx_Ancestor;    //
    Float_t Px;            //
    Float_t Py;            //
    Float_t Pz;            //
    Float_t Xv;            // origin x-vertex
    Float_t Yv;            // origin y-vertex
    Float_t Zv;            // origin z-vertex
    UInt_t Status;         //
    Bool_t IsOOBPileup;    //
    Short_t Generator;     // 0: HIJING, 1: anti-neutron injector, 2: anti-sexaquark reaction
    Bool_t IsPrimary;      //
    Bool_t IsSecFromMat;   //
    Bool_t IsSecFromWeak;  //
    Int_t ReactionID;      //
};

/*
 * Simple Tree name: "Tracks"
 */
struct Track_tt {
    UInt_t Idx;             //
    Float_t Px;             // inner parametrization
    Float_t Py;             // inner parametrization
    Float_t Pz;             // inner parametrization
    Float_t X;              //
    Float_t Y;              //
    Float_t Z;              //
    Int_t Charge;           //
    Float_t Alpha;          //
    Float_t Snp;            // local sine of the track momentum azimuthal angle
    Float_t Tgl;            // tangent of the track momentum dip angle
    Float_t Signed1Pt;      // 1/pt
    Float_t CovMatrix[15];  // covariance matrix
    // Float_t NSigmaPion;     //
    // Float_t NSigmaKaon;     //
    // Float_t NSigmaProton;   //
    // Float_t DCAxy;          // pre-calculated DCA wrt PV
    // Float_t DCAz;           // pre-calculated DCA wrt PV
    // UShort_t NTPCClusters;       //
    // Float_t NCrossedRows;        //
    // UShort_t NFindableClusters;  //
    // UShort_t NSharedClusters;    //
    // Float_t Chi2overNcls;   //
    // Bool_t IsKinkDaughter;  //
    // TBits TPCFitMap;             //
    // TBits TPCClusterMap;         //
    // TBits TPCSharedMap;          //
    Int_t Idx_True;  //
};

#endif  // T2S_TREES_STRUCTURE_HXX
