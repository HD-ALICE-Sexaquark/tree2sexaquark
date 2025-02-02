/***********************************************************************************************************/
/* !!                                                                                                   !! */
/* !! The variables of this header file should be equivalent to the ones in `AliAnalysisTaskEsd2Tree.h` !! */
/* !!                                                                                                   !! */
/***********************************************************************************************************/

#ifndef T2S_INPUT_TREES_STRUCTURE_HXX
#define T2S_INPUT_TREES_STRUCTURE_HXX

#include "Rtypes.h"

/*
 * Simple Tree name: "Events"
 */
struct Event_tt {
    UInt_t RunNumber;         //
    UInt_t DirNumber;         //
    UInt_t DirNumberB;        //
    UInt_t EventNumber;       //
    Float_t Centrality;       //
    Float_t MagneticField;    //
    Float_t PV_TrueXv;        //
    Float_t PV_TrueYv;        //
    Float_t PV_TrueZv;        //
    Bool_t IsGenPileup;       //
    Bool_t IsSBCPileup;       //
    Int_t PV_NContributors;   //
    Float_t PV_Dispersion;    //
    Float_t PV_Xv;            //
    Float_t PV_Yv;            //
    Float_t PV_Zv;            //
    Float_t PV_CovMatrix[6];  //
    Float_t SPD_PV_Zv;        //
    Float_t SPD_PV_ZvErr;     //
    UInt_t NTracks;           //
    Int_t NTPCClusters;       //
    Bool_t IsMB;              //
    Bool_t IsHighMultV0;      //
    Bool_t IsHighMultSPD;     //
    Bool_t IsCentral;         //
    Bool_t IsSemiCentral;     //
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
    inline Bool_t HasMother() { return Idx_Mother != -1; }
    inline Bool_t IsSignal() { return Generator == 2; }
    inline Bool_t IsFirstGenSignal() { return IsSignal() && !HasMother(); }
    inline Bool_t IsFinalStateSignal() { return Generator == 2 && HasMother(); }
    inline Bool_t IsSecondary() { return IsSecFromMat || IsSecFromWeak || IsSignal(); }
};

/*
 * Simple Tree name: "Tracks"
 */
struct Track_tt {
    UInt_t Idx;                  //
    Float_t Px;                  // inner parametrization
    Float_t Py;                  // inner parametrization
    Float_t Pz;                  // inner parametrization
    Float_t X;                   //
    Float_t Y;                   //
    Float_t Z;                   //
    Short_t Charge;              //
    Float_t Alpha;               //
    Float_t Snp;                 // local sine of the track momentum azimuthal angle
    Float_t Tgl;                 // tangent of the track momentum dip angle
    Float_t Signed1Pt;           // 1/pt
    Float_t CovMatrix[15];       // covariance matrix
    Float_t NSigmaPion;          //
    Float_t NSigmaKaon;          //
    Float_t NSigmaProton;        //
    Float_t DCAxy;               // pre-calculated DCA wrt PV
    Float_t DCAz;                // pre-calculated DCA wrt PV
    UShort_t NTPCClusters;       //
    Float_t NCrossedRows;        //
    UShort_t NFindableClusters;  //
    UShort_t NSharedClusters;    //
    Float_t Chi2overNcls;        //
    Bool_t IsKinkDaughter;       //
    // TBits TPCFitMap;             //
    // TBits TPCClusterMap;         //
    // TBits TPCSharedMap;          //
    UInt_t Idx_True;  //
};

#endif  // T2S_INPUT_TREES_STRUCTURE_HXX
