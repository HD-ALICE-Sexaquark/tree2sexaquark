#ifndef T2S_OUTPUT_TREES_STRUCTURE_HXX
#define T2S_OUTPUT_TREES_STRUCTURE_HXX

#include "Rtypes.h"

/*
 *
 */
struct V0_tt {
    UInt_t Idx;
    UInt_t Idx_Neg;
    UInt_t Idx_Pos;
    Int_t PID;
    Float_t Px;
    Float_t Py;
    Float_t Pz;
    Float_t E;
    Float_t Xv;  //! V0 x-vertex
    Float_t Yv;  //! V0 y-vertex
    Float_t Zv;  //! V0 z-vertex
    Float_t Neg_Px;
    Float_t Neg_Py;
    Float_t Neg_Pz;
    Float_t Pos_Px;
    Float_t Pos_Py;
    Float_t Pos_Pz;
    /* True Information */
    Int_t Idx_True;
    Int_t True_PdgCode;
    Bool_t IsSecondary;
    Bool_t IsSignal;
    Int_t ReactionID;
    Bool_t IsHybrid;
};

/*
 * Common Properties of (Anti-)Sexaquark Candidates
 */
struct Sexaquark_tt {
    Float_t Px;
    Float_t Py;
    Float_t Pz;
    Float_t E;
    Float_t E_asDecay;
    Float_t Xv;  //! secondary x-vertex
    Float_t Yv;  //! secondary y-vertex
    Float_t Zv;  //! secondary z-vertex
    Float_t DistFromPV;
    Float_t CPAwrtPV;
    Float_t DCAwrtPV;
    Float_t Chi2ndf;
    /* True Information */
    Bool_t IsSignal;
    Int_t ReactionID;
    Bool_t IsHybrid;
    Bool_t IsNonCombBkg;
    Int_t AncestorIdx;
};

/*
 *
 */
struct SexaquarkA_tt : public Sexaquark_tt {
    /* Shared with Channels "A"+"D"+"E" */
    Int_t Idx_Lambda;
    Int_t Idx_Lambda_Neg;
    Int_t Idx_Lambda_Pos;
    Float_t Lambda_DecayLength;
    Float_t DCALaSV;
    Float_t DCALaNegSV;
    Float_t DCALaPosSV;
    /* Shared with Channels "A"+"D"+"H" */
    Float_t OpeningAngle;
    /* Specific to Channel "A" */
    Int_t Idx_K0S;
    Int_t Idx_K0S_Neg;
    Int_t Idx_K0S_Pos;
    Float_t K0S_DecayLength;
    Float_t DCAK0SV;
    Float_t DCAK0NegSV;
    Float_t DCAK0PosSV;
    Float_t DCAbtwV0s;
};

/*
 *
 */
struct SexaquarkD_tt : public Sexaquark_tt {
    /* Shared with Channels "A"+"D"+"E" */
    Int_t Idx_Lambda;
    Int_t Idx_Lambda_Neg;
    Int_t Idx_Lambda_Pos;
    Float_t Lambda_DecayLength;
    Float_t DCALaSV;
    Float_t DCALaNegSV;
    Float_t DCALaPosSV;
    /* Shared with Channels "A"+"D"+"H" */
    Float_t OpeningAngle;
    /* Specific to Channel "D" */
    Int_t Idx_Kaon;
    Float_t DCAKaSV;
    Float_t DCAKaLa;
    Float_t DCALaNegKa;
    Float_t DCALaPosKa;
};

/*
 *
 */
struct SexaquarkE_tt : public Sexaquark_tt {
    /* Shared with Channels "A"+"D"+"E" */
    Int_t Idx_Lambda;
    Int_t Idx_Lambda_Neg;
    Int_t Idx_Lambda_Pos;
    Float_t Lambda_DecayLength;
    Float_t DCALaSV;
    Float_t DCALaNegSV;
    Float_t DCALaPosSV;
    /* Specific to Channel "E" */
    Int_t Idx_Kaon;
    Int_t Idx_PP;
    Int_t Idx_PiMinus;
    Int_t Idx_PiPlus;
    Float_t DCAKaSV;
    Float_t DCAKaLa;
    Float_t DCApmSV;
    Float_t DCAppSV;
    Float_t DCApmLa;
    Float_t DCAppLa;
    Float_t DCApmKa;
    Float_t DCAppKa;
};

/*
 *
 */
struct KaonPair_tt : public Sexaquark_tt {
    /* Shared with Channels "A"+"D"+"H" */
    Float_t OpeningAngle;
    /* Specific to Channel "H" */
    UInt_t Idx_KaonA;
    UInt_t Idx_KaonB;
    Float_t DCAbtwKK;
    Float_t DCAkaSV;
    Float_t DCAkbSV;
};

#endif  // T2S_OUTPUT_TREES_STRUCTURE_HXX
