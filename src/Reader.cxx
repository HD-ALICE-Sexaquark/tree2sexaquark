#include "Trees/Reader.hxx"

/*
 *
 */
void Reader::ConnectEventBranches(Bool_t IsMC) {
    fTree_Events->SetBranchAddress("RunNumber", &Event.RunNumber);
    fTree_Events->SetBranchAddress("DirNumber", &Event.DirNumber);
    if (!IsMC) fTree_Events->SetBranchAddress("DirNumberB", &Event.DirNumberB);  // PENDING
    fTree_Events->SetBranchAddress("EventNumber", &Event.EventNumber);
    fTree_Events->SetBranchAddress("Centrality", &Event.Centrality);
    fTree_Events->SetBranchAddress("PV_TrueXv", &Event.PV_TrueXv);
    fTree_Events->SetBranchAddress("PV_TrueYv", &Event.PV_TrueYv);
    fTree_Events->SetBranchAddress("PV_TrueZv", &Event.PV_TrueZv);
    fTree_Events->SetBranchAddress("IsGenPileup", &Event.IsGenPileup);
    fTree_Events->SetBranchAddress("IsSBCPileup", &Event.IsSBCPileup);
    fTree_Events->SetBranchAddress("PV_RecXv", &Event.PV_RecXv);
    fTree_Events->SetBranchAddress("PV_RecYv", &Event.PV_RecYv);
    fTree_Events->SetBranchAddress("PV_RecZv", &Event.PV_RecZv);
    fTree_Events->SetBranchAddress("PV_NContributors", &Event.PV_NContributors);
    fTree_Events->SetBranchAddress("PV_ZvErr_FromSPD", &Event.PV_ZvErr_FromSPD);
    fTree_Events->SetBranchAddress("PV_ZvErr_FromTracks", &Event.PV_ZvErr_FromTracks);
    fTree_Events->SetBranchAddress("PV_Zv_FromSPD", &Event.PV_Zv_FromSPD);
    fTree_Events->SetBranchAddress("PV_Zv_FromTracks", &Event.PV_Zv_FromTracks);
    fTree_Events->SetBranchAddress("PV_Dispersion", &Event.PV_Dispersion);
    fTree_Events->SetBranchAddress("NTracks", &Event.NTracks);
    fTree_Events->SetBranchAddress("NTPCClusters", &Event.NTPCClusters);
    fTree_Events->SetBranchAddress("IsMB", &Event.IsMB);
    fTree_Events->SetBranchAddress("IsHighMultV0", &Event.IsHighMultV0);
    fTree_Events->SetBranchAddress("IsHighMultSPD", &Event.IsHighMultSPD);
    fTree_Events->SetBranchAddress("IsCentral", &Event.IsCentral);
    fTree_Events->SetBranchAddress("IsSemiCentral", &Event.IsSemiCentral);
}

/*
 *
 */
void Reader::ConnectInjectedBranches() {
    fTree_Injected->SetBranchAddress("ReactionID", &Injected.ReactionID);
    fTree_Injected->SetBranchAddress("Px", &Injected.Px);
    fTree_Injected->SetBranchAddress("Py", &Injected.Py);
    fTree_Injected->SetBranchAddress("Pz", &Injected.Pz);
    fTree_Injected->SetBranchAddress("Nucleon_Px", &Injected.Nucleon_Px);
    fTree_Injected->SetBranchAddress("Nucleon_Py", &Injected.Nucleon_Py);
    fTree_Injected->SetBranchAddress("Nucleon_Pz", &Injected.Nucleon_Pz);
}

/*
 *
 */
void Reader::ConnectMCBranches() {
    fTree_MC->SetBranchAddress("Idx", &MC.Idx);
    fTree_MC->SetBranchAddress("PdgCode", &MC.PdgCode);
    fTree_MC->SetBranchAddress("Idx_Mother", &MC.Idx_Mother);
    fTree_MC->SetBranchAddress("Idx_Ancestor", &MC.Idx_Ancestor);
    fTree_MC->SetBranchAddress("Px", &MC.Px);
    fTree_MC->SetBranchAddress("Py", &MC.Py);
    fTree_MC->SetBranchAddress("Pz", &MC.Pz);
    fTree_MC->SetBranchAddress("Xv", &MC.Xv);
    fTree_MC->SetBranchAddress("Yv", &MC.Yv);
    fTree_MC->SetBranchAddress("Zv", &MC.Zv);
    fTree_MC->SetBranchAddress("Status", &MC.Status);
    fTree_MC->SetBranchAddress("IsOOBPileup", &MC.IsOOBPileup);
    fTree_MC->SetBranchAddress("Generator", &MC.Generator);
    fTree_MC->SetBranchAddress("IsPrimary", &MC.IsPrimary);
    fTree_MC->SetBranchAddress("IsSecFromMat", &MC.IsSecFromMat);
    fTree_MC->SetBranchAddress("IsSecFromWeak", &MC.IsSecFromWeak);
    fTree_MC->SetBranchAddress("ReactionID", &MC.ReactionID);
}

/*
 *
 */
void Reader::ConnectTracksBranches() {
    fTree_Tracks->SetBranchAddress("Idx", &Track.Idx);
    fTree_Tracks->SetBranchAddress("Px", &Track.Px);
    fTree_Tracks->SetBranchAddress("Py", &Track.Py);
    fTree_Tracks->SetBranchAddress("Pz", &Track.Pz);
    fTree_Tracks->SetBranchAddress("X", &Track.X);
    fTree_Tracks->SetBranchAddress("Y", &Track.Y);
    fTree_Tracks->SetBranchAddress("Z", &Track.Z);
    // fTree_Tracks->SetBranchAddress("Charge", &Track.Charge);
    // fTree_Tracks->SetBranchAddress("Alpha", &Track.Alpha);
    // fTree_Tracks->SetBranchAddress("Snp", &Track.Snp);
    // fTree_Tracks->SetBranchAddress("Tgl", &Track.Tgl);
    // fTree_Tracks->SetBranchAddress("Signed1Pt", &Track.Signed1Pt);
    // fTree_Tracks->SetBranchAddress("CovMatrix", &Track.CovMatrix[15]);
    // fTree_Tracks->SetBranchAddress("NSigmaPion", &Track.NSigmaPion);
    // fTree_Tracks->SetBranchAddress("NSigmaKaon", &Track.NSigmaKaon);
    fTree_Tracks->SetBranchAddress("NSigmaProton", &Track.NSigmaProton);
    // fTree_Tracks->SetBranchAddress("DCAxy", &Track.DCAxy);
    // fTree_Tracks->SetBranchAddress("DCAz", &Track.DCAz);
    // fTree_Tracks->SetBranchAddress("NTPCClusters", &Track.NTPCClusters);
    // fTree_Tracks->SetBranchAddress("NCrossedRows", &Track.NCrossedRows);
    // fTree_Tracks->SetBranchAddress("NFindableClusters", &Track.NFindableClusters);
    // fTree_Tracks->SetBranchAddress("NSharedClusters", &Track.NSharedClusters);
    // fTree_Tracks->SetBranchAddress("Chi2overNcls", &Track.Chi2overNcls);
    // fTree_Tracks->SetBranchAddress("IsKinkDaughter", &Track.IsKinkDaughter);
    fTree_Tracks->SetBranchAddress("Idx_True", &Track.Idx_True);
}
