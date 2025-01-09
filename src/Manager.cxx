#include "Trees/Manager.hxx"

Manager* Manager::Instance = nullptr;

/*
 * Get Trees Shepherd instance
 */
Manager* Manager::GetManager() {
    if (Instance == nullptr) Instance = new Manager;
    return Instance;
}

/*
 * Delete Trees Shepherd instance
 */
void Manager::DeleteManager() {
    if (Instance != nullptr) {
        delete Instance;
        Instance = nullptr;
    }
}

/*
 * Default private constructor
 */
Manager::Manager()
    : fTree_Events(nullptr),  //
      Event() {
    if (Instance) delete Instance;
    Instance = this;
}

/*
 * Private destructor
 */
Manager::~Manager() {
    if (fTree_Events) delete fTree_Events;
    fTree_Events = nullptr;
    Instance = nullptr;
}

/*
 *
 */
void Manager::Speak(Int_t a) {
    //
    InfoF("What's up? %i", a);
}

void Manager::ConnectEventBranches() {
    fTree_Events->SetBranchAddress("RunNumber", &Event.RunNumber);
    fTree_Events->SetBranchAddress("DirNumber", &Event.DirNumber);
    // if (!IsMC) fTree_Events->SetBranchAddress("DirNumberB", &Event.DirNumberB); // PENDING
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
