#include "TDirectoryFile.h"
#include "TFile.h"

#include "Utilities/Logger.hxx"

#include "Trees/Manager.hxx"

int main(int argc, char *argv[]) {

    Bool_t IsMC = kTRUE;
    TString Path_InputFile = "/home/ceres/borquez/some/esd2tree/task/attempts/local_signalMC_A1.8_18qr_test/SimpleTrees.root";

    /* Initialize Logger and TreeManager */

    Logger *ThisLogger = Logger::GetLogger();
    Manager *TreeManager = Manager::GetManager();

    /* Input File */

    TFile *file = new TFile(Path_InputFile, "READ");

    /* Events Simple Tree */

    TreeManager->SetEventsTree(file->Get<TTree>("Events"));
    TreeManager->ConnectEventBranches();

    /* The Event Loop */
    TString event_uid;
    // TDirectoryFile *event_dir;
    for (UInt_t evt_idx = 0; evt_idx < TreeManager->GetNEvents(); evt_idx++) {
        TreeManager->GetEvent(evt_idx);
        event_uid = TString::Format("A18_%i_%04u_%03u", TreeManager->Event.RunNumber, TreeManager->Event.EventNumber, TreeManager->Event.DirNumber);
        InfoF("Processing Event # %u (UID=%s)", evt_idx, event_uid.Data());
        // event_dir = file->Get<TDirectoryFile>(event_uid);
        // if (!event_dir) {
        // InfoF("TDirectoryFile %s couldn't be found, moving on...", event_uid.Data());
        // continue;
        // }
    }

    /*
        TList *list = file->GetListOfKeys();
        TIter next(list);
        TKey *key;

        TObject *obj;
        while ((key = (TKey *)next())) {
            obj = key->ReadObj();
            if (!obj->InheritsFrom(TDirectory::Class())) continue;
            event_dir = (TDirectoryFile *)obj;
            std::cout << "Event : " << obj->GetName() << std::endl;

            // Loop over MC Tree

            MC_tt mc;
            TTree *mc_tree = event_dir->Get<TTree>("MC");
            if (!mc_tree) continue;

            mc_tree->SetBranchAddress("Idx", &mc.Idx);
            mc_tree->SetBranchAddress("PdgCode", &mc.PdgCode);
            mc_tree->SetBranchAddress("Idx_Mother", &mc.Idx_Mother);
            mc_tree->SetBranchAddress("Idx_Ancestor", &mc.Idx_Ancestor);
            mc_tree->SetBranchAddress("Px", &mc.Px);
            mc_tree->SetBranchAddress("Py", &mc.Py);
            mc_tree->SetBranchAddress("Pz", &mc.Pz);
            mc_tree->SetBranchAddress("Xv", &mc.Xv);
            mc_tree->SetBranchAddress("Yv", &mc.Yv);
            mc_tree->SetBranchAddress("Zv", &mc.Zv);
            mc_tree->SetBranchAddress("Status", &mc.Status);
            mc_tree->SetBranchAddress("IsOOBPileup", &mc.IsOOBPileup);
            mc_tree->SetBranchAddress("Generator", &mc.Generator);
            mc_tree->SetBranchAddress("IsPrimary", &mc.IsPrimary);
            mc_tree->SetBranchAddress("IsSecFromMat", &mc.IsSecFromMat);
            mc_tree->SetBranchAddress("IsSecFromWeak", &mc.IsSecFromWeak);
            mc_tree->SetBranchAddress("ReactionID", &mc.ReactionID);

            for (Int_t mc_idx = 0; mc_idx < mc_tree->GetEntries(); mc_idx++) {
                // if (mc_idx > 0) break;
                mc_tree->GetEntry(mc_idx);
                std::cout << "MC Entry : " << mc_idx << std::endl;
                std::cout << TString::Format("%i, %i, %f, %f, %f", mc.PdgCode, mc.Idx_Mother, mc.Px, mc.Py, mc.Pz) << std::endl;
            }
        }
        */
}
