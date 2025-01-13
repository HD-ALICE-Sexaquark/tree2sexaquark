#include "Utilities/Logger.hxx"

#include "Analysis/Manager.hxx"

int main(int argc, char *argv[]) {

    Logger *ThisLogger = Logger::GetInstance();

    AnalysisManager *ThisAnalysis = AnalysisManager::Init();
    if (!ThisAnalysis->Configure(argc, argv)) return 1;
    if (!ThisAnalysis->OpenInputFile()) return 1;

    for (UInt_t evt_idx = 0; evt_idx < ThisAnalysis->GetN_Events(); evt_idx++) {

        if (!ThisAnalysis->GetEvent(evt_idx)) continue;

        if (ThisAnalysis->IsMC()) {
            ThisAnalysis->ProcessInjected();
            ThisAnalysis->ProcessMCParticles();
        }

        ThisAnalysis->ProcessTracks();
    }

    ThisAnalysis->DeleteInstance();
    ThisLogger->DeleteInstance();
}
