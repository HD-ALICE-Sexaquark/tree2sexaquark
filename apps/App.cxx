#include "Analysis/Manager.hxx"
#include "Utilities/Parser.hxx"

int main(int argc, char *argv[]) {

    Logger *LoggerInstance = Logger::GetInstance();

    std::unique_ptr<Parser> ParserInstance =
        std::make_unique<Parser>("Tree2Sexaquark -- Read SimpleTrees.root and create Sexaquark Analysis Results!");
    ParserInstance->Parse(argc, argv);
    if (ParserInstance->HelpOrError) return ParserInstance->ExitCode;

    std::unique_ptr<AnalysisManager> ThisAnalysis = std::make_unique<AnalysisManager>(ParserInstance->GetSettings());
    ThisAnalysis->Print();
    if (!ThisAnalysis->OpenInputFile()) return 1;
    if (!ThisAnalysis->PrepareOutputFile()) return 1;

    for (Long64_t evt_entry = 0; evt_entry < ThisAnalysis->GetN_Events(); evt_entry++) {
        if (!ThisAnalysis->GetEvent(evt_entry)) continue;
        if (ThisAnalysis->IsMC()) {
            if (ThisAnalysis->IsSignalMC()) ThisAnalysis->ProcessInjected();
            ThisAnalysis->ProcessMCParticles();
        }
        /* Tracks */
        ThisAnalysis->ProcessTracks();
        /* Findables */
        if (ThisAnalysis->IsMC()) {
            ThisAnalysis->ProcessFindableV0s();
            if (ThisAnalysis->IsSignalMC()) ThisAnalysis->ProcessFindableSexaquarks();
        }
        /* V0s */
        ThisAnalysis->KalmanV0Finder(-2212, 211, -3122);
        ThisAnalysis->KalmanV0Finder(-211, 2212, 3122);
        ThisAnalysis->KalmanV0Finder(-211, 211);
        /* Sexaquarks */
        /*
        ThisAnalysis->KalmanSexaquarkFinder(2112, {-3122, 310});              // `AntiSexaquark,Neutron -> AntiLambda,K0S`
        ThisAnalysis->KalmanSexaquarkFinder(-2112, {3122, 310});              // `Sexaquark,AntiNeutron -> Lambda,K0S`
        ThisAnalysis->KalmanSexaquarkFinder(2212, {-3122, 321, -211, 211});   // `AntiSexaquark,Proton -> AntiLambda,K+,(pi-,pi+)`
        ThisAnalysis->KalmanSexaquarkFinder(-2212, {3122, -321, -211, 211});  // `Sexaquark,AntiProton -> Lambda,K-,(pi-,pi+)`
        ThisAnalysis->KalmanSexaquarkFinder(2212, {321, 321});                // `AntiSexaquark,Proton -> K+,K+,X`
        ThisAnalysis->KalmanSexaquarkFinder(-2212, {-321, -321});             // `Sexaquark,AntiProton -> K-,K-,X`
        */
        /* End of Event */
        ThisAnalysis->EndOfEvent();
    }
    ThisAnalysis->EndOfAnalysis();

    Logger::DeleteInstance();

    return 0;
}
