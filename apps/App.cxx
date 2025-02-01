#include "Analysis/Manager.hxx"
#include "Utilities/Parser.hxx"

using namespace Tree2Sexaquark;

int main(int argc, char *argv[]) {

    Logger *LoggerInstance = Logger::GetInstance();
    Analysis::Settings *SettingsInstance = Analysis::Settings::GetInstance();

    std::unique_ptr<Parser> ParserInstance =
        std::make_unique<Parser>("Tree2Sexaquark -- Read SimpleTrees.root and create Sexaquark Analysis Results!");
    ParserInstance->Parse(argc, argv);
    if (ParserInstance->HelpOrError) return ParserInstance->ExitCode;

    Analysis::Settings::Print();

    std::unique_ptr<Analysis::Manager> ThisAnalysis = std::make_unique<Analysis::Manager>();
    if (!ThisAnalysis->OpenInputFile()) return 1;
    if (!ThisAnalysis->PrepareOutputFile()) return 1;

    for (Long64_t evt_entry = 0; evt_entry < ThisAnalysis->GetN_Events(); evt_entry++) {
        if (!ThisAnalysis->GetEvent(evt_entry)) continue;
        if (Analysis::Settings::IsMC) {
            if (Analysis::Settings::IsSignalMC) ThisAnalysis->ProcessInjected();
            ThisAnalysis->ProcessMCParticles();
        }
        /* Tracks */
        ThisAnalysis->ProcessTracks();
        /* Findables */
        if (Analysis::Settings::IsMC) {
            ThisAnalysis->ProcessFindableV0s();
            if (Analysis::Settings::IsSignalMC) ThisAnalysis->ProcessFindableSexaquarks();
        }
        /* V0s */
        ThisAnalysis->KalmanV0Finder(-2212, 211, -3122);
        // ThisAnalysis->KalmanV0Finder(-211, 2212, 3122);
        ThisAnalysis->KalmanV0Finder(-211, 211, 310);
        // ThisAnalysis->KalmanV0Finder(-211, 211, 422);
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
    }  // end of loop over events
    ThisAnalysis->EndOfAnalysis();

    Analysis::Settings::DeleteInstance();
    Logger::DeleteInstance();

    return 0;
}
