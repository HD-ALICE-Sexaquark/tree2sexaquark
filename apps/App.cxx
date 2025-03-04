#include "Analysis/Manager.hxx"
#include "Analysis/Settings.hxx"
#include "Utilities/Parser.hxx"

using namespace Tree2Sexaquark;

int main(int argc, char *argv[]) {

    Logger *LoggerInstance = Logger::GetInstance();

    std::unique_ptr<Parser> ParserInstance =
        std::make_unique<Parser>("Tree2Sexaquark -- Read AnalysisResults.root and create Sexaquark Analysis Results!");
    ParserInstance->Parse(argc, argv);
    if (ParserInstance->HelpOrError) return ParserInstance->ExitCode;

    Analysis::Settings::Print();

    std::unique_ptr<Analysis::Manager> ThisAnalysis = std::make_unique<Analysis::Manager>();
    ThisAnalysis->Init();
    /* Open Input File */
    RDataFrame DF_Input("Events", Analysis::Settings::PathInputFile);
    /* Start */
    RNode DF_Main = DF_Input;
    if (Analysis::Settings::LimitToNEvents > 0) DF_Main = DF_Input.Range(Analysis::Settings::LimitToNEvents);
    /* Events */
    DF_Main = ThisAnalysis->ProcessEvent(DF_Main);
    /* MC */
    if (Analysis::Settings::IsMC) {
        DF_Main = ThisAnalysis->ProcessMCParticles(DF_Main);
        if (Analysis::Settings::IsSignalMC) DF_Main = ThisAnalysis->ProcessInjected(DF_Main);
    }
    /* Tracks */
    DF_Main = ThisAnalysis->ProcessTracks(DF_Main);
    /* V0s */
    DF_Main = ThisAnalysis->FindV0s(DF_Main, PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
    DF_Main = ThisAnalysis->FindV0s(DF_Main, PdgCode::Lambda, PdgCode::PiPlus, PdgCode::Proton);
    DF_Main = ThisAnalysis->FindV0s(DF_Main, PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);
    /* Sexaquarks */
    DF_Main = ThisAnalysis->FindSexaquarks(DF_Main, PdgCode::Neutron,
                                           {PdgCode::AntiLambda, PdgCode::KaonZeroShort});  // `AntiSexaquark,Neutron -> AntiLambda,K0S`
    DF_Main = ThisAnalysis->FindSexaquarks(DF_Main, PdgCode::AntiNeutron, {PdgCode::Lambda, PdgCode::KaonZeroShort});  // `X -> Lambda,K0S`
    /*
    ThisAnalysis->KalmanSexaquarkFinder(2212, {-3122, 321, -211, 211});   // `AntiSexaquark,Proton -> AntiLambda,K+,(pi-,pi+)`
    ThisAnalysis->KalmanSexaquarkFinder(-2212, {3122, -321, -211, 211});  // `Sexaquark,AntiProton -> Lambda,K-,(pi-,pi+)`
    ThisAnalysis->KalmanSexaquarkFinder(2212, {321, 321});                // `AntiSexaquark,Proton -> K+,K+,X`
    ThisAnalysis->KalmanSexaquarkFinder(-2212, {-321, -321});             // `Sexaquark,AntiProton -> K-,K-,X`
    */
    // ThisAnalysis->PrintAll(DF_Main);       // DEBUG
    ThisAnalysis->EndOfAnalysis(DF_Main);

    Logger::DeleteInstance();

    return 0;
}
