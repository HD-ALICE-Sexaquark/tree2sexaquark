#include "Utilities/Logger.hxx"
#include "Utilities/Parser.hxx"

#include "Analysis/Manager.hxx"

int main(int argc, char *argv[]) {

    Logger *LoggerInstance = Logger::GetInstance();
    Parser *ParserInstance = new Parser("Tree2Sexaquark -- Read SimpleTrees.root and create Sexaquark Analysis Results!");
    ParserInstance->Parse(argc, argv);
    if (ParserInstance->HelpOrError) return ParserInstance->ExitCode;

    AnalysisManager *ThisAnalysis = new AnalysisManager(ParserInstance->GetSettings());
    ThisAnalysis->Print();
    if (!ThisAnalysis->OpenInputFile()) return 1;

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
        /* End of Event */
        ThisAnalysis->EndOfEvent();
    }
    ThisAnalysis->EndOfAnalysis();

    return 0;
}
