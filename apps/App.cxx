#include "Utilities/Logger.hxx"
#include "Utilities/Parser.hxx"

#include "Analysis/Manager.hxx"

int main(int argc, char *argv[]) {

    Logger *LoggerInstance = Logger::GetInstance();
    Parser *ParserInstance = new Parser("Tree2Sexaquark -- Read SimpleTrees.root and create Sexaquark Analysis Results!");
    ParserInstance->Parse(argc, argv);
    if (ParserInstance->HelpOrError) return ParserInstance->ExitCode;

    AnalysisManager *ThisAnalysis = new AnalysisManager(ParserInstance->InputFile, ParserInstance->IsMC, ParserInstance->LimitToNEvents);
    ThisAnalysis->Print();
    if (!ThisAnalysis->OpenInputFile()) return 1;

    for (Long64_t evt_idx = 0; evt_idx < ThisAnalysis->GetN_Events(); evt_idx++) {
        if (!ThisAnalysis->GetEvent(evt_idx)) continue;
        if (ThisAnalysis->IsMC()) {
            ThisAnalysis->ProcessInjected();
            ThisAnalysis->ProcessMCParticles();
        }
        ThisAnalysis->ProcessTracks();
    }

    return 0;
}
