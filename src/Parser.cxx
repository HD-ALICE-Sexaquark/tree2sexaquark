#include "Utilities/Parser.hxx"

#include "Analysis/Settings.hxx"
#include "Utilities/CL1.hpp"

namespace Tree2Sexaquark {

/*
 * Default constructor
 */
Parser::Parser()  //
    : ExitCode(0), HelpOrError(false), CLI_APP("") {
    //
    AddOptions();
}

/*
 * Default constructor (with arguments)
 */
Parser::Parser(std::string app_description)  //
    : ExitCode(0), HelpOrError(false), CLI_APP(app_description) {
    //
    AddOptions();
}

/*
 *
 */
void Parser::AddOptions() {
    //
    CLI_APP
        .add_option("-i,--input", Analysis::Settings::PathInputFile, "Path of input file")  //
        ->required()
        ->check(CLI::ExistingFile);
    CLI_APP.add_option("-o,--output", Analysis::Settings::PathOutputFile, "Path of output file");
    CLI_APP.add_flag("-m,--mc", Analysis::Settings::IsMC, "Flag to process MC");
    CLI_APP
        .add_flag("-s,--signal", Analysis::Settings::IsSignalMC, "Flag to process Signal MC")  //
        ->needs("-m");
    CLI_APP
        .add_option("-n,--nevents", Analysis::Settings::LimitToNEvents, "Limit to N events")  //
        ->check(CLI::NonNegativeNumber & CLI::TypeValidator<unsigned long>());
    CLI_APP
        .add_option("-j,--nthreads", Analysis::Settings::NThreads, "Number of threads")  //
        ->check(CLI::PositiveNumber & CLI::TypeValidator<unsigned int>())
        ->excludes("-n");
}

/*
 *
 */
int Parser::Parse(int argc, char* argv[]) {
    argv = CLI_APP.ensure_utf8(argv);
    try {
        CLI_APP.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        ExitCode = e.get_exit_code();
        HelpOrError = e.get_name() == "CallForHelp" || ExitCode;
        return CLI_APP.exit(e);
    }
    /* Default options */
    if (Analysis::Settings::PathOutputFile.empty()) {
        if (Analysis::Settings::IsMC) {
            if (Analysis::Settings::IsSignalMC)
                Analysis::Settings::PathOutputFile = "./SexaquarkResults_SignalMC.root";
            else
                Analysis::Settings::PathOutputFile = "./SexaquarkResults_MC.root";
        } else {
            Analysis::Settings::PathOutputFile = "./SexaquarkResults_Data.root";
        }
    }
    return 0;
}

}  // namespace Tree2Sexaquark
