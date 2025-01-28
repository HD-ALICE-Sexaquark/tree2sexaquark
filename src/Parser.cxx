#include "Utilities/Parser.hxx"

#include "Analysis/Settings.hxx"

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
        ->check(CLI::ExistingPath);
    Analysis::Settings::PathOutputFile = "./SexaquarkResults.root";  // default
    CLI_APP.add_option("-o,--output", Analysis::Settings::PathOutputFile, "Path of output file");
    CLI_APP.add_flag("-m,--mc", Analysis::Settings::IsMC, "Flag to process MC");
    CLI_APP.add_flag("-s,--signal", Analysis::Settings::IsSignalMC, "Flag to process Signal MC");
    CLI_APP
        .add_option("-n,--nevents", Analysis::Settings::LimitToNEvents, "Limit to N events")  //
        ->check(CLI::NonNegativeNumber & CLI::TypeValidator<unsigned int>());
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
    return 0;
}

}  // namespace Tree2Sexaquark
