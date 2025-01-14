#include "Utilities/Parser.hxx"

/*
 * Default constructor
 */
Parser::Parser()  //
    : InputFile(""), IsMC(true), LimitToNEvents(0), CLI_APP(""), ExitCode(0), HelpOrError(false) {
    //
    AddOptions();
}

/*
 * Default constructor (with arguments)
 */
Parser::Parser(std::string app_description)  //
    : InputFile(""), IsMC(true), LimitToNEvents(0), CLI_APP(app_description), ExitCode(0), HelpOrError(false) {
    //
    AddOptions();
}

/*
 *
 */
void Parser::AddOptions() {

    CLI_APP
        .add_option("-f,--file", InputFile, "Path of input file")  //
        ->required()
        ->check(CLI::ExistingPath);

    CLI_APP.add_flag("--mc", IsMC, "Flag to process MC");

    CLI_APP
        .add_option("-n,--nevents", LimitToNEvents, "Limit to N events")  //
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
