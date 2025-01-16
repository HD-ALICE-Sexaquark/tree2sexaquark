#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include "Utilities/CL1.hpp"

#include "Analysis/Settings.hxx"

/*
 *
 */
class Parser {
   public:
    Parser();
    Parser(std::string app_description);
    ~Parser() = default;
    int Parse(int argc, char* argv[]);
    int ExitCode;
    bool HelpOrError;
    Settings_tt GetSettings() { return Settings; }

   private:
    void AddOptions();
    CLI::App CLI_APP;
    Settings_tt Settings;
};

#endif  // T2S_PARSER_HXX
