#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include "Utilities/CL1.hpp"

class Parser {
   public:
    Parser();
    Parser(std::string app_description);
    ~Parser() = default;
    int Parse(int argc, char* argv[]);
    int ExitCode;
    bool HelpOrError;

    std::string InputFile;
    bool IsMC;
    int LimitToNEvents;

   private:
    void AddOptions();
    CLI::App CLI_APP;
};

#endif  // T2S_PARSER_HXX
