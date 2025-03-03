#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include "Utilities/CL1.hpp"

namespace Tree2Sexaquark {

class Parser {
   public:
    Parser();
    explicit Parser(std::string app_description);
    ~Parser() = default;
    int Parse(int argc, char* argv[]);
    int ExitCode;
    bool HelpOrError;

   private:
    void AddOptions();
    CLI::App CLI_APP;
};

}  // namespace Tree2Sexaquark

#endif  // T2S_PARSER_HXX
