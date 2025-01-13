#ifndef UTILITIES_LOGGER_HXX
#define UTILITIES_LOGGER_HXX

#include <cstdlib>
#include <fstream>
#include <iostream>

#include <Varargs.h>  // platform independent definition of va_copy
#include <strings.h>

#include "TArrayC.h"
#include "TString.h"

using std::ostream;

class Logger {
   public:
    enum EType_t { kError = 0, kWarning, kInfo, kDebug, kMaxType };
    static void Message(UInt_t level, const char* message, const char* function, const char* file, Int_t line);
    static Logger* GetInstance();
    static void DeleteInstance();

   private:
    Logger();
    virtual ~Logger();

    void CloseFile(Int_t type);
    FILE* GetOutputStream(Int_t type);

    void PrintMessage(UInt_t type, const char* message, const char* function, const char* file, Int_t line);
    void PrintString(Int_t type, FILE* stream, const char* format, ...);

    static Logger* Instance;     //! pointer to current instance
    static Bool_t DebugEnabled;  // flag for debug en-/disabling

    Int_t fOutputTypes[kMaxType];       // types of output streams
    TString fFileNames[kMaxType];       // file names
    FILE* fOutputFiles[kMaxType];       //! log output files
    ostream* fOutputStreams[kMaxType];  //! log output streams

    Bool_t fPrintType[kMaxType];      // print type on/off
    Bool_t fPrintModule[kMaxType];    // print module on/off
    Bool_t fPrintScope[kMaxType];     // print scope/class name on/off
    Bool_t fPrintLocation[kMaxType];  // print file and line on/off
};

/* Source: https://stackoverflow.com/a/15775519 */
inline const char* MethodName(const std::string& prettyFunction) {
    static std::string output;
    size_t colons = prettyFunction.find("::") != std::string::npos ? prettyFunction.rfind("::") : prettyFunction.size();
    size_t begin = prettyFunction.substr(0, colons).find(" ") + 1;
    size_t end = prettyFunction.rfind("(") - begin;
    output = prettyFunction.substr(begin, end) + "()";
    return output.c_str();
}
#define __METHOD_NAME__ MethodName(__PRETTY_FUNCTION__)

#define MessageF(type, format, ...)                                    \
    do {                                                               \
        TString m = TString::Format(format, __VA_ARGS__);              \
        Logger::Message(type, m, __METHOD_NAME__, __FILE__, __LINE__); \
    } while (false)

#define ErrorF(message, ...) MessageF(Logger::kError, message, __VA_ARGS__)
#define WarningF(message, ...) MessageF(Logger::kWarning, message, __VA_ARGS__)
#define InfoF(message, ...) MessageF(Logger::kInfo, message, __VA_ARGS__)
#define DebugF(message, ...) MessageF(Logger::kDebug, message, __VA_ARGS__)

#endif  // UTILITIES_LOGGER_HXX
