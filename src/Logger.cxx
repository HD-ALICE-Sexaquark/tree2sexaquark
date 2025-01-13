#include "Utilities/Logger.hxx"

using std::ios;
using std::ofstream;
using std::ostream;

// implementation of a singleton here
Logger* Logger::Instance = NULL;
Bool_t Logger::DebugEnabled = kTRUE;

/*
 * Get Logger singleton instance
 */
Logger* Logger::GetInstance() {
    if (Instance == NULL) Instance = new Logger;
    return Instance;
}

/*
 * Delete Logger singleton instance
 */
void Logger::DeleteInstance() {
    if (Instance != NULL) {
        delete Instance;
        Instance = NULL;
    }
}

/*
 * Default private constructor
 */
Logger::Logger() {

    for (Int_t iType = kError; iType < kMaxType; iType++) {
        fOutputTypes[iType] = 0;
        fFileNames[iType] = "";
        fOutputFiles[iType] = NULL;
        fOutputStreams[iType] = NULL;

        fPrintType[iType] = kTRUE;
        fPrintModule[iType] = kFALSE;
        fPrintScope[iType] = kTRUE;
        fPrintLocation[iType] = (iType == kDebug);
    }

    // replace previous instance by this one
    if (Instance) delete Instance;
    Instance = this;
}

/*
 * Private destructor
 */
Logger::~Logger() {
    for (Int_t iType = kError; iType < kMaxType; iType++) {
        CloseFile(iType);
    }

    fflush(stderr);
    fflush(stdout);

    Instance = NULL;
}

/*
 * General method to print a log message using variadac args
 * to the FILE* like (C-like) streams, e.g. stdout, stderr, or files opened by fopen.
 * Only in case of an external C++ ostream type output, the message is
 * written to that stream and the notification callback is called.
 * The message is printed by a normal vfprintf function otherwise
 */
void Logger::PrintString(Int_t type, FILE* stream, const char* format, ...) {

    if (format == NULL) return;

    va_list ap;
    va_start(ap, format);
    if (fOutputTypes[type] != 3) {
        if (stream != NULL) {
            vfprintf(stream, format, ap);
        }
    } else {
        // build the string and write everthing to the corresponding ostream
        TString fmt(format);
        TArrayC tgt(fmt.Length() * 10);  // just take a number
#ifdef R__VA_COPY
        va_list bap;
        R__VA_COPY(bap, ap);
#else
#warning definition of R__VA_COPY has disappeared
#endif  // R__VA_COPY
        Int_t iResult = 0;
        while (1) {
            iResult = vsnprintf(tgt.GetArray(), tgt.GetSize(), format, ap);
            if (iResult == -1) {
                iResult = tgt.GetSize() * 2;
            } else if (iResult < tgt.GetSize()) {
                break;
            }
#ifdef R__VA_COPY
            if (iResult < 10000) {
                tgt.Set(iResult + 1);
                va_end(ap);
                R__VA_COPY(ap, bap);
            } else
#endif  // R__VA_COPY
            {
                tgt[tgt.GetSize() - 1] = 0;
                break;
            }
        }
#ifdef R__VA_COPY
        va_end(bap);
#endif  // R__VA_COPY
        if (fOutputStreams[type]) {
            *(fOutputStreams[type]) << tgt.GetArray();
        }
    }
    va_end(ap);
}

/*
 *  Close the file for the given type if needed
 */
void Logger::CloseFile(Int_t type) {

    if ((fOutputTypes[type] == 2) && fOutputFiles[type]) {
        Bool_t closeFile = kTRUE;
        for (Int_t iType = kError; iType < kMaxType; iType++) {
            if ((iType != type) && (fOutputFiles[iType] == fOutputFiles[type])) {
                closeFile = kFALSE;
            }
        }
        if (closeFile) {
            fclose(fOutputFiles[type]);
            ofstream* stream = reinterpret_cast<ofstream*>(fOutputStreams[type]);
            stream->close();
            delete fOutputStreams[type];
        }
    }
    fOutputFiles[type] = NULL;
    fOutputStreams[type] = NULL;
    fFileNames[type] = "";
    fOutputTypes[type] = 0;
}

/*
 * Get the output stream for the given type of messages
 */
FILE* Logger::GetOutputStream(Int_t type) {

    if (type > kDebug) type = kDebug;
    if (fOutputTypes[type] == 0)
        return stdout;
    else if (fOutputTypes[type] == 1)
        return stderr;
    else if (fOutputTypes[type] == 2) {
        if (!fOutputFiles[type]) {
            FILE* file = NULL;
            ostream* stream = NULL;
            if (!fFileNames[type].IsNull()) {
                for (Int_t iType = kError; iType < kMaxType; iType++) {
                    if ((iType != type) && (fFileNames[iType].CompareTo(fFileNames[type]) == 0) && fOutputFiles[iType]) {
                        file = fOutputFiles[iType];
                        stream = fOutputStreams[iType];
                        break;
                    }
                }
                if (!file) {
                    file = fopen(fFileNames[type], "a");
                    stream = new ofstream(fFileNames[type], ios::app);
                }
            }
            fOutputFiles[type] = file;
            fOutputStreams[type] = stream;
            if (!file) CloseFile(type);
        }
        if (fOutputFiles[type]) return fOutputFiles[type];
    }

    return stdout;
}

void Logger::PrintMessage(UInt_t type, const char* message, const char* pretty_function, const char* file, Int_t line) {
    FILE* stream = GetOutputStream(type);
    static const char* typeNames[kMaxType] = {"ERROR", "WARNING", "INFO", "DEBUG"};

    if (fPrintType[type]) {
        PrintString(type, stream, "%c::", typeNames[type][0]);
    }
    if (message) {
        PrintString(type, stream, "%s: %s", pretty_function, message);
    } else {
        PrintString(type, stream, "%s", pretty_function);
    }
    if (fPrintLocation[type] && file) {
        PrintString(type, stream, " (%s:%.0d)", file, line);
    }
    if (message) {
        PrintString(type, stream, "\n");
    } else {
        PrintString(type, stream, ": ");
    }
}

void Logger::Message(UInt_t type, const char* message, const char* pretty_function, const char* file, Int_t line) {
    Instance->PrintMessage(type, message, pretty_function, file, line);
}
