#ifndef T2S_ANALYSIS_SETTINGS_HXX
#define T2S_ANALYSIS_SETTINGS_HXX

#include <string>

#include "Utilities/Logger.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

/*
 * Singleton.
 */
class Settings {
   public:
    Settings(Settings &other) = delete;
    void operator=(const Settings &) = delete;
    static Settings *GetInstance();
    static void DeleteInstance();

    static void Print() {
        InfoF("IsMC           = %i", IsMC);
        InfoF("IsSignalMC     = %i", IsSignalMC);
        InfoF("InputFile      = %s", PathInputFile.c_str());
        InfoF("OutputFile     = %s", PathOutputFile.c_str());
        InfoF("LimitToNEvents = %lu", LimitToNEvents);
        InfoF("NThreads       = %u", NThreads);
    }
    static std::string PathInputFile;
    static std::string PathOutputFile;
    static bool IsMC;
    static bool IsSignalMC;
    static unsigned long LimitToNEvents;
    static unsigned int NThreads;

   private:
    Settings() = default;
    ~Settings() = default;
    static Settings *Instance;
};

}  // namespace Analysis
}  // namespace Tree2Sexaquark

#endif  // T2S_ANALYSIS_SETTINGS_HXX
