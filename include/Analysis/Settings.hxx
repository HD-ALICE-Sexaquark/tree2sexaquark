#ifndef T2S_ANALYSIS_SETTINGS_HXX
#define T2S_ANALYSIS_SETTINGS_HXX

#include <iostream>
#include <string>

namespace Tree2Sexaquark::Analysis {

/*
 * Singleton.
 */
class Settings {
   public:
    Settings(const Settings& other) = delete;
    Settings& operator=(const Settings& other) = delete;

    static Settings& Instance();

    static void Print() {
        std::cout << "IsMC           = " << IsMC << '\n';
        std::cout << "IsSignalMC     = " << IsSignalMC << '\n';
        std::cout << "InputFile      = " << PathInputFile << '\n';
        std::cout << "OutputFile     = " << PathOutputFile << '\n';
        std::cout << "LimitToNEvents = " << LimitToNEvents << '\n';
        std::cout << "NThreads       = " << NThreads << '\n';
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
};

}  // namespace Tree2Sexaquark::Analysis

#endif  // T2S_ANALYSIS_SETTINGS_HXX
