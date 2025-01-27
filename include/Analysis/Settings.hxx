#ifndef T2S_ANALYSIS_SETTINGS_HXX
#define T2S_ANALYSIS_SETTINGS_HXX

#include <string>

struct Settings_tt {
    std::string PathInputFile;
    std::string PathOutputFile;
    bool IsMC;
    bool IsSignalMC;
    long long LimitToNEvents;
};

#endif  // T2S_ANALYSIS_SETTINGS_HXX
