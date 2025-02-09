#include "Analysis/Settings.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

/* default values */
std::string Settings::PathInputFile = "";
std::string Settings::PathOutputFile = "";
bool Settings::IsMC = false;
bool Settings::IsSignalMC = false;
unsigned long Settings::LimitToNEvents = 0;
unsigned int Settings::NThreads = 1;

Settings* Settings::Instance = nullptr;

Settings* Settings::GetInstance() {
    if (!Instance) Instance = new Settings();
    return Instance;
}

void Settings::DeleteInstance() {
    if (!Instance) {
        delete Instance;
        Instance = nullptr;
    }
}

}  // namespace Analysis
}  // namespace Tree2Sexaquark
