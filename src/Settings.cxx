#include "Analysis/Settings.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

// Define the static member variables
std::string Settings::PathInputFile = "";
std::string Settings::PathOutputFile = "";
bool Settings::IsMC = false;
bool Settings::IsSignalMC = false;
long long Settings::LimitToNEvents = 0;

Settings* Settings::Instance = nullptr;

/*
 *
 */
Settings* Settings::GetInstance() {
    if (!Instance) Instance = new Settings();
    return Instance;
}

/*
 *
 */
void Settings::DeleteInstance() {
    if (!Instance) {
        delete Instance;
        Instance = nullptr;
    }
}

}  // namespace Analysis
}  // namespace Tree2Sexaquark
