#include "Analysis/Settings.hxx"

namespace Tree2Sexaquark::Analysis {

/* default values */
std::string Settings::PathInputFile;
std::string Settings::PathOutputFile;
bool Settings::IsMC = false;
bool Settings::IsSignalMC = false;
unsigned long Settings::LimitToNEvents = 0;
unsigned int Settings::NThreads = 1;

Settings& Settings::Instance() {
    static Settings instance;
    return instance;
}

}  // namespace Tree2Sexaquark::Analysis
