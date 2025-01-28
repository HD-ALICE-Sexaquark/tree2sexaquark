#ifndef T2S_CUTS_HXX
#define T2S_CUTS_HXX

#include <functional>

#include "Rtypes.h"

#include "Utilities/Logger.hxx"

#include "Particles/V0.hxx"

namespace Tree2Sexaquark {

class Cut {
   public:
    Cut() = default;
    ~Cut() = default;

    void SetExpression(std::function<Float_t(Particle::V0)> expression) { fExpression = expression; }
    void SetMaximum(Float_t value) {
        fMaxValue = value;
        fIsMaximumSet = kTRUE;
    }
    void SetMinimum(Float_t value) {
        fMinValue = value;
        fIsMinimumSet = kTRUE;
    }
    Bool_t LessThanMaximum(Particle::V0 thisV0) { return fExpression(thisV0) < fMaxValue; }
    Bool_t MoreThanMinimum(Particle::V0 thisV0) { return fExpression(thisV0) > fMinValue; }
    Bool_t Check(Particle::V0 thisV0) {
        InfoF("Expression result: %f -- Min Value: %f -- Max Value: %f",  //
              fExpression(thisV0), fIsMinimumSet ? fMinValue : -99999., fIsMaximumSet ? fMaxValue : 99999.);
        if (fIsMaximumSet && fIsMinimumSet) return LessThanMaximum(thisV0) && MoreThanMinimum(thisV0);
        if (fIsMaximumSet) return LessThanMaximum(thisV0);
        return MoreThanMinimum(thisV0);
    }

   private:
    Float_t fMaxValue;
    Bool_t fIsMaximumSet;
    Float_t fMinValue;
    Bool_t fIsMinimumSet;
    std::function<Float_t(Particle::V0)> fExpression;
};

}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_HXX
