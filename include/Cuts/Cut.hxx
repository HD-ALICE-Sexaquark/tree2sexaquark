#ifndef T2S_CUTS_HXX
#define T2S_CUTS_HXX

#include <functional>

#include "Rtypes.h"

#include "Utilities/Logger.hxx"

#include "Particles/V0.hxx"

namespace Tree2Sexaquark {

/*
 * Container for a cut expression and its limits.
 */
class Cut {
   public:
    enum Option_t {  //
        kMinimum,
        kMaximum,
        kAbsoluteMax
    };

    Cut() = default;
    ~Cut() = default; /* This constructor requires expression and limits to be set later */

    /*
     * Set expression, but limits need to be set later
     */
    Cut(std::function<Float_t(Particle::V0)> expression)
        : fExpression(expression),  //
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {}

    /*
     * Set expression, min and max
     */
    Cut(std::function<Float_t(Particle::V0)> expression, Float_t min, Float_t max)
        : fExpression(expression),  //
          fMinimum(min),
          fMinDefined(kTRUE),
          fMaximum(max),
          fMaxDefined(kTRUE) {}

    /*
     * Sets expression and single value, which can be min., max. or abs. max.
     */
    Cut(std::function<Float_t(Particle::V0)> expression, Float_t value, Option_t option)
        : fExpression(expression),  //
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {
        //
        SetLimit(value, option);
    }

    void SetExpression(std::function<Float_t(Particle::V0)> expression) { fExpression = expression; }
    void SetLimit(Float_t value, Option_t option) {
        if (option == kAbsoluteMax) SetAbsMax(value);
        if (option == kMinimum) SetMin(value);
        if (option == kMaximum) SetMax(value);
    }
    void SetMinMax(Float_t min, Float_t max) {
        SetMin(min);
        SetMax(max);
    }
    void SetAbsMax(Float_t value) {
        SetMin(-value);
        SetMax(value);
    }
    void SetMin(Float_t value) {
        fMinimum = value;
        fMinDefined = kTRUE;
    }
    void SetMax(Float_t value) {
        fMaximum = value;
        fMaxDefined = kTRUE;
    }

    Bool_t Check(Particle::V0 thisV0) {
        if (!fExpression) {
            ErrorF("Expression not defined! %s", "");
            return kFALSE;
        }
        if (!fMinDefined && !fMaxDefined) {
            ErrorF("No limits defined! %s", "");
            return kFALSE;
        }
        /*  */
        Float_t eval = fExpression(thisV0);
        Bool_t result = kTRUE;
        if (fMinDefined) result = result && eval > fMinimum;
        if (fMaxDefined) result = result && eval < fMaximum;
        return result;
    }

   private:
    std::function<Float_t(Particle::V0)> fExpression;
    Float_t fMinimum;
    Bool_t fMinDefined;
    Float_t fMaximum;
    Bool_t fMaxDefined;
};

}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_HXX
