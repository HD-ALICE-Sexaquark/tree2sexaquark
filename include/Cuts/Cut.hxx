#ifndef T2S_CUTS_HXX
#define T2S_CUTS_HXX

#include <functional>

#include "Rtypes.h"

#include "Utilities/Constants.hxx"
#include "Utilities/Logger.hxx"

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*
 * Container for a cut expression and its limits.
 */
template <typename T>
class Cut {
   public:
    Cut() = default; /* This constructor requires expression and limits to be set later */
    ~Cut() = default;

    /*
     * Set expression, but limits need to be set later
     */
    Cut(TString name, std::function<Float_t(T)> expression)
        : fName(name),  //
          fExpression(expression),
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {}

    /*
     * Set expression, min and max
     */
    Cut(TString name, std::function<Float_t(T)> expression, Float_t min, Float_t max)
        : fName(name),  //
          fExpression(expression),
          fMinimum(min),
          fMinDefined(kTRUE),
          fMaximum(max),
          fMaxDefined(kTRUE) {}

    /*
     * Sets expression and single value, which can be min., max. or abs. max.
     */
    Cut(TString name, std::function<Float_t(T)> expression, Float_t value, Option_t option)
        : fName(name),  //
          fExpression(expression),
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {
        //
        SetLimit(value, option);
    }

    void SetExpression(std::function<Float_t(T)> expression) { fExpression = expression; }

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

    Bool_t Check(T particle) const {
        //
        if (!fExpression) {
            ErrorF("Expression not defined! %s", "");
            return kFALSE;
        }
        if (!fMinDefined && !fMaxDefined) {
            ErrorF("No limits defined! %s", "");
            return kFALSE;
        }
        /*  */
        Float_t eval = fExpression(particle);
        Bool_t result = kTRUE;
        if (fMinDefined) result = result && eval > fMinimum;
        if (fMaxDefined) result = result && eval < fMaximum;
        return result;
    }

    void Print() {
        //
        if (fMinDefined && fMaxDefined)
            InfoF("%s BETWEEN %f AND %f", fName.Data(), fMinimum, fMaximum);
        else if (fMinDefined)
            InfoF("%s GREATER THAN %f", fName.Data(), fMinimum);
        else if (fMaxDefined)
            InfoF("%s LESS THAN %f", fName.Data(), fMaximum);
        else
            ErrorF("No limits defined! %s", "");
    }

   private:
    TString fName;
    std::function<Float_t(T)> fExpression;
    Float_t fMinimum;
    Bool_t fMinDefined;
    Float_t fMaximum;
    Bool_t fMaxDefined;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_HXX
