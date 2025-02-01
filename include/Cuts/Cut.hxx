#ifndef T2S_CUTS_HXX
#define T2S_CUTS_HXX

#include <functional>

#include "Rtypes.h"

#include "Utilities/Logger.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

enum class Limit {  //
    Minimum,
    Maximum,
    AbsoluteMax
};

/*
 * Container for a cut expression and its limits.
 */
template <typename Particle>
class Cut {
   public:
    using MemFn = typename Particle::MemFn;
    Cut() = default; /* This constructor requires expression and limits to be set later */
    ~Cut() = default;

    /*
     * Set name, but expression and limits need to be set later
     */
    Cut(TString name)
        : fName(name),  //
          fExpression(nullptr),
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {}

    /*
     * Set expression, but limits need to be set later
     */
    Cut(TString name, MemFn expression)
        : fName(name),  //
          fExpression(expression),
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {}

    /*
     * Set expression, min and max
     */
    Cut(TString name, MemFn expression, Double_t min, Double_t max)
        : fName(name),  //
          fExpression(expression),
          fMinimum(min),
          fMinDefined(kTRUE),
          fMaximum(max),
          fMaxDefined(kTRUE) {}

    /*
     * Sets expression and single value, which can be min, max or abs max
     */
    Cut(TString name, MemFn expression, Double_t value, Limit limit_type)
        : fName(name),  //
          fExpression(expression),
          fMinimum(0),
          fMinDefined(kFALSE),
          fMaximum(0),
          fMaxDefined(kFALSE) {
        //
        SetLimit(value, limit_type);
    }

    void SetName(TString name) { fName = name; }
    void SetExpression(MemFn expression) { fExpression = expression; }

    void SetLimit(Double_t value, Limit limit_type) {
        if (limit_type == Limit::AbsoluteMax) SetAbsMax(value);
        if (limit_type == Limit::Minimum) SetMin(value);
        if (limit_type == Limit::Maximum) SetMax(value);
    }

    void SetMinMax(Double_t min, Double_t max) {
        SetMin(min);
        SetMax(max);
    }

    void SetAbsMax(Double_t value) {
        SetMin(-value);
        SetMax(value);
    }

    void SetMin(Double_t value) {
        fMinimum = value;
        fMinDefined = kTRUE;
    }

    void SetMax(Double_t value) {
        fMaximum = value;
        fMaxDefined = kTRUE;
    }

    Double_t Evaluate(Particle particle) { return std::invoke(fExpression, particle); }

    Bool_t Check(Particle particle) const {
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
        Double_t eval = std::invoke(fExpression, particle);
        Bool_t result = kTRUE;
        if (fMinDefined) result = result && eval > fMinimum;
        if (fMaxDefined) result = result && eval < fMaximum;
        return result;
    }

    TString GetInfo() {
        //
        TString info;
        if (fMinDefined && fMaxDefined)
            info = Form("%s BETWEEN %f AND %f", fName.Data(), fMinimum, fMaximum);
        else if (fMinDefined)
            info = Form("%s GREATER THAN %f", fName.Data(), fMinimum);
        else if (fMaxDefined)
            info = Form("%s LESS THAN %f", fName.Data(), fMaximum);
        else
            info = "NO LIMITS DEFINED!";
        return info;
    }
    void Print() { InfoF("%s", GetInfo().Data()); }

   private:
    TString fName;
    MemFn fExpression;
    Double_t fMinimum;
    Bool_t fMinDefined;
    Double_t fMaximum;
    Bool_t fMaxDefined;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_HXX
