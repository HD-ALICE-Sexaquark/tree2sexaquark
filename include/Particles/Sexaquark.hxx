#ifndef T2S_SEXAQUARK_HXX
#define T2S_SEXAQUARK_HXX

#include "Math/Vector4D.h"

#include "Utilities/Constants.hxx"

#include "Particles/Base.hxx"

namespace Tree2Sexaquark {
namespace Particle {

/*
 *
 */
class Sexaquark : public Base {
   public:
    inline Double_t LambdaPx() { return lvLambda.Px(); }
    inline Double_t LambdaPy() { return lvLambda.Py(); }
    inline Double_t LambdaPz() { return lvLambda.Pz(); }
    void Test();

    Channel_t Channel;

   private:
    ROOT::Math::PxPyPzEVector lvLambda;
};

}  // namespace Particle
}  // namespace Tree2Sexaquark

#endif  // T2S_SEXAQUARK_HXX
