#ifndef T2S_CUTS_INSPECTOR_HXX
#define T2S_CUTS_INSPECTOR_HXX

#include <vector>

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"

#include "Utilities/Logger.hxx"

#include "Cuts/Cut.hxx"
// #include "Cuts/Default.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*
 * Store and apply cuts.
 */
class Inspector {
   public:
    Inspector() = default;
    ~Inspector() = default;

    void Init();
    inline void AddCut(Cut cut) { fCutsCollection.push_back(cut); }
    Bool_t Approve(Particle::V0& thisV0);

   private:
    std::vector<Cut> fCutsCollection;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_INSPECTOR_HXX
