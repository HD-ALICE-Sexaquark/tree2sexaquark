#ifndef T2S_CUTS_INSPECTOR_HXX
#define T2S_CUTS_INSPECTOR_HXX

#include <vector>

#ifndef HomogeneousField
#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#endif
#include "KFParticle.h"

#include "Utilities/Constants.hxx"
#include "Utilities/Logger.hxx"

#include "Cuts/Cut.hxx"
#include "Particles/Sexaquark.hxx"
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

    void SetLambdaDefaultCuts();
    void SetKaonZeroDefaultCuts();
    void SetPionPairDefaultCuts();

    void SetSexaquarkDefaultCuts_ChannelA();
    void SetSexaquarkDefaultCuts_ChannelD();
    void SetSexaquarkDefaultCuts_ChannelE();
    void SetKaonPairDefaultCuts();

    void AddCut(Species_t species, TString name, std::function<Float_t(Particle::V0)> expression, Float_t value, Option_t option) {
        //
        Cut<Particle::V0> cut(name, expression, value, option);
        if (species == kLambda) fCuts_Lambda.push_back(cut);
        if (species == kKaonZero) fCuts_KaonZero.push_back(cut);
        if (species == kPionPair) fCuts_PionPair.push_back(cut);
    }

    void AddCut(Species_t species, TString name, std::function<Float_t(Particle::V0)> expression, Float_t min, Float_t max) {
        //
        Cut<Particle::V0> cut(name, expression, min, max);
        if (species == kLambda) fCuts_Lambda.push_back(cut);
        if (species == kKaonZero) fCuts_KaonZero.push_back(cut);
        if (species == kPionPair) fCuts_PionPair.push_back(cut);
    }

    void AddCut(Channel_t channel, TString name, std::function<Float_t(Particle::Sexaquark)> expression, Float_t value, Option_t option) {
        //
        Cut<Particle::Sexaquark> cut(name, expression, value, option);
        if (channel == kChannelA) fCuts_ChannelA.push_back(cut);
        if (channel == kChannelD) fCuts_ChannelD.push_back(cut);
        if (channel == kChannelE) fCuts_ChannelE.push_back(cut);
        if (channel == kChannelH) fCuts_KaonPair.push_back(cut);
    }

    Bool_t Approve(Particle::V0& thisV0);
    Bool_t Approve(Particle::Sexaquark& candidate);

    void PrintAllCuts();

   private:
    std::vector<Cut<Particle::V0>> fCuts_Lambda;
    std::vector<Cut<Particle::V0>> fCuts_KaonZero;
    std::vector<Cut<Particle::V0>> fCuts_PionPair;

    std::vector<Cut<Particle::Sexaquark>> fCuts_ChannelA;
    std::vector<Cut<Particle::Sexaquark>> fCuts_ChannelD;
    std::vector<Cut<Particle::Sexaquark>> fCuts_ChannelE;
    std::vector<Cut<Particle::Sexaquark>> fCuts_KaonPair;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_INSPECTOR_HXX
