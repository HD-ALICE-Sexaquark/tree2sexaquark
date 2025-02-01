#ifndef T2S_CUTS_INSPECTOR_HXX
#define T2S_CUTS_INSPECTOR_HXX

#include <vector>

#include "Utilities/Logger.hxx"

#include "Cuts/Cut.hxx"
#include "Particles/V0.hxx"

#include "Particles/ChannelA.hxx"
#include "Particles/ChannelD.hxx"
#include "Particles/ChannelE.hxx"
#include "Particles/KaonPair.hxx"

namespace Tree2Sexaquark {
namespace Cuts {

/*
 * Store and apply cuts.
 */
class Inspector {
   public:
    Inspector() = default;
    ~Inspector() = default;

    /* V0s */

    void SetLambdaDefaultCuts();
    void SetKaonZeroDefaultCuts();
    void SetPionPairDefaultCuts();

    void AddCut(V0::Species species, TString name, Particle::V0::MemFn expression, Double_t value, Limit limit_type) {
        Cut<Particle::V0> cut(name, expression, value, limit_type);
        if (species == V0::Species::Lambda) fCuts_Lambda.push_back(cut);
        if (species == V0::Species::KaonZeroShort) fCuts_KaonZero.push_back(cut);
        if (species == V0::Species::PionPair) fCuts_PionPair.push_back(cut);
    }

    void AddCut(V0::Species species, TString name, Particle::V0::MemFn expression, Double_t min, Double_t max) {
        Cut<Particle::V0> cut(name, expression, min, max);
        if (species == V0::Species::Lambda) fCuts_Lambda.push_back(cut);
        if (species == V0::Species::KaonZeroShort) fCuts_KaonZero.push_back(cut);
        if (species == V0::Species::PionPair) fCuts_PionPair.push_back(cut);
    }

    Bool_t Approve(Particle::V0& candidate) {
        std::vector<Cut<Particle::V0>> CutsCollection;
        if (TMath::Abs(candidate.PdgCode) == 3122)
            CutsCollection = fCuts_Lambda;
        else if (TMath::Abs(candidate.PdgCode) == 310)
            CutsCollection = fCuts_KaonZero;
        else if (TMath::Abs(candidate.PdgCode) == 422)
            CutsCollection = fCuts_PionPair;
        return ApproveParticle(candidate, CutsCollection);
    }

    /* Sexaquark Candidates */

    void SetSexaquarkDefaultCuts_ChannelA();
    void SetSexaquarkDefaultCuts_ChannelD();
    void SetSexaquarkDefaultCuts_ChannelE();
    void SetKaonPairDefaultCuts();

    template <typename T>
    void AddCut(TString name, T expression, Double_t value, Limit limit_type) {}

    template <typename T>
    void AddCut(TString name, T expression, Double_t min, Double_t max) {}

    Bool_t Approve(Sexaquark::ChannelA& candidate) { return ApproveParticle(candidate, fCuts_ChannelA); }
    Bool_t Approve(Sexaquark::ChannelD& candidate) { return ApproveParticle(candidate, fCuts_ChannelD); }
    Bool_t Approve(Sexaquark::ChannelE& candidate) { return ApproveParticle(candidate, fCuts_ChannelE); }
    Bool_t Approve(Sexaquark::KaonPair& candidate) { return ApproveParticle(candidate, fCuts_KaonPair); }

    void PrintAllCuts();

   private:
    template <typename Particle>
    Bool_t ApproveParticle(Particle& candidate, std::vector<Cut<Particle>>& CutsCollection) {
        for (auto cut : CutsCollection) {
            if (!cut.Check(candidate)) {
                // InfoF("Particle rejected. It didn't satisfy: %s (=%f)", cut.GetInfo().Data(), cut.Evaluate(candidate));
                return kFALSE;
            }
        }
        return kTRUE;
    }

    std::vector<Cut<Particle::V0>> fCuts_Lambda;
    std::vector<Cut<Particle::V0>> fCuts_KaonZero;
    std::vector<Cut<Particle::V0>> fCuts_PionPair;

    std::vector<Cut<Sexaquark::ChannelA>> fCuts_ChannelA;
    std::vector<Cut<Sexaquark::ChannelD>> fCuts_ChannelD;
    std::vector<Cut<Sexaquark::ChannelE>> fCuts_ChannelE;
    std::vector<Cut<Sexaquark::KaonPair>> fCuts_KaonPair;
};

}  // namespace Cuts
}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_INSPECTOR_HXX
