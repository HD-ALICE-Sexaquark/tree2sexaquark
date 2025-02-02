#ifndef T2S_CUTS_INSPECTOR_HXX
#define T2S_CUTS_INSPECTOR_HXX

#include <unordered_map>

#include "Particles/ChannelA.hxx"
#include "Particles/ChannelD.hxx"
#include "Particles/ChannelE.hxx"
#include "Particles/KaonPair.hxx"
#include "Particles/V0.hxx"

namespace Tree2Sexaquark {
class Cut {
   public:
    enum class Limit : Short_t {  //
        Minimum,
        Maximum,
        AbsoluteMax
    };

    Cut() = default;
    ~Cut() = default;
    Cut(Double_t minimum, Double_t maximum)
        : fIsOn(kTRUE),  //
          fMinimum(minimum),
          fMinimumDefined(kTRUE),
          fMaximum(maximum),
          fMaximumDefined(kTRUE) {}
    Cut(Double_t value, Limit limit_type)
        : fIsOn(kTRUE),  //
          fMinimum(0.),
          fMinimumDefined(kFALSE),
          fMaximum(0.),
          fMaximumDefined(kFALSE) {
        ModifyTo(value, limit_type);
    }

    void ModifyTo(Double_t minimum, Double_t maximum) {
        fMinimum = minimum;
        fMinimumDefined = kTRUE;
        fMaximum = maximum;
        fMaximumDefined = kTRUE;
    }
    void ModifyTo(Double_t value, Limit limit_type) {
        if (limit_type == Limit::Minimum) {
            fMinimum = value;
            fMinimumDefined = kTRUE;
        } else if (limit_type == Limit::Maximum) {
            fMaximum = value;
            fMaximumDefined = kTRUE;
        } else if (limit_type == Limit::AbsoluteMax) {
            fMinimum = -value;
            fMinimumDefined = kTRUE;
            fMaximum = value;
            fMaximumDefined = kTRUE;
        }
    }

    inline Bool_t Apply(Double_t value) const {
        if (!fIsOn) return kTRUE;
        if (fMinimumDefined && value < fMinimum) return kFALSE;
        if (fMaximumDefined && value > fMaximum) return kFALSE;
        return kTRUE;
    }

    inline Bool_t IsOn() const { return fIsOn; }
    void TurnOff() { fIsOn = kFALSE; }

    std::string ToString() const {
        std::string str = "";
        if (fMinimumDefined) {
            str += "Min: " + std::to_string(fMinimum);
        }
        if (fMaximumDefined) {
            str += " Max: " + std::to_string(fMaximum);
        }
        return str;
    }

   private:
    Bool_t fIsOn;
    Double_t fMinimum;
    Bool_t fMinimumDefined;
    Double_t fMaximum;
    Bool_t fMaximumDefined;
};

namespace Cuts {
class Inspector {
    using V0 = Candidate::V0;
    using Sexaquark = Candidate::Sexaquark;

   public:
    Inspector() = default;
    ~Inspector() = default;

    void InitDefaultCuts_V0();
    void InitDefaultCuts_Sexaquark();

    void AddCut(V0::Species species, std::string cut_name, Double_t minimum, Double_t maximum) {
        if (species == V0::Species::Lambda) {
            Lambda_[cut_name] = Cut(minimum, maximum);
        } else if (species == V0::Species::KaonZeroShort) {
            KaonZeroShort_[cut_name] = Cut(minimum, maximum);
        } else {
            PionPair_[cut_name] = Cut(minimum, maximum);
        }
    }

    void AddCut(V0::Species species, std::string cut_name, Double_t value, Cut::Limit limit_type) {
        if (species == V0::Species::Lambda) {
            Lambda_[cut_name] = Cut(value, limit_type);
        } else if (species == V0::Species::KaonZeroShort) {
            KaonZeroShort_[cut_name] = Cut(value, limit_type);
        } else {
            PionPair_[cut_name] = Cut(value, limit_type);
        }
    }

    void AddCut(Sexaquark::Channel channel, std::string cut_name, Double_t minimum, Double_t maximum) {
        if (channel == Sexaquark::Channel::A) {
            ChannelA_[cut_name] = Cut(minimum, maximum);
        } else if (channel == Sexaquark::Channel::D) {
            ChannelD_[cut_name] = Cut(minimum, maximum);
        } else if (channel == Sexaquark::Channel::E) {
            ChannelE_[cut_name] = Cut(minimum, maximum);
        } else {
            KaonPair_[cut_name] = Cut(minimum, maximum);
        }
    }

    void AddCut(Sexaquark::Channel channel, std::string cut_name, Double_t value, Cut::Limit limit_type) {
        if (channel == Sexaquark::Channel::A) {
            ChannelA_[cut_name] = Cut(value, limit_type);
        } else if (channel == Sexaquark::Channel::D) {
            ChannelD_[cut_name] = Cut(value, limit_type);
        } else if (channel == Sexaquark::Channel::E) {
            ChannelE_[cut_name] = Cut(value, limit_type);
        } else {
            KaonPair_[cut_name] = Cut(value, limit_type);
        }
    }

    Cut GetCut(V0::Species species, std::string cut_name) {
        if (species == V0::Species::Lambda) {
            return Lambda_[cut_name];
        } else if (species == V0::Species::KaonZeroShort) {
            return KaonZeroShort_[cut_name];
        } else {
            return PionPair_[cut_name];
        }
    }

    Cut GetCut(Sexaquark::Channel channel, std::string cut_name) {
        if (channel == Sexaquark::Channel::A) {
            return ChannelA_[cut_name];
        } else if (channel == Sexaquark::Channel::D) {
            return ChannelD_[cut_name];
        } else if (channel == Sexaquark::Channel::E) {
            return ChannelE_[cut_name];
        } else {
            return KaonPair_[cut_name];
        }
    }

    inline Bool_t Check(Candidate::V0 v0, std::string cut_name, Double_t value) {
        //
        return GetCut(v0.GetSpecies(), cut_name).Apply(value);
    }
    inline Bool_t Check(Candidate::ChannelA, std::string cut_name, Double_t value) { return ChannelA_[cut_name].Apply(value); }
    inline Bool_t Check(Candidate::ChannelD, std::string cut_name, Double_t value) { return ChannelD_[cut_name].Apply(value); }
    inline Bool_t Check(Candidate::ChannelE, std::string cut_name, Double_t value) { return ChannelE_[cut_name].Apply(value); }
    inline Bool_t Check(Candidate::KaonPair, std::string cut_name, Double_t value) { return KaonPair_[cut_name].Apply(value); }

    Bool_t Approve(Candidate::V0 v0);
    Bool_t Approve(Candidate::ChannelA sexaquark);
    Bool_t Approve(Candidate::ChannelD sexaquark);
    Bool_t Approve(Candidate::ChannelE sexaquark);
    Bool_t Approve(Candidate::KaonPair sexaquark);

    void PrintAllCuts();

   private:
    std::unordered_map<std::string, Cut> Lambda_;
    std::unordered_map<std::string, Cut> KaonZeroShort_;
    std::unordered_map<std::string, Cut> PionPair_;
    std::unordered_map<std::string, Cut> ChannelA_;
    std::unordered_map<std::string, Cut> ChannelD_;
    std::unordered_map<std::string, Cut> ChannelE_;
    std::unordered_map<std::string, Cut> KaonPair_;
};
}  // namespace Cuts

}  // namespace Tree2Sexaquark

#endif  // T2S_CUTS_INSPECTOR_HXX
