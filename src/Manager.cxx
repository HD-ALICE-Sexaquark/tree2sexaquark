#include "Analysis/Manager.hxx"

#include <set>

#include "ROOT/RResultPtr.hxx"
#include "RtypesCore.h"
#include "TROOT.h"

#include "KFParticle.h"

#include "Analysis/Settings.hxx"
#include "Cuts/Default.hxx"
#include "Math/Common.hxx"
#include "Math/KalmanFilter.hxx"

namespace Tree2Sexaquark {
namespace Analysis {

/*                   */
/* Vector Gymnastics */
/* ================= */

// mask a vector, e.g. {0,2} -> {1,0,1}
RVecI Manager::Mask(cRVecUL entries, size_t reference_size) {
    RVecI output(reference_size);
    for (auto entry : entries) output[entry] = 1;
    return output;
}

template <typename T>
RVec<T> Manager::Extract(const RVec<T> &property, cRVecUL link_) {
    return Map(link_, [&](size_t linked_entry) { return property[linked_entry]; });
}

template <typename T>
RVec<T> Manager::ExtractIf(const RVec<T> &property, cRVecL link_, cRVecI link_protection, T default_value) {
    RVec<T> output;
    output.reserve(link_protection.size());
    for (size_t i = 0; i < link_protection.size(); i++) {
        output.emplace_back(link_protection[i] ? property[link_[i]] : default_value);
    }
    return output;
}

template <typename T>
RVec<RVec<T>> Manager::ExtractVector(const RVec<T> &property, const RVec<RVecUL> &entries_) {
    return Map(entries_, [&](cRVecUL v) { return Take(property, v); });
}

template <typename T>
RVec<T> Manager::ExtractVector_First(const RVec<T> &property, const RVec<RVecUL> &entries_) {
    return Map(entries_, [&](cRVecUL v) { return property[v[0]]; });
}

template <typename T>
RVec<T> Manager::ExtractVector_Sum(const RVec<T> &property, const RVec<RVecUL> &entries_) {
    return Map(entries_, [&](cRVecUL v) { return Sum(Take(property, v)); });
}

void Manager::Init() {
    //
    ROOT::EnableThreadSafety();
    if (Settings::NThreads > 1) ROOT::EnableImplicitMT(Settings::NThreads);
}

RNode Manager::ProcessEvent(RNode df) {
    //
    /* Set Event Properties */
    /* -- Magnetic Field */
    // KFParticle::SetField(Event.MagneticField);
    /* -- Primary Vertex */
    // Float_t XYZ[3] = {Event.PV_Xv, Event.PV_Yv, Event.PV_Zv};
    // Float_t CovMatrix[6] = {Event.PV_CovMatrix[0], Event.PV_CovMatrix[1], Event.PV_CovMatrix[2],
    // Event.PV_CovMatrix[3], Event.PV_CovMatrix[4], Event.PV_CovMatrix[5]};
    // kfPrimaryVertex = Math::CreateKFVertex(XYZ, CovMatrix);
    auto CallCreateKFVertex = [](Float_t x, Float_t y, Float_t z, cRVecF cov_matrix) -> KFVertex {
        Float_t XYZ[3] = {x, y, z};
        Float_t CovMatrix[6] = {cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3], cov_matrix[4], cov_matrix[5]};
        return Math::CreateKFVertex(XYZ, CovMatrix);
    };

    RNode DF_Events = df.Define("Event_KF_PV", CallCreateKFVertex, {"PV_Xv", "PV_Yv", "PV_Zv", "PV_CovMatrix"});

    return DF_Events;
}

RNode Manager::ProcessMCParticles(RNode df) {
    //
    auto GetCharge = [](cRVecI pdg_code) -> RVecI {
        RVecI mc_charge;
        mc_charge.reserve(pdg_code.size());
        for (size_t i = 0; i < pdg_code.size(); i++) {
            auto particle = TDatabasePDG::Instance()->GetParticle(pdg_code[i]);
            if (!particle)
                mc_charge.emplace_back(0);
            else
                mc_charge.emplace_back(Int_t(particle->Charge() / 3.));  // convert units from |e|/3 to |e|
        }
        return mc_charge;
    };
    //
    auto GetDaughter = [](cRVecI mc_mask, cRVecL mother_mc_entry_) -> RVecL {
        RVecL output(mc_mask.size(), -1);
        auto daughters_entries = Nonzero(mc_mask);
        auto mothers_entries = mother_mc_entry_[mc_mask];
        for (size_t i = 0; i < mothers_entries.size(); i++) {
            output[mothers_entries[i]] = daughters_entries[i];
        }
        return output;
    };
    //
    auto GetAncestor = [](cRVecL mother_mc_entry_) -> RVecL {
        RVecL output(mother_mc_entry_.size(), -1);
        Int_t curr_entry;
        for (size_t i = 0; i < mother_mc_entry_.size(); i++) {
            // if no mother, the particle is its own ancestor
            if (mother_mc_entry_[i] < 0) {
                output[i] = i;
                continue;
            }
            // find the last ancestor
            curr_entry = mother_mc_entry_[i];
            while (curr_entry >= 0) {
                output[i] = curr_entry;
                curr_entry = mother_mc_entry_[curr_entry];
            }
        }
        return output;
    };
    auto CallExtractInt_IfNotZero = [](cRVecI property, cRVecL link_, cRVecI link_protection) -> RVecI {
        return ExtractIf<Int_t>(property, link_, link_protection, 0);
    };
    auto CallExtractUInt_IfNotZero = [](cRVecU property, cRVecL link_, cRVecI link_protection) -> RVecU {
        return ExtractIf<UInt_t>(property, link_, link_protection, 0);
    };

    RNode df_MC =
        df.Alias("MC_Mother_McEntry_", "MC_Mother_McEntry")
            .Alias("Track_McEntry_", "Track_McEntry")
            /* # Trivial Definitions */
            .Define("N_MC", "MC_Px.size()")
            .Define("MC_HasMother", "MC_Mother_McEntry_ >= 0")
            .Define("MC_IsSignal", "MC_Generator == 2")
            .Define("MC_IsFirstGenSignal", "!MC_HasMother && MC_IsSignal")
            .Define("MC_IsSecondary", "MC_IsSecFromMat || MC_IsSecFromWeak || MC_IsSignal")
            .Define("MC_IsProtonKaonPion", "abs(MC_PdgCode) == 2212 || abs(MC_PdgCode) == 321 || abs(MC_PdgCode) == 211")
            .Define("MC_Radius", "sqrt(MC_Xv * MC_Xv + MC_Yv * MC_Yv)")
            .Define("MC_Charge", GetCharge, {"MC_PdgCode"})
            /* # Mother */
            .Define("MC_Mother_PdgCode", CallExtractInt_IfNotZero, {"MC_PdgCode", "MC_Mother_McEntry_", "MC_HasMother"})
            /* # Daughters */
            .Define("MC_IsNegDau", "MC_HasMother && MC_Charge < 0")
            .Define("MC_IsPosDau", "MC_HasMother && MC_Charge > 0")
            .Define("MC_NegDau_McEntry_", GetDaughter, {"MC_IsNegDau", "MC_Mother_McEntry_"})
            .Define("MC_PosDau_McEntry_", GetDaughter, {"MC_IsPosDau", "MC_Mother_McEntry_"})
            .Define("MC_HasNegDau", "MC_NegDau_McEntry_ >= 0")
            .Define("MC_HasPosDau", "MC_PosDau_McEntry_ >= 0")
            .Define("MC_NegDau_PdgCode", CallExtractInt_IfNotZero, {"MC_PdgCode", "MC_NegDau_McEntry_", "MC_IsNegDau"})
            .Define("MC_PosDau_PdgCode", CallExtractInt_IfNotZero, {"MC_PdgCode", "MC_PosDau_McEntry_", "MC_IsPosDau"})
            /* # Ancestor */
            .Define("MC_Ancestor_McEntry_", GetAncestor, {"MC_Mother_McEntry_"})
            .Define("MC_ReactionID", CallExtractUInt_IfNotZero, {"MC_Status", "MC_Ancestor_McEntry_", "MC_IsSignal"})
            /* # Tracked */
            .Define("MC_GotTracked", Mask, {"Track_McEntry_", "N_MC"})
            .Define("MC_NegDau_GotTracked", CallExtractInt_IfNotZero, {"MC_GotTracked", "MC_NegDau_McEntry_", "MC_HasNegDau"})
            .Define("MC_PosDau_GotTracked", CallExtractInt_IfNotZero, {"MC_GotTracked", "MC_PosDau_McEntry_", "MC_HasPosDau"})
            /* # True V0s */
            .Define("MC_TrueV0_IsFindable", "(abs(MC_PdgCode) == 3122 || MC_PdgCode == 310) && MC_NegDau_GotTracked && MC_PosDau_GotTracked");

    return df_MC;
}

RNode Manager::ProcessInjected(RNode df) {
    //
    auto GetReactionProducts = [](cRVecU inj_reaction_id, cRVecU mc_reaction_id, cRVecI mc_is_firstgen_signal) -> RVec<RVecUL> {
        RVec<RVecUL> output;
        output.reserve(inj_reaction_id.size());
        RVecUL products_vec;
        products_vec.reserve(4);
        for (auto reaction_id : inj_reaction_id) {
            for (size_t mc_entry = 0; mc_entry < mc_reaction_id.size(); mc_entry++) {
                if (mc_reaction_id[mc_entry] != reaction_id || !mc_is_firstgen_signal[mc_entry]) continue;
                products_vec.emplace_back(mc_entry);
            }
            output.emplace_back(products_vec);
            products_vec.clear();
        }
        return output;
    };
    auto IsFindable_InjectedSexaquark = [](cRVecI pdg_code, cRVecI mc_got_tracked, cRVecI mc_true_v0_is_findable,
                                           const RVec<RVecUL> &injected_mc_entries_) -> RVecI {
        RVecI output;
        output.reserve(injected_mc_entries_.size());
        Bool_t all_findable, mc_is_findable;
        for (auto products_vec : injected_mc_entries_) {
            all_findable = true;
            for (auto mc_entry : products_vec) {
                if (TMath::Abs(pdg_code[mc_entry]) == 3122 || pdg_code[mc_entry] == 310)
                    mc_is_findable = mc_true_v0_is_findable[mc_entry];
                else
                    mc_is_findable = mc_got_tracked[mc_entry];
                all_findable &= mc_is_findable;
            }
            output.emplace_back(all_findable);
        }
        return output;
    };

    RNode df_Injected =  //
        df.Define("N_Injected", "Sexaquark_Px.size()")
            .Define("Injected_Pt", "sqrt(Sexaquark_Px * Sexaquark_Px + Sexaquark_Py * Sexaquark_Py)")
            /* # Vector Operations */
            .Define("Injected_Products_McEntries_", GetReactionProducts, {"ReactionID", "MC_ReactionID", "MC_IsFirstGenSignal"})
            .Define("Injected_Products_PdgCode", ExtractVector<Int_t>, {"MC_PdgCode", "Injected_Products_McEntries_"})
            .Define("Injected_Products_GotTracked", ExtractVector<Int_t>, {"MC_GotTracked", "Injected_Products_McEntries_"})
            .Define("Injected_Products_IsFindable", ExtractVector<Int_t>, {"MC_TrueV0_IsFindable", "Injected_Products_McEntries_"})
            /* # Aggregated Properties */
            .Define("Injected_Radius", ExtractVector_First<Float_t>, {"MC_Radius", "Injected_Products_McEntries_"})
            .Define("Injected_Products_SumPx", ExtractVector_Sum<Float_t>, {"MC_Px", "Injected_Products_McEntries_"})
            .Define("Injected_Products_SumPy", ExtractVector_Sum<Float_t>, {"MC_Py", "Injected_Products_McEntries_"})
            .Define("Injected_Pt_AfterReaction",
                    "sqrt(Injected_Products_SumPx * Injected_Products_SumPx + Injected_Products_SumPy * Injected_Products_SumPy)")
            .Define("Injected_IsFindable", IsFindable_InjectedSexaquark,  //
                    {"MC_PdgCode", "MC_GotTracked", "MC_TrueV0_IsFindable", "Injected_Products_McEntries_"});

    return df_Injected;
}

RNode Manager::ProcessTracks(RNode df) {
    //
    auto GetTrackCovMatrix = [](cRVecF track_sigma_y2, cRVecF track_sigma_zy, cRVecF track_sigma_z2,             //
                                cRVecF track_sigma_snp_y, cRVecF track_sigma_snp_z, cRVecF track_sigma_snp2,     //
                                cRVecF track_sigma_tgl_y, cRVecF track_sigma_tgl_z, cRVecF track_sigma_tgl_snp,  //
                                cRVecF track_sigma_tgl2, cRVecF track_sigma_1pt_y, cRVecF track_sigma_1pt_z,     //
                                cRVecF track_sigma_1pt_snp, cRVecF track_sigma_1pt_tgl, cRVecF track_sigma_1pt2) -> RVec<Float_t *> {
        RVec<Float_t *> output;
        output.reserve(track_sigma_y2.size());
        Float_t cov_matrix[15];
        for (size_t track_entry = 0; track_entry < track_sigma_y2.size(); track_entry++) {
            cov_matrix[0] = track_sigma_y2[track_entry];
            cov_matrix[1] = track_sigma_zy[track_entry];
            cov_matrix[2] = track_sigma_z2[track_entry];
            cov_matrix[3] = track_sigma_snp_y[track_entry];
            cov_matrix[4] = track_sigma_snp_z[track_entry];
            cov_matrix[5] = track_sigma_snp2[track_entry];
            cov_matrix[6] = track_sigma_tgl_y[track_entry];
            cov_matrix[7] = track_sigma_tgl_z[track_entry];
            cov_matrix[8] = track_sigma_tgl_snp[track_entry];
            cov_matrix[9] = track_sigma_tgl2[track_entry];
            cov_matrix[10] = track_sigma_1pt_y[track_entry];
            cov_matrix[11] = track_sigma_1pt_z[track_entry];
            cov_matrix[12] = track_sigma_1pt_snp[track_entry];
            cov_matrix[13] = track_sigma_1pt_tgl[track_entry];
            cov_matrix[14] = track_sigma_1pt2[track_entry];
            output.emplace_back(cov_matrix);
        }
        return output;
    };

    RNode df_Tracks =  //
        df.Define("N_Tracks", "Track_Px.size()")
            /* # PID */
            .Define("Track_PID_AsAntiProton", "abs(Track_NSigmaProton) < 3. && Track_Charge < 0")
            .Define("Track_PID_AsProton", "abs(Track_NSigmaProton) < 3. && Track_Charge > 0")
            .Define("Track_PID_AsNegKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge < 0")
            .Define("Track_PID_AsPosKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge > 0")
            .Define("Track_PID_AsPiMinus", "abs(Track_NSigmaPion) < 3. && Track_Charge < 0")
            .Define("Track_PID_AsPiPlus", "abs(Track_NSigmaPion) < 3. && Track_Charge > 0")
            .Define("AntiProton_TrackEntry_", "Nonzero(Track_PID_AsAntiProton)")
            .Define("Proton_TrackEntry_", "Nonzero(Track_PID_AsProton)")
            .Define("PiMinus_TrackEntry_", "Nonzero(Track_PID_AsPiMinus)")
            .Define("PiPlus_TrackEntry_", "Nonzero(Track_PID_AsPiPlus)")
            /* # Utilities */
            .Define("Track_CovMatrix", GetTrackCovMatrix,                                                              //
                    {"Track_SigmaY2", "Track_SigmaZY", "Track_SigmaZ2", "Track_SigmaSnpY", "Track_SigmaSnpZ",          //
                     "Track_SigmaSnp2", "Track_SigmaTglY", "Track_SigmaTglZ", "Track_SigmaTglSnp", "Track_SigmaTgl2",  //
                     "Track_Sigma1PtY", "Track_Sigma1PtZ", "Track_Sigma1PtSnp", "Track_Sigma1PtTgl", "Track_Sigma1Pt2"});

    if (Settings::IsMC) {
        df_Tracks = df_Tracks  //
                        .Define("Track_PdgCode", Extract<Int_t>, {"MC_PdgCode", "Track_McEntry"})
                        .Define("Track_IsSignal", Extract<Int_t>, {"MC_IsSignal", "Track_McEntry"});
    }

    return df_Tracks;
}

/*         */
/**  V0s  **/
/*** === ***/

/*
 * Find all V0s via Kalman Filter.
 */
RNode Manager::FindV0s(RNode df, Int_t pdg_v0, Int_t pdg_neg, Int_t pdg_pos) {
    //
    const Double_t neg_mass = TDatabasePDG::Instance()->GetParticle(pdg_neg)->Mass();
    const Double_t pos_mass = TDatabasePDG::Instance()->GetParticle(pdg_pos)->Mass();
    //
    std::string v0_name, neg_name, pos_name;
    if (pdg_v0 == 3122) {
        v0_name = "Lambdas";
        neg_name = "PiMinus";
        pos_name = "Proton";
    } else if (pdg_v0 == -3122) {
        v0_name = "AntiLambdas";
        neg_name = "AntiProton";
        pos_name = "PiPlus";
    } else if (pdg_v0 == 310) {
        v0_name = "KaonsZeroShort";
        neg_name = "PiMinus";
        pos_name = "PiPlus";
    }
    //
    auto KF_V0Finder = [pdg_v0, neg_mass, pos_mass](cRVecUL entries_neg, cRVecUL entries_pos,                                     //
                                                    cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z, RVecI charge,  //
                                                    cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt, const RVec<Float_t *> &cov_matrix,
                                                    Float_t magnetic_field, KFVertex kf_pv) -> RVec<V0_tt> {
        KFParticle::SetField(magnetic_field);
        RVec<V0_tt> output;
        V0_tt this_v0;
        ULong64_t neg, pos;
        KFParticle kf_v0, kf_neg, kf_pos;
        PxPyPzMVector lv_v0, lv_neg, lv_pos;
        XYZPoint v3_v0;
        XYZPoint v3_pv(kf_pv.GetX(), kf_pv.GetY(), kf_pv.GetZ());
        Double_t cpa_wrt_pv, dca_wrt_pv, dca_btw_dau, dca_neg_v0, dca_pos_v0;
        Double_t arm_qt, arm_alpha;
        for (auto neg : entries_neg) {
            for (auto pos : entries_pos) {
                if (neg == pos) continue;
                /* geometry */
                kf_neg = Math::CreateKFParticle(px[neg], py[neg], pz[neg], x[neg], y[neg], z[neg], charge[neg], alpha[neg], snp[neg], tgl[neg],
                                                signed1pt[neg], cov_matrix[neg], neg_mass);
                kf_pos = Math::CreateKFParticle(px[pos], py[pos], pz[pos], x[pos], y[pos], z[pos], charge[pos], alpha[pos], snp[pos], tgl[pos],
                                                signed1pt[pos], cov_matrix[pos], pos_mass);
                KFParticle kf_v0;
                kf_v0.AddDaughter(kf_neg);
                kf_v0.AddDaughter(kf_pos);
                v3_v0.SetCoordinates(kf_v0.GetX(), kf_v0.GetY(), kf_v0.GetZ());
                /* kinematics */
                lv_neg.SetCoordinates(kf_neg.Px(), kf_neg.Py(), kf_neg.Pz(), neg_mass);
                lv_pos.SetCoordinates(kf_pos.Px(), kf_pos.Py(), kf_pos.Pz(), pos_mass);
                lv_v0 = lv_neg + lv_pos;
                /* cuts */
                /* -- (anti)lambda cuts */
                if (TMath::Abs(pdg_v0) == 3122) {
                    if (TMath::Abs(lv_v0.Eta()) > Default::Lambda::AbsMaxEta) continue;
                    if (lv_v0.Pt() < Default::Lambda::MinPt) continue;
                    if (lv_v0.M() < Default::Lambda::MinMass || lv_v0.M() > Default::Lambda::MaxMass) continue;
                    if (v3_v0.Rho() < Default::Lambda::MinRadius) continue;
                    if ((v3_v0 - v3_pv).R() < Default::Lambda::MinDistFromPV) continue;
                    /*  */ cpa_wrt_pv = Math::CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
                    if (cpa_wrt_pv < Default::Lambda::MinCPAwrtPV || cpa_wrt_pv > Default::Lambda::MaxCPAwrtPV) continue;
                    /*  */ dca_wrt_pv = TMath::Abs((Double_t)kf_v0.GetDistanceFromVertex(kf_pv));
                    if (dca_wrt_pv < Default::Lambda::MinDCAwrtPV) continue;
                    /*  */ dca_btw_dau = TMath::Abs((Double_t)kf_neg.GetDistanceFromParticle(kf_pos));
                    if (dca_btw_dau > Default::Lambda::MaxDCAbtwDau) continue;
                    /*  */ dca_neg_v0 = TMath::Abs((Double_t)kf_neg.GetDistanceFromVertex(kf_v0));
                    if (dca_neg_v0 > Default::Lambda::MaxDCAnegV0) continue;
                    /*  */ dca_pos_v0 = TMath::Abs((Double_t)kf_pos.GetDistanceFromVertex(kf_v0));
                    if (dca_pos_v0 > Default::Lambda::MaxDCAposV0) continue;
                    /*  */ arm_qt = Math::ArmenterosQt(lv_v0.Vect(), lv_neg.Vect());
                    /*  */ arm_alpha = Math::ArmenterosAlpha(lv_v0.Vect(), lv_neg.Vect(), lv_pos.Vect());
                    if (arm_qt / TMath::Abs(arm_alpha) > Default::Lambda::MaxArmQtOverAlpha) continue;
                }
                /* -- k0s cuts */
                if (pdg_v0 == 310) {
                    if (TMath::Abs(lv_v0.Eta()) > Default::KaonZeroShort::AbsMaxEta) continue;
                    if (lv_v0.Pt() < Default::KaonZeroShort::MinPt) continue;
                    if (lv_v0.M() < Default::KaonZeroShort::MinMass || lv_v0.M() > Default::KaonZeroShort::MaxMass) continue;
                    if (v3_v0.Rho() < Default::KaonZeroShort::MinRadius) continue;
                    if ((v3_v0 - v3_pv).R() < Default::KaonZeroShort::MinDistFromPV ||  //
                        (v3_v0 - v3_pv).R() > Default::KaonZeroShort::MaxDistFromPV)
                        continue;
                    /*  */ cpa_wrt_pv = Math::CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
                    if (cpa_wrt_pv < Default::KaonZeroShort::MinCPAwrtPV || cpa_wrt_pv > Default::KaonZeroShort::MaxCPAwrtPV) continue;
                    /*  */ dca_wrt_pv = TMath::Abs((Double_t)kf_v0.GetDistanceFromVertex(kf_pv));
                    if (dca_wrt_pv < Default::KaonZeroShort::MinDCAwrtPV) continue;
                    /*  */ dca_btw_dau = TMath::Abs((Double_t)kf_neg.GetDistanceFromParticle(kf_pos));
                    if (dca_btw_dau > Default::KaonZeroShort::MaxDCAbtwDau) continue;
                    /*  */ dca_neg_v0 = TMath::Abs((Double_t)kf_neg.GetDistanceFromVertex(kf_v0));
                    if (dca_neg_v0 > Default::KaonZeroShort::MaxDCAnegV0) continue;
                    /*  */ dca_pos_v0 = TMath::Abs((Double_t)kf_pos.GetDistanceFromVertex(kf_v0));
                    if (dca_pos_v0 > Default::KaonZeroShort::MaxDCAposV0) continue;
                }
                /* store */
                this_v0.neg = neg;
                this_v0.pos = pos;
                this_v0.kf = kf_v0;
                this_v0.kf_neg = kf_neg;
                this_v0.kf_pos = kf_pos;
                this_v0.lv = lv_v0;
                this_v0.lv_neg = lv_neg;
                this_v0.lv_pos = lv_pos;
                output.emplace_back(this_v0);
            }  // end of loop over pos
        }      // end of loop over neg
        return output;
    };
    //
    auto GetTrackIndicesFromStruct = [](const RVec<V0_tt> &found_v0s) -> RVec<RVecUL> {
        RVec<RVecUL> output;
        output.reserve(found_v0s.size());
        RVecUL a_pair(2);
        for (auto v0 : found_v0s) {
            a_pair[0] = v0.neg;
            a_pair[1] = v0.pos;
            output.emplace_back(a_pair);
        }
        return output;
    };
    //
    auto GetTrueV0 = [pdg_v0, pdg_neg, pdg_pos](RVecI pdg_code, const RVec<RVecUL> &v0_track_entries_, cRVecUL track_mc_entry_,
                                                cRVecL mother_mc_entry_) -> RVecL {
        RVecL output;
        output.reserve(v0_track_entries_.size());
        size_t mc_neg, mc_pos;
        size_t neg_mother, pos_mother;
        for (size_t v0_entry = 0; v0_entry < v0_track_entries_.size(); v0_entry++) {
            mc_neg = track_mc_entry_[v0_track_entries_[v0_entry][0]];
            mc_pos = track_mc_entry_[v0_track_entries_[v0_entry][1]];
            neg_mother = mother_mc_entry_[mc_neg];
            pos_mother = mother_mc_entry_[mc_pos];
            if (pdg_code[mc_neg] == pdg_neg && pdg_code[mc_pos] == pdg_pos && neg_mother == pos_mother && pdg_code[neg_mother] == pdg_v0)
                output.emplace_back(neg_mother);
            else
                output.emplace_back(-1);
        }
        return output;
    };
    //
    auto IsHybrid = [](cRVecI v0_is_signal, cRVecI mc_is_signal, const RVec<RVecUL> &v0_track_entries_, cRVecUL track_mc_entry_) -> RVecI {
        RVecI output;
        output.reserve(v0_is_signal.size());
        size_t mc_neg, mc_pos;
        for (size_t v0_entry = 0; v0_entry < v0_is_signal.size(); v0_entry++) {
            mc_neg = track_mc_entry_[v0_track_entries_[v0_entry][0]];
            mc_pos = track_mc_entry_[v0_track_entries_[v0_entry][1]];
            output.emplace_back(!v0_is_signal[v0_entry] &&
                                ((mc_is_signal[mc_neg] && !mc_is_signal[mc_pos]) || (!mc_is_signal[mc_neg] && mc_is_signal[mc_pos])));
        }
        return output;
    };
    //
    auto CallExtractInt_IfNotZero = [](cRVecI property, cRVecL link_, cRVecI link_protection) -> RVecI {
        return ExtractIf<Int_t>(property, link_, link_protection, 0);
    };
    //
    auto CallExtractUInt_IfNotZero = [](cRVecU property, cRVecL link_, cRVecI link_protection) -> RVecU {
        return ExtractIf<UInt_t>(property, link_, link_protection, 0);
    };

    RNode df_V0s =
        df.Define("Found_" + v0_name, KF_V0Finder,  //
                  {neg_name + "_TrackEntry_", pos_name + "_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",
                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt", "Track_CovMatrix", "MagneticField", "Event_KF_PV"})
            .Define(v0_name + "_TrackEntries_", GetTrackIndicesFromStruct, {"Found_" + v0_name})
            .Define(v0_name + "_M", [](const RVec<V0_tt> &found_v0s) { return Map(found_v0s, [](const V0_tt &v0) { return v0.lv.M(); }); },
                    {"Found_" + v0_name});

    if (Settings::IsMC) {
        df_V0s = df_V0s  //
                     .Define(v0_name + "_Daughters_PdgCode", ExtractVector<Int_t>, {"Track_PdgCode", v0_name + "_TrackEntries_"})
                     .Define(v0_name + "_Daughters_IsSignal", ExtractVector<Int_t>, {"Track_IsSignal", v0_name + "_TrackEntries_"})
                     .Define(v0_name + "_McEntry_", GetTrueV0, {"MC_PdgCode", v0_name + "_TrackEntries_", "Track_McEntry_", "MC_Mother_McEntry_"})
                     .Define(v0_name + "_IsTrue", v0_name + "_McEntry_ != -1")
                     .Define(v0_name + "_IsSignal", CallExtractInt_IfNotZero, {"MC_IsSignal", v0_name + "_McEntry_", v0_name + "_IsTrue"})
                     .Define(v0_name + "_IsSecondary", CallExtractInt_IfNotZero, {"MC_IsSecondary", v0_name + "_McEntry_", v0_name + "_IsTrue"})
                     .Define(v0_name + "_IsHybrid", IsHybrid, {v0_name + "_IsSignal", "MC_IsSignal", v0_name + "_TrackEntries_", "Track_McEntry_"})
                     .Define(v0_name + "_PdgCode", CallExtractInt_IfNotZero, {"MC_PdgCode", v0_name + "_McEntry_", v0_name + "_IsTrue"})
                     .Define(v0_name + "_ReactionID", CallExtractUInt_IfNotZero, {"MC_ReactionID", v0_name + "_McEntry_", v0_name + "_IsTrue"});
    }

    return df_V0s;
}

/*
 *
 */
// void Manager::CollectTrueInfo_V0(Int_t pdg_neg, Int_t pdg_pos, Int_t pdg_hypothesis, Long64_t neg_entry, Long64_t pos_entry,  //
//  Bool_t& is_true, Int_t& mc_idx_v0, Int_t& mc_pdg_code, Bool_t& is_secondary, Bool_t& is_signal, Int_t& reaction_id,
//  Bool_t& is_hybrid) {
// /* Default values */
// is_true = kFALSE;
// mc_idx_v0 = -1;
// mc_pdg_code = 0;
// is_secondary = kFALSE;
// is_signal = kFALSE;
// reaction_id = -1;
// is_hybrid = kFALSE;
// /* Get MC particles linked to tracks */
// Long64_t Neg_McEntry = GetMcEntryFromTrackEntry(neg_entry);
// Long64_t Pos_McEntry = GetMcEntryFromTrackEntry(pos_entry);
// MC_tt MC_Neg, MC_Pos;
// if (!CopyMCParticle(Neg_McEntry, MC_Neg)) return;
// if (!CopyMCParticle(Pos_McEntry, MC_Pos)) return;
// /* Pion Pairs : motherless */
// if (pdg_hypothesis == 422) {
// is_true = !MC_Neg.HasMother() && !MC_Pos.HasMother() && MC_Neg.PdgCode == pdg_neg && MC_Pos.PdgCode == pdg_pos;
// if (is_true) {
// is_secondary = MC_Neg.IsSecondary() && MC_Pos.IsSecondary();
// is_signal = MC_Neg.IsFirstGenSignal() && MC_Pos.IsFirstGenSignal() && MC_Neg.ReactionID == MC_Pos.ReactionID;
// if (is_signal)
// reaction_id = MC_Neg.ReactionID;
// else
// is_hybrid = (MC_Neg.IsFirstGenSignal() && !MC_Pos.IsFirstGenSignal()) || (!MC_Neg.IsFirstGenSignal() && MC_Pos.IsFirstGenSignal());
// }
// return;
// }
// /* Get mother of MC particles */
// Bool_t SameMother = MC_Neg.HasMother() && MC_Pos.HasMother() && MC_Neg.Idx_Mother == MC_Pos.Idx_Mother;
// if (!SameMother) return;
// Long64_t V0_McEntry = GetMcEntry(MC_Neg.Idx_Mother);
// MC_tt MC_V0;
// if (!CopyMCParticle(V0_McEntry, MC_V0)) return;
// /*  */
// is_true = MC_Neg.PdgCode == pdg_neg && MC_Pos.PdgCode == pdg_pos && MC_V0.PdgCode == pdg_hypothesis;
// if (is_true) {
// mc_idx_v0 = MC_V0.Idx;
// mc_pdg_code = MC_V0.PdgCode;
// is_secondary = MC_V0.IsSecondary();
// is_signal = MC_V0.IsFirstGenSignal();
// if (is_signal)
// reaction_id = MC_V0.ReactionID;
// else
// is_hybrid = (MC_Neg.IsSignal() && !MC_Pos.IsSignal()) || (!MC_Neg.IsSignal() && MC_Pos.IsSignal());
// }
// return;
// }

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/*
 * Using Kalman Filter, reconstruct (anti-)sexaquark candidates via the reaction channel `AntiSexaquark Neutron -> AntiLambda K0S`
 */
void Manager::KalmanSexaquarkFinder_TypeA(Bool_t anti_channel) {
    //
    // Long64_t GenLambda_Neg_TrackEntry, GenLambda_Pos_TrackEntry;
    // Long64_t GenLambda_Neg_McEntry, GenLambda_Pos_McEntry;
    // Long64_t KaonZeroShort_Neg_TrackEntry, KaonZeroShort_Pos_TrackEntry;
    // Long64_t KaonZeroShort_Neg_McEntry, KaonZeroShort_Pos_McEntry;
    // /* Declare vectors */
    // ROOT::Math::PxPyPzEVector lvGenLambda;
    // ROOT::Math::PxPyPzEVector lvGenLambda_Neg, lvGenLambda_Pos;
    // ROOT::Math::PxPyPzEVector lvKaonZeroShort;
    // ROOT::Math::PxPyPzEVector lvKaonZeroShort_Neg, lvKaonZeroShort_Pos;
    // ROOT::Math::PxPyPzMVector lvStruckNucleon(0., 0., 0., GetMass(2112));
    // ROOT::Math::PxPyPzEVector lvSexaquark;
    // /* Candidate object */
    // UInt_t Counter = 0;
    // Candidate::ChannelA newSexaquark;
    // newSexaquark.SetPrimaryVertex(kfPrimaryVertex);
    // /* Declare KFParticle objects */
    // KFParticle kfKaonZeroShort;
    // KFParticle kfGenLambda;
    // /* Information from MC */
    // Bool_t IsSignal;
    // Int_t ReactionID;
    // Bool_t IsHybrid;
    // Int_t NonCombBkg_PdgCode;
    // /* Choose between `AntiSexaquark Neutron -> AntiLambda K0S` or `Sexaquark AntiNeutron -> Lambda K0S` */
    // std::vector<Candidate::V0>* GenericLambdas = anti_channel ? &Lambdas : &AntiLambdas;
    // TreeName Sexaquarks_TreeName = anti_channel ? TreeName::Sexaquarks_ALK0 : TreeName::Sexaquarks_LK0;
    // /* Loop over (Anti-)Lambda,KaonZeroShort pairs */
    // for (Candidate::V0 Lambda : *GenericLambdas) {
    // for (Candidate::V0 KaonZeroShort : KaonsZeroShort) {
    // /* Get track entries */
    // GenLambda_Neg_TrackEntry = GetTrackEntry(Lambda.EsdIdxNeg);
    // GenLambda_Pos_TrackEntry = GetTrackEntry(Lambda.EsdIdxPos);
    // KaonZeroShort_Neg_TrackEntry = GetTrackEntry(KaonZeroShort.EsdIdxNeg);
    // KaonZeroShort_Pos_TrackEntry = GetTrackEntry(KaonZeroShort.EsdIdxPos);
    // /** Sanity check: prevent any track from being repeated **/
    // std::set<Long64_t> unique_indices = {GenLambda_Neg_TrackEntry, GenLambda_Pos_TrackEntry, KaonZeroShort_Neg_TrackEntry,
    //  KaonZeroShort_Pos_TrackEntry};
    // if ((Int_t)unique_indices.size() < 4) continue;
    // /** Kalman Filter **/
    // kfGenLambda = Lambda.GetKf();
    // kfKaonZeroShort = KaonZeroShort.GetKf();
    // KFParticle kfSexaquark;
    // kfSexaquark.AddDaughter(kfGenLambda);
    // kfSexaquark.AddDaughter(kfKaonZeroShort);
    // /** Transport V0s to the secondary vertex **/
    // kfSexaquark.TransportToDecayVertex();
    // /** -- PENDING: does this even work? **/
    // kfGenLambda.SetProductionVertex(kfSexaquark);
    // kfGenLambda.TransportToProductionVertex();
    // kfKaonZeroShort.SetProductionVertex(kfSexaquark);
    // kfKaonZeroShort.TransportToProductionVertex();
    // /** Get and transport tracks **/
    // /** -- PENDING: study this further... **/
    // // kfGenLambda_Neg = Math::TransportKFParticle(kfGenLambda_Neg, kfGenLambda_Pos, GetMass(2212),  //
    // // (Int_t)Track_GenLambda_Neg.Charge);
    // // kfGenLambda_Pos = Math::TransportKFParticle(kfGenLambda_Pos, kfGenLambda_Neg, GetMass(211),  //
    // // (Int_t)Track_GenLambda_Pos.Charge);
    // // kfKaonZeroShort_Neg = Math::TransportKFParticle(kfKaonZeroShort_Neg, kfKaonZeroShort_Pos, GetMass(211),  //
    // // (Int_t)Track_KaonZeroShort_Neg.Charge);
    // // kfKaonZeroShort_Pos = Math::TransportKFParticle(kfKaonZeroShort_Pos, kfKaonZeroShort_Neg, GetMass(211),  //
    // // (Int_t)Track_KaonZeroShort_Pos.Charge);
    // // lvGenLambda_Neg =  //
    // // ROOT::Math::PxPyPzMVector(kfGenLambda_Neg.Px(), kfGenLambda_Neg.Py(), kfGenLambda_Neg.Pz(), GetMass(2212));
    // // lvGenLambda_Pos =  //
    // // ROOT::Math::PxPyPzMVector(kfGenLambda_Pos.Px(), kfGenLambda_Pos.Py(), kfGenLambda_Pos.Pz(), GetMass(211));
    // // lvKaonZeroShort_Neg =
    // // ROOT::Math::PxPyPzMVector(kfKaonZeroShort_Neg.Px(), kfKaonZeroShort_Neg.Py(), kfKaonZeroShort_Neg.Pz(), GetMass(211));
    // // lvKaonZeroShort_Pos =
    // // ROOT::Math::PxPyPzMVector(kfKaonZeroShort_Pos.Px(), kfKaonZeroShort_Pos.Py(), kfKaonZeroShort_Pos.Pz(), GetMass(211));
    // // lvGenLambda = lvGenLambda_Neg + lvGenLambda_Pos;
    // // lvKaonZeroShort = lvKaonZeroShort_Neg + lvKaonZeroShort_Pos;
    // /** Reconstruct (anti-)sexaquark candidate **/
    // lvGenLambda = ROOT::Math::PxPyPzEVector(kfGenLambda.Px(), kfGenLambda.Py(), kfGenLambda.Pz(), kfGenLambda.E());
    // lvKaonZeroShort = ROOT::Math::PxPyPzEVector(kfKaonZeroShort.Px(), kfKaonZeroShort.Py(), kfKaonZeroShort.Pz(), kfKaonZeroShort.E());
    // lvSexaquark = lvGenLambda + lvKaonZeroShort - lvStruckNucleon;
    // /** Prepare candidate object **/
    // // newSexaquark.SetSexaquarkInfo(); // PENDING
    // newSexaquark.SetKinematics(lvSexaquark, lvGenLambda, lvKaonZeroShort);
    // newSexaquark.SetGeometry(kfSexaquark, kfGenLambda, kfKaonZeroShort,  //
    //  Lambda.GetKfNeg(), Lambda.GetKfPos(), KaonZeroShort.GetKfNeg(), KaonZeroShort.GetKfPos());
    // /** Collect true information **/
    // /*
    // if (Settings::IsMC) {
    // CollectTrueInfo_ChannelA(); // PENDING
    // newSexaquark.SetTrueInfo(IsSignal, ReactionID, IsHybrid, NonCombBkg_PdgCode);
    // }
    // */
    // /** Apply cuts **/
    // if (!Inspector.Approve(newSexaquark)) continue;
    // /** Store **/
    // FillSexaquark(Sexaquarks_TreeName, newSexaquark);
    // Counter++;
    // }  // end of loop over K0S
    // }      // end of loop over (anti-)lambdas
    // InfoF("N Found Channel A Candidates: %u", Counter);
}

void Manager::KalmanSexaquarkFinder_TypeDE(Bool_t anti_channel) {
    //
}

/*
 *
 */
void Manager::KalmanSexaquarkFinder_TypeH(Bool_t anti_channel) {
    //
}

void Manager::CollectTrueInfo_ChannelA() {
    //
    // mcIdxNeg_GenLambda = getMcIdx_fromEsdIdx[esdIdxNeg_GenLambda];
    // mcIdxPos_GenLambda = getMcIdx_fromEsdIdx[esdIdxPos_GenLambda];
    // mcIdxNeg_KaonZeroShort = getMcIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort];
    // mcIdxPos_KaonZeroShort = getMcIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];
    // is_signal = isMcIdxSignal[mcIdxNeg_GenLambda] && isMcIdxSignal[mcIdxPos_GenLambda] && isMcIdxSignal[mcIdxNeg_KaonZeroShort] &&
    // isMcIdxSignal[mcIdxPos_KaonZeroShort] &&  //
    // getPdgCode_fromMcIdx[mcIdxNeg_GenLambda] == getNegPdgCode_fromV0PdgCode[pdgCodeV0] &&
    // getPdgCode_fromMcIdx[mcIdxPos_GenLambda] == getPosPdgCode_fromV0PdgCode[pdgCodeV0] &&
    // getPdgCode_fromMcIdx[mcIdxNeg_KaonZeroShort] == getNegPdgCode_fromV0PdgCode[pdgCodeKaon] &&
    // getPdgCode_fromMcIdx[mcIdxPos_KaonZeroShort] == getPosPdgCode_fromV0PdgCode[pdgCodeKaon] &&  //
    // getReactionID_fromMcIdx[mcIdxNeg_GenLambda] == getReactionID_fromMcIdx[mcIdxPos_GenLambda] &&
    // getReactionID_fromMcIdx[mcIdxPos_GenLambda] == getReactionID_fromMcIdx[mcIdxNeg_KaonZeroShort] &&
    // getReactionID_fromMcIdx[mcIdxNeg_KaonZeroShort] == getReactionID_fromMcIdx[mcIdxPos_KaonZeroShort];
    // reaction_id = getReactionID_fromMcIdx[mcIdxNeg_GenLambda] == getReactionID_fromMcIdx[mcIdxPos_GenLambda] &&
    //   getReactionID_fromMcIdx[mcIdxPos_GenLambda] == getReactionID_fromMcIdx[mcIdxNeg_KaonZeroShort] &&
    //   getReactionID_fromMcIdx[mcIdxNeg_KaonZeroShort] == getReactionID_fromMcIdx[mcIdxPos_KaonZeroShort]
    //   ? getReactionID_fromMcIdx[mcIdxNeg_GenLambda]
    //   : -1;
    // is_hybrid = (isMcIdxSignal[mcIdxNeg_GenLambda] || isMcIdxSignal[mcIdxPos_GenLambda] || isMcIdxSignal[mcIdxNeg_KaonZeroShort] ||
    //  isMcIdxSignal[mcIdxPos_KaonZeroShort]) &&
    // !is_signal;
    // same_ancestor = getAncestorMcIdx_fromMcIdx[mcIdxNeg_GenLambda] != -1 &&
    // getAncestorMcIdx_fromMcIdx[mcIdxNeg_GenLambda] == getAncestorMcIdx_fromMcIdx[mcIdxPos_GenLambda] &&
    // getAncestorMcIdx_fromMcIdx[mcIdxPos_GenLambda] == getAncestorMcIdx_fromMcIdx[mcIdxNeg_KaonZeroShort] &&
    // getAncestorMcIdx_fromMcIdx[mcIdxNeg_KaonZeroShort] == getAncestorMcIdx_fromMcIdx[mcIdxPos_KaonZeroShort];
    // is_noncomb_bkg = !is_signal && !is_hybrid && same_ancestor;
    // ancestor_idx = is_noncomb_bkg ? getAncestorMcIdx_fromMcIdx[mcIdxNeg_GenLambda] : -1;
}

void Manager::CollectTrueInfo_ChannelD() {
    //
}

void Manager::CollectTrueInfo_ChannelE() {
    //
}

void Manager::CollectTrueInfo_ChannelH() {
    //
}

/*           */
/* Utilities */
/* ========= */

void Manager::PrintAll(RNode df) {
    /* # MC */
    auto n_mc = df.Take<size_t>("N_MC");
    auto ev_mc_pdgcode = df.Take<RVecI>("MC_PdgCode");
    auto ev_mc_is_charged = df.Take<RVecI>("MC_IsProtonKaonPion");
    auto ev_mc_is_signal = df.Take<RVecI>("MC_IsSignal");
    auto ev_mc_xv = df.Take<RVecF>("MC_Xv");
    /* -- Mother */
    auto ev_mc_mother_entry = df.Take<RVecL>("MC_Mother_McEntry_");
    auto ev_mc_mother_pdgcode = df.Take<RVecI>("MC_Mother_PdgCode");
    /* -- Daughter */
    auto ev_mc_neg_entry = df.Take<RVecL>("MC_NegDau_McEntry_");
    auto ev_mc_pos_entry = df.Take<RVecL>("MC_PosDau_McEntry_");
    auto ev_mc_neg_pdgcode = df.Take<RVecI>("MC_NegDau_PdgCode");
    auto ev_mc_pos_pdgcode = df.Take<RVecI>("MC_PosDau_PdgCode");
    /* -- Ancestor */
    auto ev_mc_ancestor_entry = df.Take<RVecL>("MC_Ancestor_McEntry_");
    auto ev_mc_reaction_id = df.Take<RVecU>("MC_ReactionID");
    /* -- Linked Info */
    auto ev_mc_got_tracked = df.Take<RVecI>("MC_GotTracked");
    auto ev_mc_neg_gottracked = df.Take<RVecI>("MC_NegDau_GotTracked");
    auto ev_mc_pos_gottracked = df.Take<RVecI>("MC_PosDau_GotTracked");
    auto ev_mc_truev0_isfindable = df.Take<RVecI>("MC_TrueV0_IsFindable");
    /* # Injected */
    auto n_injected = df.Take<ULong_t>("N_Injected");
    auto ev_injected_reaction_id = df.Take<RVecU>("ReactionID");
    auto ev_injected_mc_entries = df.Take<RVec<RVecUL>>("Injected_Products_McEntries_");
    auto ev_injected_mc_pdgcode = df.Take<RVec<RVecI>>("Injected_Products_PdgCode");
    auto ev_injected_radius = df.Take<RVecF>("Injected_Radius");
    auto ev_injected_pt_after_reaction = df.Take<RVecF>("Injected_Pt_AfterReaction");
    auto ev_injected_is_findable = df.Take<RVecI>("Injected_IsFindable");
    /* # Tracks*/
    auto n_tracks = df.Take<ULong_t>("N_Tracks");
    auto ev_track_mc_entry = df.Take<RVecUL>("Track_McEntry");
    auto ev_track_nsigmaproton = df.Take<RVecF>("Track_NSigmaProton");
    auto ev_track_nsigmapion = df.Take<RVecF>("Track_NSigmaPion");
    auto ev_track_pid_as_antiproton = df.Take<RVecI>("Track_PID_AsAntiProton");
    auto ev_track_pid_as_piminus = df.Take<RVecI>("Track_PID_AsPiMinus");
    auto ev_track_pid_as_piplus = df.Take<RVecI>("Track_PID_AsPiPlus");
    auto ev_track_pdgcode = df.Take<RVecI>("Track_PdgCode");
    auto ev_track_issignal = df.Take<RVecI>("Track_IsSignal");
    /* # Found V0s */
    /* -- AntiLambdas */
    auto ev_found_antilambdas = df.Take<RVec<RVecUL>>("AntiLambdas_TrackEntries_");
    auto ev_found_antilambdas_pdgcodes = df.Take<RVec<RVecI>>("AntiLambdas_Daughters_PdgCode");
    auto ev_found_antilambdas_mc_entry = df.Take<RVecL>("AntiLambdas_McEntry_");
    auto ev_found_antilambdas_istrue = df.Take<RVecI>("AntiLambdas_IsTrue");
    auto ev_found_antilambdas_issignal = df.Take<RVecI>("AntiLambdas_IsSignal");
    auto ev_found_antilambdas_issecondary = df.Take<RVecI>("AntiLambdas_IsSecondary");
    auto ev_found_antilambdas_ishybrid = df.Take<RVecI>("AntiLambdas_IsHybrid");
    auto ev_found_antilambdas_mc_pdgcode = df.Take<RVecI>("AntiLambdas_PdgCode");
    auto ev_found_antilambdas_reaction_id = df.Take<RVecU>("AntiLambdas_ReactionID");
    /* -- Lambdas */
    /*
    auto ev_found_lambdas = df.Take<RVec<RVecUL>>("Lambdas_TrackEntries_");
    auto ev_found_lambdas_pdgcodes = df.Take<RVec<RVecI>>("Lambdas_Daughters_PdgCode");
    auto ev_found_lambdas_mc_entry = df.Take<RVecL>("Lambdas_McEntry_");
    auto ev_found_lambdas_istrue = df.Take<RVecI>("Lambdas_IsTrue");
    auto ev_found_lambdas_issignal = df.Take<RVecI>("Lambdas_IsSignal");
    auto ev_found_lambdas_issecondary = df.Take<RVecI>("Lambdas_IsSecondary");
    auto ev_found_lambdas_ishybrid = df.Take<RVecI>("Lambdas_IsHybrid");
    auto ev_found_lambdas_mc_pdgcode = df.Take<RVecI>("Lambdas_PdgCode");
    auto ev_found_lambdas_reaction_id = df.Take<RVecU>("Lambdas_ReactionID");
    */
    /* -- KaonsZeroShort */
    /*
        auto ev_found_k0s = df.Take<RVec<RVecUL>>("KaonsZeroShort_TrackEntries_");
        auto ev_found_k0s_pdgcodes = df.Take<RVec<RVecI>>("KaonsZeroShort_Daughters_PdgCode");
        auto ev_found_k0s_mc_entry = df.Take<RVecL>("KaonsZeroShort_McEntry_");
        auto ev_found_k0s_istrue = df.Take<RVecI>("KaonsZeroShort_IsTrue");
        auto ev_found_k0s_issignal = df.Take<RVecI>("KaonsZeroShort_IsSignal");
        auto ev_found_k0s_issecondary = df.Take<RVecI>("KaonsZeroShort_IsSecondary");
        auto ev_found_k0s_ishybrid = df.Take<RVecI>("KaonsZeroShort_IsHybrid");
        auto ev_found_k0s_mc_pdgcode = df.Take<RVecI>("KaonsZeroShort_PdgCode");
        auto ev_found_k0s_reaction_id = df.Take<RVecU>("KaonsZeroShort_ReactionID");
     */
    /* # Sexaquarks */
    /* -- Channel A */
    /*
    auto ev_found_channel_a = df.Take<RVec<RVecUL>>("ChannelA_V0Entries_");
    auto ev_found_channel_a_issignal = df.Take<RVecI>("ChannelA_IsSignal");
    */
    // --- //
    size_t n_events = n_tracks->size();
    for (size_t event_i = 0; event_i < n_events; event_i++) {
        std::cout << "# Event " << event_i << std::endl;
        // --- //
        std::cout << "## Summary" << std::endl;
        // printf(">> N Injected = %lu\n", (*n_injected)[event_i]);
        // printf(">> N Findable Sexaquarks = %d\n", Sum((*ev_injected_is_findable)[event_i]));
        // printf(">> N MC = %lu\n", n_mc->at(event_i));
        // printf(">> N Findable V0s = %d\n", Sum((*ev_mc_truev0_isfindable)[event_i]));
        printf(">> N Found AntiLambdas = %lu\n", (*ev_found_antilambdas)[event_i].size());
        printf(">> N Found AntiLambdas (True) = %d\n", Sum((*ev_found_antilambdas_istrue)[event_i]));
        printf(">> N Found AntiLambdas (Signal) = %d\n", Sum((*ev_found_antilambdas_issignal)[event_i]));
        printf(">> N Found AntiLambdas (Secondary) = %d\n", Sum((*ev_found_antilambdas_issecondary)[event_i]));
        printf(">> N Found AntiLambdas (Hybrid) = %d\n", Sum((*ev_found_antilambdas_ishybrid)[event_i]));
        // printf(">> N Found Lambdas = %lu\n", ev_found_lambdas[event_i].size());
        // printf(">> N Found Lambdas (True) = %d\n", Sum(ev_found_lambdas_istrue[event_i]));
        // printf(">> N Found Lambdas (Signal) = %d\n", Sum(ev_found_lambdas_issignal[event_i]));
        // printf(">> N Found Lambdas (Secondary) = %d\n", Sum(ev_found_lambdas_issecondary[event_i]));
        // printf(">> N Found KaonsZeroShort = %lu\n", (*ev_found_k0s)[event_i].size());
        // printf(">> N Found KaonsZeroShort (True) = %d\n", Sum((*ev_found_k0s_istrue)[event_i]));
        // printf(">> N Found KaonsZeroShort (Signal) = %d\n", Sum((*ev_found_k0s_issignal)[event_i]));
        // printf(">> N Found KaonsZeroShort (Secondary) = %d\n", Sum((*ev_found_k0s_issecondary)[event_i]));
        /*
        printf(">> N Found AntiSexaquarks (Channel A) = %lu\n", ev_found_channel_a[event_i].size());
        printf(">> N Found AntiSexaquarks (Channel A, Signal) = %d\n", Sum(ev_found_channel_a_issignal[event_i]));
        */
        // --- //
        std::cout << "## (MC) Particles" << std::endl;
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            printf("mc_entry=%lu, pdg_code=%d, ancestor_entry=%ld, mother_entry=%ld, neg_entry=%ld, pos_entry=%ld\n",  //
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_ancestor_entry)[event_i][mc_entry],
                   (*ev_mc_mother_entry)[event_i][mc_entry], (*ev_mc_neg_entry)[event_i][mc_entry], (*ev_mc_pos_entry)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## (MC) V0s" << std::endl;
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            if (!((*ev_mc_pdgcode)[event_i][mc_entry] == 310 || TMath::Abs((*ev_mc_pdgcode)[event_i][mc_entry]) == 3122)) continue;
            printf("mc_entry=%lu, pdg_code=%d, is_signal=%d, neg_entry=%ld, pos_entry=%ld, neg_gt=%d, pos_gt=%d, neg_pdg=%d, pos_pdg=%d, xv=%f\n",
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_is_signal)[event_i][mc_entry], (*ev_mc_neg_entry)[event_i][mc_entry],
                   (*ev_mc_pos_entry)[event_i][mc_entry], (*ev_mc_neg_gottracked)[event_i][mc_entry], (*ev_mc_pos_gottracked)[event_i][mc_entry],
                   (*ev_mc_neg_pdgcode)[event_i][mc_entry], (*ev_mc_pos_pdgcode)[event_i][mc_entry], (*ev_mc_xv)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## (MC) Charged Particles" << std::endl;
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            if (!(*ev_mc_is_charged)[event_i][mc_entry]) continue;
            printf("mc_entry=%ld, pdg_code=%d, is_signal=%d, mother_entry=%ld, mother_pdg=%d, got_tracked=%d\n",  //
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_is_signal)[event_i][mc_entry], (*ev_mc_mother_entry)[event_i][mc_entry],
                   (*ev_mc_mother_pdgcode)[event_i][mc_entry], (*ev_mc_got_tracked)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## Injected" << std::endl;
        for (size_t inj_entry = 0; inj_entry < (*n_injected)[event_i]; inj_entry++) {
            printf("reaction_id=%d, radius=%f, pt_after_reaction=%f, is_findable=%d\n",  //
                   (*ev_injected_reaction_id)[event_i][inj_entry], (*ev_injected_radius)[event_i][inj_entry],
                   (*ev_injected_pt_after_reaction)[event_i][inj_entry], (*ev_injected_is_findable)[event_i][inj_entry]);
            for (size_t prod_entry = 0; prod_entry < (*ev_injected_mc_entries)[event_i][inj_entry].size(); prod_entry++) {
                printf(">> prod_entry=%lu, pdg_code=%d\n",  //
                       (*ev_injected_mc_entries)[event_i][inj_entry][prod_entry], (*ev_injected_mc_pdgcode)[event_i][inj_entry][prod_entry]);
            }
        }
        // --- //
        std::cout << "## Tracks" << std::endl;
        for (size_t track_entry = 0; track_entry < (*n_tracks)[event_i]; track_entry++) {
            printf("track_entry=%ld, mc_entry=%lu, nsigmaproton=%f, nsigmapion=%f, pdg_code=%d, is_signal=%d\n",  //
                   track_entry, (*ev_track_mc_entry)[event_i][track_entry], (*ev_track_nsigmaproton)[event_i][track_entry],
                   (*ev_track_nsigmapion)[event_i][track_entry], (*ev_track_pdgcode)[event_i][track_entry],
                   (*ev_track_issignal)[event_i][track_entry]);
        }
        // --- //
        std::cout << "## V0s" << std::endl;
        // --- //
        std::cout << "## >> (Found) AntiLambdas" << std::endl;
        for (size_t v0_entry = 0; v0_entry < (*ev_found_antilambdas)[event_i].size(); v0_entry++) {
            printf(
                "v0_entry=%lu, neg=%lu, pos=%lu, neg_pdg=%d, pos_pdg=%d, is_true=%d, is_signal=%d, is_secondary=%d, is_hybrid=%d, mc_entry=%ld, "
                "mc_pdg_code=%d, reaction_id=%d\n",
                v0_entry, (*ev_found_antilambdas)[event_i][v0_entry][0], (*ev_found_antilambdas)[event_i][v0_entry][1],
                (*ev_found_antilambdas_pdgcodes)[event_i][v0_entry][0], (*ev_found_antilambdas_pdgcodes)[event_i][v0_entry][1],
                (*ev_found_antilambdas_istrue)[event_i][v0_entry], (*ev_found_antilambdas_issignal)[event_i][v0_entry],
                (*ev_found_antilambdas_issecondary)[event_i][v0_entry], (*ev_found_antilambdas_ishybrid)[event_i][v0_entry],
                (*ev_found_antilambdas_mc_entry)[event_i][v0_entry], (*ev_found_antilambdas_mc_pdgcode)[event_i][v0_entry],
                (*ev_found_antilambdas_reaction_id)[event_i][v0_entry]);
        }
        /*
        std::cout << "## >> (Found) Lambdas" << std::endl;
        for (size_t v0_entry = 0; v0_entry < (*ev_found_lambdas)[event_i].size(); v0_entry++) {
            printf(
                "v0_entry=%lu, neg=%lu, pos=%lu, neg_pdg=%d, pos_pdg=%d, is_true=%d, is_signal=%d, is_secondary=%d, is_hybrid=%d, mc_entry=%ld, "
                "mc_pdg_code=%d, reaction_id=%d\n",
                v0_entry, (*ev_found_lambdas)[event_i][v0_entry][0], (*ev_found_lambdas)[event_i][v0_entry][1],
                (*ev_found_lambdas_pdgcodes)[event_i][v0_entry][0], (*ev_found_lambdas_pdgcodes)[event_i][v0_entry][1],
                (*ev_found_lambdas_istrue)[event_i][v0_entry], (*ev_found_lambdas_issignal)[event_i][v0_entry],
                (*ev_found_lambdas_issecondary)[event_i][v0_entry], (*ev_found_lambdas_ishybrid)[event_i][v0_entry],
                (*ev_found_lambdas_mc_entry)[event_i][v0_entry], (*ev_found_lambdas_mc_pdgcode)[event_i][v0_entry],
                (*ev_found_lambdas_reaction_id)[event_i][v0_entry]);
            }
        */
        /*
                std::cout << "## >> (Found) KaonsZeroShort" << std::endl;
                for (size_t v0_entry = 0; v0_entry < (*ev_found_k0s)[event_i].size(); v0_entry++) {
                    printf(
                        "v0_entry=%lu, neg=%lu, pos=%lu, neg_pdg=%d, pos_pdg=%d, is_true=%d, is_signal=%d, is_secondary=%d, is_hybrid=%d,
           mc_entry=%ld, " "mc_pdg_code=%d, reaction_id=%d\n", v0_entry, (*ev_found_k0s)[event_i][v0_entry][0], (*ev_found_k0s)[event_i][v0_entry][1],
                        (*ev_found_k0s_pdgcodes)[event_i][v0_entry][0], (*ev_found_k0s_pdgcodes)[event_i][v0_entry][1],
                        (*ev_found_k0s_istrue)[event_i][v0_entry], (*ev_found_k0s_issignal)[event_i][v0_entry],
                        (*ev_found_k0s_issecondary)[event_i][v0_entry], (*ev_found_k0s_ishybrid)[event_i][v0_entry],
                        (*ev_found_k0s_mc_entry)[event_i][v0_entry], (*ev_found_k0s_mc_pdgcode)[event_i][v0_entry],
                        (*ev_found_k0s_reaction_id)[event_i][v0_entry]);
                }
        */
        // --- //
        /*
        std::cout << "## Sexaquarks" << std::endl;
        */
        // --- //
        /*
        std::cout << "## (Found) AntiSexaquarks (Channel A)" << std::endl;
        for (size_t sexa_entry = 0; sexa_entry < ev_found_channel_a[event_i].size(); sexa_entry++) {
            printf("v0a_entry=%lu, v0b_entry=%lu, is_signal=%d\n", ev_found_channel_a[event_i][sexa_entry][0],
                    ev_found_channel_a[event_i][sexa_entry][1], ev_found_channel_a_issignal[event_i][sexa_entry]);
        }
        */
    }  // end of loop over events
    /* Print NRuns */
    InfoF("NRuns = %u", df.GetNRuns());
}

void Manager::EndOfAnalysis(RNode df) {
    df.Snapshot("V0s", Settings::PathOutputFile, {"AntiLambdas_M"});
    InfoF("TFile %s has been written", Settings::PathOutputFile.c_str());
    /* Print NRuns */
    InfoF("NRuns = %u", df.GetNRuns());
}

}  // namespace Analysis
}  // namespace Tree2Sexaquark
