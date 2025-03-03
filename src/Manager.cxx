#include "Analysis/Manager.hxx"

#include <cstddef>
#include <set>

#include "KFVertex.h"
#include "RtypesCore.h"
#include "TDatabasePDG.h"
#include "TROOT.h"

#include "Analysis/Settings.hxx"
#include "Cuts/Default.hxx"
#include "Math/Common.hxx"
#include "Math/KalmanFilter.hxx"

namespace Tree2Sexaquark::Analysis {

/*                   */
/* Vector Gymnastics */
/* ================= */

// mask a vector, e.g. {0,2} -> {1,0,1}
RVecI Manager::Mask(cRVecUL entries, size_t reference_size) {
    RVecI output(reference_size);
    for (auto entry : entries) {
        output[entry] = 1;
    }
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
    if (Settings::NThreads > 1) {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT(Settings::NThreads);
    }
}

RNode Manager::ProcessEvent(RNode df) {
    //
    auto CallCreateKFVertex = [](Float_t x, Float_t y, Float_t z, cRVecF cov_matrix) -> KFVertex {
        Float_t XYZ[3] = {x, y, z};
        Float_t CovMatrix[6] = {cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3], cov_matrix[4], cov_matrix[5]};
        return Math::CreateKFVertex(XYZ, CovMatrix);
    };

    RNode DF_Events = df  //
                          .Define("Event_PV", [](Float_t x, Float_t y, Float_t z) { return XYZPoint(x, y, z); }, {"PV_Xv", "PV_Yv", "PV_Zv"})
                          .Define("Event_KF_PV", CallCreateKFVertex, {"PV_Xv", "PV_Yv", "PV_Zv", "PV_CovMatrix"});

    return DF_Events;
}

RNode Manager::ProcessMCParticles(RNode df) {
    //
    auto GetCharge = [](cRVecI pdg_codes) -> RVecI {
        RVecI mc_charge;
        mc_charge.reserve(pdg_codes.size());
        for (auto pdg_code : pdg_codes) {
            auto *particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
            if (!particle) {
                mc_charge.emplace_back(0);
            } else {
                mc_charge.emplace_back(Int_t(particle->Charge() / 3.));  // convert units from |e|/3 to |e|
            }
        }
        return mc_charge;
    };
    //
    auto GetDaughter = [](cRVecI mc_mask, cRVecL mother_mc_entry_) -> RVecL {
        RVecL output(mc_mask.size(), -1);
        auto daughters_entries = Nonzero(mc_mask);
        auto mothers_entries = mother_mc_entry_[mc_mask];
        for (size_t i = 0; i < mothers_entries.size(); i++) {
            output[mothers_entries[i]] = (Long_t)daughters_entries[i];
        }
        return output;
    };
    //
    auto GetAncestor = [](cRVecL mother_mc_entry_) -> RVecL {
        RVecL output(mother_mc_entry_.size(), -1);
        Long_t curr_entry = 0;
        for (size_t i = 0; i < mother_mc_entry_.size(); i++) {
            // if no mother, the particle is its own ancestor
            if (mother_mc_entry_[i] < 0) {
                output[i] = (Long_t)i;
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
        for (const auto &products_vec : injected_mc_entries_) {
            all_findable = true;
            for (auto mc_entry : products_vec) {
                if (TMath::Abs(pdg_code[mc_entry]) == PdgCode::Lambda || pdg_code[mc_entry] == PdgCode::KaonZeroShort) {
                    mc_is_findable = mc_true_v0_is_findable[mc_entry];
                } else {
                    mc_is_findable = mc_got_tracked[mc_entry];
                }
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
    const Double_t mass_proton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    const Double_t mass_kaon = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    const Double_t mass_pion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    //
    auto GetTrackCovMatrix = [](cRVecF track_sigma_y2, cRVecF track_sigma_zy, cRVecF track_sigma_z2,             //
                                cRVecF track_sigma_snp_y, cRVecF track_sigma_snp_z, cRVecF track_sigma_snp2,     //
                                cRVecF track_sigma_tgl_y, cRVecF track_sigma_tgl_z, cRVecF track_sigma_tgl_snp,  //
                                cRVecF track_sigma_tgl2, cRVecF track_sigma_1pt_y, cRVecF track_sigma_1pt_z,     //
                                cRVecF track_sigma_1pt_snp, cRVecF track_sigma_1pt_tgl, cRVecF track_sigma_1pt2) -> RVec<RVecF> {
        RVec<RVecF> output;
        output.reserve(track_sigma_y2.size());
        RVecF cov_matrix(15);
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

    RNode df_Tracks = df.Define("N_Tracks", "Track_Px.size()")
                          /* # PID */
                          .Define("Track_PID_AsAntiProton", "abs(Track_NSigmaProton) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsProton", "abs(Track_NSigmaProton) < 3. && Track_Charge > 0")
                          .Define("Track_PID_AsNegKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsPosKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge > 0")
                          .Define("Track_PID_AsPiMinus", "abs(Track_NSigmaPion) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsPiPlus", "abs(Track_NSigmaPion) < 3. && Track_Charge > 0")
                          .Define("AntiProton_TrackEntry_", "Nonzero(Track_PID_AsAntiProton)")
                          .Define("Proton_TrackEntry_", "Nonzero(Track_PID_AsProton)")
                          .Define("NegKaon_TrackEntry_", "Nonzero(Track_PID_AsNegKaon)")
                          .Define("PosKaon_TrackEntry_", "Nonzero(Track_PID_AsPosKaon)")
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
RNode Manager::FindV0s(RNode df, Int_t pdg_code_v0, Int_t pdg_code_neg, Int_t pdg_code_pos) {
    //
    const Double_t neg_mass = TDatabasePDG::Instance()->GetParticle(pdg_code_neg)->Mass();
    const Double_t pos_mass = TDatabasePDG::Instance()->GetParticle(pdg_code_pos)->Mass();
    //
    std::string v0_name, neg_name, pos_name;
    if (pdg_code_v0 == PdgCode::Lambda) {
        v0_name = "L";
        neg_name = "PiMinus";
        pos_name = "Proton";
    } else if (pdg_code_v0 == PdgCode::AntiLambda) {
        v0_name = "AL";
        neg_name = "AntiProton";
        pos_name = "PiPlus";
    } else if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        v0_name = "K0S";
        neg_name = "PiMinus";
        pos_name = "PiPlus";
    }
    fAnalyzed_V0sNames.push_back(v0_name);
    //
    auto KF_V0Finder = [pdg_code_v0, neg_mass, pos_mass](cRVecUL entries_neg_, cRVecUL entries_pos_,                            //
                                                         cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,         //
                                                         RVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,  //
                                                         const RVec<RVecF> &cov_matrix,                                         //
                                                         const Float_t &magnetic_field, const XYZPoint &v3_pv) -> RVec<FoundV0> {
        RVec<FoundV0> output;
        FoundV0 this_v0;
        /* protection */
        if (!entries_neg_.size() || !entries_pos_.size()) return output;
        /* one-time definitions */
        KFParticle::SetField(magnetic_field);
        /* main vars */
        XYZPoint v3_v0;
        PxPyPzMVector lv_neg;
        PxPyPzMVector lv_pos;
        PxPyPzMVector lv_v0;
        /* auxiliary vars */
        Double_t cpa_wrt_pv, dca_wrt_pv;
        Double_t dca_btw_dau, dca_neg_v0, dca_pos_v0;
        Double_t arm_qt, arm_alpha, qt_over_alpha;
        /* loop over all pairs */
        for (auto neg : entries_neg_) {
            for (auto pos : entries_pos_) {
                /* sanity check */
                if (neg == pos) continue;
                /* geometry */
                KFParticle kf_neg = Math::CreateKFParticle(px[neg], py[neg], pz[neg], x[neg], y[neg], z[neg],            //
                                                           charge[neg], alpha[neg], snp[neg], tgl[neg], signed1pt[neg],  //
                                                           cov_matrix[neg], neg_mass);
                KFParticle kf_pos = Math::CreateKFParticle(px[pos], py[pos], pz[pos], x[pos], y[pos], z[pos],            //
                                                           charge[pos], alpha[pos], snp[pos], tgl[pos], signed1pt[pos],  //
                                                           cov_matrix[pos], pos_mass);
                KFParticle kf_v0;
                kf_v0.AddDaughter(kf_neg);
                kf_v0.AddDaughter(kf_pos);
                kf_v0.TransportToDecayVertex();
                v3_v0.SetCoordinates(kf_v0.GetX(), kf_v0.GetY(), kf_v0.GetZ());
                /* transport */
                kf_neg.SetProductionVertex(kf_v0);     // TEST
                kf_pos.SetProductionVertex(kf_v0);     // TEST
                kf_neg.TransportToProductionVertex();  // TEST
                kf_pos.TransportToProductionVertex();  // TEST
                /* kinematics */
                lv_neg.SetCoordinates(kf_neg.Px(), kf_neg.Py(), kf_neg.Pz(), neg_mass);
                lv_pos.SetCoordinates(kf_pos.Px(), kf_pos.Py(), kf_pos.Pz(), pos_mass);
                lv_v0 = lv_neg + lv_pos;
                /* calculate vars */
                dca_btw_dau = TMath::Abs((Double_t)kf_neg.GetDistanceFromParticle(kf_pos));
                dca_neg_v0 = TMath::Abs((Double_t)kf_neg.GetDistanceFromVertex(kf_v0));
                dca_pos_v0 = TMath::Abs((Double_t)kf_pos.GetDistanceFromVertex(kf_v0));
                /* -- (anti)lambda cuts */
                if (TMath::Abs(pdg_code_v0) == PdgCode::Lambda) {
                    if (TMath::Abs(lv_v0.Eta()) > Default::Lambda::AbsMaxEta) continue;
                    if (lv_v0.Pt() < Default::Lambda::MinPt) continue;
                    if (lv_v0.M() < Default::Lambda::MinMass || lv_v0.M() > Default::Lambda::MaxMass) continue;
                    if (v3_v0.Rho() < Default::Lambda::MinRadius) continue;  // NOTE: could set upper limit
                    if ((v3_v0 - v3_pv).R() < Default::Lambda::MinDistFromPV) continue;
                    /*  */ cpa_wrt_pv = Math::CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
                    if (cpa_wrt_pv < Default::Lambda::MinCPAwrtPV || cpa_wrt_pv > Default::Lambda::MaxCPAwrtPV) continue;
                    /*  */ dca_wrt_pv = Math::LinePointDCA(lv_v0.Vect(), v3_v0, v3_pv);
                    if (dca_wrt_pv < Default::Lambda::MinDCAwrtPV) continue;
                    if (dca_btw_dau > Default::Lambda::MaxDCAbtwDau) continue;
                    if (dca_neg_v0 > Default::Lambda::MaxDCAnegV0) continue;
                    if (dca_pos_v0 > Default::Lambda::MaxDCAposV0) continue;
                    /*  */ arm_qt = Math::ArmenterosQt(lv_v0.Vect(), lv_neg.Vect());
                    /*  */ arm_alpha = Math::ArmenterosAlpha(lv_v0.Vect(), lv_neg.Vect(), lv_pos.Vect());
                    /*  */ qt_over_alpha = arm_qt / TMath::Abs(arm_alpha);
                    if (qt_over_alpha > Default::Lambda::AbsMaxArmQtOverAlpha) continue;
                }
                /* -- k0s cuts */
                if (pdg_code_v0 == PdgCode::KaonZeroShort) {
                    if (TMath::Abs(lv_v0.Eta()) > Default::KaonZeroShort::AbsMaxEta) continue;
                    if (lv_v0.Pt() < Default::KaonZeroShort::MinPt) continue;
                    if (lv_v0.M() < Default::KaonZeroShort::MinMass || lv_v0.M() > Default::KaonZeroShort::MaxMass) continue;
                    if (v3_v0.Rho() < Default::KaonZeroShort::MinRadius) continue;      // NOTE: could set upper limit
                    if ((v3_v0 - v3_pv).R() < Default::KaonZeroShort::MinDistFromPV ||  //
                        (v3_v0 - v3_pv).R() > Default::KaonZeroShort::MaxDistFromPV)
                        continue;
                    /*  */ cpa_wrt_pv = Math::CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
                    if (cpa_wrt_pv < Default::KaonZeroShort::MinCPAwrtPV || cpa_wrt_pv > Default::KaonZeroShort::MaxCPAwrtPV) continue;
                    /*  */ dca_wrt_pv = Math::LinePointDCA(lv_v0.Vect(), v3_v0, v3_pv);
                    if (dca_wrt_pv < Default::KaonZeroShort::MinDCAwrtPV) continue;
                    if (dca_btw_dau > Default::KaonZeroShort::MaxDCAbtwDau) continue;
                    if (dca_neg_v0 > Default::KaonZeroShort::MaxDCAnegV0) continue;
                    if (dca_pos_v0 > Default::KaonZeroShort::MaxDCAposV0) continue;
                }
                /* store */
                this_v0.idx = output.size();
                this_v0.neg = neg;
                this_v0.pos = pos;
                this_v0.v3 = v3_v0;
                this_v0.lv = lv_v0;
                this_v0.lv_neg = lv_neg;
                this_v0.lv_pos = lv_pos;
                output.push_back(this_v0);
            }  // end of loop over pos
        }      // end of loop over neg
        return output;
    };
    //
    auto CollectTrueInfo = [pdg_code_v0, pdg_code_neg, pdg_code_pos](const RVec<FoundV0> &found_v0s, cRVecUL track_mc_entry_, cRVecL mother_mc_entry_,
                                                                     RVecI pdg_code, cRVecI mc_is_secondary, cRVecI mc_is_signal,
                                                                     cRVecU mc_reaction_id) -> RVec<FoundV0_TrueInfo> {
        RVec<FoundV0_TrueInfo> output;
        output.reserve(found_v0s.size());
        FoundV0_TrueInfo this_v0_mc{};
        for (const auto &found_v0 : found_v0s) {
            /* defaults */
            this_v0_mc.entry = -1;
            this_v0_mc.pdg_code = 0;
            this_v0_mc.is_secondary = false;
            this_v0_mc.is_true = false;
            this_v0_mc.is_signal = false;
            this_v0_mc.reaction_id = 0;
            /* fill values */
            this_v0_mc.neg = track_mc_entry_[found_v0.neg];
            this_v0_mc.pos = track_mc_entry_[found_v0.pos];
            this_v0_mc.neg_pdg_code = pdg_code[this_v0_mc.neg];
            this_v0_mc.pos_pdg_code = pdg_code[this_v0_mc.pos];
            this_v0_mc.same_mother = mother_mc_entry_[this_v0_mc.neg] == mother_mc_entry_[this_v0_mc.pos];
            if (this_v0_mc.same_mother) {
                this_v0_mc.entry = mother_mc_entry_[this_v0_mc.neg];
                this_v0_mc.pdg_code = pdg_code[this_v0_mc.entry];
                this_v0_mc.is_secondary = mc_is_secondary[this_v0_mc.entry];
                this_v0_mc.is_true =
                    this_v0_mc.neg_pdg_code == pdg_code_neg && this_v0_mc.pos_pdg_code == pdg_code_pos && this_v0_mc.pdg_code == pdg_code_v0;
                if (this_v0_mc.is_true) {
                    this_v0_mc.is_signal = mc_is_signal[this_v0_mc.entry];
                    this_v0_mc.reaction_id = mc_reaction_id[this_v0_mc.entry];
                }
            }
            this_v0_mc.is_hybrid = !this_v0_mc.is_signal && ((mc_is_signal[this_v0_mc.neg] && !mc_is_signal[this_v0_mc.pos]) ||
                                                             (!mc_is_signal[this_v0_mc.neg] && mc_is_signal[this_v0_mc.pos]));
            output.emplace_back(this_v0_mc);
        }
        return output;
    };

    RNode df_V0s =
        df.Define("Found_" + v0_name, KF_V0Finder,                                              //
                  {neg_name + "_TrackEntry_", pos_name + "_TrackEntry_",                        //
                   "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",         //
                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",  //
                   "Track_CovMatrix", "MagneticField", "Event_PV"})
            .Define(v0_name + "_Neg_TrackEntry", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.neg; }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_TrackEntry", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.pos; }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_EsdIdx",
                    [](const RVec<FoundV0> &v0s, cRVecL track_esd_idx) {
                        return Map(v0s, [&track_esd_idx](const FoundV0 &v0) { return track_esd_idx[v0.neg]; });
                    },
                    {"Found_" + v0_name, "Track_EsdIdx"})
            .Define(v0_name + "_Pos_EsdIdx",
                    [](const RVec<FoundV0> &v0s, cRVecL track_esd_idx) {
                        return Map(v0s, [&track_esd_idx](const FoundV0 &v0) { return track_esd_idx[v0.pos]; });
                    },
                    {"Found_" + v0_name, "Track_EsdIdx"})
            .Define(v0_name + "_Px", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Py", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pz", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Px", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_neg.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Py", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_neg.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Pz", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_neg.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Px", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_pos.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Py", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_pos.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Pz", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv_pos.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Xv", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.v3.X(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Yv", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.v3.Y(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Zv", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.v3.Z(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Mass", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.lv.M(); }); },
                    {"Found_" + v0_name})  // DEBUG
            .Define(v0_name + "_Radius", [](const RVec<FoundV0> &v0s) { return Map(v0s, [](const FoundV0 &v0) { return v0.v3.Rho(); }); },
                    {"Found_" + v0_name});  // DEBUG

    if (Settings::IsMC) {
        df_V0s =
            df_V0s  //
                .Define(v0_name + "_TrueInfo", CollectTrueInfo,
                        {"Found_" + v0_name, "Track_McEntry_", "MC_Mother_McEntry_", "MC_PdgCode", "MC_IsSecondary", "MC_IsSignal", "MC_ReactionID"})
                .Define(v0_name + "_Neg_McEntry",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.neg; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_Pos_McEntry",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.pos; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_McEntry",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.entry; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_PdgCode",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.pdg_code; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_Neg_PdgCode",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.neg_pdg_code; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_Pos_PdgCode",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.pos_pdg_code; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_IsTrue",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.is_true; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_IsSignal",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.is_signal; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_IsSecondary",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.is_secondary; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_IsHybrid",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.is_hybrid; }); },
                        {v0_name + "_TrueInfo"})
                .Define(v0_name + "_ReactionID",
                        [](const RVec<FoundV0_TrueInfo> &v0s) { return Map(v0s, [](const FoundV0_TrueInfo &v0_mc) { return v0_mc.reaction_id; }); },
                        {v0_name + "_TrueInfo"});
    }

    return df_V0s;
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/*
 * Using Kalman Filter, reconstruct (anti)sexaquark candidates via the reaction channels
 * - `AntiSexaquark Neutron -> AntiLambda K0S`
 * - `??? -> Lambda K0S`
 */
RNode Manager::FindSexaquarks_TypeA(RNode df, Bool_t anti_channel) {

    const Double_t neutron_mass = TDatabasePDG::Instance()->GetParticle(2112)->Mass();
    const Double_t proton_mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    const Double_t pion_mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    std::string sexaquark_name;
    std::string v0a_name;
    Double_t v0a_neg_mass, v0a_pos_mass;
    if (!anti_channel) {
        sexaquark_name = "ASA";      // antisexaquark type a
        v0a_name = "AL";             // antilambda
        v0a_neg_mass = proton_mass;  // antiproton
        v0a_pos_mass = pion_mass;    // piplus
    } else {
        sexaquark_name = "BA";       // background type a
        v0a_name = "L";              // lambda
        v0a_neg_mass = pion_mass;    // piminus
        v0a_pos_mass = proton_mass;  // proton
    }
    std::string v0b_name = "K0S";  // kaonzeroshort

    auto KF_SexaquarkFinder = [neutron_mass, v0a_neg_mass, v0a_pos_mass, pion_mass](const RVec<FoundV0> &found_v0a, const RVec<FoundV0> &found_v0b,
                                                                                    cRVecI charge, const Float_t &magnetic_field,
                                                                                    const KFVertex &kf_pv) -> RVec<TypeA> {
        RVec<TypeA> output;
        TypeA this_sexa;
        /* protection */
        if (!found_v0a.size() || !found_v0b.size()) return output;
        /* one-time definitions */
        KFParticle::SetField(magnetic_field);
        XYZPoint v3_pv(kf_pv.GetX(), kf_pv.GetY(), kf_pv.GetZ());
        PxPyPzMVector lv_struck_neutron(0., 0., 0., neutron_mass);
        /* useful vars */
        KFParticle kf_v0a, kf_v0a_neg, kf_v0a_pos;
        KFParticle kf_v0b, kf_v0b_neg, kf_v0b_pos;
        PxPyPzMVector lv_v0a, lv_v0a_neg, lv_v0a_pos;
        PxPyPzMVector lv_v0b, lv_v0b_neg, lv_v0b_pos;
        PxPyPzMVector lv_sexa;
        PxPyPzMVector lv_sexa_asdecay;
        XYZPoint v3_sexa;
        /* auxiliary vars */
        Double_t radius, rapidity;
        Double_t cpa_wrt_pv, dca_la_sv, dca_k0s_sv;
        Double_t dca_la_neg_sv, dca_la_pos_sv, dca_k0s_neg_sv, dca_k0s_pos_sv;
        Double_t dca_btw_v0s;
        for (auto v0a : found_v0a) {
            for (auto v0b : found_v0b) {
                /* sanity check */
                std::set<ULong64_t> unique_track_entries = {v0a.neg, v0a.pos, v0b.neg, v0b.pos};
                if (unique_track_entries.size() != 4) continue;
                /* (debug) */
                // std::cout << "v0a = " << v0a.kf.GetX() << " " << v0a.kf.GetY() << " " << v0a.kf.GetZ() << std::endl;
                // std::cout << "v0b = " << v0b.kf.GetX() << " " << v0b.kf.GetY() << " " << v0b.kf.GetZ() << std::endl;
                /* geometry */
                /* -- copy KFParticle objects... */
                // kf_v0a = v0a.kf;
                // kf_v0a_neg = v0a.kf_neg;
                // kf_v0a_pos = v0a.kf_pos;
                // kf_v0b = v0b.kf;
                // kf_v0b_neg = v0b.kf_neg;
                // kf_v0b_pos = v0b.kf_pos;
                /* (debug) */
                // std::cout << "v0a (after copy) = " << kf_v0a.GetX() << " " << kf_v0a.GetY() << " " << kf_v0a.GetZ() << std::endl;
                // std::cout << "v0b (after copy) = " << kf_v0b.GetX() << " " << kf_v0b.GetY() << " " << kf_v0b.GetZ() << std::endl;
                //
                KFParticle kf_sexa;
                kf_sexa.AddDaughter(kf_v0a);
                kf_sexa.AddDaughter(kf_v0b);
                kf_sexa.TransportToDecayVertex();
                if (TMath::IsNaN(kf_sexa.GetX()) || TMath::IsNaN(kf_sexa.GetY()) || TMath::IsNaN(kf_sexa.GetZ())) continue;
                v3_sexa.SetCoordinates(kf_sexa.GetX(), kf_sexa.GetY(), kf_sexa.GetZ());
                /* (debug) */
                std::cout << "v0a (after added) = " << kf_v0a.GetX() << " " << kf_v0a.GetY() << " " << kf_v0a.GetZ() << '\n';
                std::cout << "v0b (after added) = " << kf_v0b.GetX() << " " << kf_v0b.GetY() << " " << kf_v0b.GetZ() << '\n';
                std::cout << "sexa = " << kf_sexa.GetX() << " " << kf_sexa.GetY() << " " << kf_sexa.GetZ() << '\n';
                /* transport V0s */
                kf_v0a.SetProductionVertex(kf_sexa);
                kf_v0b.SetProductionVertex(kf_sexa);
                kf_v0a.TransportToProductionVertex();
                kf_v0b.TransportToProductionVertex();
                /* (debug) */
                // std::cout << "v0a (after supposed transport) = " << kf_v0a.GetX() << " " << kf_v0a.GetY() << " " << kf_v0a.GetZ() << std::endl;
                // std::cout << "v0b (after supposed transport) = " << kf_v0b.GetX() << " " << kf_v0b.GetY() << " " << kf_v0b.GetZ() << std::endl;
                /* transport V0s daughters */
                // kf_v0a_neg = Math::TransportKFParticle(magnetic_field, kf_v0a_neg, kf_v0a_pos, v0a_neg_mass, charge[v0a.neg]);
                // kf_v0a_pos = Math::TransportKFParticle(magnetic_field, kf_v0a_pos, kf_v0a_neg, v0a_pos_mass, charge[v0a.pos]);
                // kf_v0b_neg = Math::TransportKFParticle(magnetic_field, kf_v0b_neg, kf_v0b_pos, pion_mass, charge[v0b.neg]);
                // kf_v0b_pos = Math::TransportKFParticle(magnetic_field, kf_v0b_pos, kf_v0b_neg, pion_mass, charge[v0b.pos]);
                /* kinematics */
                lv_v0a_neg.SetCoordinates(kf_v0a_neg.Px(), kf_v0a_neg.Py(), kf_v0a_neg.Pz(), v0a_neg_mass);
                lv_v0a_pos.SetCoordinates(kf_v0a_pos.Px(), kf_v0a_pos.Py(), kf_v0a_pos.Pz(), v0a_pos_mass);
                lv_v0b_neg.SetCoordinates(kf_v0b_neg.Px(), kf_v0b_neg.Py(), kf_v0b_neg.Pz(), pion_mass);
                lv_v0b_pos.SetCoordinates(kf_v0b_pos.Px(), kf_v0b_pos.Py(), kf_v0b_pos.Pz(), pion_mass);
                lv_v0a = lv_v0a_neg + lv_v0a_pos;
                lv_v0b = lv_v0b_neg + lv_v0b_pos;
                lv_sexa = lv_v0a + lv_v0b - lv_struck_neutron;
                lv_sexa_asdecay = lv_v0a + lv_v0b;
                /* calculate properties */
                radius = v3_sexa.Rho();
                rapidity = lv_sexa.Rapidity();
                cpa_wrt_pv = Math::CosinePointingAngle(lv_sexa.Vect(), v3_sexa, v3_pv);
                dca_la_sv = TMath::Abs(kf_v0a.GetDistanceFromVertex(kf_sexa));
                dca_la_neg_sv = TMath::Abs(kf_v0a_neg.GetDistanceFromVertex(kf_sexa));
                dca_la_pos_sv = TMath::Abs(kf_v0a_pos.GetDistanceFromVertex(kf_sexa));
                dca_k0s_sv = TMath::Abs(kf_v0b.GetDistanceFromVertex(kf_sexa));
                dca_k0s_neg_sv = TMath::Abs(kf_v0b_neg.GetDistanceFromVertex(kf_sexa));
                dca_k0s_pos_sv = TMath::Abs(kf_v0b_pos.GetDistanceFromVertex(kf_sexa));
                dca_btw_v0s = TMath::Abs(kf_v0a.GetDistanceFromParticle(kf_v0b));
                /* (debug) */
                std::cout << "radius=" << radius << " rapidity=" << rapidity << " cpa_wrt_pv=" << cpa_wrt_pv << " dca_la_sv=" << dca_la_sv
                          << " dca_la_neg_sv=" << dca_la_neg_sv << " dca_la_pos_sv=" << dca_la_pos_sv << " dca_k0s_sv=" << dca_k0s_sv
                          << " dca_k0s_neg_sv=" << dca_k0s_neg_sv << " dca_k0s_pos_sv=" << dca_k0s_pos_sv << " dca_btw_v0s=" << dca_btw_v0s << '\n';
                /* cuts */
                if (radius < Default::ChannelA::MinRadius || radius > 200.) continue;
                if (TMath::Abs(rapidity) > Default::ChannelA::AbsMaxRapidity) continue;
                if (cpa_wrt_pv < Default::ChannelA::MinCPAwrtPV || cpa_wrt_pv > Default::ChannelA::MaxCPAwrtPV) continue;
                if (dca_la_sv > Default::ChannelA::MaxDCALaSV) continue;
                if (dca_la_neg_sv > Default::ChannelA::MaxDCALaNegSV) continue;
                if (dca_la_pos_sv > Default::ChannelA::MaxDCALaPosSV) continue;
                if (dca_k0s_sv > Default::ChannelA::MaxDCAK0SV) continue;
                if (dca_k0s_neg_sv > Default::ChannelA::MaxDCAK0NegSV) continue;
                if (dca_k0s_pos_sv > Default::ChannelA::MaxDCAK0PosSV) continue;
                if (dca_btw_v0s > Default::ChannelA::MaxDCAbtwV0s) continue;
                std::cout << "a candidate passed all cuts!" << '\n';
                /* store */
                this_sexa.v0a = v0a;
                this_sexa.v0b = v0b;
                this_sexa.kf = kf_sexa;
                this_sexa.lv = lv_sexa;
                this_sexa.lv_asdecay = lv_sexa_asdecay;
                output.push_back(this_sexa);
            }  // end of loop over k0s
        }      // end of loop over (anti)lambdas
        return output;
    };
    //
    auto CollectTrueInfo = [](const RVec<TypeA> &found, const RVec<FoundV0_TrueInfo> &v0a_mc,
                              const RVec<FoundV0_TrueInfo> &v0b_mc) -> RVec<TypeA_TrueInfo> {
        RVec<TypeA_TrueInfo> output;
        output.reserve(found.size());
        TypeA_TrueInfo this_sexa_mc{};
        /* auxiliar vars */
        Bool_t same_reaction_id;
        for (const auto &found_sexa : found) {
            /* defaults */
            this_sexa_mc.reaction_id = 0;
            this_sexa_mc.is_signal = false;
            /* fill values */
            this_sexa_mc.v0a_mc = v0a_mc[found_sexa.v0a.idx];
            this_sexa_mc.v0b_mc = v0b_mc[found_sexa.v0b.idx];
            same_reaction_id = this_sexa_mc.v0a_mc.reaction_id == this_sexa_mc.v0b_mc.reaction_id;
            if (same_reaction_id) {
                this_sexa_mc.reaction_id = this_sexa_mc.v0a_mc.reaction_id;
                this_sexa_mc.is_signal = this_sexa_mc.v0a_mc.is_signal && this_sexa_mc.v0b_mc.is_signal;
            }
            this_sexa_mc.is_hybrid = (!this_sexa_mc.is_signal && ((this_sexa_mc.v0a_mc.is_signal && !this_sexa_mc.v0b_mc.is_signal) ||
                                                                  (!this_sexa_mc.v0a_mc.is_signal && this_sexa_mc.v0b_mc.is_signal))) ||
                                     this_sexa_mc.v0a_mc.is_hybrid || this_sexa_mc.v0b_mc.is_hybrid;
            output.emplace_back(this_sexa_mc);
        }
        return output;
    };

    RNode df_Sexaquarks = df.Define("Found_" + sexaquark_name, KF_SexaquarkFinder,  //
                                    {"Found_" + v0a_name, "Found_" + v0b_name, "Track_Charge", "MagneticField", "Event_KF_PV"})
                              .Define(sexaquark_name + "_Mass",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.lv.M(); }); },
                                      {"Found_" + sexaquark_name})
                              .Define(sexaquark_name + "_Pt",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.lv.Pt(); }); },
                                      {"Found_" + sexaquark_name})
                              .Define(sexaquark_name + "_Xv",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.kf.GetX(); }); },
                                      {"Found_" + sexaquark_name})
                              .Define(sexaquark_name + "_Yv",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.kf.GetY(); }); },
                                      {"Found_" + sexaquark_name})
                              .Define(sexaquark_name + "_Radius",
                                      [](const RVec<TypeA> &sexaquarks) {
                                          return Map(sexaquarks, [](const TypeA &sexa) {
                                              return TMath::Sqrt(sexa.kf.GetX() * sexa.kf.GetX() + sexa.kf.GetY() * sexa.kf.GetY());
                                          });
                                      },
                                      {"Found_" + sexaquark_name});

    if (Settings::IsMC) {
        df_Sexaquarks =
            df_Sexaquarks                                               //
                .Define(sexaquark_name + "_TrueInfo", CollectTrueInfo,  //
                        {"Found_" + sexaquark_name, v0a_name + "_TrueInfo", v0b_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0a_Neg_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0a_mc.neg; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0a_Pos_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0a_mc.pos; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0b_Neg_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0b_mc.neg; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0b_Pos_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0b_mc.pos; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0a_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0a_mc.entry; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_V0b_McEntry",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.v0b_mc.entry; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(
                    sexaquark_name + "_IsSignal",
                    [](const RVec<TypeA_TrueInfo> &sexaquarks) { return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.is_signal; }); },
                    {sexaquark_name + "_TrueInfo"})
                .Define(sexaquark_name + "_ReactionID",
                        [](const RVec<TypeA_TrueInfo> &sexaquarks) {
                            return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.reaction_id; });
                        },
                        {sexaquark_name + "_TrueInfo"})
                .Define(
                    sexaquark_name + "_IsHybrid",
                    [](const RVec<TypeA_TrueInfo> &sexaquarks) { return Map(sexaquarks, [](const TypeA_TrueInfo &sexa) { return sexa.is_hybrid; }); },
                    {sexaquark_name + "_TrueInfo"});
    }

    return df_Sexaquarks;
}

RNode Manager::FindSexaquarks_TypeDE(RNode df, Bool_t anti_channel) {
    //
    return df;
}

RNode Manager::FindSexaquarks_TypeH(RNode df, Bool_t anti_channel) {
    //
    return df;
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
    // --- //
    size_t n_events = n_tracks->size();
    for (size_t event_i = 0; event_i < n_events; event_i++) {
        std::cout << "# Event " << event_i << '\n';
        // --- //
        std::cout << "## Summary" << '\n';
        printf(">> N Injected = %lu\n", (*n_injected)[event_i]);
        printf(">> N Findable Sexaquarks = %d\n", Sum((*ev_injected_is_findable)[event_i]));
        printf(">> N MC = %lu\n", n_mc->at(event_i));
        printf(">> N Findable V0s = %d\n", Sum((*ev_mc_truev0_isfindable)[event_i]));
        // --- //
        std::cout << "## (MC) Particles" << '\n';
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            printf("mc_entry=%lu, pdg_code=%d, ancestor_entry=%ld, mother_entry=%ld, neg_entry=%ld, pos_entry=%ld\n",  //
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_ancestor_entry)[event_i][mc_entry],
                   (*ev_mc_mother_entry)[event_i][mc_entry], (*ev_mc_neg_entry)[event_i][mc_entry], (*ev_mc_pos_entry)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## (MC) V0s" << '\n';
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            if (((*ev_mc_pdgcode)[event_i][mc_entry] != PdgCode::KaonZeroShort && TMath::Abs((*ev_mc_pdgcode)[event_i][mc_entry]) != PdgCode::Lambda))
                continue;
            printf("mc_entry=%lu, pdg_code=%d, is_signal=%d, neg_entry=%ld, pos_entry=%ld, neg_gt=%d, pos_gt=%d, neg_pdg=%d, pos_pdg=%d, xv=%f\n",
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_is_signal)[event_i][mc_entry], (*ev_mc_neg_entry)[event_i][mc_entry],
                   (*ev_mc_pos_entry)[event_i][mc_entry], (*ev_mc_neg_gottracked)[event_i][mc_entry], (*ev_mc_pos_gottracked)[event_i][mc_entry],
                   (*ev_mc_neg_pdgcode)[event_i][mc_entry], (*ev_mc_pos_pdgcode)[event_i][mc_entry], (*ev_mc_xv)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## (MC) Charged Particles" << '\n';
        for (size_t mc_entry = 0; mc_entry < (*n_mc)[event_i]; mc_entry++) {
            if ((*ev_mc_is_charged)[event_i][mc_entry] == 0) continue;
            printf("mc_entry=%ld, pdg_code=%d, is_signal=%d, mother_entry=%ld, mother_pdg=%d, got_tracked=%d\n",  //
                   mc_entry, (*ev_mc_pdgcode)[event_i][mc_entry], (*ev_mc_is_signal)[event_i][mc_entry], (*ev_mc_mother_entry)[event_i][mc_entry],
                   (*ev_mc_mother_pdgcode)[event_i][mc_entry], (*ev_mc_got_tracked)[event_i][mc_entry]);
        }
        // --- //
        std::cout << "## Injected" << '\n';
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
        std::cout << "## Tracks" << '\n';
        for (size_t track_entry = 0; track_entry < (*n_tracks)[event_i]; track_entry++) {
            printf("track_entry=%ld, mc_entry=%lu, nsigmaproton=%f, nsigmapion=%f, pdg_code=%d, is_signal=%d\n",  //
                   track_entry, (*ev_track_mc_entry)[event_i][track_entry], (*ev_track_nsigmaproton)[event_i][track_entry],
                   (*ev_track_nsigmapion)[event_i][track_entry], (*ev_track_pdgcode)[event_i][track_entry],
                   (*ev_track_issignal)[event_i][track_entry]);
        }
    }  // end of loop over events
    InfoF("NRuns = %u", df.GetNRuns());
}

void Manager::EndOfAnalysis(RNode df) {
    std::vector<std::string> column_list;
    /* Add V0s Properties */
    for (const auto &v0_name : fAnalyzed_V0sNames) {
        column_list.insert(column_list.end(), {v0_name + "_Neg_TrackEntry", v0_name + "_Pos_TrackEntry",       //
                                               v0_name + "_Neg_EsdIdx", v0_name + "_Pos_EsdIdx",               //
                                               v0_name + "_Px", v0_name + "_Py", v0_name + "_Pz",              //
                                               v0_name + "_Neg_Px", v0_name + "_Neg_Py", v0_name + "_Neg_Pz",  //
                                               v0_name + "_Pos_Px", v0_name + "_Pos_Py", v0_name + "_Pos_Pz",  //
                                               v0_name + "_Xv", v0_name + "_Yv", v0_name + "_Zv",              //
                                               v0_name + "_Mass", v0_name + "_Radius"});
        if (Settings::IsMC) {
            column_list.insert(column_list.end(), {v0_name + "_Neg_McEntry", v0_name + "_Pos_McEntry", v0_name + "_McEntry", v0_name + "_PdgCode",
                                                   v0_name + "_Neg_PdgCode", v0_name + "_Pos_PdgCode", v0_name + "_IsTrue", v0_name + "_IsSignal",
                                                   v0_name + "_IsSecondary", v0_name + "_IsHybrid", v0_name + "_ReactionID"});
        }
    }
    /* Add (Anti)Sexaquarks Properties */
    /*
    std::string sexaquark_name = "ASA";
    column_list.insert(column_list.end(), {sexaquark_name + "_Mass", sexaquark_name + "_Pt", sexaquark_name + "_Xv", sexaquark_name + "_Yv",
                                           sexaquark_name + "_Radius"});
    if (Settings::IsMC) {
        column_list.insert(column_list.end(),
                           {sexaquark_name + "_V0a_Neg_McEntry", sexaquark_name + "_V0a_Pos_McEntry", sexaquark_name + "_V0b_Neg_McEntry",
                            sexaquark_name + "_V0b_Pos_McEntry", sexaquark_name + "_V0a_McEntry", sexaquark_name + "_V0b_McEntry",
                            sexaquark_name + "_IsSignal", sexaquark_name + "_ReactionID", sexaquark_name + "_IsHybrid"});
    }
    */
    /* Write to disk */
    df.Snapshot("Events", Settings::PathOutputFile, column_list);
    InfoF("TFile %s has been written", Settings::PathOutputFile.c_str());
    /* Print info */
    InfoF("NRuns = %u", df.GetNRuns());
}

}  // namespace Tree2Sexaquark::Analysis
