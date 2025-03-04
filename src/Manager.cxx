#include "Analysis/Manager.hxx"

#include <cstddef>
#include <set>

#include "RtypesCore.h"
#include "TDatabasePDG.h"

#include "Analysis/Settings.hxx"
#include "Cuts/Default.hxx"
#include "Math/Common.hxx"
#include "Math/KalmanFilter.hxx"

#define READ_ESD_INDICES 0

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
    /* init map */
    fParticleName_[PdgCode::PiMinus] = "PiMinus";
    fParticleName_[PdgCode::PiPlus] = "PiPlus";
    fParticleName_[PdgCode::NegKaon] = "NegKaon";
    fParticleName_[PdgCode::PosKaon] = "PosKaon";
    fParticleName_[PdgCode::AntiProton] = "AntiProton";
    fParticleName_[PdgCode::Proton] = "Proton";
    fParticleName_[PdgCode::AntiNeutron] = "AntiNeutron";
    fParticleName_[PdgCode::Neutron] = "Neutron";
    fParticleName_[PdgCode::AntiLambda] = "AL";
    fParticleName_[PdgCode::Lambda] = "L";
    fParticleName_[PdgCode::KaonZeroShort] = "K0S";
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

/*                  */
/**  MC Particles  **/
/*** ============ ***/

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
            /* # Trivial definitions */
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

/*            */
/**  Tracks  **/
/*** ====== ***/

RVec<KF_Track> Manager::Tracks_KF_Creator(Double_t pdg_mass, cRVecUL entries_,                                    //
                                          cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,          //
                                          cRVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,  //
                                          const RVec<RVecF> &cov_matrix, const Float_t &magnetic_field) {
    KFParticle::SetField(magnetic_field);
    RVec<KF_Track> output;
    output.reserve(entries_.size());
    KF_Track track{};
    for (auto entry : entries_) {
        track.entry = entry;
        track.kf = Math::CreateKFParticle(px[entry], py[entry], pz[entry], x[entry], y[entry], z[entry],          //
                                          charge[entry], alpha[entry], snp[entry], tgl[entry], signed1pt[entry],  //
                                          cov_matrix[entry], pdg_mass);
        track.lv = PxPyPzMVector(px[entry], py[entry], pz[entry], pdg_mass);
        output.emplace_back(track);
    }
    return output;
}

RVec<MC_Track> Manager::Linked_MC_Creator(cRVecUL linked_mc_entries_, cRVecL mother_mc_entries_,  //
                                          cRVecI pdg_code, cRVecI is_secondary, cRVecI is_signal, cRVecU reaction_id) {
    //
    RVec<MC_Track> output;
    output.reserve(linked_mc_entries_.size());
    MC_Track particle{};
    for (auto mc_entry : linked_mc_entries_) {
        particle.mother_entry = mother_mc_entries_[mc_entry];
        particle.pdg_code = pdg_code[mc_entry];
        particle.is_secondary = is_secondary[mc_entry];
        particle.is_signal = is_signal[mc_entry];
        particle.reaction_id = reaction_id[mc_entry];
        output.emplace_back(particle);
    }
    return output;
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
    //
    auto CallKFProtonCreator = [mass_proton](cRVecUL entries_, cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,  //
                                             cRVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,            //
                                             const RVec<RVecF> &cov_matrix, const Float_t &magnetic_field) -> RVec<KF_Track> {
        return Tracks_KF_Creator(mass_proton, entries_, px, py, pz, x, y, z, charge, alpha, snp, tgl, signed1pt, cov_matrix, magnetic_field);
    };
    //
    auto CallKFKaonCreator = [mass_kaon](cRVecUL entries_, cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,  //
                                         cRVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,            //
                                         const RVec<RVecF> &cov_matrix, const Float_t &magnetic_field) -> RVec<KF_Track> {
        return Tracks_KF_Creator(mass_kaon, entries_, px, py, pz, x, y, z, charge, alpha, snp, tgl, signed1pt, cov_matrix, magnetic_field);
    };
    //
    auto CallKFPionCreator = [mass_pion](cRVecUL entries_, cRVecF px, cRVecF py, cRVecF pz, cRVecF x, cRVecF y, cRVecF z,  //
                                         cRVecI charge, cRVecF alpha, cRVecF snp, cRVecF tgl, cRVecF signed1pt,            //
                                         const RVec<RVecF> &cov_matrix, const Float_t &magnetic_field) -> RVec<KF_Track> {
        return Tracks_KF_Creator(mass_pion, entries_, px, py, pz, x, y, z, charge, alpha, snp, tgl, signed1pt, cov_matrix, magnetic_field);
    };

    RNode df_Tracks = df.Define("N_Tracks", "Track_Px.size()")
                          /* # PID */
                          .Define("Track_PID_AsAntiProton", "abs(Track_NSigmaProton) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsProton", "abs(Track_NSigmaProton) < 3. && Track_Charge > 0")
                          .Define("Track_PID_AsNegKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsPosKaon", "abs(Track_NSigmaKaon) < 3. && Track_Charge > 0")
                          .Define("Track_PID_AsPiMinus", "abs(Track_NSigmaPion) < 3. && Track_Charge < 0")
                          .Define("Track_PID_AsPiPlus", "abs(Track_NSigmaPion) < 3. && Track_Charge > 0")
                          /*  */
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
                                   "Track_Sigma1PtY", "Track_Sigma1PtZ", "Track_Sigma1PtSnp", "Track_Sigma1PtTgl", "Track_Sigma1Pt2"})
                          /* # KF */
                          .Define("KF_AntiProton", CallKFProtonCreator,                                         //
                                  {"AntiProton_TrackEntry_",                                                    //
                                   "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",         //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",  //
                                   "Track_CovMatrix", "MagneticField"})
                          .Define("KF_Proton", CallKFProtonCreator,
                                  {"Proton_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",  //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",                 //
                                   "Track_CovMatrix", "MagneticField"})
                          .Define("KF_NegKaon", CallKFKaonCreator,
                                  {"NegKaon_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",  //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",                  //
                                   "Track_CovMatrix", "MagneticField"})
                          .Define("KF_PosKaon", CallKFKaonCreator,
                                  {"PosKaon_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",  //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",                  //
                                   "Track_CovMatrix", "MagneticField"})
                          .Define("KF_PiMinus", CallKFPionCreator,
                                  {"PiMinus_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",  //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",                  //
                                   "Track_CovMatrix", "MagneticField"})
                          .Define("KF_PiPlus", CallKFPionCreator,
                                  {"PiPlus_TrackEntry_", "Track_Px", "Track_Py", "Track_Pz", "Track_X", "Track_Y", "Track_Z",  //
                                   "Track_Charge", "Track_Alpha", "Track_Snp", "Track_Tgl", "Track_Signed1Pt",                 //
                                   "Track_CovMatrix", "MagneticField"});

    if (Settings::IsMC) {
        df_Tracks = df_Tracks                                    //
                        .Define("MC_Tracks", Linked_MC_Creator,  //
                                {"Track_McEntry_", "MC_Mother_McEntry_", "MC_PdgCode", "MC_IsSecondary", "MC_IsSignal", "MC_ReactionID"});
    }

    return df_Tracks;
}

/*         */
/**  V0s  **/
/*** === ***/

/*
 * Using Kalman Filter, find V0s.
 */
RVec<KF_V0> Manager::V0s_KF_Finder(PdgCode pdg_code_v0, Double_t neg_mass, Double_t pos_mass,           //
                                   const RVec<KF_Track> &neg_tracks, const RVec<KF_Track> &pos_tracks,  //
                                   const Float_t &magnetic_field, const XYZPoint &v3_pv) {
    RVec<KF_V0> output;
    KF_V0 this_v0;
    KF_Track cp_neg, cp_pos;
    /* protection */
    if (!neg_tracks.size() || !pos_tracks.size()) return output;
    /* one-time definitions */
    KFParticle::SetField(magnetic_field);
    /* loop over all pairs */
    for (const auto &neg : neg_tracks) {
        for (const auto &pos : pos_tracks) {
            /* sanity check */
            if (neg.entry == pos.entry) continue;
            cp_neg = neg;
            cp_pos = pos;
            /* fit */
            KFParticle kf_v0(cp_neg.kf, cp_pos.kf);
            //
            kf_v0.TransportToDecayVertex();
            cp_neg.kf.SetProductionVertex(kf_v0);
            cp_pos.kf.SetProductionVertex(kf_v0);
            /* transport daughters to V0 vertex */
            cp_neg.kf.TransportToProductionVertex();
            cp_pos.kf.TransportToProductionVertex();
            cp_neg.lv.SetCoordinates(cp_neg.kf.Px(), cp_neg.kf.Py(), cp_neg.kf.Pz(), neg_mass);
            cp_pos.lv.SetCoordinates(cp_pos.kf.Px(), cp_pos.kf.Py(), cp_pos.kf.Pz(), pos_mass);
            /* fill struct */
            this_v0.idx = output.size();
            this_v0.neg = cp_neg;
            this_v0.pos = cp_pos;
            this_v0.kf = kf_v0;
            this_v0.lv = cp_neg.lv + cp_pos.lv;
            /* apply cuts and store */
            if (!V0s_PassesCuts(pdg_code_v0, this_v0, v3_pv)) continue;
            output.push_back(this_v0);
        }  // end of loop over pos
    }      // end of loop over neg
    return output;
};

RVec<MC_V0> Manager::V0s_TrueInfoCollector(PdgCode pdg_code_v0, PdgCode pdg_code_neg, PdgCode pdg_code_pos,  //
                                           const RVec<KF_V0> &found_v0s, const RVec<MC_Track> &linked_mc,    //
                                           cRVecI pdg_code, cRVecI is_signal, cRVecU reaction_id) {
    RVec<MC_V0> output;
    output.reserve(found_v0s.size());
    MC_V0 mc_v0{};
    for (const auto &found_v0 : found_v0s) {
        mc_v0.mc_entry = -1;
        mc_v0.pdg_code = 0;
        mc_v0.is_true = false;
        mc_v0.is_signal = false;
        mc_v0.reaction_id = 0;
        /* fill */
        mc_v0.neg = linked_mc[found_v0.neg.entry];
        mc_v0.pos = linked_mc[found_v0.pos.entry];
        mc_v0.has_mc = mc_v0.neg.mother_entry != -1 && mc_v0.neg.mother_entry == mc_v0.pos.mother_entry;
        if (mc_v0.has_mc) {
            mc_v0.mc_entry = (Long64_t)mc_v0.neg.mother_entry;
            mc_v0.pdg_code = pdg_code[mc_v0.mc_entry];
            mc_v0.is_true = mc_v0.neg.pdg_code == pdg_code_neg && mc_v0.pos.pdg_code == pdg_code_pos && mc_v0.pdg_code == pdg_code_v0;
            if (mc_v0.is_true) {
                mc_v0.is_signal = is_signal[mc_v0.mc_entry];
                mc_v0.reaction_id = reaction_id[mc_v0.mc_entry];
            }
        }
        mc_v0.is_hybrid = !mc_v0.is_signal &&  //
                          ((mc_v0.neg.is_signal && !mc_v0.pos.is_signal) || (!mc_v0.neg.is_signal && mc_v0.pos.is_signal));
        output.emplace_back(mc_v0);
    }
    return output;
};

Bool_t Manager::V0s_PassesCuts(PdgCode pdg_code_v0, const KF_V0 &v0, const XYZPoint &v3_pv) {
    //
    XYZPoint v3_v0(v0.kf.GetX(), v0.kf.GetY(), v0.kf.GetZ());
    /* -- (anti)lambda cuts */
    if (TMath::Abs(pdg_code_v0) == PdgCode::Lambda) {
        if (TMath::Abs(v0.lv.Eta()) > Default::Lambda::AbsMaxEta) return kFALSE;
        if (v0.lv.Pt() < Default::Lambda::MinPt) return kFALSE;
        if (v0.lv.M() < Default::Lambda::MinMass || v0.lv.M() > Default::Lambda::MaxMass) return kFALSE;
        if (v3_v0.Rho() < Default::Lambda::MinRadius) return kFALSE;  // NOTE: could set upper limit
        if ((v3_v0 - v3_pv).R() < Default::Lambda::MinDistFromPV) return kFALSE;
        /*  */ Double_t cpa_wrt_pv = Math::CosinePointingAngle(v0.lv.Vect(), v3_v0, v3_pv);
        if (cpa_wrt_pv < Default::Lambda::MinCPAwrtPV || cpa_wrt_pv > Default::Lambda::MaxCPAwrtPV) return kFALSE;
        /*  */ Double_t dca_wrt_pv = Math::LinePointDCA(v0.lv.Vect(), v3_v0, v3_pv);
        if (dca_wrt_pv < Default::Lambda::MinDCAwrtPV) return kFALSE;
        /*  */ Double_t dca_btw_dau = TMath::Abs((Double_t)v0.neg.kf.GetDistanceFromParticle(v0.pos.kf));
        if (dca_btw_dau > Default::Lambda::MaxDCAbtwDau) return kFALSE;
        /*  */ Double_t dca_neg_v0 = TMath::Abs((Double_t)v0.neg.kf.GetDistanceFromVertex(v0.kf));
        if (dca_neg_v0 > Default::Lambda::MaxDCAnegV0) return kFALSE;
        /*  */ Double_t dca_pos_v0 = TMath::Abs((Double_t)v0.pos.kf.GetDistanceFromVertex(v0.kf));
        if (dca_pos_v0 > Default::Lambda::MaxDCAposV0) return kFALSE;
        /*  */ Double_t arm_qt = Math::ArmenterosQt(v0.lv.Vect(), v0.neg.lv.Vect());
        /*  */ Double_t arm_alpha = Math::ArmenterosAlpha(v0.lv.Vect(), v0.neg.lv.Vect(), v0.pos.lv.Vect());
        /*  */ Double_t qt_over_alpha = arm_qt / TMath::Abs(arm_alpha);
        if (qt_over_alpha > Default::Lambda::AbsMaxArmQtOverAlpha) return kFALSE;
    }
    /* -- k0s cuts */
    if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        if (TMath::Abs(v0.lv.Eta()) > Default::KaonZeroShort::AbsMaxEta) return kFALSE;
        if (v0.lv.Pt() < Default::KaonZeroShort::MinPt) return kFALSE;
        if (v0.lv.M() < Default::KaonZeroShort::MinMass || v0.lv.M() > Default::KaonZeroShort::MaxMass) return kFALSE;
        if (v3_v0.Rho() < Default::KaonZeroShort::MinRadius) return kFALSE;  // NOTE: could set upper limit
        if ((v3_v0 - v3_pv).R() < Default::KaonZeroShort::MinDistFromPV ||   //
            (v3_v0 - v3_pv).R() > Default::KaonZeroShort::MaxDistFromPV)
            return kFALSE;
        /*  */ Double_t cpa_wrt_pv = Math::CosinePointingAngle(v0.lv.Vect(), v3_v0, v3_pv);
        if (cpa_wrt_pv < Default::KaonZeroShort::MinCPAwrtPV || cpa_wrt_pv > Default::KaonZeroShort::MaxCPAwrtPV) return kFALSE;
        /*  */ Double_t dca_wrt_pv = Math::LinePointDCA(v0.lv.Vect(), v3_v0, v3_pv);
        if (dca_wrt_pv < Default::KaonZeroShort::MinDCAwrtPV) return kFALSE;
        /*  */ Double_t dca_btw_dau = TMath::Abs((Double_t)v0.neg.kf.GetDistanceFromParticle(v0.pos.kf));
        if (dca_btw_dau > Default::KaonZeroShort::MaxDCAbtwDau) return kFALSE;
        /*  */ Double_t dca_neg_v0 = TMath::Abs((Double_t)v0.neg.kf.GetDistanceFromVertex(v0.kf));
        if (dca_neg_v0 > Default::KaonZeroShort::MaxDCAnegV0) return kFALSE;
        /*  */ Double_t dca_pos_v0 = TMath::Abs((Double_t)v0.pos.kf.GetDistanceFromVertex(v0.kf));
        if (dca_pos_v0 > Default::KaonZeroShort::MaxDCAposV0) return kFALSE;
    }
    /*  */
    return kTRUE;
}

RNode Manager::FindV0s(RNode df, PdgCode pdg_code_v0, PdgCode pdg_code_neg, PdgCode pdg_code_pos) {
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
    auto CallV0Finder = [pdg_code_v0, neg_mass, pos_mass](const RVec<KF_Track> &neg_tracks, const RVec<KF_Track> &pos_tracks,  //
                                                          const Float_t &magnetic_field, const XYZPoint &v3_pv) -> RVec<KF_V0> {
        return V0s_KF_Finder(pdg_code_v0, neg_mass, pos_mass, neg_tracks, pos_tracks, magnetic_field, v3_pv);
    };
    //
    auto CallTrueCollector = [pdg_code_v0, pdg_code_neg, pdg_code_pos](const RVec<KF_V0> &found_v0s, const RVec<MC_Track> &linked_mc,  //
                                                                       cRVecI pdg_code, cRVecI is_signal, cRVecU reaction_id) -> RVec<MC_V0> {
        return V0s_TrueInfoCollector(pdg_code_v0, pdg_code_neg, pdg_code_pos,  //
                                     found_v0s, linked_mc,                     //
                                     pdg_code, is_signal, reaction_id);
    };

    RNode df_V0s =
        df.Define("Found_" + v0_name, CallV0Finder, {"KF_" + neg_name, "KF_" + pos_name, "MagneticField", "Event_PV"})
            .Define(v0_name + "_Neg_TrackEntry", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.neg.entry; }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_TrackEntry", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.pos.entry; }); },
                    {"Found_" + v0_name})
#if READ_ESD_INDICES
            .Define(v0_name + "_Neg_EsdIdx",
                    [](const RVec<KF_V0> &v0s, cRVecL esd_idx) { return Map(v0s, [&esd_idx](const KF_V0 &v0) { return esd_idx[v0.neg.entry]; }); },
                    {"Found_" + v0_name, "Track_EsdIdx"})
            .Define(v0_name + "_Pos_EsdIdx",
                    [](const RVec<KF_V0> &v0s, cRVecL esd_idx) { return Map(v0s, [&esd_idx](const KF_V0 &v0) { return esd_idx[v0.pos.entry]; }); },
                    {"Found_" + v0_name, "Track_EsdIdx"})
#endif
            .Define(v0_name + "_Px", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.lv.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Py", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.lv.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pz", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.lv.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Px", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.neg.lv.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Py", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.neg.lv.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Neg_Pz", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.neg.lv.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Px", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.pos.lv.Px(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Py", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.pos.lv.Py(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Pos_Pz", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.pos.lv.Pz(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Xv", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.kf.GetX(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Yv", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.kf.GetY(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Zv", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.kf.GetZ(); }); },
                    {"Found_" + v0_name})
            .Define(v0_name + "_Mass", [](const RVec<KF_V0> &v0s) { return Map(v0s, [](const KF_V0 &v0) { return v0.lv.M(); }); },
                    {"Found_" + v0_name})  // DEBUG
            .Define(v0_name + "_Radius",
                    [](const RVec<KF_V0> &v0s) {
                        return Map(v0s, [](const KF_V0 &v0) { return TMath::Sqrt(v0.kf.GetX() * v0.kf.GetX() + v0.kf.GetY() * v0.kf.GetY()); });
                    },
                    {"Found_" + v0_name});  // DEBUG

    if (Settings::IsMC) {
        df_V0s =
            df_V0s                                           //
                .Define("MC_" + v0_name, CallTrueCollector,  //
                        {"Found_" + v0_name, "MC_Tracks", "MC_PdgCode", "MC_IsSignal", "MC_ReactionID"})
                .Define(v0_name + "_Neg_McEntry", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.neg.entry; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_Pos_McEntry", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.pos.entry; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_McEntry", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.mc_entry; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_PdgCode", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.pdg_code; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_Neg_PdgCode",
                        [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.neg.pdg_code; }); }, {"MC_" + v0_name})
                .Define(v0_name + "_Pos_PdgCode",
                        [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.pos.pdg_code; }); }, {"MC_" + v0_name})
                .Define(v0_name + "_IsTrue", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.is_true; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_IsSignal", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.is_signal; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_IsHybrid", [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.is_hybrid; }); },
                        {"MC_" + v0_name})
                .Define(v0_name + "_ReactionID",
                        [](const RVec<MC_V0> &v0s) { return Map(v0s, [](const MC_V0 &mc_v0) { return mc_v0.reaction_id; }); }, {"MC_" + v0_name});
    }

    return df_V0s;
}

/*                */
/**  Sexaquarks  **/
/*** ========== ***/

/** # Channel A **/

/*
 * Using Kalman Filter, reconstruct (anti)sexaquark candidates via the reaction channels
 * - `AntiSexaquark Neutron -> AntiLambda K0S`
 * - `X -> Lambda K0S`
 */
RVec<TypeA> Manager::TypeA_KF_Finder(const RVec<KF_V0> &found_v0a, const RVec<KF_V0> &found_v0b,  //
                                     const Float_t &magnetic_field, const KFVertex &kf_pv) {
    RVec<TypeA> output;
    TypeA this_sexa;
    KF_V0 cp_v0a, cp_v0b;
    /* protection */
    if (!found_v0a.size() || !found_v0b.size()) return output;
    /* one-time definitions */
    KFParticle::SetField(magnetic_field);
    const Double_t neutron_mass = TDatabasePDG::Instance()->GetParticle(PdgCode::Neutron)->Mass();
    /* loop over all pairs */
    for (const auto &v0a : found_v0a) {
        for (const auto &v0b : found_v0b) {
            /* sanity check */
            std::set<ULong64_t> unique_track_entries = {v0a.neg.entry, v0a.pos.entry, v0b.neg.entry, v0b.pos.entry};
            if (unique_track_entries.size() < 4) continue;
            cp_v0a = v0a;
            cp_v0b = v0b;
            /* fit */
            KFParticle kf_sexa(cp_v0a.kf, cp_v0b.kf);
            kf_sexa.SetProductionVertex(kf_pv);
            //
            kf_sexa.TransportToDecayVertex();
            cp_v0a.kf.SetProductionVertex(kf_sexa);
            cp_v0b.kf.SetProductionVertex(kf_sexa);
            //
            cp_v0a.kf.TransportToDecayVertex();
            cp_v0a.neg.kf.SetProductionVertex(cp_v0a.kf);
            cp_v0a.pos.kf.SetProductionVertex(cp_v0a.kf);
            //
            cp_v0b.kf.TransportToDecayVertex();
            cp_v0b.neg.kf.SetProductionVertex(cp_v0b.kf);
            cp_v0b.pos.kf.SetProductionVertex(cp_v0b.kf);
            /* transport tracks to V0s vertices */
            cp_v0a.neg.kf.TransportToProductionVertex();
            cp_v0a.pos.kf.TransportToProductionVertex();
            cp_v0b.neg.kf.TransportToProductionVertex();
            cp_v0b.pos.kf.TransportToProductionVertex();
            cp_v0a.neg.lv.SetCoordinates(cp_v0a.neg.kf.Px(), cp_v0a.neg.kf.Py(), cp_v0a.neg.kf.Pz(), cp_v0a.neg.lv.M());
            cp_v0a.pos.lv.SetCoordinates(cp_v0a.pos.kf.Px(), cp_v0a.pos.kf.Py(), cp_v0a.pos.kf.Pz(), cp_v0a.pos.lv.M());
            cp_v0b.neg.lv.SetCoordinates(cp_v0b.neg.kf.Px(), cp_v0b.neg.kf.Py(), cp_v0b.neg.kf.Pz(), cp_v0b.neg.lv.M());
            cp_v0b.pos.lv.SetCoordinates(cp_v0b.pos.kf.Px(), cp_v0b.pos.kf.Py(), cp_v0b.pos.kf.Pz(), cp_v0b.pos.lv.M());
            /* transport V0s to secondary vertex */
            cp_v0a.kf.TransportToProductionVertex();
            cp_v0b.kf.TransportToProductionVertex();
            cp_v0a.lv = PxPyPzEVector(cp_v0a.kf.Px(), cp_v0a.kf.Py(), cp_v0a.kf.Pz(), cp_v0a.kf.E());
            cp_v0b.lv = PxPyPzEVector(cp_v0b.kf.Px(), cp_v0b.kf.Py(), cp_v0b.kf.Pz(), cp_v0b.kf.E());
            /* fill struct */
            this_sexa.v0a = cp_v0a;
            this_sexa.v0b = cp_v0b;
            this_sexa.kf = kf_sexa;
            this_sexa.lv = PxPyPzEVector(cp_v0a.lv.Px() + cp_v0b.lv.Px(), cp_v0a.lv.Py() + cp_v0b.lv.Py(), cp_v0a.lv.Pz() + cp_v0b.lv.Pz(),
                                         cp_v0a.lv.E() + cp_v0b.lv.E() - neutron_mass);
            this_sexa.lv_asdecay = cp_v0a.lv + cp_v0b.lv;
            /* apply cuts and store */
            if (!TypeA_PassesCuts(this_sexa, kf_pv)) continue;
            output.push_back(this_sexa);
        }  // end of loop over k0s
    }      // end of loop over (anti)lambdas
    return output;
}

RVec<MC_TypeA> Manager::TypeA_TrueInfoCollector(const RVec<TypeA> &found,  //
                                                const RVec<MC_V0> &mc_v0a, const RVec<MC_V0> &mc_v0b) {
    RVec<MC_TypeA> output;
    output.reserve(found.size());
    MC_TypeA mc_sexa{};
    for (const auto &found_sexa : found) {
        mc_sexa.reaction_id = 0;
        mc_sexa.is_signal = false;
        /* fill values */
        mc_sexa.mc_v0a = mc_v0a[found_sexa.v0a.idx];
        mc_sexa.mc_v0b = mc_v0b[found_sexa.v0b.idx];
        if (mc_sexa.mc_v0a.reaction_id == mc_sexa.mc_v0b.reaction_id) {
            mc_sexa.reaction_id = mc_sexa.mc_v0a.reaction_id;
            mc_sexa.is_signal = mc_sexa.mc_v0a.is_signal && mc_sexa.mc_v0b.is_signal;
        }
        mc_sexa.is_hybrid = !mc_sexa.is_signal &&
                            ((mc_sexa.mc_v0a.is_signal && !mc_sexa.mc_v0b.is_signal) || (!mc_sexa.mc_v0a.is_signal && mc_sexa.mc_v0b.is_signal) ||  //
                             mc_sexa.mc_v0a.is_hybrid || mc_sexa.mc_v0b.is_hybrid);
        output.emplace_back(mc_sexa);
    }
    return output;
}

Bool_t Manager::TypeA_PassesCuts(const TypeA &sexa, const KFVertex &kf_pv) {
    //
    XYZPoint v3_sexa(sexa.kf.GetX(), sexa.kf.GetY(), sexa.kf.GetZ());
    XYZPoint v3_pv(kf_pv.GetX(), kf_pv.GetY(), kf_pv.GetZ());
    /*  */
    if (v3_sexa.Rho() < Default::ChannelA::MinRadius) return kFALSE;
    if (TMath::Abs(sexa.lv.Rapidity()) > Default::ChannelA::AbsMaxRapidity) return kFALSE;
    /*  */ Double_t cpa_wrt_pv = Math::CosinePointingAngle(sexa.lv.Vect(), v3_sexa, v3_pv);
    if (cpa_wrt_pv < Default::ChannelA::MinCPAwrtPV || cpa_wrt_pv > Default::ChannelA::MaxCPAwrtPV) return kFALSE;
    /*  */ Double_t dca_la_sv = TMath::Abs((Double_t)sexa.v0a.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_sv > Default::ChannelA::MaxDCALaSV) return kFALSE;
    /*  */ Double_t dca_la_neg_sv = TMath::Abs((Double_t)sexa.v0a.neg.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_neg_sv > Default::ChannelA::MaxDCALaNegSV) return kFALSE;
    /*  */ Double_t dca_la_pos_sv = TMath::Abs((Double_t)sexa.v0a.pos.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_pos_sv > Default::ChannelA::MaxDCALaPosSV) return kFALSE;
    /*  */ Double_t dca_k0s_sv = TMath::Abs((Double_t)sexa.v0b.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_k0s_sv > Default::ChannelA::MaxDCAK0SV) return kFALSE;
    /*  */ Double_t dca_k0s_neg_sv = TMath::Abs((Double_t)sexa.v0b.neg.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_k0s_neg_sv > Default::ChannelA::MaxDCAK0NegSV) return kFALSE;
    /*  */ Double_t dca_k0s_pos_sv = TMath::Abs((Double_t)sexa.v0b.pos.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_k0s_pos_sv > Default::ChannelA::MaxDCAK0PosSV) return kFALSE;
    /*  */ Double_t dca_btw_v0s = TMath::Abs((Double_t)sexa.v0a.kf.GetDistanceFromParticle(sexa.v0b.kf));
    if (dca_btw_v0s > Default::ChannelA::MaxDCAbtwV0s) return kFALSE;
    return kTRUE;
}

RNode Manager::FindSexaquarks_TypeA(RNode df, PdgCode pdg_struck_nucleon, const std::vector<PdgCode> &pdg_reaction_products) {
    //
    std::string sexaquark_name;
    if (pdg_struck_nucleon > 0)
        sexaquark_name = "ASA";  // antisexaquark type a
    else
        sexaquark_name = "BA";  // background type a
    std::string v0a_name = fParticleName_[pdg_reaction_products[0]];
    std::string v0b_name = fParticleName_[pdg_reaction_products[1]];
    fAnalyzed_Channels.push_back(sexaquark_name);

    RNode df_Sexaquarks = df.Define("Found_" + sexaquark_name, TypeA_KF_Finder,  //
                                    {"Found_" + v0a_name, "Found_" + v0b_name, "MagneticField", "Event_KF_PV"})
                              .Define(sexaquark_name + "_V0a_Idx",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.v0a.idx; }); },
                                      {"Found_" + sexaquark_name})
                              .Define(sexaquark_name + "_V0b_Idx",
                                      [](const RVec<TypeA> &sexaquarks) { return Map(sexaquarks, [](const TypeA &sexa) { return sexa.v0b.idx; }); },
                                      {"Found_" + sexaquark_name})
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
            df_Sexaquarks.Define("MC_" + sexaquark_name, TypeA_TrueInfoCollector, {"Found_" + sexaquark_name, "MC_" + v0a_name, "MC_" + v0b_name})
                .Define(sexaquark_name + "_V0a_Neg_McEntry",
                        [](const RVec<MC_TypeA> &sexaquarks) {
                            return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0a.neg.entry; });
                        },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_V0a_Pos_McEntry",
                        [](const RVec<MC_TypeA> &sexaquarks) {
                            return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0a.pos.entry; });
                        },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_V0b_Neg_McEntry",
                        [](const RVec<MC_TypeA> &sexaquarks) {
                            return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0b.neg.entry; });
                        },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_V0b_Pos_McEntry",
                        [](const RVec<MC_TypeA> &sexaquarks) {
                            return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0b.pos.entry; });
                        },
                        {"MC_" + sexaquark_name})
                .Define(
                    sexaquark_name + "_V0a_McEntry",
                    [](const RVec<MC_TypeA> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0a.mc_entry; }); },
                    {"MC_" + sexaquark_name})
                .Define(
                    sexaquark_name + "_V0b_McEntry",
                    [](const RVec<MC_TypeA> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.mc_v0b.mc_entry; }); },
                    {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_IsSignal",
                        [](const RVec<MC_TypeA> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.is_signal; }); },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_ReactionID",
                        [](const RVec<MC_TypeA> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.reaction_id; }); },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_IsHybrid",
                        [](const RVec<MC_TypeA> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeA &mc_sexa) { return mc_sexa.is_hybrid; }); },
                        {"MC_" + sexaquark_name});
    }

    return df_Sexaquarks;
}

/** # Channel D **/

/*
 * Using Kalman Filter, reconstruct (anti)sexaquark candidates via the reaction channels
 * - `AntiSexaquark Proton -> AntiLambda K+`
 * - `X -> Lambda K-`
 */
RVec<TypeD> Manager::TypeD_KF_Finder(const RVec<KF_V0> &found_v0s, const RVec<KF_Track> &bach_tracks,  //
                                     const Float_t &magnetic_field, const KFVertex &kf_pv) {
    RVec<TypeD> output;
    TypeD this_sexa;
    KF_V0 cp_v0;
    KF_Track cp_ba;
    /* protection */
    if (!found_v0s.size() || !bach_tracks.size()) return output;
    /* one-time definitions */
    KFParticle::SetField(magnetic_field);
    const Double_t proton_mass = TDatabasePDG::Instance()->GetParticle(PdgCode::Proton)->Mass();
    /* loop over all pairs */
    for (const auto &v0 : found_v0s) {
        for (const auto &ba : bach_tracks) {
            /* sanity check */
            std::set<ULong64_t> unique_track_entries = {v0.neg.entry, v0.pos.entry, ba.entry};
            if (unique_track_entries.size() < 3) continue;
            cp_v0 = v0;
            cp_ba = ba;
            /* fit */
            KFParticle kf_sexa(cp_v0.kf, cp_ba.kf);
            kf_sexa.SetProductionVertex(kf_pv);
            /* assign secondary vertex as production vertex for V0 and bachelor */
            kf_sexa.TransportToDecayVertex();
            cp_v0.kf.SetProductionVertex(kf_sexa);
            cp_ba.kf.SetProductionVertex(kf_sexa);
            /* assign V0 as production vertex for V0 daughters */
            cp_v0.kf.TransportToDecayVertex();
            cp_v0.neg.kf.SetProductionVertex(cp_v0.kf);
            cp_v0.pos.kf.SetProductionVertex(cp_v0.kf);
            /* transport daughters to V0 vertex */
            cp_v0.neg.kf.TransportToProductionVertex();
            cp_v0.pos.kf.TransportToProductionVertex();
            cp_v0.neg.lv.SetCoordinates(cp_v0.neg.kf.Px(), cp_v0.neg.kf.Py(), cp_v0.neg.kf.Pz(), cp_v0.neg.lv.M());
            cp_v0.pos.lv.SetCoordinates(cp_v0.pos.kf.Px(), cp_v0.pos.kf.Py(), cp_v0.pos.kf.Pz(), cp_v0.pos.lv.M());
            /* transport V0 and bachelor to secondary vertex */
            cp_v0.kf.TransportToProductionVertex();
            cp_ba.kf.TransportToProductionVertex();
            cp_v0.lv = PxPyPzEVector(cp_v0.kf.Px(), cp_v0.kf.Py(), cp_v0.kf.Pz(), cp_v0.kf.E());
            cp_ba.lv.SetCoordinates(cp_ba.kf.Px(), cp_ba.kf.Py(), cp_ba.kf.Pz(), cp_ba.lv.M());
            /* fill struct */
            this_sexa.v0 = cp_v0;
            this_sexa.ba = cp_ba;
            this_sexa.kf = kf_sexa;
            this_sexa.lv = PxPyPzEVector(cp_v0.lv.Px() + cp_ba.lv.Px(), cp_v0.lv.Py() + cp_ba.lv.Py(), cp_v0.lv.Pz() + cp_ba.lv.Pz(),
                                         cp_v0.lv.E() + cp_ba.lv.E() - proton_mass);
            this_sexa.lv_asdecay = cp_v0.lv + cp_ba.lv;
            /* apply cuts and store */
            if (!TypeD_PassesCuts(this_sexa, kf_pv)) continue;
            output.push_back(this_sexa);
        }  // end of loop over charged kaons
    }      // end of loop over (anti)lambdas
    return output;
}

RVec<MC_TypeD> Manager::TypeD_TrueInfoCollector(PdgCode pdg_code_ba, const RVec<TypeD> &found,  //
                                                const RVec<MC_V0> &mc_v0, const RVec<MC_Track> &mc_ba) {
    RVec<MC_TypeD> output;
    output.reserve(found.size());
    MC_TypeD mc_sexa{};
    for (const auto &found_sexa : found) {
        mc_sexa.reaction_id = 0;
        mc_sexa.is_signal = false;
        /* fill values */
        mc_sexa.mc_v0 = mc_v0[found_sexa.v0.idx];
        mc_sexa.mc_ba = mc_ba[found_sexa.ba.entry];
        if (mc_sexa.mc_v0.reaction_id == mc_sexa.mc_ba.reaction_id) {
            mc_sexa.reaction_id = mc_sexa.mc_v0.reaction_id;
            mc_sexa.is_signal = mc_sexa.mc_v0.is_signal && mc_sexa.mc_ba.is_signal && mc_sexa.mc_ba.pdg_code == pdg_code_ba;
        }
        mc_sexa.is_hybrid = !mc_sexa.is_signal &&  //
                            ((mc_sexa.mc_v0.is_signal && !mc_sexa.mc_ba.is_signal) || (!mc_sexa.mc_v0.is_signal && mc_sexa.mc_ba.is_signal) ||
                             mc_sexa.mc_v0.is_hybrid);
        output.emplace_back(mc_sexa);
    }
    return output;
}

Bool_t Manager::TypeD_PassesCuts(const TypeD &sexa, const KFVertex &kf_pv) {
    //
    XYZPoint v3_sexa(sexa.kf.GetX(), sexa.kf.GetY(), sexa.kf.GetZ());
    XYZPoint v3_pv(kf_pv.GetX(), kf_pv.GetY(), kf_pv.GetZ());
    /*  */
    if (v3_sexa.Rho() < Default::ChannelD::MinRadius) return kFALSE;
    if (TMath::Abs(sexa.lv.Rapidity()) > Default::ChannelD::AbsMaxRapidity) return kFALSE;
    /*  */ Double_t cpa_wrt_pv = Math::CosinePointingAngle(sexa.lv.Vect(), v3_sexa, v3_pv);
    if (cpa_wrt_pv < Default::ChannelD::MinCPAwrtPV || cpa_wrt_pv > Default::ChannelD::MaxCPAwrtPV) return kFALSE;
    /*  */ Double_t dca_la_sv = TMath::Abs((Double_t)sexa.v0.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_sv > Default::ChannelD::MaxDCALaSV) return kFALSE;
    /*  */ Double_t dca_la_neg_sv = TMath::Abs((Double_t)sexa.v0.neg.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_neg_sv > Default::ChannelD::MaxDCALaNegSV) return kFALSE;
    /*  */ Double_t dca_la_pos_sv = TMath::Abs((Double_t)sexa.v0.pos.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_la_pos_sv > Default::ChannelD::MaxDCALaPosSV) return kFALSE;
    /*  */ Double_t dca_ka_sv = TMath::Abs((Double_t)sexa.ba.kf.GetDistanceFromVertex(sexa.kf));
    if (dca_ka_sv > Default::ChannelD::MaxDCAKaSV) return kFALSE;
    /*  */ Double_t dca_ka_la = TMath::Abs((Double_t)sexa.ba.kf.GetDistanceFromVertex(sexa.v0.kf));
    if (dca_ka_la > Default::ChannelD::MaxDCAKaLa) return kFALSE;
    return kTRUE;
}

RNode Manager::FindSexaquarks_TypeD(RNode df, PdgCode pdg_struck_nucleon, const std::vector<PdgCode> &pdg_reaction_products) {
    //
    std::string sexaquark_name;
    if (pdg_struck_nucleon > 0)
        sexaquark_name = "ASD";  // antisexaquark type d
    else
        sexaquark_name = "BD";  // background type d
    fAnalyzed_Channels.push_back(sexaquark_name);
    PdgCode pdg_code_v0 = pdg_reaction_products[0];
    PdgCode pdg_code_ba = pdg_reaction_products[1];
    std::string v0_name = fParticleName_[pdg_code_v0];
    std::string bach_name = fParticleName_[pdg_code_ba];
    //
    auto CallTrueCollector = [pdg_code_ba](const RVec<TypeD> &found,  //
                                           const RVec<MC_V0> &mc_v0, const RVec<MC_Track> &mc_ba) -> RVec<MC_TypeD> {
        return TypeD_TrueInfoCollector(pdg_code_ba, found, mc_v0, mc_ba);
    };

    RNode df_Sexaquarks =
        df.Define("Found_" + sexaquark_name, TypeD_KF_Finder,  //
                  {"Found_" + v0_name, "KF_" + bach_name, "MagneticField", "Event_KF_PV"})
            .Define(sexaquark_name + "_V0_Idx",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.v0.idx; }); },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_Ba_Entry",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.ba.entry; }); },
                    {"Found_" + sexaquark_name})
#if READ_ESD_INDICES
            .Define(sexaquark_name + "_V0_Neg_EsdIdx",
                    [](const RVec<TypeD> &sexaquarks, cRVecL esd_idx) {
                        return Map(sexaquarks, [&esd_idx](const TypeD &sexa) { return esd_idx[sexa.v0.neg.entry]; });
                    },
                    {"Found_" + sexaquark_name, "Track_EsdIdx"})
            .Define(sexaquark_name + "_V0_Pos_EsdIdx",
                    [](const RVec<TypeD> &sexaquarks, cRVecL esd_idx) {
                        return Map(sexaquarks, [&esd_idx](const TypeD &sexa) { return esd_idx[sexa.v0.pos.entry]; });
                    },
                    {"Found_" + sexaquark_name, "Track_EsdIdx"})
            .Define(sexaquark_name + "_Ba_EsdIdx",
                    [](const RVec<TypeD> &sexaquarks, cRVecL esd_idx) {
                        return Map(sexaquarks, [&esd_idx](const TypeD &sexa) { return esd_idx[sexa.ba.entry]; });
                    },
                    {"Found_" + sexaquark_name, "Track_EsdIdx"})
#endif
            .Define(sexaquark_name + "_Xv",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.kf.GetX(); }); },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_Yv",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.kf.GetY(); }); },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_Mass",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.lv.M(); }); },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_Pt",
                    [](const RVec<TypeD> &sexaquarks) { return Map(sexaquarks, [](const TypeD &sexa) { return sexa.lv.Pt(); }); },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_Radius",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks,
                                   [](const TypeD &sexa) { return TMath::Sqrt(sexa.kf.GetX() * sexa.kf.GetX() + sexa.kf.GetY() * sexa.kf.GetY()); });
                    },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_CPAwrtPV",
                    [](const RVec<TypeD> &sexaquarks, const XYZPoint &pv) {
                        return Map(sexaquarks, [&pv](const TypeD &sexa) {
                            return Math::CosinePointingAngle(sexa.lv.Vect(), XYZPoint(sexa.kf.GetX(), sexa.kf.GetY(), sexa.kf.GetZ()), pv);
                        });
                    },
                    {"Found_" + sexaquark_name, "Event_PV"})
            .Define(sexaquark_name + "_DCALaSV",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks, [](const TypeD &sexa) { return TMath::Abs((Double_t)sexa.v0.kf.GetDistanceFromVertex(sexa.kf)); });
                    },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_DCALaNegSV",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks, [](const TypeD &sexa) { return TMath::Abs((Double_t)sexa.v0.neg.kf.GetDistanceFromVertex(sexa.kf)); });
                    },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_DCALaPosSV",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks, [](const TypeD &sexa) { return TMath::Abs((Double_t)sexa.v0.pos.kf.GetDistanceFromVertex(sexa.kf)); });
                    },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_DCAKaSV",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks, [](const TypeD &sexa) { return TMath::Abs((Double_t)sexa.ba.kf.GetDistanceFromVertex(sexa.kf)); });
                    },
                    {"Found_" + sexaquark_name})
            .Define(sexaquark_name + "_DCAKaLa",
                    [](const RVec<TypeD> &sexaquarks) {
                        return Map(sexaquarks, [](const TypeD &sexa) { return TMath::Abs((Double_t)sexa.ba.kf.GetDistanceFromVertex(sexa.v0.kf)); });
                    },
                    {"Found_" + sexaquark_name});

    if (Settings::IsMC) {
        df_Sexaquarks =
            df_Sexaquarks
                .Define("MC_" + sexaquark_name, CallTrueCollector,  //
                        {"Found_" + sexaquark_name, "MC_" + v0_name, "MC_Tracks"})
                .Define(
                    sexaquark_name + "_V0_Neg_McEntry",
                    [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.mc_v0.neg.entry; }); },
                    {"MC_" + sexaquark_name})
                .Define(
                    sexaquark_name + "_V0_Pos_McEntry",
                    [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.mc_v0.pos.entry; }); },
                    {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_Ba_McEntry",
                        [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.mc_ba.entry; }); },
                        {"MC_" + sexaquark_name})
                .Define(
                    sexaquark_name + "_V0_McEntry",
                    [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.mc_v0.mc_entry; }); },
                    {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_IsSignal",
                        [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.is_signal; }); },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_ReactionID",
                        [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.reaction_id; }); },
                        {"MC_" + sexaquark_name})
                .Define(sexaquark_name + "_IsHybrid",
                        [](const RVec<MC_TypeD> &sexaquarks) { return Map(sexaquarks, [](const MC_TypeD &mc_sexa) { return mc_sexa.is_hybrid; }); },
                        {"MC_" + sexaquark_name});
    }

    return df_Sexaquarks;
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
    std::cout << "NRuns = " << df.GetNRuns() << '\n';
}

void Manager::EndOfAnalysis(RNode df) {
    //
    std::vector<std::string> column_list = {"EventNumber"};
    /* Add V0s Properties */
    for (const auto &v0_name : fAnalyzed_V0sNames) {
        column_list.insert(column_list.end(), {
            v0_name + "_Neg_TrackEntry", v0_name + "_Pos_TrackEntry",  //
#if READ_ESD_INDICES
                v0_name + "_Neg_EsdIdx", v0_name + "_Pos_EsdIdx",  //
#endif
                v0_name + "_Px", v0_name + "_Py", v0_name + "_Pz",              //
                v0_name + "_Neg_Px", v0_name + "_Neg_Py", v0_name + "_Neg_Pz",  //
                v0_name + "_Pos_Px", v0_name + "_Pos_Py", v0_name + "_Pos_Pz",  //
                v0_name + "_Xv", v0_name + "_Yv", v0_name + "_Zv",              //
                v0_name + "_Mass", v0_name + "_Radius"
        });
        if (Settings::IsMC) {
            column_list.insert(column_list.end(),
                               {v0_name + "_Neg_McEntry", v0_name + "_Pos_McEntry", v0_name + "_Neg_PdgCode", v0_name + "_Pos_PdgCode",  //
                                v0_name + "_McEntry", v0_name + "_PdgCode", v0_name + "_IsTrue", v0_name + "_IsSignal", v0_name + "_IsHybrid",
                                v0_name + "_ReactionID"});
        }
    }
    /* Add (Anti)Sexaquarks Properties */
    for (const auto &channel_name : fAnalyzed_Channels) {
        /* -- Channel A */
        if (channel_name.back() == 'A') {
            column_list.insert(column_list.end(), {channel_name + "_V0a_Idx", channel_name + "_V0b_Idx",  //
                                                   channel_name + "_Xv", channel_name + "_Yv",            //
                                                   channel_name + "_Mass", channel_name + "_Pt", channel_name + "_Radius"});
            if (Settings::IsMC) {
                column_list.insert(column_list.end(),
                                   {channel_name + "_V0a_Neg_McEntry", channel_name + "_V0a_Pos_McEntry", channel_name + "_V0b_Neg_McEntry",
                                    channel_name + "_V0b_Pos_McEntry", channel_name + "_V0a_McEntry", channel_name + "_V0b_McEntry",
                                    channel_name + "_IsSignal", channel_name + "_ReactionID", channel_name + "_IsHybrid"});
            }
        }
        /* -- Channel D */
        if (channel_name.back() == 'D') {
            column_list.insert(column_list.end(), {
                channel_name + "_V0_Idx", channel_name + "_Ba_Entry",
#if READ_ESD_INDICES
                    channel_name + "_V0_Neg_EsdIdx", channel_name + "_V0_Pos_EsdIdx", channel_name + "_Ba_EsdIdx",  //
#endif
                    channel_name + "_Xv", channel_name + "_Yv",                              //
                    channel_name + "_Mass", channel_name + "_Pt", channel_name + "_Radius",  //
                    channel_name + "_CPAwrtPV", channel_name + "_DCALaSV",                   //
                    channel_name + "_DCALaNegSV", channel_name + "_DCALaPosSV",              //
                    channel_name + "_DCAKaSV", channel_name + "_DCAKaLa"
            });
            if (Settings::IsMC) {
                column_list.insert(column_list.end(), {channel_name + "_V0_Neg_McEntry", channel_name + "_V0_Pos_McEntry",
                                                       channel_name + "_Ba_McEntry", channel_name + "_V0_McEntry", channel_name + "_IsSignal",
                                                       channel_name + "_ReactionID", channel_name + "_IsHybrid"});
            }
        }
    }
    /* Write to disk */
    df.Snapshot("Events", Settings::PathOutputFile, column_list);
    std::cout << "TFile " << Settings::PathOutputFile << " has been written" << '\n';
    /* Print info */
    std::cout << "NRuns = " << df.GetNRuns() << '\n';
}

}  // namespace Tree2Sexaquark::Analysis
