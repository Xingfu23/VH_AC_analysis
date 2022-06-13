#ifndef __PRESELECTION_K_H__
#define __PRESELECTION_K_H__

#endif  // __PRESELECTION_K_H__

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiE4D.h"

// using namespace std;
using Position4_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>;

float DeltaPhi(float phi_a, float phi_b) {
    float dphi = phi_a - phi_b;
    if (dphi > M_PI) {
        dphi -= 2. * M_PI;
    } else if (dphi < -M_PI) {
        dphi += 2. * M_PI;
    }
    return dphi;
};

bool Event_Preselection_sig(Position4_t a1, Position4_t a2, float met_Pt, float met_Phi) {
    int indicator = 0;
    Position4_t diphoton = a1 + a2;
    if ((a1.Pt() / diphoton.M() < 0.333) || (a2.Pt() / diphoton.M() < 0.2)) indicator += 1;
    if (diphoton.M() < 120. || diphoton.M() > 130.) indicator += 1;
    if (a1.Pt() < 30. || a2.Pt() < 25. || met_Pt < 50.) indicator += 1;
    if (fabs(a1.Eta() > 2.5) || (fabs(a1.Eta()) < 1.57 && fabs(a1.Eta()) > 1.44)) indicator += 1;
    if (fabs(a2.Eta() > 2.5) || (fabs(a2.Eta()) < 1.57 && fabs(a2.Eta()) > 1.44)) indicator += 1;
    if (fabs(DeltaPhi(diphoton.Phi(), met_Phi)) < 2.) indicator += 1;

    bool result = (indicator == 0) ? true : false;
    return result;
}

bool Event_Preselection_bkg(Position4_t a1, Position4_t a2, float met_Pt, float met_Phi) {
    int indicator = 0;
    Position4_t diphoton = a1 + a2;
    if ((a1.Pt() / diphoton.M() < 0.333) || (a2.Pt() / diphoton.M() < 0.2)) indicator += 1;
    if (diphoton.M() > 120. && diphoton.M() < 130.) indicator += 1;
    if (a1.Pt() < 30. || a2.Pt() < 25. || met_Pt < 50.) indicator += 1;
    if (fabs(a1.Eta() > 2.5) || (fabs(a1.Eta()) < 1.57 && fabs(a1.Eta()) > 1.44)) indicator += 1;
    if (fabs(a2.Eta() > 2.5) || (fabs(a2.Eta()) < 1.57 && fabs(a2.Eta()) > 1.44)) indicator += 1;
    if (fabs(DeltaPhi(diphoton.Phi(), met_Phi)) < 2.) indicator += 1;

    bool result = (indicator == 0) ? true : false;
    return result;
}

bool Event_Preselection_jet(Position4_t jet){
    int indicator = 0;
    if (jet.Pt() < 15. || fabs(jet.Eta() > 2.0)) indicator += 1;

    bool result = (indicator == 0) ? true : false;
    return result;
}