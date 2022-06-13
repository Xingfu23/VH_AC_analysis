#ifndef __PRESELECTION_H__
#define __PRESELECTION_H__

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiE4D.h"

using namespace std;
using Position4_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>;


bool Event_Preselection(Position4_t a1, Position4_t a2, float met_Pt, float met_Phi) {
    // For background study, here we need to erase the signal region (diphoton mass between 120 and 130 GeV)
    // And we ask missing energy's pT needs bigger than 50 GeV
    // The leading photon pT bigger than 30 GeV, and the subleading one need to larger than 25GeV
    // Ratio between leading photon and diphoton mass is bigger than 1/3 and the same ratio but for subleading photon is required bigger than 1/5
    // the phi difference between diphoton should over 2 and adding the eta restriction
    int indicator = 0;
    Position4_t diphoton = a1 + a2;

    if ((a1.Pt() / diphoton.M() < 0.333) || (a2.Pt() / diphoton.M() < 0.2)) indicator += 1;
    if (diphoton.M() > 120. && diphoton.M() < 130.) indicator += 1;
    if (diphoton.M() == 0.) indicator += 1;
    if (a1.Pt() < 30. || a2.Pt() < 25. || met_Pt < 15.) indicator += 1;
    if (fabs(a1.Eta() > 2.5) || (fabs(a1.Eta()) < 1.57 && fabs(a1.Eta()) > 1.44)) indicator += 1;
    if (fabs(a2.Eta() > 2.5) || (fabs(a2.Eta()) < 1.57 && fabs(a2.Eta()) > 1.44)) indicator += 1;
    // if (fabs(diphoton.Phi() - met_Phi) < 2.) indicator += 1;
    if (fabs(DeltaPhi(diphoton.Phi(), met_Phi)) < 2.) indicator += 1;
    if (met_Pt < 50.) indicator += 1;

    if (indicator == 0) {
        return true;
    } else {
        return false;
    }
};

bool Event_Preselection_sig(Position4_t a1, Position4_t a2, float met_Pt, float met_Phi) {
    int indicator = 0;
    Position4_t diphoton = a1 + a2;

    if ((a1.Pt() / diphoton.M() < 0.333) || (a2.Pt() / diphoton.M() < 0.2)) indicator += 1;
    if (diphoton.M() < 120. || diphoton.M() > 130.) indicator += 1;
    if (diphoton.M() == 0.) indicator += 1;
    if (a1.Pt() < 30. || a2.Pt() < 25. || met_Pt < 15.) indicator += 1;
    if (fabs(a1.Eta() > 2.5) || (fabs(a1.Eta()) < 1.57 && fabs(a1.Eta()) > 1.44)) indicator += 1;
    if (fabs(a2.Eta() > 2.5) || (fabs(a2.Eta()) < 1.57 && fabs(a2.Eta()) > 1.44)) indicator += 1;
    if (fabs(diphoton.Phi() - met_Phi) < 2.) indicator += 1;
    if (met_Pt < 50.) indicator += 1;

    if (indicator == 0) {
        return true;
    } else {
        return false;
    }
};

#endif  // __PRESELECTION_H__