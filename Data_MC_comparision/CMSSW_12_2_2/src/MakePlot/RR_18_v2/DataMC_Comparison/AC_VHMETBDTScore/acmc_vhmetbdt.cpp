#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObject.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "tools/drawingoptions.h"
#include "tools/preselection_k.h"

using Position4_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>;
using std::cin;
using std::cout;
using std::string;
using std::unique_ptr;
using std::vector;
using namespace TMVA;

float RR18_Lumi = 56.69;
float scale_factor_GJets = 0.97;
float scale_factor_DiPhoton = 2.39;

// common variables
float weight;
float dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE, dipho_leadIDMVA;
float dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE, dipho_subleadIDMVA;
float met_Pt, met_Phi, met_sumEt;
float n_jets, jet_Pt[15], jet_Eta[15], jet_Phi[15], jet_E[15], jet_deepbtag[15];

// *file path and names of SM background and histograms' colours.
// data
string AbsRoute = "/wk_cms2/iceeric02/Workspace/VH_AC_analysis/Data_MC_comparision/CMSSW_12_2_2/src/DownloadSample/RR_18_v2";
vector<string> filetype_data = {"/Output_file_RR_18_Data/output_EGamma_spigazzi-Era2018_RR-17Sep2018_v2-legacyRun2FullV2-v0-Run2018.root"};

// *Input AC, sig and bkg file route
// SM signal and AC signal
// vector<string> filetype_list_zh = {"/Output_file_RR_18_MCSet1/output_ZHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"};
// vector<string> filetype_list_wh = {"/Output_file_RR_18_MCSet1/output_WHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"};
vector<string> filetype_list_zh = {"/output_ZHiggs0PMToGG_M125_13TeV_JHUGenV7011_pythia8.root"};
vector<string> filetype_list_wh = {"/output_WHiggs0PMToGG_M125_13TeV_JHUGenV7011_pythia8.root"};
vector<string> filetype_list_zhJHUGenV7 = {"/output_ZHiggs0PHToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_ZHiggs0MToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_ZHiggs0L1ToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_ZHiggs0PMToGG_M125_13TeV_JHUGenV7011_pythia8.root"};
vector<string> filetype_list_whJHUGenV7 = {"/output_WHiggs0PHToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_WHiggs0MToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_WHiggs0L1ToGG_M125_13TeV_JHUGenV7011_pythia8.root",
                                           "/output_WHiggs0PMToGG_M125_13TeV_JHUGenV7011_pythia8.root"};
vector<string> mc_routename_sig = {"tagsDumper/trees/zh_125_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/wh_125_13TeV_RECO_VH_MET_Tag0"};
vector<string> mc_routename_JHUGenV7 = {"tagsDumper/trees/zh_125GenV7_13TeV_RECO_VH_MET_Tag0",
                                        "tagsDumper/trees/wh_125GenV7_13TeV_RECO_VH_MET_Tag0"};

// SM background
vector<string> filetype_list_Diboson = {"/Output_file_RR_18_MCSet3/output_WW_TuneCP5_13TeV-pythia8.root",
                                        "/Output_file_RR_18_MCSet3/output_WZ_TuneCP5_13TeV-pythia8.root",
                                        "/Output_file_RR_18_MCSet3/output_ZZ_TuneCP5_13TeV-pythia8.root",
                                        "/Output_file_RR_18_MCSet3/output_WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.root",
                                        "/Output_file_RR_18_MCSet3/output_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.root"};
vector<string> filetype_list_Diphoton = {"/Output_file_RR_18_MCSet3/output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root"};
vector<string> filetype_list_ggh = {"/Output_file_RR_18_MCSet2/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root",
                                    "/Output_file_RR_18_MCSet2/output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root"};
vector<string> filetype_list_tth = {"/Output_file_RR_18_MCSet2/output_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root"};
vector<string> filetype_list_vbf = {"/Output_file_RR_18_MCSet2/output_VBFHToGG_M-125_TuneCP5_13TeV_powheg_pythia8_withDipoleRecoil.root",
                                    "/Output_file_RR_18_MCSet2/output_VBFHToGG_M125_13TeV_amcatnlo_pythia8.root"};
vector<string> filetype_list_DY = {"/Output_file_RR_18_MCSet3/output_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root"};
vector<string> filetype_list_Topbkg = {"/Output_file_RR_18_MCSet3/output_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
                                       "/Output_file_RR_18_MCSet3/output_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root"};
vector<string> filetype_list_gplusjets = {"/Output_file_RR_18_MCSet2/output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
                                          "/Output_file_RR_18_MCSet2/output_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8.root",
                                          "/Output_file_RR_18_MCSet2/output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root"};

vector<string> mc_routename_bkg = {"tagsDumper/trees/diboson_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/diphoton_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/ggh_125_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/tth_125_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/vbf_125_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/DYjets_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/Top_bkg_13TeV_RECO_VH_MET_Tag0",
                                   "tagsDumper/trees/gjet_fakephoton_13TeV_RECO_VH_MET_Tag0"};

// *Input weight file
string Weightfile_route = "/wk_cms2/iceeric02/Workspace/VH_AC_analysis/VH_BDT_train/VH_AC_BDTClassifier/BDTClassifier/Output_xml";
vector<string> Weightfiletypezh_list = {"/output_ZHiggs0PHToGG_M125.xml", "/output_ZHiggs0MToGG_M125.xml", "/output_ZHiggs0L1ToGG_M125.xml"};
vector<string> Weightfiletypewh_list = {"/output_WHiggs0PHToGG_M125.xml", "/output_WHiggs0MToGG_M125.xml", "/output_WHiggs0L1ToGG_M125.xml"};

// customise functions
void SettingBranch_Photon(TTree *tree);
void SettingBranch_MET(TTree *tree);
void SettingBranch_Jet(TTree *tree);
float ConvertProbability(float bdtscore);

vector<string> file_naming_zh = {"ZHiggs0PHToGG_M125_VHBDTScore", "ZHiggs0MToGG_M125_VHBDTScore", "ZHiggs0L1ToGG_M125_VHBDTScore"};
vector<string> file_naming_wh = {"WHiggs0PHToGG_M125_VHBDTScore", "WHiggs0MToGG_M125_VHBDTScore", "WHiggs0L1ToGG_M125_VHBDTScore"};
// vector<Color_t> ac_histo_colour = {kRed - 7, kGreen + 2, kAzure + 2};
vector<string> ac_histo_ledgendname = {"(f_{a2}=1.0) #times 500", "(f_{a3}=1.0) #times 500", "(f_{#Lambda1}=1.0) #times 500"};

vector<string> filetype_legname = {"Dibosn", "#gamma#gamma", "ggh", "tth", "vbf", "Drell-Yen", "Top", "#gamma + Jets"};
vector<Color_t> histo_colours = {kOrange + 1, kAzure, kPink, kPink - 4, kPink + 6, kViolet + 8, kCyan + 1, kAzure - 4};

int main() {
    // sig_type == 0 means zh part, compare to sig_type == 1, which means wh part
    for (int sig_type = 0; sig_type < 2; sig_type++) {
        // Using AC sample as loop iterator, one AC sample output one graph
        // 0 means fa2 = 1., 1 means fa3 = 1. and 2 means fLambda1 = 1.
        for (int AC_type = 0; AC_type < 3; AC_type++) {
            //* Booking Canvas for vhbdt score
            auto c_bdtscore = new TCanvas("c_bdtscore", " ", 600, 650);
            auto c_probability = new TCanvas("c_probability", " ", 600, 650);
            auto leg_bdtscore = new TLegend(0.63, 0.6, 0.88, 0.88);
            auto leg_probability = new TLegend(0.63, 0.6, 0.88, 0.88);
            c_bdtscore->SetLogy();
            c_probability->SetLogy();

            // *Booking stack histogram for continuous background
            THStack *bkg_bdtscore_hs = new THStack("bkg_bdtscore_hs", " ");
            THStack *bkg_probability_hs = new THStack("bkg_probability_hs", " ");

            // *********************************************************************************************************** //
            // * Making signal Plot, and using the plots as main plots                                                     //
            // *********************************************************************************************************** //
            // Booking Histogram
            TH1F *sig_bdtscore_histo = new TH1F("sig_bdtscore_histo", " ", 50, -1., 1.);
            TH1F *sig_probability_histo = new TH1F("sig_probability_histo", " ", 50, 0., 1.);

            // Import file
            string filename_main, weight_file;
            if (sig_type == 0) {
                filename_main = AbsRoute + filetype_list_zh[0];
                weight_file = Weightfile_route + Weightfiletypezh_list[AC_type];
            } else {
                filename_main = AbsRoute + filetype_list_wh[0];
                weight_file = Weightfile_route + Weightfiletypewh_list[AC_type];
            }
            auto file_main = std::unique_ptr<TFile>{new TFile(filename_main.c_str())};
            cout << "Input File: " << filename_main << "\n";
            cout << "Input Weight File: " << weight_file << "\n";
            TTree *tree_main;
            // file_main->GetObject(mc_routename_sig[sig_type].c_str(), tree_main);
            file_main->GetObject(mc_routename_JHUGenV7[sig_type].c_str(), tree_main);

            SettingBranch_Photon(tree_main);
            SettingBranch_MET(tree_main);
            SettingBranch_Jet(tree_main);

            // Setup TMVA Reader
            std::unique_ptr<TMVA::Reader> reader;

            float _pho1_eta, _pho2_eta;
            float _pho1_ptoM, _pho2_ptoM;
            float _min_phoId, _max_phoId;
            float _dipho_cosphi;
            float _met, _met_sumEt;
            float _dphi_dipho_met;
            float _pt_balance;
            float _njet;
            float _max_jet_pt;
            float _max_jet_dCSV;
            float _min_dphi_jet_met;

            _pho1_eta = -999.;
            _pho2_eta = -999.;
            _pho1_ptoM = -999.;
            _pho2_ptoM = -999.;
            _min_phoId = -999.;
            _max_phoId = -999.;
            _dipho_cosphi = -999.;
            _met = -999.;
            _met_sumEt = -999.;
            _dphi_dipho_met = -999.;
            _pt_balance = -999.;
            _njet = -999.;
            _max_jet_pt = -999.;
            _max_jet_dCSV = -999.;
            _min_dphi_jet_met = -999.;

            reader.reset(new TMVA::Reader("!Color:!Silent"));
            reader->AddVariable("pho1_eta", &_pho1_eta);
            reader->AddVariable("pho2_eta", &_pho2_eta);
            reader->AddVariable("pho1_ptoM", &_pho1_ptoM);
            reader->AddVariable("pho2_ptoM", &_pho2_ptoM);
            reader->AddVariable("min_phoId", &_min_phoId);
            reader->AddVariable("max_phoId", &_max_phoId);
            reader->AddVariable("dipho_cosphi", &_dipho_cosphi);
            reader->AddVariable("met", &_met);
            reader->AddVariable("met_sumEt", &_met_sumEt);
            reader->AddVariable("dphi_dipho_met", &_dphi_dipho_met);
            reader->AddVariable("pt_balance", &_pt_balance);
            reader->AddVariable("njet", &_njet);
            reader->AddVariable("max_jet_pt", &_max_jet_pt);
            reader->AddVariable("max_jet_dCSV", &_max_jet_dCSV);
            reader->AddVariable("min_dphi_jet_met", &_min_dphi_jet_met);
            reader->BookMVA("BDT", weight_file);

            // Loop over Event, and calculate VH BDTScore
            for (int evt_entry = 0; evt_entry < tree_main->GetEntries(); ++evt_entry) {
                tree_main->GetEvent(evt_entry);

                Position4_t photon1, photon2, diphoton;
                photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                diphoton = photon1 + photon2;

                // *Put event preselection here.
                if (!(Event_Preselection_sig(photon1, photon2, met_Pt, met_Phi))) continue;
                float Max_PhotonID_Origin = std::max(dipho_leadIDMVA, dipho_subleadIDMVA);
                float Min_PhotonID_Origin = std::min(dipho_leadIDMVA, dipho_subleadIDMVA);

                if (dipho_leadIDMVA < -0.7 || dipho_subleadIDMVA < -0.7) continue;

                vector<Position4_t> vec_jet;
                for (int jetarray_entry = 0; jetarray_entry < n_jets; ++jetarray_entry) {
                    Position4_t jet;
                    jet.SetCoordinates(jet_Pt[jetarray_entry], jet_Eta[jetarray_entry], jet_Phi[jetarray_entry], jet_E[jetarray_entry]);
                    if (jet.Pt() < 20.) continue;
                    if (fabs(jet.Eta()) > 2.4) continue;
                    vec_jet.push_back(jet);
                }

                vector<float> jetpt, dphijetmet;
                float max_jet_pt, max_jet_dCSV, min_dphi_jet_met;
                if (vec_jet.size() == 0) {
                    max_jet_pt = -1.;
                    max_jet_dCSV = -2.;
                    min_dphi_jet_met = 4.;
                } else {
                    for (int jetvecentry = 0; jetvecentry < n_jets; ++jetvecentry) {
                        jetpt.push_back(vec_jet[jetvecentry].Pt());
                        dphijetmet.push_back(fabs(DeltaPhi(vec_jet[jetvecentry].Phi(), met_Phi)));
                    }
                    max_jet_pt = *std::max_element(jetpt.begin(), jetpt.end());
                    max_jet_dCSV = *std::max_element(jet_deepbtag, jet_deepbtag + 15);
                    min_dphi_jet_met = *std::min_element(dphijetmet.begin(), dphijetmet.end());
                }

                // assign TMVA featuring variables
                _pho1_eta = photon1.Eta();
                _pho2_eta = photon2.Eta();
                _pho1_ptoM = photon1.Pt() / diphoton.M();
                _pho2_ptoM = photon2.Pt() / diphoton.M();
                _min_phoId = Min_PhotonID_Origin;
                _max_phoId = Max_PhotonID_Origin;
                _dipho_cosphi = cos(DeltaPhi(photon1.Phi(), photon2.Phi()));
                _njet = n_jets;
                _max_jet_pt = max_jet_pt;
                _max_jet_dCSV = max_jet_dCSV;
                _met = met_Pt;
                _met_sumEt = met_sumEt;
                _min_dphi_jet_met = min_dphi_jet_met;
                _dphi_dipho_met = fabs(DeltaPhi(met_Phi, diphoton.Phi()));
                _pt_balance = (diphoton.Pt() - met_Pt) / diphoton.Pt();

                // Calculate the BDTScore
                float bdtScore_sig = reader->EvaluateMVA("BDT");
                float prob_sig = ConvertProbability(bdtScore_sig);
                sig_bdtscore_histo->Fill(bdtScore_sig, weight);
                sig_probability_histo->Fill(prob_sig, weight);
            }

            // *********************************************************************************************************** //
            // * Drawing anomoulous plot                                                                                   //
            // *********************************************************************************************************** //
            // Booking Histogram
            TH1F *sigac_bdtscore_histo = new TH1F("sigac_bdtscore_histo", " ", 50, -1., 1.);
            TH1F *sigac_probability_histo = new TH1F("sigac_probability_histo", " ", 50, 0., 1.);

            // Import file
            string filename_ac;
            if (sig_type == 0) {
                filename_ac = AbsRoute + filetype_list_zhJHUGenV7[AC_type];
            } else {
                filename_ac = AbsRoute + filetype_list_whJHUGenV7[AC_type];
            }
            auto file_ac = std::unique_ptr<TFile>(new TFile(filename_ac.c_str()));
            cout << "Input File: " << filename_ac << "\n";
            cout << "Input Weight File: " << weight_file << "\n";
            TTree *tree_ac;
            file_ac->GetObject(mc_routename_JHUGenV7[sig_type].c_str(), tree_ac);

            SettingBranch_Photon(tree_ac);
            SettingBranch_MET(tree_ac);
            SettingBranch_Jet(tree_ac);

            // Setup TMVA Reader
            std::unique_ptr<TMVA::Reader> reader_ac;

            float _pho1_eta_ac, _pho2_eta_ac;
            float _pho1_ptoM_ac, _pho2_ptoM_ac;
            float _min_phoId_ac, _max_phoId_ac;
            float _dipho_cosphi_ac;
            float _met_ac, _met_sumEt_ac;
            float _dphi_dipho_met_ac;
            float _pt_balance_ac;
            float _njet_ac;
            float _max_jet_pt_ac;
            float _max_jet_dCSV_ac;
            float _min_dphi_jet_met_ac;

            _pho1_eta_ac = -999.;
            _pho2_eta_ac = -999.;
            _pho1_ptoM_ac = -999.;
            _pho2_ptoM_ac = -999.;
            _min_phoId_ac = -999.;
            _max_phoId_ac = -999.;
            _dipho_cosphi_ac = -999.;
            _met_ac = -999.;
            _met_sumEt_ac = -999.;
            _dphi_dipho_met_ac = -999.;
            _pt_balance_ac = -999.;
            _njet_ac = -999.;
            _max_jet_pt_ac = -999.;
            _max_jet_dCSV_ac = -999.;
            _min_dphi_jet_met_ac = -999.;

            reader_ac.reset(new TMVA::Reader("!Color:!Silent"));
            reader_ac->AddVariable("pho1_eta", &_pho1_eta_ac);
            reader_ac->AddVariable("pho2_eta", &_pho2_eta_ac);
            reader_ac->AddVariable("pho1_ptoM", &_pho1_ptoM_ac);
            reader_ac->AddVariable("pho2_ptoM", &_pho2_ptoM_ac);
            reader_ac->AddVariable("min_phoId", &_min_phoId_ac);
            reader_ac->AddVariable("max_phoId", &_max_phoId_ac);
            reader_ac->AddVariable("dipho_cosphi", &_dipho_cosphi_ac);
            reader_ac->AddVariable("met", &_met_ac);
            reader_ac->AddVariable("met_sumEt", &_met_sumEt_ac);
            reader_ac->AddVariable("dphi_dipho_met", &_dphi_dipho_met_ac);
            reader_ac->AddVariable("pt_balance", &_pt_balance_ac);
            reader_ac->AddVariable("njet", &_njet_ac);
            reader_ac->AddVariable("max_jet_pt", &_max_jet_pt_ac);
            reader_ac->AddVariable("max_jet_dCSV", &_max_jet_dCSV_ac);
            reader_ac->AddVariable("min_dphi_jet_met", &_min_dphi_jet_met_ac);
            reader_ac->BookMVA("BDT", weight_file);

            // Loop over Event, and calculate VH BDTScore
            for (int evt_entry = 0; evt_entry < tree_ac->GetEntries(); ++evt_entry) {
                tree_ac->GetEvent(evt_entry);

                Position4_t photon1, photon2, diphoton;
                photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                diphoton = photon1 + photon2;

                // *Put event preselection here.
                if (!(Event_Preselection_sig(photon1, photon2, met_Pt, met_Phi))) continue;
                float Max_PhotonID_Origin = std::max(dipho_leadIDMVA, dipho_subleadIDMVA);
                float Min_PhotonID_Origin = std::min(dipho_leadIDMVA, dipho_subleadIDMVA);

                if (dipho_leadIDMVA < -0.7 || dipho_subleadIDMVA < -0.7) continue;

                vector<Position4_t> vec_jet;
                for (int jetarray_entry = 0; jetarray_entry < n_jets; ++jetarray_entry) {
                    Position4_t jet;
                    jet.SetCoordinates(jet_Pt[jetarray_entry], jet_Eta[jetarray_entry], jet_Phi[jetarray_entry], jet_E[jetarray_entry]);
                    if (jet.Pt() < 20.) continue;
                    if (fabs(jet.Eta()) > 2.4) continue;
                    vec_jet.push_back(jet);
                }

                vector<float> jetpt, dphijetmet;
                float max_jet_pt, max_jet_dCSV, min_dphi_jet_met;
                if (vec_jet.size() == 0) {
                    max_jet_pt = -1.;
                    max_jet_dCSV = -2.;
                    min_dphi_jet_met = 4.;
                } else {
                    for (int jetvecentry = 0; jetvecentry < n_jets; ++jetvecentry) {
                        jetpt.push_back(vec_jet[jetvecentry].Pt());
                        dphijetmet.push_back(fabs(DeltaPhi(vec_jet[jetvecentry].Phi(), met_Phi)));
                    }
                    max_jet_pt = *std::max_element(jetpt.begin(), jetpt.end());
                    max_jet_dCSV = *std::max_element(jet_deepbtag, jet_deepbtag + 15);
                    min_dphi_jet_met = *std::min_element(dphijetmet.begin(), dphijetmet.end());
                }

                // assign TMVA featuring variables
                _pho1_eta_ac = photon1.Eta();
                _pho2_eta_ac = photon2.Eta();
                _pho1_ptoM_ac = photon1.Pt() / diphoton.M();
                _pho2_ptoM_ac = photon2.Pt() / diphoton.M();
                _min_phoId_ac = Min_PhotonID_Origin;
                _max_phoId_ac = Max_PhotonID_Origin;
                _dipho_cosphi_ac = cos(DeltaPhi(photon1.Phi(), photon2.Phi()));
                _njet_ac = n_jets;
                _max_jet_pt_ac = max_jet_pt;
                _max_jet_dCSV_ac = max_jet_dCSV;
                _met_ac = met_Pt;
                _met_sumEt_ac = met_sumEt;
                _min_dphi_jet_met_ac = min_dphi_jet_met;
                _dphi_dipho_met_ac = fabs(DeltaPhi(met_Phi, diphoton.Phi()));
                _pt_balance_ac = (diphoton.Pt() - met_Pt) / diphoton.Pt();

                // Calculate the BDTScore
                float bdtScore_ac = reader_ac->EvaluateMVA("BDT");
                float prob_ac = ConvertProbability(bdtScore_ac);

                sigac_bdtscore_histo->Fill(bdtScore_ac, weight);
                sigac_probability_histo->Fill(prob_ac, weight);
            }

            // *************************************************************************************************************** //
            // * Drawing SM background VHMET BDT score plots (without gamma plus jet background)                               //
            // *************************************************************************************************************** //
            // Booking Histograms
            TH1F *Diboson_bdtscore_histo = new TH1F("Diboson_bdtscore_histo", " ", 50, -1, 1);
            TH1F *Diphoton_bdtscore_histo = new TH1F("Diphoton_bdtscore_histo", " ", 50, -1, 1);
            TH1F *ggh_bdtscore_histo = new TH1F("ggh_bdtscore_histo", " ", 50, -1, 1);
            TH1F *tth_bdtscore_histo = new TH1F("tth_bdtscore_histo", " ", 50, -1, 1);
            TH1F *vbf_bdtscore_histo = new TH1F("vbf_bdtscore_histo", " ", 50, -1, 1);
            TH1F *DY_bdtscore_histo = new TH1F("DY_bdtscore_histo", " ", 50, -1, 1);
            TH1F *Top_bdtscore_histo = new TH1F("Top_bdtscore_histo", " ", 50, -1, 1);
            TH1F *gplusjets_bdtscore_histo = new TH1F("gplusjets_bdtscore_histo", " ", 50, -1, 1);

            TH1F *Diboson_probability_histo = new TH1F("Diboson_probability_histo", " ", 50, 0, 1);
            TH1F *Diphoton_probability_histo = new TH1F("Diphoton_probability_histo", " ", 50, 0, 1);
            TH1F *ggh_probability_histo = new TH1F("ggh_probability_histo", " ", 50, 0, 1);
            TH1F *tth_probability_histo = new TH1F("tth_probability_histo", " ", 50, 0, 1);
            TH1F *vbf_probability_histo = new TH1F("vbf_probability_histo", " ", 50, 0, 1);
            TH1F *DY_probability_histo = new TH1F("DY_probability_histo", " ", 50, 0, 1);
            TH1F *Top_probability_histo = new TH1F("Top_probability_histo", " ", 50, 0, 1);
            TH1F *gplusjets_probability_histo = new TH1F("gplusjets_probability_histo", " ", 50, 0, 1);

            vector<TH1 *> histograms_bdtscore_fillorder = {Diboson_bdtscore_histo, Diphoton_bdtscore_histo, ggh_bdtscore_histo,
                                                           tth_bdtscore_histo, vbf_bdtscore_histo, DY_bdtscore_histo,
                                                           Top_bdtscore_histo, gplusjets_bdtscore_histo};
            vector<TH1 *> histograms_probability_fillorder = {Diboson_probability_histo, Diphoton_probability_histo, ggh_probability_histo,
                                                              tth_probability_histo, vbf_probability_histo, DY_probability_histo,
                                                              Top_probability_histo, gplusjets_probability_histo};

            vector<vector<string>> bkg_list;
            bkg_list.push_back(filetype_list_Diboson);
            bkg_list.push_back(filetype_list_Diphoton);
            bkg_list.push_back(filetype_list_ggh);
            bkg_list.push_back(filetype_list_tth);
            bkg_list.push_back(filetype_list_vbf);
            bkg_list.push_back(filetype_list_DY);
            bkg_list.push_back(filetype_list_Topbkg);

            // Looping the filetype_list and reading all the background file
            for (int bkg_type = 0; bkg_type < bkg_list.size(); bkg_type++) {
                // Booking subhistogram
                TH1F *bkg_bdtscore_subhisto[bkg_list[bkg_type].size()];
                TH1F *bkg_probability_subhisto[bkg_list[bkg_type].size()];
                for (int filelist_entry = 0; filelist_entry < bkg_list[bkg_type].size(); filelist_entry++) {
                    // Import file
                    string filename = AbsRoute + bkg_list[bkg_type][filelist_entry];
                    auto file = std::unique_ptr<TFile>{new TFile(filename.c_str())};
                    cout << "Input Background File: " << filename << "\n";
                    TTree *tree;
                    file->GetObject(mc_routename_bkg[bkg_type].c_str(), tree);
                    SettingBranch_Photon(tree);
                    SettingBranch_MET(tree);
                    SettingBranch_Jet(tree);

                    bkg_bdtscore_subhisto[filelist_entry] = new TH1F(Form("bkg_bdtscore_subhisto%d", filelist_entry), " ", 50, -1, 1);
                    bkg_probability_subhisto[filelist_entry] = new TH1F(Form("bkg_probability_subhisto%d", filelist_entry), " ", 50, 0, 1);

                    // Setup TMVA Reader
                    std::unique_ptr<TMVA::Reader> reader;

                    float _pho1_eta, _pho2_eta;
                    float _pho1_ptoM, _pho2_ptoM;
                    float _min_phoId, _max_phoId;
                    float _dipho_cosphi;
                    float _met, _met_sumEt;
                    float _dphi_dipho_met;
                    float _pt_balance;
                    float _njet;
                    float _max_jet_pt;
                    float _max_jet_dCSV;
                    float _min_dphi_jet_met;

                    _pho1_eta = -999.;
                    _pho2_eta = -999.;
                    _pho1_ptoM = -999.;
                    _pho2_ptoM = -999.;
                    _min_phoId = -999.;
                    _max_phoId = -999.;
                    _dipho_cosphi = -999.;
                    _met = -999.;
                    _met_sumEt = -999.;
                    _dphi_dipho_met = -999.;
                    _pt_balance = -999.;
                    _njet = -999.;
                    _max_jet_pt = -999.;
                    _max_jet_dCSV = -999.;
                    _min_dphi_jet_met = -999.;

                    reader.reset(new TMVA::Reader("!Color:!Silent"));
                    reader->AddVariable("pho1_eta", &_pho1_eta);
                    reader->AddVariable("pho2_eta", &_pho2_eta);
                    reader->AddVariable("pho1_ptoM", &_pho1_ptoM);
                    reader->AddVariable("pho2_ptoM", &_pho2_ptoM);
                    reader->AddVariable("min_phoId", &_min_phoId);
                    reader->AddVariable("max_phoId", &_max_phoId);
                    reader->AddVariable("dipho_cosphi", &_dipho_cosphi);
                    reader->AddVariable("met", &_met);
                    reader->AddVariable("met_sumEt", &_met_sumEt);
                    reader->AddVariable("dphi_dipho_met", &_dphi_dipho_met);
                    reader->AddVariable("pt_balance", &_pt_balance);
                    reader->AddVariable("njet", &_njet);
                    reader->AddVariable("max_jet_pt", &_max_jet_pt);
                    reader->AddVariable("max_jet_dCSV", &_max_jet_dCSV);
                    reader->AddVariable("min_dphi_jet_met", &_min_dphi_jet_met);
                    reader->BookMVA("BDT", weight_file);

                    // Loop over Event, and calculate VH BDTScore
                    for (int evt_entry = 0; evt_entry < tree->GetEntries(); evt_entry++) {
                        tree->GetEntry(evt_entry);

                        Position4_t photon1, photon2, diphoton;
                        photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                        photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                        diphoton = photon1 + photon2;

                        // *Put event preselection here.
                        if (!(Event_Preselection_bkg(photon1, photon2, met_Pt, met_Phi))) continue;
                        if (dipho_leadIDMVA < -0.7 || dipho_subleadIDMVA < -0.7) continue;

                        float Max_PhotonID_Origin, Min_PhotonID_Origin;
                        Max_PhotonID_Origin = std::max(dipho_leadIDMVA, dipho_subleadIDMVA);
                        Min_PhotonID_Origin = std::min(dipho_leadIDMVA, dipho_subleadIDMVA);

                        vector<Position4_t> vec_jet;
                        for (int jetarray_entry = 0; jetarray_entry < n_jets; ++jetarray_entry) {
                            Position4_t jet;
                            jet.SetCoordinates(jet_Pt[jetarray_entry], jet_Eta[jetarray_entry], jet_Phi[jetarray_entry], jet_E[jetarray_entry]);
                            if (jet.Pt() < 20.) continue;
                            if (fabs(jet.Eta()) > 2.4) continue;
                            vec_jet.push_back(jet);
                        }

                        vector<float> jetpt, dphijetmet;
                        float max_jet_pt, max_jet_dCSV, min_dphi_jet_met;
                        if (vec_jet.size() == 0) {
                            max_jet_pt = -1.;
                            max_jet_dCSV = -2.;
                            min_dphi_jet_met = 4.;
                        } else {
                            for (int jetvecentry = 0; jetvecentry < n_jets; ++jetvecentry) {
                                jetpt.push_back(vec_jet[jetvecentry].Pt());
                                dphijetmet.push_back(fabs(DeltaPhi(vec_jet[jetvecentry].Phi(), met_Phi)));
                            }
                            max_jet_pt = *std::max_element(jetpt.begin(), jetpt.end());
                            max_jet_dCSV = *std::max_element(jet_deepbtag, jet_deepbtag + 15);
                            min_dphi_jet_met = *std::min_element(dphijetmet.begin(), dphijetmet.end());
                        }

                        // assign TMVA featuring variables
                        _pho1_eta = photon1.Eta();
                        _pho2_eta = photon2.Eta();
                        _pho1_ptoM = photon1.Pt() / diphoton.M();
                        _pho2_ptoM = photon2.Pt() / diphoton.M();
                        _min_phoId = Min_PhotonID_Origin;
                        _max_phoId = Max_PhotonID_Origin;
                        _dipho_cosphi = cos(DeltaPhi(photon1.Phi(), photon2.Phi()));
                        _njet = n_jets;
                        _max_jet_pt = max_jet_pt;
                        _max_jet_dCSV = max_jet_dCSV;
                        _met = met_Pt;
                        _met_sumEt = met_sumEt;
                        _min_dphi_jet_met = min_dphi_jet_met;
                        _dphi_dipho_met = fabs(DeltaPhi(met_Phi, diphoton.Phi()));
                        _pt_balance = (diphoton.Pt() - met_Pt) / diphoton.Pt();

                        // output bdt score and fill the histogram
                        float bdtScore_bkg = reader->EvaluateMVA("BDT");
                        float prob_bkg = ConvertProbability(bdtScore_bkg);
                        bkg_bdtscore_subhisto[filelist_entry]->Fill(bdtScore_bkg, weight * RR18_Lumi);
                        bkg_probability_subhisto[filelist_entry]->Fill(prob_bkg, weight * RR18_Lumi);
                    }
                    // Add the subhistogram into main histograms
                    histograms_bdtscore_fillorder[bkg_type]->Add(bkg_bdtscore_subhisto[filelist_entry]);
                    histograms_probability_fillorder[bkg_type]->Add(bkg_probability_subhisto[filelist_entry]);
                    delete bkg_bdtscore_subhisto[filelist_entry];
                    delete bkg_probability_subhisto[filelist_entry];
                }
            }

            // *************************************************************************************************************** //
            // * Drawing SM background VHMET BDT score plots: gamma plus jet background(consider data-driven scale factor)     //
            // *************************************************************************************************************** //
            // * Store the information of data's max photon ID in sideband region
            // Import file
            cout << "Store data's max photon ID in sideband region. \n";
            std::string filename_data = AbsRoute + filetype_data[0];
            auto file_data = std::unique_ptr<TFile>{new TFile(filename_data.c_str())};
            cout << "Input Data File: " << filename_data << "\n";
            TTree *tree_data;
            file_data->GetObject("tagsDumper/trees/Data_13TeV_RECO_VH_MET_Tag0", tree_data);

            SettingBranch_Photon(tree_data);
            // *Store information Max photon ID locates in sideband region of data in vector: MaxPhoID_Data_SidebandRegion
            // *This information will be used in generating gamma plus jet background (data-dirven method)
            vector<int> MaxPhoID_Data_SidebandRegion_No = {};
            vector<float> MaxPhoID_Data_SidebandRegion = {};

            // Loop over Event, and calculate VH BDTScore
            for (int evt_entry = 0; evt_entry < tree_data->GetEntries(); ++evt_entry) {
                tree_data->GetEvent(evt_entry);

                Position4_t photon1, photon2, diphoton;
                photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                diphoton = photon1 + photon2;

                // *Put event preselection here.
                if (!(Event_Preselection_bkg(photon1, photon2, met_Pt, met_Phi))) continue;
                float Max_PhotonID_Origin = std::max(dipho_leadIDMVA, dipho_subleadIDMVA);
                float Min_PhotonID_Origin = std::min(dipho_leadIDMVA, dipho_subleadIDMVA);

                // Store the information of data's max photon ID in sideband region
                if ((Min_PhotonID_Origin > -1. && Min_PhotonID_Origin < -0.7) && (Max_PhotonID_Origin > -0.7 && Max_PhotonID_Origin < 1.)) {
                    MaxPhoID_Data_SidebandRegion.push_back(Max_PhotonID_Origin);
                    MaxPhoID_Data_SidebandRegion_No.push_back(evt_entry);
                }
            }

            cout << "Photon Sideband Region Size: " << MaxPhoID_Data_SidebandRegion.size() << "\n";
            cout << "Photon Sideband Region Number Size: " << MaxPhoID_Data_SidebandRegion_No.size() << "\n";

            bkg_list.push_back(filetype_list_gplusjets);

            // Booking histogram and realtive subhistogram for fake photon ID from gamma plus jets
            TH1F *FakePhotonGJets_histo = new TH1F("FakePhotonGJets_histo", " ", 40, -0.9, 1.);
            TH1F *FakePhotonGJets_subhisto[filetype_list_gplusjets.size()];

            // Reading the gamma plus jets files to get fake photon ID PDF
            for (int filelist_entry = 0; filelist_entry < filetype_list_gplusjets.size(); ++filelist_entry) {
                string filename_fake = AbsRoute + filetype_list_gplusjets[filelist_entry];
                auto file_fake = std::unique_ptr<TFile>{new TFile(filename_fake.c_str())};
                TTree *tree_fake;
                file_fake->GetObject("tagsDumper/trees/gjet_fakephoton_13TeV_RECO_VH_MET_Tag0", tree_fake);
                SettingBranch_Photon(tree_fake);
                SettingBranch_MET(tree_fake);
                float dipho_lead_prompt, dipho_sublead_prompt;
                tree_fake->SetBranchAddress("dipho_lead_prompt", &dipho_lead_prompt);
                tree_fake->SetBranchAddress("dipho_sublead_prompt", &dipho_sublead_prompt);

                // Booking the subhistogram
                FakePhotonGJets_subhisto[filelist_entry] = new TH1F(Form("FakePhotonGJets_subhisto%d", filelist_entry), " ", 40, -0.9, 1.);

                // Loop over the events
                for (int evt_entry = 0; evt_entry < tree_fake->GetEntries(); ++evt_entry) {
                    tree_fake->GetEntry(evt_entry);

                    // Put offline selection for background.
                    Position4_t photon1, photon2;
                    photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                    photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                    Event_Preselection_bkg(photon1, photon2, met_Pt, met_Phi);
                    if (dipho_leadIDMVA < -0.9 || dipho_subleadIDMVA < -0.9) continue;

                    // Selecting the fake photon candidates
                    if (dipho_lead_prompt != 1. && dipho_sublead_prompt != 1.) {
                        FakePhotonGJets_subhisto[filelist_entry]->Fill(std::min(dipho_leadIDMVA, dipho_subleadIDMVA), weight * RR18_Lumi);
                    } else if (dipho_lead_prompt != 1. && dipho_sublead_prompt == 1.) {
                        FakePhotonGJets_subhisto[filelist_entry]->Fill(dipho_leadIDMVA, weight * RR18_Lumi);
                    } else if (dipho_lead_prompt == 1. && dipho_sublead_prompt != 1.) {
                        FakePhotonGJets_subhisto[filelist_entry]->Fill(dipho_subleadIDMVA, weight * RR18_Lumi);
                    } else
                        continue;
                }
                FakePhotonGJets_histo->Add(FakePhotonGJets_subhisto[filelist_entry]);
            }
            FakePhotonGJets_histo->Scale(1. / FakePhotonGJets_histo->Integral());

            // Now we try to use 7th order polynomial to fit histogram to get fake photon PDF.
            TF1 *FakePhotonPDF_fit = new TF1("FakePhotonPDF_fit", "pol7", -0.9, 1.);
            FakePhotonGJets_histo->Fit(FakePhotonPDF_fit);
            // *Here we using max photon id in data sideband region with pre-event weight to calculate gamma plus jets background's vh bdt score
            // *And the min photon is generate from fake photon PDF within region [-0.7, max photon id]
            // Import file
            std::string filename_main2 = AbsRoute + filetype_data[0];
            auto file_main2 = std::unique_ptr<TFile>{new TFile(filename_main2.c_str())};
            cout << "Input File: " << filename_main2 << "\n";
            TTree *tree_main2;
            file_main2->GetObject("tagsDumper/trees/Data_13TeV_RECO_VH_MET_Tag0", tree_main2);

            SettingBranch_Photon(tree_main2);
            SettingBranch_MET(tree_main2);
            SettingBranch_Jet(tree_main2);

            // Setup TMVA Reader
            std::unique_ptr<TMVA::Reader> reader2;
            /// TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

            float _pho1_eta_dd, _pho2_eta_dd;
            float _pho1_ptoM_dd, _pho2_ptoM_dd;
            float _min_phoId_dd, _max_phoId_dd;
            float _dipho_cosphi_dd;
            float _met_dd, _met_sumEt_dd;
            float _dphi_dipho_met_dd;
            float _pt_balance_dd;
            float _njet_dd;
            float _max_jet_pt_dd;
            float _max_jet_dCSV_dd;
            float _min_dphi_jet_met_dd;

            _pho1_eta_dd = -999.;
            _pho2_eta_dd = -999.;
            _pho1_ptoM_dd = -999.;
            _pho2_ptoM_dd = -999.;
            _min_phoId_dd = -999.;
            _max_phoId_dd = -999.;
            _dipho_cosphi_dd = -999.;
            _met_dd = -999.;
            _met_sumEt_dd = -999.;
            _dphi_dipho_met_dd = -999.;
            _pt_balance_dd = -999.;
            _njet_dd = -999.;
            _max_jet_pt_dd = -999.;
            _max_jet_dCSV_dd = -999.;
            _min_dphi_jet_met_dd = -999.;

            reader2.reset(new TMVA::Reader("!Color:!Silent"));
            reader2->AddVariable("pho1_eta", &_pho1_eta_dd);
            reader2->AddVariable("pho2_eta", &_pho2_eta_dd);
            reader2->AddVariable("pho1_ptoM", &_pho1_ptoM_dd);
            reader2->AddVariable("pho2_ptoM", &_pho2_ptoM_dd);
            reader2->AddVariable("min_phoId", &_min_phoId_dd);
            reader2->AddVariable("max_phoId", &_max_phoId_dd);
            reader2->AddVariable("dipho_cosphi", &_dipho_cosphi_dd);
            reader2->AddVariable("met", &_met_dd);
            reader2->AddVariable("met_sumEt", &_met_sumEt_dd);
            reader2->AddVariable("dphi_dipho_met", &_dphi_dipho_met_dd);
            reader2->AddVariable("pt_balance", &_pt_balance_dd);
            reader2->AddVariable("njet", &_njet_dd);
            reader2->AddVariable("max_jet_pt", &_max_jet_pt_dd);
            reader2->AddVariable("max_jet_dCSV", &_max_jet_dCSV_dd);
            reader2->AddVariable("min_dphi_jet_met", &_min_dphi_jet_met_dd);
            reader2->BookMVA("BDT", weight_file);

            // Loop over Event, and calculate VH BDTScore
            for (int evt_entry = 0; evt_entry < MaxPhoID_Data_SidebandRegion_No.size(); ++evt_entry) {
                tree_main2->GetEntry(MaxPhoID_Data_SidebandRegion_No[evt_entry]);

                Position4_t photon1, photon2, diphoton;
                photon1.SetCoordinates(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadE);
                photon2.SetCoordinates(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadE);
                diphoton = photon1 + photon2;
                float MaxPhotonID_Origin, pre_event_weight;
                float MaxPhotonID, MinPhotonID;
                MaxPhotonID_Origin = MaxPhoID_Data_SidebandRegion[evt_entry];
                pre_event_weight = FakePhotonPDF_fit->Integral(-0.7, MaxPhotonID_Origin) / FakePhotonPDF_fit->Integral(-0.9, -0.7);
                MaxPhotonID = MaxPhotonID_Origin * pre_event_weight;
                MinPhotonID = FakePhotonPDF_fit->GetRandom(-0.7, MaxPhotonID_Origin);

                // *Put event preselection here.
                if (!(Event_Preselection_bkg(photon1, photon2, met_Pt, met_Phi))) continue;

                vector<Position4_t> vec_jet;
                for (int jetarray_entry = 0; jetarray_entry < n_jets; ++jetarray_entry) {
                    Position4_t jet;
                    jet.SetCoordinates(jet_Pt[jetarray_entry], jet_Eta[jetarray_entry], jet_Phi[jetarray_entry], jet_E[jetarray_entry]);
                    if (jet.Pt() < 20.) continue;
                    if (fabs(jet.Eta()) > 2.4) continue;
                    vec_jet.push_back(jet);
                }

                vector<float> jetpt, dphijetmet;
                float max_jet_pt, max_jet_dCSV, min_dphi_jet_met;
                if (vec_jet.size() == 0) {
                    max_jet_pt = -1.;
                    max_jet_dCSV = -2.;
                    min_dphi_jet_met = 4.;
                } else {
                    for (int jetvecentry = 0; jetvecentry < n_jets; ++jetvecentry) {
                        jetpt.push_back(vec_jet[jetvecentry].Pt());
                        dphijetmet.push_back(fabs(DeltaPhi(vec_jet[jetvecentry].Phi(), met_Phi)));
                    }
                    max_jet_pt = *std::max_element(jetpt.begin(), jetpt.end());
                    max_jet_dCSV = *std::max_element(jet_deepbtag, jet_deepbtag + 15);
                    min_dphi_jet_met = *std::min_element(dphijetmet.begin(), dphijetmet.end());
                }

                // assign TMVA featuring variables
                _pho1_eta_dd = photon1.Eta();
                _pho2_eta_dd = photon2.Eta();
                _pho1_ptoM_dd = photon1.Pt() / diphoton.M();
                _pho2_ptoM_dd = photon2.Pt() / diphoton.M();
                _min_phoId_dd = MinPhotonID;
                _max_phoId_dd = MaxPhotonID;
                _dipho_cosphi_dd = cos(DeltaPhi(photon1.Phi(), photon2.Phi()));
                _njet_dd = n_jets;
                _max_jet_pt_dd = max_jet_pt;
                _max_jet_dCSV_dd = max_jet_dCSV;
                _met_dd = met_Pt;
                _met_sumEt_dd = met_sumEt;
                _min_dphi_jet_met_dd = min_dphi_jet_met;
                _dphi_dipho_met_dd = fabs(DeltaPhi(met_Phi, diphoton.Phi()));
                _pt_balance_dd = (diphoton.Pt() - met_Pt) / diphoton.Pt();

                // output bdt score and fill the histogram
                float bdtScore_gjets = reader2->EvaluateMVA("BDT");
                float prob_gjets = ConvertProbability(bdtScore_gjets);
                gplusjets_bdtscore_histo->Fill(bdtScore_gjets, pre_event_weight);
                gplusjets_probability_histo->Fill(prob_gjets, pre_event_weight);
            }

            // *********************************************************************************************************** //
            // * Plotting                                                                                                  //
            // *********************************************************************************************************** //
            // *BDT score: sig
            c_bdtscore->cd();
            sig_bdtscore_histo->SetStats(0);
            sig_bdtscore_histo->SetTitle(" ; VHMET BDTScore; Events");
            sig_bdtscore_histo->GetYaxis()->SetTitleSize(0.04);
            sig_bdtscore_histo->Scale(500.);
            HistogramLabelOption(sig_bdtscore_histo, 0.025, 1.2, 1.05);
            HistogramLimitOption(sig_bdtscore_histo, 0.05, 500000);
            HistogramLineOption(sig_bdtscore_histo, 2, kAzure, 1, 0.5);
            sig_bdtscore_histo->Draw("hist");
            leg_bdtscore->AddEntry(sig_bdtscore_histo, "(f_{a1}=1.0 (SM)) #times 500", "f");
            LegendOption(leg_bdtscore, 0.025);

            // *BDT score: ac
            HistogramLineOption(sigac_bdtscore_histo, 2, kRed + 1, 1, 0.5);
            sigac_bdtscore_histo->Scale(500.);
            leg_bdtscore->AddEntry(sigac_bdtscore_histo, ac_histo_ledgendname[AC_type].c_str(), "f");

            // *BDT score: SM background
            Diphoton_bdtscore_histo->Scale(scale_factor_DiPhoton);
            gplusjets_bdtscore_histo->Scale(scale_factor_GJets);
            histograms_bdtscore_fillorder[1] = Diphoton_bdtscore_histo;
            histograms_bdtscore_fillorder[7] = gplusjets_bdtscore_histo;
            Make_Stack(bkg_bdtscore_hs, histograms_bdtscore_fillorder, histo_colours, leg_bdtscore, filetype_legname);
            bkg_bdtscore_hs->Draw("hist, same");
            sig_bdtscore_histo->Draw("hist, same");
            sigac_bdtscore_histo->Draw("hist, same");

            leg_bdtscore->Draw();
            if (sig_type == 0) {
                c_bdtscore->SaveAs((file_naming_zh[AC_type] + "_bdtscore.pdf").c_str());
            } else {
                c_bdtscore->SaveAs((file_naming_wh[AC_type] + "_bdtscore.pdf").c_str());
            }

            // *Probability: sig
            c_probability->cd();
            sig_probability_histo->SetStats(0);
            sig_probability_histo->SetTitle(" ; Probability; Events");
            sig_probability_histo->GetYaxis()->SetTitleSize(0.04);
            sig_probability_histo->Scale(500.);
            HistogramLabelOption(sig_probability_histo, 0.025, 1.2, 1.05);
            HistogramLimitOption(sig_probability_histo, 0.05, 500000);
            HistogramLineOption(sig_probability_histo, 2, kAzure, 1, 0.5);
            sig_probability_histo->Draw("hist");
            leg_probability->AddEntry(sig_probability_histo, "(f_{a1}=1.0 (SM)) #times 500", "f");
            LegendOption(leg_probability, 0.025);

            // *Probability: ac
            HistogramLineOption(sigac_probability_histo, 2, kRed + 1, 1, 0.5);
            sigac_probability_histo->Scale(500.);
            leg_probability->AddEntry(sigac_probability_histo, ac_histo_ledgendname[AC_type].c_str(), "f");

            // *Probability: SM background
            Diphoton_probability_histo->Scale(scale_factor_DiPhoton);
            gplusjets_probability_histo->Scale(scale_factor_GJets);
            histograms_probability_fillorder[1] = Diphoton_probability_histo;
            histograms_probability_fillorder[7] = gplusjets_probability_histo;
            Make_Stack(bkg_probability_hs, histograms_probability_fillorder, histo_colours, leg_probability, filetype_legname);
            bkg_probability_hs->Draw("hist, same");
            sig_probability_histo->Draw("hist, same");
            sigac_probability_histo->Draw("hist, same");

            leg_probability->Draw();
            if (sig_type == 0) {
                c_probability->SaveAs((file_naming_zh[AC_type] + "_probability.pdf").c_str());
            } else {
                c_probability->SaveAs((file_naming_wh[AC_type] + "_probability.pdf").c_str());
            }

            // Delete the plots and histograms in the end of loop
            delete c_bdtscore;
            delete leg_bdtscore;
            delete bkg_bdtscore_hs;
            delete sig_bdtscore_histo;
            delete sigac_bdtscore_histo;
            delete c_probability;
            delete leg_probability;
            delete bkg_probability_hs;
            delete sig_probability_histo;
            delete sigac_probability_histo;
        }
    }
}

void SettingBranch_Photon(TTree *tree) {
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("dipho_leadPt", &dipho_leadPt);
    tree->SetBranchAddress("dipho_leadEta", &dipho_leadEta);
    tree->SetBranchAddress("dipho_leadPhi", &dipho_leadPhi);
    tree->SetBranchAddress("dipho_leadE", &dipho_leadE);
    tree->SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA);

    tree->SetBranchAddress("dipho_subleadPt", &dipho_subleadPt);
    tree->SetBranchAddress("dipho_subleadEta", &dipho_subleadEta);
    tree->SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi);
    tree->SetBranchAddress("dipho_subleadE", &dipho_subleadE);
    tree->SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA);
}

void SettingBranch_MET(TTree *tree) {
    tree->SetBranchAddress("met_Pt", &met_Pt);
    tree->SetBranchAddress("met_Phi", &met_Phi);
    tree->SetBranchAddress("met_sumEt", &met_sumEt);
}

void SettingBranch_Jet(TTree *tree) {
    tree->SetBranchAddress("n_jets", &n_jets);
    for (int jet_num = 0; jet_num < 15; jet_num++) {
        tree->SetBranchAddress(Form("jet%d_Pt", jet_num + 1), &jet_Pt[jet_num]);
        tree->SetBranchAddress(Form("jet%d_Eta", jet_num + 1), &jet_Eta[jet_num]);
        tree->SetBranchAddress(Form("jet%d_Phi", jet_num + 1), &jet_Phi[jet_num]);
        tree->SetBranchAddress(Form("jet%d_E", jet_num + 1), &jet_E[jet_num]);
        tree->SetBranchAddress(Form("jet%d_deepbtag", jet_num + 1), &jet_deepbtag[jet_num]);
    }
}

float ConvertProbability(float bdtscore) {
    return 1. / (1. + exp(0.5 * log(2. / (bdtscore + 1.) - 1)));
}