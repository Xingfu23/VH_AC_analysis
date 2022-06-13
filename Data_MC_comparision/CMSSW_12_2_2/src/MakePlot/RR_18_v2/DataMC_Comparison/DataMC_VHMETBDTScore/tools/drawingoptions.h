#ifndef __DRAWINGOPTIONS_H__
#define __DRAWINGOPTIONS_H__

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TObject.h>
#include <TPad.h>
#include <TROOT.h>
#include <TTree.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

void HistogramLineOption(TH1 *histogram, int linewidth, Color_t linecolour, int linestyle, int fillcolour) {
    histogram->SetStats(0);
    histogram->SetLineWidth(linewidth);
    histogram->SetLineColor(linecolour);
    histogram->SetLineStyle(linestyle);
    histogram->SetFillColor(fillcolour);
};

void HistogramLimitOption(TH1 *histogram, float min, float max) {
    histogram->SetMinimum(min);
    histogram->SetMaximum(max);
};

void HistogramLabelOption(TH1 *histogram, float labelsize, float xlabeloffset, float ylabeloffset) {
    histogram->SetLabelSize(labelsize, "xy");
    histogram->GetXaxis()->SetTitleOffset(xlabeloffset);
    histogram->GetYaxis()->SetTitleOffset(ylabeloffset);
};

void HistogramDotOption(TH1 *histogram, float marker_style, float marker_size, float marker_colour = 1) {
    histogram->SetMarkerStyle(marker_style);
    histogram->SetMarkerSize(marker_size);
    histogram->SetMarkerColor(marker_colour);
}

void DrawingSMHiggs(TH1 *histogram_zh, TH1 *histogram_wh, TLegend *legend) {
    // index is 1 means drawing maximum photon id and index equals to 2 means minmum photon id
    HistogramLineOption(histogram_zh, 2, kRed + 1, 1, 0);
    HistogramLineOption(histogram_wh, 2, kRed + 1, 2, 0);
    histogram_zh->Draw("hist, same");
    histogram_wh->Draw("hist, same");
    legend->AddEntry(histogram_zh, "Z(#rightarrow#nu#nu)H#times 500", "f");
    legend->AddEntry(histogram_wh, "WH#times 500", "f");
};

void Make_Stack(THStack *hs, std::vector<TH1 *> v, std::vector<Color_t> colours, TLegend *legend, std::vector<std::string> label) {
    if (v.size() != colours.size() || v.size() != label.size()) {
        std::cout << "ERROR, please check the length of histograms list with colours or label list. They should be the same." << std::endl;
    } else {
        for (int iter = 0; iter < v.size(); ++iter) {
            HistogramLineOption(v[iter], 1, kBlack, 1, colours[iter]);
            legend->AddEntry(v[iter], label[iter].c_str(), "f");
            hs->Add(v[iter]);
        }
    }
}

void LegendOption(TLegend *legend, float textsize, float bordersize = 0.) {
    legend->SetTextSize(textsize);
    legend->SetBorderSize(bordersize);
};

void RatioPadSetting(TPad *pad) {
    pad->SetTickx(1);
    pad->SetTicky(1);
    pad->Draw();
    pad->cd();
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.2);
};

void RatioPlotSettings(TH1 *ratio_histogram, TH1 *total_histogram) {
    ratio_histogram->Sumw2();
    ratio_histogram->Divide(total_histogram);
    ratio_histogram->GetXaxis()->SetTitleSize(0.09);
    ratio_histogram->GetYaxis()->SetTitleSize(0.1);
    ratio_histogram->GetYaxis()->SetNdivisions(5);
    HistogramLineOption(ratio_histogram, 1, kPink, 1, 0);
    HistogramLabelOption(ratio_histogram, 0.1, 0.92, 0.4);
    HistogramLimitOption(ratio_histogram, 0., 2.1);
    ratio_histogram->SetMarkerStyle(8);
    ratio_histogram->Draw("p");
};

void AddText() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(.7, .92, "Era2018_RR_v2");
};

void AddText_Fake() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(.70, .92, "Era2018_RR_v2");
    latex.SetTextSize(0.04);
    latex.DrawLatex(.75, .85, "#gamma + jets");
};

void HistogramLabelOption2(std::shared_ptr<TH1> histogram, float labelsize, float xlabeloffset, float ylabeloffset) {
    histogram->SetLabelSize(labelsize, "xy");
    histogram->GetXaxis()->SetTitleOffset(xlabeloffset);
    histogram->GetYaxis()->SetTitleOffset(ylabeloffset);
};

#endif  // __DRAWINGOPTIONS_H__