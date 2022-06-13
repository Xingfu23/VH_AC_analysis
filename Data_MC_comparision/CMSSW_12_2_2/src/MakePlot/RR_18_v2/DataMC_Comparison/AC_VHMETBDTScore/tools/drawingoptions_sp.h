#ifndef __DRAWINGOPTIONS_SP_H__
#define __DRAWINGOPTIONS_SP_H__

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

void HistogramLineOption(std::shared_ptr<TH1> histogram, int linewidth, Color_t linecolour, int linestyle, int fillcolour) {
    histogram->SetStats(0);
    histogram->SetLineWidth(linewidth);
    histogram->SetLineColor(linecolour);
    histogram->SetLineStyle(linestyle);
    histogram->SetFillColor(fillcolour);
};

void HistogramLimitOption(std::shared_ptr<TH1> histogram, float min, float max) {
    histogram->SetMinimum(min);
    histogram->SetMaximum(max);
};

void HistogramLabelOption(std::shared_ptr<TH1> histogram, float labelsize, float xlabeloffset, float ylabeloffset) {
    histogram->SetLabelSize(labelsize, "xy");
    histogram->GetXaxis()->SetTitleOffset(xlabeloffset);
    histogram->GetYaxis()->SetTitleOffset(ylabeloffset);
};

void HistogramDotOption(std::shared_ptr<TH1> histogram, float marker_style, float marker_size, float marker_colour = 1) {
    histogram->SetMarkerStyle(marker_style);
    histogram->SetMarkerSize(marker_size);
    histogram->SetMarkerColor(marker_colour);
}

void LegendOption(std::shared_ptr<TLegend> legend, float textsize, float bordersize = 0.) {
    legend->SetTextSize(textsize);
    legend->SetBorderSize(bordersize);
};

#endif  // __DRAWINGOPTIONS_SP_H__