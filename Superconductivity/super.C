#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <TString.h>
#include "TMath.h"
#include <TPaveStats.h>
#include <TText.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TGraphErrors* readData(TString);
void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);
void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);

TGraphErrors *graphLarge, *graphSmall, *allPeaks;

/*
 * MAIN FUNCTION
 */

void super(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    TGraphErrors *data = readData("magneticField");

    TF1 *f1 = new TF1("fitRight", "[0]+[1]*TMath::Abs(TMath::Sin( [2]*(x-[3]) )/( [2]*(x-[3]) ) )", -200, 200);
    f1->SetParameter(0, 0.);
    f1->SetParameter(1, 0.36);
    f1->SetParameter(2, 0.05);
    f1->SetParameter(3, -15);
    f1->SetParName(0, "C");
    f1->SetParName(1, "N");
    f1->SetParName(2, "#pi#Phi/#Phi_{0}I");
    f1->SetParName(3, "x_{0}");
    f1->SetLineColor(kBlue);

    TF1 *f2 = new TF1("fitLeft", "[0]+[1]*TMath::Abs(TMath::Sin( [2]*(x-[3]) )/( [2]*(x-[3]) ) )");
    f2->SetParameter(0, 50./2000.);
    f2->SetParameter(1, 700./2000.);
    f2->SetParameter(2, 0.06283);
    f2->SetParameter(3, -15);
    f2->SetParName(0, "C");
    f2->SetParName(1, "N");
    f2->SetParName(2, "A");
    f2->SetParName(3, "A");
    f2->SetLineColor(kBlue);

    // data->Fit(f1, "ME", "", 18, 60);
    // data->Fit(f2, "ME", "", -200, -33);

    data->Fit(f1, "ME", "", -90, 70);


    Double_t pi = 3.1415926538, L = 7.5e-6, l = 1.75e-9, lambda = 39.78e-9, conv = 0.54e-4, t=2*lambda+l;
    Double_t fQauntaRight = pi*L*t*conv/f1->GetParameter(2);//, fQauntaLeft = pi*L*t*conv/f2->GetParameter(2);
    Double_t fQuantaError = f1->GetParError(2)*pi*L*t*conv/f1->GetParameter(2)/f1->GetParameter(2);

    cout << "Flux quanta from right is " << fQauntaRight << " +- " << fQuantaError << endl;//<< " and from left is " << fQauntaLeft << endl;

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    graph(data, "", can, "Solenoid current (mA)", "Zero-voltage current (mA)", "magField");


}

TGraphErrors* readData(TString dataset){

    TString fileName = "./Data/" + dataset + ".txt";

    Int_t nEntries = 0;
    Double_t current = 0, voltage = 0;
    Double_t currentArr[1000], currentError[1000], voltageArr[1000], voltageError[1000];

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}
    
    while(ifs >> current >> voltage){

        if(!(current > -85 && current < -60) && !(current > 45 && current < 65)) continue;

        currentArr[nEntries] = current;
        currentError[nEntries] = 0.8;
        voltageArr[nEntries] = voltage/2000.;
        voltageError[nEntries] = 5./2000.;
        nEntries++;

    }

    ifs.close();

    TGraphErrors *gr = new TGraphErrors(nEntries, currentArr, voltageArr, currentError, voltageError);

    return gr;

}


/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );
    // gStyle->SetOptFit(kFALSE);
    gStyle->SetOptFit(1100);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("4.3f");

    if(!histoData) return;

    // TF1 *f = histoData->GetFunction("Peak_0");

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.6);
        histoData->SetMarkerColor(kRed);
        // histoData->SetMarkerColor(kBlue);
        histoData->SetLineColor(kBlack);
        // f->SetLineColor(kBlue);
        // f->SetLineWidth(4);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        // f->SetLineColor(kRed);
        // f->SetLineWidth(2);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.755,0.675,0.985,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data", (type.Contains("dot") ? "lep" : "l"));
    // leg->AddEntry(f, "Voigtian fits","l");
    
    histoData->Draw((type.Contains("dot") ? "E1" : ""));
    // leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.84);
    histoData->GetXaxis()->SetLabelOffset(-0.001);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(1.0);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
}

void histogram(TH1D *histoData1, TH1D *histoData2, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoData1->Sumw2();
    histoData1->Scale(10.0/histoData1->Integral());
    histoData1->SetLineWidth(2);
    histoData1->SetMarkerStyle(20);
    histoData1->SetMarkerSize(1.5);
    histoData1->SetMarkerColor(kRed);
    histoData1->SetLineColor(kBlack);
    // TF1 *f1 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f1->SetParameter(0, 200);
    // f1->SetParameter(1, histoData1->GetMaximum()/10.);
    // f1->SetParameter(2, histoData1->GetMean());
    // f1->SetParameter(3, .05);
    // f1->SetParameter(4, .05);
    // f1->SetParName(0, "Offset");
    // f1->SetParName(1, "Normalization");
    // f1->SetParName(2, "Mean");
    // f1->SetParName(3, "Sigma");
    // f1->SetParName(4, "gamma");
    // f1->SetLineColor(kBlue);
    // f1->SetLineWidth(4);
    // histoData1->Fit(f1, "ME");
    histoData1->Fit("gaus", "0");

    histoData2->Sumw2();
    histoData2->Scale(10.0/histoData2->Integral());
    histoData2->SetLineWidth(2);
    histoData2->SetMarkerStyle(21);
    histoData2->SetMarkerSize(1.5);
    histoData2->SetMarkerColor(kBlue);
    histoData2->SetLineColor(kBlack);
    // TF1 *f2 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f2->SetParameter(0, 200);
    // f2->SetParameter(1, histoData1->GetMaximum()/10.);
    // f2->SetParameter(2, histoData1->GetMean());
    // f2->SetParameter(3, .05);
    // f2->SetParameter(4, .05);
    // f2->SetParName(0, "Offset");
    // f2->SetParName(1, "Normalization");
    // f2->SetParName(2, "Mean");
    // f2->SetParName(3, "Sigma");
    // f2->SetParName(4, "gamma");
    // f2->SetLineColor(kBlue);
    // f2->SetLineWidth(4);
    // histoData2->Fit(f2, "ME");
    histoData2->Fit("gaus", "0");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.15,0.6,0.385,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData1, "Hydrogen","lep");
    leg->AddEntry(histoData2, "Deuterium","lep");
    // leg->AddEntry(f, "Voigt fit","l");

    Double_t max = histoData1->GetMaximum();
    if(histoData2->GetMaximum() > max) max = histoData2->GetMaximum();
    max*=1.2;
    histoData1->SetMaximum(max);
    histoData1->SetMinimum(0.);
    histoData2->SetMaximum(max);
    histoData2->SetMinimum(0.);
    
    histoData1->Draw("E1");
    histoData2->Draw("same E1");
    leg->Draw("same");

    // add axis labels
    histoData1->GetXaxis()->SetTitle(xTitle);
    histoData1->GetXaxis()->CenterTitle();
    histoData1->GetXaxis()->SetTitleSize(0.055);
    histoData1->GetXaxis()->SetTitleOffset(0.90);
    histoData1->GetXaxis()->SetLabelOffset(0.010);
    histoData1->GetXaxis()->SetLabelSize(0.05);
    histoData1->GetYaxis()->SetTitle(yTitle);
    histoData1->GetYaxis()->CenterTitle();
    histoData1->GetYaxis()->SetTitleSize(0.055);
    histoData1->GetYaxis()->SetTitleOffset(1.0);
    histoData1->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData1->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");


}

void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    gStyle->SetFitFormat("4.3f");

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    // TF1 *f = new TF1(graphName, "[0]+[1]*TMath::Power(TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])),0.5)");
    // TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    // f->SetParameter(0, 250);
    // f->SetParameter(1, 100);
    // f->SetParameter(2, 0.395);
    // f->SetParameter(3, 1);
    // f->SetParameter(4, 10000);
    // f->SetParameter(5, 100);
    // f->SetParName(0, "Offset");
    // f->SetParName(1, "Peak intensity");
    // f->SetParName(2, "Aperture");
    // f->SetParName(3, "Center value");
    // f->SetParName(4, "Amp");
    // f->SetParLimits(0, 100, 300);
    // f->SetParLimits(1, 10, 100);
    // f->SetParLimits(2, 0.1, 1);
    // f->SetParLimits(3, -0.2, 2);
    // f->SetParLimits(4, 50, 10000);
    // f->SetParLimits(5, 1, 10000);
    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    TF1 *f1 = graph->GetFunction("fitRight");
    // TF1 *f2 = graph->GetFunction("linFit2");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    leg->AddEntry(f1, "Fit","l");
    // leg->AddEntry(f2, "Anti-Stokes linear fit","l");
    
    graph->Draw("ap");
    leg->Draw("same");

    // add axis labels
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetYaxis()->SetTitle(yTitle);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(0.9);
    graph->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    graph->SetTitle(graphName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}

void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph1 || !graph2){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph1->SetLineWidth(2);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(1.5);
    graph1->SetMarkerColor(kRed);
    graph1->SetLineColor(kBlack);
    graph2->SetLineWidth(2);
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(1.5);
    graph2->SetMarkerColor(kRed);
    graph2->SetLineColor(kBlack);

    TF1 *f1 = graph1->GetFunction("fitSmall");
    TF1 *f2 = graph2->GetFunction("fitLarge");

    f1->SetLineColor(kBlue);
    f1->SetLineWidth(4);
    f2->SetLineColor(kGreen+1);
    f2->SetLineWidth(4);
    f2->SetLineStyle(9);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph1, "Data","lep");
    leg->AddEntry(f1, "Small Stokes fit","l");
    leg->AddEntry(f2, "Large Stokes fit","l");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Nitrogen");

    mg->Add(graph1);
    mg->Add(graph2);
    
    // graph1->Draw("ap");
    // graph2->Draw("same ap");
    mg->Draw("ap");
    leg->Draw("same");

    // add axis labels
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitle(yTitle);
    mg->GetYaxis()->CenterTitle();
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(0.9);
    mg->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    mg->SetTitle(graphName);
    // can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}
