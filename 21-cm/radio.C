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
#include <TF2.h>
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
#include <TGraphErrors.h>
#include <TProfile.h>
#include "TTimer.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void graph(TGraphErrors*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(TH2D*, const TString, TCanvas*, const TString, const TString, const TString, const TString, Bool_t anim = kFALSE);
TGraphErrors* timeScan(TString, TString);
TGraphErrors* linScan(TString, Double_t, Double_t, TString);
void npScan(TH2D*, TString, Double_t, Double_t);
Bool_t badData(TProfile*);
Double_t getAverage(TProfile*);
Double_t getError(TProfile*);

UInt_t h = 0;

/*
 * MAIN FUNCTION
 */

void radio(){

    TH1::StatOverflows(kTRUE);

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);

    TH2D *npscan1 = new TH2D("npscan1", "npscan1", 11, -5.5, 5.5, 11, -5.5, 5.5);
    TH2D *npscan2 = new TH2D("npscan2", "npscan2", 11, -5.5, 5.5, 11, -5.5, 5.5);
    TH2D *npscan3 = new TH2D("npscan3", "npscan3", 11, -5.5, 5.5, 11, -5.5, 5.5);
    TH2D *npscan4 = new TH2D("npscan4", "npscan4", 11, -5.5, 5.5, 11, -5.5, 5.5);
    TH2D *npscan5 = new TH2D("npscan5", "npscan5", 11, -11, 11, 11, -11, 11);
    TH2D *npscan6 = new TH2D("npscan6", "npscan6", 11, -11, 11, 11, -11, 11);
    TH2D *npscan7 = new TH2D("npscan7", "npscan7", 11, -11, 11, 11, -11, 11);
    TH2D *npscan8 = new TH2D("npscan8", "npscan8", 11, -16.5, 16.5, 16, -16.5, 16.5);

    // TGraphErrors *azscan1 = linScan("azscan_nov8_1", 1, 30, "az");
    // TGraphErrors *azscan2 = linScan("azscan_nov8_2", 1, 30, "az");
    //TGraphErrors *azscan3 = linScan("azscan_nov9_1", 1, 30, "az");
    // TGraphErrors *azscan4 = linScan("azscan_nov10_1", 0.5, 20, "az");
    // TGraphErrors *azscan5 = linScan("azscan_nov11_1", 0.5, 20, "az");
    // TGraphErrors *azscan6 = linScan("azscan_nov11_2", 0.5, 20, "az");

    // TGraphErrors *elscan1 = linScan("elscan_nov8_1", 1, 30, "el");
    // TGraphErrors *elscan2 = linScan("elscan_nov9_1", 1, 30, "el");
    // TGraphErrors *elscan3 = linScan("elscan_nov11_1", 0.5, 20, "el");
    //TGraphErrors *elscan4 = linScan("elscan_nov11_2", 0.5, 20, "el");

    // TGraphErrors *sunDrift1 = timeScan("sundrift_nov8_1", "Sun");
    // TGraphErrors *sunDrift2 = timeScan("sundrift_nov8_2", "Sun");
    // TGraphErrors *sunDrift3 = timeScan("sundrift_nov8_3", "Sun");
    TGraphErrors *sunDrift4 = timeScan("sundrift_nov9_1", "Sun");
    // TGraphErrors *sunDrift5 = timeScan("sundrift_nov11_1", "Sun");

    // npScan(npscan1, "npscan_1_nov9_1", 5, 1);
    // npScan(npscan2, "npscan_1_nov10_1", 5, 1);
    npScan(npscan3, "npscan_1_nov11_1", 5, 1);
    // npScan(npscan4, "npscan_1_nov11_2", 5, 1);
    // npScan(npscan5, "npscan_2_nov8_1", 10, 2);
    // npScan(npscan6, "npscan_2_nov10_1", 10, 2);
    npScan(npscan7, "npscan_2_nov11_1", 10, 2);
    // npScan(npscan8, "npscan_3_nov8_1", 15, 3);

    // graph(azscan1, "Azimuthal scan 1", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan1");
    // graph(azscan2, "Azimuthal scan 2", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan2");
    //graph(azscan3, "Azimuthal scan", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan");
    // graph(azscan4, "Azimuthal scan 4", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan4");
    // graph(azscan5, "Azimuthal scan 5", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan5");
    // graph(azscan6, "Azimuthal scan 6", can, "Offset (degrees)", "Average brightness temperature (K)", "azscan6");

    // graph(elscan1, "Elevation scan 1", can, "Offset (degrees)", "Average brightness temperature (K)", "elscan1");
    // graph(elscan2, "Elevation scan 2", can, "Offset (degrees)", "Average brightness temperature (K)", "elscan2");
    // graph(elscan3, "Elevation scan 3", can, "Offset (degrees)", "Average brightness temperature (K)", "elscan3");
    //graph(elscan4, "Elevation scan", can, "Offset (degrees)", "Average brightness temperature (K)", "elscan");

    // graph(sunDrift1, "Sun drift 1", can, "Time (min)", "Average brightness temperature (K)", "sundrift1");
    // graph(sunDrift2, "Sun drift 2", can, "Time (min)", "Average brightness temperature (K)", "sundrift2");
    // graph(sunDrift3, "Sun drift 3", can, "Time (min)", "Average brightness temperature (K)", "sundrift3");
    graph(sunDrift4, "Solar drift (Corrected)", can, "Time (min)", "Average brightness temperature (K)", "sundrift");
    // graph(sunDrift5, "Sun drift 5", can, "Time (min)", "Average brightness temperature (K)", "sundrift5");

    // histogram(npscan1, "Sun npoint scan 1", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan1");
    // histogram(npscan2, "Sun npoint scan 2", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan2");
    //histogram(npscan3, "Solar scan (1 deg)", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "Brightness temperature (K)", "npscan1");
    // histogram(npscan4, "Sun npoint scan 4", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan4");
    // histogram(npscan5, "Sun npoint scan 5", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan5");
    // histogram(npscan6, "Sun npoint scan 6", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan6");
    //histogram(npscan7, "Solar scan (2 deg)", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "Brightness temperature (K)", "npscan2");
    // histogram(npscan8, "Sun npoint scan 8", can, "Azimuthal offset (deg)", "Elevation offset (deg)", "npscan8");

}

TGraphErrors* linScan(TString inputfile, Double_t step, Double_t maxOffset, TString type){

    if(type != "az" && type != "el"){
        cout << "Error: Type of linear scan \"" << type << "\" not defined" << endl;
        return NULL;
    }

    TFile file("./Data/" + inputfile + ".root");
    if(!file.IsOpen()){
        cout << "Error: File \"" << inputfile << ".root\" does not exist" << endl;
        return NULL;
    }

    Double_t offset[100], intLumi[100], offsetError[100], intLumiError[100];
    Int_t badPoints = 0;

    for(Double_t x = -maxOffset; x <= maxOffset; x += step){
        TProfile *p = (TProfile*) file.Get((TString::Format( (type == "az" ? "offset %1.1f 0" : "offset 0 %1.1f"), x)).Data());
        if(!p) p = (TProfile*) file.Get((TString::Format( (type == "az" ? "offset %1i 0" : "offset 0 %1i"), Nint(x))).Data());
        if(!p){
            cout << "Error: File seems to be corrupted or the scan parameters are wrong: " << (TString::Format( (type == "az" ? "offset %1i 0" : "offset 0 %1i"), x)).Data() << " not found" << endl;
            return NULL;
        }
        if(getAverage(p) == 0.){
            badPoints++;
            continue;
        }
        offset[Nint(x+maxOffset-badPoints)] = x;
        offsetError[Nint(x+maxOffset-badPoints)] = 0.1;
        intLumi[Nint(x+maxOffset-badPoints)] = getAverage(p);
        intLumiError[Nint(x+maxOffset-badPoints)] = getError(p);

        //intLumi[x+maxOffset] = p->GetMaximum();
        //intLumi[x+maxOffset] = p->GetBinContent(p->GetXaxis()->GetNbins()/2+1);
    }

    file.Close();

    TGraphErrors *graph = new TGraphErrors(2*maxOffset+1-badPoints, offset, intLumi, offsetError, intLumiError);

    return graph;

}

TGraphErrors* timeScan(TString inputfile, TString scanName){

    TFile file("./Data/" + inputfile + ".root");
    if(!file.IsOpen()){
        cout << "Error: File \"" << inputfile << ".root\" does not exist" << endl;
        return NULL;
    }

    Double_t t[1000], intLumi[1000], tError[1000], intLumiError[1000];
    Int_t nEntries=0, badPoints=0;

    for(Int_t x = 0; x+badPoints <= 1000; x++){
        TProfile *p = (TProfile*) file.Get((TString::Format(" %s_%1i", scanName.Data(), x+badPoints)).Data());
        if(!p && x == 0){
            cout << "Error: File seems to be corrupted or the scan parameters are wrong: " << (TString::Format("%s_%1i", scanName.Data(), x)).Data() << " not found" << endl;
            return NULL;
        }
        else if(!p){
            nEntries = x;
            break;
        }
        if(badData(p)){
            x--;
            badPoints++;
            continue;
        }
        t[x] = (x+badPoints)*7./60.;
        tError[x] = 0.5/60.;
        intLumi[x] = getAverage(p);
        intLumiError[x] = getError(p);
        //intLumi[x] = p->GetMaximum();
        //intLumi[x] = p->GetBinContent(p->GetXaxis()->GetNbins()/2+1);
    }

    file.Close();

    TGraphErrors *graph = new TGraphErrors(nEntries, t, intLumi, tError, intLumiError);

    return graph;

}

void npScan(TH2D *histo, TString inputfile, Double_t maxOffset, Double_t step){

    TFile file("./Data/" + inputfile + ".root");
    if(!file.IsOpen()){
        cout << "Error: File \"" << inputfile << ".root\" does not exist" << endl;
        return;
    }

    //TH2D *histo = new TH2D(inputfile, inputfile, 2*maxOffset/step + 1, -maxOffset-step/2., maxOffset+step/2., 2*maxOffset/step + 1, -maxOffset-step/2., maxOffset+step/2.);
    //TH2D *histo1 = new TH2D("test", "test", 11, -20, 20, 11, -20, 20);

    for(Int_t x = -maxOffset; x <= maxOffset; x+=step){
        for(Int_t y = -maxOffset; y <= maxOffset; y+=step){
            TProfile *p = (TProfile*) file.Get((TString::Format("offset %1i %1i", x, y)).Data());
            if(!p){
                cout << "Error: File seems to be corrupted or the scan parameters are wrong: " << (TString::Format("offset %1i %1i", x, y)).Data() << " not found" << endl;
                return;
            }
            histo->Fill(x, y, getAverage(p));
            //histo->Fill(x, y, p->GetMaximum());
            //histo->Fill(x, y, p->GetBinContent(p->GetXaxis()->GetNbins()/2+1));
        }
    }

    file.Close();

}

void graph(TGraphErrors *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    gStyle->SetOptStat(2210);
    //gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1100);

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    // TF1 *f = new TF1(graphName, "[0]+[1]*TMath::Power(TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])),0.5)");
    TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    f->SetParameter(0, 250);
    f->SetParameter(1, 100);
    f->SetParameter(2, 0.395);
    f->SetParameter(3, 1);
    f->SetParameter(4, 10000);
    f->SetParameter(5, 100);
    f->SetParName(0, "Offset");
    f->SetParName(1, "Peak intensity");
    f->SetParName(2, "Aperture");
    f->SetParName(3, "Center value");
    f->SetParName(4, "Amp");
    f->SetParLimits(0, 100, 300);
    f->SetParLimits(1, 10, 100);
    f->SetParLimits(2, 0.1, 1);
    f->SetParLimits(3, -0.2, 2);
    f->SetParLimits(4, 50, 10000);
    f->SetParLimits(5, 1, 10000);
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    //leg->AddEntry(f, "Airy disk fit","l");
    
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

void histogram(TH2D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString zTitle, const TString name, Bool_t anim){

    gStyle->SetOptStat(2210);
    //gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);

    if(!histoData){
        cout << "Error: Histogram \"" << histName << "\" not defined" << endl;
        return;
    }

    TF2 *f2 = new TF2("func","[0]+[1]*TMath::Gaus(x,[2],[3])*TMath::Gaus(y,[4],[5])");
    f2->SetParameter(0, 200);
    f2->SetParameter(1, 300);
    f2->SetParameter(2, 2);
    f2->SetParameter(3, 4);
    f2->SetParameter(4, 1);
    f2->SetParameter(5, 4);
    f2->SetParName(0, "Constant");
    f2->SetParName(1, "Amplitude");
    f2->SetParName(2, "Offset az");
    f2->SetParName(3, "Sigma az");
    f2->SetParName(4, "Offset el");
    f2->SetParName(5, "Sigma el");
    f2->SetLineColor(kBlack);
    f2->SetLineWidth(1);
    histoData->Fit(f2, "ME");

    //cout << histoData << endl;
    //gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.755,0.87, 0.999, 0.999);
    leg->SetTextSize(0.055);
    leg->AddEntry(f2, "Gaussian fit","l");
    
    histoData->Draw("lego2");
    leg->Draw("same");
    

    //add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(1.25);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(1.25);
    histoData->GetYaxis()->SetLabelSize(0.05);
    histoData->GetZaxis()->SetTitle(zTitle);
    histoData->GetZaxis()->CenterTitle();
    histoData->GetZaxis()->SetTitleSize(0.055);
    histoData->GetZaxis()->SetTitleOffset(0.9);
    histoData->GetZaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");

    if(anim){
        Double_t phi = 30;
        gStyle->SetCanvasPreferGL(true);
        //gStyle->SetFrameFillColor(42);
        Double_t steps = 360;
        for(Int_t x = 1; x <= steps; x++){
            phi += 360./steps;
            gPad->SetPhi(phi);
            gPad->Modified();
            gPad->Update();
            if(x == steps) gPad->Print("testanim.gif++3++");
            else gPad->Print(name + ".gif+3");
            cout << "saved " << x << " of " << steps << endl;
        }
    }
}

Double_t getAverage(TProfile *p)
{
    if(!p) return 0;

    Double_t total = 0;

    for(UInt_t bin = 11; bin <= p->GetXaxis()->GetNbins()-2-10; bin++){
        total += p->GetBinContent(bin);
    }

    Double_t average = total/Double_t(p->GetXaxis()->GetNbins()-2-20);

    return average;
}

Double_t getError(TProfile *p)
{
    if(!p) return 0;

    Double_t total = 0;

    for(UInt_t bin = 11; bin <= p->GetXaxis()->GetNbins()-2-10; bin++){
        total += p->GetBinError(bin);
    }

    Double_t error = total/Double_t(p->GetXaxis()->GetNbins()-2-20);

    return error;
}


/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    // histoData->SetLineWidth(2);
    // histoData->SetMarkerStyle(20);
    // histoData->SetMarkerSize(1.5);
    // histoData->SetMarkerColor(kRed);
    // histoData->SetLineColor(kBlack);
    // TF1 *f = new TF1(histName, "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f->SetParameter(0, 200);
    // f->SetParameter(1, histoData->GetMaximum()/10.);
    // f->SetParameter(2, histoData->GetMean());
    // f->SetParameter(3, .05);
    // f->SetParameter(4, .05);
    // f->SetParName(0, "Offset");
    // f->SetParName(1, "Normalization");
    // f->SetParName(2, "Mean");
    // f->SetParName(3, "Sigma");
    // f->SetParName(4, "gamma");
    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    // histoData->Fit(f, "ME");
    //TF1 *f = histoData->GetFunction("fit");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data","lep");
    //leg->AddEntry(f, "Double Voigtian fit","l");
    
    histoData->Draw("E1");
    leg->Draw("same");

    // // add axis labels
    // histoData->GetXaxis()->SetTitle(xTitle);
    // histoData->GetXaxis()->CenterTitle();
    // histoData->GetXaxis()->SetTitleSize(0.055);
    // histoData->GetXaxis()->SetTitleOffset(0.90);
    // histoData->GetXaxis()->SetLabelOffset(0.010);
    // histoData->GetXaxis()->SetLabelSize(0.05);
    // histoData->GetYaxis()->SetTitle(yTitle);
    // histoData->GetYaxis()->CenterTitle();
    // histoData->GetYaxis()->SetTitleSize(0.055);
    // histoData->GetYaxis()->SetTitleOffset(1.0);
    // histoData->GetYaxis()->SetLabelSize(0.05);
    // gStyle->SetTitleSize(0.08, "t");
    // histoData->SetTitle(histName);
    // can->Update();

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

Bool_t badData(TProfile *p){

    if(!p) return kTRUE;

    h++;

    TH1D *histo = new TH1D((TString::Format("tempHisto_%1u", h)).Data(), (TString::Format("tempHisto_%1u", h)).Data(), 1000,
        p->GetBinContent(p->GetMinimumBin()), p->GetBinContent(p->GetMaximumBin()));

    for(UInt_t bin = 11; bin <= p->GetXaxis()->GetNbins()-2-10; bin++){
        histo->Fill(p->GetBinContent(bin));
    }

    Double_t stdDev = histo->GetStdDev();

    if(stdDev > 200) return kTRUE;

    return kFALSE;
}
