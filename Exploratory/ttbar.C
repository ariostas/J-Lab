#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
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
#include <TGraph.h>

#include <TLorentzVector.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type);
void histogram(vector<TH1D*> histos, vector<TString> histNames, TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);

/*
 * MAIN FUNCTION
 */

void ttbar(TString inputFile = "root://cmsxrootd.fnal.gov//store/user/rbi/merged/HIEWQExo-HIRun2015-PromptReco-v1-forest_v2/0.root"){

    TH1D *leadingMuPt = new TH1D("leadingMuPt", "leadingMuPt", 100, 0, 80);
    TH1D *leadingElePt = new TH1D("leadingElePt", "leadingElePt", 100, 0, 80);
    TH1D *leadingJetPt = new TH1D("leadingJetPt", "leadingJetPt", 100, 0, 80);
    
    TFile* infile = TFile::Open(inputFile);

    if(!infile){cout << "File not opened correctly. Exiting..." << endl; return;}
    
    TTree* leptons = (TTree*) infile->Get("ggHiNtuplizer/EventTree");
    TTree* jets = (TTree*) infile->Get("akPu4PFJetAnalyzer/t");

    if(!leptons || !jets){cout << "Trees not found. Tree: " << leptons << endl; return;}

    Long64_t totalEntries = leptons->GetEntries();

    vector<float> *muPt = 0, *elePt = 0;
    Int_t *nref = 0;
    Float_t jtpt[10000];

    leptons->SetBranchAddress("muPt",  &muPt);
    leptons->SetBranchAddress("elePt",  &elePt);
    jets->SetBranchAddress("nref",  &nref);
    jets->SetBranchAddress("jtpt",  &jtpt);
    // TClonesArray *jtpt = new TClonesArray("Float_t", 100);
    // jets->SetBranchAddress("jtpt", &jtpt);

    for(int entry = 0; entry < totalEntries; entry++){

        if(entry%10000 == 0) cout << entry << "/" << totalEntries << endl;

        leptons->GetEntry(entry);
        jets->GetEntry(entry);

        if(muPt->size() > 0) leadingMuPt->Fill(muPt->at(0));
        if(elePt->size() > 0) leadingElePt->Fill(elePt->at(0));
        //if(jtpt->size() > 0) 
        leadingJetPt->Fill(jtpt[0]);

        // if(jtpt->GetEntries() > 0) leadingJetPt->Fill(*(Float_t*) (*jtpt)[0]);

    }

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    histogram(leadingMuPt, "Leading #mu p_{T}", can, "Transverse momentum (GeV)", "Count", "muPt", "dot");
    histogram(leadingElePt, "Leading e^{-} p_{T}", can, "Transverse momentum (GeV)", "Count", "elePt", "dot");
    histogram(leadingJetPt, "Leading jet p_{T}", can, "Transverse momentum (GeV)", "Count", "jetPt", "dot");

    vector<TH1D*> histos = {leadingElePt, leadingMuPt, leadingJetPt};
    vector<TString> names = {"Electrons", "Muons", "Jets"};
    histogram(histos, names, "Leading objects p_{T}", can, "Transverse momentum (GeV)", "Normalized count", "leadingPt");

}

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

void histogram(vector<TH1D*> histos, vector<TString> histNames, TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){   
    
    if(histos.size()!=histNames.size()){ cout << "Number of histograms and names don't match." << endl; return;}

    Double_t max=0, min=1;

    for(Int_t i=0; i<histos.size(); i++){
        Double_t integral = histos.at(i)->Integral();
        if(integral != 0.){
            histos.at(i)->Scale(Double_t(1)/integral);
            if(histos.at(i)->GetMaximum()>max) max=histos.at(i)->GetMaximum();
            if(histos.at(i)->GetMinimum()<min) min=histos.at(i)->GetMinimum();
        }
    }

    max*=1.25;
    min*=0.75;
    //if(name == "Theta"){max=.1; min=0;}

    vector<Int_t> colors, lines;
    colors.push_back(kBlue); colors.push_back(kRed); colors.push_back(kGreen+2); colors.push_back(kBlack); colors.push_back(kRed); colors.push_back(kBlack); 
    lines.push_back(1); lines.push_back(9); lines.push_back(10); lines.push_back(7); lines.push_back(9);

    for(Int_t i=0; i<histos.size(); i++){
        if(histNames.at(i) == "Backgrounds") continue;
        histos.at(i)->SetMaximum(max);
        histos.at(i)->SetMinimum(min);
        histos.at(i)->SetLineWidth(4);
        histos.at(i)->SetLineColor(colors.at(i%6));
        histos.at(i)->SetLineStyle(lines.at(i%5));
        if(i==0){
            histos.at(i)->Draw("][");
        }
        else histos.at(i)->Draw("same ][");
        
    }

    histos.at(0)->GetXaxis()->SetTitle(xTitle);
    histos.at(0)->GetXaxis()->CenterTitle();
    histos.at(0)->GetXaxis()->SetTitleSize(0.055);
    histos.at(0)->GetXaxis()->SetTitleOffset(0.87);
    histos.at(0)->GetXaxis()->SetLabelOffset(0.010);
    histos.at(0)->GetXaxis()->SetLabelSize(0.05);
    histos.at(0)->GetYaxis()->SetTitle(yTitle);
    histos.at(0)->GetYaxis()->CenterTitle();
    histos.at(0)->GetYaxis()->SetTitleSize(0.055);
    histos.at(0)->GetYaxis()->SetTitleOffset(0.95);
    histos.at(0)->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histos.at(0)->SetTitle(histName);
    can->Update();

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.70,0.675,0.90,0.875);
    //leg->SetTextFont(72);
    leg->SetTextSize(0.04);
    leg->SetFillColor(kWhite);
    for(Int_t i=0; i<histos.size(); i++){
        if(histNames.at(i) == "Backgrounds") continue;
        leg->AddEntry(histos.at(i),histNames.at(i),"l");
    }
    leg->Draw("same");

    TString filename=name; filename+=".png";

    can->SaveAs(filename);
}
