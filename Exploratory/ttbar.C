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
void histogram(vector<TH1D*> histosData, vector<TH1D*> histosMC, vector<TString> histNames, TCanvas *can, vector<TString> xTitle, vector<TString> yTitle, vector<TString> names);
vector<TH1D*> analize(TString input = "data");

/*
 * MAIN FUNCTION
 */

 void ttbar(){

    TH1::StatOverflows(kTRUE);
    gStyle->SetOptStat(kFALSE);

    vector<TH1D*> histosData = analize("data"), histosMC = analize("MC");
    vector<TString> histNames = {"Leading muons", "Leading electrons", "Leading muons", "Leading electrons", "Leading muons", "Leading electrons", "Dilepton system", "Dilepton system", "Dilepton system", "Missing transverse energy"};
    vector<TString> xTitles = {"Transverse momentum (GeV)", "Transverse momentum (GeV)", "#eta", "#eta", "#phi", "#phi", "Mass (GeV)", "Transverse momentum (GeV)", "#DeltaR", "MET (GeV)"};
    vector<TString> yTitles = {"Normalized count", "Normalized count", "Normalized count", "Normalized count", "Normalized count","Normalized count", "Normalized count", "Normalized count", "Normalized count", "Normalized count"};
    vector<TString> filenames = {"leadPt", "subleadPt", "leadEta", "subleadEta", "leadPhi", "subleadPhi", "dilepMass", "dilepPt", "dilepDR", "met"};

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);

    histogram(histosData, histosMC, histNames, can, xTitles, yTitles, filenames);

    // histogram(leadingMuPt, "Leading #mu p_{T}", can, "Transverse momentum (GeV)", "Count", "muPt", "dot");
    // histogram(leadingElePt, "Leading e^{-} p_{T}", can, "Transverse momentum (GeV)", "Count", "elePt", "dot");
    // histogram(leadingJetPt, "Leading jet p_{T}", can, "Transverse momentum (GeV)", "Count", "jetPt", "dot");

    // histogram(dileptonPt, "Dilepton system", can, "Transverse momentum (GeV)", "Count", "dilepPt", "line");
    // histogram(dileptonMass, "Dilepton system", can, "Mass (GeV)", "Count", "dilepMass", "line");

    // vector<TH1D*> histos = {leadingElePt, leadingMuPt, leadingJetPt};
    // vector<TString> names = {"Electrons", "Muons", "Jets"};
    // histogram(histos, names, "Leading objects p_{T}", can, "Transverse momentum (GeV)", "Normalized count", "leadingPt");

 }

vector<TH1D*> analize(TString input){

    TH1D *leadingLepPt = new TH1D("leadingLepPt_"+input, "leadingLepPt_"+input, 100, 20, 140);
    TH1D *subleadingLepPt = new TH1D("subleadingLepPt_"+input, "subleadingLepPt_"+input, 100, 20, 140);
    TH1D *leadingLepEta = new TH1D("leadingLepEta_"+input, "leadingLepEta_"+input, 100, -2.5, 2.5);
    TH1D *subleadingLepEta = new TH1D("subleadingLepEta_"+input, "subleadingLepEta_"+input, 100, -2.5, 2.5);
    TH1D *leadingLepPhi = new TH1D("leadingLepPhi_"+input, "leadingLepPhi_"+input, 100, -3.1416, 3.1416);
    TH1D *subleadingLepPhi = new TH1D("subleadingLepPhi_"+input, "subleadingLepPhi_"+input, 100, -3.1416, 2.1416);
    TH1D *subleadingLepIso = new TH1D("subleadingLepIso_"+input, "subleadingLepIso_"+input, 100, 0, 5);

    TH1D *dileptonMass = new TH1D("dileptonMass_"+input, "dileptonMass_"+input, 100, 0, 300);
    TH1D *dileptonPt = new TH1D("dileptonPt_"+input, "dileptonPt_"+input, 100, 0, 300);
    TH1D *dileptonDR = new TH1D("dileptonDR_"+input, "dileptonDR_"+input, 100, 0, 6);
    TH1D *met = new TH1D("met_"+input, "met_"+input, 100, 0, 100);

    vector<TH1D*> histograms = {leadingLepPt, subleadingLepPt, leadingLepEta, subleadingLepEta, leadingLepPhi, subleadingLepPhi, dileptonMass, dileptonPt, dileptonDR, met};

    TString inputFile;
    // if(input == "data") inputFile = "root://cmsxrootd.fnal.gov//store/user/rbi/merged/HIEWQExo-HIRun2015-PromptReco-v1-forest_v2/0.root";
    if(input == "data") inputFile = "root://cmsxrootd.fnal.gov//store/user/rbi/merged/rbi-HIEWQExo-HIRun2015-PromptReco-emu-skim-forest/0.root";
    else if(input == "MC") inputFile = "root://cmsxrootd.fnal.gov//store/user/rbi/merged/rbi-TTbar_5TeV_TuneCUETP8M1_75X_GEN_SIM-forest_2/0.root";
    else{cout << "Wrong input" << endl; return histograms;}
    
    TFile* infile = TFile::Open(inputFile);

    if(!infile){cout << "File not opened correctly. Exiting..." << endl; return histograms;}
    
    TTree* leptons = (TTree*) infile->Get("ggHiNtuplizer/EventTree");
    TTree* jets = (TTree*) infile->Get("akPu3PFJetAnalyzer/t");
    TTree* pftree = (TTree*) infile->Get("pfcandAnalyzer/pfTree");

    if(!leptons || !jets){cout << "Trees not found. Tree: " << leptons << endl; return histograms;}

    Long64_t totalEntries = leptons->GetEntries();

    vector<float> *muPt = 0, *muEta = 0, *muPhi = 0, *muCharge = 0, *muIsoTrk = 0;
    vector<int> *muIsGood = 0;
    vector<float> *elePt = 0, *eleEta = 0, *elePhi = 0, *eleCharge = 0, *elePFChIso04 = 0, *eleSCEta = 0, *eleBrem = 0, *eleSigmaIEtaIEta = 0;
    vector<float> *pfPt = 0, *pfEta = 0, *pfPhi = 0, *pfEnergy = 0;
    Int_t nref;
    Float_t jtpt[1000], jteta[1000], jtphi[1000], jtm[1000];
    Float_t muMass = 0.105658, eleMass = 0.000510998;

    leptons->SetBranchAddress("muPt",  &muPt);
    leptons->SetBranchAddress("muEta",  &muEta);
    leptons->SetBranchAddress("muPhi",  &muPhi);
    leptons->SetBranchAddress("muCharge",  &muCharge);
    leptons->SetBranchAddress("muIsoTrk",  &muIsoTrk);
    leptons->SetBranchAddress("muIsGood",  &muIsGood);
    leptons->SetBranchAddress("elePt",  &elePt);
    leptons->SetBranchAddress("eleEta",  &eleEta);
    leptons->SetBranchAddress("elePhi",  &elePhi);
    leptons->SetBranchAddress("eleCharge",  &eleCharge);
    leptons->SetBranchAddress("elePFChIso04",  &elePFChIso04);
    leptons->SetBranchAddress("eleSCEta",  &eleSCEta);
    leptons->SetBranchAddress("eleBrem",  &eleBrem);
    leptons->SetBranchAddress("eleSigmaIEtaIEta",  &eleSigmaIEtaIEta);
    jets->SetBranchAddress("nref",  &nref);
    jets->SetBranchAddress("jtpt",  &jtpt);
    jets->SetBranchAddress("jteta",  &jteta);
    jets->SetBranchAddress("jtphi",  &jtphi);
    jets->SetBranchAddress("jtm",  &jtm);
    pftree->SetBranchAddress("pfPt",  &pfPt);
    pftree->SetBranchAddress("pfEta",  &pfEta);
    pftree->SetBranchAddress("pfPhi",  &pfPhi);
    pftree->SetBranchAddress("pfEnergy",  &pfEnergy);

    TString out = TString::Format("Reading %s...   ", input.Data());
    out.Resize(18);

    Int_t eventsRemaining = 0, onlyMuons = 0, onlyElectrons = 0;
    cout << "Total events: " << totalEntries << endl;
    TLorentzVector tempP4, tempP42;

    for(int entry = 0; entry < totalEntries; entry++){

        if(entry == 0) cout << TString::Format("%s%2.0f\% \n", out.Data(), 0.);
        else if(entry%1000 == 0) cout << TString::Format("\e[A%s%2.0f\% \n", out.Data(), Double_t(entry)/Double_t(totalEntries)*100.);
        else if(entry == totalEntries-1) cout << TString::Format("\e[A%s\033[1;32mdone\033[0m", out.Data());

        leptons->GetEntry(entry);
        jets->GetEntry(entry);
        // pftree->GetEntry(entry);

        Int_t muon = -1, electron = -1, jet1 = -1, jet2 = -1, nMuons = 0, nElectrons = 0, nJets = 0;

        TLorentzVector metP4; metP4.SetPxPyPzE(0,0,0,0);
        // for(UInt_t x = 0; x < pfPt->size(); x++){
        //     tempP4.SetPtEtaPhiE(pfPt->at(x), pfEta->at(x),pfPhi->at(x), pfEnergy->at(x));
        //     metP4 -= tempP4;
        // }

        for(UInt_t x = 0; x < muPt->size(); x++){
            if(muPt->at(x) <= 18) continue;
            if(Abs(muEta->at(x)) >= 2.4) continue;
            if(muIsGood->at(x) == 0) continue;
            // if(muIsoTrk->at(x) / muPt->at(x) > 0.4) continue;

            tempP42.SetPtEtaPhiM(muPt->at(x), muEta->at(x), muPhi->at(x), muMass);

            bool isInJet = false;
            for(UInt_t y = 0; y < nref; y++){
                if(jtpt[y] < muPt->at(x)) continue;
                tempP4.SetPtEtaPhiM(jtpt[y], jteta[y], jtphi[y], jtm[y]);
                if(tempP4.DeltaR(tempP42) < 0.3){
                    isInJet = true;
                    // cout << "Muon eta: " << muEta->at(x) << " Jet eta: " << jteta[y] << endl;
                    // cout << "Muon phi: " << muPhi->at(x) << " Jet phi: " << jtphi[y] << endl;
                    // cout << "Muon pt: " << muPt->at(x) << " Jet pt: " << jtpt[y] << endl << endl;
                    break;
                }
            }
            // if(isInJet) continue;

            nMuons++;
            if(muon == -1) muon = x;
            else if(muPt->at(x) > muPt->at(muon)) muon = x;
        }

        for(UInt_t x = 0; x < elePt->size(); x++){
            if(elePt->at(x) <= 20) continue;
            if(Abs(eleEta->at(x)) >= 2.4) continue;
            if(eleBrem->at(x) <= -0.1 || eleBrem->at(x) >= 0.2) continue;
            // if(eleSigmaIEtaIEta->at(x) >= 0.024) continue;
            // if(Abs(eleSCEta->at(x)) > 1.4442) continue;
            // if(elePFChIso04->at(x) / elePt->at(x) > 0.4) continue;
            tempP42.SetPtEtaPhiM(elePt->at(x), eleEta->at(x), elePhi->at(x), eleMass);

            bool isInJet = false;
            for(UInt_t y = 0; y < nref; y++){
                if(jtpt[y] < elePt->at(x)) continue;
                tempP4.SetPtEtaPhiM(jtpt[y], jteta[y], jtphi[y], jtm[y]);
                if(tempP4.DeltaR(tempP42) < 0.3){
                    isInJet = true;
                    // cout << "Electron eta: " << eleEta->at(x) << " Jet eta: " << jteta[y] << endl;
                    // cout << "Electron phi: " << elePhi->at(x) << " Jet phi: " << jtphi[y] << endl;
                    // cout << "Electron pt: " << elePt->at(x) << " Jet pt: " << jtpt[y] << endl << endl;
                    break;
                }
            }
            // if(isInJet) continue;

            nElectrons++;
            if(electron == -1) electron = x;
            else if(elePt->at(x) > elePt->at(electron)) electron = x;
        }

        Double_t leadPt, subleadPt, leadEta, subleadEta, leadPhi, subleadPhi;
        if(nMuons > 0){
            leadPt = muPt->at(muon);
            leadEta = muEta->at(muon);
            leadPhi = muPhi->at(muon);
            onlyMuons++;
            leadingLepPt->Fill(leadPt);
            leadingLepEta->Fill(leadEta);
            leadingLepPhi->Fill(leadPhi);
        }
        if(nElectrons > 0){
            subleadPt = elePt->at(electron);
            subleadEta = eleEta->at(electron);
            subleadPhi = elePhi->at(electron);
            onlyElectrons++;
            subleadingLepPt->Fill(subleadPt);
            subleadingLepEta->Fill(subleadEta);
            subleadingLepPhi->Fill(subleadPhi);
            subleadingLepIso->Fill(elePFChIso04->at(electron) / elePt->at(electron));
        }

        if(nMuons < 1 || nElectrons < 1) continue;
        if(muCharge->at(muon) == eleCharge->at(electron)) continue;

        // for(Int_t x = 0; x < nref; x++){
        //     if(jtpt[x] < 20) continue;

        //     nJets++;
        //     if(jet1 == -1) jet1 = x;
        //     else if(jtpt[x] > jtpt[jet1]) jet1 = x;
        //     else if(jet2 == -1) jet2 = x;
        //     else if(jtpt[x] > jtpt[jet2]) jet2 = x;
        // }

        // if(nJets < 2) continue;

        TLorentzVector muP4; muP4.SetPtEtaPhiM(muPt->at(muon), muEta->at(muon), muPhi->at(muon), muMass);
        TLorentzVector eleP4; eleP4.SetPtEtaPhiM(elePt->at(electron), eleEta->at(electron), elePhi->at(electron), eleMass);

        TLorentzVector dilepP4 = muP4 + eleP4;

        // Double_t leadPt, subleadPt, leadEta, subleadEta, leadPhi, subleadPhi;
        // if(muPt->at(muon) > elePt->at(electron)){
            // leadPt = muPt->at(muon);
            // subleadPt = elePt->at(electron);
            // leadEta = muEta->at(muon);
            // subleadEta = eleEta->at(electron);
            // leadPhi = muPhi->at(muon);
            // subleadPhi = elePhi->at(electron);
        // }
        // else{
        //     subleadPt = muPt->at(muon);
        //     leadPt = elePt->at(electron);
        //     subleadEta = muEta->at(muon);
        //     leadEta = eleEta->at(electron);
        //     subleadPhi = muPhi->at(muon);
        //     leadPhi = elePhi->at(electron);
        // }

        // leadingLepPt->Fill(leadPt);
        // subleadingLepPt->Fill(subleadPt);
        // leadingLepEta->Fill(leadEta);
        // subleadingLepEta->Fill(subleadEta);
        // leadingLepPhi->Fill(leadPhi);
        // subleadingLepPhi->Fill(subleadPhi);

        dileptonPt->Fill(dilepP4.Pt());
        dileptonMass->Fill(dilepP4.M());
        dileptonDR->Fill(muP4.DeltaR(eleP4));
        met->Fill(metP4.Pt());


        eventsRemaining++;

    }

    cout << TString::Format("     Remaining events: %i. Events with a muon: %i. Events with an electron: %i\n", eventsRemaining, onlyMuons, onlyElectrons);

    return histograms;

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

    for(UInt_t i=0; i<histos.size(); i++){
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

    for(UInt_t i=0; i<histos.size(); i++){
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
    for(UInt_t i=0; i<histos.size(); i++){
        if(histNames.at(i) == "Backgrounds") continue;
        leg->AddEntry(histos.at(i),histNames.at(i),"l");
    }
    leg->Draw("same");

    TString filename=name; filename+=".png";

    can->SaveAs(filename);
}

void histogram(vector<TH1D*> histosData, vector<TH1D*> histosMC, vector<TString> histNames, TCanvas *can, vector<TString> xTitle, vector<TString> yTitle, vector<TString> names){   
    
    if(histosData.size()!=histNames.size() || histosData.size()!=histosMC.size()){ cout << "Number of histograms don't match." << endl; return;}

    for(UInt_t x = 0; x < histosData.size(); x++){

        histosData.at(x)->Scale(1.0/histosData.at(x)->Integral());
        histosMC.at(x)->Scale(1.0/histosMC.at(x)->Integral());

        Double_t max=-999, min=999;

        if(histosData.at(x)->GetMaximum() > max) max = histosData.at(x)->GetMaximum();
        if(histosMC.at(x)->GetMaximum() > max) max = histosMC.at(x)->GetMaximum();
        if(histosData.at(x)->GetMinimum() < min) min = histosData.at(x)->GetMinimum();
        if(histosMC.at(x)->GetMinimum() < min) min = histosMC.at(x)->GetMinimum();

        max*=1.25;
        min*=0.75;

        histosData.at(x)->SetMaximum(max);
        histosData.at(x)->SetMinimum(min);
        histosMC.at(x)->SetMaximum(max);
        histosMC.at(x)->SetMinimum(min);

        histosData.at(x)->SetLineWidth(3);
        histosMC.at(x)->SetLineWidth(3);

        histosData.at(x)->SetLineStyle(1);
        histosMC.at(x)->SetLineStyle(9);
        histosData.at(x)->SetLineColor(kBlue);
        histosMC.at(x)->SetLineColor(kRed);

        histosData.at(x)->Draw("hist");
        histosMC.at(x)->Draw("hist same");

        histosData.at(x)->GetXaxis()->SetTitle(xTitle.at(x));
        histosData.at(x)->GetXaxis()->CenterTitle();
        histosData.at(x)->GetXaxis()->SetTitleSize(0.055);
        histosData.at(x)->GetXaxis()->SetTitleOffset(0.87);
        histosData.at(x)->GetXaxis()->SetLabelOffset(0.010);
        histosData.at(x)->GetXaxis()->SetLabelSize(0.05);
        histosData.at(x)->GetYaxis()->SetTitle(yTitle.at(x));
        histosData.at(x)->GetYaxis()->CenterTitle();
        histosData.at(x)->GetYaxis()->SetTitleSize(0.055);
        histosData.at(x)->GetYaxis()->SetTitleOffset(0.95);
        histosData.at(x)->GetYaxis()->SetLabelSize(0.05);
        gStyle->SetTitleSize(0.08, "t");
        histosData.at(x)->SetTitle(histNames.at(x));
        can->Update();

        gStyle->SetLegendBorderSize(0);
        TLegend *leg = new TLegend(0.70,0.675,0.90,0.875);
        //leg->SetTextFont(72);
        leg->SetTextSize(0.04);
        leg->SetFillColor(kWhite);
        leg->AddEntry(histosData.at(x), "Data", "l");
        leg->AddEntry(histosMC.at(x), "MC", "l");
        leg->Draw("same");

        TString filename=names.at(x); filename+=".png";

        can->SaveAs(filename);
    }
}
