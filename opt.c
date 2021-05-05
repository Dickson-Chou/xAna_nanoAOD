// draw the graph after the selection
// example code to run 2016 NSSM MC X->Y+H samples
// .L xAna_nano_nssm.C++
// xAna_nano_nssm("test.root") or xAna_nano_nssm("input.txt")
// example root file is at /afs/cern.ch/work/s/syu/public/forTiKai/nssm_nano.root
//==========================================
//updated = 2020/04/28
//==========================================

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TLorentzVector.h>
//#include "setTDRStyle.C"

const double L2016 = 35.9*1000; // unit pb
// xs unit pb & number 

using namespace std;

//define punzi_eq
  	double punzi(double sigeff, double bg){
  		return sigeff/(1+TMath::Power(bg,0.5));
  	}


void opt(std::string outputFileName = "opt.root"){

//setTDRStyle();
gStyle->SetOptStat(0);
const int LINEWIDTH=3;
const int NXDIVISIONS=5;
const int NYDIVISIONS=6;
const double LATEXSIZE = 0.070;// x axis // 0.070
const double LABELSIZE = 37.5; //numbers
Float_t eta[24]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};//delta eta value
float pun[24]={0.00116173,0.00117028,0.00118103,0.00118863,0.00119593,0.00120309,0.00120745,0.0012083,0.00121738,0.00121898,0.00121791,0.00122125,0.00122419,0.00122309,0.00122211,0.00122424,0.00122289,0.00122305,0.00122273,0.00122078,0.00122441,0.00122331,0.00122301,0.00122264};
float M90[24]={0.00135201,0.00135252,0.00135876,0.00135637,0.00135984,0.00136463,0.00136149,0.00135972,0.00136724,0.00136694,0.001365,0.00136867,0.00137119,0.00136754,0.00136406,0.00136895,0.00136407,0.0013631,0.00136395,0.00135993,0.00136504,0.00136198,0.0013638,0.00136094};
float M100[24]={0.00140182,0.00140448,0.00140951,0.00141073,0.00141175,0.00141044,0.00141342,0.00141123,0.0014206,0.00142128,0.00141566,0.00141867,0.00142331,0.00142047,0.00142012,0.00142037,0.00141703,0.00141794,0.00141741,0.0014145,0.00141785,0.00141609,0.00141519,0.00141558};
float M125[24]={0.00135682,0.00135646,0.0013641,0.00136308,0.00136592,0.00136911,0.0013705,0.00136575,0.0013737,0.00137524,0.00137282,0.0013752,0.0013762,0.00137455,0.00137514,0.00137606,0.0013734,0.00137377,0.00137257,0.00137111,0.00137444,0.00137236,0.00137204,0.00137148};
float M150[24]={0.00126358,0.00126973,0.00127065,0.0012722,0.00127896,0.00128073,0.00128048,0.00127922,0.00128012,0.00128519,0.00128381,0.00128286,0.00128829,0.00128697,0.00128289,0.00128472,0.00128308,0.00128252,0.0012821,0.00128022,0.00128368,0.00128297,0.00128097,0.00128085};
float M200[24]={0.00129508,0.00129497,0.00130059,0.00130007,0.00130601,0.00130991,0.00130842,0.00130884,0.00131395,0.00131114,0.00131164,0.00131461,0.00131656,0.00131538,0.00131335,0.00131526,0.00131357,0.00131198,0.00131088,0.00130828,0.00131,0.00131325,0.00131016,0.00131142};
float M250[24]={0.00128209,0.00128072,0.00128635,0.00129174,0.00129508,0.00129468,0.00129738,0.00129468,0.00130555,0.00130267,0.00130082,0.00130466,0.00130551,0.00130288,0.00130215,0.00130365,0.00130385,0.00130373,0.00130317,0.00130231,0.00130698,0.00130812,0.00130743,0.00130566};
float D13[16]={0.00137119,0.00142331,0.0013762,0.00128829,0.00131656,0.00130551,0.00129219,0.0012782,0.00126358,0.00125187,0.00120751,0.00115938,0.00104727,0.000858562,0.000425336,0.000428151};
//open file
  TFile* signal = TFile::Open("signal.root");
  TFile* back_500_700 = TFile::Open("back_500_700.root");
  TFile* back_700_1000 = TFile::Open("back_700_1000.root");
  TFile* back_1000_1500 = TFile::Open("back_1000_1500.root");
  TFile* back_1500_2000 = TFile::Open("back_1500_2000.root");
  TFile* back_2000_inf = TFile::Open("back_2000_inf.root");
//  TFile* ttbar = TFile::Open("ttbar.root");


  TH1F* HBB = (TH1F*) signal->Get("signal150_ 1.3");
  
  TH1F* HBB_5 = (TH1F*) back_500_700->Get("signal_ 1.3");
  TH1F* HBB_7 = (TH1F*) back_700_1000->Get("signal_ 1.3");
  TH1F* HBB_10 = (TH1F*) back_1000_1500->Get("signal_ 1.3");
  TH1F* HBB_15 = (TH1F*) back_1500_2000->Get("signal_ 1.3");
  TH1F* HBB_20 = (TH1F*) back_2000_inf->Get("signal_ 1.3");
//  TH1F* HBB_tt = (TH1F*) ttbar->Get("signal_ 1.3");


  TH1F* heve_s = (TH1F*) signal->Get("heve");
  TH1F* heve_5 = (TH1F*) back_500_700->Get("heve");
  TH1F* heve_7 = (TH1F*) back_700_1000->Get("heve");
  TH1F* heve_10 = (TH1F*) back_1000_1500->Get("heve");
  TH1F* heve_15 = (TH1F*) back_1500_2000->Get("heve");
  TH1F* heve_20 = (TH1F*) back_2000_inf->Get("heve");

// merge back_HBB
  TH1F* bHBB = new TH1F("bHBB","",100,0,1);
  TH1F* PVE = new TH1F("Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* M_90 = new TH1F("M90_Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* M_100 = new TH1F("M100_Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* M_125 = new TH1F("M125_Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* M_150 = new TH1F("M150_Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* M_200 = new TH1F("M200_Punzi_vs_delta_eta","",24,0,2.4);
  TH1F* D_13 = new TH1F("Punzi_vs_YMASS","",16,0,16);
  TH1F* DDE = (TH1F*) signal->Get("DDE");

  TH1F* punzi_Net = new TH1F("PN","",100,0,1);
  
  TH1F* effs_b = new TH1F("S vs B","",100,0,1);
  	std::vector<double> xs;
 		xs.push_back(29980);
 		xs.push_back(6334);
 		xs.push_back(1088);
 		xs.push_back(99.11);
 		xs.push_back(20.23);

  	std::vector<double> number;
  		number.push_back(heve_5->GetBinContent(1));
  		number.push_back(heve_7->GetBinContent(1));
  		number.push_back(heve_10->GetBinContent(1));
  		number.push_back(heve_15->GetBinContent(1));
  		number.push_back(heve_20->GetBinContent(1));


  	std::vector<double> weight;
  	double W_B = 0;
  	for (int i = 0 ; i < 5 ; i++){
  		W_B = xs[i]*35900/number[i];
  		weight.push_back(W_B);
  		}
  	for(int i = 0 ; i < 5 ; i++){
  		cout << "weight[i]" << weight[i] << endl;
  		}
  
  	bHBB->Add(HBB_5,weight[0]);
  	bHBB->Add(HBB_7,weight[1]);
  	bHBB->Add(HBB_10,weight[2]);
  	bHBB->Add(HBB_15,weight[3]);
  	bHBB->Add(HBB_20,weight[4]);

//======================================================================
//punzi
//======================================================================
  	double nSigEventTotal = heve_s->GetBinContent(2);
  	double nSigEvent = HBB->Integral();
  	double nBgEvent = bHBB->Integral();
  	int nBin = HBB->GetSize();
  	cout << nBin << endl;
  	vector <vector<float>> effs(2);//2-D vector
  	float Teffs[100], Teffb[100];
  	vector <vector<double>> effs_t(2);//with total
  	vector <vector<float>> effb(2);//bg
  	vector <vector<double>> effb_t(2);//with total
  	vector <vector<double>> punziList(2);
  	double event[2][2] = {0};

  	cout << nSigEvent << "," << nBin << " , " << nBgEvent << endl; 
  	for(int i = 0 ; i < nBin ; i++){ 
  	  	event[0][0] += HBB->GetBinContent(i);//from 0 to end for sig
  	  	event[1][1] += HBB->GetBinContent(nBin-i-2);//from end to 0 for sig
  	  	event[0][1] += bHBB->GetBinContent(i);
  	  	event[1][0] += bHBB->GetBinContent(nBin-i-2);
  	  	effs[0].push_back(event[0][0]/nSigEvent);//(sig/other)
  	  	effs[1].push_back(event[1][1]/nSigEvent);
  	  	effs_t[0].push_back(event[0][0]/nSigEventTotal);//(sig/total)
  	  	effs_t[1].push_back(event[1][1]/nSigEventTotal);
  	  	effb[0].push_back(event[0][1]/nBgEvent);
  	  	effb[1].push_back(event[1][0]/nBgEvent);
		cout << ">> event " << event[1][0] << ": nBgEvent " << nBgEvent << ": EFFB " << effb[1][i] << " : " << i <<  endl;
  	  	punziList[0].push_back(punzi(event[0][0]/nSigEvent,event[0][1]));
  	  	punziList[1].push_back(punzi(event[1][1]/nSigEvent,event[1][0]));
  	  	if(effb[1][i] > 1){
  	  		effb[1][i] = 1;
  	  	}
  	  }
  	  	reverse(effs[1].begin(),effs[1].end());
  		reverse(effb[1].begin(),effb[1].end());
  		reverse(punziList[1].begin(),punziList[1].end());
  	cout << "finished1" << endl;
  	
  	for(int i = 1 ; i < nBin-1 ; i++){
		
  		punzi_Net->SetBinContent(i,punziList[1][i]);
  		Teffs[i] = effs[1][i];
  		Teffb[i] = 1 - effb[1][i];
  		cout << Teffb[i] << " || " << Teffs[i] << " || " << punziList[1][i] << " : " << i << endl; 
  	}
	for(int  i = 1 ; i < 25 ; i++){
		PVE->SetBinContent(i,pun[i-1]);
		M_90->SetBinContent(i,M90[i-1]);
		M_100->SetBinContent(i,M100[i-1]);
		M_125->SetBinContent(i,M125[i-1]);
		M_150->SetBinContent(i,M150[i-1]);
		M_200->SetBinContent(i,M200[i-1]);
	}
	for(int i = 1 ; i < 17 ; i++){
		D_13->SetBinContent(i,D13[i-1]);
	}
	
	float_t punz = 0 ; 
	Int_t max_bin = 0;
	punz = punzi_Net->GetMaximum();
	max_bin = punzi_Net->GetMaximumBin();
	cout << "punzi_maximum = "  << punz << endl;
	cout << "punzi_max_bin = "  << max_bin-1 << endl;
	auto c2 = new TCanvas();
  	TGraph *sigEff_vs_bkgEff = new TGraph(101, Teffs, Teffb);

  	sigEff_vs_bkgEff->SetMaximum(1.2);
  	sigEff_vs_bkgEff->Draw("ap");
  	sigEff_vs_bkgEff->SetTitle("");
  	sigEff_vs_bkgEff->GetXaxis()->SetTitle("Sig Efficiency");
  	sigEff_vs_bkgEff->GetXaxis()->SetLimits(0 , 1.0);
  	sigEff_vs_bkgEff->GetXaxis()->SetTickSize(0.03);
  	sigEff_vs_bkgEff->GetXaxis()->SetTitleOffset(1.2); 
  	sigEff_vs_bkgEff->GetXaxis()->SetLabelOffset(0.015); 
  	sigEff_vs_bkgEff->GetYaxis()->SetTitle("Bkg rejection Efficiency");
  	sigEff_vs_bkgEff->GetYaxis()->SetTitleOffset(1.3);
  	//sigEff_vs_bkgEff->GetYaxis()->SetNdivisions(505);
 	sigEff_vs_bkgEff->GetYaxis()->SetTickSize(0.03);
  	sigEff_vs_bkgEff->GetYaxis()->SetLabelOffset(0.005);
  	sigEff_vs_bkgEff->GetHistogram()->SetMaximum(1.2); // along
  	sigEff_vs_bkgEff->GetHistogram()->SetMinimum(0.);

  	sigEff_vs_bkgEff->SetLineColor(kRed);
  	sigEff_vs_bkgEff->SetLineWidth(3);
  	sigEff_vs_bkgEff->SetMarkerColor(kBlue);
  	sigEff_vs_bkgEff->SetMarkerSize(1);
  	sigEff_vs_bkgEff->SetMarkerStyle(20);
  	sigEff_vs_bkgEff->Draw("p");
	c2->Print("sgrejection_particle_1.3.png");
	
    auto c1 = new TCanvas();
		
  		punzi_Net->GetXaxis()->SetTitle("particle net");
  		punzi_Net->GetYaxis()->SetTitle("Punzi_significance");
  		punzi_Net->GetYaxis()->SetTitleOffset(1.4);
		punzi_Net->SetLineColor(kRed);
  		punzi_Net->SetLineWidth(2);
  		punzi_Net->Draw("hist");
		c1->Print("punzi_particle_m150_1.3.png");
 	
	 auto c3 = new TCanvas();
	 	HBB->SetLineColor(kRed);
	 	HBB->Draw("hist");
		HBB->GetXaxis()->SetTitle("particle net");
		HBB->SetTitle("QCD:sig_vs_bg");
		bHBB->Draw("hist same");
		c3->Print("bg+sig_m150_1.3.png");
	auto c4 = new TCanvas();
		D_13->GetYaxis()->SetTitle("Punzi_significace");
		D_13->GetXaxis()->SetTitle("YMASS_range");
		D_13->GetYaxis()->SetTitleOffset(1.4);
		D_13->Draw("hist");
		c4->Print("punz_vs_ymass.png");
	auto c5 = new TCanvas();
		DDE->GetXaxis()->SetTitle("delta_eta_value");
		DDE->Draw("hist");
		c5->Print("delta_eta_value.png");
/*
	auto c4 = new TCanvas();
		PVE->GetYaxis()->SetTitle("Punzi_significace");
  		PVE->GetXaxis()->SetTitle("delta_eta_value");
		PVE->Draw("hist");
		c4->Print("punz_vs_delta.png");
	auto c5 = new TCanvas();
		M_90->GetYaxis()->SetTitle("Punzi_significace");
  		M_90->GetXaxis()->SetTitle("delta_eta_value");
		M_90->Draw("hist");
		c5->Print("punz_vs_delta_90.png");
	auto c6 = new TCanvas();
		M_100->GetYaxis()->SetTitle("Punzi_significace");
  		M_100->GetXaxis()->SetTitle("delta_eta_value");
		M_100->Draw("hist");
		c6->Print("punz_vs_delta_100.png");
	auto c7 = new TCanvas();
		M_125->GetYaxis()->SetTitle("Punzi_significace");
  		M_125->GetXaxis()->SetTitle("delta_eta_value");
		M_125->Draw("hist");
		c7->Print("punz_vs_delta_125.png");
*/
  	
/*
    TCanvas* c1 = new TCanvas("c1","",700,1000);  
    
    THStack *hs = new THStack("hs","histograms");
  	bHBB->SetFillColor(kBlue);
  	hs->Add(bHBB);
  	HBB_tt->SetFillColor(kGreen);
  	hs->Add(HBB_tt);

    hs->Draw();
*/

//  	auto c3 = new TCanvas();
//  		sigEff_vs_bkgEff->Draw("l");

}

