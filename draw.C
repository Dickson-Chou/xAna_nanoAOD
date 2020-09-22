// draw the graph after the selection
// example code to run 2016 NSSM MC X->Y+H samples
// .L xAna_nano_nssm.C++
// xAna_nano_nssm("test.root") or xAna_nano_nssm("input.txt")
// example root file is at /afs/cern.ch/work/s/syu/public/forTiKai/nssm_nano.root

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TLorentzVector.h>
#include "setNCUStyle.C"

using namespace std;

void draw(std::string filename="*.root", std::string outputFileName="result.root"){
           
  std::vector<std::string> inputFiles;

  // check first if this is a root file
  if(filename.find(".root")!=string::npos)
    {
      cout << "This is a single input root file" << endl;
      inputFiles.push_back(filename);
    }
  else // assume this is a text file
    {
      cout << "This is a text file " << endl;
      ifstream fin;
      fin.open(filename.data());
      string temp;
      fin >> temp;
      while(!fin.eof())
	{
	  inputFiles.push_back(temp);
	  fin >> temp;
	}
      cout << "There are " << inputFiles.size() << " files" << endl;
      for(unsigned int ifile=0; ifile < inputFiles.size(); ifile++)
	cout << "Input file " << ifile << " is " << inputFiles[ifile] << endl;
    }
  cout << "Output file name is " << outputFileName << endl;

  setNCUStyle(true);
  //get TTree from file ...  
  TreeReader data(inputFiles,"variable");
  TTree* thisTree = data.GetTree();

 
//  Float_t*  FatJet_ParticleNetMD_ProbXbb = data.GetPtrFloat("FatJet_ParticleNetMD_ProbXbb");
//  Float_t*  I_FatJet_ParticleNetMD_ProbXbb = data.GetPtrFloat("I_FatJet_ParticleNetMD_ProbXbb");


  //get TH1F data from file ...
  TFile* integral = new TFile("histo.root","READ");

  TH1F* Ydistri = (TH1F*) integral->Get("Ydistri");
  TH1F* Ydistri_P1 = (TH1F*) integral->Get("Ydistri_P1");
  TH1F* Ydistri_P2 = (TH1F*) integral->Get("Ydistri_P2");
  TH1F* Ydistri_P3 = (TH1F*) integral->Get("Ydistri_P3");
  TH1F* Ydistri_P4 = (TH1F*) integral->Get("Ydistri_P4");
  TH1F* Ydistri_P5 = (TH1F*) integral->Get("Ydistri_P5");
  TH1F* Ydistri_P6 = (TH1F*) integral->Get("Ydistri_P6");
  TH1F* Ydistri_P7 = (TH1F*) integral->Get("Ydistri_P7");
  TH1F* Ydistri_P8 = (TH1F*) integral->Get("Ydistri_P8");

  TH1F* Unequal_Ydistri_initial = (TH1F*) integral->Get("Unequal_Ydistri_initial");
  TH1F* Unequal_Ydistri_P1 = (TH1F*) integral->Get("Unequal_Ydistri_P1");
  TH1F* Unequal_Ydistri_P2 = (TH1F*) integral->Get("Unequal_Ydistri_P2");
  TH1F* Unequal_Ydistri_P3 = (TH1F*) integral->Get("Unequal_Ydistri_P3");
  TH1F* Unequal_Ydistri_P4 = (TH1F*) integral->Get("Unequal_Ydistri_P4");
  TH1F* Unequal_Ydistri_P5 = (TH1F*) integral->Get("Unequal_Ydistri_P5");
  TH1F* Unequal_Ydistri_P6 = (TH1F*) integral->Get("Unequal_Ydistri_P6");
  TH1F* Unequal_Ydistri_P7 = (TH1F*) integral->Get("Unequal_Ydistri_P7");
  TH1F* Unequal_Ydistri_P8 = (TH1F*) integral->Get("Unequal_Ydistri_P8");

  TH1F* HBB = (TH1F*) integral->Get("HBB");
  TH1F* hbb = (TH1F*) integral->Get("hbb");
  TH1F* HDeep = (TH1F*) integral->Get("HDeep");
  TH1F* Hdeep = (TH1F*) integral->Get("Hdeep");
 // TH1F* HBB = new TH1F("HBB"," ",100,0,1);


//  TCanvas *c1 = new TCanvas("c1","c1",3);
  
  TH1F* Total = (TH1F*) integral->Get("heve");
  TH1F* eff_final = new TH1F("Eff_final","",16,0,16);
  TH1F* eff1 = new TH1F("EFF1","",16,0,16);
  TH1F* eff2 = new TH1F("EFF2","",16,0,16);
  TH1F* eff3 = new TH1F("EFF3","",16,0,16);
  TH1F* eff4 = new TH1F("EFF4","",16,0,16);
  TH1F* eff5 = new TH1F("EFF5","",16,0,16);
  TH1F* eff6 = new TH1F("EFF6","",16,0,16);
  TH1F* eff7 = new TH1F("EFF7","",16,0,16);
  TH1F* eff8 = new TH1F("EFF8","",16,0,16);

  TH1F* eff_final_UN = new TH1F("Eff_final_UN","",16,0,16);
  TH1F* eff1_UN = new TH1F("EFF1_UN","",16,0,16);
  TH1F* eff2_UN = new TH1F("EFF2_UN","",16,0,16);
  TH1F* eff3_UN = new TH1F("EFF3_UN","",16,0,16);
  TH1F* eff4_UN = new TH1F("EFF4_UN","",16,0,16);
  TH1F* eff5_UN = new TH1F("EFF5_UN","",16,0,16);
  TH1F* eff6_UN = new TH1F("EFF6_UN","",16,0,16);
  TH1F* eff7_UN = new TH1F("EFF7_UN","",16,0,16);
  TH1F* eff8_UN = new TH1F("EFF8_UN","",16,0,16);


  
  		//definition
  Float_t YM[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};//vector of ymass points
  double ym[17]={70,95,120,145,195,245,295,395,495,595,695,795,895,995,1195,1395,1595};
  std::vector<float> YDistri_initial;
  std::vector<float> YDistri;
  std::vector<float> total;
  std::vector<float> Ryy(16,0);
  std::vector<float> Ryt(16,0);
  

  

  for(int i = 1 ; i < 17 ; i++ ){
//  	YDistri_initial.push_back(Ydistri_initial->GetBinContent(i));
  	YDistri.push_back(Ydistri->GetBinContent(i));
  	total.push_back(Total->GetBinContent(3));
  }
 
  //th1f distri draw
 
/*
  for(int i = 0 ; i < 16 ; i++){

  	Ryy[i] = YDistri[i]/YDistri_initial[i];
  	Ryt[i] = YDistri[i]/total[i];
  	effyy->Fill(i,Ryy[i]);
  	effyy->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	effyt->Fill(i,Ryt[i]);
  	effyt->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  }
 */ 
//  cout << "initial = " << YDistri_initial[2] << endl;
//  cout << "Ydistri = " << YDistri[2] << endl;
//  cout << "npass[1] = " << total[2] << endl;

//divide::
  auto c1 = new TCanvas();
/*
  eff1_UN->Divide(Unequal_Ydistri_P1,Unequal_Ydistri_initial,1,1,"B");
  eff2_UN->Divide(Unequal_Ydistri_P2,Unequal_Ydistri_initial,1,1,"B");
  eff3_UN->Divide(Unequal_Ydistri_P3,Unequal_Ydistri_initial,1,1,"B");
  eff4_UN->Divide(Unequal_Ydistri_P4,Unequal_Ydistri_initial,1,1,"B");
  eff5_UN->Divide(Unequal_Ydistri_P5,Unequal_Ydistri_initial,1,1,"B");
  eff6_UN->Divide(Unequal_Ydistri_P6,Unequal_Ydistri_initial,1,1,"B");
  eff7_UN->Divide(Unequal_Ydistri_P7,Unequal_Ydistri_initial,1,1,"B");
  eff8_UN->Divide(Unequal_Ydistri_P8,Unequal_Ydistri_initial,1,1,"B");

 
  eff1->Divide(Ydistri_P1,Ydistri_initial,1,1,"B");
  eff2->Divide(Ydistri_P2,Ydistri_initial,1,1,"B");
  eff3->Divide(Ydistri_P3,Ydistri_initial,1,1,"B");
  eff4->Divide(Ydistri_P4,Ydistri_initial,1,1,"B");
  eff5->Divide(Ydistri_P5,Ydistri_initial,1,1,"B");
  eff6->Divide(Ydistri_P6,Ydistri_initial,1,1,"B");
  eff7->Divide(Ydistri_P7,Ydistri_initial,1,1,"B");
  eff8->Divide(Ydistri_P8,Ydistri_initial,1,1,"B");



  for(int i = 0 ; i < 16 ; i++){
  	
  	eff1->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff2->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
 	eff3->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
 	eff4->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i])); 
  	eff5->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff6->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff7->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff8->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	
  	eff3_UN->Fill(YM[i]);
  	eff3_UN->SetBins(16,ym);
  }

  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat("");//cancel the label
  	
 // eff->GetMaximum()
 // subleading cut eff
  eff1->SetMinimum(0);
  eff1->SetLineColor(kRed);
  eff1->SetMarkerColor(kRed);
  eff1->SetMarkerStyle(kFullCircle);
  eff2->SetLineColor(kBlue+2);
  eff2->SetMarkerColor(kBlue+2);
  eff2->SetMarkerStyle(kFullSquare);
  eff3->SetLineColor(kGreen+3);
  eff3->SetMarkerColor(kGreen+3);
  eff3->SetMarkerStyle(kFullTriangleUp);
  eff4->SetLineColor(kRed+3);
  eff4->SetMarkerColor(kRed+3);
  eff4->SetMarkerStyle(kFullTriangleDown);

  eff1->Draw();
  eff2->Draw("SAME");
  eff3->Draw("SAME");
  eff4->Draw("SAME");

  auto legend = new TLegend(0.2,0.2,0.4,0.4);//(x1,y1,x2,y2)
  legend->SetHeader("cut2 to cut 5","C"); // option "C" allows to center the header
  legend->SetNColumns(2);
  legend->AddEntry(eff1,"cut2","lep");
  legend->AddEntry(eff2,"cut3","lep");
  legend->AddEntry(eff3,"cut4","lep");
  legend->AddEntry(eff4,"cut5","lep");
  legend->Draw();

  c1->Modified();
  c1->Print("Divide_eff_sub.pdf");

  c1->Update();
  eff5->SetMinimum(0);
  eff5->SetLineColor(kRed);
  eff5->SetMarkerColor(kRed);
  eff5->SetMarkerStyle(kFullCircle);
  eff6->SetLineColor(kBlue+2);
  eff6->SetMarkerColor(kBlue+2);
  eff6->SetMarkerStyle(kFullSquare);
  eff7->SetLineColor(kGreen+3);
  eff7->SetMarkerColor(kGreen+3);
  eff7->SetMarkerStyle(kFullTriangleUp);
  eff8->SetLineColor(kRed+3);
  eff8->SetMarkerColor(kRed+3);
  eff8->SetMarkerStyle(kFullTriangleDown);


//  eff5->Draw();
//  eff6->Draw("SAME");
  eff7->Draw();
  eff8->Draw("SAME");

  auto legend_1 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
  legend_1->SetHeader("cut6 to cut 10","C"); // option "C" allows to center the header
  legend_1->SetNColumns(2);
  legend_1->AddEntry(eff5,"cut6","lep");
  legend_1->AddEntry(eff6,"cut7","lep");
  legend_1->AddEntry(eff7,"cut8","lep");
  legend_1->AddEntry(eff8,"cut9","lep");

  legend_1->Draw();

  c1->Modified();
  c1->Print("Divide_eff_lead.pdf");

  

  auto c2 = new TCanvas();

  
  c2->Divide(2,1);
//  auto h1 = (TH1F*)Ydistri_initial->Clone();
//  h1->Reset();
  c2->cd(1);
  Ydistri_P3->SetLineColor(kBlue+3);
  Ydistri_P4->SetLineColor(kOrange+9);

  
  Ydistri_P3->SetFillColor(kBlue+3);
  Ydistri_P4->SetFillColor(kOrange+9);
//  Ydistri_P3->SetTitle("Y Distribution of N-1");
  Ydistri_P3->GetXaxis()->SetTitle("Y Mass points");
  Ydistri_P3->GetYaxis()->SetTitle("Number of Y");
  Ydistri_P3->Draw("hist");
  Ydistri_P4->Draw("hist SAME");

  auto legend_2 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_2->SetHeader("Distri of ","C"); // option "C" allows to center the header
//  legend_2->SetNColumns(2);
  legend_2->AddEntry(Ydistri_P3," no particle net","f");
  legend_2->AddEntry(Ydistri_P4,"particle net","f");
  legend_2->Draw();

  c2->cd(2);
  eff3->SetLineColor(kBlue+3);
  eff3->SetMarkerColor(kBlue+3);
  eff3->SetMarkerStyle(kFullCircle);
  eff4->SetLineColor(kOrange+9);
  eff4->SetMarkerColor(kOrange+9);
  eff4->SetMarkerStyle(kFullCircle);
//  eff3->GetXaxis()->SetTitle("Y Mass points");
//  eff3->GetYaxis()->SetTitle("efficiency");
//  eff3->GetXaxis()->SetTitleSize(1);
//  eff3->GetYaxis()->SetTitleSize(1);

  eff3->Draw();
  eff4->Draw("SAME");

  auto legend_3 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_3->SetHeader("cut6 to cut 10","C"); // option "C" allows to center the header
//  legend_1->SetNColumns(2);
  legend_3->AddEntry(eff3,"no particle net","lep");
  legend_3->AddEntry(eff4,"particle net","lep");
  legend_3->Draw();

  c2->Modified();
  c2->Print("N-1_Net.pdf");

  auto c3 = new TCanvas();
  Unequal_Ydistri_P3->SetLineColor(kBlue+3);
  Unequal_Ydistri_P4->SetLineColor(kOrange+9);

  Unequal_Ydistri_P3->SetFillColor(kBlue+3);
  Unequal_Ydistri_P4->SetFillColor(kOrange+9);

  Unequal_Ydistri_P3->Draw("hist");
  Unequal_Ydistri_P4->Draw("hist SAME");

  auto legend_4 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_2->SetHeader("Distri of ","C"); // option "C" allows to center the header
//  legend_2->SetNColumns(2);
  legend_4->AddEntry(Unequal_Ydistri_P3," no particle net","f");
  legend_4->AddEntry(Unequal_Ydistri_P4,"particle net","f");
  legend_4->Draw();

  c3->Modified();
  c3->Print("N-1_Net_un.pdf");
  
*/
  HBB->SetLineColor(kBlue+3);
  HBB->SetFillColor(kBlue+3);
  hbb->SetLineColor(kOrange);
  hbb->SetFillColor(kOrange);
  HBB->Draw();
  hbb->Draw("SAME");

  auto l = new TLegend(0.15,0.7,0.35,0.9);
  l->AddEntry(HBB,"No particle net","f");
  l->AddEntry(hbb,"paritcle net","f");
  l->Draw();
  c1->Draw();

  auto c2 = new TCanvas();

  HDeep->SetLineColor(kGreen+3);
  HDeep->SetFillColor(kGreen+3);
  Hdeep->SetLineColor(kBlue-2);
  Hdeep->SetFillColor(kBlue-2);
  HDeep->Draw();
  Hdeep->Draw("SAME");

  auto l2 = new TLegend(0.15,0.7,0.35,0.9);
  l2->AddEntry(HDeep,"No DeepAKB","f");
  l2->AddEntry(Hdeep,"DeepAKB","f");
  l2->Draw();
  c2->Draw();
//TGraphAsymmErrors::

//  c1->Update();

//  auto gr = new TGraphAsymmErrors(Ydistri , Ydistri_initial , "cl = 0.683 b(1,1) mode");

//  Ydistri_P1->Draw("PLC PMC");
//  Ydistri_P2->Draw("SAME PLC PMC");
//  c1->Modified();
//  c1->Print("TGraphAsymmErrors_eff.pdf");
  
  //store 
  TFile* outFile = new TFile(outputFileName.data(),"recreate");
//  Ydistri_initial->Write();
//  eff1->Write();
//  eff2->Write();
//  eff3->Write();
//  eff3_UN->Draw();
//  gr->Write();
  outFile->Close();
}






