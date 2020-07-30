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

using namespace std;
void xAna_nano_nssm(std::string filename="devdatta_nanoAOD_nssm.root", std::string outputFileName="histo.root"){
           
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
  
  //get TTree from file ...  
  TreeReader data(inputFiles,"Events");
  TTree* thisTree = data.GetTree();

  // check if this tree has the required branches for trigger
  std::vector<std::string> trigNames;
  trigNames.push_back("HLT_PFHT800");
  trigNames.push_back("HLT_PFHT900");
  trigNames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5");
  trigNames.push_back("HLT_AK8PFJet360_TrimMass30");
  trigNames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
  trigNames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50");
  trigNames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");

  std::vector<std::string> tempNames = trigNames;

  // first check if the tirgger-branch exists or not
  for(unsigned int itrig=0; itrig < tempNames.size(); itrig++)
    {
      TBranch* thisBranch = thisTree->FindBranch(tempNames[itrig].data());
      if(thisBranch==NULL){
	cerr << "Branch: " << tempNames[itrig] << " is not present in the tree! " << endl;
	trigNames.erase(trigNames.begin()+itrig);
      }
    }

  unsigned int nTrigs= trigNames.size();
  cout << "Available number of trigger paths = " << nTrigs << endl;
  
  // second check if the Ymass-brance exists or not
  std::vector<std::string> YmassNames;
  YmassNames.push_back("GenModel_YMass_90");
  YmassNames.push_back("GenModel_YMass_100");
  YmassNames.push_back("GenModel_YMass_125");
  YmassNames.push_back("GenModel_YMass_150");

  std::vector<std::string> TYmassNames = YmassNames;

  for(unsigned int iymass=0; iymass < TYmassNames.size(); iymass++)
  {
    TBranch* massBranch = thisTree->FindBranch(TYmassNames[iymass].data());
    if(massBranch==NULL){
  cerr << "Branch: " << TYmassNames[iymass] << " is not present in the tree! " << endl;
  YmassNames.erase(YmassNames.begin()+iymass);
    }
  }

  unsigned int nYmass = YmassNames.size();
  cout << "Available number of Ymass value = " << nYmass << endl;



  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  Int_t Ymassp[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};
  double_t YM[2]={90,100};
  Float_t eff = 0;
  double npassym_90 = 0 ;
  double npassym_100 = 0;
  int npassym_125 = 0;
  int npassym_150 = 0;
  int npassym_200 = 0;
  int npassym_250 = 0;
  int npassym_300 = 0;
  int npassym_400 = 0;
  int npassym_500 = 0;
  int npassym_600 = 0;
  int npassym_700 = 0;
  int npassym_800 = 0;
  int npassym_900 = 0;
  int npassym_1000 = 0;
  int npassym_1200 = 0;
  int npassym_1400 = 0;
  
  const unsigned int nLabels=21;
  TH1F* heve=new TH1F("heve","",nLabels,-0.5,20.5);
  TH1F* FatJetPt = new TH1F("FatPt","",400,0,400 );
  TH1F* YMass = new TH1F("YMass","",2,0,2);
  


  heve->SetYTitle("Number of Events");
  heve->LabelsOption("v");
  YMass->SetYTitle("Number of YMass");
  YMass->LabelsOption("v");
  //store the variable into the tree called data
  TTree* variable = new TTree("variable", "variable");
  
  Float_t FatJet_Pt[2],FatJet_Eta[2],FatJet_Tau1[2],FatJet_Tau2[2],FatJet_Tau21[2],FatJet_Mass[2],FatJet_Msoftdrop[2],FatJet_BtagHbb[2],FatJet_ParticleNetMD_ProbXbb[2],FatJet_DeepTagMD_ZHbbvsQCD[2];
  Bool_t GenModel_YMAss_150,GenModel_YMAss_90,GenModel_YMAss_100,GenModel_YMAss_125,GenModel_YMAss_200,GenModel_YMAss_250,GenModel_YMAss_300,GenModel_YMAss_400;
  Bool_t GenModel_YMAss_500,GenModel_YMAss_600,GenModel_YMAss_700,GenModel_YMAss_800,GenModel_YMAss_900,GenModel_YMAss_1000,GenModel_YMAss_1200,GenModel_YMAss_1400;
  
  variable->Branch("FatJet_Pt",&FatJet_Pt,"FatJet_Pt[2]/F");
  variable->Branch("FatJet_Eta",&FatJet_Eta,"FatJet_Eta[2]/F");
  variable->Branch("FatJet_Tau1",&FatJet_Tau1,"FatJet_Tau1[2]/F");
  variable->Branch("FatJet_Tau2",&FatJet_Tau2,"FatJet_Tau2[2]/F");
  variable->Branch("FatJet_Tau21",&FatJet_Tau21,"FatJet_Tau21[2]/F");
  variable->Branch("FatJet_Mass",&FatJet_Mass,"FatJet_Mass[2]/F");
  variable->Branch("FatJet_Msoftdrop",&FatJet_Msoftdrop,"FatJet_Msoftdrop[2]/F");
  variable->Branch("FatJet_BtagHbb",&FatJet_BtagHbb,"FatJet_BtagHbb[2]/F");
  variable->Branch("FatJet_ParticleNetMD_ProbXbb",&FatJet_ParticleNetMD_ProbXbb,"FatJet_ParticleNetMD_ProbXbb[2]/F");
  variable->Branch("FatJet_DeepTagMD_ZHbbvsQCD",&FatJet_DeepTagMD_ZHbbvsQCD,"FatJet_DeepTagMD_ZHbbvsQCD[2]/F");
  variable->Branch("GenModel_YMAss_90",&GenModel_YMAss_90,"GenModel_YMAss_90/O");
  variable->Branch("GenModel_YMAss_100",&GenModel_YMAss_100,"GenModel_YMAss_100/O");
  variable->Branch("GenModel_YMAss_125",&GenModel_YMAss_125,"GenModel_YMAss_125/O");
  variable->Branch("GenModel_YMAss_150",&GenModel_YMAss_150,"GenModel_YMAss_150/O");
  variable->Branch("GenModel_YMAss_200",&GenModel_YMAss_200,"GenModel_YMAss_200/O");
  variable->Branch("GenModel_YMAss_250",&GenModel_YMAss_250,"GenModel_YMAss_250/O");
  variable->Branch("GenModel_YMAss_300",&GenModel_YMAss_300,"GenModel_YMAss_300/O"); 
  variable->Branch("GenModel_YMAss_400",&GenModel_YMAss_400,"GenModel_YMAss_400/O");
  variable->Branch("GenModel_YMAss_500",&GenModel_YMAss_500,"GenModel_YMAss_500/O");
  variable->Branch("GenModel_YMAss_600",&GenModel_YMAss_600,"GenModel_YMAss_600/O");
  variable->Branch("GenModel_YMAss_700",&GenModel_YMAss_700,"GenModel_YMAss_700/O");
  variable->Branch("GenModel_YMAss_800",&GenModel_YMAss_800,"GenModel_YMAss_800/O");
  variable->Branch("GenModel_YMAss_900",&GenModel_YMAss_900,"GenModel_YMAss_900/O");
  variable->Branch("GenModel_YMAss_1000",&GenModel_YMAss_1000,"GenModel_YMAss_1000/O");
  variable->Branch("GenModel_YMAss_1200",&GenModel_YMAss_1200,"GenModel_YMAss_1200/O");
  variable->Branch("GenModel_YMAss_1400",&GenModel_YMAss_1400,"GenModel_YMAss_1400/O");

  const char *label[nLabels];
  label[0]="Total";
  for(unsigned int i=1; i< nLabels; i++){
    label[i] = Form("Cut %d",i);
  }

  // start looping over events
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;
    heve->Fill(label[0],1.);
    //1. trigger 

    Bool_t passTrigger=false;
      
    for(unsigned int itrig=0; itrig< nTrigs; itrig++)
      {
	if(data.GetBool(trigNames[itrig].data())==true)
	  {
	    passTrigger=true;
	    break;
	  }
      } // end of loop over required trigger paths

    if(!passTrigger)continue;
    nPass[0]++;
    heve->Fill(label[1],1.);

    // loop over fatjets
    Float_t*  FatJet_pt = data.GetPtrFloat("FatJet_pt");
    Float_t*  FatJet_eta = data.GetPtrFloat("FatJet_eta");
    Float_t*  FatJet_tau1 = data.GetPtrFloat("FatJet_tau1");
    Float_t*  FatJet_tau2 = data.GetPtrFloat("FatJet_tau2");
    Float_t*  FatJet_mass = data.GetPtrFloat("FatJet_mass");
    Float_t*  FatJet_msoftdrop = data.GetPtrFloat("FatJet_msoftdrop");
    Float_t*  FatJet_btagHbb = data.GetPtrFloat("FatJet_btagHbb");
    Float_t*  FatJet_deepTagMD_ZHbbvsQCD = data.GetPtrFloat("FatJet_deepTagMD_ZHbbvsQCD");
    Float_t*  FatJet_ParticleNetMD_probXbb = data.GetPtrFloat("FatJet_ParticleNetMD_probXbb");
/*    Bool_t   GenModel_YMass_150 = data.GetBool("GenModel_YMass_150");
    Bool_t   GenModel_YMass_90 = data.GetBool("GenModel_YMass_90");
    Bool_t   GenModel_YMass_100 = data.GetBool("GenModel_YMass_100");
    Bool_t   GenModel_YMass_125 = data.GetBool("GenModel_YMass_125");
    Bool_t   GenModel_YMass_200 = data.GetBool("GenModel_YMass_200");
    Bool_t   GenModel_YMass_250 = data.GetBool("GenModel_YMass_250");
    Bool_t   GenModel_YMass_300 = data.GetBool("GenModel_YMass_300");
    Bool_t   GenModel_YMass_400 = data.GetBool("GenModel_YMass_400");
    Bool_t   GenModel_YMass_500 = data.GetBool("GenModel_YMass_500");
    Bool_t   GenModel_YMass_600 = data.GetBool("GenModel_YMass_600");
    Bool_t   GenModel_YMass_700 = data.GetBool("GenModel_YMass_700");
    Bool_t   GenModel_YMass_800 = data.GetBool("GenModel_YMass_800");
    Bool_t   GenModel_YMass_900 = data.GetBool("GenModel_YMass_900");
    Bool_t   GenModel_YMass_1000 = data.GetBool("GenModel_YMass_1000");
    Bool_t   GenModel_YMass_1200 = data.GetBool("GenModel_YMass_1200");
    Bool_t   GenModel_YMass_1400 = data.GetBool("GenModel_YMass_1400"); */



    UInt_t nFatJet = data.GetInt("nFatJet");
    UInt_t nGoodPair = 0;
    Float_t Max_pt = 0;
    Float_t sub_Max_pt = 0;
    Int_t leadingID = 0;
    Int_t subleadingID = 0;
    Float_t FatJet_tau21_i = 0;
    Float_t FatJet_tau21_j = 0;

    Bool_t masspass = false;
    

    for(UInt_t ij=0; ij < nFatJet; ij++){

      if(FatJet_pt[ij] < 300)continue;
      if(fabs(FatJet_eta[ij]) > 2.4)continue;
      Float_t tau21_i = FatJet_tau2[ij]/FatJet_tau1[ij];
      if(tau21_i > 0.55)continue;
      if(!(FatJet_deepTagMD_ZHbbvsQCD[ij] > 0.8 && FatJet_ParticleNetMD_probXbb[ij] > 0.85))continue;



      for(UInt_t jj=0; jj < ij; jj++){

	     if(FatJet_pt[jj] < 300)continue;
	     if(fabs(FatJet_eta[jj]) > 2.4)continue;
	     Float_t tau21_j = FatJet_tau2[jj]/FatJet_tau1[jj];
	     if(tau21_j > 0.55)continue;
       if(!(FatJet_deepTagMD_ZHbbvsQCD[jj] > 0.8 && FatJet_ParticleNetMD_probXbb[jj] > 0.85))continue;


	// eta difference
	     if(fabs(FatJet_eta[ij]-FatJet_eta[jj]) > 1.3)continue;
       
       if(FatJet_pt[jj] > Max_pt){
        Max_pt = FatJet_pt[jj];leadingID = jj;
      }
       if(jj == leadingID){
        FatJet_tau21_j = tau21_j;
      }



	     FatJetPt->Fill(FatJet_pt[ij]);
	
	     nGoodPair++;
	
      } // end of inner jet loop

      if(FatJet_pt[ij] > sub_Max_pt && sub_Max_pt < Max_pt){
        sub_Max_pt = FatJet_pt[ij] ;
        subleadingID = ij;
      
      }
      if(ij == subleadingID){
        FatJet_tau21_i = tau21_i;
      }
      
    } // end of outer jet loop

// store the leading & subleading information into the vector

    FatJet_Pt[0] = FatJet_pt[leadingID];
    FatJet_Pt[1] = FatJet_pt[subleadingID];
    FatJet_Eta[0] = FatJet_eta[leadingID];
    FatJet_Eta[1] = FatJet_eta[subleadingID];
    FatJet_Tau1[0] = FatJet_tau1[leadingID];
    FatJet_Tau1[1] = FatJet_tau1[subleadingID];
    FatJet_Tau2[0] = FatJet_tau2[leadingID];
    FatJet_Tau2[1] = FatJet_tau2[subleadingID];
    FatJet_Tau21[0] = FatJet_tau21_j ;
    FatJet_Tau21[1] = FatJet_tau21_i ;
    FatJet_Mass[0] = FatJet_mass[leadingID];
    FatJet_Mass[1] = FatJet_mass[subleadingID];
    FatJet_Msoftdrop[0] = FatJet_msoftdrop[leadingID];
    FatJet_Msoftdrop[1] = FatJet_msoftdrop[subleadingID];
    FatJet_BtagHbb[0] = FatJet_btagHbb[leadingID];
    FatJet_BtagHbb[1] = FatJet_btagHbb[subleadingID];
    FatJet_DeepTagMD_ZHbbvsQCD[0] = FatJet_deepTagMD_ZHbbvsQCD[leadingID];
    FatJet_DeepTagMD_ZHbbvsQCD[1] = FatJet_deepTagMD_ZHbbvsQCD[subleadingID];
    FatJet_ParticleNetMD_ProbXbb[0] = FatJet_ParticleNetMD_probXbb[leadingID];
    FatJet_ParticleNetMD_ProbXbb[1] = FatJet_ParticleNetMD_probXbb[subleadingID];
    
    
//    cout << npassym_200 << " , " << nPass[1] << endl;
    
    if(nGoodPair<1)continue;
    nPass[1]++;
    heve->Fill(label[2],1.);
    // get number of each value of Ymass
       if(GenModel_YMAss_90){
        npassym_90+=1;
        
       }
       if(GenModel_YMAss_100){
        npassym_100+=1;
       }
       if(GenModel_YMAss_125){
        npassym_125++;
       }
       if(GenModel_YMAss_150){
        npassym_150++;
       }
       if(GenModel_YMAss_200){
        npassym_200++;
       }
       if(GenModel_YMAss_250){
        npassym_250++;
       }
       if(GenModel_YMAss_300){
        npassym_300++;
       }
       if(GenModel_YMAss_400){
        npassym_400++;
       }
       if(GenModel_YMAss_500){
        npassym_500++;
       }
       if(GenModel_YMAss_600){
        npassym_600++;
       }
       if(GenModel_YMAss_700){
        npassym_700++;
       }
       if(GenModel_YMAss_800){
        npassym_800++;
       }
       if(GenModel_YMAss_900){
        npassym_900++;
       }
       if(GenModel_YMAss_1000){
        npassym_1000++;
       }
       if(GenModel_YMAss_1200){
        npassym_1200++;
       }
       if(GenModel_YMAss_1400){
        npassym_1400++;
       }

    GenModel_YMAss_90  = GenModel_YMass_90;
    GenModel_YMAss_150 = GenModel_YMass_150;
    GenModel_YMAss_100 = GenModel_YMass_100;
    GenModel_YMAss_125 = GenModel_YMass_125;
    GenModel_YMAss_200 = GenModel_YMass_200;
    GenModel_YMAss_250 = GenModel_YMass_250;
    GenModel_YMAss_300 = GenModel_YMass_300;
    GenModel_YMAss_400 = GenModel_YMass_400;
    GenModel_YMAss_500 = GenModel_YMass_500;
    GenModel_YMAss_600 = GenModel_YMass_600;
    GenModel_YMAss_700 = GenModel_YMass_700;
    GenModel_YMAss_800 = GenModel_YMass_800;
    GenModel_YMAss_900 = GenModel_YMass_900;
    GenModel_YMAss_1000 = GenModel_YMass_1000;
    GenModel_YMAss_1200 = GenModel_YMass_1200;
    GenModel_YMAss_1400 = GenModel_YMass_1400;
    YMass->Fill(GenModel_YMAss_150);
    variable->Fill();


    

    

  } // end of loop over events
  Int_t NY[16]={npassym_125,npassym_150,npassym_200,npassym_300,npassym_400,npassym_500,npassym_600,npassym_700,npassym_800,npassym_900,npassym_1000,npassym_1200,npassym_1400};
  double_t ny[2]={npassym_90,npassym_100};
  cout << "Number of total events = " << nTotal << endl;

  for(int i=0;i<20;i++)
    if(nPass[i]>0)cout << "nPass["<<i<<"]= " << nPass[i] << endl;

  Float_t a = 0 ;
  Float_t b = 0 ; 
  a = nTotal;
  b = nPass[1];

  cout << " Efficiency = " << ((b)/a)*100 << "%" << endl;
  cout << " Npassym_90 = " << npassym_90 << endl;
  cout << " Npassym_100 = " << npassym_100 << endl;
  cout << " Npassym_125 = " << npassym_125 << endl;
  cout << " Npassym_150 = " << npassym_150 << endl;
  cout << " Npassym_200 = " << npassym_200 << endl;
  cout << " Npassym_250 = " << npassym_250 << endl;
  cout << " Npassym_300 = " << npassym_300 << endl;
  cout << " Npassym_400 = " << npassym_400 << endl;
  cout << " Npassym_500 = " << npassym_500 << endl;
  cout << " Npassym_600 = " << npassym_600 << endl;
  cout << " Npassym_700 = " << npassym_700 << endl;
  cout << " Npassym_800 = " << npassym_800 << endl;
  cout << " Npassym_900 = " << npassym_900 << endl;
  cout << " Npassym_1000 = " << npassym_1000 << endl;
  cout << " Npassym_1200 = " << npassym_1200 << endl;
  cout << " Npassym_1400 = " << npassym_1400 << endl;
  cout << " Ny_total = " << npassym_1400+npassym_1200+npassym_1000+npassym_900+npassym_800+npassym_700+npassym_600+npassym_500+npassym_400+npassym_300+npassym_250+npassym_200+npassym_150+npassym_100+npassym_125+npassym_90 << endl;
  TGraph* gr = new TGraph(16,Ymassp,NY);
  gr->SetTitle("Ymass_distribution");
  gr->SetMarkerStyle(20);
  TExec *ex = new TExec("ex","drawtext();");
  gr->GetListOfFunctions()->Add(ex);
  gr->Draw("ACP");
  // writing example output file
  TFile* outFile = new TFile(outputFileName.data(),"recreate");
  gr->Write();
  heve->Write();
  FatJetPt->Write();
  YMass->Write();
  variable->Write();
  outFile->Close();
  
  
} // end of macro

void drawtext()
  {
   Int_t i;
   double_t YM,ny;
   TLatex *l;

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   for (i=1; i<16; i++) {
      g->GetPoint(i,YM,ny);
      l = new TLatex(YM,ny+0.2,Form("%4.2f",ny));
      l->SetTextSize(0.025);
      l->SetTextFont(42);
      l->SetTextAlign(21);
      l->Paint();
   }
}