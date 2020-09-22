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
  YmassNames.push_back("GenModel_YMass_200");
  YmassNames.push_back("GenModel_YMass_250");
  YmassNames.push_back("GenModel_YMass_300");
  YmassNames.push_back("GenModel_YMass_400");
  YmassNames.push_back("GenModel_YMass_500");
  YmassNames.push_back("GenModel_YMass_600");
  YmassNames.push_back("GenModel_YMass_700");
  YmassNames.push_back("GenModel_YMass_800");
  YmassNames.push_back("GenModel_YMass_900");
  YmassNames.push_back("GenModel_YMass_1000");
  YmassNames.push_back("GenModel_YMass_1200");
  YmassNames.push_back("GenModel_YMass_1400");


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
  double ym[17]={70,95,120,145,195,245,295,395,495,595,695,795,895,995,1195,1395,1595};
  Float_t YM[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};//vector of ymass points
  std::vector<float> npassym(16,0);//vector of the number of ym in each mass.
  std::vector<float> npassym_initial(16,0);//vector of the number of ym without selection
  std::vector<float> npassym_P1(16,0);
  std::vector<float> npassym_P2(16,0);
  std::vector<float> npassym_P3(16.0);
  std::vector<float> npassym_P4(16,0);
  std::vector<float> npassym_P5(16,0);
  std::vector<float> npassym_P6(16,0);
  std::vector<float> npassym_P7(16,0);
  std::vector<float> npassym_P8(16,0);
  std::vector<float> npassym_P9(16,0);
  std::vector<float> npassym_P10(16,0);
  std::vector<float> npassym_P11(16,0);
  std::vector<float> npassym_P12(16,0);
  

  std::vector<float> M_H;
  std::vector<float> M_Y;

  

  const unsigned int nLabels=10;
//  TCanvas *c1 = new TCanvas("c1","c1",3);
  
  TH1F* heve=new TH1F("heve","",nLabels,-0.5,20.5);
  TH1F* m_H = new TH1F("m_H"," ",80,100,140);
  TH1F* m_Y = new TH1F("m_Y"," ",90,140,320);
 // TH1F* H_FatJet_ParticleNetMD_probXbb = new TH1F("H_FatJet_ParticleNetMD_probXbb"," ",100,0,1);

//  TH1F* m_Y = new TH1F("m_Y"," ",)
  TH1F* YMass = new TH1F("YMass","",2,0,2);
  TH1F* Ydistri_initial = new TH1F("Ydistri_initial","",16,0,16);
  TH1F* Ydistri_P1 = new TH1F("Ydistri_P1","",16,0,16);
  TH1F* Ydistri_P2 = new TH1F("Ydistri_P2","",16,0,16);
  TH1F* Ydistri_P3 = new TH1F("Ydistri_P3","",16,0,16);
  TH1F* Ydistri_P4 = new TH1F("Ydistri_P4","",16,0,16);
  TH1F* Ydistri_P5 = new TH1F("Ydistri_P5","",16,0,16);
  TH1F* Ydistri_P6 = new TH1F("Ydistri_P6","",16,0,16);
  TH1F* Ydistri_P7 = new TH1F("Ydistri_P7","",16,0,16);
  TH1F* Ydistri_P8 = new TH1F("Ydistri_P8","",16,0,16);
  TH1F* Ydistri_P9 = new TH1F("Ydistri_P9","",16,0,16);
  TH1F* Ydistri_P10 = new TH1F("Ydistri_P10","",16,0,16);
  TH1F* Ydistri_P11 = new TH1F("Ydistri_P11","",16,0,16);
  TH1F* Ydistri_P12 = new TH1F("Ydistri_P12","",16,0,16);
  

  TH1F* Unequal_Ydistri_initial = new TH1F("Unequal_Ydistri_initial","",16,0,1600);
  TH1F* Unequal_Ydistri_P1 = new TH1F("Unequal_Ydistri_P1","",16,0,16);
  TH1F* Unequal_Ydistri_P2 = new TH1F("Unequal_Ydistri_P2","",16,0,16);
  TH1F* Unequal_Ydistri_P3 = new TH1F("Unequal_Ydistri_P3","",16,0,16);
  TH1F* Unequal_Ydistri_P4 = new TH1F("Unequal_Ydistri_P4","",16,0,16);
  TH1F* Unequal_Ydistri_P5 = new TH1F("Unequal_Ydistri_P5","",16,0,16);
  TH1F* Unequal_Ydistri_P6 = new TH1F("Unequal_Ydistri_P6","",16,0,16);
  TH1F* Unequal_Ydistri_P7 = new TH1F("Unequal_Ydistri_P7","",16,0,16);
  TH1F* Unequal_Ydistri_P8 = new TH1F("Unequal_Ydistri_P8","",16,0,16);
  TH1F* Unequal_Ydistri_P9 = new TH1F("Unequal_Ydistri_P9","",16,0,16);
  TH1F* Unequal_Ydistri_P10 = new TH1F("Unequal_Ydistri_P10","",16,0,16);
  TH1F* Unequal_Ydistri_P11 = new TH1F("Unequal_Ydistri_P11","",16,0,16);
  TH1F* Unequal_Ydistri_P12 = new TH1F("Unequal_Ydistri_P12","",16,0,16);
  
  TH1F* HBB = new TH1F("HBB","",100,0,1);
  TH1F* hbb = new TH1F("hbb","",100,0,1);

  TH1F* HDeep = new TH1F("HDeep","",100,0,1);
  TH1F* Hdeep = new TH1F("Hdeep"," ",100,0,1);

  heve->SetYTitle("Number of Events");
  heve->LabelsOption("v");
  YMass->SetYTitle("Number of YMass");
  YMass->LabelsOption("v");
  


  //store the variable into the tree called data
  TTree* variable = new TTree("variable", "variable");
  
  Float_t FatJet_Pt[2],FatJet_Eta[2],FatJet_Mass[2],FatJet_Msoftdrop[2],FatJet_BtagHbb[2],FatJet_ParticleNetMD_ProbXbb,FatJet_DeepTagMD_ZHbbvsQCD;
  Float_t FatJet_Phi[2],I_FatJet_ParticleNetMD_ProbXbb,I_FatJet_DeepTagMD_ZHbbvsQCD;
  Bool_t GenModel_YMAss_150,GenModel_YMAss_90,GenModel_YMAss_100,GenModel_YMAss_125,GenModel_YMAss_200,GenModel_YMAss_250,GenModel_YMAss_300,GenModel_YMAss_400;
  Bool_t GenModel_YMAss_500,GenModel_YMAss_600,GenModel_YMAss_700,GenModel_YMAss_800,GenModel_YMAss_900,GenModel_YMAss_1000,GenModel_YMAss_1200,GenModel_YMAss_1400;
  
  variable->Branch("FatJet_Pt",&FatJet_Pt,"FatJet_Pt[2]/F");
  variable->Branch("FatJet_Eta",&FatJet_Eta,"FatJet_Eta[2]/F");
  variable->Branch("FatJet_Phi",&FatJet_Phi,"FatJet_Phi[2]/F");
  variable->Branch("FatJet_Mass",&FatJet_Mass,"FatJet_Mass[2]/F");
  variable->Branch("FatJet_Msoftdrop",&FatJet_Msoftdrop,"FatJet_Msoftdrop[2]/F");
  variable->Branch("FatJet_BtagHbb",&FatJet_BtagHbb,"FatJet_BtagHbb[2]/F");
  variable->Branch("FatJet_ParticleNetMD_ProbXbb",&FatJet_ParticleNetMD_ProbXbb,"FatJet_ParticleNetMD_ProbXbb/F");
  variable->Branch("I_FatJet_ParticleNetMD_ProbXbb",&I_FatJet_ParticleNetMD_ProbXbb,"I_FatJet_ParticleNetMD_ProbXbb/F");
  variable->Branch("I_FatJet_DeepTagMD_ZHbbvsQCD",&I_FatJet_DeepTagMD_ZHbbvsQCD,"I_FatJet_DeepTagMD_ZHbbvsQCD/F");
  variable->Branch("FatJet_DeepTagMD_ZHbbvsQCD",&FatJet_DeepTagMD_ZHbbvsQCD,"FatJet_DeepTagMD_ZHbbvsQCD/F");
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

    if (jEntry % 5000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;
    heve->Fill(label[0],1.);
  // store the initial number of the y into the vector 
  for(unsigned int i = 0 ; i < nYmass ; i++ ){   
    if(data.GetBool(YmassNames[i].data())==true){
      npassym_initial[i] += 1;
     }
    }
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
    Float_t*  FatJet_phi = data.GetPtrFloat("FatJet_phi");
    Float_t*  FatJet_tau1 = data.GetPtrFloat("FatJet_tau1");
    Float_t*  FatJet_tau2 = data.GetPtrFloat("FatJet_tau2");
    Float_t*  FatJet_mass = data.GetPtrFloat("FatJet_mass");
    Float_t*  FatJet_msoftdrop = data.GetPtrFloat("FatJet_msoftdrop");
    Float_t*  FatJet_btagHbb = data.GetPtrFloat("FatJet_btagHbb");
    Float_t*  FatJet_deepTagMD_ZHbbvsQCD = data.GetPtrFloat("FatJet_deepTagMD_ZHbbvsQCD");
    Float_t*  FatJet_ParticleNetMD_probXbb = data.GetPtrFloat("FatJet_ParticleNetMD_probXbb");

    Bool_t GenModel_YMass_90 = data.GetBool("GenModel_YMass_90");
    Bool_t GenModel_YMass_100 = data.GetBool("GenModel_YMass_100");
    Bool_t GenModel_YMass_125 = data.GetBool("GenModel_YMass_125");
    Bool_t GenModel_YMass_150 = data.GetBool("GenModel_YMass_150");
    Bool_t GenModel_YMass_200 = data.GetBool("GenModel_YMass_200");
    Bool_t GenModel_YMass_250 = data.GetBool("GenModel_YMass_250");
    Bool_t GenModel_YMass_300 = data.GetBool("GenModel_YMass_300");
    Bool_t GenModel_YMass_400 = data.GetBool("GenModel_YMass_400");
    Bool_t GenModel_YMass_500 = data.GetBool("GenModel_YMass_500");
    Bool_t GenModel_YMass_600 = data.GetBool("GenModel_YMass_600");
    Bool_t GenModel_YMass_700 = data.GetBool("GenModel_YMass_700");
    Bool_t GenModel_YMass_800 = data.GetBool("GenModel_YMass_800");
    Bool_t GenModel_YMass_900 = data.GetBool("GenModel_YMass_900");
    Bool_t GenModel_YMass_1000 = data.GetBool("GenModel_YMass_1000");
    Bool_t GenModel_YMass_1200 = data.GetBool("GenModel_YMass_1200");
    Bool_t GenModel_YMass_1400 = data.GetBool("GenModel_YMass_1400");

    // bool_t .. = data.GetBool("");


    UInt_t nFatJet = data.GetInt("nFatJet");
    UInt_t nGoodPair = 0;
    UInt_t nGoodPairJet = 0;
    Float_t Max_pt = 0;
    Float_t sub_Max_pt = 0;
    Int_t leadingID = 0;
    Int_t subleadingID = 0;
    
  //  Float_t m_H = 0;

    Bool_t masspass = false;
    Bool_t pt_etapass_sub = false;
    Bool_t pt_etapass = false;
    Bool_t deep = false;
    Bool_t net = false;

    

    for(UInt_t ij=0; ij < nFatJet; ij++){
      
      if(FatJet_pt[ij] < 300)continue;
      if(fabs(FatJet_eta[ij]) > 2.4)continue;
      pt_etapass_sub = true;
//      Float_t tau21_i = FatJet_tau2[ij]/FatJet_tau1[ij];
//      if(tau21_i < 0.55){
//      tau21pass_sub = true;
//    }
//      if(!(FatJet_deepTagMD_ZHbbvsQCD[ij] > 0.8 && FatJet_ParticleNetMD_probXbb[ij] > 0.85))continue;
//      if(FatJet_deepTagMD_ZHbbvsQCD[ij] > 0.8){
//      deep_sub = true;
//    }
//      if(FatJet_ParticleNetMD_probXbb[ij] > 0.85){
//      net_sub = true;
//    }

      for(UInt_t jj=0; jj < ij; jj++){

	     if(FatJet_pt[jj] < 300)continue;
	     if(fabs(FatJet_eta[jj]) > 2.4)continue;
       pt_etapass = true;
//	     Float_t tau21_j = FatJet_tau2[jj]/FatJet_tau1[jj];
//	     if(tau21_j < 0.55){
//       tau21pass = true;
//     }
//       if(!(FatJet_deepTagMD_ZHbbvsQCD[jj] > 0.8 && FatJet_ParticleNetMD_probXbb[jj] > 0.85))continue;
//       if(FatJet_deepTagMD_ZHbbvsQCD[jj] > 0.8){
//       deep = true;
//     }
//       if(FatJet_ParticleNetMD_probXbb[jj] > 0.85){
//       net = true;
//     }
       

	// eta difference
	     if(fabs(FatJet_eta[ij]-FatJet_eta[jj]) > 1.3)continue;
       
       if(FatJet_pt[jj] > Max_pt){
        Max_pt = FatJet_pt[jj];
        leadingID = jj;
      }
/*       if(jj == leadingID){
        FatJet_tau21_j = tau21_j;
      }
*/	
	     nGoodPair++;
	
      } // end of inner jet loop

      if(FatJet_pt[ij] > sub_Max_pt && sub_Max_pt < Max_pt){
        sub_Max_pt = FatJet_pt[ij] ;
        subleadingID = ij;
      
      }
         
    } // end of outer jet loop


    if(pt_etapass_sub){
      nPass[1]++;
      heve->Fill(label[2],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P1[iymass]+=1;
      }
    }
  }  
    if(pt_etapass){
      nPass[2]++;
      heve->Fill(label[3],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P2[iymass]+=1;
      }
    }
  }
    if(nGoodPair<1)continue;
    nPass[3]++;
    heve->Fill(label[4],1.);
    for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P3[iymass]+=1;
      }
    }

    if(105 < FatJet_msoftdrop[leadingID] && FatJet_msoftdrop[leadingID] < 135 && 150 < FatJet_msoftdrop[subleadingID]){
      M_H.push_back(FatJet_msoftdrop[leadingID]);
      M_Y.push_back(FatJet_msoftdrop[subleadingID]);
      nGoodPairJet++;
      masspass = true;
      I_FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[subleadingID] ; 
      HBB->Fill(FatJet_ParticleNetMD_probXbb[subleadingID]);
      I_FatJet_DeepTagMD_ZHbbvsQCD = FatJet_deepTagMD_ZHbbvsQCD[subleadingID];
      HDeep->Fill(FatJet_deepTagMD_ZHbbvsQCD[subleadingID]);

      if(FatJet_deepTagMD_ZHbbvsQCD[subleadingID] > 0.8){
        deep = true;
        FatJet_DeepTagMD_ZHbbvsQCD = FatJet_deepTagMD_ZHbbvsQCD[subleadingID];
        Hdeep->Fill(FatJet_deepTagMD_ZHbbvsQCD[subleadingID]);

      }
      if(FatJet_ParticleNetMD_probXbb[subleadingID] > 0.85){
        net = true;
        FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[subleadingID];
        hbb->Fill(FatJet_ParticleNetMD_probXbb[subleadingID]);
      }
    }

    if(105 < FatJet_msoftdrop[subleadingID] && FatJet_msoftdrop[subleadingID] < 135 && 150 < FatJet_msoftdrop[leadingID]){
      M_H.push_back(FatJet_msoftdrop[subleadingID]);
      M_Y.push_back(FatJet_msoftdrop[leadingID]);
      nGoodPairJet++;
      masspass = true;
      I_FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[leadingID] ; 
      HBB->Fill(FatJet_ParticleNetMD_probXbb[leadingID]);
      I_FatJet_DeepTagMD_ZHbbvsQCD = FatJet_deepTagMD_ZHbbvsQCD[leadingID];
      HDeep->Fill(FatJet_deepTagMD_ZHbbvsQCD[leadingID]);

      if(FatJet_deepTagMD_ZHbbvsQCD[leadingID] > 0.8){
        deep = true;
        FatJet_DeepTagMD_ZHbbvsQCD = FatJet_deepTagMD_ZHbbvsQCD[leadingID];
        Hdeep->Fill(FatJet_deepTagMD_ZHbbvsQCD[leadingID]);
      }
      if(FatJet_ParticleNetMD_probXbb[leadingID] > 0.85){
        net = true;
        FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[leadingID];
        hbb->Fill(FatJet_ParticleNetMD_probXbb[leadingID]);
      }
    }



// store the leading & subleading information into the vector(tree)

    FatJet_Pt[0] = FatJet_pt[leadingID];
    FatJet_Pt[1] = FatJet_pt[subleadingID];
    FatJet_Eta[0] = FatJet_eta[leadingID];
    FatJet_Eta[1] = FatJet_eta[subleadingID];
    FatJet_Phi[0]=FatJet_phi[leadingID];
    FatJet_Phi[1]=FatJet_phi[subleadingID];  
    
    FatJet_Mass[0] = FatJet_mass[leadingID];
    FatJet_Mass[1] = FatJet_mass[subleadingID];

    FatJet_Msoftdrop[0] = FatJet_msoftdrop[leadingID];
    FatJet_Msoftdrop[1] = FatJet_msoftdrop[subleadingID];
    FatJet_BtagHbb[0] = FatJet_btagHbb[leadingID];
    FatJet_BtagHbb[1] = FatJet_btagHbb[subleadingID];
    
    
    if(masspass == true){
      nPass[4]++;
      heve->Fill(label[5],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P4[iymass]+=1;
        
      }
    }
  }
    
    if(deep == true){
      nPass[5]++;
      heve->Fill(label[6],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P5[iymass]+=1;
      }
    }
  }
    if(net == true){
      nPass[6]++;
      heve->Fill(label[7],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P6[iymass]+=1;
      }
    }
  }
    if(deep == true && net == true){
      nPass[7]++;
      heve->Fill(label[8],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P7[iymass]+=1;
      }
    }
  }
  
    

/*
    // get number of each value of Ymass 
    int pairY = 0 ;

    for(unsigned int iymass=0 ; iymass < nYmass; iymass++ ){
      if(data.GetBool(YmassNames[iymass].data())==true){
        npassym[iymass]+=1;
        pairY +=1;
      }
    }

//    cout << pairY <<endl;

    Ypair->Fill(pairY);

    if(pairY>1) YY +=1;

*/
    // store the data(bool) into the branch: in these data it might not pass thought the mass selection

    GenModel_YMAss_90  = GenModel_YMass_90;
    GenModel_YMAss_100 = GenModel_YMass_100;
    GenModel_YMAss_125 = GenModel_YMass_125;
    GenModel_YMAss_150 = GenModel_YMass_150;
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

  
 // Float_t ny[16] = {npassym[0],npassym[1],npassym[2],npassym[3],npassym[4],npassym[5],npassym[6],npassym[7],npassym[8],npassym[9],npassym[10],npassym[11],npassym[12],npassym[13],npassym[14],npassym[15]};
  
  
  Float_t nytotal = 0 ;
  for(int i = 0 ; i < 16 ; i++){
    nytotal += npassym_P7[i];
  }
  cout << "nytotal = " << nytotal << endl;
  cout << "Number of total events = " << nTotal << endl;

  for(int i=0;i<20;i++)
    if(nPass[i]>0)cout << "nPass["<<i<<"]= " << nPass[i] << endl;
/* 
  for(auto x : M_Y){
    cout << x << endl;
  }

*/  

  for(int i = 0 ; i < 16 ; i++){
    for(int k = 0 ; k < npassym_initial[i] ; k++ ){
      Ydistri_initial->Fill(i);
      Ydistri_initial->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_initial->Fill(YM[i]);
      Unequal_Ydistri_initial->SetBins(16,ym);
    }
    for(int j = 0 ; j < npassym_P1[i] ; j++){
      Ydistri_P1->Fill(i); 
      Ydistri_P1->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P1->Fill(YM[i]);
      Unequal_Ydistri_P1->SetBins(16,ym);
//    Ydistri_P1->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
    }

    for(int z = 0 ; z < npassym_P2[i] ; z++ ){
      Ydistri_P2->Fill(i);
      Ydistri_P2->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P2->Fill(YM[i]);
      Unequal_Ydistri_P2->SetBins(16,ym);
    }
    for(int a = 0 ; a < npassym_P3[i] ; a++){
      Ydistri_P3->Fill(i);
      Ydistri_P3->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P3->Fill(YM[i]);
      Unequal_Ydistri_P3->SetBins(16,ym);
    }
    for(int b = 0 ; b < npassym_P4[i] ; b++){
      Ydistri_P4->Fill(i);
      Ydistri_P4->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P4->Fill(YM[i]);
      Unequal_Ydistri_P4->SetBins(16,ym);
      m_H->Fill( M_H[i] );
      m_Y->Fill( M_Y[i] );
    }
    for(int c = 0 ; c < npassym_P5[i] ; c++){
      Ydistri_P5->Fill(i);
      Ydistri_P5->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P5->Fill(YM[i]);
      Unequal_Ydistri_P5->SetBins(16,ym);
    }
    for(int d = 0 ; d < npassym_P6[i] ; d++){
      Ydistri_P6->Fill(i);
      Ydistri_P6->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P6->Fill(YM[i]);
      Unequal_Ydistri_P6->SetBins(16,ym);
    }
    for(int e = 0 ; e < npassym_P7[i] ; e++){
      Ydistri_P7->Fill(i);
      Ydistri_P7->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P7->Fill(YM[i]);
      Unequal_Ydistri_P7->SetBins(16,ym);
    }

  /*  for(int f = 0 ; f < npassym_P8[i] ; f++){
      Ydistri_P8->Fill(i);
      Ydistri_P8->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P8->Fill(YM[i]);
      Unequal_Ydistri_P8->SetBins(16,ym);
    }
    for(int g = 0 ; g < npassym_P9[i] ; g++){
      Ydistri_P9->Fill(i);
      Ydistri_P9->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P9->Fill(YM[i]);
      Unequal_Ydistri_P9->SetBins(16,ym);
    }
    for(int h = 0 ; h < npassym_P10[i] ; h++){
      Ydistri_P10->Fill(i);
      Ydistri_P10->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P10->Fill(YM[i]);
      Unequal_Ydistri_P10->SetBins(16,ym);
    }
    for(int l = 0 ; l < npassym_P11[i] ; l++){
      Ydistri_P11->Fill(i);
      Ydistri_P11->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P11->Fill(YM[i]);
      Unequal_Ydistri_P11->SetBins(16,ym);
    }
    for(int m = 0 ; m < npassym_P12[i] ; m++){
      Ydistri_P12->Fill(i);
      Ydistri_P12->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
      Unequal_Ydistri_P12->Fill(YM[i]);
      Unequal_Ydistri_P12->SetBins(16,ym);
    }
 */   
  } 




  for(int j=0 ; j<16; j++){
    cout << "npassym_P7["<<j<<"]= " << npassym_P7[j] << endl;

  }

  //tlatex at heve 
/*
  heve->Draw("histo text 0");
  heve->Draw("same");
  c1->Update();
  TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  TLatex *myt = new TLatex(0,0,Form("Eff(cut2/total) = %1.4f ",(b/a)));
  myt->SetTextFont(42);
  myt->SetTextSize(0.03);
  myt->SetTextColor(kRed);
  listOfLines->Add(myt);
  heve->SetStats(0);
  c1->Modified();
  //c1->Update();
  c1->Print("mypdf.pdf");

  //tlatex at ny distri

  Ydistri->Draw("histo text 0");
  Ydistri->Draw("same");
  c1->Update();
  TPaveStats *ps1 = (TPaveStats*)c1->GetPrimitive("stats");
  ps1->SetName("mystats");
  TList *listOfLines1 = ps1->GetListOfLines();x
  TLatex *myt1 = new TLatex(0,0,Form("Eff(totalny/nPass[1]) = %1.4f ", nytotal/b));
  myt1->SetTextFont(42);
  myt1->SetTextSize(0.025);
  myt1->SetTextColor(kRed);
  listOfLines1->Add(myt1);
  Ydistri->SetStats(0);
  c1->Modified();
  //c1->Update();
  c1->Print("mypdf1.pdf"); 

*/
/*  R_yt->Draw("same");
  c1->Update();
  TPaveStats *ps2 = (TPaveStats*)c1->GetPrimitive("stats");
  ps2->SetName("mystats");
  TList *listOfLines2 = ps2->GetListOfLines();
  TLatex *myt2 = new TLatex(0,0,Form("Eff(totalny/nPass[1]) = %1.4f ", nytotal/b));
  myt1->SetTextFont(42);
  myt1->SetTextSize(0.025);
  myt1->SetTextColor(kRed);
  listOfLines1->Add(myt1);
  Ydistri->SetStats(0);
  c1->Modified();
*/

 /* //graph
  TGraph* gr = new TGraph(16,YM,ny);
  gr->SetTitle("Ymass_distribution");
  gr->SetMarkerStyle(20);
  TExec *ex = new TExec("ex","drawtext();");
  gr->GetListOfFunctions()->Add(ex);
  gr->Draw("ACP");
*/ 

  //heve->Draw("hist text 0");
  // writing example output file
  TFile* outFile = new TFile(outputFileName.data(),"recreate");
//  gr->Write();
  Ydistri_initial->Write();
  Ydistri_P1->Write();
  Ydistri_P2->Write();
  Ydistri_P3->Write();
  Ydistri_P4->Write();
  Ydistri_P5->Write();
  Ydistri_P6->Write();
  Ydistri_P7->Write();
  Ydistri_P8->Write();


  Unequal_Ydistri_initial->Write();
  Unequal_Ydistri_P1->Write();
  Unequal_Ydistri_P2->Write();
  Unequal_Ydistri_P3->Write();
  Unequal_Ydistri_P4->Write();
  Unequal_Ydistri_P5->Write();
  Unequal_Ydistri_P6->Write();
  Unequal_Ydistri_P7->Write();
  Unequal_Ydistri_P8->Write();
  
  heve->Write();
  m_H->Write();
  m_Y->Write();
  HBB->Write();
  hbb->Write();
  HDeep->Write();
  Hdeep->Write();
  YMass->Write();
  variable->Write();
  outFile->Close();
  
  
  
  
 } // end of macro
 
  


/*
//text of Tgraph
void drawtext()
  {
   Int_t i;
   double_t YM,ny;
   TLatex *l;

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   for (i=0; i<16; i++) {
      g->GetPoint(i,YM,ny);
      l = new TLatex(YM,ny+0.6,Form("%4.2f",ny));
      l->SetTextSize(0.025);
      l->SetTextFont(42);
      l->SetTextAlign(21);
      l->Paint();
   }
} */
