#ifndef Analyze_k100sim_H
#define Analyze_k100sim_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TGraph.h"

class Analyze_k100sim : public NtupleVariables{

 public:
  Analyze_k100sim(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~Analyze_k100sim();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(int nEvents);
  void     BookHistogram(const char *);
  
  TFile *oFile;
  TH1F* h_nParticles;
  TH1F* h_parent;
  TH1F* h_gamma_z;
  TH1F* h_parent_energy_log10;
  TH1F* h_parent_energy_log10_on_nCap;
  TH1F* h_Egamma_on_nCap;
  TH1F* h_Egamma;
  TH1F* h_e_dep_on_nCap;
  TH1F* h_Ngamma_on_nCap;

  TH1F* h_Ngamma_on_nCap_I127;
  
  TH1F* h_Egamma_on_NON_nCap;

  TH2F* h_2D_Ngamma_EGamma_nCap_NaI;

  TH1F* h_DeltaTime_log10;

  TH1F* h_k100sim_bkg_si_Edep;
  TH1F* h_k100sim_bkg_si_Edep_withTimeCut;
  TH1F* h_k100sim_bkg_NaI_Edep;
  TH1F* h_k100sim_bkg_NaI_Edep_withTimeCut;
  TH1F* h_k100sim_bkg_NaI_Edep_NaiNcap;
  TH1F* h_tile_multiplicity;
  TH1F* h_tile_energy;
  TH1F* h_neutron_NaI;
  TH1F* h_neutron_NaI_nCap;
  TH1F* h_edep_NaI_time;
  TH2F* h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity;
  TH1F* h_neutron_NaI_nCap_time;
  
  TGraph* rate_in_window;

  TH1F* h_preStepTime_log10;
  TH1F* h_postStepTime_log10;

  TH1F* h_timeGamma_after_nCap;
  TH1F* h_edep_NaI_time_after_nCap;

  TH1F* h_tile_energy_reso;
  TH1F* h_total_Egamma_on_nCap;

  TH2F* h_total_Egamma_on_nCap_vs_log10_Eneutron;
  TH1F* h_tile_energy_allTime;

  TH1F* h_tile_energy_onNcap;
  TH1F* h_tile_energy_non_Ncap;

  TH1F* h_tile_energy_vertical;
  TH1F* h_tile_energy_horizontal;

  TH1F* h_E_gamma_nCap_photons[9];
  TH1F* h_E_gamma_nCap_photons_inclusive;

  

  TH1F* h_thermalNeutronEnergy;
  TH1F* h_E_gamma_nCap_photons_thermalNeutron[9];
  TH1F* h_E_gamma_nCap_photons_inclusive_thermalNeutron;




  // TDirectory *d_Layer1;
  // TH1F *h_ADChg[128];
};
#endif

#ifdef Analyze_k100sim_cxx

void Analyze_k100sim::BookHistogram(const char *outFileName) {

  char* hname = new char[200];

//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  h_nParticles = new TH1F("h_nParticles","h_nParticles",100,0,100);
  h_parent = new TH1F("h_parent","h_parent",20,0,20);
  h_gamma_z = new TH1F("h_gamma_z","h_gamma_z",550,0,55);
  h_parent_energy_log10 = new TH1F("h_parent_energy_log10","h_parent_energy_log10",150,-10,5);
  h_parent_energy_log10->GetXaxis()->SetTitle("log10(Energy [MeV])");
  h_parent_energy_log10_on_nCap = new TH1F("h_parent_energy_log10_on_nCap","h_parent_energy_log10_on_nCap",150,-10,5);
  h_parent_energy_log10_on_nCap->GetXaxis()->SetTitle("log10(Energy [MeV])");
  h_Egamma = new TH1F("h_Egamma","h_Egamma",1000,0,10000);
  h_Egamma->GetXaxis()->SetTitle("Energy [keV]");
  h_Egamma_on_nCap = new TH1F("h_Egamma_on_nCap","h_Egamma_on_nCap",150,0,15);
  h_Egamma_on_nCap->GetXaxis()->SetTitle("Energy [MeV]");

  
  for(int i = 0; i < 9; i++) {
    sprintf(hname,"h_E_gamma_nCap_photons_%d",i);
    h_E_gamma_nCap_photons_thermalNeutron[i] = new TH1F(hname,hname,1500,0,15000);
  }
  
  h_E_gamma_nCap_photons_inclusive = new TH1F("h_E_gamma_nCap_photons_inclusive","h_E_gamma_nCap_photons_inclusive",1500,0,15000);



  h_thermalNeutronEnergy = new TH1F("h_thermalNeutronEnergy","h_thermalNeutronEnergy",150,-10,5);
  h_thermalNeutronEnergy->GetXaxis()->SetTitle("log10(Energy [MeV])");

  for(int i = 0; i < 9; i++) {
    sprintf(hname,"h_E_gamma_nCap_photons_thermalNeutron_%d",i);
    h_E_gamma_nCap_photons[i] = new TH1F(hname,hname,1500,0,15000);
  }

  h_E_gamma_nCap_photons_inclusive_thermalNeutron = new TH1F("h_E_gamma_nCap_photons_inclusive_thermalNeutron","h_E_gamma_nCap_photons_inclusive_thermalNeutron",1500,0,15000);

  h_total_Egamma_on_nCap = new TH1F("h_total_Egamma_on_nCap","h_Egamma_on_nCap",500,0,50);
  h_total_Egamma_on_nCap->GetXaxis()->SetTitle("Energy [MeV]");

  h_e_dep_on_nCap = new TH1F("h_e_dep_on_nCap","h_e_dep_on_nCap",150,-15,5);
  h_e_dep_on_nCap->GetXaxis()->SetTitle("log10(Energy [MeV])");

  h_Ngamma_on_nCap = new TH1F("h_Ngamma_on_nCap","h_Ngamma_on_nCap",50,0,50);
  h_Ngamma_on_nCap->GetXaxis()->SetTitle("# of photons on nCap in NaI");

  h_Ngamma_on_nCap_I127 = new TH1F("h_Ngamma_on_nCap_I127","h_Ngamma_on_nCap_I127",50,0,50);
  h_Ngamma_on_nCap_I127->GetXaxis()->SetTitle("# of photons on nCap in NaI");

  h_2D_Ngamma_EGamma_nCap_NaI = new TH2F("h_2D_Ngamma_EGamma_nCap_NaI","nGamma vs energy on nCap in NaI",50,0,50,2000,0,20000);
  h_2D_Ngamma_EGamma_nCap_NaI->GetXaxis()->SetTitle("# of photons");
  h_2D_Ngamma_EGamma_nCap_NaI->GetYaxis()->SetTitle("Energy of photon [keV]");

  h_Egamma_on_NON_nCap = new TH1F("h_Egamma_on_NON_nCap","Gamma energy on NON nCap events",500,0,50);
  h_Egamma_on_NON_nCap->GetXaxis()->SetTitle("Energy [MeV]");

  h_k100sim_bkg_si_Edep = new TH1F("h_k100sim_bkg_si_Edep","Energy deposited in Si on non-nCap events",500,0,50);
  h_k100sim_bkg_si_Edep->GetXaxis()->SetTitle("Energy [MeV]");

  h_k100sim_bkg_si_Edep_withTimeCut = new TH1F("h_k100sim_bkg_si_Edep_withTimeCut","Energy deposited in Si on non-nCap events < 1ms",500,0,50);
  h_k100sim_bkg_si_Edep_withTimeCut->GetXaxis()->SetTitle("Energy [MeV]");

  h_k100sim_bkg_NaI_Edep = new TH1F("h_k100sim_bkg_NaI_Edep","Energy deposited in NaI on non-nCap events",500,0,50);
  h_k100sim_bkg_NaI_Edep->GetXaxis()->SetTitle("Energy [MeV]");


  h_k100sim_bkg_NaI_Edep_withTimeCut = new TH1F("h_k100sim_bkg_NaI_Edep_withTimeCut","Energy deposited in NaI on non-nCap events  < 1ms",500,0,50);
  h_k100sim_bkg_NaI_Edep_withTimeCut->GetXaxis()->SetTitle("Energy [MeV]");

  h_k100sim_bkg_NaI_Edep_NaiNcap = new TH1F("h_k100sim_bkg_NaI_Edep_NaiNcap","Energy deposited in NaI on non-nCap events (Ncap in NaI)",500,0,50);
  h_k100sim_bkg_NaI_Edep_NaiNcap->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_multiplicity = new TH1F("h_tile_multiplicity","Numbers of tiles with non-zero Edep",25,0,25);
  h_tile_multiplicity->GetXaxis()->SetTitle("Multiplicity");

  h_tile_energy = new TH1F("h_tile_energy","Edep in individual tiles",200,0,20);
  h_tile_energy->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_onNcap = new TH1F("h_tile_energy_onNcap","Edep in individual tiles on NaI-nCap",200,0,20);
  h_tile_energy_onNcap->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_non_Ncap = new TH1F("h_tile_energy_non_Ncap","Edep in individual tiles on NON NaI-nCap",200,0,20);
  h_tile_energy_non_Ncap->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_vertical = new TH1F("h_tile_energy_vertical","Edep in individual tiles [vertical]",200,0,20);
  h_tile_energy_vertical->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_horizontal = new TH1F("h_tile_energy_horizontal","Edep in individual tiles [horizontal]",200,0,20);
  h_tile_energy_horizontal->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_allTime = new TH1F("h_tile_energy_allTime","Edep in individual tiles for all time",200,0,20);
  h_tile_energy_allTime->GetXaxis()->SetTitle("Energy [MeV]");

  h_tile_energy_reso = new TH1F("h_tile_energy_reso","Edep in individual tiles",200,0,20);
  h_tile_energy_reso->GetXaxis()->SetTitle("Energy [MeV]");

  h_neutron_NaI = new TH1F("h_neutron_NaI","Energy of neutrons flux in NaI",1500,-10,5);
  h_neutron_NaI->GetXaxis()->SetTitle("log10(Neutron [MeV])");

  h_neutron_NaI_nCap = new TH1F("h_neutron_NaI_nCap","Energy of neutrons flux in NaI on nCap",1500,-10,5);
  h_neutron_NaI_nCap->GetXaxis()->SetTitle("log10(Neutron [MeV])");

  h_neutron_NaI_nCap_time = new TH1F("h_neutron_NaI_nCap_time","Neutron capture time in NaI",355,-0.5,35);
  h_neutron_NaI_nCap_time->GetXaxis()->SetTitle("log10(time3 [ns])");

  h_edep_NaI_time = new TH1F("h_edep_NaI_time","time3",355,-0.5,35);
  h_edep_NaI_time->GetXaxis()->SetTitle("log10(time3 [ns])");

  h_edep_NaI_time_after_nCap = new TH1F("h_edep_NaI_time_after_nCap","time3",355,-0.5,35);
  h_edep_NaI_time_after_nCap->GetXaxis()->SetTitle("log10(time3 [ns])");

  h_timeGamma_after_nCap = new TH1F("h_timeGamma_after_nCap","time3",355,-0.5,35);
  h_timeGamma_after_nCap->GetXaxis()->SetTitle("log10(time3 [ns])");

  h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity = new TH2F("h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity","Nai Edep vs multiplicity",500,0,50,25,0,25);
  h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity->GetXaxis()->SetTitle("Energy [MeV]");
  h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity->GetYaxis()->SetTitle("Multiplicity");


  h_DeltaTime_log10 = new TH1F("h_DeltaTime_log10","time3 - time1",405,-0.5,40);
  h_DeltaTime_log10->GetXaxis()->SetTitle("log10(time3 - time1) [ns]");

  h_postStepTime_log10 = new TH1F("h_postStepTime_log10","postStep global time",405,-0.5,40);
  h_postStepTime_log10->GetXaxis()->SetTitle("log10(postStepGlobalTime) [ns]");

  h_preStepTime_log10 = new TH1F("h_preStepTime_log10","preStep global time",405,-0.5,40);
  h_preStepTime_log10->GetXaxis()->SetTitle("log10(preStepGlobalTime) [ns]");

  h_total_Egamma_on_nCap_vs_log10_Eneutron = new TH2F("h_total_Egamma_on_nCap_vs_log10_Eneutron","Neutron energy vs total gamma energy",1500,-10,5,3000,0,30000);
  h_total_Egamma_on_nCap_vs_log10_Eneutron->GetXaxis()->SetTitle("log10(neutron Energy [MeV])");
  h_total_Egamma_on_nCap_vs_log10_Eneutron->GetYaxis()->SetTitle("total gamma [keV]");

  // rate_in_window = new TGraph();
  // rate_in_window->SetMarkerStyle(kFullCircle);
  // rate_in_window->SetLineColor(kBlack);
  // rate_in_window->GetXaxis()->SetTitle("Upper Threshold [MeV]");
  // rate_in_window->GetYaxis()->SetTitle("Trigger rate [Hz]");



  // d_Layer1 = oFile->mkdir("Layer1");
  // char hname[200];
  // d_Layer1-> cd();
  // for(int ii=0; ii<128; ii++) {
  //   sprintf(hname, "ADC_HG_%d", ii+1);
  //   h_ADChg[ii] = new TH1F(hname, hname, 100, 0, 400);
  // }
}


Analyze_k100sim::Analyze_k100sim(const TString &inputFileList, const char *outFileName, const char* dataset) {

  TChain *tree = new TChain("simtree");
  

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  NtupleVariables::Init(tree);

  BookHistogram(outFileName);
  
}

Bool_t Analyze_k100sim::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
    
  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

Long64_t Analyze_k100sim::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  return centry;

  // if (centry!=0)
  //   return centry;
  // else return -1;
}

Analyze_k100sim::~Analyze_k100sim() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif
