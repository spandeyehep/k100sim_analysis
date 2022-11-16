#define Analyze_k100sim_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include "TVector3.h"
#include "TRandom3.h"
#include "TF1.h"
#include "Analyze_k100sim.h"
#include <stdlib.h>

using namespace std;


bool isFe56(vector<int> *particles) {
  for(int i = 0; i < (int) particles->size(); i++) {
    if(particles->at(i) == 260570)
      return true;
  }
  return false;
}

int main(int argc, char* argv[])
{

  if (argc < 3) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" 
    << " " << " number_of_events"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *event_arg   = argv[3];
  

  Analyze_k100sim k100sim(inputFileList, outFileName);
  

  k100sim.EventLoop(atoi(event_arg));

  return 0;
}

void Analyze_k100sim::EventLoop(int nEvents) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  

  bool DEBUG = false;

  if(DEBUG) cout<<"DEBUG ==>  DEBUG flag is ON."<<endl;

  Long64_t nbytes = 0, nb = 0;
 
  int decade = 0;
  long neutron_capture = 0;

  if(DEBUG) cout<<"DEBUG ==> Entering into event loop."<<endl;

  bool tile_detID_bug = false;
  double bottom_tile_position = 0.5*(406.0 + 20.0);
  
  if(tile_detID_bug) {
    cout<<"Tile detectorID bug is in effect!!!!"<<endl;
  }
  

  long total_nCap_events_Si = 0;
  long total_nCap_events_NaI = 0;

  long capture_on_I = 0;
  long capture_on_Na = 0;
  long capture_on_others = 0;
  long capture_on_Fe = 0;

  long total_trigger_events = 0;
  long nCap_NaI_trigger_events = 0;
  vector<long> trig_events_in_windows;
  vector<long> trig_events_in_windows_nCap_NaI;
  vector<long> trig_events_in_windows_non_nCap_NaI;
  for(float j = 3.0; j <= 14.0; j+=0.2) {
    trig_events_in_windows.push_back(0);
    trig_events_in_windows_nCap_NaI.push_back(0);
    trig_events_in_windows_non_nCap_NaI.push_back(0);

  }
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    

    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;
    
    if(DEBUG) cout<<"DEBUG ==>  Reading entry : "<<jentry<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(DEBUG) cout<<"DEBUG ==>  Analysis code starts now."<<endl;

    //cout<<"***** jentry : "<<jentry<<endl;
    int nPart = (int)P->size();
    h_nParticles->Fill(nPart);
    

    int parent = 0;
    long trackID_pre = 0;
    long stepID_pre = 0;
    double Edep = 0;
    bool Ncap = false;
    bool foundParent = false;
    vector<long> tracks;
    vector<double> gamma_En;
    vector<double> gamma_En_nonNcap;
    tracks.clear();
    gamma_En.clear();
    gamma_En_nonNcap.clear();
    double E_neutron_nCap = 0.0;
    int double_capture_count = 0;
    long TS_debug = -1;
    bool ncap_si = false;
    double edep_si = 0.;
    double edep_nai = 0.;
    double edep_si_tCut = 0.;
    double edep_nai_tCut = 0.;
    bool ncap_nai = false;
    double neutron_energy = -1.;
    double neutron_energy_nCap = 0.;
    bool first_step_in_detector = true;
    double nCap_NaI_time = 0.;
    double total_Egamma_nCap = 0.0;
    int nPhotons = 0;
    vector<unsigned int> daughters;
    daughters.clear();
    //cout<<"Entry = "<<jentry<<endl;
    //cout<<"P size = "<<P->size()<<endl;
    //cout<<"DT size = "<<DT->size()<<endl;
    // cout<<"Type size = "<<Type->size()<<endl;
    // cout<<"E1 size = "<<E1->size()<<endl;
    // cout<<"X1 size = "<<X1->size()<<endl;
    // cout<<"X3 size = "<<X3->size()<<endl;
    // cout<<"nCap size = "<<nCap->size()<<endl;
    // if(jentry > 1000)
    //   break;

    bool FeEvent = isFe56(Type);
    if(FeEvent) capture_on_Fe++;
    
    if((int)P->size() != (int)DT->size()) {
      cout<<"Unqual sizes for P and DT for event = "<<jentry<<endl;
      cout<<"P size = "<<P->size()<<endl;
      cout<<"DT size = "<<DT->size()<<endl;
      exit(0);
    }
    

    double Edep_tiles_temp[23];
    double Edep_tiles_temp_alltime[23];
    for(int i = 0; i < 23; i++) {
      Edep_tiles_temp[i] = 0.;
      Edep_tiles_temp_alltime[i] = 0.;
    }

    bool time_cut_1ms = true;
    bool time_check = false;

    for(int i = 0; i < nPart; i++) {

      h_DeltaTime_log10->Fill(log10(time1->at(i)  - time3->at(i)));
      h_preStepTime_log10->Fill(log10(time3->at(i)));
      h_postStepTime_log10->Fill(log10(time1->at(i)));
      //cout<<"time1->at(i)  - time3->at(i) = "<<(time1->at(i)  - time3->at(i))<<endl;

      // if(time1->at(i) > 1.e30 && jentry != 7) {
      //   cout<<"jentry : i : Type : TS : energy : time1 :: " << jentry << " : " << i <<" : "<< Type->at(i) <<" : "<< TS->at(i) <<" : "<< E1->at(i) <<" : "<< time1->at(i) <<endl;
      //   time_check = true;
      // }

      if(ncap_nai) {
        if(D3->at(i) > 0 && DT->at(i) != 1) {
          h_edep_NaI_time_after_nCap->Fill(log10(time1->at(i)));
        }
        if(P->at(i) == trackID_pre*1.e5) {
          long trackID_current = TS->at(i)/1.e5;
          long stepID_current = TS->at(i) - trackID_current*1.e5;
          if(stepID_current == 1) {
            daughters.push_back(Type->at(i));
            if(Type->at(i) == 22) {
              h_Egamma_on_nCap->Fill(E1->at(i));
              h_timeGamma_after_nCap->Fill(log10(time1->at(i)));
              total_Egamma_nCap += E1->at(i);
              nPhotons++;
              gamma_En.push_back(E1->at(i)*1.e3);
            }
            
          }
        }


        // if(Type->at(i) == 22) {
        //   if(P->at(i) == trackID_pre*1.e5) {
        //     long trackID_current = TS->at(i)/1.e5;
        //     long stepID_current = TS->at(i) - trackID_current*1.e5;
        //     if(stepID_current == 1) {
        //       h_Egamma_on_nCap->Fill(E1->at(i));
        //       h_timeGamma_after_nCap->Fill(log10(time1->at(i)));
        //       total_Egamma_nCap += E1->at(i);
        //       nPhotons++;
        //       gamma_En.push_back(E1->at(i)*1.e3);

        //     }
        //   }
        // }
        
      }

      if(jentry == 111) {
        cout<<" i : Type : TS : P : X : Y : Z : DT : E1 : D3 : nCap : time1 :: "<< i << " : " << Type->at(i) << " : " << TS->at(i) << " : " << P->at(i) << " : " << X1->at(i) << " : " << Y1->at(i) << " : " << Z1->at(i) << " : " << DT->at(i) << " : " << E1->at(i) << " : " << D3->at(i) << " : " << nCap->at(i) << " : " <<time1->at(i)  <<endl;
      }

      
      if(DT->at(i) > 2000) { //vertical tiles : 2001 - 2020
        Edep_tiles_temp_alltime[DT->at(i) - 2000 - 1] += D3->at(i);

      }
      else { //horizontal tiles : 1001 - 1003
        Edep_tiles_temp_alltime[DT->at(i) - 1000 - 1 + 20] += D3->at(i);
        // cout<<"DT->at(i) : DT->at(i) - 1000 - 1 :: "<<(DT->at(i))<<" : "<<(DT->at(i) - 1000 - 1 + 20)<<endl;
        // return;
      }

      if(time1->at(i) > 1.e11) continue;
      //if(time1->at(i) > 1.e6) time_cut_1ms = false;

      if(nCap->at(i) == true && DT->at(i) == 1) {
        ncap_si = true;
        total_nCap_events_Si++;
      }

      if(Type->at(i) == 2112 && DT->at(i) != 1 && first_step_in_detector){
        neutron_energy = E1->at(i);
        TS_debug = TS->at(i);
        first_step_in_detector = false;
      }

      

      if(nCap->at(i) == true && DT->at(i) != 1) {
        ncap_nai = true;
        total_nCap_events_NaI++;
        trackID_pre = TS->at(i)/1.e5;
        if(Type->at(i) != 2112) {
          cout<<"Neutron capture without neutron?? Check event = "<<jentry<<endl;
          exit(0);
        }
        else {
          neutron_energy_nCap = E1->at(i);
          nCap_NaI_time = time1->at(i);
        }
        
      }

      
      if(DT->at(i) > 2000 && abs(X1->at(i)) < bottom_tile_position && tile_detID_bug) {
        switch(DT->at(i)) {
          case 2001: DT->at(i) = 1001; break;
          case 2002: DT->at(i) = 1002; break;
          case 2003: DT->at(i) = 1003; break;
          default: break;
        }
        //cout<<"DetID bug: DT : X1 :: "<<DT->at(i)<<" : "<<X1->at(i)<<endl;
      } 
        

      

      if(time_cut_1ms) {
        if(DT->at(i) == 1)
          edep_si_tCut += D3->at(i);
        else {
          edep_nai += D3->at(i);
          h_edep_NaI_time->Fill(log10(time1->at(i)));
          if(DT->at(i) > 2000) { //vertical tiles : 2001 - 2020
            Edep_tiles_temp[DT->at(i) - 2000 - 1] += D3->at(i); // after calculation: vertical tile indices 0 to 19

          }

          else { //horizontal tiles : 1001 - 1003
            Edep_tiles_temp[DT->at(i) - 1000 - 1 + 20] += D3->at(i); // after calculation: horizontal tile indices 20 to 22
            // cout<<"DT->at(i) : DT->at(i) - 1000 - 1 :: "<<(DT->at(i))<<" : "<<(DT->at(i) - 1000 - 1 + 20)<<endl;
            // return;
          }
        }  
      }
      
      long trackID_temp = TS->at(i)/1.e5;
      //long stepID_current = TS->at(i) - trackID_current*1.e5;
      if(IsUnique(tracks, trackID_temp) && Type->at(i) == 22) {
        gamma_En_nonNcap.push_back(E1->at(i));
      }

    }

    if(time_check) return;
    
    if(total_Egamma_nCap) h_total_Egamma_on_nCap->Fill(total_Egamma_nCap);
    //if(total_Egamma_nCap > 25.0) cout<<"jentry : total_Egamma_nCap : nPhotons :: "<<jentry<<" : " <<total_Egamma_nCap << " : " << nPhotons<<endl;
      

    if(neutron_energy_nCap > 0. && ncap_nai) {
      h_total_Egamma_on_nCap_vs_log10_Eneutron->Fill(log10(neutron_energy_nCap),total_Egamma_nCap*1.e3);
    }
      

    if(DEBUG) cout<<"DEBUG ==>  End of particle track loop."<<endl;



      // if(Type->at(i) == 22) h_gamma_z->Fill(Z1->at(i));
      // if(P->at(i) == 0 && !foundParent)  {
      //   parent++;
      //   h_parent_energy_log10->Fill(log10(E1->at(i)));
      //   E_neutron_nCap = E1->at(i);
      //   foundParent = true;
      // }

      // if(nCap->at(i)) {
      //   neutron_capture++;
      //   Ncap = true;
      //   double_capture_count++;
      // }
      
      //Edep += D3->at(i);
      // // cout<<"i : TS : P : type : (x,y,z) : KE : Edep :: "<<i<<" : "<<TS->at(i)<<" : "<<P->at(i)<<" : "<<Type->at(i)<<" : ("<<X1->at(i)<<", "<<Y1->at(i)<<" ,"<<Z1->at(i)<<") : "<<E1->at(i)<<" : "<<D3->at(i)<<endl;

      
      // // if(i > 0 && false) {
      // //   trackID_pre = TS->at(i-1)/1.e5;
      // //   stepID_pre = TS->at(i-1) - trackID_pre*1.e5;
      // //   if(trackID_pre == trackID_current && stepID_current == stepID_pre+1) {
      // //     TVector3 v_current(X1->at(i),Y1->at(i),Z1->at(i));
      // //     TVector3 v_pre(X1->at(i-1),Y1->at(i-1),Z1->at(i-1));
      // //     cout<<"DeltaR = "<<v_current.DeltaR(v_pre)<<endl;
      // //   }
      // // }


    



    // if(ncap_si)
    //   continue;

    if(ncap_nai) {
      h_Ngamma_on_nCap->Fill(nPhotons);
      // if(nPhotons == 0)  {
      //   cout<<"jentry for nPhotons = 0 : "<<jentry<<endl;
      //   return;
      // }
      //cout<<"Hello 1"<<endl;
      for (int i = 0; i < (int)gamma_En.size();i++) {
        h_2D_Ngamma_EGamma_nCap_NaI->Fill(gamma_En.size(),gamma_En[i]);
        h_E_gamma_nCap_photons_inclusive->Fill(gamma_En[i]);
        if((int)gamma_En.size() < 1) h_E_gamma_nCap_photons[0]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 2) h_E_gamma_nCap_photons[1]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 3) h_E_gamma_nCap_photons[2]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 4) h_E_gamma_nCap_photons[3]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 5) h_E_gamma_nCap_photons[4]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 6) h_E_gamma_nCap_photons[5]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 7) h_E_gamma_nCap_photons[6]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 8) h_E_gamma_nCap_photons[7]->Fill(gamma_En[i]);
        else if((int)gamma_En.size() < 9) h_E_gamma_nCap_photons[8]->Fill(gamma_En[i]);


        if(neutron_energy_nCap < 1.e-7) {
          h_thermalNeutronEnergy->Fill(log10(neutron_energy_nCap));
          h_E_gamma_nCap_photons_inclusive_thermalNeutron->Fill(gamma_En[i]);
          if((int)gamma_En.size() < 1) h_E_gamma_nCap_photons[0]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 2) h_E_gamma_nCap_photons_thermalNeutron[1]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 3) h_E_gamma_nCap_photons_thermalNeutron[2]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 4) h_E_gamma_nCap_photons_thermalNeutron[3]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 5) h_E_gamma_nCap_photons_thermalNeutron[4]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 6) h_E_gamma_nCap_photons_thermalNeutron[5]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 7) h_E_gamma_nCap_photons_thermalNeutron[6]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 8) h_E_gamma_nCap_photons_thermalNeutron[7]->Fill(gamma_En[i]);
          else if((int)gamma_En.size() < 9) h_E_gamma_nCap_photons_thermalNeutron[8]->Fill(gamma_En[i]);
        }
      }
      //cout<<"Hello 2"<<endl;
      bool cap_iodine = false;
      bool cap_sodium = false;
      bool cap_others = false;
      for(int i = 0 ; i < (int)daughters.size(); i++) {
        if(daughters[i] < 1000) continue;
        if(daughters[i]%100 == 53) cap_iodine = true;
        else if(daughters[i]%100 == 11) cap_sodium = true;
        else cap_others = true;

        
      }
      //cout<<"Hello 3"<<endl;
      if(cap_iodine) {
        capture_on_I++;
        h_Ngamma_on_nCap_I127->Fill(nPhotons);
      }
      else if(cap_sodium) capture_on_Na++;
      else capture_on_others++;

      // if((cap_iodine && cap_sodium) || (cap_iodine && cap_others) || (cap_others && cap_sodium)) {
      //   // cout<<"Multiple daughters in jentry (see below) = "<<jentry<<endl;
      //   // for(int i = 0 ; i < (int)daughters.size(); i++) {
      //   //   cout<<i<<"\t"<<daughters[i]<<endl;
      //   // }
      //   // return;
      // }
      // else {
      //   if(cap_iodine) capture_on_I++;
      //   else if(cap_sodium) capture_on_Na++;
      //   else capture_on_others++;
      // }
      //cout<<"Hello 4"<<endl;
    }

    else {
      //cout<<"Hello 5"<<endl;
      for (int i = 0; i < (int)gamma_En_nonNcap.size();i++) {
        h_Egamma_on_NON_nCap->Fill(gamma_En_nonNcap[i]);
      }
    }


    if(DEBUG) cout<<"DEBUG ==> End of NaI capture condition."<<endl;

    if(neutron_energy > 0.)
      h_neutron_NaI->Fill(log10(neutron_energy));
    if(neutron_energy_nCap > 0.) {
      h_neutron_NaI_nCap->Fill(log10(neutron_energy_nCap));
      h_neutron_NaI_nCap_time->Fill(log10(nCap_NaI_time));
    }

    h_k100sim_bkg_si_Edep->Fill(edep_si);
    h_k100sim_bkg_NaI_Edep->Fill(edep_nai);
    h_k100sim_bkg_si_Edep_withTimeCut->Fill(edep_si_tCut);

    if(edep_nai_tCut > 0.) {
      h_k100sim_bkg_NaI_Edep_withTimeCut->Fill(edep_nai_tCut);
      if(ncap_nai)
        h_k100sim_bkg_NaI_Edep_NaiNcap->Fill(edep_nai_tCut);
    }

    int n_tiles = 0;
    double edep_tiles_total = 0;

    // There are two vetos for NaI Edep: 
    // 1) Reject the event if Edep in any of the tile exceeds 10 MeV,
    // 2) consider only those tiles for which  3 MeV < Edep < 10 MeV

    bool tile_10MeV_veto = false;
    bool trigger_window = false;

    vector<bool> trig_window;
    vector<bool> trig_veto;
    vector<double> upper_threshold;

    if(DEBUG) cout<<"DEBUG ==>  Start calculating triggered events."<<endl;
    
    for(float j = 3.0; j <= 14.0; j+=0.2) {
      trig_window.push_back(false);
      trig_veto.push_back(false);
      upper_threshold.push_back(j);
    }

    TF1* reso = new TF1("reso","171.48*x^(-0.504)");
    TRandom3* rand = new TRandom3();

    for(int i = 0; i < 23; i++) {
      double saturated_tile_energy = Edep_tiles_temp[i] < 12.0 ? Edep_tiles_temp[i] : 12.0;
      double resolution = (reso->Eval(saturated_tile_energy*1.e3)/100.0)*(saturated_tile_energy*1.e3);
      double smeared_tile_energy = rand->Gaus(saturated_tile_energy,(resolution/1.e3));
      h_tile_energy->Fill(Edep_tiles_temp[i]);
      h_tile_energy_allTime->Fill(Edep_tiles_temp_alltime[i]);

      if(ncap_nai) {
        h_tile_energy_onNcap->Fill(Edep_tiles_temp[i]);
      } else {
        h_tile_energy_non_Ncap->Fill(Edep_tiles_temp[i]);
      }

      if(i < 20) h_tile_energy_vertical->Fill(Edep_tiles_temp[i]);
      else h_tile_energy_horizontal->Fill(Edep_tiles_temp[i]);
      double tile_energy = Edep_tiles_temp[i];
      //double tile_energy = smeared_tile_energy > 0.? smeared_tile_energy : 0.;

      if(smeared_tile_energy > 0.)
        h_tile_energy_reso->Fill(smeared_tile_energy);
      int l = 0;
      for(float j = 3.0; j <= 14.0; j+=0.2) {
        if(tile_energy > upper_threshold[l])
          trig_veto[l] = true;
        else if(tile_energy > 3.0)
          trig_window[l] = true;
        l++;
      }
      if(tile_energy > 10.)
        tile_10MeV_veto = true;
      else if(tile_energy > 3.0)
        trigger_window = true;

      if(tile_energy > 3.)
        n_tiles++;
      edep_tiles_total += tile_energy;

    }

    if(DEBUG) cout<<"DEBUG ==>  Triggered events calculation done."<<endl;

    if(n_tiles > 0)
      h_tile_multiplicity->Fill(n_tiles);

    if(trigger_window && !tile_10MeV_veto) {
      total_trigger_events++;
      if(ncap_nai)
        nCap_NaI_trigger_events++;
    }


    // if(edep_tiles_total > 4.) {
    //   cout<<"jentry : EV :: "<<jentry<<" : "<<EV<<endl;
    // }


    int l = 0;
    for(float j = 3.0; j <= 14.0; j+=0.2) {
      if(trig_window[l] && !trig_veto[l]) {
        trig_events_in_windows[l]++;
        if(ncap_nai)
          trig_events_in_windows_nCap_NaI[l]++;
        else
          trig_events_in_windows_non_nCap_NaI[l]++;
      }

      l++;
    }

    if(DEBUG) cout<<"DEBUG ==>  Total trigger counting done."<<endl;
    //if(Ncap) cout<<"jentry = "<<jentry<<endl;
    // if(Ncap ) {
    //   cout<<"found neutron capture for event = "<<jentry<<endl;
    //   for(int i =0; i < nPart; i++) {
    //     cout<<"i : nCap->at(i) :: "<<i <<" : "<<nCap->at(i)<<endl;
    //   }
    // }
    // if(double_capture_count > 1) {
    //   cout<<double_capture_count<<" neutron captures found for event = "<<jentry<<endl;
    //   break;
    // }

    // h_parent->Fill(parent);
    // vector<double> sorted_en;
    // sorted_en.clear();
    // sorted_en = sort(gamma_En);
    // if(Ncap) { 
    //   h_e_dep_on_nCap->Fill(log10(Edep));
    //   if(foundParent) {
    //     h_parent_energy_log10_on_nCap->Fill(log10(E_neutron_nCap));
    //     //cout<<"neutron energy = "<<E_neutron_nCap<<endl;  
    //   }
    // }
    // if(sorted_en.size() != 0) {
    //   if(foundParent) h_Egamma->Fill(1.e3 * sorted_en[0]);
    //   if(Ncap) h_Egamma_on_nCap->Fill(1.e3 * sorted_en[0]);;
    // }
    
    ////////////////////////////////
    /////// For debugging //////////
    ////////////////////////////////

    // if(false) {
    //   if(!foundParent && Ncap ) {
    //     if(sorted_en.size() == 0) continue;
    //     if((sorted_en[0] <= 2.22 || sorted_en[0] > 2.23)) continue; 
    //     cout<<"######## Parent not found in event = "<<jentry<<endl;
    //     cout<<"Number of unique photons = "<<sorted_en.size()<<endl;
    //     if(sorted_en.size() > 0) cout<<"Energy of leading photon [keV] = "<<(1.e3 * sorted_en[0])<<endl;
    //     cout<<"----- Entering particle loop -----"<<endl;
    //     for(int i = 0; i < nPart; i++) {
    //       cout<<"i : TS : P : type : (x,y,z) : KE : Edep :: "<<i<<" : "<<TS->at(i)<<" : "<<P->at(i)<<" : "<<Type->at(i)<<" : ("<<X1->at(i)<<", "<<Y1->at(i)<<" ,"<<Z1->at(i)<<") : "<<E1->at(i)<<" : "<<D3->at(i)<<endl;  
          
    //     }
    //     cout<<"----- End of particle loop -----"<<endl;
    //     break;
    //   }  
    // }
    
    ////////////////////////////////

    if(P->size() != X1->size()) {
      cout<<"Unequal size. Breaking"<<endl;
      break;
    }
    
    //if(jentry > 500000) break;
    
    if(DEBUG) cout<<"DEBUG ==>  End of event loop for jentry = "<<jentry<<endl;

  } // loop over entries

  //cout<<"Neutron capture events = "<<neutron_capture<<endl;
  cout<<" Livetime for PuBe 21.8 seconds for 50 M events."<<endl;
  cout<<" //////////////////////////////// "<<endl;
  cout<<" //////// Run summary /////////// "<<endl;
  cout<<" //////////////////////////////// "<<endl;

  cout<<"Neutron capture events in silicon = "<<total_nCap_events_Si<<endl;
  cout<<"Neutron capture events in NaI = "<<total_nCap_events_NaI<<endl;
  cout<<"nCap NaI : capture on I = "<<capture_on_I<<endl;
  cout<<"nCap NaI : capture on Na = "<<capture_on_Na<<endl;
  cout<<"nCap NaI : capture on others = "<<capture_on_others<<endl;
  cout<<"Total triggerd events = "<<total_trigger_events<<endl;
  cout<<"Trigger events on nCap in NaI = "<<nCap_NaI_trigger_events<<endl;
  cout<<endl<<endl;
  cout<<"-------Rates-------"<<endl;
  double Sim_Events = (double)nEvents;
  double livetime = 21.8 * (Sim_Events/50);
  livetime = 0.218;
  cout<<"livetime for this run ("<< nEvents <<" M events) =  "<<livetime<<" seconds"<<endl;
  cout<<"Total triggerd rate = "<<(double)total_trigger_events/livetime<<" Hz"<<endl;
  cout<<"Neutron capture rate in silicon = "<<(double)total_nCap_events_Si/livetime<<" Hz"<<endl;
  cout<<"Neutron capture rate in NaI = "<<(double)total_nCap_events_NaI/livetime<<" Hz"<<endl;

  cout<<"---- For different upper thresholds ----"<<endl;
  int l = 0;
  double x[56], y[56];
  cout<<"Upper threshold [MeV] Rate [Hz]"<<endl;
  for(float j = 3.0; j <= 14.0; j+=0.2) {
    x[l] = j;
    y[l] = (double)trig_events_in_windows[l]/livetime;
    //rate_in_window->SetPoint(l,x[l],y[l]);
    //cout<<"Upper threshold [MeV] : Rate [Hz] :: "<<x[l]<<" : "<<y[l]<<endl;
    cout<<x[l]<<" "<<y[l]<<endl;
    l++;
  }

  rate_in_window = new TGraph(56,x,y);
  rate_in_window->SetMarkerStyle(kFullCircle);
  rate_in_window->SetLineColor(kBlack);
  rate_in_window->GetXaxis()->SetTitle("Upper Threshold [MeV]");
  rate_in_window->GetYaxis()->SetTitle("Trigger rate [Hz]");
  rate_in_window->Write();

  l = 0;
  cout<<"Upper threshold [MeV] Rate [Hz] on nCap in NaI"<<endl;
  for(float j = 3.0; j <= 14.0; j+=0.2) {
    x[l] = j;
    y[l] = (double)trig_events_in_windows_nCap_NaI[l]/livetime;
    //rate_in_window->SetPoint(l,x[l],y[l]);
    //cout<<"Upper threshold [MeV] : Rate [Hz] :: "<<x[l]<<" : "<<y[l]<<endl;
    cout<<x[l]<<" "<<y[l]<<endl;
    l++;
  }


  l = 0;
  cout<<"Upper threshold [MeV] Rate [Hz] on Non-nCap in NaI"<<endl;
  for(float j = 3.0; j <= 14.0; j+=0.2) {
    x[l] = j;
    y[l] = (double)trig_events_in_windows_non_nCap_NaI[l]/livetime;
    //rate_in_window->SetPoint(l,x[l],y[l]);
    //cout<<"Upper threshold [MeV] : Rate [Hz] :: "<<x[l]<<" : "<<y[l]<<endl;
    cout<<x[l]<<" "<<y[l]<<endl;
    l++;
  }
}

