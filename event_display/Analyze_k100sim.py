#!/usr/bin/env python
# coding: utf-8


import os, sys
import uproot, awkward
import ROOT as rt
import numpy as np
from array import array
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm





#inFile_k100 = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_2M_neutrons/sim_2M_events_PuBe.root"
#inFile_k100 = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_NaI_51M/sim_NaI_51M.root"
#inFile_k100 = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_100k_nCap_cylinder/sim_100k_nCap_cylinder.root"
inFile_k100 = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_50M_PuBe_sourceAndshields/sim_50M_PuBe_sourceAndshields_0001.root"
#inFile_k100 = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_50M_PuBe_sourceAndshields_boron_shields/sim_50M_PuBe_sourceAndshields_boronShield_V1H2.root"
file_k100 = uproot.open(inFile_k100)
if(not file_k100):
    print("could not open file: %s"%(inFile_k100))
    sys.exit(0)
    
tree_k100 = file_k100["simtree"]

if(not tree_k100):
    print("simtree does not exist in file: %s"%(file_k100))
    sys.exit(0)

outFileName = "output_"+inFile_k100.split("/")[-1]
outFile = rt.TFile.Open(outFileName,"recreate")
outFile.cd()
save_ = False
# textOnTop = rt.TLatex()
# textOnTop.SetTextSize(0.05)

print ("Reading file: %s"%(inFile_k100))

# From k100sim
EV = tree_k100["EV"].array(library="np")
P = tree_k100["P"].array(library="np")
Type = tree_k100["Type"].array(library="np")
DT = tree_k100["DT"].array(library="np")
TS = tree_k100["TS"].array(library="np")
E1 = tree_k100["E1"].array(library="np")
D3 = tree_k100["D3"].array(library="np")
X1 = tree_k100["X1"].array(library="np")
Y1 = tree_k100["Y1"].array(library="np")
Z1 = tree_k100["Z1"].array(library="np")
X3 = tree_k100["X3"].array(library="np")
Y3 = tree_k100["Y3"].array(library="np")
Z3 = tree_k100["Z3"].array(library="np")
time1 = tree_k100["time1"].array(library="np")
time3 = tree_k100["time3"].array(library="np")
nCap = tree_k100["nCap"].array(library="np")




h_k100sim_Nphotons_all = rt.TH1F("h_Nphotons_all","Number of photons from k100Sim",50,0,50)
h_pdgid_parent = rt.TH1F("h_pdgid","PDGID of parent particle",3000,0,3000)

h_parent_energy_log10 = rt.TH1F("h_parent_energy_log10","Energy of neutrons",150,-10,5)
h_parent_energy_log10.GetXaxis().SetTitle("log10(Parent Energy [MeV])")

h_parent_energy_log10_inSi = rt.TH1F("h_parent_energy_log10_inSi","Energy of neutrons traversing silicon",150,-10,5)
h_parent_energy_log10_inSi.GetXaxis().SetTitle("log10(Parent Energy [MeV])")

h_parent_energy_log10_on_nCap = rt.TH1F("h_parent_energy_log10_on_nCap","Energy of neutrons on nCap",150,-10,5)
h_parent_energy_log10_on_nCap.GetXaxis().SetTitle("log10(Parent Energy [MeV])")

h_energy_neutron_nCap = rt.TH1F("h_energy_neutron_nCap","Energy of neutrons for nCap",150,0,1.5)
h_energy_neutron_nCap.GetXaxis().SetTitle("Energy [eV]")

h_Edep_neutron = rt.TH1F("h_Edep_neutron","Energy deposited by neutrons before nCap",150,-10,5)
h_Edep_neutron.GetXaxis().SetTitle("log10(Neutron Edep [MeV])")

h_neutron_x_y = rt.TH2F("h_neutron_x_y","Position of neutrons registered in detector",120,-60,60,120,-60,60)
h_neutron_x_y.GetXaxis().SetTitle("x [mm]")
h_neutron_x_y.GetYaxis().SetTitle("y [mm]")



h_neutron_x_y_nCap = rt.TH2F("h_neutron_x_y_nCap","Position of neutron captures registered in Si crystal",120,-60,60,120,-60,60)
h_neutron_x_y_nCap.GetXaxis().SetTitle("x [mm]")
h_neutron_x_y_nCap.GetYaxis().SetTitle("y [mm]")

h_neutron_x_y_z_nCap = rt.TH3F("h_neutron_x_y_z_nCap","Position of neutron captures registered in Si crystal",120,-60,60,120,-60,60,60,10.0, 70.0)
h_neutron_x_y_z_nCap.GetXaxis().SetTitle("x [mm]")
h_neutron_x_y_z_nCap.GetYaxis().SetTitle("y [mm]")
h_neutron_x_y_z_nCap.GetYaxis().SetTitle("z [mm]")


h_neutron_x = rt.TH1F("h_neutron_x","X position of neutrons registered in the detector",120,-60,60)
h_neutron_x.GetXaxis().SetTitle("x [mm]")

h_neutron_y = rt.TH1F("h_neutron_y","y position of neutrons registered in the detector",120,-60,60)
h_neutron_y.GetXaxis().SetTitle("y [mm]")

h_neutron_z = rt.TH1F("h_neutron_z","Z position of neutrons registered in the detector",600, 10.0, 70.0)
h_neutron_z.GetXaxis().SetTitle("z [mm]")

h_neutron_dR_nCap = rt.TH1F("h_neutron_dR_nCap","dR(preSetp, postStep) for neutron capture",200, 0.0, 200.0)
h_neutron_dR_nCap.GetXaxis().SetTitle("dR [mm]")

h_photon_energy_nCapInNaI = rt.TH1F("h_photon_energy_nCapInNaI","Energy of photons produced on nCap in NaI as per Geant4",150,0,15)
h_photon_energy_nCapInNaI.GetXaxis().SetTitle("Energy [MeV]")


total_capture_events = 0
capture_in_si = 0
capture_in_NaI = 0
events_with_no_neutron = []
events_with_neutron = []
events_with_nCap_in_NaI = []
events_with_nCap_in_Si = []

for i in tqdm(range(len(EV))):
    tracks = TS[i]/100000
    tracks = tracks.astype('int32')
    neutron_tracks = tracks[Type[i] == 2112]
    neutron_parents = P[i][Type[i] == 2112]
    
    #neutron_tracks = neutron_tracks[neutron_parents == 0]
    
    
    if(len(neutron_tracks) == 0): # No neutron in the event
        events_with_no_neutron.append(i)
        continue
        
    unique_tracks, indices = np.unique(neutron_tracks, return_index=True)
    
#     if(len(unique_tracks) !=1): # Sanity check
#         print("More than one unique tracks for event = ",i)
#         break
    
    events_with_neutron.append(i) 
    neutron_energy = E1[i][Type[i] == 2112]
    neutron_parents = P[i][Type[i] == 2112]
#     if(len(neutron_energy[neutron_parents == 0]) == 0):
#         continue
    neutron_DT = DT[i][Type[i] == 2112]
    neutron_D3 = D3[i][Type[i] == 2112]
    
    neutron_X1 = X1[i][Type[i] == 2112]
    neutron_Y1 = Y1[i][Type[i] == 2112]
    neutron_Z1 = Z1[i][Type[i] == 2112]
    neutron_X3 = X3[i][Type[i] == 2112]
    neutron_Y3 = Y3[i][Type[i] == 2112]
    neutron_Z3 = Z3[i][Type[i] == 2112]
    neutron_capture = nCap[i][Type[i] == 2112]
    
    parent_neutron_en = (neutron_energy)[0]
    parent_neutron_X3 = (neutron_X3)[0]
    parent_neutron_Y3 = (neutron_Y3)[0]
    parent_neutron_Z3 = (neutron_Z3)[0]
    parent_neutron_DT = (neutron_DT)[0]
    
    h_parent_energy_log10.Fill(np.log10(parent_neutron_en)) # Use first neutron resgiterd in Si
    h_neutron_x_y.Fill(parent_neutron_X3,parent_neutron_Y3)
    h_neutron_x.Fill(parent_neutron_X3)
    h_neutron_y.Fill(parent_neutron_Y3)
    h_neutron_z.Fill(parent_neutron_Z3)
    
    if(parent_neutron_DT == 1):
        h_parent_energy_log10_inSi.Fill(np.log10(parent_neutron_en))
        
    for i_ncap in (np.where(neutron_capture == 1))[0]:
        Edep_ncap = neutron_D3[i_ncap]
        if(Edep_ncap > 0.):
            print("non zero Edep %f on nCap by neutron for event = %d"%(Edep_ncap,i))
#         if(np.sum(neutron_D3) > )
        total_capture_events += 1
#         h_parent_energy_log10_on_nCap.Fill(np.log10(neutron_energy[i_ncap]))
#         h_energy_neutron_nCap.Fill(neutron_energy[i_ncap]*1.e6)
        dx = neutron_X3[i_ncap] - neutron_X1[i_ncap]
        dy = neutron_Y3[i_ncap] - neutron_Y1[i_ncap]
        dz = neutron_Z3[i_ncap] - neutron_Z1[i_ncap]
        h_Edep_neutron.Fill(np.sum(neutron_D3))
        h_neutron_dR_nCap.Fill(np.sqrt(dx*dx + dy*dy + dz*dz))
        if(neutron_DT[i_ncap] == 1):
            capture_in_si += 1
            h_neutron_x_y_nCap.Fill(neutron_X3[i_ncap],neutron_Y3[i_ncap])
            h_neutron_x_y_z_nCap.Fill(neutron_X3[i_ncap],neutron_Y3[i_ncap],neutron_Z3[i_ncap])
            h_parent_energy_log10_on_nCap.Fill(np.log10(neutron_energy[i_ncap]))
            h_energy_neutron_nCap.Fill(neutron_energy[i_ncap]*1.e6)
            events_with_nCap_in_Si.append(i)
        else:
            capture_in_NaI += 1
            events_with_nCap_in_NaI.append(i)
    
        
print ("Total events = %d"%(len(EV)))
print ("Events with no registered neutrons in detector = %d"%(len(events_with_no_neutron)))
print ("Events with registered neutrons in NaI+Si = %d"%(len(events_with_neutron)))
#print (events_with_no_neutron)
print ("Neutrons registered in Silicon = %d"%(h_parent_energy_log10_inSi.GetEntries()))
print ("Neutrons registered in NaI = %d"%(len(events_with_neutron) - h_parent_energy_log10_inSi.GetEntries()))
print ("Total neutron capture events = %d out of %d events."%(total_capture_events,len(events_with_neutron)))
print ("Total neutron capture events in Si = %d out of %d events."%(capture_in_si,total_capture_events))
print ("Total neutron capture events in NaI = %d out of %d events."%(capture_in_NaI,total_capture_events))





###################################
####### Background study  #########
###################################

for i in events_with_nCap_in_NaI:
    tracks = TS[i]/100000
    # #type(tracks)
    ncap_index = (np.where(nCap[i] == 1))[0][0]
    tracks = tracks.astype('int32')
    parent_track = tracks[ncap_index]*100000
    parent_type = Type[i][ncap_index]
    photon_tracks = tracks[Type[i] == 22]
    photon_energy = E1[i][Type[i] == 22]
    photon_parents = P[i][Type[i] == 22]    
    photon_DT = DT[i][Type[i] == 22]
    unique_tracks, indices = np.unique(photon_tracks, return_index=True)
    photon_parents = photon_parents[indices]
    #parent_type = Type[i][indices]
    photon_DT = photon_DT[indices]
    #print (photon_energy[indices])
#     print (parent_track)
#     print (parent_type)
#     print (photon_parents)
    #break
    for energy in photon_energy[indices]:
        h_photon_energy_nCapInNaI.Fill(energy)






#background study

h_k100sim_bkg_si_Edep = rt.TH1F("h_k100sim_bkg_si_Edep","Energy deposited in Si on non-nCap events",500,0,50)
h_k100sim_bkg_si_Edep.GetXaxis().SetTitle("Energy [MeV]")

h_k100sim_bkg_si_Edep_withTimeCut = rt.TH1F("h_k100sim_bkg_si_Edep_withTimeCut","Energy deposited in Si on non-nCap events < 1ms",500,0,50)
h_k100sim_bkg_si_Edep_withTimeCut.GetXaxis().SetTitle("Energy [MeV]")

h_k100sim_bkg_NaI_Edep = rt.TH1F("h_k100sim_bkg_NaI_Edep","Energy deposited in NaI on non-nCap events",500,0,50)
h_k100sim_bkg_NaI_Edep.GetXaxis().SetTitle("Energy [MeV]")

h_k100sim_bkg_NaI_Edep_withTimeCut = rt.TH1F("h_k100sim_bkg_NaI_Edep_withTimeCut","Energy deposited in NaI on non-nCap events  < 1ms",500,0,50)
h_k100sim_bkg_NaI_Edep_withTimeCut.GetXaxis().SetTitle("Energy [MeV]")

h_k100sim_bkg_NaI_Edep_NaiNcap = rt.TH1F("h_k100sim_bkg_NaI_Edep_NaiNcap","Energy deposited in NaI on non-nCap events (Ncap in NaI)",500,0,50)
h_k100sim_bkg_NaI_Edep_NaiNcap.GetXaxis().SetTitle("Energy [MeV]")

h_tile_multiplicity = rt.TH1F("h_tile_multiplicity","Numbers of tiles with non-zero Edep",25,0,25)
h_tile_multiplicity.GetXaxis().SetTitle("Multiplicity")

h_tile_energy = rt.TH1F("h_tile_energy","Edep in individual tiles",2500,0,50)
h_tile_energy.GetXaxis().SetTitle("Energy [MeV]")

h_neutron_NaI = rt.TH1F("h_neutron_NaI","Energy of neutrons flux in NaI",1500,-10,5)
h_neutron_NaI.GetXaxis().SetTitle("log10(Neutron [MeV])")

h_neutron_NaI_nCap = rt.TH1F("h_neutron_NaI","Energy of neutrons flux in NaI on nCap",1500,-10,5)
h_neutron_NaI_nCap.GetXaxis().SetTitle("log10(Neutron [MeV])")

h_edep_NaI_time = rt.TH1F("h_edep_NaI_time","time3",355,-0.5,35)
h_edep_NaI_time.GetXaxis().SetTitle("log10(time3 [ns])")


h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity = rt.TH2F("h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity","Nai Edep vs multiplicity",500,0,50,25,0,25)
h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity.GetXaxis().SetTitle("Energy [MeV]")
h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity.GetYaxis().SetTitle("Multiplicity")



debug_count = 0
Edep_tiles = np.zeros(23)

total_event_count = 0
event_count_3_10MeV = 0
total_event_count_in_trig_window = 0
ncap_nai_event_count_in_trig_window = 0

neutron_in_NaI = 0
neutron_in_NaI_nCap = 0

events_with_NaI_Edep_10MeV_or_more = []

for i in tqdm(range(len(EV))):
    ncap_si = False
    edep_si = 0.
    edep_nai = 0.
    edep_si_tCut = 0.
    edep_nai_tCut = 0.
    ncap_nai = False
    Edep_tiles_temp = np.zeros(23)
    neutron_energy = 0.
    neutron_energy_nCap = 0.
    

    
    first_step_in_detector = True
    TS_debug = -1
    for j,capture in enumerate(nCap[i]):
        if(capture == True and DT[i][j] == 1):
            ncap_si = True
            break
        if(Type[i][j] == 2112 and DT[i][j] != 1 and first_step_in_detector):
            neutron_energy = E1[i][j]
            TS_debug = TS[i][j]
            first_step_in_detector = False
            

        if(capture == True and DT[i][j] != 1):
            ncap_nai = True
#             if(TS[i][j] != TS_debug):
#                 print("something is worng with neuton in NaI @ event %d"%(i))
#                 sys.exit(0)
            if(Type[i][j] != 2112):
                print("Neutron capture without neutron?? Check event = %d"%(i))
                sys.exit()
            else:
                neutron_energy_nCap = E1[i][j]
            
        if(DT[i][j] == 1):
            edep_si += D3[i][j]
            if(time3[i][j] < 1.e6):
                edep_si_tCut += D3[i][j]
        else:
            edep_nai += D3[i][j]
            h_edep_NaI_time.Fill(np.log10(time3[i][j]))
            if(time3[i][j] < 1.e6):
                edep_nai_tCut += D3[i][j]
                if(DT[i][j] > 2000):
                    Edep_tiles_temp[DT[i][j] - 2000 - 1] += D3[i][j]
                else:
                    Edep_tiles_temp[DT[i][j] - 1000 - 1] += D3[i][j]
                
    if(ncap_si):
        continue
    
    h_k100sim_bkg_si_Edep.Fill(edep_si)
    h_k100sim_bkg_NaI_Edep.Fill(edep_nai)
    h_k100sim_bkg_si_Edep_withTimeCut.Fill(edep_si_tCut)
    if(edep_nai_tCut > 0.):
        h_k100sim_bkg_NaI_Edep_withTimeCut.Fill(edep_nai_tCut)
        if(ncap_nai):
            h_k100sim_bkg_NaI_Edep_NaiNcap.Fill(edep_nai_tCut)
        
    Edep_tiles += Edep_tiles_temp
    
    
    
    n_tiles = np.count_nonzero(Edep_tiles_temp)
    h_tile_multiplicity.Fill(n_tiles)
    h_k100sim_bkg_NaI_Edep_withTimeCut_vs_tile_multiplicity.Fill(edep_nai_tCut,n_tiles)
    
    tile_10MeV_veto = False
    energy_in_trig_window = 0.  # This is the energy sum of panels with energies between 3 and 10 MeV each panel
    for j,energy in enumerate(Edep_tiles_temp):
        h_tile_energy.Fill(energy)
        if(energy > 10.):
            tile_10MeV_veto = True
        if(energy > 3. and energy < 10. ):
            energy_in_trig_window += energy
    
    if(edep_nai > 0.):
        total_event_count += 1
        
    if(energy_in_trig_window > .0):
        event_count_3_10MeV += 1
        if(not tile_10MeV_veto):
            total_event_count_in_trig_window += 1
            if(ncap_nai):
                ncap_nai_event_count_in_trig_window += 1
    
    if(neutron_energy > 0.):
        h_neutron_NaI.Fill(np.log10(neutron_energy))
        neutron_in_NaI += 1
    if(neutron_energy_nCap > 0.):
        h_neutron_NaI_nCap.Fill(np.log10(neutron_energy_nCap))
        neutron_in_NaI_nCap += 1
        

       
    if(edep_nai > 10.):
        events_with_NaI_Edep_10MeV_or_more.append(i)
        if(debug_count < 10):
            print("event number with Edep_NaI > 10 MeV = %d"%(i))
            debug_count += 1
    if(edep_nai_tCut > 1. and debug_count < 2 and False):
        for j,energy in enumerate(Edep_tiles_temp):
            print ("Tile : Edep :: %d : %f"%(j,energy))
        print("Event : total esum :: %d : %f "%(i,edep_nai_tCut))
        debug_count += 1


print ("Events with neutrons passing through NaI = %d"%(neutron_in_NaI))
print ("Events with neutrons passing through NaI and nCap = %d"%(neutron_in_NaI_nCap))
print ("nCap in NaI event fraction = %0.2f%%"%(100*neutron_in_NaI_nCap/neutron_in_NaI))
print ("Total event count with Edep NaI > 0. = %d"%(total_event_count))
print ("Total Event count between 3 and 10 MeV = %d"%(event_count_3_10MeV))

print ("Total Event count for (3 < Edep <  10 MeV) && (Edep_tile < 10 MeV) = %d"%(total_event_count_in_trig_window))
print ("Ncap-in-NaI Event count for (3 < Edep <  10 MeV) && (Edep_tile < 10 MeV) = %d"%(ncap_nai_event_count_in_trig_window))



# ## For Debug




i = 0

for j in range(len(TS[i])):
    if(D3[i][j] > 0 or True):
        print("i:Type:TS:P:X:Y:Z:DT:E1:D3:nCap :: %d : %d : %d : %d : %0.2f : %0.2f : %0.2f : %d : %0.3f : %0.3f : %d"%(j,Type[i][j],TS[i][j],P[i][j],X1[i][j],Y1[i][j],Z1[i][j],DT[i][j],E1[i][j],D3[i][j],nCap[i][j]))


   
tracks = TS[i]/100000
tracks = tracks.astype('int32')
neutron_tracks = tracks[Type[i] == 2112]
neutron_parents = P[i][Type[i] == 2112]
#neutron_tracks = neutron_tracks[neutron_parents == 0]
neutron_energy = E1[i][Type[i] == 2112]
neutron_DT = DT[i][Type[i] == 2112]

neutron_X1 = X1[i][Type[i] == 2112]
neutron_Y1 = Y1[i][Type[i] == 2112]
neutron_Z1 = Z1[i][Type[i] == 2112]
neutron_capture = nCap[i][Type[i] == 2112]

unique_tracks, indices = np.unique(neutron_tracks, return_index=True)

neutron_energy = neutron_energy[indices]
neutron_X1 = neutron_X1[indices]
neutron_Y1 = neutron_Y1[indices]
neutron_Z1 = neutron_Z1[indices]
neutron_capture = neutron_capture[indices]

print(Type[i])
print(P[i])
print(D3[i])
print(neutron_capture)
print(neutron_energy)


# In[41]:


debug_count = 0
h_nCap_NaI_Edep_10MeV_or_more = rt.TH1F("h_nCap_NaI_Edep_10MeV_or_more","nCap in NaI when Edep_NaI > 10 MeV",4,-1,3)
h_secondary_PID_nCap_NaI = rt.TH1F("h_secondary_PID_nCap_NaI","PID of secondaries generated on NaI nCap",400000,-200000,200000)

h_energy_gamma_nCap_NaI = rt.TH1F("h_energy_gamma_nCap_NaI","Energy of gammas on nCap NaI",1000,0,50)
h_energy_gamma_nCap_NaI.GetXaxis().SetTitle("Energy [MeV]")

h_energySum_gamma_nCap_NaI = rt.TH1F("h_energySum_gamma_nCap_NaI","EnergySum of gammas on nCap NaI",1000,0,50)
h_energySum_gamma_nCap_NaI.GetXaxis().SetTitle("Energy [MeV]")

h_energy_other_nCap_NaI = rt.TH1F("h_energy_other_nCap_NaI","Energy of all other PIDs on nCap NaI",1000,0,50)
h_energy_other_nCap_NaI.GetXaxis().SetTitle("Energy [MeV]")

h_energySum_other_nCap_NaI = rt.TH1F("h_energySum_other_nCap_NaI","Energy of all other PIDs on nCap NaI",1000,0,50)
h_energySum_other_nCap_NaI.GetXaxis().SetTitle("Energy [MeV]")

h_nGammas_nCap_NaI = rt.TH1F("h_nGammas_nCap_NaI","Nb of gammas on nCap in NaI",20,0,20)
h_nGammas_nCap_NaI.GetXaxis().SetTitle("Nb of gammas")

print ("events_with_NaI_Edep_10MeV_or_more = ",len(events_with_NaI_Edep_10MeV_or_more))
debug = True
for i in events_with_NaI_Edep_10MeV_or_more:
    ncap_si = False
    edep_si = 0.
    edep_nai = 0.
    edep_si_tCut = 0.
    edep_nai_tCut = 0.
    ncap_nai = False
    Edep_tiles_temp = np.zeros(23)
    neutron_energy = 0.
    neutron_energy_nCap = 0.
    
    first_step_in_detector = True
    TS_debug = -1
    neutron_track = 0 
    for j,capture in enumerate(nCap[i]):
        if(D3[i][j] != 1 and time3[i][j] < 1.e6):
                edep_nai_tCut += D3[i][j]
#         if(capture == True and DT[i][j] == 1):
#             ncap_si = True
#             break
#         if(Type[i][j] == 2112 and DT[i][j] != 1 and first_step_in_detector):
#             neutron_energy = E1[i][j]
#             TS_debug = TS[i][j]
#             first_step_in_detector = False
         
        if(capture == True and DT[i][j] != 1):
            ncap_nai = True
            track = int(TS[i][j]/100000)*100000
            step = TS[i][j] - (track)
            neutron_track = track
#             print ("event : TS : track : step :: %d : %d : %d : %d"%(i,TS[i][j],track,step))
#             sys.exit(0)

            if(Type[i][j] != 2112):
                print("Neutron capture without neutron?? Check event = %d"%(i))
                sys.exit()
    
    if(edep_nai_tCut > 10.):
        h_nCap_NaI_Edep_10MeV_or_more.Fill(int(ncap_nai))
        
    if(ncap_nai):
        esum_gammas = 0.
        esum_others = 0.
        n_gammas = 0
        for j,capture in enumerate(nCap[i]):
            track = int(TS[i][j]/100000)*100000
            step = TS[i][j] - (track)
            
            if(P[i][j] == neutron_track and step == 1):
                #print("event : TS : Type : E1 :: %d : %d : %d : %f"%(i, TS[i][j], Type[i][j], E1[i][j]))
                h_secondary_PID_nCap_NaI.Fill(Type[i][j])
                if(Type[i][j] == 22):
                    h_energy_gamma_nCap_NaI.Fill(E1[i][j])
                    esum_gammas += E1[i][j]
                    n_gammas += 1
                else:
                    h_energy_other_nCap_NaI.Fill(E1[i][j])
                    esum_others += E1[i][j]
        h_nGammas_nCap_NaI.Fill(n_gammas)
        if(debug == True and n_gammas == 7):
            print("Nb gamma = 7 for event %d"%(i))
            debug = False
        if(esum_gammas > 0.):
            h_energySum_gamma_nCap_NaI.Fill(esum_gammas)
        if(esum_others > 0.):
            h_energySum_other_nCap_NaI.Fill(esum_others)
            
    
            
outFile.Write()
outFile.Close()
print ("%s written."%(outFileName))






