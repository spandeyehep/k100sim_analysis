import os, sys
import uproot, awkward
#import ROOT
import numpy as np
from array import array
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

import utils.randomize_nrCascadeSim as nrcRand
import utils.helper as helper


parser = OptionParser()
parser.add_option("-d", "--basedir", dest="basedir",action="store",
                      help="base directory containing k100sim step1 root file",default="/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_2M_neutrons")
parser.add_option("-f", "--inFileName", dest="inFileName",action="store",
                      help="k100sim step1 root file path",default="sim_2M_events_PuBe.root")

(options, args) = parser.parse_args()
print(options)

base_dir = options.basedir

inFile_k100 = base_dir+"/"+options.inFileName
file_k100 = uproot.open(inFile_k100)
if(not file_k100):
	print("could not open file: %s"%(inFile_k100))
	sys.exit(0)
tree_k100 = file_k100["simtree"]
if(not tree_k100):
	print("simtree does not exist in file: %s"%(file_k100))
	sys.exit(0)

event = tree_k100['EV'].array(library="np")
x = tree_k100['X3'].array(library="np") 
y = tree_k100['Y3'].array(library="np") 
z = tree_k100['Z3'].array(library="np") 
pdgid = tree_k100['Type'].array(library="np")
nCap = tree_k100['nCap'].array(library="np") 
Edep = tree_k100['D3'].array(library="np") 
DT = tree_k100['DT'].array(library="np") 

print ("] Calculating neutron capture events in silicon...")
neutron_capture_events = helper.GetNcapEvents(nCap,DT)
print ("] Found %d events out of %d events in k100sim step1 file."%(neutron_capture_events,len(nCap)))

#sys.exit(0)


inFile_nrCascade = inFile_k100.strip(".root")+"_nrCascade_randomized.root"
if(not os.path.exists(inFile_nrCascade)):
	nrCascaseFileName = inFile_k100.strip(".root")+"_nrCascade.root"
	if(not os.path.exists(nrCascaseFileName)):
		print ("] Corresponding nrCascadeSim file does not exist")
		print("] Generating nrCascadeSim file...")
		helper.generate_nrCascadeSim(neutron_capture_events,nrCascaseFileName)
		print ("] Initiating randomization...")
		nrcRand.randomize(nrCascaseFileName)
	else:
		print("] Unrandomized nrCascadeSim file found.")
		print ("] Initiating randomization...")
		nrcRand.randomize(nrCascaseFileName)
	




file_nrCascade = uproot.open(inFile_nrCascade)
if(not file_nrCascade):
	print("could not open file: %s"%(inFile_nrCascade))
	sys.exit(0)



tree_nrCascade = file_nrCascade["cascade"]
if(not tree_nrCascade):
	print("simtree does not exist in file: %s"%(file_nrCascade))
	sys.exit(0)


n = tree_nrCascade["n"].array(library="np")   # one value per event
Eg = tree_nrCascade["Eg"].array(library="np") # vector type per event

print ("] Found %d events in nrCascadeSim file."%(len(n)))
# print("len(n) = %d"%(len(n)))
# print("len(Eg) = %d"%(len(Eg)))

# randomization of nrCascadeSim events
# shuffler = np.random.permutation(len(n)) 
# n = n[shuffler]
# Eg = Eg[shuffler]

#sys.exit(0)

nCap_flag = array('h')
nGamma = array('h')
x_cap = array('d')
y_cap = array('d')
z_cap = array('d')
EGamma_builder = awkward.ArrayBuilder()
counter = 0

print ("] Combining k100sim step1 and nrCascadeSim files...")

#print (nCap[124])
# print (np.where(nCap[108] == 1))
# events = array('L')
# for i in range(len(event)):
# 	for j,it in enumerate((nCap[i])):
# 		if(it == 1):
# 			#print("Found nCap in event",i)
# 			events.append(i)

# for i in events:
# 	print(i," : ",len(np.where(nCap[i])))
print("len(n) = ",len(n))
for i in range(len(event)):
	#print ("Event = ",i)
	EGamma_temp = array('d')
	ncap_index = helper.getNcapIndex(nCap[i])	
	if(ncap_index >= 0 and DT[i][ncap_index] == 1 and counter < len(n)):
		nCap_flag.append(1)
		nGamma.append(n[counter])
		x_cap.append(x[i][ncap_index])
		y_cap.append(y[i][ncap_index])
		z_cap.append(z[i][ncap_index])
		#print ("i : counter : n : energy :: ",i,counter,n[counter],Eg[counter])
		for j,energy in enumerate(Eg[counter]):
			EGamma_temp.append(energy)
		counter += 1
		EGamma_builder.append(EGamma_temp)
	# else:
	# 	nCap_flag.append(0)
	# 	nGamma.append(0)
	# 	x_cap.append(-999.0)
	# 	y_cap.append(-999.0)
	# 	z_cap.append(-999.0)
	# EGamma_builder.append(EGamma_temp)
	




#print ("counter = ",counter)

EGamma_array = EGamma_builder.snapshot()

print("----Lengths----")
print("nCap_flag : nGamma : x_cap : y_cap : z_cap : EGamma_array :: %d : %d : %d : %d : %d : %d"%(len(nCap_flag), len(nGamma), len(x_cap), len(y_cap), len(z_cap), len(EGamma_array)))
print("---------------")

# for i,it in enumerate(EGamma_array):
# 	if(nGamma[i]>0):
# 		print ("i : n : Eg : x_cap : y_cap : z_cap :: ",i,nGamma[i],it,x_cap[i],y_cap[i],z_cap[i] )

outFileName = inFile_k100.strip(".root")+"_nrCascade_step2.root"

if(os.path.exists(outFileName)):
	os.system("rm -f %s"%(outFileName))

outFile = uproot.recreate(outFileName)
outFile["combined"] = { 
						"Event": np.arange(len(nCap_flag)),
 						"nCap_flag": nCap_flag,
 						#"n_gamma": nGamma,
 						"E_gamma": EGamma_array,
 						"x_cap": x_cap,
 						"y_cap": y_cap,
 						"z_cap": z_cap

 					   }

print ("] Step2 done.")
print ("] Output File : %s"%(outFileName))