import os, sys
import uproot, awkward
#import ROOT
import numpy as np
from array import array
from mpl_toolkits import mplot3d
import numpy as np
from tqdm import tqdm
from optparse import OptionParser

import utils.helper as helper


parser = OptionParser()
parser.add_option("-d", "--basedir", dest="basedir",action="store",
                      help="base directory containing k100sim step1 root file",
                      default="/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_2M_neutrons")
parser.add_option("-f", "--step1FileName", dest="step1file",action="store",
                      help="k100sim step1 root file",default="inin.root")
parser.add_option("-c", "--step3FileName", dest="step3file",action="store",
                      help="k100sim root generated at step3",default="inin.root")


(options, args) = parser.parse_args()
print(options)

base_dir = options.basedir

inFile_k100_step1 = base_dir+"/"+options.step1file
inFile_nrCascade = inFile_k100_step1.strip(".root")+"_nrCascade_randomized.root"
inFile_k100_step3 = options.step3file


file_k100_step1 = uproot.open(inFile_k100_step1)
if(not file_k100_step1):
	print("could not open file: %s"%(inFile_k100_step1))
	sys.exit(0)

file_k100_step3 = uproot.open(inFile_k100_step3)
if(not file_k100_step3):
	print("could not open file: %s"%(inFile_k100_step3))
	sys.exit(0)

file_nrCascade = uproot.open(inFile_nrCascade)
if(not file_nrCascade):
	print("could not open file: %s"%(inFile_nrCascade))
	sys.exit(0)



tree_k100_step1 = file_k100_step1["simtree"]
if(not tree_k100_step1):
	print("simtree does not exist in file: %s"%(file_k100_step1))
	sys.exit(0)

EV_step1 = tree_k100_step1['EV'].array(library="np")
DT_step1 = tree_k100_step1['DT'].array(library="np")
TS_step1 = tree_k100_step1['TS'].array(library="np")
P_step1 = tree_k100_step1['P'].array(library="np")
Type_step1 = tree_k100_step1['Type'].array(library="np")
E1_step1 = tree_k100_step1['E1'].array(library="np")
D3_step1 = tree_k100_step1['D3'].array(library="np")
X1_step1 = tree_k100_step1['X1'].array(library="np") 
Y1_step1 = tree_k100_step1['Y1'].array(library="np") 
Z1_step1 = tree_k100_step1['Z1'].array(library="np") 
X3_step1 = tree_k100_step1['X3'].array(library="np") 
Y3_step1 = tree_k100_step1['Y3'].array(library="np") 
Z3_step1 = tree_k100_step1['Z3'].array(library="np") 
PX1_step1 = tree_k100_step1['PX1'].array(library="np") 
PY1_step1 = tree_k100_step1['PY1'].array(library="np") 
PZ1_step1 = tree_k100_step1['PZ1'].array(library="np") 
PX3_step1 = tree_k100_step1['PX3'].array(library="np") 
PY3_step1 = tree_k100_step1['PY3'].array(library="np") 
PZ3_step1 = tree_k100_step1['PZ3'].array(library="np") 
time1_step1 = tree_k100_step1['time1'].array(library="np")
time3_step1 = tree_k100_step1['time3'].array(library="np")
nCap_step1 = tree_k100_step1['nCap'].array(library="np")


tree_nrCascade = file_nrCascade["cascade"]
if(not tree_nrCascade):
	print("cascade does not exist in file: %s"%(file_nrCascade))
	sys.exit(0)

n = tree_nrCascade["n"].array(library="np")   # one value per event
cid = tree_nrCascade["cid"].array(library="np")   # array  per event
Elev = tree_nrCascade["Elev"].array(library="np")   # array per event
taus = tree_nrCascade["taus"].array(library="np")   # array value per event
E = tree_nrCascade["E"].array(library="np")   # array value per event
delE = tree_nrCascade["delE"].array(library="np")   # array value per event
I = tree_nrCascade["I"].array(library="np")   # array value per event
Ei = tree_nrCascade["Ei"].array(library="np")   # array value per event
time = tree_nrCascade["time"].array(library="np")   # array value per event
Eg = tree_nrCascade["Eg"].array(library="np") # vector type per event

tree_k100_step3 = file_k100_step3["simtree"]
if(not tree_k100_step3):
	print("simtree does not exist in file: %s"%(file_k100_step3))
	sys.exit(0)

EV_step3 = tree_k100_step3['EV'].array(library="np")
DT_step3 = tree_k100_step3['DT'].array(library="np")
TS_step3 = tree_k100_step3['TS'].array(library="np")
P_step3 = tree_k100_step3['P'].array(library="np")
Type_step3 = tree_k100_step3['Type'].array(library="np")
E1_step3 = tree_k100_step3['E1'].array(library="np")
D3_step3 = tree_k100_step3['D3'].array(library="np")
X1_step3 = tree_k100_step3['X1'].array(library="np") 
Y1_step3 = tree_k100_step3['Y1'].array(library="np") 
Z1_step3 = tree_k100_step3['Z1'].array(library="np") 
X3_step3 = tree_k100_step3['X3'].array(library="np") 
Y3_step3 = tree_k100_step3['Y3'].array(library="np") 
Z3_step3 = tree_k100_step3['Z3'].array(library="np") 
PX1_step3 = tree_k100_step3['PX1'].array(library="np") 
PY1_step3 = tree_k100_step3['PY1'].array(library="np") 
PZ1_step3 = tree_k100_step3['PZ1'].array(library="np") 
PX3_step3 = tree_k100_step3['PX3'].array(library="np") 
PY3_step3 = tree_k100_step3['PY3'].array(library="np") 
PZ3_step3 = tree_k100_step3['PZ3'].array(library="np") 
time1_step3 = tree_k100_step3['time1'].array(library="np")
time3_step3 = tree_k100_step3['time3'].array(library="np")
nCap_step3 = tree_k100_step3['nCap'].array(library="np")

print ("] All files and trees loaded successfully.")
print ("] Now combining files...")

counter = 0
step3_evnt_counter = 0


nCapFlag_bigtree = array('h')
eventid_bigtree = array('L')
n_bigtree = array('H')
cid_bigtree = array('L')
Elev_builder = awkward.ArrayBuilder()
taus_builder = awkward.ArrayBuilder()
E_builder = awkward.ArrayBuilder()
delE_builder = awkward.ArrayBuilder()
I_builder = awkward.ArrayBuilder()
Ei_builder = awkward.ArrayBuilder()
time_builder = awkward.ArrayBuilder()
Eg_builder = awkward.ArrayBuilder()


for i in tqdm(range(len(EV_step1))):
	eventid_bigtree.append(i)
	ncap_index = helper.getNcapIndex(nCap_step1[i])
	if(ncap_index >= 0 and DT_step1[i][ncap_index] == 1 and counter < len(n)): # Condition to check for neutron capture in Si in step1 file
		if(counter+1 != EV_step3[step3_evnt_counter]): # To comensate for event count where photons escaped sensitive material step3 geant4 simulation 
			# do nothing
			# save zeros
			print ("counter+1 : EV_step3[step3_evnt_counter] :: %d : %d"%(counter+1,EV_step3[step3_evnt_counter]))
			sys.exit(0)
			empty_array = array('d')
			nCapFlag_bigtree.append(0)
			n_bigtree.append(0)
			cid_bigtree.append(0)
			Elev_builder.append(empty_array)
			taus_builder.append(empty_array)
			E_builder.append(empty_array)
			delE_builder.append(empty_array)
			I_builder.append(empty_array)
			Ei_builder.append(empty_array)
			time_builder.append(empty_array)
			Eg_builder.append(empty_array)
			counter += 1
		else:
			########################################################################
			### do sanity check with EV_step3 counter, i.e. step3_evnt_counter  ####
			########################################################################

			for j,TS in enumerate(TS_step3[step3_evnt_counter]):
				if(int(P_step3[step3_evnt_counter][j]) != 0):
					continue

				track = int(TS/1.e5)
				step = int(TS) - track*100000
				if(step != 1):
					continue
				if((X1_step3[step3_evnt_counter][j] != X3_step1[i][ncap_index]) or (Y1_step3[step3_evnt_counter][j] != Y3_step1[i][ncap_index]) or (Z1_step3[step3_evnt_counter][j] != Z3_step1[i][ncap_index])):
					print("ERROR : Sanity check failed. Neutron capture at step1 and photon production at step3 are not same.")
					print("-------See below for details---------")
					print("step1 entry = ",i)
					print("step3 entry = ",step3_evnt_counter)
					print("step1 neutron capture location = ({0},{1},{2})".format(X3_step1[i][ncap_index], Y3_step1[i][ncap_index], Z3_step1[i][ncap_index]))
					print("step3 gamma production location = ({0},{1},{2})".format(X1_step3[step3_evnt_counter][j],Y1_step3[step3_evnt_counter][j],Z1_step3[step3_evnt_counter][j]))
					print("-------------------------------------")
					sys.exit(0)
			
			######################## sanity check passed ################################
			#############################################################################


			# do whatever
			#print ( " step 1 event : nrCascade counter : step3_evnt_counter :: {0} : {1} : {2}".format(i,counter,step3_evnt_counter))
			# print ("From step 1:")
			# print (X1_step1[i])
			# print ("From step 3:")
			# print (X1_step3[step3_evnt_counter])

			DT_step1[i] = DT_step3[step3_evnt_counter]
			TS_step1[i] = TS_step3[step3_evnt_counter]
			P_step1[i] = P_step3[step3_evnt_counter]
			Type_step1[i] = Type_step3[step3_evnt_counter]
			E1_step1[i] = E1_step3[step3_evnt_counter]
			D3_step1[i] = D3_step3[step3_evnt_counter]
			X1_step1[i] = X1_step3[step3_evnt_counter]
			Y1_step1[i] = Y1_step3[step3_evnt_counter]
			Z1_step1[i] = Z1_step3[step3_evnt_counter]
			X3_step1[i] = X3_step3[step3_evnt_counter]
			Y3_step1[i] = Y3_step3[step3_evnt_counter]
			Z3_step1[i] = Z3_step3[step3_evnt_counter]
			PX1_step1[i] = PX1_step3[step3_evnt_counter]
			PY1_step1[i] = PY1_step3[step3_evnt_counter]
			PZ1_step1[i] = PZ1_step3[step3_evnt_counter]
			PX3_step1[i] = PX3_step3[step3_evnt_counter]
			PY3_step1[i] = PY3_step3[step3_evnt_counter]
			PZ3_step1[i] = PZ3_step3[step3_evnt_counter]
			time1_step1[i] = time1_step3[step3_evnt_counter]
			time3_step1[i] = time3_step3[step3_evnt_counter]
			nCapFlag_bigtree.append(1)
			nrcsim_evnt = int(EV_step3[step3_evnt_counter]) - 1
			n_bigtree.append(n[nrcsim_evnt])
			cid_bigtree.append(cid[nrcsim_evnt])
			Elev_builder.append(Elev[nrcsim_evnt])
			taus_builder.append(taus[nrcsim_evnt])
			E_builder.append(E[nrcsim_evnt])
			delE_builder.append(delE[nrcsim_evnt])
			I_builder.append(I[nrcsim_evnt])
			Ei_builder.append(Ei[nrcsim_evnt])
			time_builder.append(time[nrcsim_evnt])
			Eg_builder.append(Eg[nrcsim_evnt])
			
			# X1_step1[i] = X1_step3[step3_evnt_counter]
			# print ("updated step 1 values:")
			# print (X1_step1[i])
			# sys.exit(0)
			# .
			# .
			step3_evnt_counter += 1
			counter += 1
	else:
		empty_array = array('d')
		nCapFlag_bigtree.append(0)
		n_bigtree.append(0)
		cid_bigtree.append(0)
		Elev_builder.append(empty_array)
		taus_builder.append(empty_array)
		E_builder.append(empty_array)
		delE_builder.append(empty_array)
		I_builder.append(empty_array)
		Ei_builder.append(empty_array)
		time_builder.append(empty_array)
		Eg_builder.append(empty_array)
		#sanity check




Elev_bigtree = Elev_builder.snapshot()
taus_bigtree = taus_builder.snapshot()
E_bigtree = E_builder.snapshot()
delE_bigtree = delE_builder.snapshot()
I_bigtree = I_builder.snapshot()
Ei_bigtree = Ei_builder.snapshot()
time_bigtree = time_builder.snapshot()
Eg_bigtree = Eg_builder.snapshot()


# print ((eventid_bigtree))
# sys.exit(0)
# print ("length of step1 file = ",len(EV_step1))
# print ("lenth of DT_step1 = ",len(DT_step1))
# print ("length of n_bigtree = ",len(n_bigtree))
# print ("length of Eg_bigtree = ",len(Eg_bigtree))		
# print ("counter = ",counter)
# print ("step3_evnt_counter = ",step3_evnt_counter)
# sys.exit(0)
# event = 153
# print( "event = ",event)
# print ("capture index = ",getNcapIndex(nCap_step1[event]))
# print ("nCapFlag_bigtree[event] = ",nCapFlag_bigtree[event])
# print ("n_bigtree[event] = ",n_bigtree[event])
# print ("Eg_bigtree[event] = ",Eg_bigtree[event])


print ("] Files combined.")
print ("] Writing output file...")

outFileName = inFile_k100_step1.strip(".root")+"_nrCascade_step4.root"
if(os.path.exists(outFileName)):
	os.system("rm -f %s"%(outFileName))

outFile = uproot.recreate(outFileName)
outFile["simtree"] = { 
						"eventid" : eventid_bigtree,
						"nCap_flag" : nCapFlag_bigtree,
						"DT" : DT_step1,
						"TS" : TS_step1,
						"P" : P_step1,
						"Type" : Type_step1,
						"E1" : E1_step1,
						"D3" : D3_step1,
						"X1" : X1_step1,
						"Y1" : Y1_step1,
						"Z1" : Z1_step1,
						"X3" : X3_step1,
						"Y3" : Y3_step1,
						"Z3" : Z3_step1,
						"PX1" : PX1_step1,
						"PY1" : PY1_step1,
						"PZ1" : PZ1_step1,
						"PX3" : PX3_step1,
						"PY3" : PY3_step1,
						"PZ3" : PZ3_step1,
						"time1" : time1_step1,
						"time3" : time3_step1,
						"n" : n_bigtree,
 						"cid": cid_bigtree,
 						"Elev": Elev_bigtree,
 						"taus": taus_bigtree,
 						"E": E_bigtree,
 						"delE": delE_bigtree,
 						"I": I_bigtree,
 						"Ei": Ei_bigtree,
 						"time": time_bigtree,
 						"Eg": Eg_bigtree

 					   }


print ("] Step 4 done.")
print ("] output file: %s"%(outFileName))
		# print ("From step 1:")
		# print (X1_step1[i])
		# print ("From step 3:")
		# print (X1_step3[counter])
		
		# for j,TS in enumerate(TS_step3[counter]):
		# 	track = int(TS/1.e5)
		# 	step = int(TS) - track*100000
		# 	print ("TS : track : step : P : x : y :z :: {0} : {1} : {2} : {3} : {4} : {5} : {6}".format(TS,track,step,int(P_step3[counter][j]),X1_step3[counter][j],Y1_step3[counter][j],Z1_step3[counter][j]) )
		# break
		
