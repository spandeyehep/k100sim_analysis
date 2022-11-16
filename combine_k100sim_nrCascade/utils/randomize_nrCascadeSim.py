import os, sys
import uproot, awkward
#import ROOT
import numpy as np
from array import array
#from mpl_toolkits import mplot3d
import numpy as np
#import matplotlib.pyplot as plt



def randomize(inFile_nrCascade):
	file_nrCascade = uproot.open(inFile_nrCascade)
	if(not file_nrCascade):
		print("could not open file: %s"%(inFile_nrCascade))
		sys.exit(0)

	print ("] File opened succssfully...")
	tree_nrCascade = file_nrCascade["cascade"]
	if(not tree_nrCascade):
		print("simtree does not exist in file: %s"%(file_nrCascade))
		sys.exit(0)

	print ("] Tree loaded succssfully...")

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

	print ("] Randomizing...")

	shuffler = np.random.permutation(len(n)) 
	n = n[shuffler]
	cid = cid[shuffler]
	Elev = Elev[shuffler]
	taus = taus[shuffler]
	E = E[shuffler]
	delE = delE[shuffler]
	I = I[shuffler]
	Ei = Ei[shuffler]
	time = time[shuffler]
	Eg = Eg[shuffler]

	print ("] Randomized succssfully...")

	OnlyInRootFile = inFile_nrCascade.split("/")[-1]
	outFile_nrCascade = inFile_nrCascade.strip(OnlyInRootFile) + OnlyInRootFile.strip(".root") + "_randomized.root"
	#outFile_nrCascade = inFile_nrCascade.split("/")[-1].strip(".root") + "_randomized.root"

	outFile = uproot.recreate(outFile_nrCascade)
	outFile["cascade"] = { 
							"original_Eid": shuffler,
							"n" : n,
	 						"cid": cid,
	 						"Elev": Elev,
	 						"taus": taus,
	 						"E": E,
	 						"delE": delE,
	 						"I": I,
	 						"Ei": Ei,
	 						"time": time,
	 						"Eg": Eg
	 					   }


	print ("] File written succssfully...")
	print ("] Output file: %s"%(outFile_nrCascade))