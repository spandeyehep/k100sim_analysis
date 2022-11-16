import os, sys
import uproot, awkward
#import ROOT
import numpy as np
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from tqdm import tqdm


inFile = "/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_50M_PuBe_sourceAndshields/sim_50M_PuBe_sourceAndshields.root"
file = uproot.open(inFile)
tree = file["simtree"]
event = tree['EV'].array(library="np")
x_prestep = tree['X1'].array(library="np") 
y_prestep = tree['Y1'].array(library="np") 
z_prestep = tree['Z1'].array(library="np") 
pdgid = tree['Type'].array(library="np") 
Edep = tree['D3'].array(library="np") 
TrackStep = tree['TS'].array(library="np") 
Time1 = tree['time1'].array(library="np") 
# print (x_prestep[0])
# print (Edep[0])
save_ = False

print("x_prestep : Edep :: %d : %d"%(len(x_prestep[321]), len(Edep[321])))


def data_for_cylinder_along_z(center_x,center_y,radius,z_center,height_z):
    z = np.linspace(z_center - height_z/2, z_center + height_z/2, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def data_for_rectangle_along_z(center_x,center_y,center_z, width, length, height):
    x_grid = np.linspace(center_x - width , center_x + width, 100)
    y_grid = np.linspace(center_y - length , center_y + length, 100)
    z_grid = np.linspace(center_z - height , center_z + height, 100)
    z = np.linspace(center_z - height , center_z + height, 100)
    x = np.linspace(center_x - width , center_x + width, 100)
    y = np.linspace(center_y - length , center_y + length, 100)
    x_grid, y_grid = np.meshgrid(x, y)
    for i in z:
    	


    # theta = np.linspace(0, 2*np.pi, 50)
    # theta_grid, z_grid=np.meshgrid(theta, z)
    # x_grid = radius*np.cos(theta_grid) + center_x
    # y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def show_event(event_, time_, x_, y_, z_, pid_, elevation, azimuth):
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	
	for i in range(len(pid_)):
		if(pid_[i] == 11):
			ax.scatter(x_[i], y_[i], z_[i],marker="v",c="r")
		elif(pid_[i] == 22):
			ax.scatter(x_[i], y_[i], z_[i],marker="^",c="b")
		else:
			ax.scatter(x_[i], y_[i], z_[i],marker="o",c="g")
	
	# Draw particle source (Am60)
	# Xc,Yc,Zc = data_for_cylinder_along_z(0.0,0.0,1.143 , 58.055899, 2)
	# ax.plot_surface(Xc, Yc, Zc, color="r")

	# Draw detector crystal
	#Xc,Yc,Zc = data_for_cylinder_along_z(0.0,0.0,50.0,20.0 + 33.3/2,33.3)
	# center = (-405, 0, 100.248) mm
	# size/2 = (61, 213 ,637.5) mm 
	# Xc,Yc,Zc = data_for_rectangle_along_z(-405., 0., 100.2, 61., 213. ,637.5)
	# ax.plot_surface(Xc, Yc, Zc, alpha=0.5) 

	
	ax.set_xlabel('X1 [mm]')
	ax.set_ylabel('Y1 [mm]')
	ax.set_zlabel('Z1 [mm]')
	ax.set_title("Global time = %0.1f ps"%(time_*10000))
	print("Global time ",(time_*1000))
	fig.suptitle("Am60 14 keV gamma, Event : %d"%(event_), fontsize=16)
	#ax.set_xlim(-60.0,60.0)
	#ax.set_ylim(-60.0,60.0)
	#ax.set_zlim(15.0,55.0)
	ax.view_init(elev=elevation, azim=azimuth)
	#savename = "./plots/gamma_100kev_PuBeGammaOff/event%d/event%d_time%d.png"%(event_,event_,int(time_*10000))
	savename = "./plots/event%d_time%d.png"%(event_,int(time_*10000))
	if(save_):
		plt.savefig(savename)
	print("plot saved : %s"%(savename))
	#sys.exit(0)
	#plt.show()
	#input("raw")

for i in tqdm(range(1)):
	i = 321
	#print ("event = %d"%(i))
	if( len(x_prestep[i]) != len(Edep[i]) ):
		print("Incompatible length for event %d"%(i))
		break

	Track = (TrackStep[i]/100000)
	Track = Track.astype('int32')

	# create 2D array from two 1D array
	combined = np.vstack((Time1[i],Track)).T

	# now add columns to 2D array
	combined = np.c_[combined,pdgid[i]]
	combined = np.c_[combined,x_prestep[i]]
	combined = np.c_[combined,y_prestep[i]]
	combined = np.c_[combined,z_prestep[i]]
	# 0 : time
	# 1 : track
	# 2 : pid
	# 3 : x
	# 4 : y
	# 5 : z


	#print(i, " : ",combined)

	# now sort 2D array accroding to time (0th column)
	combined = combined[combined[:,0].argsort()]
	print("After sorting")
	#print(i, " : ",combined)
	#print("length of combined: ", len(combined))
	x_ = []
	y_ = []
	z_ = []
	pid_ = []
	
	time_= 0.0
	unique_time = len(np.unique(Time1[i]))
	print("unique_time = ",unique_time)
	elv = np.linspace(30,1,unique_time)
	azim = np.linspace(60,1,unique_time)
	counter = 0
	for j in tqdm(range(len(combined))):
		pid_.append(combined[j][2])
		x_.append(combined[j][3])
		y_.append(combined[j][4])
		z_.append(combined[j][5])

		if(time_ != combined[j][0]):
			# print ("**** time_ = ",time_," *****")
			# print (x_)
			# print (y_)
			# print (z_)
			show_event(i, combined[j][0], x_, y_, z_, pid_,elv[counter],azim[counter])
			# x_ = []
			# y_ = []
			# z_ = []
			# pid_ = []
			counter += 1
			time_ = combined[j][0]
			#unique_time.append(time_)

		# else:
		# 	pid_.append(combined[j][2])
		# 	x_.append(combined[j][3])
		# 	y_.append(combined[j][4])
		# 	z_.append(combined[j][5])
	# print("unique_time = ",unique_time)
	# print("count = ",len(unique_time))



print("Done!!!")