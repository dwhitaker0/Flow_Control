import numpy as np
import pandas as pd
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

from Machines import PHDSyringePump as PHDPump
import time
import seabreeze.spectrometers as sb
import datetime
import logging
import os
import sys
import numpy as np
import serial
import Python_ChemFuncts as CF

###################################
############USER INPUT###############
###################################
experiment_name = input("Enter Experiment Name: ")
experiment_time = float(input("How long should the experiment run? (Minutes) "))
experiment_time = float(experiment_time * 60)
pump_flow1 =   0 #float(input("Pump1 Flow Rate (ml/min): "))
pump_flow2 =   0 #float(input("Pump2 Flow Rate (ml/min): "))
spectral_int_time =  float(input("Spectral Integration Time (Seconds): "))
start_time = float(time.time())
today = datetime.date.today()  
todaystr = today.isoformat() 
#experiment_name = "./data/" + todaystr + "/" + experiment_name
root_experiment_name = "./data/" + todaystr + "/" + str(experiment_name)
spectral_int_time = float(spectral_int_time * 1000000)
experiment_number = 1

initial_x = input("Enter Initial X: ") #0.5
initial_x = float(initial_x)
initial_y = input("Enter Initial Y: ") #0.5
initial_y = float(initial_y)
initial_step_x = input("Initial Step X: ") #0.1
initial_step_x = float(initial_step_x)
initial_step_y = input("Initial Step Y: ") #0.1
initial_step_y = float(initial_step_y)

col_names =['Flow 1', 'Flow 2', 'Ratio']


###Algorithm Variables #####

max_iterations = 10
n = 2 #Number of variables
no_improv_thresh = 1
no_improv_break = 3

a = 1		#NM algorithm factor
b = 2		#NM algorithm factor
g = 0.5		#NM algorithm factor
k = 1.5 	#NM algorithm factor
d = 0.5		#NM algorithm factor

## Constraints

x_max = 1
x_min = 0.05

y_max = 1
y_min = 0.05





if os.path.exists(root_experiment_name):
		print("Experiment already exists . . . . exiting")
		time.sleep (5)
		sys.exit()
else:
	os.makedirs(root_experiment_name)

###################################
###################################
###################################

#############	
##Set up logging##
#############

log_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%d/%m/%Y %H:%M:%S')

#File to log to
logFile = os.path.join(root_experiment_name) + "/logfile.txt"

#Setup File handler
file_handler = logging.FileHandler(logFile)
file_handler.setFormatter(log_formatter)
file_handler.setLevel(logging.INFO)

#Setup Stream Handler (console)
#stream_handler = logging.StreamHandler()
#stream_handler.setFormatter(log_formatter)
#stream_handler.setLevel(logging.INFO)

#Get logger
logger = logging.getLogger('root')
logger.setLevel(logging.INFO)

#Add both Handlers
logger.addHandler(file_handler)
#logger.addHandler(stream_handler)

#################
##Connect to machines##
#################

#Pump

Pump1 = PHDPump.sPump("COM6") ##Change COM to where pump is connected
Pump2 = PHDPump.sPump("COM8")
#Dilution_pump = connect_here
logger.info("Connected to Pumps")

#Spectrometer
specdevs = sb.list_devices()
spec = sb.Spectrometer(specdevs[0])
logger.info("Connected to: " + str(specdevs[0]))

#####################################
#############Set up plotting ##############
#####################################

#QtGui.QApplication.setGraphicsSystem('raster')
app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle("Flow Experimental System")
#mw.resize(600,400)
mw.setGeometry(6,150, 650,450 )
cw = QtGui.QWidget()
mw.setCentralWidget(cw)
l = QtGui.QVBoxLayout()
cw.setLayout(l)


pw = pg.PlotWidget(name="UV Spectrum") 
l.addWidget(pw)
pw2 = pg.PlotWidget(name="Lambda Max Vs Time")
l.addWidget(pw2)

mw.show()
## Create an empty plot curve to be filled later, set its pen
p1 = pw.plot(pen=2, title="UV Spectrum" )
pw.enableAutoRange(enable=True)
pw.setYRange(-2,2,padding=0)
pw.setLabel('left', 'Absorbance', units='Arbitr. Units')
pw.setLabel('bottom', 'Wavelength', units='nm')
p2 = pw2.plot(pen=2, symbol = 'o', symbolPen = 2, title="Ratio Vs Time")
pw2.enableAutoRange(enable=True)
pw2.setLabel('left', 'Ratio', units='Arbitr. Units')
pw2.setLabel('bottom', 'Time', units='s')

######################################
######################################
######################################



#To initialise simplex a starting co-ordinate is required (Xin), n+1 points are generated with step_size (initial_step_x & y) away from Xin

#####################################

def start_dark_reference():
	global dark_ref
	spec.integration_time_micros(spectral_int_time)
	dark_ref = spec.intensities()


def start_reference():
	global reference
	#Dilution_pump.start()
	#time.sleep(10)
	spec.integration_time_micros(spectral_int_time)
	reference = spec.intensities()
	#time.sleep(1)
	#Dilution_pump.stop()
	
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



def start_experiment():
	global wl
	global raw_spectral_data
	global refd_spectral_data
	global start_time
	global spec_time
	global lambda_max
	global smooth
	#global SGderiv
	global peak_ratio
	global iteration_result
	global experiment_number
	
	stored_exeption = None
	
	experiment_name = root_experiment_name + "/" + str(experiment_number)
	
	if os.path.exists(experiment_name):
		print("Experiment already exists . . . . exiting")
		time.sleep (5)
		sys.exit()
	else:
		os.makedirs(experiment_name)

	wl = spec.wavelengths()
	raw_spectral_data = spec.wavelengths()
	Abs_spectral_data = spec.wavelengths()
	#SGderiv = spec.wavelengths()
	spec_time = np.array(0) 
	start_time = time.time()
	lambda_max = np.array(0)
	ratio = np.array(0)
	logger.info ("Experiment started at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
	
	#set-up and start pump
	logger.info("Setting Pump 1 @ " + str(pump_flow1))
	Pump1.set_flow_rate(pump_flow1)
	logger.info("Setting Pump 2 @ " +str( pump_flow2))
	Pump2.set_flow_rate(pump_flow2)
	Pump1.start()
	Pump2.start()

	
	while (float(time.time()) < start_time+float(experiment_time)):
		try:	
		
			spectrum = spec.intensities()
			abs_spectrum = np.log10((reference - dark_ref)/(spectrum - dark_ref))
			Where_Nan = np.isnan(abs_spectrum)
			Where_Inf = np.isinf(abs_spectrum)
			abs_spectrum[Where_Nan] = 0
			abs_spectrum[Where_Inf] = 0
			smooth = CF.whitsm(abs_spectrum, lmda = 100)
			deriv_spectrum = CF.savitzky_golay(smooth, window_size = 15, order = 1, deriv=1, rate=1)
			sp_time = time.time() 
			raw_spectral_data = np.vstack((raw_spectral_data, spectrum))
			Abs_spectral_data = np.vstack((Abs_spectral_data, abs_spectrum))
			spec_time =  np.append(spec_time, [sp_time-start_time])
			curr_lambda_max = np.amax(spectrum - reference)
			lambda_max = np.append(lambda_max, [curr_lambda_max])			
			#SGderiv = np.vstack(SGderiv, deriv_spectrum)
			ratio = np.append(ratio, (deriv_spectrum[wl == find_nearest(wl,260)] + deriv_spectrum[wl == find_nearest(wl,290)] ) / deriv_spectrum[wl == find_nearest(wl,310)])
			
			
			#Update Plot
			
			#p1.setData(wl,np.log10((reference - dark_ref)/(spectrum - dark_ref)))
			p1.setData(wl,smooth)
			p2.setData(spec_time, ratio)
			pg.QtGui.QApplication.processEvents()
			
			time.sleep((spectral_int_time/1000000)+0.05)
			
			if stored_exeption:
				break
			
		except KeyboardInterrupt:
			logger.info("Experiment Interrupted at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
			stored_exeption = sys.exc_info()
		
	logger.info("Experiment complete at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
	Pump1.stop()
	Pump2.stop()
	final_data = raw_spectral_data
	final_Abs = Abs_spectral_data
	Where_Nan = np.isnan(Abs_spectral_data)
	Where_Inf = np.isinf(Abs_spectral_data)
	final_Abs[Where_Nan] = 0
	final_Abs[Where_Inf] = 0
	np.savetxt(os.path.join(experiment_name) + "/spectral_results.csv", final_data , fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/Abs_spectral_results.csv", final_Abs, fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/reference.csv", reference, fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/dark_reference.csv", dark_ref, fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/wl.csv", wl, fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/times.csv", spec_time, fmt="%s", delimiter=",")
	np.savetxt(os.path.join(experiment_name) + "/ratio.csv", ratio, fmt="%s", delimiter=",")
	logger.info("Spectral data saved as: " + os.path.join(experiment_name) + "/spectral_results.csv")
	iteration_result = np.median(ratio[-len(ratio * 0.2): ])
	experiment_number = experiment_number + 1

#############
####Begin#####
#############

iterations =  0
Results	 = pd.DataFrame()

input("Close light shutter and push enter . . . ")
time.sleep(1)
start_dark_reference()
input("Open light shutter and push enter . . . ")
input("Start dilution pump and push enter when ready to take reference . . .")
start_reference()
logger.info ("Reference Acquired at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
input("Push Enter to Start Experiment. . . ")

pump_flow1 = initial_x 
pump_flow2 = initial_y 

start_experiment()
Results = pd.DataFrame([(pump_flow1,pump_flow2,iteration_result)], columns=col_names)

pump_flow1 = initial_x + initial_step_x
pump_flow2 = initial_y

start_experiment()
Results = Results.append(pd.DataFrame([(pump_flow1, pump_flow2, iteration_result)], columns=col_names), ignore_index=True)

pump_flow1 = initial_x 
pump_flow2 = initial_y + initial_step_y

start_experiment()
Results = Results.append(pd.DataFrame([(pump_flow1, pump_flow2, iteration_result)], columns=col_names), ignore_index=True)
print(Results)

#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

##############

###########OPTIMISATION ALGORITHM###############

no_improv = 0 


while iterations < max_iterations:
	
	sorted = Results.sort_values(by='Ratio', ascending = False)
	yield_h = sorted.iloc[n,n]
	yield_s = sorted.iloc[n-1,n]
	yield_l = sorted.iloc[n-2,n]
	 
	if yield_l - yield_s >= no_improv_thresh:
		no_improv = 0
	else:
		no_improv = no_improv + 1
		
	if no_improv == no_improv_break:
		print ("No Further Improvement, exiting......")
		logger.info("No Further Improvement, experiment stopped" + time.strftime("%H:%M:%S | %d-%m-%Y"))
		Pump1.disconnect()
		Pump2.disconnect()
		exit()
		print (Results)
		np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
		break
	 
	 
	h = sorted.iloc[n,:n]
	s = sorted.iloc[n-1,:n]
	l = sorted.iloc[n-2,:n]

	best_vertex = [(s.iloc[0],l.iloc[0]),(s.iloc[1],l.iloc[1])]
	centroid = (sum(best_vertex[0])/len(best_vertex[0]),sum(best_vertex[1])/len(best_vertex[1]))



####Reflection####

	Xr = centroid + (centroid - h)
	print ("Next Point = " + "\n" + str(Xr))
	
	Xr.iloc[0] = np.clip(Xr.iloc[0], x_min, x_max)
	Xr.iloc[1] = np.clip(Xr.iloc[1], y_min, y_max)
	
	pump_flow1 = Xr.iloc[0]
	pump_flow2 = Xr.iloc[1]
	
	start_experiment()

	Res_Xr =  iteration_result 

	Results = Results.append(pd.DataFrame([(pump_flow1, pump_flow2, iteration_result)], columns=col_names), ignore_index=True)
	print(Results)
	np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
	#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))  




####NO MOVE####

	if Res_Xr < yield_l and Res_Xr >= yield_s:
		iterations = iterations +1
		print ("No Move")
		logger.info("No Move")
		continue



####Expand####

	if Res_Xr > yield_l:
		Xe = centroid + b*(centroid - h)
		print ("Next Point = " + "\n" + str(Xe))
		
		Xe.iloc[0] = np.clip(Xe.iloc[0], x_min, x_max)
		Xe.iloc[1] = np.clip(Xe.iloc[1], y_min, y_max)
		
		pump_flow1 = Xe.iloc[0]
		pump_flow2 = Xe.iloc[1]
	
		start_experiment()
		#Res_Xe = float(raw_input("Xe: "))
		Res_Xe = iteration_result

#Final condition#

		if Res_Xe >= yield_l:
			Results = Results.append(pd.DataFrame([(Xe.iloc[0], Xe.iloc[1], Res_Xe)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

			iterations = iterations + 1
			continue

		else: #doing nothing will just keep R as the new point
			Results = Results.append(pd.DataFrame([(Xe.iloc[0], Xe.iloc[1], Res_Xe)], columns=col_names), ignore_index=True)
			iterations = iterations + 1
			continue

	  

####Inside Contraction####

	if Res_Xr <= yield_h:

		#Xic = g*centroid + g*h
		Xic = centroid - g*(Xr - centroid)
		print ("Next Point = " + "\n" + str(Xic))
		
		Xic.iloc[0] = np.clip(Xic.iloc[0], x_min, x_max)
		Xic.iloc[1] = np.clip(Xic.iloc[1], y_min, y_max)
		
		pump_flow1 = Xic.iloc[0]
		pump_flow2 = Xic.iloc[1]
	
		start_experiment()

		Res_Xic = iteration_result

	#Final Condition#

		if Res_Xic >= yield_h:

			Results = Results.append(pd.DataFrame([(Xic.iloc[0], Xic.iloc[1], Res_Xic)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

			iterations = iterations + 1
			continue

		else: ##Shrink##

			Xredh = l + d*(h - l)
			print ("Next Set Points = " + "\n" + str(Xredh))

			Xredh.iloc[0] = np.clip(Xredh.iloc[0], x_min, x_max)
			Xredh.iloc[1] = np.clip(Xredh.iloc[1], y_min, y_max)
			
			pump_flow1 = Xredh.iloc[0]
			pump_flow2 = Xredh.iloc[1]
	
			start_experiment()
			Res_Xredh = iteration_result
			
			Results = Results.append(pd.DataFrame([(Xic.iloc[0], Xic.iloc[1], Res_Xic)], columns=col_names), ignore_index=True)
			Results = Results.append(pd.DataFrame([(Xredh.iloc[0], Xredh.iloc[1], Res_Xredh)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))


			Xreds = l + d*(s - l)
			print ("Next Set Points = " + "\n" + str(Xreds))
			
			Xreds.iloc[0] = np.clip(Xreds.iloc[0], x_min, x_max)
			Xreds.iloc[1] = np.clip(Xreds.iloc[1], y_min, y_max)

			pump_flow1 = Xreds.iloc[0]
			pump_flow2 = Xreds.iloc[1]
	
			start_experiment()
			Res_Xreds = iteration_result


			Results = Results.append(pd.DataFrame([(Xreds.iloc[0], Xreds.iloc[1], Res_Xreds)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

			
			iterations = iterations + 1
			continue

						

####Outside Contraction####

	if Res_Xr > yield_h and Res_Xr <= yield_s:

		#Xoc = k*centroid - g*h
		Xoc = centroid + g*(Xr - centroid)
		print ("Next Point = " + "\n" + str(Xoc))
		
		Xoc.iloc[0] = np.clip(Xoc.iloc[0], x_min, x_max)
		Xoc.iloc[1] = np.clip(Xoc.iloc[1], y_min, y_max)

		pump_flow1 = Xoc.iloc[0]
		pump_flow2 = Xoc.iloc[1]
	
		start_experiment()
		Res_Xoc = iteration_result
	  


#Final Condition#

		if Res_Xoc > Res_Xr:

			Results = Results.append(pd.DataFrame([(Xoc.iloc[0], Xoc.iloc[1], Res_Xoc)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

			iterations = iterations + 1
			continue

		else: ##Shrink##

			Xredh = l + d*(h - l)
			print ("Next Set Points = "  + "\n" + str(Xredh))
			
			Xredh.iloc[0] = np.clip(Xredh.iloc[0], x_min, x_max)
			Xredh.iloc[1] = np.clip(Xredh.iloc[1], y_min, y_max)

			pump_flow1 = Xredh.iloc[0]
			pump_flow2 = Xredh.iloc[1]
	
			start_experiment()
			Res_Xredh = iteration_result

			Results = Results.append(pd.DataFrame([(Xoc.iloc[0], Xoc.iloc[1], Res_Xoc)], columns=col_names), ignore_index=True)
			Results = Results.append(pd.DataFrame([(Xredh.iloc[0], Xredh.iloc[1], Res_Xredh)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))


			Xreds = l + d*(s - l)
			print ("Next Set Points = "+ "\n" + str(Xreds))
			
			Xreds.iloc[0] = np.clip(Xreds.iloc[0], x_min, x_max)
			Xreds.iloc[1] = np.clip(Xreds.iloc[1], y_min, y_max)
		
			pump_flow1 = Xreds.iloc[0]
			pump_flow2 = Xreds.iloc[1]
	
			start_experiment()
			Res_Xreds = iteration_result

		
			Results = Results.append(pd.DataFrame([(Xreds.iloc[0], Xreds.iloc[1], Res_Xreds)], columns=col_names), ignore_index=True)
			print(Results)
			np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
			#p1.setData(y=np.array(Results.iloc[:,0]), x=np.array(Results.iloc[:,1]))

					
			iterations = iterations + 1
			continue
			
np.savetxt(os.path.join(root_experiment_name) + "/results.csv", Results , fmt="%s", delimiter=",")
Pump1.disconnect()
Pump1.disconnect()
exit()

