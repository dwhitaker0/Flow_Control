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


import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from PySide import QtCore


experiment_name = input("Enter Experiment Name: ")
experiment_time = float(input("How long should the experiment run? (Minutes) "))
experiment_time = float(experiment_time * 60)
pump_flow1 =  float(input("Pump1 Flow Rate (ml/min): "))
pump_flow2 =  float(input("Pump2 Flow Rate (ml/min): "))
spectral_int_time =  float(input("Spectral Integration Time (Seconds): "))
start_time = float(time.time())
today = datetime.date.today()  
todaystr = today.isoformat() 
experiment_name = "./data/" + todaystr + "/" + experiment_name
spectral_int_time = float(spectral_int_time * 1000000)


#create directory

#if not os.path.exists(experiment_name):
 #   os.makedirs(experiment_name)

if os.path.exists(experiment_name):
	print("Experiment already exists . . . . exiting")
	time.sleep (5)
	sys.exit()
else:
	os.makedirs(experiment_name)


#############	
##Set up logging##
#############

log_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%d/%m/%Y %H:%M:%S')

#File to log to
logFile = os.path.join(experiment_name) + "/logfile.txt"

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

Pump1 = PHDPump.sPump("COM6")
Pump2 = PHDPump.sPump("COM8")
#Dilution_pump = connect_here
logger.info("Connected to Pumps")

#Spectrometer
specdevs = sb.list_devices()
spec = sb.Spectrometer(specdevs[0])
logger.info("Connected to: " + str(specdevs[0]))


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
	
	stored_exeption = None

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
	Pump1.set_flow_rate(pump_flow1)
	Pump2.set_flow_rate(pump_flow2)
	Pump1.start()
	Pump2.start()
	
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
	
	
	
	
	while (float(time.time()) < start_time+float(experiment_time)):
		try:	
		
			spectrum = spec.intensities()
			abs_spectrum = np.log10((reference - dark_ref)/(spectrum - dark_ref))
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
			
			p1.setData(wl,np.log10((reference - dark_ref)/(spectrum - dark_ref)))
			p2.setData(spec_time, ratio)
			pg.QtGui.QApplication.processEvents()
			
			time.sleep((spectral_int_time/1000000)+0.05)
			
			if stored_exeption:
				break
			
		except KeyboardInterrupt:
			logger.info("Experiment Interrupted at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
			stored_exeption = sys.exec_info()
		
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
	
###############
## RUNNING CODE ##
###############
input("Close light shutter and push enter . . . ")
time.sleep(5)
start_dark_reference()
input("Open light shutter and push enter . . . ")
input("Start dilution pump and push enter when ready to take reference . . .")
start_reference()
logger.info ("Reference Acquired at: " + time.strftime("%H:%M:%S | %d-%m-%Y"))
input("Push Enter to Start Experiment. . . ")
start_experiment()
shutdown_and_save()
Pump1.disconnect()
Pump2.disconnect()
exit()


	
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()





