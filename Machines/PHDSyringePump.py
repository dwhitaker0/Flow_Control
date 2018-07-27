import os
import sys
import traceback
import time
import serial



# User defined specifications:

#SyringeDiameter = 20.0 #specified in mm (up to 3 d.p)

#FlowRate = 8 #specified in ml/hr (up to 5 d.p)

#ComPort = 4

#DEBUG = True  #Debug flag



#Function to write to serial port

def send_cmd(string,serial_port):
	cmd_send = string + chr(13) + chr(10)
	serial_port.write(cmd_send.encode('ascii'))
	time.sleep(0.3)
	

 

class sPump:	
	def __init__(self, port): #usage pump = FusionSyringePump.sPump('COM8')
		self.port = serial.Serial(port = port,
		baudrate = 115200,
		parity	 = serial.PARITY_NONE, 
		stopbits = serial.STOPBITS_TWO,
		bytesize = serial.EIGHTBITS,
		timeout	 = 0.3)
		#xonxoff = False,#disable software flow 
		#controlrtscts = False) #disable hardware flow control

	
	def connect(self): 
		#Open the serial port
		self.port.open()
		print('Pump Connected')

	def disconnect(self): 
		#Close the serial port
		self.port.close() 
		print('Serial port closed')
		   

	def status(self):
		send_cmd('status',self.port)
		data_raw   = self.port.readline()
		data_pump  = data_raw.decode('ascii')
		print('Response to status: ' +str(data_pump)) #0:stopped, 1:running, 2:paused, 3:delayed, 4:stalled	   

	#Configure the pump: units as ml/hr (1), transfer mode as infusion, syringe diameter as xx mm----------------
	def configure(self, SyringeDiameter):
		send_cmd('hexw2 1 0 '+ str(SyringeDiameter) + ' 20',self.port)
		

	#Read the pump operating limits ------------------------------------------------------------
	
	def limits(self):
		send_cmd('read limit parameter',self.port)
		data_raw   = self.port.readline()  #pump should respond with: max rate, min rate, max volume, min volume
		data_pump  = data_raw.decode('ascii')
		print('Response to read limit parameter: ' +str(data_pump))

	
	#Set the pump flow rate ------------------------------------------------------------
	def set_flow_rate(self, FlowRate):
		send_cmd('irat '+str(FlowRate) + " ml/m",self.port)
		

	#Read the current pump parameters ------------------------------------------------------------
	def read_params(self):
		send_cmd('view parameter',self.port)
		data_raw   = self.port.readline()  #pump should respond with: rate units, syringe inner diameter, transfer rate, priming rate, time, transfer volume, time delay
		data_pump  = data_raw.decode('ascii')
		print('Response to view parameter: ' +str(data_pump))


	#Start the pump ------------------------------------------------------------
	def start(self):
		send_cmd('irun',self.port)
		

   
	#Pause the pump ------------------------------------------------------------
	def pause(self):
		send_cmd('pause',self.port)
		

	
	#Stop the pump ------------------------------------------------------------
	def stop(self):
		send_cmd('stop',self.port)	
		
