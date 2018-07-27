import Python_ChemFuncts as CF
import numpy as np
import matplotlib.pyplot as plt


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



wl = np.loadtxt("./wl.csv", delimiter = ",")
spectrum = np.loadtxt("./spectrum.csv", delimiter = ",")

smooth = CF.whitsm(spectrum, lmda = 100) #Smooth Spectrum

derivative = CF.savitzky_golay(smooth, window_size = 15, order = 1, deriv=1, rate=1) #SG derivative of smoothed spectrum


plt.subplot(3,1,1)
plt.plot(wl, spectrum, "-")

plt.subplot(3,1,2)
plt.plot(wl, smooth, "-")

plt.subplot(3,1,3)
plt.plot(wl, derivative, "-")

plt.show()



ratio = (derivative[wl == find_nearest(wl,260)] + derivative[wl == find_nearest(wl,290)] ) / derivative[wl == find_nearest(wl,310)]




