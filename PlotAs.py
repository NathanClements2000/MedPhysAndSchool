import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


Data = pd.read_csv('ML3-Glider-M09.csv')

plt.plot(Data['TRun1'],Data['M2Run1'])
plt.show()

Times = Data['TRun1'].copy()

#print(Times)

#Fourier = np.cos(2*np.pi*0.1*Times)
#print(Fourier)

Freqs = np.linspace(0.1, 1.5, 300)
FCosTrans1 = np.zeros(300)
FSinTrans1 = np.zeros(300)
FCosTrans2 = np.zeros(300)
FSinTrans2 = np.zeros(300)
FCosTrans3 = np.zeros(300)
FSinTrans3 = np.zeros(300)
#print(Freqs)

for i in range(len(Freqs)):
    FourierC = np.cos(2*np.pi*Freqs[i]*Times)
    ConvoluteC = FourierC*Data['M1Run1']
    FourierS = np.sin(2*np.pi*Freqs[i]*Times)
    ConvoluteS = FourierS*Data['M1Run1']    
    integral = ConvoluteC.sum()*(Times[1]-Times[0])
    FCosTrans1[i] = integral
    integral = ConvoluteS.sum()*(Times[1]-Times[0])
    FSinTrans1[i] = integral
    
plt.plot(Freqs, FCosTrans1)
plt.plot(Freqs, FSinTrans1)
plt.plot(Freqs, np.sqrt(FCosTrans1**2 + FSinTrans1**2))
plt.show()


for i in range(len(Freqs)):
    FourierC = np.cos(2*np.pi*Freqs[i]*Times)
    ConvoluteC = FourierC*Data['M2Run1']
    FourierS = np.sin(2*np.pi*Freqs[i]*Times)
    ConvoluteS = FourierS*Data['M2Run1']    
    integral = ConvoluteC.sum()*(Times[1]-Times[0])
    FCosTrans2[i] = integral
    integral = ConvoluteS.sum()*(Times[1]-Times[0])
    FSinTrans2[i] = integral
    
plt.plot(Freqs, FCosTrans2)
plt.plot(Freqs, FSinTrans2)
plt.plot(Freqs, np.sqrt(FCosTrans2**2 + FSinTrans2**2))
plt.show()

for i in range(len(Freqs)):
    FourierC = np.cos(2*np.pi*Freqs[i]*Times)
    ConvoluteC = FourierC*Data['M3Run1']
    FourierS = np.sin(2*np.pi*Freqs[i]*Times)
    ConvoluteS = FourierS*Data['M3Run1']    
    integral = ConvoluteC.sum()*(Times[1]-Times[0])
    FCosTrans3[i] = integral
    integral = ConvoluteS.sum()*(Times[1]-Times[0])
    FSinTrans3[i] = integral
    
plt.plot(Freqs, FCosTrans3)
plt.plot(Freqs, FSinTrans3)
plt.plot(Freqs, np.sqrt(FCosTrans3**2 + FSinTrans3**2))
plt.show()

plt.plot(Freqs, np.sqrt(FCosTrans1**2 + FSinTrans1**2)+2,label='Mass 1')
plt.plot(Freqs, np.sqrt(FCosTrans2**2 + FSinTrans2**2),label = 'Mass 2')
plt.plot(Freqs, np.sqrt(FCosTrans3**2 + FSinTrans3**2)-2, label = 'Mass 3')
plt.legend(loc = 'upper right')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Offset Acceleration Amplitude ($cm/s^2$)')
plt.title('Fourier Spectrum for run 1 - Absolute values')
plt.savefig('SpectrumRun1.pdf')
plt.show()

plt.plot(Freqs, FCosTrans1+2,label='Mass 1')
plt.plot(Freqs, FCosTrans2,label = 'Mass 2')
plt.plot(Freqs, FCosTrans3-2, label = 'Mass 3')
plt.legend(loc = 'upper right')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Offset Acceleration Amplitude ($cm/s^2$)')
plt.title('Fourier cosine Spectrum for run 1')
plt.savefig('SpectrumRun1Cos.pdf')
plt.show()

plt.plot(Freqs, FCosTrans1,label='Mass 1')
plt.plot(Freqs, FCosTrans2,label = 'Mass 2')
plt.plot(Freqs, FCosTrans3, label = 'Mass 3')
plt.legend(loc = 'upper right')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Offset Acceleration Amplitude ($cm/s^2$)')
plt.title('Fourier cosine Spectrum for run 1 no offset')
plt.savefig('SpectrumRun1CosNO.pdf')
plt.show()


