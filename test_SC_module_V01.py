# This scipt is for testing the SC correction module v.0.1

import matplotlib.pyplot as plot1
import numpy as np
# from math import pi

import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module

# test the object

BetaObject = SC.BETA('betaAndTc') # us class BETA in mudule SC to create BetaObject

# p = 32 # pressure, bar
pressureArray = np.arange(0.0, 34.0+2.0, 2.0)

for p in pressureArray:
    
    BetaObject.c1_function(SC.P,SC.c1,p);c1p = BetaObject.c1p
    BetaObject.c2_function(SC.P,SC.c2,p);c2p = BetaObject.c2p
    BetaObject.c3_function(SC.P,SC.c3,p);c3p = BetaObject.c3p
    BetaObject.c4_function(SC.P,SC.c4,p);c4p = BetaObject.c4p
    BetaObject.c5_function(SC.P,SC.c5,p);c5p = BetaObject.c5p
    BetaObject.tc_function(SC.P,SC.Tc,p);Tcp = BetaObject.tcp

    BetaObject.mStar_function(SC.P,SC.Ms,p);mEffective = BetaObject.ms
    BetaObject.vFermi_function(SC.P,SC.VF,p);vFermi = BetaObject.vf

    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,' ,effective mass is ', mEffective,' ,Fermi velocity is ', vFermi,'\n\n\n')

    print('\npressure is ',p,' ,effective mass is ', mEffective,' ,Fermi velocity is ', vFermi,'\n\n\n')


# test the vectorization, SC_Container is container list for all coefficents arraies,

SC_Container = []
BetaObject.c1_function(SC.P,SC.c1,pressureArray);SC_Container.append(BetaObject.c1p)
BetaObject.c2_function(SC.P,SC.c2,pressureArray);SC_Container.append(BetaObject.c2p)
BetaObject.c3_function(SC.P,SC.c3,pressureArray);SC_Container.append(BetaObject.c3p)
BetaObject.c4_function(SC.P,SC.c4,pressureArray);SC_Container.append(BetaObject.c4p)
BetaObject.c5_function(SC.P,SC.c5,pressureArray);SC_Container.append(BetaObject.c5p)
BetaObject.tc_function(SC.P,SC.Tc,pressureArray);SC_Container.append(BetaObject.tcp)

BetaObject.mStar_function(SC.P,SC.Ms,pressureArray);SC_Container.append(BetaObject.ms)
BetaObject.vFermi_function(SC.P,SC.VF,pressureArray);SC_Container.append(BetaObject.vf)

print(SC_Container)

plot1.plot(pressureArray,SC_Container[0],'-x',pressureArray,SC.c1,'or')
plot1.plot(pressureArray,SC_Container[1],'-x',pressureArray,SC.c2,'or')
plot1.plot(pressureArray,SC_Container[2],'-x',pressureArray,SC.c3,'or')
plot1.plot(pressureArray,SC_Container[3],'-x',pressureArray,SC.c4,'or')
plot1.plot(pressureArray,SC_Container[4],'-x',pressureArray,SC.c5,'or')
plot1.plot(pressureArray,SC_Container[5],'-x',pressureArray,SC.Tc,'or')
plot1.plot(pressureArray,SC_Container[6],'-x',pressureArray,SC.Ms,'or')
plot1.plot(pressureArray,SC_Container[7],'-x',pressureArray,SC.VF,'or')

plot1.show()

