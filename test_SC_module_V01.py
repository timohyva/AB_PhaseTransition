# This scipt is for testing the SC correction module

# import matplotlib.pyplot as plot1
import numpy as np
# from math import pi


import Module_SC_CorrectionObject as SC # strong coupling correction module

# test the object

BetaObject = SC.BETA('betaAndTc') # us class BETA in mudule SC to create BetaObject

# p = 32 # pressure, bar
pressureArray = np.arange(0.0, 34.0, 2.0)

for p in pressureArray:
    
    judgementList = SC.region_judgement(p) # region_judgement is in SC module

    print('judgementList is', judgementList)

    print('low pressure is: ',judgementList[0], ',high pressure is: ',judgementList[1], ',interpolation region is: ', judgementList[2])

    pk = judgementList[0]; pk1 = judgementList[1]; k = judgementList[2]

    BetaObject.c1_function(SC.c1,p,pk,pk1,k);c1p = BetaObject.c1p
    BetaObject.c2_function(SC.c2,p,pk,pk1,k);c2p = BetaObject.c2p
    BetaObject.c3_function(SC.c3,p,pk,pk1,k);c3p = BetaObject.c3p
    BetaObject.c4_function(SC.c4,p,pk,pk1,k);c4p = BetaObject.c4p
    BetaObject.c5_function(SC.c5,p,pk,pk1,k);c5p = BetaObject.c5p
    BetaObject.tc_function(SC.Tc,p,pk,pk1,k);Tcp = BetaObject.tcp

    BetaObject.mStar_function(SC.Ms,p,pk,pk1,k);mEffective = BetaObject.ms
    BetaObject.vFermi_function(SC.VF,p,pk,pk1,k);vFermi = BetaObject.vf

    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,' ,effective mass is ', mEffective,' ,Fermi velocity is ', vFermi,'\n\n\n')

    print('\npressure is ',p,' ,effective mass is ', mEffective,' ,Fermi velocity is ', vFermi,'\n\n\n')
