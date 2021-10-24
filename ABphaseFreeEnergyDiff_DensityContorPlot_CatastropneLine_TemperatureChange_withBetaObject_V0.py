# This script is used for calculating the difference between
# equilibruium free energy of A and B phase.
# the pressure region is 21 - 34 bar.

# and also it use the eigne values of curvature matrix (coresponding to the gap/mass of higgs modes)
# to find the catastrophne line, on which curvatures vanished.
# The eigen values of curvature matrix came from Mathematica script.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version.3 of ebergy difference code, in which the pico-Joule(pJ 10^-12) unit is used
# J.m^-3 = 10^3 pJ.mm^-3

# And this script also uses the Module_SC_CorrectionObject_V01 to calculate the SC coefficients.

# author: Quang. Zhang (github@hyvatimo)
# version: 0.0

# zeta3 = 1.2020569;


import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module
import matplotlib.pyplot as plot1
import numpy as np

from math import pi
import math

# main code starts from here

m = 1;s = 1; J = 1; Kelvin = 1 # Length unit, Time unit, Energy unit
bar = 1
zeta3 = 1.2020569;
# kb = 8.617333262145*(10**(-5)) #Boltzmann ev.K^-1
kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1
c = 2.99792458*(10**(8)) # speed of light, m.s^-1

# hbar = 6.582119569*(10**(-16)) #plank constant, eV.s
hbar = 1.054571817*(10**(-34)) # planck constant, J.s
# u = 9.3149410242*(10**(8))*eV*(c**(-2)) # atomic mass unit, Dalton, eV.c^-2
u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom


###########################################
# build object, and check the intepolatinon 

BetaObject = SC.BETA('betaAndTc')
    
###################################################
# calculate free energy under different temperature     

stepT = 0.1*(10**-3) 
Temperature = np.arange(0.0*(10**-3), 2.5*(10**-3)+stepT, stepT) #Kelvin

stepPressure = 0.01*bar
pressure = np.arange(20.0, 34.0*bar+stepPressure, stepPressure)

print('Temperature is', Temperature, '\n length of Temperature is ', len(Temperature))
lengthT = len(Temperature)

print('Pressure is',pressure,'\n length of Delta is ', len(pressure))
lengthPressure = len(pressure)

fBGL_array = np.zeros((lengthPressure,lengthT))
fAGL_array = np.zeros((lengthPressure,lengthT))
DiffFABGL = np.zeros((lengthPressure,lengthT)) # save data of the energy difference

DiffFABGLScaled = np.zeros((lengthPressure,lengthT)) # save data of the energy differene, which be scaled by |fBGL|

for iP in range(0, lengthPressure, 1):
    print('\n\n Now P is:', pressure[iP], '\n\n')
   
    indexP = iP
    print('indexP is ',indexP)

    p = pressure[iP]
    
    BetaObject.c1_function(SC.P,SC.c1,p);c1p = BetaObject.c1p
    BetaObject.c2_function(SC.P,SC.c2,p);c2p = BetaObject.c2p
    BetaObject.c3_function(SC.P,SC.c3,p);c3p = BetaObject.c3p
    BetaObject.c4_function(SC.P,SC.c4,p);c4p = BetaObject.c4p
    BetaObject.c5_function(SC.P,SC.c5,p);c5p = BetaObject.c5p
    BetaObject.tc_function(SC.P,SC.Tc,p);Tcp = (BetaObject.tcp)*(10**(-3))
    
    BetaObject.mStar_function(SC.P,SC.Ms,p); mEffective = (BetaObject.ms)*m3
    BetaObject.vFermi_function(SC.P,SC.VF,p); vFermi = BetaObject.vf
    

    c245p = c2p + c4p + c5p;c12p = c1p + c2p;c345p = c3p + c4p + c5p

#    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,'\n\n')

    N0 = ((mEffective**(2))*vFermi)/((2*pi*pi)*(hbar**(3))) # energy density of Fermi surface

    print('\npressure is, ',p,' effective mass is, ', mEffective, ' Fermi velocity is,', vFermi, ' N(0) is ',N0,'\n\n')
    
    for iT in range(0, lengthT, 1):
        #indexDelta = math.floor(delta/stepDelta)
        indexT = iT
#        print('indexT is',indexT)
       
        t = Temperature[indexT]/Tcp
#        print('temperatureis:, ',t)

        if t > 1:

           print(" bro, we just got temperature at Tc, save a np.nan. ")

           DiffFABGL[indexP,indexT] = np.nan
           DiffFABGLScaled[indexP,indexT] = np.nan

           fBGL_array[indexP,indexT] = fBGLRed
           fAGL_array[indexP,indexT] = fAGLRed

           
        else:    

           alphaRed = (1/3)*(t-1)

           beta245Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(2+t*c245p) # A Phase
           beta12345Red = ((7*zeta3)/(240*pi*pi*kb*kb*Tcp*Tcp))*(3.0*(1+t*c12p) + (2+t*c345p)) # B phase

           deltaA = math.sqrt((-alphaRed)/(4*beta245Red)) # A -> B

           deltaB = math.sqrt((-alphaRed)/(2*beta12345Red)) # A -> B
        
           fAGLRed = alphaRed*2*(deltaA**2) + 4*beta245Red*(deltaA**4)
 
#       fBGLRed = alphaRed*3*(delta**2) + 3*beta12345Red*(delta**4)
           fBGLRed = alphaRed*3*(deltaB**2) + 3*beta12345Red*(deltaB**4)


           
           DiffFABGLScaled[indexP,indexT] = (fAGLRed - fBGLRed)/abs(fBGLRed)
           DiffFABGL[indexP,indexT] = DiffFABGLScaled[indexP,indexT]*abs(fBGLRed*N0)
        
           fBGL_array[indexP,indexT] = fBGLRed
           fAGL_array[indexP,indexT] = fAGLRed


##################################################################################################
##################################################################################################
# the following codes are for catastrophe lines

beta4beta5_array = np.zeros((lengthPressure,lengthT))
beta3beta4beta5_array = np.zeros((lengthPressure,lengthT))
beta5_array = np.zeros((lengthPressure,lengthT))
beta1beta3_array = np.zeros((lengthPressure,lengthT))
beta2beta4beta5_array = np.zeros((lengthPressure,lengthT))

for iP in range(0, lengthPressure, 1):
    print('\n\n Now P is:', pressure[iP], '\n\n')
   
    indexP = iP
    print('indexP is ',indexP)

    p = pressure[iP]
    
    BetaObject.c1_function(SC.P,SC.c1,p);c1p = BetaObject.c1p
    BetaObject.c2_function(SC.P,SC.c2,p);c2p = BetaObject.c2p
    BetaObject.c3_function(SC.P,SC.c3,p);c3p = BetaObject.c3p
    BetaObject.c4_function(SC.P,SC.c4,p);c4p = BetaObject.c4p
    BetaObject.c5_function(SC.P,SC.c5,p);c5p = BetaObject.c5p
    BetaObject.tc_function(SC.P,SC.Tc,p);Tcp = (BetaObject.tcp)*(10**(-3))
    
    BetaObject.mStar_function(SC.P,SC.Ms,p); mEffective = (BetaObject.ms)*m3
    BetaObject.vFermi_function(SC.P,SC.VF,p); vFermi = BetaObject.vf
    

    c245p = c2p + c4p + c5p;c12p = c1p + c2p;c345p = c3p + c4p + c5p

#    print('\npressure is ',p,' ,c1p is ',c1p,' ,c2p is ',c2p,' ,c3p is ',c3p,' ,c4p is ',c4p,' ,c4p ',' ,c5p ',c5p,' ,tcp is ',Tcp,'\n\n')

    N0 = ((mEffective**(2))*vFermi)/((2*pi*pi)*(hbar**(3))) # energy density of Fermi surface

    print('\npressure is, ',p,' effective mass is, ', mEffective, ' Fermi velocity is,', vFermi, ' N(0) is ',N0,'\n\n')
    
    for iT in range(0, lengthT, 1):
        #indexDelta = math.floor(delta/stepDelta)
        indexT = iT
#        print('indexT is',indexT)
       
        t = Temperature[indexT]/Tcp
#        print('temperatureis:, ',t)

        if t > 1:

           print(" bro, we just got temperature at Tc, do nothing ")

                     
        else:    

           deltaA = math.sqrt((-alphaRed)/(4*beta245Red)) # A -> B
            
           betatilde = (7*zeta3)/(240*pi*pi*kb*kb)

           beta1WC = (-1/(Tcp*Tcp))*betatilde
           beta1SC = (c1p/(Tcp*Tcp))*betatilde
           beta1 = beta1WC +beta1SC

           beta2WC = (2/(Tcp*Tcp))*betatilde
           beta2SC = (c2p/(Tcp*Tcp))*betatilde
           beta2 = beta2WC + beta2SC
           
           beta3WC = (2/(Tcp*Tcp))*betatilde
           beta3SC = (c3p/(Tcp*Tcp))*betatilde
           beta3 = beta3WC + beta3SC

           beta4WC = (2/(Tcp*Tcp))*betatilde
           beta4SC = (c4p/(Tcp*Tcp))*betatilde
           beta4 = beta4WC + t*beta4SC

           beta5WC = (-2/(Tcp*Tcp))*betatilde
           beta5SC = (c5p/(Tcp*Tcp))*betatilde
           beta5 = beta5WC + t*beta5SC

           beta4beta5_array[indexP,indexT] = -4*(beta4 + beta5)
           beta3beta4beta5_array[indexP,indexT] = 4*(beta3-beta4-beta5)
           beta5_array[indexP,indexT] = -8*beta5
           beta1beta3_array[indexP,indexT] = 8*(beta1 + beta3)
           beta2beta4beta5_array[indexP,indexT] = 8*(beta2+beta4+beta5)


    # print('fBGLRed_Delta is:', fBGL_array[indexP,:])
    # print('fAGLRed_Delta is:', fAGL_array[indexP,:])
    # print('difference of fGLRed_Delta is:', DiffFABGL[indexP,:])
    # print('difference between A & B fGLRed_Delta is:', DiffFABGLScaled[indexP,:])
#    plo
#    plot1.plot(Delta, fGL_array[indexT,:], 'o-'); plot1.ylabel('ev^2'); plot1.xlabel('ev')
#    plot1.show()
    

# # Plot Free energy difference
# for iP in range(0, lengthPressure, 1):
#     #indexT = math.floor(T/stepT)
#     indexP = iP
#     plot1.plot(Temperature*1000, DiffFABGL[indexP,:], '-.'); plot1.ylabel(r'$(f{A}_{GL}-f^{B}_{GL}).N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')
    
# plot1.show()    
# plot1.clf()
# plot1.cla()
# plot1.close()    

# # Plot Free Energies for A and B 
# for iP in range(0, lengthPressure, 1):
#     #indexT = math.floor(T/stepT)
#     indexP = iP
   
#     plot1.plot(Temperature*1000, fAGL_array[indexP,:], '--'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')
#     plot1.plot(Temperature*1000, fBGL_array[indexP,:], '-'); plot1.ylabel(r'$f_{GL}.N(0)^{-1}/ev^{2}$'); plot1.xlabel(r'$T$/mK')

# plot1.show()    
 
# plot1.clf()
# plot1.cla()
# plot1.close()    

# density and contour plot of the fAGL - fBGL
# DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, DiffFABGL*(10**(3)));plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)
#Cplot = plot1.contour(Temperature*1000, pressure, DiffFABGL*(10**(3)));plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
#CplotCatas = plot1.contour(Temperature*1000, pressure, beta5_array);plot1.clabel(CplotCatas, inline=True, fontsize=8.5, colors='k')
#plot1.savefig('DensityPlot_FreeEnergyDiff_SI_unit_moduleV01.pdf');

DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, beta2beta4beta5_array);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)
plot1.show()

plot1.clf()
plot1.cla()
plot1.close()    

# density and contour plot of (fAGL - fBGL)/|fBGL|
#DensityPlot = plot1.pcolormesh(Temperature*1000, pressure, DiffFABGLScaled);plot1.ylabel(r'$p/bar$'); plot1.xlabel(r'$T$/mK');plot1.colorbar(DensityPlot)
#Cplot = plot1.contour(Temperature*1000, pressure, DiffFABGLScaled);plot1.clabel(Cplot, inline=True, fontsize=8.5, colors='r')
# plot1.savefig('DensityPlot_FreeEnergyDiff_Scaled_SI_unit_moduleV01.pdf');
#plot1.show()

#plot1.clf()
#plot1.cla()
#plot1.close()


# for iP in [70, 80, 90, 100, 110, 120]:
#     print(" pressure is: ",pressure[iP])
#     #indexT = math.floor(T/stepT)
#     indexP = iP
   
#     plot1.plot(Temperature*1000, DiffFABGLScaled[indexP,:], '-.'); plot1.ylabel(r'$(f{A}_{GL}-f^{B}_{GL})$/$|f^{B}|$'); plot1.xlabel(r'$T$/mK')

# plot1.savefig('FreeEnergyDiff_Scaled.pdf');
# plot1.show()

