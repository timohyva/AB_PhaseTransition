############################################################
####        Important notations and descriptions       #####
############################################################


# This script is used for calculating and plotting the ratio f(\phi_{-})/f_{\phi_{+}} of free energy
# deduced from original GL free energy by setting \phi \times D_{\alpha i} = A_{B} - A_{A}. Here A_{A} and A_{B}
# are particular order parameters of A and B phases. The explicit expression see the Mathematica note book.
# The resulted contour plot and density plot show \bar{\lambda} > 0.93 in whole p-T region with H = 0.

# This script uses SC_Beta module to caltulate the M^{2}, \delta and \lambda parameter in the p-T plane,
# then plugging them into the expressions of f(\phi_{-}) and f(\phi_{+}) to get the ratio f(\phi_{-})/f_{\phi_{+}}

# A new version of this code with Lotynk't experimentla parameters will be developed soon based on this code

#log time: 10th. January. 2022


###########################################################################################################
##############         Significant Developments & Processes of WP2 project         ########################
###########################################################################################################

# Two new modules are created based on the primary SC_CorrectionObject module.
# One is SC_Beta module and other one is plot_TAB_line module.

# The basic idea is using A. J. Leggtt's simply considerations when f_A - f_B is small and
# R_c / \xi, R_c/t are both very tiny, and close to T_AB
# For details, checking his Journal of Low Temperature vol.87 571 1992 and Yip & Leggtt review
# in Helium Three 1990.

# the critical radius are plotted as countours, and compared with the experiment data from Lotynk's experiment.
# the result is amazing, and suggests R_c ~ 400 \xi

# Mark suggested another evaluation about R_c in 19th Nov 2021, see the email in details
# the basic idea is R_c ~ ((f_A \xi_GL)/ \Deltaf_AB)/\xi_GL = f_A/ \Deltaf_AB
# temperature region from 0.0 mk to 2.4 mK

# the pressure region is 21 - 34 bar.

# This script uses the SC_Beta module and plot_TAB module, both of them use SC_CorrectionObject module.

# strong coupling correction coefficients of \beta comes from PRB. 101. 024517
# the free energy difference is rescaled by absolute value of equilibrium free energy of B phase

# This is version.2 of ebergy difference code, in which the pico-Joule(pJ 10^-12) unit is used 

# author: Quang. Zhang (github@hyvatimo)

# zeta3 = 1.2020569;

############################################################################
######                       Modules & Constants                     #######
############################################################################

import Module_SC_Beta_V05 as SCC
import Module_AB_wall_V00 as WB

# Tc, TAB data digitized from Parpia's manuscript
import Module_p_Tc_TAB_Parpia as pTcTAB

# import Module_plot_TAB_line as TAB_line
import Module_Lotynk_pT_parameter as Lotynk
import Module_New_Parpia_pT_parameter as Parpia

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from math import pi
import math

# main code starts from here

# Length unit, Time unit, Energy unit, mass unit, pressure unit 
zeta3 = 1.2020569;
m = 1;s = 1; J = 1; Kelvin = 1; kg =1; bar = 1

kb = 1.380649*(10**(-23)) # Boltzmann constant J.K^-1

hbar = 1.054571817*(10**(-34)) # planck constant, J.s

u = 1.66053906660*(10**(-27)) # atomic mass unit, Dalton, kg

m3 = 3.016293*u #mass of helium3 atom
              

#############################################################################
####                 calculate ratio for Lotynk parameter               #####
#############################################################################

# import the Lotynk experimental p-T parameters
# Temperature_Lotynk = np.asarray(Lotynk.Temperature) # mili Kelvin
# pressure_Lotynk = np.asarray(Lotynk.pressure)

# # saving array 
# fB_Lotynk = np.array([])
# fBA_Lotynk = np.array([])
# xiGL_Lotynk = np.array([])
# tensionAB_Lotynk = np.array([])

##############################################################################
#####                    calculate the Rc array                          #####
##############################################################################

Pressure = pTcTAB.p

Tc_pressure = (pTcTAB.Tc)*(10**(-3)) # Kelvin

temperature = np.arange(0.5,1.1,0.01) # scaled Temperature t=T/Tc(p)

t_array, P_array = np.meshgrid(temperature, Pressure)

T_array = np.zeros((len(Pressure), len(temperature))); print("\n\n",T_array)
# T_array = np.zeros((np.size(t_array))); print("\n\n",T_array)

# print("\n P_array looks like\n", P_array,"\n t_array looks like\n", t_array, " \n ", t_array[3,:], "\n\n",P_array[:,12])

print("\n P_array looks like\n", P_array,"\n \n t_array looks like \n", t_array)

for ii in range(0, len(Pressure), 1):
  T_array[ii,:] = t_array[ii, :]*Tc_pressure[ii]

print(" \n T_array looks like\n ", T_array)

##############################################################################
#####                   calcuate sigmaAB for pTcTAB.p                    #####
##############################################################################

tensionAB_p = np.array([])
for p in Pressure:
  tensionAB_p = np.append(tensionAB_p, WB.get_wall_n_surfaceEnergy(p))

print("\n tensionAB_p looks like ", tensionAB_p)   

##############################################################################


# # saving array 
# fB = np.array([])
# fBA = np.array([])
# xiGL = np.array([])

Rc = np.zeros((len(Pressure),len(temperature)))


# #############################################################################

for iP in range(0, len(Pressure), 1):
    print('\n\n Now P is:', Pressure[iP], '\n\n')

    p = Pressure[iP]

    Tcp = Tc_pressure[iP]

    sigma_p = tensionAB_p[iP]

    # T = Temperature_Lotynk[iP]*(10**(-3))*Kelvin
    # t = Temperature[indexT]/Tcp
    # t = T/SCC.Tcp(p)

    for it in range(0, len(temperature), 1):

       t = temperature[it]

       print('\n scaled temperature is:, ',t)
    
       if t >= 1:

        print(" bro, we just got temperature at Tc, save a np.nan. ")

        # mask data from t > 1 region
        # fB = np.append(fB, np.nan)
    
        # fBA = np.append(fBA, np.nan)

        # xiGL = np.append(xiGL, np.nan)

        Rc[iP, it] = np.nan

                 
       else:

        # fB = np.append(fB, SCC.alpha_td(p,T)/SCC.betaB_td(p,T))
     
        # fBA = np.append(fBA, (-SCC.alpha_td(p,T)/SCC.betaB_td(p,T)) + (SCC.alpha_td(p,T)/SCC.betaA_td(p,T)))

        # xiGL = np.append(xiGL, SCC.xiGL_OC(p, T))

        # tensionAB = np.append(tensionAB, WB.get_wall_n_surfaceEnergy(p))

        fB = SCC.alpha_td_Parpia(p,t)/SCC.betaB_td_Parpia(p,t)
     
        fBA = (-SCC.alpha_td_Parpia(p,t)/SCC.betaB_td_Parpia(p,t)) + (SCC.alpha_td_Parpia(p,t)/SCC.betaA_td_Parpia(p,t))
        
        xi_tuple = SCC.xiGL_OC_Parpia(p, Tcp, t)

        xiGL = xi_tuple[0] # t - dependent xiGL

        print(" \n fB is ", fB, " fBA is ", fBA, " xiGL is ", xiGL, " sigmaAB is", sigma_p, " xi_tuple looks like :", xi_tuple)

        Rc[iP, it] = -1.*((sigma_p*fB)/fBA)*xiGL
        print("\n sigmaAB/fBA looks like : ", Rc[iP, it])

# -2sigma/fBA
# print("\n tensionAB looks like",tensionAB)
# Rc = -1.*((tensionAB*fB)/fBA)*xiGL; print("\n sigmaAB/fBA looks like : ", Rc)

print(" \n Rc array looks like ",Rc)


#############################################################################

#############################################################################
####                 calculate ratio for Parpia new parameter           #####
#############################################################################

# saving array 
# fB_IC_Parpia = np.array([])
# fBA_IC_Parpia = np.array([])
# xiGL_IC_Parpia = np.array([])
# tensionAB_IC_Parpia = np.array([])

# fB_HEC_Parpia = np.array([])
# fBA_HEC_Parpia = np.array([])
# xiGL_HEC_Parpia = np.array([])
# tensionAB_HEC_Parpia = np.array([])

#############################################################################

# for iP in range(0, len(Parpia.pressure), 1):
#     print('\n\n Now Parpia P is:', Parpia.pressure[iP], '\n\n')

#     p = Parpia.pressure[iP]

#     TIC = Parpia.T_IC[iP]*(10**(-3))*Kelvin

#     THEC = Parpia.T_HEC[iP]*(10**(-3))*Kelvin

#     # t = Temperature/Tcp
#     t_IC = TIC/SCC.Tcp(p)

#     t_HEC = THEC/SCC.Tcp(p)

#     print("\n T_IC is :", TIC, ' scaled temperature is: ', t_IC, " T_HEC is : ", THEC, " scaled temperature is: ", t_HEC)
    
#     if t_IC >= 1:

#       print(" bro, we just got temperature at Tc, save a np.nan. ")

#       # mask data from t > 1 region
#       fB_IC_Parpia = np.append(fB_IC_Parpia, np.nan)
    
#       fBA_IC_Parpia = np.append(fBA_IC_Parpia, np.nan)

#       xiGL_IC_Parpia = np.append(xiGL_IC_Parpia, np.nan)

#       tensionAB_IC_Parpia = np.append(tensionAB_IC_Parpia, np.nan)

                 
#     else:

#       fB_IC_Parpia = np.append(fB_IC_Parpia, SCC.alpha_td(p,TIC)/SCC.betaB_td(p,TIC))
     
#       fBA_IC_Parpia = np.append(fBA_IC_Parpia, (-SCC.alpha_td(p,TIC)/SCC.betaB_td(p,TIC)) + (SCC.alpha_td(p,TIC)/SCC.betaA_td(p,TIC)))

#       xiGL_IC_Parpia = np.append(xiGL_IC_Parpia, SCC.xiGL_OC(p, TIC))

#       tensionAB_IC_Parpia = np.append( tensionAB_IC_Parpia, WB.get_wall_n_surfaceEnergy(p))

#       print(" \n Parpia fB_IC is ", fB_IC_Parpia[iP], " Parpia fBA_IC is ", fBA_IC_Parpia[iP], " Parpia xiGL is ", xiGL_IC_Parpia[iP], " Parpia tension IC is ",tensionAB_IC_Parpia[iP])

#     ###########################################################################
#     ###########################################################################

#     if t_HEC >= 1:

#       print(" bro, we just got temperature at Tc, save a np.nan. ")

#       # mask data from t > 1 region
#       fB_HEC_Parpia = np.append(fB_HEC_Parpia, np.nan)
    
#       fBA_HEC_Parpia = np.append(fBA_HEC_Parpia, np.nan)

#       xiGL_HEC_Parpia = np.append(xiGL_HEC_Parpia, np.nan)

#       tensionAB_HEC_Parpia = np.append(tensionAB_HEC_Parpia, np.nan)

                 
#     else:

#       fB_HEC_Parpia = np.append(fB_HEC_Parpia, SCC.alpha_td(p,THEC)/SCC.betaB_td(p,THEC))
     
#       fBA_HEC_Parpia = np.append(fBA_HEC_Parpia, (-SCC.alpha_td(p,THEC)/SCC.betaB_td(p,THEC)) + (SCC.alpha_td(p,THEC)/SCC.betaA_td(p,THEC)))

#       xiGL_HEC_Parpia = np.append(xiGL_HEC_Parpia, SCC.xiGL_OC(p, THEC))

#       tensionAB_HEC_Parpia = np.append(tensionAB_HEC_Parpia, WB.get_wall_n_surfaceEnergy(p)) 

#       print(" \n Parpia fB_HEC is ", fB_HEC_Parpia[iP], " Parpia fBA_HEC is ", fBA_HEC_Parpia[iP], " xiGL_HEC is ", xiGL_HEC_Parpia[iP], " Parpia tension HEC is ",tensionAB_HEC_Parpia[iP])  

# # ratio -2fBxiGL/fBA
# print("\n xiGL_OC looks like ", xiGL_IC_Parpia)
# print("\n tensionAB_IC_Parpia looks like",tensionAB_IC_Parpia)
# # Rc_IC_Parpia = -2.*((tensionAB_IC_Parpia*fB_IC_Parpia)/fBA_IC_Parpia)*xiGL_IC_Parpia; print("\n 2sigma/fBA IC looks like : ", Rc_IC_Parpia) # sphere evaluation
# Rc_IC_Parpia = -1.*((tensionAB_IC_Parpia*fB_IC_Parpia)/fBA_IC_Parpia)*xiGL_IC_Parpia; print("\n 2sigma/fBA IC looks like : ", Rc_IC_Parpia) # Celiydill evaluation

# print("\n xiGL_OC looks like ", xiGL_IC_Parpia)
# print("\n tensionAB_IC_Parpia looks like",tensionAB_IC_Parpia)
# # Rc_HEC_Parpia = -2.*((tensionAB_HEC_Parpia*fB_HEC_Parpia)/fBA_HEC_Parpia)*xiGL_HEC_Parpia; print("\n fBxiGL/fBA HEC looks like : ", Rc_HEC_Parpia) # sphere evaluation
# Rc_HEC_Parpia = -1.*((tensionAB_HEC_Parpia*fB_HEC_Parpia)/fBA_HEC_Parpia)*xiGL_HEC_Parpia; print("\n fBxiGL/fBA HEC looks like : ", Rc_HEC_Parpia) # Celydrell evaluation

#############################################################################
#####           scatter plot the ratio verse the pressure            ########
#############################################################################

fig, ax = plt.subplots(1,1)
Levels = [0.5, 1., 1.2, 1.5, 2., 2.5, 3, 4., 6., 8., 10.]

c1 = ax.contourf(T_array*(10**(3)), P_array, Rc*(10**(6)), cmap=cm.PuBu_r, levels=Levels);
ax.set_ylabel(r'$p/bar$');
ax.set_xlabel(r'$T$/mK');

fig.colorbar(c1, ax=ax)

c2 = ax.contour(T_array*(10**(3)), P_array, Rc*(10**(6)), levels=Levels, colors='black');
plt.clabel(c2, inline=True, fontsize=8.5, colors='r')

ax.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow")
ax.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x")

ax.set_title(r"contour plot of $(\sigma_{AB}(p)/f_{AB}(p,T))/{\mu}m$")
ax.set_xlim([1.8, 2.3])
ax.grid(True)

fig.savefig("Rc_contour_RWS_Tc_digitised_from_ParpiaManuscript.pdf")

plt.show()


# fig, ax = plt.subplots(2,1)
# ax[0].scatter(Pressure, temperature, Rc/(10**(-6)), color = "green"); 
# # ax[0].scatter(Parpia.pressure, Rc_IC_Parpia/(10**(-6)), color = "blue");
# # ax[0].scatter(Parpia.pressure, Rc_HEC_Parpia/(10**(-6)), color = "red", marker ="x")

# # ax[0].legend(labels=("PRL 2021 slow cooling events","p-T_IC", "p-T_HEC"), loc="upper right")
# ax[0].legend(labels=("p-T_IC", "p-T_HEC"), loc="upper right")

# ax[0].set_xlabel(r"p/bar");
# ax[0].set_ylabel(r"$R_{c}=(\frac{-\sigma_{AB}(p)}{f_{B}-f_{A}})/{\mu}m$")
# ax[0].grid(True)
# ax[0].set_title(r"$R_{c}$ of half-cylinder-bubble (1st fig) and (p, T) plot of events in HEC, IC")

# ax[1].scatter(Parpia.pressure, Parpia.T_IC, color = "blue")
# ax[1].scatter(Parpia.pressure, Parpia.T_HEC, color = "red", marker = "x")
# ax[1].legend(labels=("p-T_IC", "p-T_HEC"), loc="upper right")
# ax[1].set_xlabel(r"p/bar");
# ax[1].set_ylabel(r"T/mK");
# ax[1].grid(True)
# ax[1].set_title(r" $(p, T)$ plot of Transition events ")

# fig.savefig("CriticalRadius_MarkDate_halfShpere_Parpia_AndEvents.pdf")
# plt.show()

# fig1, ax1 = plt.subplots(1,1)
# ax1.scatter(Parpia.pressure, tensionAB_IC_Parpia, color = "green")
# ax1.set_xlabel(r"p/bar")
# ax1.set_ylabel(r"$\sigma_{AB}(p)/\xi_{GL}^{OC}(T)|f_{B}(p, T)|$")
# ax1.grid(True)

# ax1.set_title(r"domain wall tension $\sigma_{AB}/\xi_{GL}^{OC}(p,T)|f_{B}(p,T)|$")
# fig1.savefig("surface_tension_IC_Parpia.pdf")
# plt.show()               


