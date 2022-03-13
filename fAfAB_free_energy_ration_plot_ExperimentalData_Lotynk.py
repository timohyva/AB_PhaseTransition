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


# import Module_SC_CorrectionObject_V01 as SC # strong coupling correction module
import Module_SC_Beta_V02 as SCC 
# import Module_plot_TAB_line as TAB_line
import Module_Lotynk_pT_parameter as Lotynk

import matplotlib.pyplot as plt
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
Temperature_Lotynk = np.asarray(Lotynk.Temperature) # mili Kelvin
pressure_Lotynk = np.asarray(Lotynk.pressure)

# saving array 
fA_Lotynk = np.array([])
fBA_Lotynk = np.array([])

#############################################################################

for iP in range(0, len(pressure_Lotynk), 1):
    print('\n\n Now Lotynk P is:', pressure_Lotynk[iP], '\n\n')

    p = pressure_Lotynk[iP]

    T = Temperature_Lotynk[iP]*(10**(-3))*Kelvin

    # t = Temperature[indexT]/Tcp
    t = T/SCC.Tcp(p)   

    print("\n Temperature is :", T, ' scaled temperature is:, ',t)
    
    if t >= 1:

      print(" bro, we just got temperature at Tc, save a np.nan. ")

      # mask data from t > 1 region
      fA_Lotynk = np.append(fA_Lotynk, np.nan)
    
      fBA_Lotynk = np.append(fBA_lotynk, np.nan)

                 
    else:

      fA_Lotynk = np.append(fA_Lotynk, SCC.alpha_td(p,T)/SCC.betaA_td(p,T))
     
      fBA_Lotynk = np.append(fBA_Lotynk, (SCC.alpha_td(p,T)/SCC.betaB_td(p,T)) - (SCC.alpha_td(p,T)/SCC.betaA_td(p,T)))

      print(" \n Lotynk fA is ", fA_Lotynk[iP], " Lotynk fBA is ", fBA_Lotynk[iP])

# ratio fA/fBA
ratio_Lotynk = fA_Lotynk/fBA_Lotynk; print("\n fA/fBA looks like : ", ratio_Lotynk)

#############################################################################
#####           scatter plot the ratio verse the pressure            ########
#############################################################################

fig, ax = plt.subplots(1,1)
ax.scatter(pressure_Lotynk, ratio_Lotynk, color = "green"); ax.set_xlabel(r"p/bar"); ax.set_ylabel(r"$\frac{f_{A}}{f_{B}-f_{A}}$")
ax.grid(True)
plt.show()



