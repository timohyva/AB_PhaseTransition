
'''
contour plot of p-T GL coherent length

author: Quang, (timohyva@github)
'''


import Module_SC_Beta_V05 as SCC
import he3_tools_Vn01 as h


import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

import Module_New_Parpia_pT_parameter as Parpia
import Module_Parpia_pT_ConsQ as ConsQ

SCC.turn_on_PLTS()

Pressure = np.arange(18.0,30.1,0.01) # bar
# Pressure = np.arange(20.0,24.01,0.01) # bar

Temperature = (np.arange(1.8, 2.4, 0.005))*(10**(-3)) # Kelvin
# Temperature = (np.arange(2.0, 2.4, 0.005))*(10**(-3)) # Kelvin

T_arr, P_arr = np.meshgrid(Temperature, Pressure)
print("\n T_arr looks like \n", T_arr,"\n P_arr looks like \n", P_arr )

# xiGL_RWS_arr = np.zeros(T_arr.shape)

xiGL_RWS_arr = SCC.xiGL_JWS(P_arr, T_arr)*(10**6)

print(" \n xiGL_RWS_arr looks like ",xiGL_RWS_arr)


###########################################################################
####                         mask TAB tails                          ######
###########################################################################

TAB_PLTS_masked = np.array([])

for p in Pressure:

    if p >= h.p_pcp_bar:
       TAB_PLTS_masked = np.append(TAB_PLTS_masked, (10**3)*SCC.TAB_RWSco(p))

    elif p < h.p_pcp_bar:
       TAB_PLTS_masked = np.append(TAB_PLTS_masked, np.nan)        

# TAB_arr = SCC.TAB_RWSco(Pressure)*(10**3)

Tc_arr = SCC.Tcp(Pressure)*(10**3)

    
###########################################################################
#####                         plot contours                          ######
###########################################################################
    
# Levels1 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3]
# Levels1 = np.arange(0.005, 1.02, 0.005) # micro meter
Levels1 = [0.309, 0.62, 1.1]

mSize = 100 # scatter marker size
            
fig1, ax = plt.subplots(1,1)

# cf1 = ax.contourf(T_arr*(10**3), P_arr, 10*xiGL_RWS_arr, cmap=cm.PuBu_r, levels=Levels1)

# fig1.colorbar(cf1, ax=ax, location = 'left')

c1 = ax.contour(T_arr*(10**3), P_arr, 10*xiGL_RWS_arr, levels=Levels1, colors='magenta')
plt.clabel(c1, inline=True, fontsize=18, colors='r')

# # plot Greywall T_AB line in axx[0]
# TABGreywall1, = ax.plot(h.TAB_poly_Greywall(Pressure-h.p_pcp_bar), Pressure, color = "orange")
# TcGreywall1, = ax.plot(h.Tc_poly_Greywall(Pressure), Pressure, color = "red")

# plot PLTS T_AB lien
# TABPlts1, = ax.plot(h.TAB_poly_PLTS(Pressure), Pressure, color = "orange")
TABPlts1, = ax.plot(TAB_PLTS_masked, Pressure, color = "orange")

# plot PLTS Tc
# TCPlts, = ax.plot(h.Tc_mK(Pressure, "PLTS"), Pressure, color = "purple")
TCPlts1, = ax.plot(Tc_arr, Pressure, color = "red")


# scatter plot of parpia's constant pressure data
sIC1 = ax.scatter(Parpia.T_IC, Parpia.pressure, color = "yellow", marker = "o", s = mSize, label=" Const-P, IC")
sHEC1 = ax.scatter(Parpia.T_HEC, Parpia.pressure, color = "red", marker = "x", s = mSize, label="Const-P, HEC")

# scatter plot of parpia's constant Q data

sHEC_CQ1 = ax.scatter(ConsQ.T_HEC, ConsQ.p_HEC, color = "purple", marker = "s", s = mSize, label="Cosnt-Q, HEC")
sIC_CQ1 = ax.scatter(ConsQ.T_IC, ConsQ.p_IC, color = "cyan", marker = "1", s = mSize, label="Const-Q, IC")


# leg1 = ax.legend([TABGreywall1,TcGreywall1, sIC1, sHEC1, sIC_CQ1, sHEC_CQ1],[r"$T_{AB}^{Greywall}$",r"$T_{c}^{Greywall}$",r"$p-T-IC$",r"$p-T-HEC$",r"$p-T-IC-CQ$",r"$p-T-HEC-CQ$"],  fontsize=15.0, loc='lower right')
leg1 = ax.legend([TABPlts1,TCPlts1, sIC1, sHEC1, sIC_CQ1, sHEC_CQ1],[r"$T_{AB}^{PLTS}$",r"$T_{c}^{PLTS}$",r"$Const-p, IC$",r"$Const-p, HEC$",r"$Const-Q, IC$",r"$Const-Q, HEC$"],  fontsize=15.0, loc='lower right')

ax.set_ylabel(r'$p/bar$', fontsize=18.0);
ax.set_xlabel(r'$T$/mK', fontsize=18.0);
# ax.set_xlim([2.10, 2.3])
# ax.set_ylim([18.0, 22.0])
ax.set_title(r"$\xi_{GL}(p, T)/{\mu}m$", fontsize=16.0)
ax.grid(True)

text_kwargs1 = dict(ha='center', va='center', fontsize=28, color='C1')
text_kwargs2 = dict(ha='center', va='center', fontsize=28, color='blue')
plt.text(1.95, 28., 'thin wall region', **text_kwargs1)
plt.text(2.32, 24., 'exteme thick wall region', **text_kwargs2)

plt.show()
