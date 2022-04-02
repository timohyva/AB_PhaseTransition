

#########################################################################################
# this script shows the t_AB in he3_tools is same as my implement in Module_plot_TAB_line
#########################################################################################

import he3_tools as h
import Module_plot_TAB_line as TAB_RWS

import matplotlib.pyplot as plt
import numpy as np

Pressure = np.arange(21.5,30.1,0.1)


TAB1 = h.T_mK(h.t_AB(Pressure), Pressure)

fig, ax = plt.subplots(1,1)

TABrws = ax.contour(TAB_RWS.X*(10**3), TAB_RWS.Y, TAB_RWS.EnergyDensity_Difference_fABGL, levels=[0.0], colors='red')

ax.plot(TAB1, Pressure, color="blue", marker="x")

plt.show()
