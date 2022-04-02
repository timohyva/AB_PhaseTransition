


import he3_tools as h
import numpy as np
import matplotlib.pyplot as plt

Pressure = np.arange(21.1,30.1,0.1) # bar
TAB_RWS = h.T_mK(h.t_AB(Pressure), Pressure)

print("\n the first 15 elemnts of TAB_RWS looks like,", TAB_RWS[:16])

fig, ax = plt.subplots(1,1)

Tc_G, = ax.plot(h.Tc_mK(Pressure), Pressure, color="purple")

TAB_greywall_poly, = ax.plot(h.T_AB_Greywall_poly(Pressure), Pressure, color="orange")

TAB_greywall_poly_abs, = ax.plot(h.T_AB_Greywall_poly(Pressure, "abs"), Pressure, color="orange")

ax.legend([Tc_G, TAB_greywall_poly, TAB_greywall_poly_abs],[r"$T_{c}^{G}$",r"$T_{AB}^{G-poly}$",r"$T_{AB}^{G-poy-abs}$"])

ax.set_xlabel(r"$T/mK$")
ax.set_ylabel(r"$p/bar$")

plt.show()
