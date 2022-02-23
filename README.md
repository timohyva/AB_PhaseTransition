# AB_PhaseTransition wiki page

Repsositry of modules and codes devoted on A-B transtion, bubble nucleation and some confinements physics 

![WP2][QUEST_DMC_WP2]

## Descriptions
This repository provides modules and codes about many calculations based on strong coupling (SC) corrected Ginzburg-Landau model of Superfluid Helium-3. It is developed to implement different kinds of calculations involving SC-Ginzburg-Landau free energy and dydic tensor field/order parameter of symmetry-breking phases of Helium-3 in mK level temperature. The strong coupling correction data sheet comes from Physical Review B papers by Robert.C.Regan, J.J.Wiman and James. A. Sauls:
1. *"Vortex phase diagram of rotating superfluid  3He−B"*
    [Vortex phase diagram of rotating superfluid 3He−B](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.024517)
2. *"Superfluid phases of 3He in nanoscale channels"* 
    [Superfluid phases of 3He in nanoscale channels](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.144515)


The basic block of this repo consists of three modules:
<p align="center">
    Module-SCCO-V0*, Module-SC-Beta-V0*, Module-SCC-bB-V0*
</p>

[SCCO](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SCCO_V02.py) offers an object, which directly talks with strong coupling data sheet and generates required SC corrected data. The calling requirments is set as coming from [SC-Beta](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SC_Beta_V03.py). This module provides the interfaces, though which the pre-indices, the dimensionless munerical values of SC data and the characteristic dimensional qualities (e.g., Tc, N(0), Kb * Tc, V_{Fermi} and coherent length) can be independently called as object attributes. It allows most of interested qualities in Helium-3 can be called though the handy module [SCC-bB](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SCC_bB_V00.py). For the detailed informations about how these features are implemented, please read the head comments of module files and How To Use)


(* The implement of physics introduced by confiment geometry will be considered/developed late *)

## Languages (current)
* python3


## License
[GPL] (https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.en.html)

## How To Use

### Python ![Python](https://github.com/timohyva/AB_PhaseTransition/blob/master/logo_languge11.png)
The defult way is launching from bash. And it's ok for runing in IDE.
```shell
# import SCC-bB module. 
~$ python3 import Module_SCC_bB_V0* as SCCB 
```
After importing, the phyical qualities can be invoked as object attributes:
```python
SCCB.alpha_bar(p, T)
SCCB.beta1_bar(p, T)
SCCB.beta2_bar(p, T)
```
or
```python
SCCB.betaA_bar(p, T)
SCCB.betaB_bar(p, T)
```

and
```python
SCCB.DeltaA2_bar(p, T)
SCCB.DeltaB2_bar(p, T)
```
...

## Module fphi-AtoB

Module fphi-AtoB-V0* is an example about using SCCB objects. This module provides barrier location (\phi_{-}), barrier height (f(\phi_{-})), free energy (potential line) etc along A and B phases:
```python
import Module_fphi_AtoB_V01 as AB
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1)
phi = np.arange(0., 6.5, 0.1, dtype=np.float64)
for p in range(0, 34, 2): ax.plot(phi, AB.fphi_bar(p, 0.5*(10**(-3)), phi))
ax.set_xlabel(r"$\phi/k_{b} \dot T_{c}$");ax.set_ylabel(r"$f(\phi)/(N(0) \times (kb \dot Tc)^{2})$");ax.grid(True)

plt.show()
```

##
<p align="center">
<img width="450" src="https://github.com/timohyva/AB_PhaseTransition/blob/master/Contour_And_Density_Plot_Of_1st_EigenvalueOfCurvatureMatirx.png" alt="curvature">
</p>

<p align="center">
<img width="450" src="https://github.com/timohyva/AB_PhaseTransition/blob/master/Lotynk_Data_Plot_fafab.png" alt="Lotynk">
</p>

<p align="center">
<img width="450" src="https://github.com/timohyva/AB_PhaseTransition/blob/master/fphiAB.png" alt="fphiAB">
</p>


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

[plot1]: https://github.com/timohyva/AB_PhaseTransition/blob/master/Contour_And_Density_Plot_Of_1st_EigenvalueOfCurvatureMatirx.png

[plot2]: https://github.com/timohyva/AB_PhaseTransition/blob/master/Lotynk_Data_Plot_fafab.png 

[QUEST_DMC_WP2]: https://github.com/timohyva/AB_PhaseTransition/blob/master/QUEST_DMC1.png
