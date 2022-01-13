# AB_PhaseTransition wiki page

Repsositry and codes devoted on A-B transtion, Bubble Nucleation and Some Confinements Physics 

![WP2][QUEST_DMC_WP2]

## Descriptions
This repositoty contaitains a coefficients obeject ([Module_SC](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SC_CorrectionObject_V01.py)), which calculates the forth order coefficients of Ginzburg-Landau with Strong coupling correction. An other modeule named [Module_SC_Beta](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SC_Beta_V01.py) uses this coefficient object to provides all GL paramters, gaps energies and coherent length over p-T plane.
The strong coupling correection model comes from Physical Review B papers by Robert.C.Regan, J.J.Wiman and James. A. Sauls:
1. *"Vortex phase diagram of rotating superfluid  3He−B"*
    [Vortex phase diagram of rotating superfluid 3He−B](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.024517)
2. *"Superfluid phases of 3He in nanoscale channels"* 
    [Superfluid phases of 3He in nanoscale channels](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.144515)

The codes are developed to implement different kinds of calculations involving SC-Ginzburg-Landau free energy and dydic tensor field/order parameter of symmetry-breking phases of liquid Helium-3 in mK level temperature. 

(* The implement of physics introduced by confiment geometry will be considered/developed late *)

## Languages (current)
* python3


## License
![license](https://github.com/timohyva/AB_PhaseTransition/blob/master/GPL3logo2.png)
[GPV3] (https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.en.html)

## Running The codes

### Python ![Python](https://github.com/timohyva/AB_PhaseTransition/blob/master/logo_languge11.png)
The defult way is launching from bash. And it's ok for rungning in IDE.
```shell
# import strong coupling correction module. there is a class in this modeule
~$ python3 import Module_SC_CorrectionObject_V01 as SC 

# calculate and show the A-B equilibrium line 
~$ python3 ABphaseFreeEnergyDiff_DensityContorPlot_TemperatureChange_withBetaObject_V3.py > out.out 

# curvature countour plot
~$ python3 plot_curvatures_pressureFixed_APhase_ContourPlot_V0.py > out.out
```
![alt text][plot1]
![Lotynk_Data_Plot][plot2]

### Mathematica ![Mathmatica](https://github.com/timohyva/AB_PhaseTransition/blob/master/logo_languge3.png)
Same with any other .nb script. No powerful resourse is needed.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

[plot1]: https://github.com/timohyva/AB_PhaseTransition/blob/master/Contour_And_Density_Plot_Of_1st_EigenvalueOfCurvatureMatirx.png

[plot2]: https://github.com/timohyva/AB_PhaseTransition/blob/master/Lotynk_Data_Plot_fafab.png 

[QUEST_DMC_WP2]: https://github.com/timohyva/AB_PhaseTransition/blob/master/QUEST_DMC1.png
