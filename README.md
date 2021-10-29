# AB_PhaseTransition wiki page
Repsositry and codes devoted on A-B transtion, Bubble Nucleation and Some Confinements Physics 

## Descriptions
This repositoty contaitains a coefficients obeject ([Module_SC](https://github.com/timohyva/AB_PhaseTransition/blob/master/Module_SC_CorrectionObject_V01.py)), which calculates the forth order coefficients of Ginzburg-Landau with Strong coupling correction.
The strong coupling correection model comes from Physical Review B papers by Robert.C.Regan, J.J.Wiman and James. A. Sauls:
1. *"Vortex phase diagram of rotating superfluid  3He−B"*
    [Vortex phase diagram of rotating superfluid 3He−B](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.024517)
2. *"Superfluid phases of 3He in nanoscale channels"* 
    [Superfluid phases of 3He in nanoscale channels](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.144515)

The codes are developed to implement different kinds of calculations involving SC-Ginzburg-Landau free energy and dydic tensor field/order parameter of symmetry-breking phases of liquid Helium-3 in mK level temperature. 

(* The implemet of physics introduced by confiment geometry will be considered/developed soon *)

## Languages (current)
* python3
* mathematica

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

### Mathematica ![Mathmatica](https://github.com/timohyva/AB_PhaseTransition/blob/master/logo_languge3.png)
Same with any other .nb script. No powerful resourse is needed.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
