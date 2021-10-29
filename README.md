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
[GPV2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.en.html)

## Running The codes

### Python ![alt text](https://github.com/timohyva/AB_PhaseTransition/blob/master/logo_languge1.png)
The defult way is launching from bash. And it's ok for rungning in IDE.
```bush
# import strong coupling correction module. there is a class in this modeule
~$ python3 import Module_SC_CorrectionObject_V01 as SC 

# calculate and show the A-B equilibrium line 
~$ python3 ABphaseFreeEnergyDiff_DensityContorPlot_TemperatureChange_withBetaObject_V3.py > out.out 

# curvature countour plot
~$ python3 plot_curvatures_pressureFixed_APhase_ContourPlot_V0.py > out.out
```

### Mathematica
Same with any other .nb script. No powerful resourse is needed.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Alternatively, for H1 and H2, an underline-ish style:

Alt-H1
======

Alt-H2
------

Emphasis, aka italics, with *asterisks* or _underscores_.

Strong emphasis, aka bold, with **asterisks** or __underscores__.

Combined emphasis with **asterisks and _underscores_**.

Strikethrough uses two tildes. ~~Scratch this.~~

1. First ordered list item
2. Another item
⋅⋅* Unordered sub-list. 
1. Actual numbers don't matter, just that it's a number
⋅⋅1. Ordered sub-list
4. And another item.

⋅⋅⋅You can have properly indented paragraphs within list items. Notice the blank line above, and the leading spaces (at least one, but we'll use three here to also align the raw Markdown).

⋅⋅⋅To have a line break without a paragraph, you will need to use two trailing spaces.⋅⋅
⋅⋅⋅Note that this line is separate, but within the same paragraph.⋅⋅
⋅⋅⋅(This is contrary to the typical GFM line break behaviour, where trailing spaces are not required.)

* Unordered list can use asterisks
- Or minuses
+ Or pluses

[I'm an inline-style link](https://www.google.com)

[I'm an inline-style link with title](https://www.google.com "Google's Homepage")

[I'm a reference-style link][Arbitrary case-insensitive reference text]

[I'm a relative reference to a repository file](../blob/master/LICENSE)

[You can use numbers for reference-style link definitions][1]

Or leave it empty and use the [link text itself].

URLs and URLs in angle brackets will automatically get turned into links. 
http://www.example.com or <http://www.example.com> and sometimes 
example.com (but not on Github, for example).

Some text to show that the reference links can follow later.

[arbitrary case-insensitive reference text]: https://www.mozilla.org
[1]: http://slashdot.org
[link text itself]: http://www.reddit.com

Here's our logo (hover to see the title text):

Inline-style: 
![alt text](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 1")

Reference-style: 
![alt text][logo]

[logo]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 2"

# AB_PhaseTransition
numeric codes developments for A-B transiton models
