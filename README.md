Infotaxis Library
=================

This library implements the infotaxis strategy for autonomous search for a diffusive source without gradients <sup>\[1\]\[2\]</sup>. This library provides a class infotaxis::InfotaxisGrid which maintains a set of parameters (diffusivity, wind direction and velocity, etc...) as well as a probability field, allowing to update said field step by step. It also provides a method getOptimalMove to calculate which direction (4 directions + staying still) will likely provide the most information (maximum entropy reduction).

A main() is also provided as an example. It instanciates an InfotaxisGrid, then simulates a robot using the grid to find the diffusive source. It will output frames of each iteration to a ./pictures folder (that you'll need to create), representing the path since start on the current probability field.

Building
--------

infotaxis.cpp requires C++11 compiling but no external library. main.cpp requires png++ for frames output. Then just do

	cmake .
	make

Examples
--------

Here are some gifs compiled from aforementioned frames. Common parameters are :
Diffusivity = 1; Emission rate = 1; Expected lifetime of the particles = 400; Delta time = 1; Grid size = 100x100; resolution = 0.05;

Wind velocity = 1, angle = -Pi/2 (southbound), sensor radius = 1

![gif](http://i.imgur.com/El4AEaD.gif)

Same wind, sensor radius = 0.1

![gif](http://i.imgur.com/vyUG5eG.gif)

-------------
<sub>[1] RistiCampus c, B., Skvortsov, A., & Walker, A. (2014). Autonomous search for a diffusive source in an unknown structured environment. Entropy, 16, 789–813. doi:10.3390/e16020789</sub>

<sub>[2] Vergassola, M., Villermaux, E., & Shraiman, B. I. (2007). “Infotaxis” as a strategy for searching without gradients. Nature, 445(7126), 406–409.</sub>