Infotaxis Library
=================

This library implements the infotaxis strategy for autonomous search for a diffusive source without gradients <sup>\[1\]\[2\]</sup>. This library provides a class infotaxis::InfotaxisGrid which maintains a set of parameters (diffusivity, wind direction and velocity, etc...) as well as a probability field, allowing to update said field step by step. It also provides a method getOptimalMove to calculate which direction (4 directions + staying still) will likely provide the most information (maximum entropy reduction).

A main() is also provided as an example. It instanciates an InfotaxisGrid, then simulates a robot using the grid to find the diffusive source. It will output frames of each iteration to a ./pictures folder (that you'll need to create), representing the path since start on the current probability field.

Building
--------

The sources require C++11 compiling and png++/libpng, additionally main.cpp requires boost. Then just do

	cmake .
	make
	
CMake will check for libpng, but won't check for png++'s header's presence, and so trying to compile without it will fail in make.
	
### Without libpng

libpng is strictly required for main.cpp, but you can build the library without it. If libpng isn't found or if you manually specify -DINFOTAXIS_NOPNG to cmake, main.cpp won't be compiled, and infotaxis.cpp will be compiled without png-related methods.

True vs Fast Infotaxis
----------------------

One goal of this library was also to test several methods for driving the robot, i.e selecting the direction to go toward next.

- The "true" infotaxis method (as exposed in the referenced articles) implies calculating the expected encounter rate at the given neighboring cell, and from that, calculating the expected diminution of entropy, and going wherever this value is the highest in absolute value. This process takes a hell lot of processing, since it implies duplicating the grid, and simulating several cases of detection , in all 5 directions (4 sides, and staying in place).
- An alternative method, which is used if an InfotaxisGrid is created with `trueInfotaxis=false`, goes where detection is expected to be the most probable. So it computes the expected encounter rate like previously, and simply goes where this value is maximized.

You can test the alternative method in *main* by passing `-f` / `--fast` . Speedups of 5x~10x per iteration are observed on a 100x100 grid, however the increase in search time, as well as its robustness, are yet to experiment for, particularly in an increased grid size.

Examples
--------

Here are some gifs compiled from aforementioned frames. Common parameters are :
Diffusivity = 1; Emission rate = 1; Expected lifetime of the particles = 400; Delta time = 1; Grid size = 100x100; resolution = 0.05;

Wind velocity = 1, angle = -Pi/2 (southbound), sensor radius = 1

![gif](http://i.imgur.com/El4AEaD.gif)

Same wind, sensor radius = 0.1

![gif](http://i.imgur.com/vyUG5eG.gif)

Here's what happens when it doesn't detect anything (Artificially. With no wind)

![gif](http://i.imgur.com/QyuMeP6.gif)

-------------
<sub>[1] RistiCampus c, B., Skvortsov, A., & Walker, A. (2014). Autonomous search for a diffusive source in an unknown structured environment. Entropy, 16, 789–813. doi:10.3390/e16020789</sub>

<sub>[2] Vergassola, M., Villermaux, E., & Shraiman, B. I. (2007). “Infotaxis” as a strategy for searching without gradients. Nature, 445(7126), 406–409.</sub>