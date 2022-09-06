# Heat accumulatio in pulsed laser processing (drilling)
![](https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/overview.gif)

## Overview
- This code was intended to be used for the calculation of heat accumulation in pulsed laser processing, primarily ultrafast drilling.
- The key idea is an efficient and flexible implementation of superposition of analytical solutions to the heat equation.
- Determination and superposition of these heatkernels is quite efficient - especially in the case of equidistant timestepping - which lends itself to an implementation based on fast convolution (frequency domain convolution).
- As is typical for these calculations, the temperature is evaluated immediately before the succeeding pulse.
- The solution is radially symmetric.
- The simulation domain is the half-infinite body.
- Sources are mirrored about the surface.
- Convection is not considered.
- Fluid dynamics cannot be considered in this formulation.
- Point / Gaussian and Circular "donut" sources are implemented.

## Features
- Pulse repetition rate can be constant or variable (chirped).
- Energy distribution over depth can be variable and time-dependent.
- A special case could be burst-processing.
- Irradiated / residual pulse energy can be time-dependent.
- Gaussian / "donut" source: Lateral extent of source can be time-dependent.
- In conclusion: The implementation is flexible enough to provide (in theory) a useful and fast approximation of heat accumulation in ultrafast laser processes - in particular percussion drilling.

%% Examples
- Drilling with an energy ramp
![https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/eramp.mp4](https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/eramp.gif)

- Drilling with a variable pulse repetition rate (chirp)
![https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/chirp.mp4](https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/chirp.gif)

- Burst processing
![https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/burst.mp4](https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/burst.gif)

- Not sure what to call this
![https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/fun.mp4](https://github.com/david-brinkmeier/heatAccumulation/blob/main/resources/media/fun.gif)

## Get Started
- Simply run [main.m](main.m)
- Beyond that you're unfortunately on your own, but in general I'm happy to answer question. Contact me.
- (The code is not pretty, but I think there are some neat ideas / implementations used in [heatacc](/functions/heatacc.m) and [get2D_Distributed_Heatkernels](/functions/get2D_Distributed_Heatkernels.m).

## Why wasn't this used?
- Unfortunately both the drilling progress as well as the energy distribution (i.e. absorbed energy over the borehole depth resulting in the temperature increase) is unknown.
- Even worse: Existing models that estimate the drilling progress are based on considerations (absorbed energy distributions) which are at odds with experimental observations.
- Determination of the energy distribution is in pulsed drilling is non-trivial and until that problem is solved these calculations are, realistically, nothing more than pretty to look at.