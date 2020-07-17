A simple chemical enrichment model for galaxies including stellar population
modelling in a star-by-star fashion. A complimentary semi-analytic model
to an under-development chemical enrichment and individual star model in
the cosmological hydrodynamics code Enzo. Intended use (for now) is to quickly
and cheaply characterize the effect of modifying free parameters in the 
chemical and stellar model without running expensive hydro simulations.

Current state (Fall 2019):

Many improvements to code has been made, with some examples on how to run
many models in parallel to speed up tests of various parameters and things.
Currently working on porting the slowest chunks of the code to cython to 
get some actual speed improvements to be able to run this for more useful
science. Things are still very rough, with changes to the API at all levels
still likely (backwards compatability will not be garunteed for a while). 
Code ported to Python3 now

Prior state (Oct 2016):

As is, the star-by-star modelling in this onezone code is highly inefficient
(there is a reason why SSP's is the standard protocol!!!). However, we
only follow star-by-star at the moment to be a direct analog of the much more
expensive hydro simulations. At the moment, running many calculations is 
computationally infeasible, as a few Gyrs takes on order of a few hours 
(depending strongly on the number of particles). This can be improved immensely
with some parallelization and optimzation work, that I hope to get around to
at some point. The code currently works, though is not currently user friendly
(unless you happen to be the one who wrote it). Improving this is on my list
of to-do's.

See here for my current fork of Enzo: https://bitbucket.org/aemerick/enzo-emerick
