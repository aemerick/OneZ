A simple chemical enrichment model for galaxies including stellar population
modelling in a star-by-star fashion. A complimentary semi-analytic model
to an under-development chemical enrichment and individual star model in
the cosmological hydrodynamics code Enzo. Intended use (for now) is to quickly
and cheaply characterize the effect of modifying free parameters in the 
chemical and stellar model without running expensive hydro simulations.

Current state:

As is, the star-by-star modelling in this onezone code is highly inefficient
(there is a reason why following SSP's is the standard protocol). However, we
only allow star-by-star at the moment to be a direct analog of the much more
expensive hydro simulations. At the moment, a few Gyr of evolution may take a 
day on a single processor. This could be significantly improved with some
optimization and parallelization
