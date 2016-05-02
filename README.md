Numerical code used in the computation of absorption rates of surface polartion (SP(h)P) modes.
At the moment the code includes plane wave SP(h)P as well as OAM carrying vortex SP(h)P modes. Extension to new SP(h)P modes should be straightforward by specifying the vector potential of the modes of interest.

Implementation notes:
 * Throughout the program, the units of dimention are normalized by the Bohr radius (a0) which is the natural lengthscale of the electronic system.
 * The gradient of the wavefunction is computed numerical using the dimentionless parametar delta, which was found to be numerically robust for the calculations.
 * The limit of integration of the radial dimention is defined through the dimentionless input parameter R. Due to the exponentially decaying nature of the wavefunction, the parameter R=300 is great enough to correctly compute the overlap integral for common hydrogen wavefunctions.


The codes for LaguerreGen.m (author Matthias.Trampisch@rub.de) and legendre_newP.m are not of my authorship but are included here since they are necessary in the computation of the absorption rates.

The code for AbsorptionRate_Final.m was developed by Francisco Machado (fmachado@mit.edu) in the context of his work in [Marin Soljačić's group](http://www.mit.edu/~soljacic/) at MIT. 