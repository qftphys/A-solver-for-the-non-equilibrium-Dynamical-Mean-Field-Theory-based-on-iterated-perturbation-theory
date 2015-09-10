# NON-EQUILIBRIUM DMFT-IPT

##### *Code is still in beta-developing state*

A solver for the non-equilibrium Dynamical Mean-Field Theory based on *iterated perturbation theory*. 

The code is based on:  
* SciFortran [https://github.com/aamaricci/SciFortran]  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools]

The code structure is as follow:  
* The set of modules compile into a top layer named `NEQ_DMFT_IPT.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* 
In the driver code the user must includes the `NEQ_DMFT_IPT` module and call the necessary procedures to solve the DMFT equations.

An example, solving the Hubbard model on the Bethe lattice, is contained in the file `neqdmft_bethe_quench.f90`.


--

***COPYRIGHT & LICENSING***  
Copyright 2012 - 2017 (c), Adriano Amaricci.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

