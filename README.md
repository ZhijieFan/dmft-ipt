# dmft-ipt

An *iterative perturbation theory* (IPT) based solver for the Dynamical Mean-Field Theory. 

The code is based on:
* SciFortran [https://github.com/aamaricci/SciFortran]
* DMFT_Tools [https://github.com/aamaricci/DMFTtools]
and give  access to the solution of the equilibrium DMFT equations using either **Matsubara**, **Real-axis** and **Keldysh** formalism. 

The code structure is as follow:
* The set of modules compile into a top layer named `DMFT_IPT.f90`
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
In the driver code the user must includes the `DMFT_IPT` module and call the necessary procedures to solve the DMFT equations.

An example, solving the Hubbard model at finite temperature, is contained in the file `ipt_hm_matsubara.f90`.


***COPYRIGHT & LICENSING***  
Copyright 2012 - 2017 (c), Adriano Amaricci.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

