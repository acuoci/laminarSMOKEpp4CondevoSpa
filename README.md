# laminarSMOKE++
CFD solver for laminar reacting flows with detailed kinetic mechanisms based on OpenFOAM and [OpenSMOKE++ framework][1]

If you use laminarSMOKE++ for your publications, we kindly ask you to cite the following papers:

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> Numerical modeling of laminar flames with detailed kinetics based on the operator-splitting method
> (2013) Energy and Fuels, 27 (12), pp. 7730-7753, DOI: 10.1021/ef4016334
 
> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> A computational tool for the detailed kinetic modeling of laminar flames: Application to C2H4/CH4 coflow flames
> (2013) Combustion and Flame, 160 (5), pp. 870-886, DOI: 10.1016/j.combustflame.2013.01.011
 
> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014

Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- Boost C++ (http://www.boost.org/)
- OpenSMOKE++ (provided with the current version of laminarSMOKE++)
- Intel MKL (https://software.intel.com/en-us/intel-mkl)

Optional libraries
------------------
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)

Optional libraries (under testing)
----------------------------------
- ISATLib (mauro.bracconi@polimi.it)

Compilation
-----------------------------------------------------
1. Open the `config.mkl` file and adjust the paths to the compulsory libraries
2. Type: `source config.mkl`
3. Compile the following laminarSMOKE++ libraries included in the `libs` folder through the usual `wmake` command: `boundaryConditions`, `clusteringAlgorithms`, `materialSynthesis`, and `radiationModels`
4. Compile the desired solver(s) included in the `applications` folder: for example, from the `applications/laminarBuoyantPimpleSMOKE++` folder, type `wmake` to compile the unsteady solver for buoyant flows

[1]: https://www.opensmokepp.polimi.it/
