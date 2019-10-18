# Comparison of N-methylacetamide conformation between Classical MD and QM/MM/MD simulations

**Description:** Several QM/MM tutorials adapted to CP2K. Tutorials adapted from AMBER, GMX, NAMD and CPMD.

**Author:** Salomé Llabrés Prat, PhD

## Original tutorial from AMBER suite
Tutorial developed by Dr. Ross Walker for the [AMBER 15 software](http://ambermd.org). The goal of the tutorial is to do fast and accurate coupled potential semi-empirical QM/MM simulations with full treatment of long range electrostatics (PME) or Generalized Born implicit solvent simulations. 

**Reference:** http://ambermd.org/tutorials/advanced/tutorial2/index.htm

<br/><br/>

---

# Tutorial adapted to CP2K
In this tutorial we will run two simulations of a toy system using two different schemes (classical MD simulation and QM/MM MD simulation usign a semi-empirical PM3 Hamiltonian) and compare the results. The toy system we are going to use is a simple molecule as N-methylacetamide (**NMA**) in a periodic box of TIP3P water. 
The tutorial is structured into 4 different sections:
- **Section 1:** Setup of the system. 
- **Section 2:** Classical MD simulations
- **Section 3:** Semi-empirical QM/MM simulations
- **Section 4:** Comparison of results

<br/><br/>

---

## Section 1: Setup of the system
The first stage is to build a topology and coordinate file for N-methylacetamide (**NMA**) in TIP3P waterbox. NMA can be built  simply by joining an ACE and an NME residue. There are several options to do so. Here we propose two:
- Create a *fake* initial PDB structure containing 3 atoms of each the residues involved. Leap (xLeap or tLeap) will add the missing atoms for us.  
- Using pymol to build it using the Builder command. 

### 1.1 *Fake* PDB input 
Here we provide a *fake* PDB file that specifies the first atom of the two residues. 
**NMA_skeleton.pdb**
```javascript
ATOM      1  CH3 ACE     1      -0.410   0.146   0.154  1.00  0.00           C
ATOM      2  C   ACE     1       1.048   0.238  -0.245  1.00  0.00           C
ATOM      3  O   ACE     1       1.659  -0.755  -0.662  1.00  0.00           O
ATOM      4  N   NME     2       1.662   1.368  -0.150  1.00  0.00           N
ATOM      5  CH3 NME     2       3.061   1.468  -0.530  1.00  0.00           C
ATOM      6  H   NME     2       1.171   2.183   0.190  1.00  0.00           H
TER
```

### 1.2 Pymol commands

Use the [Pymol](https://sourceforge.net/projects/pymol/) builder tool, here is a snapshot of the pymol screenshot highlighting the Builder menu. Save the molecule under the name of **NMA_skeleton.pdb**. 

![pymol_screenshot](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/images/Pymol_screenshot.png)


Once we have the NMA_skeleton.pdb file, we prepare the **leap.in** script.
```javascript
source leaprc.protein.ff14SB
source leaprc.water.tip3p
NMA = loadpdb NMA_skeleton.pdb
solvatebox NMA TIP3PBOX 15
savepdb NMA NMA.pdb
saveamberparm NMA NMA.prmtop NMA.inpcrd
quit
```
Run the following command:
```javascript
$ tleap -f leap.in
```

You should have created 3 files:
- **NMA.pdb**: PDB file to visualise the system. 
- **NMA.prmtop**: Topology file. 
- **NMA.inpcrd**: coordinates. 

> TIP: It is good practice to visualise the PDB files created by the Leap before starting any simulations.

<br/><br/>

---

## Section 2: Classical MD simulations
Now that we have our topology and coordinate files, we are ready to carry out classical molecular dynamics. First, we need to minimize our system to remove any bad contacts created by solvation and also to allow our artificially created NMA structure to relax. For this we will run 4500 steps of minimization.

<br/><br/>

### 2.1 Minimisation

Our input file is very similar to files we have used in the previous tutorials. Since this is a periodic simulation with PME we can safely use a cut off of 8 angstroms. **WIP**

You will find the input file here: **AMBER_NMA/section2/CP2K/cp2k_mm_minimisation.inp**

** Highlighted regions of cp2k_mm_minimisation.inp**

We have to set up calculation type as `GEO_OPT` in the `&GLOBAL` section:

```javascript
&GLOBAL
  PROJECT EM          ! Name of the calculation
  PRINT_LEVEL LOW     ! Verbosity of the output
  RUN_TYPE GEO_OPT    ! Calculation type: Geometry optimisation
&END GLOBAL
```

In the `&FORCE_EVAL` section, we define the basci system definitions (such as topology and coordinates). Here we are only going to explain the most important ones:
- `METHOD FIST` : Classical MD method in CP2K
- In the `&MM` subsection, we set the forcefield type `PARMTYPE AMBER` and the topology file `PARM_FILE_NAME NMA.prmtop`. 
- In the `&SUBSYS` subsection, we define the box size in the `&CELL` subsection and both AMBER topology `CONN_FILE_NAME NMA.prmtop` and AMBER COORDINATES `COORD_FILE_NAME NMA.inpcrd`. 
have to set up calculation type as `GEO_OPT` in the `&GLOBAL` section:

```javascript
&FORCE_EVAL
  METHOD FIST                     ! FIST Method for the MM calculations
  &MM
    &FORCEFIELD                   ! Force_field for the classical calculations
      PARMTYPE AMBER              ! Forcefield format : AMBER
      PARM_FILE_NAME NMA.prmtop   ! Topology file. 
      &SPLINE                     ! Parameters for the splines used in the nonboned interactions
        EMAX_SPLINE 1.0E8         ! Maximum value of the potential up to which splines will be constructed
        RCUT_NB [angstrom] 8      ! Cutoff of non-bonding interactions
      &END SPLINE
    &END FORCEFIELD
    &POISSON                      ! Poisson resolutor
      &EWALD                      ! Ewald parameters controlling electrostatic only for CLASSICAL MM
        EWALD_TYPE SPME           ! The type of ewald: smooth particle mesh using beta-Euler splines
        ALPHA .40                 ! Alpha parameter
        GMAX 80                   ! Number of grid points
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS                                                  ! Subsystem: coordinates, topology, molecules and cell
    &CELL                                                  ! Input parameters needed to set up the CELL
      ABC [angstrom] 40.1760410  39.3025410  39.5609130.   ! Set box dimensions here
      ALPHA_BETA_GAMMA 90 90 90                            ! Box shape
    &END CELL
    &TOPOLOGY                                              ! Section specifying topology for classical runs.
      CONN_FILE_FORMAT AMBER                               ! Connectivity file format : AMBER
      CONN_FILE_NAME NMA.prmtop                            ! Connectivity file name
      COORD_FILE_FORMAT CRD                                ! Coordinates file format : CRD
      COORD_FILE_NAME NMA.inpcrd                           ! Coordinates file name
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

> TIP: CP2K dooes not read the box size from the AMBER CRD file. You **MUST** copy the box size from the last line of the CRD file. 


In the `&MOTION` section, we define the parameters for the geometry optimisation. 
```javascript
&MOTION
  &GEO_OPT                             ! Geometry optimisation
    OPTIMIZER LBFGS                    ! Geometry optimisation algorithm
    MAX_ITER 4500                      ! Maximum number of iterations
  &END
  &PRINT                               ! Printing properties during a geometry optimization run
    &TRAJECTORY                        ! Controls the output of the trajectory
      FORMAT XYZ                       ! Format of the output trajectory:  XYZ
      &EACH                            ! New trajectory frame will be printed each 500 steps
          GEO_OPT 500
      &END EACH
    &END TRAJECTORY
    &RESTART                           ! This section controls the printing of restart files
      &EACH                            ! A restart file will be printed every 500 steps
        GEO_OPT 500
      &END EACH
    &END RESTART
    &RESTART_HISTORY                   ! This section controls dumping of unique restart files during the run keeping all of them. Most useful if recovery is needed at a later point.
      &EACH                            ! A new restart file will be printed every 500 md steps
        GEO_OPT 500
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

Then run the following command: 

```javascript
$ cp2k.popt cp2k_mm_minimisation.inp > cp2k_mm_minimisation.out
```

Once the job has finished, we need to check that the minimisation has converged. For each optimisation step, we have this section: 

```javascript
  --------------------------
 OPTIMIZATION STEP:    368
 --------------------------

 ENERGY| Total FORCE_EVAL ( FIST ) energy (a.u.):            -28.450804265384360


 --------  Informations at step =   368 ------------
  Optimization Method        =                LBFGS
  Total Energy               =       -28.4508042654
  Real energy change         =        -0.0039438240
  Decrease in energy         =                  YES
  Used time                  =                0.174

  Convergence check :
  Max. step size             =         0.0099827301
  Conv. limit for step size  =         0.0100000000
  Convergence in step size   =                  YES
  RMS step size              =         0.0007891138
  Conv. limit for RMS step   =         0.0050000000
  Convergence in RMS step    =                  YES
  Max. gradient              =         0.0061057290
  Conv. limit for gradients  =         0.0100000000
  Conv. in gradients         =                  YES
  RMS gradient               =         0.0002605894
  Conv. limit for RMS grad.  =         0.0050000000
  Conv. in RMS gradients     =                  YES
 ---------------------------------------------------
```

Here we can see that it took 368 steps to converge this system. CP2K evaluates the convergence using 4 different criteria:

- The calculated displacement for the next step must be essentially 0:

    **Conv. limit for step size**  =        4.50000000E-004
  
- The root-mean-square (RMS) of the calculated displacement for the next step must be also essentially 0:

    **Conv. limit for RMS step**   =        1.50000000E-003

- The forces/gradient must be essentially 0: 

    **Conv. limit for gradients**  =        3.00000000E-003

- The RMS of the forces/gradient must be also essentially 0:

    **Conv. limit for RMS grad.**  =        3.00000000E-004
  

Once all four criteria are met, the optimisation phase stops and CP2K reevaluates the energy at the minimum and finishes the calculation:

```javascript
 *******************************************************************************
 ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
 *******************************************************************************

                    Reevaluating energy at the minimum

 ENERGY| Total FORCE_EVAL ( FIST ) energy (a.u.):            -28.450804265384360
```

<br/><br/>

### 2.2 MD (without equilibration)
Once we have minimized the system, we will run 1,000 fs of molecular dynamics at a constant temperature of 300K. Since our toy system is so simple we can do so. However, normally we would slowly heat our system to 300K and run a short MD both at NPT and NVT to make sure the volume and pressure are correct. 

You will find the input file here: **AMBER_NMA/section2/CP2K/cp2k_mm_npt.inp**

**Highlighted regions of cp2k_mm_npt.inp**
Here we will explain how the CP2K input changes compared to the minimisation input file. 

We have to set up calculation type as `MD` in the `&GLOBAL` section:

```javascript
&GLOBAL
  PROJECT NPT          ! Name of the calculation
  PRINT_LEVEL LOW      ! Verbosity of the output
  RUN_TYPE MD          ! Calculation type: Molecular Dynamics
&END GLOBAL
```

The main changes are located in the `&MOTION` section:
- `&MD` subsection : Define all the parameter to run a MD. Including:
  - `ENSEMBLE` : NPT_I, NVT...
  - `&BAROSTAT` subsection 
  - `&THERMOSTAT` subsection   

```javascript
&MOTION
  &MD                       ! Parameters needed perform an MD run.
    ENSEMBLE NPT_I          ! Ensemble/integrator to use for MD propagation
    TIMESTEP [fs] 0.5       ! Integration step
    STEPS    20000          ! The number of MD steps to perform
    TEMPERATURE 298         ! Temperature (K)
    &BAROSTAT 
      TIMECON [fs] 100      ! Barostat time constant
      PRESSURE [bar] 1.0    ! Initial pressure
    &END BAROSTAT
    &THERMOSTAT
      REGION GLOBAL         ! Region each thermostat is attached to.
      TYPE CSVR             ! Canonical sampling through velocity rescaling thermostat
      &CSVR
        TIMECON [fs] 10.    ! Time constant of the CSVR thermostat.
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &RESTART               ! Controls the printing of restart files
      &EACH                ! A restart file will be printed every 10000 md steps
        MD 25000
      &END
    &END
    &TRAJECTORY            ! Controls the output of the trajectory
      FORMAT DCD           ! Format of the output trajectory is DCD
      &EACH                ! New trajectory frame will be printed each 100 md steps
        MD 500
      &END
    &END
    &RESTART_HISTORY       ! Controls dumping of unique restart files during the run keeping all of them.
      &EACH                ! A new restart file will be printed every 10000 md steps
        MD 5000
      &END
    &END
  &END PRINT
&END MOTION
```

We have an additional section `&EXT_RESTART` that enable to restart the simulation from the output of the minimisation step.
``` javascript
&EXT_RESTART                        ! External restart, specifies an external input file where to take positions
  RESTART_FILE_NAME EM-1.restart    ! Specifies the name of restart file
  RESTART_COUNTERS .FALSE.          ! Restarts the counters in MD schemes and optimization STE
&END
```

Then run the following command: 

```javascript
$ cp2k.popt cp2k_mm_npt.inp > cp2k_mm_npt.out
```

This simulation produces 3 output files (as we have defined in the `&MOTION/&PRINT` section of the input file:
- **NPT-1.cell**: Time evolution of the cell size.
- **NPT-1.ener**: Time evolution od the system energy.
- **NPT-pos-1.dcd**: Resulting trajectory of the simulation.
- **NPT-1.restart**: Restart file.

The CP2K output lists the system properties at each simulation step. It is worth check this information, it can inform us that something weird is going on. 

```javascript
ENERGY| Total FORCE_EVAL ( FIST ) energy (a.u.):            -21.356208918405574

*******************************************************************************
 ENSEMBLE TYPE                =                                            NPT_I
 STEP NUMBER                  =                                             2000
 TIME [fs]                    =                                      1000.000000
 CONSERVED QUANTITY [hartree] =                              -0.222957087224E+02

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                        0.20                 0.22
 ENERGY DRIFT PER ATOM [K]    =          0.534744069806E+01   0.514644830721E+01
 POTENTIAL ENERGY[hartree]    =         -0.213562089184E+02  -0.213438912680E+02
 KINETIC ENERGY [hartree]     =          0.613345909757E+01   0.611129281552E+01
 TEMPERATURE [K]              =                     294.860              293.795
 PRESSURE [bar]               =         -0.611160475868E+02  -0.135378276079E+01
 BAROSTAT TEMP[K]             =          0.449844101096E+03   0.144704296670E+04
 VOLUME[bohr^3]               =          0.391890442973E+06   0.423127981084E+06
 CELL LNTHS[bohr]             =    0.7409758E+02   0.7248656E+02   0.7296308E+02
 AVE. CELL LNTHS[bohr]        =    0.7600430E+02   0.7435183E+02   0.7484061E+02
 *******************************************************************************
```

We should also check the temperature, both the kinetic energy and potential energy of the system and the simulation box size. If we plot them against simulation time, we observe that the system has convergenced. 

Temperature | Kinetic Energy | Potential Energy | Simulation box volume
------------ | ------------- | ------------- | -------------
![MM_temperature](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section2/CP2K/Temp.png) | ![MM_Kinetic Energy](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section2/CP2K/KineticE.png) | ![MM Potential Energy](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section2/CP2K/PotentialE.png) | ![MM_box_volume](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section2/CP2K/Volume.png) 

<br/><br/>

---

## Section 3: Semi-empirical QM/MM simulations
Now, let's try repeating the simulation using a coupled QM/MM simulation. To model the NMA molecule, we are going to use a semi-empirical PM3 Hamiltonian. We are going to treat the solution classically. In this situation, we have a clear separation between the QM region and the MM region (no bonds crossing the QM/MM boundary). There is no need to define link atoms.

![QM/MM setup](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/images/QMMM.setup.png)

<br/><br/>

### 3.0 Modification of the Lennard-Jones parameters: 
We need to be aware of the limitations of the classical models and how to combine them with the QM methods. We have set up the solvent using TIP3P water model, which has no Lennard-Jones parameters defined for the hydrogen atoms. We have to set them in order to avoid unrealistic interactions with the QM region. Here in this example, we are going to use the LJ parameters in the GAFF2 forcefield. 

To modify them, there are two options: modifying the prmtop file using [parmed](https://parmed.github.io/ParmEd/html/index.html) or adding the missing LJ parameters into the CP2K input file. Here in this tutorial we are going to use the first option.

- Modify directly the LJ parameters using Parmed. Here are the PARMED commands to run:
```javascript
$ parmed NMA.prmtop
```

**parmed commands**
```javascript
changeLJSingleType :WAT@H1 0.3019 0.047
changeLJSingleType :WAT@H2 0.3019 0.047
outparm NMA.LJ_mod.prmtop
quit
```

- Adding the missing parameters into the CP2K parameters into the `&FORCEFIELD/&NONBONDED` section of the CP2K input file.

```javascript
&FORCEFIELD
  ...
  &NONBONDED
    &LENNARD-JONES
        ATOMS HW N
        EPSILON [kcalmol] 0.052
        SIGMA [angstrom]  2.42
        RCUT [angstrom] 9.0
    &END
    &LENNARD-JONES
        ATOMS HW O
        EPSILON [kcalmol] 0.058
        SIGMA [angstrom]  2.2612
        RCUT [angstrom] 9.0
    &END
  &END
  ...
&END
```

<br/><br/>

### 3.1 Minimisation
Our input file is very similar to files we have used in the previous tutorials. Since this is a periodic simulation with PME we can safely use a cut off of 8 angstroms. **WIP** 

You will find the input file here: **AMBER_NMA/section3/CP2K/cp2k_qmmm_minimisation.inp**

**Highlighted regions of cp2k_qmmm_minimisation.inp**

The most important changes on the CP2K input file are located in the `&FORCE_EVAL` section:
- `METHOD QMMM` as we are doing a QM/MM simulation.
- In the `&DFT` subsection, we define all the important parameters for the QM region:
  - `CHARGE 0`
  - `&QS` subsection : CP2K uses the QUICKSTEP method to run QM simulations.
    - `METHOD PM3` the semiempirical PM3 here. 
    - `&SE/&COULOMB` subsection we can set up tbe cutoff for Coulomb interactions.
    - `&SE/&EXCHANGE` subsection we can set up tbe cutoff for Exchange interactions.
  - `&SCF` subsection.    
- In the `&QMMM` subsection, we define all the important parameters for the QM/MM interaction:
  - `ECOUPL COULOMB` Energy Coupling
  - `&CELL` subsection defines the box for the QM region. 
  - `&QM_KIND element` We have to add a section for each element we have in the QM region and list their corresponding `MM_INDEX` indexes. 
  
```javascript
&FORCE_EVAL
  METHOD QMMM                                ! Method: QM/MM 
  STRESS_TENSOR ANALYTICAL                   ! Controls the calculation of the stress tensor. 
  &DFT                                       ! DFT programs
    CHARGE 0                                 ! The total charge of the system
    &QS         
      METHOD PM3                             ! Electronic structure method 
      &SE                                    ! Parameters needed for Semi-empirical methods
         &COULOMB                            ! Evaluation of the COULOMB term in SE calculations
           CUTOFF [angstrom] 10.0
         &END
         &EXCHANGE                           ! Evaluation of the EXCHANGE and core Hamiltonian terms in SE calculations
           CUTOFF [angstrom] 10.0
         &END
      &END
    &END QS
    &SCF                                     ! Perform an SCF run.
      MAX_SCF 30                             ! Maximum number of SCF iteration
      EPS_SCF 1.0E-6                         ! Target accuracy for the SCF convergence.
      SCF_GUESS ATOMIC                       ! Initial guess for the wavefunction (ATOMIC: atomic density using the atomic code)
      &OT                                    ! orbital transformation (OT) method.
        MINIMIZER DIIS                       ! Minimizer (DIIS: Direct inversion in the iterative subspace)
        PRECONDITIONER FULL_SINGLE_INVERSE   ! Preconditioner to be used with all minimization schemes.
      &END
      &OUTER_SCF                             ! Parameters controlling the outer SCF loop
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END
    &END SCF 
  &END DFT
  &MM
    &FORCEFIELD
      PARMTYPE AMBER
      PARM_FILE_NAME NMA.LJ_mod.prmtop
      &SPLINE
        EMAX_SPLINE 1.0E8
        RCUT_NB [angstrom] 10
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA .40
        GMAX 80
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS
    &CELL
      ABC [angstrom] 40.1760410  39.3025410  39.5609130
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &TOPOLOGY
      CONN_FILE_FORMAT AMBER
      CONN_FILE_NAME NMA.LJ_mod.prmtop
      COORD_FILE_FORMAT CRD
      COORD_FILE_NAME NMA.inpcrd
    &END TOPOLOGY
  &END SUBSYS
  &QMMM                                        ! Input for QM/MM calculations
    ECOUPL COULOMB                             ! Type of the QM - MM electrostatic coupling
    &CELL                                      ! Parameters needed to set up the CELL
      ABC 40 40 40
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &QM_KIND N                                 ! QM kind in the QM/MM scheme: element
      MM_INDEX 7                               ! The indexes of the MM atoms that have this kind
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 6
    &END QM_KIND
    &QM_KIND C
      MM_INDEX 2 5 9
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 1 3 4 8 10 11 12
    &END QM_KIND
  &END QMMM
&END FORCE_EVAL
```

Then run the following command: 

```javascript
$ cp2k.popt cp2k_qmmm_minimisation.inp > cp2k_qmmm_minimisation.out
```

Once the simulation has finished, we can analyse the results from this minimisation using the same protocol than in the **section 2.1**. It is worth noting that the ouput for every optimisation step looks like this:

```javascript
--------------------------
 OPTIMIZATION STEP:    515
 --------------------------

  Translating the system in order to center the QM fragment in the QM box.

 Number of electrons:                                                         30
 Number of occupied orbitals:                                                 15
 Number of molecular orbitals:                                                15

 Number of orbital functions:                                                 27
 Number of independent orbital functions:                                     27

 Extrapolation method: ASPC


 SCF WAVEFUNCTION OPTIMIZATION

  ----------------------------------- OT ---------------------------------------
  Minimizer      : DIIS                : direct inversion
                                         in the iterative subspace
                                         using   7 DIIS vectors
                                         safer DIIS on
  Preconditioner : FULL_SINGLE_INVERSE : inversion of 
                                         H + eS - 2*(Sc)(c^T*H*c+const)(Sc)^T
  Precond_solver : DEFAULT
  stepsize       :    0.08000000                  energy_gap     :    0.08000000
  eps_taylor     :   0.10000E-15                  max_taylor     :             4
  ----------------------------------- OT ---------------------------------------

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 OT DIIS     0.80E-01    0.0     0.00895204       -34.4682165326 -3.45E+01
...
    15 OT DIIS     0.80E-01    0.0     0.00000084       -34.4731676669 -3.18E-10

  *** SCF run converged in    15 steps ***


  Core-core repulsion energy [eV]:                           2336.86507766754858
  Core Hamiltonian energy [eV]:                             -1441.45102479507000
  Two-electron integral energy [eV]:                        -3630.00026212931107
  Electronic energy [eV]:                                   -3256.45115585972553
  QM/MM Electrostatic energy:                                  -0.67899964369650

  Total energy [eV]:                                         -938.06259813526765

  Atomic reference energy [eV]:                               917.68881409010157
  Heat of formation [kcal/mol]:                              -469.83065871492789

  outer SCF iter =    1 RMS gradient =   0.84E-06 energy =        -34.4731676669
  outer SCF loop converged in   1 iterations or   15 steps


 ENERGY| Total FORCE_EVAL ( QMMM ) energy (a.u.):            -62.312895111901369


 --------  Informations at step =   515 ------------
  Optimization Method        =                LBFGS
  Total Energy               =       -62.3128951119
  Real energy change         =        -0.0007502778
  Decrease in energy         =                  YES
  Used time                  =                0.385

  Convergence check :
  Max. step size             =         0.0497827994
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                   NO
  RMS step size              =         0.0040703138
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                   NO
  Max. gradient              =         0.0043510329
  Conv. limit for gradients  =         0.0004500000
  Conv. for gradients        =                   NO
  RMS gradient               =         0.0001596248
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. in RMS gradients     =                  YES
 ---------------------------------------------------
```

<br/><br/>

### 3.2 MD without equilibration

Once we have minimized the system, we will run 1,000 fs of molecular dynamics at a constant temperature of 300K. Since our toy system is so simple we can do so. However, normally we would slowly heat our system to 300K and run a short MD both at NPT and NVT to make sure the volume and pressure are correct. 

You will find the input file here: **AMBER_NMA/section3/CP2K/cp2k_qmmm_npt.inp**

In the **cp2k_qmmm_npt.inp** we combine the CP2K sections explained in the 2.2 section with the `&FORCE_EVAL` section explained in the 3.1 section above. 

Then run the following command:

```javascript
cp2k.popt cp2k_qmmm_npt.inp > cp2k_qmmm_npt.out
```

To ensure that the system is converged, we perform the same basic analysis explained in **Section2.2**. 

Temperature | Kinetic Energy | Potential Energy | Simulation box volume
------------ | ------------- | ------------- | -------------
![MM_temperature](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section3/CP2K/Temp.png) | ![MM_Kinetic Energy](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section3/CP2K/KineticE.png) | ![MM Potential Energy](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section3/CP2K/PotentialE.png) | ![MM_box_volume](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section3/CP2K/Volume.png) 

<br/><br/>

---

## Section 4: Comparison of results

If you use VMD to visualise the trajectory, we can observe the behaviour of the NMA molecule in water solution using classical treatment ot QM/MM ensemble. We think it is worth to monitor the following: 

- Visual analysis of the simulations
- O-C-N-H out of plane angle as a function of time.
- Cluster analysis 
- hydrogen bond analysis. 

There are many ways of doing such analysis. You can use the tools builtin AMBER suite to do this analysis. [cpptraj](http://ambermd.org/tutorials/analysis/) is able to read mdcrd and netcdf files (AMBER trajectories) and DCD (NAMD trajectories) and uses prmtop topology files (AMBER topologies). Also the same analysis can be done using other software such as [MDAnalysis](https://www.mdanalysis.org) or [MDtraj](http://mdtraj.org/1.9.3/). Some visualisation may help too: [VMD](https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node128.html).

Here we are going to use MDtraj to do so, we invite you to choose your analysis tool of preference and try to perform the same analysis. 

You will find the a jupyter notebook here: **AMBER_NMA/section4/Analysis.ipynb**

Here, you'll see the dihedral angle formed by O-C-N-H atoms for the simulations in section 2 (blue) and in section 3 (red). 

 Analysis | MM trajectory | QM/MM trajectory 
------------ | ------------ | ------------- 
Visual analysis (degrees) | ![MM_snapshot](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section4/mm.png) | ![QMMM_snapshot](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section4/qmmm.png)
Dihedral analysis (radians) | ![MM_dihedral](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section4/mm.dihedral.png) | ![QMMM_dihedral](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/AMBER_NMA/section4/qmmm.dihedral.png) 

<br/><br/>

