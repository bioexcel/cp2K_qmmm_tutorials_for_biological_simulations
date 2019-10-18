# Comparison of energy barriers for a Diels-Alder reaction in different environments. 

**Description:** Several QM/MM tutorials adapted to CP2K. Tutorials adapted from AMBER, GMX, NAMD and CPMD.

**Author:** Salomé Llabrés Prat, PhD

# Original tutorial from GMX suite
## 1. Diels-Alder Antibody Catalyst (DAA)
Tutorial developed by Dr. Gerrit Groenhof for the [GMX](http://www.gromacs.org/Documentation/Tutorials). The goal of this tutorial is to set up a pure QM simulation and QM/MM simulations that reproduce a Diels-Alder reaction in different environmnts: in vacuum, in solution and in complex with the protein. 

**Reference:** http://wwwuser.gwdg.de/%7Eggroenh/EMBO2004/html/introduction.html

<br/><br/>

---

# Tutorial adapted to CP2K
In this tutorial we will run three simulations of a Diels-Alder reaction in different environments: vaccum, in solution and in a biological environment. 

The tutorial is structured into 4 different sections:
- **Section 0:** Introduction.
- **Section 1:** Single point calculations of Diels-Alder reaction in vacuum.
- **Section 2:** Diels-Alder reaction in vacuum.
- **Section 3:** Diels-Alder reaction in solution.
- **Section 4:** Diels-Alder reaction complexed with the antibody.
- **Section 5:** Comparison of the energy barriers among the different environments.

<br/><br/>

---

## Section 0: Introduction
The essence of biological catalysis is the complementarity between the transition state of the process catalyzed and the enzyme. The transition state "fits" perfectly into the enzyme's active site pocket and is therefore stabilized, which enhances the reaction rate. 

Also, the essence of our immune system is also complementarity, but in this case it is the complementarity between the antigen, a compound that should not be inside your body, and the antibody. Upon exposure to an antigen the immune system generates immunoglobins that can bind the antigen and thereby make it harmless. This natural immuno response can be exploited to generate so-called antibody catalysts for chemical reactions. The immune system is triggered by exposing it to a compound that mimicks the transition state of the reaction of interest. Such a compound is called a transition state analogue. Even though it structurally resembles the transition state, it is a stable molecule. The generated immunoglobins will strongly bind to the analogue, and, as the analogue is structurally related to the transition state of the reaction, the immunoglobins will posses a certain degree of catalytic activity for that reaction.

The Diels Alder reaction we are going to study is the following:

![Diels Alder Reaction](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/F1.large.png)


Reference of the study: [Evolution of Shape Complementarity and Catalytic Efficiency from a Primordial Antibody Template](https://science.sciencemag.org/content/286/5448/2345.long)

**PDB id**: [1C1E](https://www.rcsb.org/structure/1c1e)

<br/><br/>

---

## Section 1: Single point calculations of Diels-Alder reaction in vacuum

First, we are going to perform geometry optimisation calculations on 3 structures: 

- Reactants (**R**), 
- Transition state (**TS**) and
- Products (**P**).

We have provided these 3 structures in XYZ format. These has been inferred from the X-ray structure (PDB id 1C1E) and previosly optimised using Gaussian09. 

<br/><br/>

### 1.1 Geometry optimisation of the chemical structures of the reactants and products

First, we are going to perform a QM geometry optimisation for the reactant and product states. We are going to use the semi-empirical PM3.

You will find the input file here: 
- **GMX_DAA/section1/CP2K/cp2k_R_geoopt.pm3.inp**
- **GMX_DAA/section1/CP2K/cp2k_P_geoopt.pm3.inp**


**Highlighted regions of cp2k_R_geoopt.inp**

We have to set up calculation type as `GEO_OPT` in the `&GLOBAL` section:

```javascript
&GLOBAL
  PROJECT DA.R        ! Name of the calculation
  PRINT_LEVEL LOW     ! Verbosity of the output
  RUN_TYPE GEO_OPT    ! Calculation type: Geometry optimisation
&END GLOBAL
```

In the `&FORCE_EVAL` section, we define the basci system definitions (such as topology and coordinates). Here we are only going to explain the most important ones:
- `METHOD QS` : Quickstep is the Quantum mechanical method in CP2K 
- In the `&SUBSYS` subsection, we define several parameters of the system. 
  - `&CELL` subsection defines the box size that contains the QM atoms.
  - `&COORD` subsection locates the coordinates of the reactant state. 
  - `&KIND` subsections lists the basis set and the potential for each element in the simulation. 
- In the `&DFT` subsection, we define several parameters of the DFT basis set
  - `&CHARGE 0` sets the charge of the system.
  - `&QS` subsection parameters. We specify here that we are going to use `METHOD PM3`
  - `&MGRID` subsection defines all the important parameters for the GPW.
  - `&SCF` subsection parameters.

```javascript
&FORCE_EVAL                              ! parameters needed to calculate energy and forces
  METHOD QS                              ! QUICKSTEP method
  &SUBSYS                                ! a subsystem: coordinates, topology, molecules and cell
    &CELL                                ! Input parameters needed to set up the simulation cell
      ABC 12.4138 12.4138 12.4138
    &END CELL
    &COORD                               ! Coordinates for simple systems specified using explicit XYZ coordinates
  ... R coordinates ...
    &END COORD
  &END SUBSYS  
  &DFT                                   ! Parameter needed by LCAO DFT programs
    CHARGE 0                             ! The total charge of the system
    &QS                                  ! parameters needed to set up the Quickstep framework
      METHOD PM3                         ! Electronic structure method 
      &SE                                ! Parameters needed to set up the Semi-empirical methods
         &COULOMB                        ! Evaluation of the COULOMB term in SE calculations
           CUTOFF [angstrom] 10.0
         &END
         &EXCHANGE                       ! Evaluation of the EXCHANGE and core Hamiltonian terms in SE calculations
           CUTOFF [angstrom] 10.0
         &END
      &END
    &END QS
    &MGRID                              ! multigrid information
      CUTOFF 200
      NGRIDS 4
      REL_CUTOFF 30
    &END MGRID
    &SCF                                ! Parameters needed to perform an SCF run.
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-05
      MAX_SCF 2000
      &OT                               ! orbital transformation (OT) method
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END
      &OUTER_SCF                        ! parameters controlling the outer SCF loop
        EPS_SCF 1.0E-6
        MAX_SCF 100
      &END
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
&END FORCE_EVAL
```

In the `&MOTION` section, we define the parameters for the geometry optimisation. 

```javascript
&MOTION                        
  &GEO_OPT                    ! Environment of the geometry optimizer
    TYPE MINIMIZATION         ! Kind of geometry optimization
      MAX_ITER 4000
    OPTIMIZER CG
    &CG                       ! Conjugate gradient optimization
      MAX_STEEP_STEPS  0
      RESTART_LIMIT 9.0E-01
    &END CG
  &END GEO_OPT
&END MOTION
```

Then run the following command: 

```javascript
$ cp2k.popt cp2k_R_geoopt.inp > cp2k_R_geoopt.out
```

Once the job has finished, we need to check that the minimisation has converged. For each optimisation step, we have this section: 

```javascript
  Core-core repulsion energy [eV]:                         111210.31809538803645
  Core Hamiltonian energy [eV]:                             -6070.07521304926922
  Two-electron integral energy [eV]:                      -218128.32380280233338
  Electronic energy [eV]:                                 -115134.23711445042863

  Total energy [eV]:                                        -3923.91901906240992

  Atomic reference energy [eV]:                              3926.78210430983836
  Heat of formation [kcal/mol]:                                66.02431952621455

  outer SCF iter =    6 RMS gradient =   0.72E-06 energy =       -144.2013768850
  outer SCF loop converged in   6 iterations or    8 steps


 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):             -144.201376885183862


 *******************************************************************************
 ***                 BRENT   - NUMBER OF ENERGY EVALUATIONS :       1        ***
 *******************************************************************************

 --------  Informations at step =   197 ------------
  Optimization Method        =                   SD
  Total Energy               =      -144.2013768852
  Real energy change         =        -0.0000000005
  Decrease in energy         =                   NO
  Used time                  =                0.857

  Convergence check :
  Max. step size             =         0.0000000000
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                  YES
  RMS step size              =         0.0000000000
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                  YES
  Max. gradient              =         0.0079886843
  Conv. limit for gradients  =         0.0004500000
  Conv. for gradients        =                   NO
  RMS gradient               =         0.0015348400
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. for gradients        =                   NO
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


> **Tip:** Use this command to monitorise the runs: grep -a12 "onvergence check" cp2k_R_geoopt.pm3.out

<br/><br/>

### 1.2 Geometry optimisation and Frequency calulations of the transition states

#### Geometry optimisation of the Transition State

Subsequently, we are going to perform a QM geometry optimisation for the transition state. We are going to use the semi-empirical PM3 again to have consistent results. The CP2K input file is similar to the geometry optimisation for the reactant and product states, but there are several differences:

You will find the input file here: **GMX_DAA/section1/CP2K/cp2k_TS_geoopt.pm3.inp**

**Highlighted regions of cp2k_TS_geoopt.inp**

To calculate the goemetry for the TS, we need to slightly modify the `&GEO_OPT` in the `&MOTION` section. To use the [DIMER method](https://aip.scitation.org/doi/10.1063/1.480097) to calculate the transition states in CP2K. 

```javascript
&MOTION
  &GEO_OPT                          ! Environment of the geometry optimizer
    TYPE TRANSITION_STATE           ! Kind of geometry optimization
    MAX_ITER 2000
    OPTIMIZER CG
    &CG                             ! Conjugate gradient optimization
      MAX_STEEP_STEPS 1000
    &END
    &TRANSITION_STATE               ! Transition state search
      METHOD DIMER                  ! Method to use for locating transition states
      &DIMER                        ! parameters for Dimer Method
        DR 0.0001
        ANGLE_TOLERANCE [deg] 0.5
        INTERPOLATE_GRADIENT
        &ROT_OPT
          OPTIMIZER CG
          MAX_ITER 1000
          &CG
            MAX_STEEP_STEPS 1000
          &END CG
        &END
      &END
    &END
  &END GEO_OPT
&END MOTION
```

<br/><br/>

#### Frequency calculation

Additionally, we need to an extra calculation to validate the TS structure. We need to make sure that there are negative frequencies on the atomic vibrations and that they correspond to the formation of the bonds. For the vibrational analysis to be correct, we need to use the same QM method used in the geometry optimisation used for the geometry optimisation. 

You will find the input file here: **GMX_DAA/section1/CP2K/cp2k_TS_freq.pm3.inp**

**Highlighted regions of cp2k_TS_freq.pm3.inp**

We have to set up calculation type as `NORMAL_MODES` in the `&GLOBAL` section:

```javascript
&GLOBAL
  PROJECT DA.TS.freq       ! Name of the calculation
  PRINT_LEVEL MEDIUM       ! Verbosity of the output
  RUN_TYPE NORMAL_MODES    ! Calculation type: Normal modes
&END GLOBAL
```

We maintain the same `&FORCE_EVAL` section than in the geometry optimisations. We omit the `&MOTION` section and add to following one:

```javascript
&VIBRATIONAL_ANALYSIS     ! Normal Modes, vibrational, or phonon analysis
  NPROC_REP 1             ! number of processors to be used per replica environment
  DX 0.01                 ! increment to construct the HESSIAN with finite difference method 
  FULLY_PERIODIC          ! Avoids to clean rotations from the Hessian matrix
  &PRINT
    &MOLDEN_VIB
    &END
    &CARTESIAN_EIGS
    &END
    &PROGRAM_RUN_INFO
      &EACH
        REPLICA_EVAL 1
      &END
    &END
  &END
&END
```

This calculation creates a series of files: 
- **TS-VIBRATIONS-1.mol** : Resulting frequencies and structures. 
- **TS-VIBRATIONS-1.eig** : binary file containing the eigenvectors and eigenstates of the frequencies.
- **TS-Hessian-1.hess** : binary file containing the Hessian matrix.

In the **TS-VIBRATIONS-1.mol** file, we can see that the top listed frequecies are negative. Ideally, we would have only one negative frequency that corresponds to the formation of the two C-C bonds in a synchonous way. 

```javascript
 [Molden Format]
 [FREQ]
     -739.565631
      -63.251123
      -36.147818
       15.614584
       30.341246  
       ...
```

The easiest way to validate the atomic vibrations is to visualise them. The default output format in CP2K is [MOLDEN](http://cheminf.cmbi.ru.nl/molden/) format (**TS-VIBRATIONS-1.mol**), you might find this software difficult to install or to use, we have provided a python script that converts this format to a multiple XYZ file for each vibration. You can visualise these multiple XYZ files using Pymol of VMD. 

This image shows the first vibration found using CP2K:

![DA TS vibrations](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section1/CP2K/vibration.png)

<br/><br/>

### 1.3 Calculation of the energy barrier.

Here we compare the single-point energy calculations between Gaussian and CP2K:

Energy / kcalmol-1 | Reactants | Transition States | Product
------------ | ------------- | ------------- | -------------
Structure | ![DA R](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section1/DA.reactants.png) | ![DA TS](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section1/DA.ts.png) | ![DA P](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section1/DA.product.png)
Gaussian09 | 0.00 | 47.81 | -8.69
CP2K  | 0.00 | 36.48 | -42.02

![Single Point Energy profile](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section1/CP2K/energy_profile.png)

**Issues found**

Since the results obtained with CP2K do not agree with the results obtained Gaussian09, we tried to fine-tune different paramaters on the CP2K input file to achieve reproducible results. Here are the most relevant findings:

- The optimiser has a very significant effect on the structure. 
  - CG gets a similar structure than Gaussian. 
  - BFGS leads to the structure after SO2 expulsion. This method cannot be used with the DIMER method: ```BFGS method not yet working with DIMER```
  
- SCF convergence: Using MAX_SCF 2000 is useful to converge the SCF in these simulations. 
- GRID CUTOFF convergence doesn't affect the accuracy of these small systems:
```javascript
# Grid cutoff vs total energy
# REL_CUTOFF = 40
# Cutoff (Ry) | Total Energy (Ha)
    150.00  -3921.4603694477
    200.00  -3921.4603694477
    250.00  -3921.4603694477
    300.00  -3921.4603694477
    350.00  -3921.4603694477
    400.00  -3921.4603694477
    450.00  -3921.4603694477
    500.00  -3921.4603694477
    550.00  -3921.4603694477
    600.00  -3921.4603694477
```


<br/><br/>

---

## Section 2: Diels-Alder reaction in vacuum

To obtain the reaction profile of the reaction in vaccum, we are going to use the [Nundged Elastic Band (NEB)](https://theory.cm.utexas.edu/henkelman/pubs/jonsson98_385.pdf) to estimate the energy profile of the Diels Alder reaction. 

We are going to use the resulting structures from the previous section as guesses for the reaction path. Since the structures were optimised using the semi-empirical PM3. 

You will find the input file here: **GMX_DAA/section2/CP2K/cp2k_neb.pm3.inp**

**Highlighted regions of cp2k_neb.pm3.inp**

We have to set up calculation type as `BAND` in the `&GLOBAL` section:

```javascript
&GLOBAL
  PROJECT DAA_NEB
  RUN_TYPE BAND          ! Band methods
  PRINT_LEVEL MEDIUM
&END GLOBAL
```

We use a similar `&FORCE_EVAL` section of the input file because we want to maintain the same QM level. However, we define a new `&MOTION` section. 

- `BAND_TYPE CI-NEB`
- `NUMBER_OF_REPLICA` 10
- `K_SPRING` 0.05
- `&REPLICA` subsection. Here we list the structures as initial guesses for NEB.

```javascript
&MOTION
  &BAND                    ! BAND run
    BAND_TYPE CI-NEB       ! Type of BAND calculation
    NUMBER_OF_REPLICA 15   ! Number of Replica to use in the BAND
    K_SPRING 0.05          ! Spring constant
    
    &CONVERGENCE_CONTROL.  ! control the convergence criteria for BAND
      MAX_FORCE 0.0010
      RMS_FORCE 0.0050
    &END
    
    ROTATE_FRAMES TRUE.    ! Compute at each BAND step the RMSD and rotate the frames in order to minimize it
    ALIGN_FRAMES TRUE      ! Alignment of the frames at the beginning of a BAND calculation
    
    &CI_NEB                ! CI-NEB type calculation only.
      NSTEPS_IT 2
    &END
    
    &OPTIMIZE_BAND         ! optimization method for the band 
      OPT_TYPE DIIS        ! DIIS based optimization procedure for BAND
      OPTIMIZE_END_POINTS TRUE
      &DIIS
        MAX_STEPS 1000
      &END
      
    &END
    &PROGRAM_RUN_INFO
    &END
    &CONVERGENCE_INFO
    &END

    &REPLICA               ! Coordinates and velocities (possibly) of the replica
      COORD_FILE_NAME ./DA.R-pos-1.xyz
    &END
    &REPLICA
      COORD_FILE_NAME ./DA.TS-pos-1.xyz
    &END
    &REPLICA
      COORD_FILE_NAME ./DA.P-pos-1.xyz
    &END
  &END BAND

&END MOTION
```

> **TIP**: The order of the atoms MUST be the same on all the structures. 

The energy profile obtained is the following: 

Nodes  |  Energy / kcalmol-1
------------ | ------------- 
R  | 0.00
1  | 0.00 
2  | -0.89 
3  | -1.59 
4  | -1.44 
5  | -1.43 
6  | -0.77 
7  | -0.02 
8  | 0.20 
9  | -0.90 
10  | 0.16 
11  | -56.15 
12  | -97.13 
13  | -97.92 
P  | -97.92 

![NEB_15replicas](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section2/CP2K/NEB_15rep.png)

**Issues found:**

These results do not match with the energies found in the **Section 1**. Therefore, we tried to use with different paramaters on the input file. Here are the most relevant findings:

- The final structure of the product is not the structure we found using either Gaussian nor CP2K (using CG optimiser). The product found has expulsed SO2 after the Diels Alder reaction. 
- Using ```OPTIMIZE_END_POINTS TRUE```, CP2K also optimises the end points and generates a smoother Energy profile.


<br/><br/>

---

## Section 3: Diels-Alder reaction in solution

The next step is to simulate the Diels Alder reaction in solution. To do so we are going to use a QM/MM approach where part of the molecule is going to be defined as the QM region and the rest of the organic molecule and the solvent are goint to be defined as MM. This setup will require us to define link atoms, we will explain how to do so in this section.

This section is divided into 4 parts:
- Setup of the system. 
- Classical minimisation and equilibration of the system.
- Monitorisation of the QM/MM set up.
- Biased simulations runs. 

<br/><br/>

### Section 3.1: Setup of the system.

**Ligand parameterisation**

We are only going to use a simple approach, we will parameterise the reaction product and model the reaction backwards ( **P --> R** ). To do so, we will use the ambertools protocol (using antechamber and prmchk2). For more details on the process: 
http://ambermd.org/tutorials/basic/tutorial4b/

We are going to use the full molecule this time:

![Full product molecule](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/full.product.png)

```javascript
$ antechamber -i product.pdb -fi pdb -o product.mol2 -fo mol2 -nc -1 -c bcc -at gaff2 -rn DAA
$ parmchk2 -i product.mol2 -f mol2 -o product.frcmod
```

Using this protocol, we are going to use GAFF2 forcefield parameters to define atomtupes and use AM1-BCC QM level to calculate the charges. 

<br/><br/>

**System preparation**

Using tleap from the ambertools suite, we are going to solvate the system in a waterbox, including a Na+ counterion to neutralise the system.

The tleap input file (**leap.in**) looks like this:
```javascript
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams product.frcmod

LIG = loadmol2 product.mol2

solvatebox LIG TIP3PBOX 15
addions LIG Na+ 1
savepdb LIG system.pdb
saveamberparm LIG system.prmtop system.inpcrd

quit
```

Finally we execute tleap like this:

```javascript
$ tleap -f leap.in
```

> TIP: Always visualise the system to check everything is OK.

<br/><br/>

### Section 3.2: Classical minimisation and equilibration of the system.

Here we are going to leverage the AMBER tools and AMBER tutorial [B0](http://ambermd.org/tutorials/basic/tutorial0/). Using AMBER I used a similar protocol to the one in AMBER forcefield to minimise & equilibrate the system MM. We are not going to explain in detail the AMBER input file, you will find a detailed description in the aforementioned tutorial. 

You will find the input file here: 
- **GMX_DAA/section3/CP2K/classical_equilibration/min_classical.in**
- **GMX_DAA/section3/CP2K/classical_equilibration/heat_classical.in**

**min_classical.in**

```javascript
Initial minimisation of our structure
 &cntrl
  imin=1, maxcyc=4500, ncyc=2000,
  cut=8.0, ntb=1, ntc=2, ntf=2
 /
```

**heat_classical.in**

```javascript
Heat
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim=10000, dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0, temp0=300.0,
  ntpr=100, ntwx=100,
  cut=8.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1, ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
```

To run the minimisation, you have to use this command:

```javascript
$ sander -O -i min_classical.in -o min_classical.out -p system.prmtop -c system.inpcrd -r system.min.r
$ sander -O -i heat_classical.in -o heat_classical.out -p system.prmtop -c system.min.r -r system.md.r -x system.md.nc
```

<br/><br/>

### Section 3.3: Monitorisation of the QM/MM set up.

There are several steps in order to set up a QM/MM system:
- Modifications of the Lennard-Jones parameters, 
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

**Modification of Lennard-Jones parameters**

We need to be aware of the limitations of the classical models and how to combine them with the QM methods. We have set up the solvent using TIP3P water model, which has no Lennard-Jones parameters defined for the hydrogen atoms. We have to set them in order to avoid unrealistic interactions with the QM region. Here in this example, we are going to use the LJ parameters in the GAFF2 forcefield.

To modify them, we are going to modify the prmtop file using parmed. Here are the PARMED commands to run:

```javascript
$ parmed system.prmtop
```

```javascript
changeLJSingleType :WAT@H1 0.3019 0.047
changeLJSingleType :WAT@H2 0.3019 0.047
outparm system.LJ_mod.prmtop
quit
```

<br/><br/>

**Fix the atomic coordinates format obtained form the equilibration step**

In order to use the restart file from the heat simulations, we need to change the format from binary netcdf to a six columns formatted restart file. CPPTRAJ (from AMBERtools) can do that for us and simultaneosuly remove periodic boundary conditions artifacts.

```javascript
$ cpptraj system.prmtop
```

```javascript
trajin system.md.r
autoimage
trajout system.md.crd restrt
go
quit
```

> **TIP**: CRD AMBER format (10 columns) and RST7 AMBER (6 columns) format differ on the number of columns. We are converting the file to restrt (RST7) format but naming the file as .crd as CP2K requires. 

For a more detailed explanation of the cpptraj commands, you can look at the [CPPTRAJ manual](https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml). 

<br/><br/>

**Definition of the QM/MM boundaries**

After that we set up the QM/MM partition as follows: 

![QM/MM partition](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section3.png)

QM region is highlighted in blue and link atoms are shown in purple. The rest of the system (water molecules and counterions) is shown in white. 

As shown in the scheme above, the QM/MM boundary divides the ligand and hence, covalent bonds. To deal with broken covalent bonds, we are going to use the [**IMOMM** scheme](https://www.semanticscholar.org/paper/IMOMM%3A-A-New-Integrated-Ab-Initio-%2B-Molecular-of-Maseras-Morokuma/5e1bb154312e6751c044c8226a9ee73b69db89f6). This scheme adds an additional atomic centre (usually a hydrogen atom) covalently bonded to the last QM atom. This additional atom is not real and its only purpose is to saturate its valency replacing the bond that has been broken. 

To set the link atoms in CP2K, we have to do the following steps.

- Delete the atoms that you want to exclude from the QM part.
- Identify the link atoms.
- Add the link atoms to the system under the QM kind part: 

```javascript
   &LINK
      MM_INDEX  5
      QM_INDEX  38
      LINK_TYPE IMOMM
    &END LINK
```

Here it is an image showing all the indexes of the atoms for this solvated system:

![QMMM indexes](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section3_indexes.png)

<br/><br/>

You will find the input file here: **GMX_DAA/section3/CP2K/monitorisation_qmmm/qmmm_md.link_atoms.inp**

We have to set up calculation type as `MD` in the `&GLOBAL` section:

```javascript
&GLOBAL 
  PROJECT QMMM
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL
```

The ```&FORCE_EVAL``` section retains part of the definitions we used in the vacuum QM calculations. We have to a lot of new parameters:
- ```METHOD QMMM``` : We will use method QMMM
- Similar ```&DFT``` section
- ```&MM``` section: Here we specify all the necessary options to run the MM region of the system. 
  - ```&FORCEFIELD``` section
  - ```&POISSON``` section
- ```&SUBSYS``` section: List all the details of the system:
  - ```&CELL``` section: Size of the system. 
  - ```&TOPOLOGY``` section: Here we specify the AMBER input files. 
  - ```&KIND NA+``` section to specify all the elements that are defined in the system with a different element name (Na -- NA+)
- ```&QMMM``` section: 
  - ```ECOUPL COULOMB``` electrostatic coupling scheme
  - ```&CELL``` section: size of the QM subsystem. 
  - ```&QM_KIND``` section: We list all the indexes of each element.  
  - ``` &LINK``` section to list the link atom indexes.

```javascript
&FORCE_EVAL
  METHOD QMMM
  STRESS_TENSOR ANALYTICAL
  &DFT
    CHARGE 0
    &QS
      METHOD PM3
    &END QS
    &SCF
      MAX_SCF 2000
      EPS_SCF 1.0E-6
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 100
      &END
    &END SCF
  &END DFT
    &MM
    &FORCEFIELD
      PARMTYPE AMBER
      PARM_FILE_NAME post.LJ_mod.prmtop
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
    !Set box dimensions here
      ABC [angstrom] 49.9302020  41.8911420  42.8511700
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &TOPOLOGY
      CONN_FILE_FORMAT AMBER
      CONN_FILE_NAME post.LJ_mod.prmtop
      COORD_FILE_FORMAT CRD
      COORD_FILE_NAME post.md.crd
    &END TOPOLOGY
    &KIND NA+
     ELEMENT Na
    &END KIND
  &END SUBSYS
    &QMMM
    ECOUPL COULOMB
    &CELL
      ABC 40 40 40
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &QM_KIND N
      MM_INDEX 33
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 16 17 35 37
    &END QM_KIND
    &QM_KIND C
      MM_INDEX 27 28 29 30 31 32 34 36 38
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 19 20 21 22
    &END QM_KIND
    &QM_KIND CL
      MM_INDEX 23 24 25 26
    &END QM_KIND
    &QM_KIND S
      MM_INDEX 18
    &END QM_KIND
    &LINK
      MM_INDEX  5
      QM_INDEX  38
      LINK_TYPE IMOMM
    &END LINK
  &END QMMM
&END FORCE_EVAL
```

The main changes are located in the ```&MOTION``` section:

- ```&MD subsection``` : Define all the parameter to run a MD. Including:
  - ```ENSEMBLE``` : NPT_I, NVT...
  - ```&BAROSTAT``` subsection
  - ```&THERMOSTAT``` subsection
  - ```&PRINT``` subsection

```javascript
&MOTION
  &MD
  ENSEMBLE NPT_I
  TIMESTEP [fs] 0.5
  STEPS    20000
  TEMPERATURE 298
  &BAROSTAT
    TIMECON [fs] 100
    PRESSURE [bar] 1.0
  &END BAROSTAT
  &THERMOSTAT
    REGION GLOBAL
    TYPE CSVR
    &CSVR
      TIMECON [fs] 10.
    &END CSVR
  &END THERMOSTAT
  &END MD
  &PRINT
    &RESTART                                    ! This section controls the printing of restart files
      &EACH                                     ! A restart file will be printed every 10000 md steps
        MD 25000
      &END
    &END
    &TRAJECTORY                                 ! Thes section Controls the output of the trajectory
      FORMAT DCD                                ! Format of the output trajectory is DCD
      &EACH                                     ! New trajectory frame will be printed each 100 md steps
        MD 4
      &END
    &END
    &RESTART_HISTORY                            ! This section controls dumping of unique restart files during the run keeping all of them.Most useful if recovery is needed at a later point.
      &EACH                                     ! A new restart file will be printed every 10000 md steps
        MD 5000
      &END
    &END
    &CELL
      &EACH
        MD 100
      &END
    &END
  &END PRINT
&END MOTION
```

<br/><br/>

### Section 3.4: QM/MM enhanced sampling (Metadynamics)

As we saw in the previous QM/MM simulations, the system remains estable and runs smoothly. However, for the Diels Alder reaction to be reversed we should extend our equilibrium simulations for a very long time and probably nothing would happen given the energy difference between reactants and products. Therefore, we need to biase the system in order to sample the Diels Alder reaction. To biase the system we need to define collective variables that describe the reaction or change we want to observe and use enhanced sampling methods to push the system onto that coordinate. 

CP2K has different [collective variables (CV)](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html) and several [enhanced sampling methods](https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/FREE_ENERGY.html) implemented. However, for this particular problem none of the implemented CV fulfil our needs. Luckily, we can used the [PLUMED2 plugin](https://www.plumed.org) to define a custom CV. To do so, you will need an additional file (plumed.dat) describing your CV to CP2K.

**plumed.dat**:

``` javascript
COM ATOMS=30,31 LABEL=com1
COM ATOMS=29,32 LABEL=com2

DISTANCE ATOMS=com1,com2 LABEL=d1
UPPER_WALLS ARG=d1 AT=6.0 KAPPA=150.0 EXP=2 EPS=1 OFFSET=0 LABEL=uwall

#
# Activate metadynamics in d1
# depositing a Gaussian every 500 time steps,
# with height equal to 1.2 kJ/mol,
# and width 0.35 rad for both CVs.
#
METAD ARG=d1 PACE=500 HEIGHT=1.2 SIGMA=0.35,0.35 FILE=HILLS LABEL=metad

# monitor the d1 distance, the upper wall and the metadynamics bias potential
PRINT STRIDE=10 ARG=d1,uwall,metad.bias FILE=COLVAR
```

**qmmm_md.link_atoms.metadynamics.inp**

In the ```&GLOBAL``` section goes as: 

```javascript
&GLOBAL
  PROJECT METAD
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL
```

In the ```&MOTION``` section add: 

```javascript
&FREE_ENERGY
  &METADYN
    USE_PLUMED .TRUE.
    PLUMED_INPUT_FILE  ./plumed.dat
  &END METADYN
&END FREE_ENERGY
  ```

You can also run metadynamics using the collective variables defined in [CP2K](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html). Here I add a sample version of the ```&FREE_ENERGY``` section:

```javascript
  &FREE_ENERGY
    METHOD METADYN
    &METADYN
      DO_HILLS  .TRUE.
      NT_HILLS 100
      WW [kcalmol] 1.5
      &METAVAR
        WIDTH 0.5 !Also known as scale
        COLVAR 1
      &END METAVAR
    &PRINT
        &COLVAR
           COMMON_ITERATION_LEVELS 3
           &EACH
             MD 10
           &END
        &END
      &END
    &END METADYN
  &END FREE_ENERGY
```

And the COLVAR definition is added to the ```&MM``` section of the ```&FORCE_EVAL```:

```javascript
    &COLVAR
       &DISTANCE_FUNCTION
          COEFFICIENT 1
          ATOMS 29 30 32 31
       &END DISTANCE_FUNCTION
    &END COLVAR
```

<br/><br/>

---

## Section 4: Diels-Alder reaction in protein

The protocol used here is very similar to the protocol described for the QM/MM system in solution. The steps are the in order to set up a QM/MM system:

- Modifications of the Lennard-Jones parameters,
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

### Section 4.1: Setup of the system. 

**Cleaning the initial structure (from PDB database)**

We start with the catalytic antibody structure (PDB id 1C1E) and we are goingto build the system shown below:

![Protein-ligand system](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section4.png)

First, you split the system into 3 parts (protein, ligand, waters). See **fix_1c1e.sh**:

```javascript
# Split water molecules, protein and ligand ENH
grep "HOH" 1c1e.pdb > 1c1e.xray_waters.pdb
grep -v "HOH" 1c1e.pdb | grep -v "ENH" > 1c1e.protein.pdb
grep "ENH" 1c1e.pdb > 1c1e.ligand.pdb
```

Secondly, you use pymol to align the Diels Alder product to the ligand ENH and save the new coordinates (Result: **DAA.aligned.pdb**). 

<br/><br/>

**Ligand parameterisation**

Again, we will use the ambertools protocol (using antechamber and prmchk2). For more details on the process: 
http://ambermd.org/tutorials/basic/tutorial4b/

```javascript
$ antechamber -i DAA.aligned.pdb -fi pdb -o product.mol2 -fo mol2 -nc -1 -c bcc -at gaff2 -rn DAA
$ parmchk2 -i product.mol2 -f mol2 -o product.frcmod
```

Using this protocol, we are going to use GAFF2 forcefield parameters to define atomtupes and use AM1-BCC QM level to calculate the charges. 

<br/><br/>

**System preparation**

Using tleap for the ambertools suite, we are going to prepare the system (which inclused the product solvated in a waterbox with a counterion to neutralise the system. 

The tleap input file (**leap.in**) looks like this:

```javascript
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams product.frcmod

PROT = loadpdb 1c1e.protein.pdb
LIG = loadmol2 product.mol2
AGUAS = loadpdb 1c1e.xray_waters.pdb

SYS = combine {PROT LIG AGUAS}

solvatebox SYS TIP3PBOX 15
addions SYS Cl- 5
savepdb SYS system.pdb
saveamberparm SYS system.prmtop system.inpcrd

quit
```

Finally we execute tleap like this:

```javascript
$ tleap -f leap.in
```

> TIP: Always visualise the system to check everything is OK.

<br/><br/>


### Section 4.2: Classical minimisation and equilibration of the system.

Here we are going to lto use the same protocol described in the **Section 3.2**. Here are the input files:

- GMX_DAA/section4/CP2K/classical_equilibration/min_classical.in
- GMX_DAA/section4/CP2K/classical_equilibration/heat_classical.in

**min_classical.in**

```javascript
Initial minimisation of our structure
 &cntrl
  imin=1, maxcyc=4500, ncyc=2000,
  cut=8.0, ntb=1, ntc=2, ntf=2
 /
 ```
 
**heat_classical.in**

```javascript
Heat
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim=10000, dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0, temp0=300.0,
  ntpr=100, ntwx=100,
  cut=8.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1, ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
```

To run the minimisation, you have to use this command:
```javascript
$ sander -O -i min_classical.in -o min_classical.out -p system.prmtop -c system.inpcrd -r system.min.r
$ sander -O -i heat_classical.in -o heat_classical.out -p system.prmtop -c system.min.r -r system.md.r -x system.md.nc
```

<br/><br/>

### Section 4.3: Monitorisation of the QM/MM set up.

There are several steps in order to set up a QM/MM system:
- Modifications of the Lennard-Jones parameters, 
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

**Modification of Lennard-Jones parameters**

Again we need to modify the ennard-Jones parameters of TIP3P water model and the hydrogen atom of the hydroxyl groups of the Ser and Tyr residues of the protein. We are going to use the LJ parameters in the GAFF2 forcefield again.

To modify them, we are going to modify the prmtop file using parmed. Here are the PARMED commands to run:

```javascript
$ parmed system.prmtop
```

```javascript
changeLJSingleType :TYR@HH 0.3019 0.047
changeLJSingleType :SER@HG 0.3019 0.047
changeLJSingleType :WAT@H1 0.3019 0.047
changeLJSingleType :WAT@H2 0.3019 0.047
outparm system.LJ_mod.prmtop
quit
```

<br/><br/>

**Fix the atomic coordinates format obtained form the equilibration step**

In order to use the output of the heat simulations, we need to convert the structure from netcdf (binary) to CRD (AMBER) format. We are going to use the CPPTRAJ tool ( provided by the AMBERtools free suite ) to conver the format and recenter the water box. 

```javascript
$ cpptraj system.prmtop
```

```javascript
trajin system.md.r
autoimage
trajout system.md.crd restrt
go
quit
```

For a more detailed explanation of the cpptraj commands, you can look at the [CPPTRAJ manual](https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml). 

<br/><br/>

**Definition of the QM/MM boundaries**

We are going to use the same QM MM partition used in **section 3** and hence, we have to set up the appropiate link atoms. You need to check the atom indexes for each system. 

```javascript
   &LINK
      MM_INDEX  6616
      QM_INDEX  6649
      LINK_TYPE IMOMM
    &END LINK
```

<br/><br/>

The CP2K input resembles the one used in the **Section 3.3**. The index numbers have changed and will need to be updated: **GMX_DAA/section4/CP2K/monitorisation_qmmm/qmmm_md.link_atoms.inp**

```javascript
    &QMMM
    ECOUPL COULOMB
    &CELL
      ABC 40 40 40
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &QM_KIND N
      MM_INDEX 6644
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 6648 6646 6628 6627
    &END QM_KIND
    &QM_KIND C
      MM_INDEX 6638 6639 6640 6641 6642 6643 6649 
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 6630 6631 6632 6633
    &END QM_KIND
    &QM_KIND CL
      MM_INDEX 6634 6635 6636 6637 
    &END QM_KIND
    &QM_KIND S
      MM_INDEX 6629
    &END QM_KIND
    &LINK
      MM_INDEX 6616
      QM_INDEX  6649
      LINK_TYPE IMOMM
    &END LINK
  &END QMMM
```

<br/><br/>


### Section 4.4: QM/MM enhanced sampling (Metadynamics)

As we saw in the previous QM/MM simulations, the system remains estable and runs smoothly. However, we need to set up biased simulations for the Diels Alder reaction to happen. We are going to use the same protocol described before in **section 3.4**

**plumed.dat**:

```javascript
COM ATOMS=6640,6643 LABEL=com1
COM ATOMS=6641,6642 LABEL=com2

DISTANCE ATOMS=com1,com2 LABEL=d1
UPPER_WALLS ARG=d1 AT=6.0 KAPPA=150.0 EXP=2 EPS=1 OFFSET=0 LABEL=uwall

#
# Activate metadynamics in d1
# depositing a Gaussian every 500 time steps,
# with height equal to 1.2 kJ/mol,
# and width 0.35 rad for both CVs.
#
METAD ARG=d1 PACE=500 HEIGHT=1.2 SIGMA=0.35,0.35 FILE=HILLS LABEL=metad

# monitor the d1 distance, the upper wall and the metadynamics bias potential
PRINT STRIDE=10 ARG=d1,uwall,metad.bias FILE=COLVAR
```

**qmmm_md.link_atoms.metadynamics.inp**

In the ```&GLOBAL``` section goes as: 

```javascript
&GLOBAL
  PROJECT METAD
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL
```

In the ```&MOTION``` section add: 

```javascript
&FREE_ENERGY
  &METADYN
    USE_PLUMED .TRUE.
    PLUMED_INPUT_FILE  ./plumed.dat
  &END METADYN
&END FREE_ENERGY
  ```

<br/><br/>
