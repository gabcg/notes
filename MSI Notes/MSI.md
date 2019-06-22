# Molecular Simulations

## Table of Contents

* [Part 1\. Classical Unbiased MD](#part-1-classical-unbiased-md)
* [Molecular Simulations](#molecular-simulations-1)
  * [1\. Introduction to Molecular Dynamics Simulations](#1-introduction-to-molecular-dynamics-simulations)
    * [1\. Introductory Concepts](#1-introductory-concepts)
      * [1\.1 Timescales](#11-timescales)
      * [1\.2 Protein Motion](#12-protein-motion)
    * [2\. Molecular Mechanics](#2-molecular-mechanics)
      * [2\.1 Types of Energy](#21-types-of-energy)
        * [Stretching Energy](#stretching-energy)
        * [Bending Energy](#bending-energy)
        * [Torsion Energy](#torsion-energy)
        * [Non\-Bonded Energy](#non-bonded-energy)
        * [Summary](#summary)
    * [3\. Molecular Simulations Algorithms and Definitions](#3-molecular-simulations-algorithms-and-definitions)
      * [3\.1 Types of Force Field](#31-types-of-force-field)
      * [3\.2 Scheme of the Algorithm](#32-scheme-of-the-algorithm)
      * [3\.3 Integration Algorithms](#33-integration-algorithms)
      * [3\.4 Simulation Environment and Water Models](#34-simulation-environment-and-water-models)
        * [Solvent and Water Models](#solvent-and-water-models)
      * [3\.4 Periodic Boundary Conditions](#34-periodic-boundary-conditions)
      * [3\.6 Short and Long Range Interactions](#36-short-and-long-range-interactions)
      * [3\.7 Temperature](#37-temperature)
      * [3\.8 Pressure and Chemical Potential](#38-pressure-and-chemical-potential)
      * [3\.9 Ensembles](#39-ensembles)
    * [4\. Typical MD Simulation](#4-typical-md-simulation)
      * [4\.1 Steps](#41-steps)
      * [4\.2 Programs](#42-programs)
      * [4\.3 Limitations](#43-limitations)
  * [2\. VMD](#2-vmd)
    * [Visualization and Basic Analysis of Simulations](#visualization-and-basic-analysis-of-simulations)
      * [Representations](#representations)
      * [Scripting with Tcl](#scripting-with-tcl)
    * [Visualization of MDs](#visualization-of-mds)
  * [Class 3: GROMACS](#class-3-gromacs)
    * [NOTES TAKEN IN CLASS](#notes-taken-in-class)
    * [Introduction](#introduction)
    * [Preparation of the Structure](#preparation-of-the-structure)
    * [Generation of the Topology File](#generation-of-the-topology-file)
    * [Box Definition and Solvation](#box-definition-and-solvation)
    * [Minimization](#minimization)
    * [Equilibration](#equilibration)
    * [Production Run](#production-run)
    * [Analysis](#analysis)
    * [NOTES BY REPEATING THE TUTORIAL](#notes-by-repeating-the-tutorial)
      * [Generate Topology](#generate-topology)
      * [Define box and Solvate](#define-box-and-solvate)
      * [Add ions](#add-ions)
      * [Energy Minimization](#energy-minimization)
      * [Equilibration](#equilibration-1)
      * [Production MD](#production-md)
      * [Analysis](#analysis-1)
  * [Class 4\. GROMACS II](#class-4-gromacs-ii)
    * [Minimization Protocol](#minimization-protocol)
      * [Integrator](#integrator)
    * [NVT Equilibration Protocol](#nvt-equilibration-protocol)
      * [Integrator](#integrator-1)
      * [Output Control](#output-control)
      * [Bond Parameters: Constraints](#bond-parameters-constraints)
      * [Nonbonded Settings](#nonbonded-settings)
      * [Electrostatics](#electrostatics)
      * [Temperature Coupling](#temperature-coupling)
      * [Periodic Boundary Conditions](#periodic-boundary-conditions)
      * [Velocities](#velocities)
      * [Analysis of the Simulation in VMD](#analysis-of-the-simulation-in-vmd)
  * [Class 5\. Application of MD to a Membrane Protein](#class-5-application-of-md-to-a-membrane-protein)
    * [Simulation Setup](#simulation-setup)
      * [Step 0: Structure Preparation](#step-0-structure-preparation)
    * [Simulation Setup](#simulation-setup-1)
    * [Equilibration and Production](#equilibration-and-production)
    * [Analysis](#analysis-2)
  * [Class 6: Markov\-Based Analysis of Biomolecular Systems](#class-6-markov-based-analysis-of-biomolecular-systems)
    * [Introduction](#introduction-1)
* [Part 2\. Beyond Classical MD](#part-2-beyond-classical-md)
  * [C8\. Enhanced Sampling Techniques](#c8-enhanced-sampling-techniques)
    * [Introduction](#introduction)
    * [CV\-Dependent Methods](#cv-dependent-methods)
      * [Metadynamics](#metadynamics)
      * [Umbrella Sampling](#umbrella-sampling)
    * [CV\-Free Methods: Replica Exchange](#cv-free-methods-replica-exchange)
    * [PLUMED](#plumed)
  * [C10\. Umbrella Sampling Simulations](#c10-umbrella-sampling-simulations)
  * [C11\. Interactive and Steered Molecular Dynamics](#c11-interactive-and-steered-molecular-dynamics)
    * [IMD](#imd)
    * [SMD](#smd)

## General

**Professors**

* Coordinator: [Jana Selent](mailto:jana.selent@upf.edu)

**Some Links and Useful Information**

* [gmx_tutorials](https://github.com/jalemkul/gmx_tutorials_livecoms/blob/master/gmx_tutorials.pdf): it consists on the tutorials we use for GROMACS, but in a paper-like PDF.
* [Tutorial on VMD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2972669/figure/F2/): there is also a PDF version, but all the figures appear at the end, which is inconvenient.
* [Salt Bridges Plugin, Version 1.1](https://www.ks.uiuc.edu/Research/vmd/plugins/saltbr/)
* [VMD Script Library](http://www.ks.uiuc.edu/Research/vmd/script_library/)

* To view the boundaries of a box, open VMD tcl console and type `pbc box`.

Part 1. Classical Unbiased MD
-----------------------------

## 1. Introduction to Molecular Dynamics Simulations

> Date: 02/04/2019

### 1. Introductory Concepts

Molecular simulations can be considered as a virtual microscope with high temporal and spacial resolution, that complements conventional experiments.

They allow studying complex, dynamic processes occuring in biological systems, such as:

* Protein stability.
* Conformational changes.
* Protein folding.
* Molecular recognition: proteins, DNA, membranes, complexes.
* Ion transport in biological systems.

#### 1.1 Timescales

Different processes occur in different time scales:

* Bond vibration: femtoseconds ($10^{-15}$ s).
* Bond rotation: picoseconds ($10^{-12}$ s).
* Ion crossing: nanoseconds ($10^{-19}$ s).

Depending on the size of the system or the computing power, different timescales can be achieved. For example, small systems can be simulated on a milisecond ($10^{-3}$ s) scale, but bigger systems would need a milisecond for this.

To arrive to the scale of seconds, we can wait for technological advances or use advanced algorithms.

#### 1.2 Protein Motion

In proteins, bond vibration happens in picoseconds to nanoseconds, while an helix movement can happen in miliseconds to seconds.

We can consider a protein conformation as a snapshot in a single point of time, but proteins are dynamic and change from one conformation to other. Some of those cluster in zones of similar conformations, called **conformational states**. Typically, we can consider active, inactive and intermediate conformational states.

<img src="msi-notes.assets/C1-1_conformational-states.png" alt=""
	title="" width="300"/>

To perform a molecular simulation, it is important to use good experimental data:

* Crystallography allows a good resolution, and shows the average of a given conformational state.
* Spectroscopy (NMR/fluorescence): resolution-wise, small proteins offer a good resolution with this technique, while big proteins do not behave so well. As an advantage, it allows to see different structures (even if they are just changes of an inserted chemical probe).

### 2. Molecular Mechanics

Molecular mechanics is a computational method to calculate geometries and energies.

Systems are simplified in order to perform the calculations with larger molecules. An example of this is that atoms are treated as rubber balls of different size (atom types) joined together by springs of varying length.

<img src="msi-notes.assets/C1-2_atoms.png" alt=""
	title="" width="150"/>


Molecular mechanics enables the calculation of the total steric energy in terms of deviations from reference “unstrained” bond length angles and torsions plus non-bonded interaction (VdW, electrostatic):

$$
E_{tot} = E_{str} + E_{bend} + E_{tor} + E_{elec} + E_{vdb} + ... 
$$

#### 2.1 Types of Energy

##### Stretching Energy

* It consists on a sum of the energy of each bond.
* There are two parameters that are assigned to each pair of bonded atoms based on their type (e.g. C-C, C-H):
	* $k_b$: spring constant, controls the stiffness of the bond: the higher, the stiffer, which means that more energy is required to deform the bond from its equilibrium value.
	* $r_0$: it comes from experiments such as nature, crystallography, ab initio; or from quantum mechanics calculations. It is collected in the force field file.
* High differences between $r$ and $r_0$ increase the energy.

<img src="msi-notes.assets/C1-3_stretching-energy.png" alt=""
	title="" width="400"/>

##### Bending Energy

* The equation estimates the energy associated with the vibration about the equilibrium bond angle.
* In this case, parameters are assigned to triplets of atoms, based on their type.
	* $\theta_0$: reference value for which energy E adopts the optimum.
	* $K_0$: controls the stiffness of the angle spring. A larger value means that more energy is required to deform or bond the angle from its equilibrium value. 

<!--The reference is still in the force field.-->

<img src="msi-notes.assets/C1-4_bending-energy.png" alt=""
	title="" width="400"/>

##### Torsion Energy

* Torsion energy can be considered as how one part rotates in reference to the other.
* It corrects the other energy terms, rather than representing a physical process. It makes the total energy agree with the experiment.
* Parameters are asssigned to quartets of atoms:
	* $A$: controls the amplitude of the curve.
	* $n$: controls its periodicity
	* ɸ: shifts the entire curve along the rotation angle axis $\tau$.
* It is modeled by a simple periodic function, which has this form because rotation can keep happening.

<img src="msi-notes.assets/C1-5_torsion-energy.png" alt=""
	title="" width="400"/>

In the example, it is not very preferred and A would be high.

##### Non-Bonded Energy

* It is the combination of two terms: van der Waals and electrostatic.
* In both terms, $r_{ij}$ represents radius.
* The fact that they are summing means that in long distance they do not feal each other.
* Van der Waals
	* Approximated by a Lennard-Jones 6-12 potential. Itself has two terms, attraction and repulsion:
		* Attraction: it occurs at a short range, and rapidly dies off when the distance between the interacting pair of atoms increases (by just a few armstrongs, as it has a 1/r^6 dependency).
		* Repulsion: occurs when the distance is slightly less than the sum of their contact radii. It is modeled by an equation that  is designed to rapidly blow up at close distances  (1/r^12 dependency).
	* There are two parameters (which are defined in the force field) that control the depth and position (interatomic distance) of the potential energy well for a given pair of atoms:
		* $A$ defines stickinness of the attraction. It has a negative sign.
		* $B$ defines the level of hardness of the atoms (marshmallow-like, billiard ball-like, etc.): sof atoms can be deformed when they are close, hard cannot. It belongs to the repulsion term.
* Electrostatic energy:
	* It is approximated by a Coulomb potential.
	* It is a function of the charge of each atom ($q_i$ and $q_j$), their interatomic distance ($r_{ij}$) and a molecular dielectric expression that accounts for the attenuation of electrostatic interaction by the environment (such as the solvent or the molecule itself).

<img src="msi-notes.assets/C1-6_non-bonded-energy.png" alt=""
	title="" width="400"/>

* Sum: in long distances they dont feel each other.

##### Summary

<img src="msi-notes.assets/C1-7_summary-energy.png" alt=""
	title="" width="400"/>

### 3. Molecular Simulations Algorithms and Definitions

#### 3.1 Types of Force Field

Force fields contain parameters of the energies we have seen and potential energy functions, which enable molecular mechanics calculations.

There are different types of force field, depending on the application. They evolve, as the longer simulation times allow to see problems, which are corrected.

In a CHARMM36 force field, we can observe:

* The mass of the atoms (even for different types of the same element).
* Bond definition.
* Spring constants for two atoms.
* Angle definition and spring constants for angles (three atoms).
* It finishes with comments, which defines the origin of the information contained in the force field.

#### 3.2 Scheme of the Algorithm

* The solvation of $E_{tot}$ is the computationally most expensive step of a simulation.
* A molecular simulation can be seen as the computation of the behaviour of a molecular system as a function of time.
* Newton's laws of motions are solved in a simulation, concretely the second ($F=ma$). This is done iteratively.
* A time step $\Delta t$ needs to be chosen. The larger it is, the faster the simulation.
* It is desirable to choose a small $\Delta t$, but it is problematic.
	* The integration step is limited by the highest frequency motions: the vibration of hydrogen bonds (typically 0.5 fs).
	* Constraining those bonds can speed up the simulations to 2 fs.
	* Biomolecular system with constraints on all bond lengths: 4 to 5 fs.
	* From a doubt:iIn the simulation, you don't really fix H from side chain and the backbone (but you can if you add more restraints).
* After choosing the initial position of atoms and $\Delta t$: get acting force.
* Then, the atom acceletation. 
* The next step is to move the atoms, calculating velocities. Those are also dependent on the temperature (which is present in the formula).
	* *It seems that for some parts, the velocity is assigned randomly* 
* Finally, the time is moved forward and the procedure is repeated.


<img src="msi-notes.assets/C1-8_summary-md-alg.png" alt=""
	title="" width="450"/>

> The number of iterations for 1 μs is = 500 x 10^6. *I guess this is with the 2 fs timestep.* Special purpose computers such as Anton allow to simulate big systems during long times.

#### 3.3 Integration Algorithms

The leap-frog algorithm and the velocity Verlet integrator are present in most molecular dynamics software. Both algorithms produce identical trajectories when used without speciall additional features:

**Leap-frog algorithm**

In this algorithm, the position is called $r$, and is defined at time $t$. The velocity is calculated between each time ($t-1/2\Delta t$).

Positions and velocities are updated using the forces $F(t)$ determined by the positions at time $t$.

<img src="msi-notes.assets/C1-9_leap-frog-1.png" alt=""
	title="" width="250"/>

<img src="msi-notes.assets/C1-10_leap-frog-2.png" alt=""
	title="" width="250"/>
 
**The velocity Verlet integrator**

In this case, positions $r$, velocities $v$ and forces $F(t)$ are determined at time $t$

<img src="msi-notes.assets/C1-11_verlet.png" alt=""
	title="" width="300"/>


*Probably "time t" is actually "time delta t" for both cases*

#### 3.4 Simulation Environment and Water Models

##### Solvent and Water Models

Simulating in vacuum is unrealistic: a solvent is needed. Some approaches put dielectric constants of solvents to solve the molecular system, but its not as good as having water.

Different water models have been created. Parameters are defined to reproduce the properties of water at room temperature and atmospheric pressure. The most widely used are rigid (its vibration is stiff).

###### Transferable Interaction Potential (TIP)

TIP3P: TIP model with 3 interaction sites centred on the atomic nuclei. Positive partial charges on the hydrogen atoms and negative on the oxygen.

TIP4P: TIP model with 4 interaction sites. The negative charge is moved 0.015 nm off the oxygen towards the hydrogens along the bisector of the HOH angle. Slightly better than TIP3P, but computationally more expensive.

###### Simple Point Charge (SPC)

As TIP, it is a series of models. All of them have 3 interaction sites centred on the atomic nuclei. Positive partial charges on the hydrogens and negative on the oxygen.

SPC with C6 Lennard-Jones parameters on the hydrogen atoms. The best model with three sites within the TIP and SPC families.

###### Model Election

The water model can be chosen independently of the biomolecular force field, but often this is not a good idea:

* AMBER, CHARMM, and OPLS protein force fields have been parameterised with TIP3P.
* The GROMOS protein force field has been parameterised with SPC.

#### 3.4 Periodic Boundary Conditions

Periodic boundary conditions are related with the amount of water to use. Adding more water requires volume, increasing the computation time.

One solution for this is the creation of identical virtual unit cells. The atoms in the surrounding virtual system interact with atoms in the real system.

These modeling conditions allow to eliminate the surface interaction of the water molecules, and the representation they create is more similar to the *in vivo* environment that a water sphere surrounded by vacuum.

#### 3.6 Short and Long Range Interactions

The electrostatic forces are separated in two groups: short and long range. They are summed in Fourier space.

<img src="msi-notes.assets/C1-12_sri-1.png" alt=""
	title="" width="250"/>

The particle-mesh-Ewald (PME) methos is the most used to resolve long range electrostatics:

* Cutoff: 12 Å (after this cutoff, calculations are not performed anymore?)
* Switching distance: 10 Å

<img src="msi-notes.assets/C1-13_sri-2.png" alt=""
	title="" width="400"/>
	
#### 3.7 Temperature

Temperature is an experimental condition that needs to be reproduced. It can be assumed as kinetic energy, and control of it is needed as velocity is lost due to cutoffs.

Control of the temperature is done by adding a thermostat. Some types are:

* Berendsen thermostat (weak coupling).
* Nosé-Hoover thermostat (extended-system Hamiltonian).
* Andersen thermostat (stochastic coupling). Not present in all simulation software.
* Langevin dynamics

#### 3.8 Pressure and Chemical Potential

Pressure can be controled by adjusting the cell size (big cells mean lower pressures and viceversa). The following barostats are used:

* Berendsen barostat (weak coupling).
* Parrinello-Rahman barostat and variants (extended-system Hamiltonian).

Regarding the chemical potential, it can be reproduced by adding or removing particles.

#### 3.9 Ensembles

Ensembles keep some parameters fixed. As a summary, the canonical ensemble has a termostat, and the isobaric also a barostat.

**Microcanonical Ensemble (NVE)**

The thermodynamic state characterized by a fixed number of atoms, $N$, a fixed volume $V$ and a fixed energy $E$. This corresponds to an isolated system.

**Canonical Ensemble (NVT)**

Collection of all systems whose thermodynamic state is characterized by a fixed number of atoms, $N$, a fixed volume $V$ and a fixed temperature $T$.

**Isobaric-isothermal Ensemble (NPT)**

Characterized by afixed number of atoms, $N$, a fixed pressure $P$ and a fixed temperature $T$.

### 4. Typical MD Simulation

#### 4.1 Steps 

A molecular dynamics simulation involves a series of steps.

1. Get a protein structure (for example, from PDB).
2. Fix missing segments, side chains and define protonation states. 
	* Some side chains do not appear because they are flexible and are not captured when the structure is being defined.
	* With respect to the protonation state, amino acids can be protonated or not depending on the environment, so it cas to be calculated if some side chains should be protonated or not.
3. Prepare input files: add the solvent, corresponding ions and choose the right parameters.
4. Energy minimization: geometry is optimized here.
5. Equilibration simulation: allows getting the right box size?
6. Run production simulation.
7. Analyze the output trajectory.

#### 4.2 Programs

There are different programs, with advantages and disadvantages. It is recommended to stick to one or two and learn them in detail, as they have a lot of options. The choosen algorithms needs to be justified when publishing.

ACEMD is a commercial program, but it is really fast, as it runs on GPUs. This allows good parallelization. GROMACS is an open source alternative which also has a GPU version it allows molecular dynamics simulations and energy minimization.

#### 4.3 Limitations

* Parameter definition for the simulation is not perfect, so it can always be further refined.
* Limited polarization effects; waters can reorient but partial charges are fixed.
* The phase space is not sampled exhaustively.

Also, seeing a single in a simulation does not mean it is relevant: when it is observed more than once it can be considered significant.

## 2. VMD

> Date: 04/04/2019

### Visualization and Basic Analysis of Simulations

Open `1UBQT.pdb` the structure with `File > New Molecule...`

The `Mouse` menu allows to change between modes, but there are also keys:

* `R` for rotation mode.
* `T` for translation mode.
* `S` for scale mode (zoom in/out).

Other useful keys:

* `1` to click on atoms.
* `2` to click on bonds.
* `3` to click on angles.
* `4` to click on dihedrals (torsion angles, between 4 atoms ⚛️).
* `=` to re-center view.

#### Representations

The `Graphics > Representation` menu controls everything about the representation, with options to manipulate it. It is divided in various tabs:

* Draw settings: includes the coloring method and the drawing method. Drawing methods such as 'Lines' allow to change thickness.

Some drawing methods:

* VDW: atoms are represented as spheres. Allows changing the Sphere Scale and Sphere Resolution
* Tube: allows to see the shape of the backbone. *You can change its side*
* NewCartoon: highlights second structure traits (*useful to count them*).

Coloring methods: *S is in yellow*

* Name: each atom has its own color
* Residue type: self-descriptive.

<!--Append the drawing methods image-->

> In our example we observe a beta sheet (formed by 4 beta strands. The arrow represents the C-terminal. A beta sheet can be parallel or antiparallel depending on the direction of the strand. As this has both, it is a mixed beta sheet.
> 
> The small helix is a little bit tighter than a normal alpha helix. <!--Add the image in slide 11 about structure of ubiquitin-->

Other tab is called **Selections**. It allows to select certain elements, to create representations about them or to hide what we are not interested in (*for example, we can have all protein in the line drawing method, but serines in VDW*). This relies on **keywords** (*i.e. selecting alpha helix and then index, shows the atom numbers, starting by 0*).

```reformat
> PDB Files. Index 1 of pdb in vmd shows the dirst atom. Keywords in vmd are equivalent to pdb files.
>
>Selet vdw and in selection type index 1: shows first atom. Then index 1 to 10
>
>also resid 1-10 (first strand of beta sheet if we put new cartoon again.
>
>in name N in VDW, we see
>
>res name in selection, keyword resname shows the residues there. It doesnt change when i say resname LYS THR but it should show them
>
>if you type
>
>use x-y-z coordinates:
>
> x>10 armstrong
```

In selections, connectors `and`, `or` and `not` can be used to define more complex selections. Examples:

* `betasheet and name N` to select nitrogens in the beta sheets. or to see both.
* `not betasheet and not helix` to select just loops.
* `water within 3 of protein` water molecules within 3 armstrongs of the protein. The bigger the armstrongs, the more water molecules we see. <!--Doubt, guess i need a specific drawing method or something to see them-->

Finally, particular representations can be saved. When we have multiple ones, we can deactivate some (to focus) or see all of them at once (to highlight traits in the context of the structure).

When we finish creating representations of a protein and we have it as we want, we can save this information for another session with the menu `File -> Save State` or create an image with `File -> Render`

<!--Tasks: put the why globulins have shape, sequence viewer extension-->

#### Scripting with Tcl

The Tcl scripting language can be used to automate actions in VMD. It is a straightforward way to work. It can be accessed through `Extensions > Console`.

Some basics:

* Printing: `puts "Hello World!"`.
* Arithmetic operations: `expr -3 * 10`.
* Variables: `set x [expr -3 * 10]`.
* Printing variables: `puts $x`.

We can select and measure atoms:

```
set sel [atomselect top "all"]
```

Top is the  molecule that is defined as "top" in the molecule list.

When a selection is done, we can get information from it:

```
$sel num
```

Commands from slide:

```
set sel [atomselect top "all"]
$sel num
$sel moveby {10 0 0}
set sel [atomselect top "hydrophobic"] $sel get resname
set sel [atomselect top "hydrophobic and name CA"] $sel get resname
$sel get {resname resid}
measure center $sel
measure minmax $sel
```

If we don't know a particular command, we can activate `File > Log Tcl Commands to Console` and then do what we want graphically. We will see in the console the equivalent command.

### Visualization of MDs

Dynamics of drug receptor interactions: the delta opioid receptor bound to naltrindole.

First we need to open the following file `structure_WT.pdb` as new molecule. Once we've done it, select `Load Data into Molecule` and choose the `.dcd` file (which contains the frames). To finish, load the visualization state (`.vmd` extension).

Perform RMSD alignment: `Extensions > Analisys > RMSD Trajectory Tool`

Steps in slide

Plot saved as rmsd_plot.ps

Natrindol rings

3: 1
5: 2 <!--Ask about the strange oxygen-->
6: 4 (2 are aromatic

so we need to do a selection of ALA LEU VAL ILE PRO PHE MET TRP that are at least 3 armstrongs

((resname ALA or resname LEU or resname VAL or resname ILE or resname PRO or resname PHE or resname MET or resname TRP) as within 3 of resname EJ4) and not hydrogen

((resname ALA or resname LEU or resname VAL or resname ILE or resname PRO or resname PHE or resname MET or resname TRP) as within 3 of resname EJ4) and not hydrogen

resname ALA as within 3 of resname EJ4 and not hydrogen

## 3: GROMACS

> Date 09/04/2019

### NOTES TAKEN IN CLASS

### Introduction

The goal of today is to learn to setup a simulation with GROMACS. An overview of the pipeline that the program follows is:

<!--Insert image from PDF-->

1. Checking the PDB file for missing parts. Certain parts of the structure are too flexible to be captured by X-rays. This usully happens with the surface of the protein, as it is interacting with the solvent.
2. Add hydrogens to the structure and generate the topology with `pdb2gmx`.
3. Build a simulation box with `editconf` and fill it with a solvent with `solvate`.
4. Add ions with `genion` to neutralize possible net charges.
5. Use `grompp` to combine different inuts into a binary `.tpr` file.
6. Run mdrun to simulate.
7. Analyse the simulation.

### Preparation of the Structure

Most of this part is already "done", as the structure we use does not have missing parts. The only thing we need to do remove the water molecules from the PDB file with the text editor of our choice.

### Generation of the Topology File

We look into the topology file. It has info ab the molecule. You get 4 each residue the atom tyoe, resnumber atom itself, charge, mass, sum of the charge.

Some extra things it says. We are going to do equilibration by constraining heavy atoms (those that are not hydrogen).

At the bottom there is a list of molecules: we will put more molecules in the next steps there.

### Box Definition and Solvation

This is done with `editconf` and `solvate`. We will use a cubic box, but using a rhombic dodecahedral box would be a better option, as it optimizes the volume. This accelerates computation because there are less solvent molecules.

> We make the periodic boundary conditions in this boxes. *I don't know where this comes from, i think it has to do with the size of the box*. Even in assymetric proteins there is always 1 nm of distance between the protein and the box wall.

Once we solvate the protein, we can observe in the topology file that a certain amount of solvent has been added to the molecules list.

The water module we use (spc model) is used as an standard by GROMACS, but there are other models.

Neutralize the system. Reasons:

* In nature, net charge is neutral, so our system also does it. But how much do we need to add? We can see in the topology file a sum of the charges: its net charge.

This step is not so straightforward. First, it requires to create a `.tpr` file (as in the simulations). It is a binary file that compiles the information of the input. This files are always generated using `grompp`.

The input needed for adding ions is the solvation box and the topology file, and also a `.mdp` file (which is provided by the tutorial).

Once the binary file is generated, run `genion`. It generates another `.gro` file, edits the topology file (to add the ions to the molecules list).

We add Na or Cl to neutralize (depends on the net charge. It is also needed to select which type of molecule will be exchanged (solvent molecules in this case, as we don't want it to exchange protein atoms). This changes are reflected both in the `.gro` and topology files.

> It seems that the output of `genion` already tells how many ions were added, so it's not really needed to open the topology file.
> 
> Anyway, we can observe the `.gro` file in VMD to check if the ions were added: open it and create a representation with a selection `ions` and draw them as VDW.

In nature net charges are zero. Maybe they concentrate somewhere, or move, but net is zero.

------

### Minimization

Solvated neutralized system. Reason: we added estimated hydrogens (maybe they need to be minim and also find their right orientation within the system). We also added water and ions.

At the beginning we have to find some positional constraits, which is only heavy tom. During this minim basically you minimize in the 1st step the H atoms. For that you need to use grompp to make a binary file with info from topology, coordinates and the mdp protocol. 

> This protocols can be obtained from publications when we want to repeat what was published, but if not they are kind of dependent on the force field you are using. Most of the values are default values, but you can do whatever you want if you can justify it. The more you imulate, you learn more what he parameters should be. Starting with cutoffs of how many minimizations, this is basically something that is learnt by doing. It is also possible to optimize protocols for types of proteins (i.e. proteases) and use this protocol for all of them. *And I guess you use it as a template, modifying it a bit for specific proteins within the type*.

> Minimization is important because if you start from high energies, the simularion algorithm explodes.

Some stuff contained in the minimization protocol:

* Integrator: steepest descent minimization.
	> When minimizing you don't explore the entire energetic landscape, but go to the next local minimum. The global might be in other place. The algorithm always looks that the structure is going down, it will never try to go up. 
* When to stop the minimization, based in what happens first:
	* A certain number of iterations (50000 in this case).
	* When the maximum force is below a certain value.
* Cutoffs for non-bonded interactions (electrostatic, VDW), divided in:
	* Short interactions: treated with the potential energy function, which is 10 A.
	* Long interctions are treated with the (...) equation, an approximation to solve and converge quickly on the long range using Fourier space.
* A list of pairs. Always calc bw 2 particles. To speed up calculations there is a step where it checks, makes a list of pairs around a certain atom and checks, puts (...) and then you only need to calculate from this list of pairs the non-bonded interactions. This is updated each timestep. *<I don't understand this part*
* Cutoffs scheme: a neighbor searching? just means to the normal cutoff. Srves as a kind of buffer that maybe a atom moves slightly from the cutodd so its buffered a little bit. *neither this one*

**Results of the minimization**

When we run the EM, the outputs are:

* em.log: ASCII-text log file of the EM process
* em.edr: Binary energy file
* em.trr: Binary full-precision trajectory
* em.gro: Energy-minimized structure

We can use `.edr` files to analyse the whole minimization (as it contains the energy of every step). We first convert it to a non-binary format with `energy`. The command always prompts what you want to compute. In this case, we type `10`.

The resulting file can be plotted with our plotting tool of choice.

> With GNUplot, the fist lines of the file need to be removed, ad the command is `plot '<filename>.xvg' with line`

We can see an animation of the minimization by opening the molecule (the last `.gro`), and then opening the `.trr` file with "Load into molecule".

> She says she is disapponted as the `.trr` only contained one frame, so it' not a real animation. We could have configured it to contain more frames.

### Equilibration

The equilibration has two steps. For both of them, a `.mdp` protocol and the generation of a binary `.tpr` is needed.

The structure input is the minimized one, with the topology. The protocol nvt.mdp: constant volume and temperature. We want to equilibrate this ensemble.

* We try not to use position restraints on the protein. That means water and ions can move and also actually the side chains of the protein, but not the backbone.
* Integrator: leap-frog intgrator.
* 50000 integrations, one integration step we use 2 fs. 2*50000 = we simulate foe 100 ps.
* Output control: how often we save a frame into the trajectory file. Step here refers to fs, so each one picosecond it writes a frame.
* Same thing for velocities and energy, and the log file.
* Bond-parameters:
	* First, say if the simulation is a continuation of other or not (in this case it's not).
	* Then, (from theory class): we can apply a timestep of 2 fs because we constrain all H bonds, which is the fastest vibration (it is below 2 fs).
		* When simulating one itration, you have to do it capturing the fastest movement in the systems (H vibration 0.5 fs). We can do the iteration faster by constraining H bonds, so they don't vibrate and are always kept at the same distance. Therefore we define this. Taking off the constraint would mean that we woud go to a lower timestamp and see the bond vibrations. 
* Nonbonded-param
* Cutoff: 1 nm
* The buffer
* List of atoms within the 1 nm cutoff but it has a buffers to 1.2 nm (so if during the simulation it moves, we still consider it as within 1 nm).
* Electrostatic
* The PME: we have to get an approx for long-range because if not it is computationally too demanding, so there are different things you can use but the one ther eis the most default one.
* A aram  for cubic box interpolation. Something about Fourier and a grid with 0.6 nm (at each point you calculate the electrostatics).
* A thermostat to control the temperature.
* A reference temperature (300 K, usual for running simulations).
* As we don't use NPT (yet) there is not barostat.
* Periodic boundary conditions
* Velocity generations: yes, because we start from the beginning, so we want to start the molecule according to the kinetic energy and the temperature (...) Maxwell distribution.

Afte doing this equilibration, we take the result, use `energy` and plot `temperature.xvg`. We want to see if the system has a constant temperature or not.

Now we **generate the NPT ensemble to continue the simulation**, by adding a barostat. There two main changes in the `.mdp` protocol. If we look it:

* It defines that we continue from a previous simulation (the NTP).
* It has pressure coupling on (it was off before.

> The input to do the plot can be found in `3.5_NPT`.

### Production Run

Generate the `.tpr` using the proper protocol, the structure obtained from NPT and the topology file.

### Analysis

For the analysis, the simulation needs to be aligned and centered on the protein. For that, you can check if the protein is stable or not by calculating the RMSD. You choose the backbone atoms. Make a plot of the RMSD. Check if the protein changes conformation or converges to a more or less stable conformation.

> To do this from the output she gives, go to the folder.

Usually RMSD changes: increases a little bit and maybe it changes again, and identify different conformations of the protein. We observe how the protein tries to converge to something.

We can also do this with VMD:

* Load `nvt.gro` (from a `3.5` folder).
* Then load inside the molecule `md_0_1.xtc` (from `3.7`). *This is the trajectory of last step. Now it is an `.xtc`, a more compressed simulation file*.
* Go to the RMSD Trajectory tool, click "Align" and then click "RMSD".

Random forgotten notes:

* take md_0_1.gro from 3.6
* class 3.7 xtc

### NOTES BY REPEATING THE TUTORIAL

#### Generate Topology

**Step One: Prepare the Topology**

* Hen egg white lysozyme (PDB code 1AKI), from [PDB](http://www.rcsb.org/pdb/home/home.do). Some preprocessing:
	* Delete waters:
		```
		grep -v HOH 1aki.pdb > 1AKI_clean.pdb
		```
		> This should not be done with all molecules: sometimes water is functional.
	* Check entries with the `MISSING` tag (atoms or whole residues not present). Absent terminal regions are not a problem. Model the molecule if needed, if not, `pdb2gmx` will fail. This module only genereates topologies.
* Next we execute `pdb2gmx`:
	* Command: `gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce`
		* It asks for a force field to process the molecule, in this case an all-atom OPLS force field (so type `15`). 
	* Generated files:
		* The topology for the molecule: `topol.top`
			* Contains all the information necessary to define the molecule within a simulation:
				* nonbonded parameters (atom types and charges)
				* bonded parameters (bonds, angles, and dihedrals) 
		* A position restraint file: `posre.itp`
			* defines a force constant used to keep atoms in place during equilibration  
		* A post-processed structure file: `1AKI_processed.gro`.
			* contains all the atoms defined within the force field (i.e., H atoms have been added to the amino acids in the protein)

There are many other options that can be passed to pdb2gmx. Some commonly used ones are listed here:

* -ignh: Ignore H atoms in the PDB file; especially useful for NMR structures. Otherwise, if H atoms are present, they must be in the named exactly how the force fields in GROMACS expect them to be. Different conventions exist, so dealing with H atoms can occasionally be a headache! If you need to preserve the initial H coordinates, but renaming is required, then the Linux sed command is your friend.
* -ter: Interactively assign charge states for N- and C-termini.
* -inter: Interactively assign charge states for Glu, Asp, Lys, Arg, and His; choose which Cys are involved in disulfide bonds.

**Step Two: Examine the Topology**

We can examinate the topology file with any text editor.

* It starts with comments (preceded by `;`).
* Then it has a line calling the force field file.
* After this, we find the protein name (Protein_A as the PDB had a chain A), ith 3 exclusions for bonden neighbors (read in manual).
* Then, there is a section that defines the atoms with several columns.
* Other sections are [ bonds ], [ pairs ], [ angles ], and [ dihedrals ]. <!--The parameters and function types associated with these sections are elaborated on in Chapter 5 of the GROMACS manual-->
* Then there is a call for the position restraints, which are in a file also created by pdb2gmx.
* Remainder: defining other molecules (solvent, ions) and providing system-level descriptions (system that will be printed in the output and a list of all molecules in the system).

> The type of water is SPC/E, as we passed `-water spce`. Other typical choices for water include SPC, TIP3P, and TIP4P.
> 
> We observe in the file that water is position-restrained, using a force constant (k_pr) of 1000 kJ mol-1 nm-2.

--

> key notes about the [ molecules ] directive:

> 1. The order of the listed molecules must exactly match the order of the molecules in the coordinate (in this case, .gro) file.
> 2. The names listed must match the [ moleculetype ] name for each species, not residue names or anything else.
> 
> If you fail to satisfy these concrete requirements at any time, you will get fatal errors from grompp (discussed later) about mismatched names, molecules not being found, or a number of others.

#### Define box and Solvate

**Step Three: Defining the Unit Cell & Adding Solvent**

* We will simulate in an aqueus system (though other solvents can be used).
* For this we define a box and fill it with solvent. This has two steps.

Define the box dimensions using the editconf module.

We will use a simple cubic box (but a rhombic dodecahedron is recommended, as its volume is ~71% of the cubic box of the same periodic distance, thus saving on the number of water molecules). The command is:

```
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```

In this command we say:

* Center the protein in the box (-c)
* Place it at least 1.0 nm from the box edge (-d 1.0)
* The box type is a cube (-bt cubic)

> Explanation of the distance:
> 
> Since we will be using periodic boundary conditions, we must satisfy the minimum image convention. That is, a protein should never see its periodic image, otherwise the forces calculated will be spurious. Specifying a solute-box distance of 1.0 nm will mean that there are at least 2.0 nm between any two periodic images of a protein. This distance will be sufficient for just about any cutoff scheme commonly used in simulations.
 
Next, we fill the box with water:

```
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```

Explanation of the command:

* The configuration of the protein (-cp) is contained in the output of the previous editconf step
* The configuration of the solvent (-cs) is part of the standard GROMACS installation.
	* We are using spc216.gro, which is a generic equilibrated 3-point solvent model. You can use spc216.gro as the solvent configuration for SPC, SPC/E, or TIP3P water, since they are all three-point water models.
* Output: 1AKI_solv.gro
* We tell the name of the topology file (topol.top) to modify it (it writes the number of waters in the [ molecules ] directive, thi behaviour doesnt happen with non-water solvents).

#### Add ions

* `pdb2gmx`, based on the amino acid composition told us that the protein has a net charge of +8e (last line of your [ atoms ] directive in topol.top).
* Ions should be added to the system with the `genion` tool.
	* It reads through the topology replaces water molecules with the ions that the user specifies.
	* Its input is called run input file (.tpr), produced by the `grompp` module. It is also used when running the simulation.
		* `grompp` processes the coordinate file and topology (which describes the molecules) to generate an atomic-level input .tpr. This file contains all the parameters for all atoms in the system.

To produce the .tpr. a .mdp (molecular dynamics parameter file) is needed. grompp assemples the arameters of this file with the coordinates and topology into the tpr file.

> mdp file is normally used to run energy minimization or an MD simulation, but in this case is simply used to generate an atomic description of the system.

-- 

> *This file is fownloaded from the tutorial page*
> 
> In reality, the .mdp file used at this step can contain any legitimate combination of parameters. I typically use an energy-minimization script, because they are very basic and do not involve any complicated parameter combinations. Please note that the files provided with this tutorial are intended only for use with the OPLS-AA force field. Settings, particularly nonbonded interaction settings, will be different for other force fields.

Building the .tpr file is done with:

gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

hen we pass the file to genion. It will prompt us, choose group 12 SOL (so it doesnt replace parts of the protein with ions):

gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

The options of this command are:


* The input is a structure/state file (-s) as input.
* The output is a a .gro file (-o).
* Process the topology (-p) to reflect the removal of water molecules and addition of ions
* Define positive and negative ion names (-pname and -nname, respectively).
* Tell genion to add only the ions necessary to neutralize the net charge on the protein by adding the correct number of negative ions (-neutral, which in this case will add 8 Cl- ions to offset the +8 charge on the protein).
	* You can also use genion to add a specified concentration of ions in addition to simply neutralizing the system by specifying the -neutral and -conc options in conjunction. Refer to the genion man page for information on how to use these options.

About ion names: they used to be dependent on the specific force field file, but now they are standard. They are specified in capital letters. This molecule type is then written to the topology file:

```
[ molecules ]
; Compound      #mols
Protein_A         1
SOL           10636
CL                8
```

> Residue or atom names may or may not append the sign of the charge (+/-), depending on the force field. Do not use atom or residue names in the genion command, or you will encounter errors in subsequent steps.

#### Energy Minimization

**Step Five: Energy Minimization**

* Ensure that the system has no steric clashes or inappropriate geometry. The structure is relaxed through a process called energy minimization (EM).
* The process is similar to adding ions: grompp generates a tpr, but instead of using genion we use mdrun (GROMACS MD engine).
* Again, the input parameter file from gromp is given by the tutorial.

```
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```

We proceed with the energy minimization:

```
gmx mdrun -v -deffnm em
```

* We run the command in verbose because it can be slow.
* The -deffnm flag will define the file names of the input and output. 
	* So, if you did not name your grompp output "em.tpr," you will have to explicitly specify its name with the mdrun -s flag. 

mdrun generates the followng files:

* em.log: ASCII-text log file of the EM process
* em.edr: Binary energy file
* em.trr: Binary full-precision trajectory
* em.gro: Energy-minimized structure

And this, even if we don't use verbose:

```
Steepest Descents converged to Fmax < 1000 in 831 steps
Potential Energy  = -5.8702450e+05
Maximum force     =  9.0797430e+02 on atom 736
Norm of force     =  2.0987490e+01
```

We can **evaluate** if the EM was successful by:

* Potential energy should be negative, and (for a simple protein in water) on the order of 10-5-10-6, depending on the system size and number of water molecules.
* The maximum force had a targe, specified at the *minim.mdp - "emtol = 1000.0"*
	* indicating a target Fmax of no greater than 1000 kJ mol-1 nm-1. It is possible to arrive at a reasonable Epot with Fmax > emtol. If this happens, your system may not be stable enough for simulation. Evaluate why it may be happening, and perhaps change your minimization parameters (integrator, emstep, etc).
* Also, the file em.edr contains all energy terms collected by GROMACS during EM. We can analyse it with the energy module:

	```
	gmx energy -f em.edr -o potential.xvg
	```
	> At the prompt, type "10 0" to select Potential (10); zero (0) terminates input
	
	* The result of this command is an average Epot and a file called potential xvg, which can be `ñptted. It will show a steady convergence of the Epot. can be plotted, and we'll see the average Epot 
	
		```
		Statistics over 831 steps [ 0.0000 through 830.0000 ps ], 1 data sets
		All statistics are over 658 points (frames)
		
		Energy                      Average   Err.Est.       RMSD  Tot-Drift
		-------------------------------------------------------------------------------
		Potential                   -564593      11000    28730.6     -75544  (kJ/mol)
		```
		<!--Insert the plot-->
		
#### Equilibration

Step Six: Equilibration

* Now, we have a reasonable starting structure in terms of geometry and solvent orientation.
* Now we should equilibrate the solvent and ions around the protein.
	* If not (unrestrained dynamics) the system may collapse, because the solvent is optimized within itself, and not necesarily with the solute.
	* We need to bring the solvent to the temperature we want to simulate. *The correct temperature is based on kinetic energies*.
	* Then, establish the proper orientation about the solute (the protein).
	* At last, we apply pressure to the system until it reaches the proper density.

We will use the posre.itp file that pdb2gmx generated:

* It applies position restraining forces on heavy atoms of the protein (non-H atoms).
	* This allows movement, but only after overcoming a substantial energy penalty.
	* Position restraints allow to equilibrate the solvent around the protein.
	* The origin of the position restraints (the coordinates at which the restraint potential is zero) is provided via a coordinate file passed to the -r option of grompp.

Equilibration is usually conducted in two phases:

* Phase 1: NVT ensemble (constant Number of particles, Volume, and Temperature), or "isothermal-isochoric" or "canonical".
	* The timeframe is dependent upon the contents of the system:
		* The temperature of the system should reach a plateau at the desired value.
		* If the temperature has not estabilized, then more time is required.
		* Typically, 100-150 ps should be enough.
	* A .mdp file is required for grompp to run (downloaded from tutorial). Some parameters are:
		* `gen_vel = yes`: Initiates velocity generation. Using different random seeds (`gen_seed`) gives different initial velocities, and thus multiple (different) simulations can be conducted from the same starting structure.
		* `tcoupl = V-rescale`: The velocity rescaling thermostat is an improvement upon the Berendsen weak coupling method, which did not reproduce a correct kinetic ensemble.
		* `pcoupl = no`: Pressure coupling is not applied.
* Phse 2: NPT or "isothermal-isobaric" ensemble, which more closely resembles experimental conditions.
	* This stabilizes the pressure (thus also the density) of the system.
	* Again, a mdp file is needed (and downloaded). It is similar to the nvt.mdp, but has some changes:
		* It has a pressure coupling section which uses the Parrinello-Rahman barostat. 
		* `continuation = yes`: We are continuing the simulation from the NVT equilibration phase
		* `gen_vel = no`: Velocities are read from the trajectory (see below).
	* Again, `grompp` is called and then we run the simulation and evaluate energies. 

The commands for the equilibration are:

```
mx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -deffnm nvt
```

Then we can anlyse the temperature progression:

```
gmx energy -f nvt.edr -o temperature.xvg
```

> Type `16 0` at the prompt to select the temperature of the system and exit.

In the pressure equilibration:

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -deffnm npt
```

> We are now including the -t flag to include the checkpoint file from the NVT equilibration; this file contains all the necessary state variables to continue our simulation. To conserve the velocities produced during NVT, we must include this file.
> 
> The coordinate file (-c) is the final output of the NVT simulation.

```
gmx energy -f npt.edr -o pressure.xvg
```

> Type `18 0` at the prompt to select the pressure of the system and exit.

The obtainded plots for temperaure are:

From the plot, it is clear that the temperature of the system quickly reaches the target value (300 K), and remains stable over the remainder of the equilibration. For this system, a shorter equilibration period (on the order of 50 ps) may have been adequate.

About the plots for pressure:

The pressure value fluctuates widely over the course of the 100-ps equilibration phase, but this behavior is not unexpected. The running average of these data are plotted as the red line in the plot. Over the course of the equilibration, the average value of the pressure is 7.5 ± 160.5 bar. Note that the reference temperature was set to 1 bar, so is this outcome acceptable? Pressure is a quantity that fluctuates widely over the course of an MD simulation, as is clear from the large root-mean-square fluctuation (160.5 bar), so statistically speaking, one cannot distinguish a difference between the obtained average (7.5 ± 160.5 bar) and the target/reference value (1 bar).

We also take a look at density:

```
gmx energy -f npt.edr -o density.xvg
```

> entering "24 0" at the prompt.

As with the pressure, the running average of the density is also plotted in red. The average value over the course of 100 ps is 1019 ± 3 kg m-3, close to the experimental value of 1000 kg m-3 and the expected density of the SPC/E model of 1008 kg m-3. The parameters for the SPC/E water model closely replicate experimental values for water. The density values are very stable over time, indicating that the system is well-equilibrated now with respect to pressure and density.

Please note: I frequently get questions about why density values obtained do not match my results. Pressure-related terms are slow to converge, and thus you may have to run NPT equilibration slightly longer than is specified here.

#### Production MD

**Step Eight: Production MD**

* The system is now well-equilibrated at the desired temperature and pressure.
* We will release the position restraints and run production MD for data collection.
* We will use of the checkpoint file which contains preserve pressure coupling information. This is used by grompp with a script called `md.mdp`:
	
	```
	gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
	```
	> `grompp` will print an estimate of how many ocessors should be dedicated to the PME calculation:
	>
	> * For a cubic box, the optimal setup will have a PME load of 0.25 (3:1 PP:PME - we're very close to optimal!)
	> * For a dodecahedral box, the optimal PME load is 0.33 (2:1 PP:PME)
* Then, we run the simulation. It is a 1-ns imulation (guess the script specifies it).

	```
	gmx mdrun -deffnm md_0_1
	```

#### Analysis

**Step Nine: Analysis**

We use some of the GROMACS tools in order to analyse the simulation:

* `trjconv`: a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc). **We'll use tis "corrected" trajectory for further analysis.**
	* In this exercise we use it to account for any periodicity in the system. The protein will diffuse through the unit cell, and may appear "broken" or may "jump" across to the other side of the box. To account for such actions, issue the following:
		
		```
		gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
		```
		> Select 1 ("Protein") as the group to be centered and 0 ("System") for output
* With the trajectory obtained from last step, we can calculate the RMSD:
	
	```
	gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
	```
	> Choose 4 ("Backbone") for both the least-squares fit and the group for RMSD calculation
	
	* The -tu flag will output the results in terms of ns, even though the trajectory was written in ps.
		* This is done for clarity of the output (especially if you have a long simulation - 1e+05 ps does not look as nice as 100 ns).
	* The output plot will show the RMSD relative to the structure present in the minimized, equilibrated system.
* We can also calculate the RMSD relative to the crystal structure:

	```
	gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
	```
	
	* Both time series show the RMSD levels off to ~0.1 nm (1 Å), indicating that the structure is very stable. Subtle differences between the plots indicate that the structure at t = 0 ns is slightly different from this crystal structure. This is to be expected, since it has been energy-minimized, and because the position restraints are not 100% perfect, as discussed previously.
* With `gyrate` we can measure the radious of fyration (a measure of its compactness).
	
	```
	gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
	```
	> Choose group 1 (Protein) for analysis.
	
	* Two things can happen:
		* Maintaina a relatively stable value for R_g: it is stably folded.
		* R_g changes over time: the protein unfolds.


## Class 4. GROMACS II

<!--Date: 11/04/2019-->

First we focus a bit on some details, and then we analyse a simulation.

Min questions were about the simulation protocols, so here we go.

### Minimization Protocol

It is contained in the file `minim.mdp`.

#### Integrator

The integrator is the algorithm that performs the minimization. They look through the energetic landscape by different iterations. In each of them, they calculate the potential energy by summing all the energies of the system:

$$
E_{tot}=E_{str}+E_{bend}+E_{tor}+E_{vdw}+E_{ele}+...
$$

The steepest descent minimizer is a first order method, which gradually changes the coordinates of atoms as moving the system closer to the minimum.

It looks for local minimums. It starts in a random point, so it might fall in the global minimum or in a local minimum.

This algorithm can be explained with a 2D simplification, in which what it does is similar to a line search:

* Draw a line and pick 3 points on it randomly (as we don't know where the minimum for this line is).
* Calculate the potential energy function for each point.
* Get the point with less energy. Decide in which segment there could be a smaller energy. *In the picture E2 is the point with less energy, and that probably between it and E3 there is a smaller energy*.
* Put three points in the selected segment and calculate the energies.

The segment selection is repeated until the lowest point is found. When this happens, a perpendicular line that crosses the first through this point is drawn, and the three point strategy is done again.

Ths strategy is fast, but is forced to make right-angled turns at each point even though it might not be the best  route to the minimum. This acceptable in narrow/steep minimums, as it can go directly to the bottom.

> Gives undesirable properties in long narrow valley <!--Doubt-->

> Each iteration depicted by arrows in slide 11.

For shallow minimums, the conjugate gradient minimization is more appropiate. It is a bit more computationally expensive, as it stores data from previous iterations (history gradient information: performs calculations following a gradient). This information is used to choose a direction, which is used to go more directly for the local minimum.

### NVT Equilibration Protocol

It is contained in the `nvt.mdp` file.

#### Integrator

The integrator in this case uses a leap-frog algorithm. During the cycle of each iteration of the simulation, we have to obtain both the coordinates and the velocities of the atoms. This algorithm calculates each at different times:

* The coordinates are calculated at each time-step (in the tutorial, 2 fs).
* The velocities are calculated in the half-time.

This is done because the velocities can be very different, so the average of the velocities is taken. This approach also allows better parallelization (probably).

#### Output Control

The protocol defines how the data is stored. 

The coordinates are written to `.gro` files at each timestep. If a coordinates file was generated for each step, there would be a lot of files, so they are combined in a trajectory file (`.trr` or `.xtc` in the case of GROMACS, `.dcd` in NAMD).

> We write out the coordinates in each file, 2fs each gro file corresponds to 1 ps, so each gro is a frame or step. The time of the ste depended on the frequencies.

Other files in the output are the energy file (binary) and the log file.

#### Bond Parameters: Constraints

The vibrations of the hydrogen bonds are the fastest of the system. By constraining them, we can use timesteps of 2 fs.

#### Nonbonded Settings

Calculating non-bonded interactions in too much detail is computationally very expensive for for big systems such as proteins or in fields such as drug design (tough it is possible to do it with small systems). Then, a cutoff scheme is used.

The cutoff scheme of this protocol is Verlet, and it performs a buffered neighbor search within certain cutoffs:

* Coulomb for short-range electrostatic interactions.
* VDW interactions.

It calculates the non-bonded interactions by combining whatever. This is speeded up by the fact that beforehand, a list of atoms is calculated, and the interactions are calculated based on this list.

> Pairs of atoms have a role here <!--Doubt-->

During simulations, one pair of atoms is above the cutoff, so is not in the list. To prevent a problem, a buffer is used, so pairs with a slightly larger radius can be included in the list.

> I think in this case it has a default buffer, you don't give a temperature.

#### Electrostatics

For long-range electrstatic interactions, calculations are performed in a Fourier space (so it is an approximation).

> In the image, positive charges are red, and you see how actually the PMD density looks like, so its an approx.
> 

> **Some notes about the cutoffs and other values:**
> 
> Cutoff values have been proved for many years, so typically the same ones are used. They are changed when the problem is very specific (this is usually related with temperatures).
> 
> When a lot of data is being generated, it is saved at bigger timestamps. But when a really fast event has to be observed, saving more frames (smaller timestamp) is adequate.
> 
> Fourier spacing and grid spacing: when you have a grid, it is required to define the space between each point that is calculated. The smaller the space, the higher the computation. The spacing values we see have been used for decades.

#### Temperature Coupling

When there are a lot of cutoffs, energy is lost and this adds inaccuracy. This has to be compensated by adjusting by whatever velocities. *More information about this and pressure coupling in the future*.

#### Periodic Boundary Conditions

They are normally used. Maybe they arent if the simulations are in vacuum.

#### Velocities

When starting the simulation atoms do not have velocities yet. They are assigned (randomly, basically).

Velocities are coupled to temperature (the higher, the faster and viceversa). This doesn't mean that at a certain temperature all molecules have the same velocity: we observe that they correspond to a distribution: maybe some molecules are slower and others way faster, but most of them are at an average velocity (or closer to it).

In the following image we observe how the maxwell distribution is for two samples at different temperatures:

* Each follows a different distribution.

In the simulation, the velocities are assigned based on the Maxwell distribution.




* Cold temp: you would have a certain distr, 2 molecules are very slow and 2 very high, the average velocity lies whatever and the higher, the better c_cold and v_warm is the average. Most molecules are on the average.
	* As the system is closed, the number of molecules is limited. There is a top to the molecules in the graph of the distr bc we work in closed systems.
	* So we assign (for xample, to the solvent the average speed of the maxwell distr). Certain temp doesnt mean that the atom has just one velocity, we have to look at the distr of all molecules.
	* Here we use it for assigning velocities, but if you go further and rescale velocities then you still use the distr to be sure that all the system follows a normal distribution.  In the next step you collect the velocities, so they are not reassigned.
	* gen_seed is random, means the vel, the distr, the average vel fulfilled, but it could be that  an atom is faster than other. This makes sure that if yu start form the same system the you can deviate into diff firection and evolve it into diff directions  which we like to see, because better for cnformation space.
	* The assignment of vel is not the same for a heavy or light atom: the average velocity is proportional to sqrt(boltzmann*T/m), so the heavier, the slower.

> If we had a system with 5 molecules, maybe 1 would be slower. Increasing the temperature would increase the velocities, but the slower one might still be slower than the others.


--------------------------
	
#### Analysis of the Simulation in VMD

<!--Audio: from 40:30-->

**Examination the Simulation**

We open VMD and we load the following files from `GROMACS_output/complex`:

* Load the structure file `npt.gro`.
* Load the trajectory `md_01.xtc` into `npt.gro`.

We observe a lot of lines. If we create a representation for *resid 9780* (a water) with CPK, we observe that one of the hydrogen bonds is normal, but the other has a long bond.

The reason is that the hydrogen lies where the box finishes. In the simulation it would be moving to the next periodic box, but in this representation we see it in the other side of the box. 

So this lines are an artifact caused by the periodic boundary conditions. If we rewrapped the molecules around the protein so they won't split, we would see the entire molecules, not the artifact (as they wouldn't split from one image to the other).

We can observe this with a representation of *resid 8429* with licorice: we will observe how it appears on different sides during the simulation.

Next, we create the following representations:

* Protein with cartoon
* Ions with VDW

If we play the simulation, we observe that some ions disappear. If we create a representation for *resid 10772* (an ion) and hide the other representations (except the protein) we will observe how it appears on different sides during the simulation. This is also caused by the boundary conditions.

> Wrapping centers the protein in the middle so it doesn't get splitted. We could also define that ions or waters shouldn't be splitted. This is specially applicable with proteins with more than one atom, as we can see half a molecule in one side and the other half in the opposite. Ins are not a good example, as they are just one atom and don't split. Just "teletransport".

**Alignment and Backbone RMSD as a first step of analysis**

Playing the simulation at a faster speed allows to see a bit of we see the difussion and rotation of the protein. If we want to analyse this, we need to align.

> We align to the backbone. Sometimes you just align to helices or whatever to focus on different things.

We open the RMSD Trajectory Tool and we make sure that the *Backbone* box is checked and that the first frame is the reference. Then, we click on align then to RMSD (the *Plot* box needs to be checked).

The RMSD goes around 1.2, which is normal for a system. It doesn't have big conformational changes, just adjusts a little bit the backbone. *It is important to see if the protein is stable or does crazy things, based on the RMSD*.

**Average RMSD per residue**

The plot we obtained is an average for all backbone atoms. Now we wold like to see *inhibitor region of the protein*, the prot has a chain of residues (from 1 to 100) we want  to see for each residue what is the RMSD. We can define regions that are more flexible.

Calculating the average rmsd is done by looking at each frame (ref frame is usually the first one).

We want to find the average RMSD over time of each residue in the protein. Again, we will consider the first frame as the reference, so the second frame is compared to the first, checks how much the residue changed and puts an RMSD. Also for the third, and so.

> It doesnt matter which frame is the reference, but we always need to compre to the reference.
> 

We use a script that stores data and contains a function that you use to calc the rmsd of the residue over time. To use it, we need to define our selection (what you want to calculate the rmsd over). It is defined in the slide.

* top is used when we have more molecules, you see the id in the window list, we can highlight one with a T, that is the top one.

The RMSD per residue. We'll see the more stable residues and which fluctuate more. See for each residue. The other was the whole structure.

protein and alpha: in the pdb the name calls. the atom type. By selecting c alpha is resume dfor alpha. You could also say name calpha We select just one atom per residue (a residue has usually 10-15 atoms, we pick one  bc we want the resid) If we picked more than one atom, we'll hace a lot of resids per residue (so the name would be repeated). Thats why we select just one atom. We can see in other examples what happens when you select more than one atom per residue.

This is achieved with tcl scripting. We need to call the script:

```tcl
source residue_rmsd.tcl
```

```tcl
cd class4/GROMACS_output/complex
source

set sel_resid [[atomselect top "protein and CA"] get resid]
```
> : gets a smaller list, just CA

An erroneous way to do the selection is:

```tcl
set sel_resid [[atomselect top "protein"] get resid]
```
> You get 20 times '1', because the residue has 20 atoms, for each atom you get the resid.

```tcl
set sel_resid [[atomselect top "protein and CA"] get resid]
```

> The list here is smaller, as we are selecting just alpha carbons.

> **View more about the residues by looking at slide 35**
> 

Once we set the selection, we pass it to the function with function name + structure + selection:

```tcl
rmsd_residue_over_time top $sel_resid
```

As a result we will obtain a file called `residue_rmsd.dat` (in the `complex` folder). We have to plot it to identify residues with large conformational changes during energy optimization.

> Seems the terminal also gives the list of the average RMSD per residue, but we use the printed file for convenience.

From the plot, we observe that some residues with high RMSD are 15, 68, 35, 114, 128 (they appear as peaks).

> gnuplot plot dat with line

Represent them with VDW and colored by name and play the simulation.

```
resid 15 or resid 35 or resid 68 or resid 114 or resid 128
```

Where are the residues with a high RMSD located?

The residues with a high RMSD are located in the surface (including the coils). This residues interact with the solvate and are very flexible. These residues are the ones that normally are problematic to crystalise. In the crysta they usually have high V values. Match experiment with analysis: you can explain x-ray struct with this.

**Showing the Flexibility**

The flexibility of a residue can be made visible by simultaneously showing multiple conformational states along the simulation.

Create a representation with the residues with high RMSD with drawing method "Licorice". Remove the hydrogens from the selection.

```
(and no h)
```

Then go to the Trajectory tab and in the field "Draw Multiple Frames" type `0:10:102`, which means write from frames 1 to 102 with a step of 10.

This can also be done with the protein. 

**Electrostatic Interactions**

Salt bridges have an important contribution to protein stability. They can be found with Extensions > Analysis > Salt Bridges.

This generates a series of files with information about the salt bridges, which can be plotted.

Substrates and ligands often use a water channel to reach the active site of the protein. To identify water channels, we can generate an occupancy map of water (Extensions > Analysis > VolMap Tool).

> We select water within 3 of protein so we don't have all the representation filled with water. The volmap type is "occupancy" because it is easy to interpret (1 = always there, 0.5 half time of simulation). 

<!--Insert screenshot-->

It will create an Isosurface representation. Color it with colorid 15 to distinguish from the protein.

Then, create a VolumeSlice Representation. Move it along an axis and identify better the axis. Try different axis (in this case, it is the Y). We will see here the water enters: that is the pocket.

We want to know if the pocket is hydrophobix or polar. To do this:

* Create a representation of the protein (lines colored by name)
* Check residues along the pocket.
* Create a representation with those residues to highlight them (use the resid, VDW and color by name).

--

> Seems some residues she selects are:
> 
> * trp28, leu56, trp 108, , 31, 8, 12
> 
> The residues I select are: trp108, leu56, met 12, ile 55, , glu 35, ala 31, trp28
> 
> ile 55 no tan cerca
> 

List of hydrophobic residues: ALA LEU VAL ILE PRO PHE MET TRP

trp108 (A), leu56 (A), met 12, ile 55, , glu 35, ala 31, trp28

The pocket is hydrophobic


> Nice pocket to target. when you see water in a hole you expect certain feature. If it is for instance hydrophobic the ligand dispaces waters and we get a nice affinity


put waters areound in a repr to see how they go in and out

water and same residue as within 3 of resid (update selection for every frame and color every frame)

Sum of the missing parts: ion binding sites search with volmap. Ions lie around or stably bind somewhere? If so, you would see it in the occupancy map for the ions. this seems to be q5

## Class 5. Application of MD to a Membrane Protein

<!--Date: 25/04/2019-->

We use [CHARMM-GUI](http://www.charmm-gui.org/?doc=input) to set up the simulation. After that, the equilibration and production steps are done with GROMACS. Finally, the simulation is analysed.

### Simulation Setup

When using CHARMM-GUI, the following steps are followed:

![](msi-notes.assets/5.1.png)

#### Step 0: Structure Preparation

VMD is used to prepare the protein ([3PBL](https://www.rcsb.org/structure/3pbl)) in this case.

When we open it, we observe that it has two receptors. Normally, the one that has less missing segments/atoms structure is taken. We will keep chain A, as there is not much difference.

> Represent chain A as NewCartoon to see which one it is.

To save the chain into another file, we open the TkConsole:

```tk
set sel [atomselect top "chain A"]
$sel writepdb 3pdbl_chainA.pdb 
```

> This script just takes the chain A of structure marked as top. The `atomselect` key is used to select parts of a structure. The selection is saved into a variable, which later is used to generate the PDB.
> 

We load the new file and observe it:

So the chain A, after loading, will look like this:

![](msi-notes.assets/5.3.png)

It still contains the T4 lysozime (introduced into the receptor to facilitate crystallisation).

In order to remove it, we use again Tk, but first we need to know which residues are the lysozime and with the chain. The residue IDs can be checked with the Sequence Viewer. In it, we can see how the residue number goes from 221 to 1002 suddently: that is the lysozyme. *Prot goes *

![](msi-notes.assets/5.4.png)

![](msi-notes.assets/5.5.png)

The selection is done with:

```tcl
set sel [atomselect top "resid < 1000 or resid > 1200"]

$sel writepdb 3pbl_chainA_cut.pdb
```
> We know thay the protein goes until residue 222 and then restarts at 319, and also, over 1200 is not the lysozime.
> 

The resulting structure is shown in the image, but if we create a VDW representatin with "not protein" as selection, we will see small molecules:

![](msi-notes.assets/5.6.png)

Those need to be removed:

```tcl
set sel [atomselect top "protein"]
$sel writepdb 3pbl_chainA_cut_protein.pdb
```

**<New chain ID's**

With the removal of the lysozyme, a broken loop appeared (221-319, 100 residues difference between them). This is a bit artificial and has potential to mess up, so we will assign a chain ID to each chain, and save it:

```tcl
set sel1 [atomselect top "resid 1 to 221"]
set sel2 [atomselect top "resid > 318"]

$sel1 set chain A
$sel2 set chain B

set sel [atomselect top "protein"]
$sel writepdb 3pbl_chainA_final.pdb
```

Result:

![](msi-notes.assets/5.7.png)

### Simulation Setup

**Upload and Selection**

The procedure starts by going to the CHARMM-GUI membrane builder and uploading the structure. Then it will allow to select the chains (not selected by default, make sure that both A and B are selected).

![](msi-notes.assets/5.charmm-1.png)

![](msi-notes.assets/5.charmm-2.png)

**PDB Manipulation**

The gap that we were talking about before generates two N and C-terminal ends. To avoid charged artifacts with them, they are capped with ACE and CT3, so they remain neutral.

The disulfide bonds need to be identified. To know where they are, we go back to VMD:

* Select all cysteine residues
* Draw method: licorice
* Color: by name

If we remove the rest of representations, we will see that some cysteines are together: those are the bonds. We can select the atoms to see the residue ID.

![](msi-notes.assets/5.charmm-3.png)

![](msi-notes.assets/5.charmm-4.png)

**Generate PDB and Orient Molecule**

Orienting the protein because...

In the case of GPCRs, we use "Align the First Principal Axis Along Z".

Lastly, the D3R has a pore were small ligands bind. This has to be filled with water. We do it with the option *Using protein geometry*.

![](msi-notes.assets/5.charmm-5.png)

<!--min 41.-->

**Calculate Cross-Sectional Area** <!--min 44-->

When the calculation is performed, two plots and some characteristics are returned, which say if the protein was well oriented or not:

* Protein Top Area: 1190.41383
* Protein Bot Area: 1139.70367
* Pore Minimum Radius: 0

![](msi-notes.assets/5.charmm-pore.png)

![](msi-notes.assets/5.charmm-cross.png)

**System Size Determination**

In this step, we choose the box type and its size:

* The box will be rectangular.
* Length of Z will be based on a water thickness of 12 (A or nm??<!--doubt-->. This is done to maintain the system as small as possible.
* Length of X and Y will be based on the ratios of lipid components, but we need to guess a size. We put 80 (80x80).

Next, there is a table with lipids, wich contains lipid types. The upperleaflet and lowerleaflet ratios refer to each side of the membrane. We choose POPC as the lipid, with a 1-1 ratio (*same amount of lipids on each side*).

Finally, we click on the "Show the system info" button. With that information, we can know if the system is to small.

![](msi-notes.assets/5.charmm-6.png)

If we click on `step3_packing.pdb`, we can see the result so far.

![](msi-notes.assets/5.charmm-7.png)

**Generate Components**

![](msi-notes.assets/5.charmm-8.png)

> There are two popular methods commonly used to build a realistic protein/membrane complex. In the first method, lipid-like pseudo atoms are first distributed around a protein and then replaced by lipid molecules one at a time [19]–[21]. Individual lipid molecules are randomly selected from a lipid library that contains various conformations of lipid molecules. This method allows one to easily control the system size and the number of lipid molecules while it generates a lipid bilayer nicely packed around the protein. In the second method, a hole is first created in a pre-equilibrated lipid bilayer and then the membrane protein is inserted into the hole [22]–[24]. In general, weak repulsive radial forces perpendicular to the membrane normal are applied to a lipid bilayer until the hole is large enough to accommodate the protein. This method provides a well-equilibrated lipid bilayer, and one might expect less equilibration time than in the first method. For the sake of convenience, the first method is called the “replacement method” and the second method is called the “insertion method” hereinafter. Although both methods are well explained in the literature [19]–[24], considerable efforts and experiences with MD simulation software are required to build a realistic protein/membrane system.
> 
> From [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000880)
> 

Q2: 



----

Check lipid penetration

The protein surface penetration check finds the lipid tails that go beyond the protein surface, and the lipid ring penetration check detects the lipid tails that pass through the cyclic groups (e.g., cholesterol ring) in the simulation systems. Energy minimization can resolve many of these bad contacts, but one might need to visually check the following lipid molecules to ensure the following contacts are resolved. The user can regenerate the lipid bilayer if necessary.
Protein surface penetration:

No protein surface penetration is found.
Lipid ring penetration:
No lipid ring penetration is found.
Building Ion and Waterbox

Membrane components are generated. Due to time constrains, we first generate the lipid bilayer then generate ions and the water box. Click "Next Step" to generate ions and the water box.

----

**Define equilibration conditions** <!--Before min 1:05-->

![](msi-notes.assets/5.charmm-9.png)

Finally, we download all files.

![](msi-notes.assets/5.charmm-10.png)

### Equilibration and Production

In the `2_equilibration_and_production` directory, there are some restraints of the system that are applied during the equilibration. They are:

* BB: force constant to keep the backbone atoms constraint.
* SC: force constant to keep the sidechain atoms constraint.
* wforce: force constant to keep the water molecules away from the hydrophobic core.
* tforce: force constant to keep the lipid tail below +/- %.
* mforce: force constant to keep the lipid head groups close to target values.

The two first retraints apply to the protein. Then, the lipids are packed around the protein (before there were some empty spaces) and water is prevented from entering the membrane.

![](msi-notes.assets/5.equilibration-production-1.png)

![](msi-notes.assets/5.equilibration-production-2.png)

Now look ar her files, the folder we said. Inside, gromacs. Interested in steps 6 to 7. Around min 1:18.

Proa_rest_itp: fc_bb is force constraint backbone.

`toppar` contains stuff about the molecules found in our simulation.

Fixed only one atom of the head. You don't allow to move it from z so it doestn move up and down: they just diffuse along x and z.

There is a slide on how to run it.

![](msi-notes.assets/5.equilibration-production-3.png)

### Analysis

The function of many GPCRs is sensitive to sodium ions. We will: 

* Detect sodium ion binding sites by computing the occupancy map.
* Generate a figure showing residues that bind the ion and label these residues (e.g. D400).
* Generate a plot that shows the ion entrance over time by plotting the Z coordinate over time.

Open the `.psf` with VMD and then load the `.xtc` into it.

> Tip: as the simulation has too many frames, we can make it less comoutationally intensive by decreasing the number of frames: right click the simulation and choose "Delete frames". In stride put two, so it reduces the number of frames by half:
> 
> ![](msi-notes.assets/5.vmd2-2.png)

A first observation is that waters are very dynamic, as they can travel a lot. But we don't care about this because we don't analyse water.

The first thing to do is aligning the backbone of the protein with the RMSD Trajectories tool (*check the backbone box*).

![](msi-notes.assets/5.vmd2-1.png)

Next, a VolMap is generated with:

* Selection: resname SOD and within 5 of protein
* Type: occupancy

![](msi-notes.assets/5.vmd2-3.png)

When this is done, go to the generated representation and plot it as isosurface with a value of 0.01. It is also a good idea to make it translucent instead of opaque.

![](msi-notes.assets/5.vmd2-4.png)

In the simulation we saw that one of the ions entered the channel, stayed at one point and then goes further inside. The VolMap confirms this behaviour: we observe two very defined binding sites.

To plot the ion entrance over time, the following steps are needed:

* To get the index, select it and go to labels. In this case, it is 56129 (this can change for each simulation).
	![](msi-notes.assets/5.vmd2-1.png)
* Using a script, select your ion and obtain the Z coordinate:
	```
	set sel [atomselect top "index 56129"]
	$sel get z
	```

The last step gives the same information as the labels window: the Z coordinate for the current frame. In order to get the coordinate for each frame, a loop is needed:

```tcl
# Open file to write
set file [open "sodium.dat" w]

# 1. Define your selection
set sel [atomselect top "index 56129"]

# Second: know how many frames:
set n [molinfo top get numframes]

for { set i 0 } { $i < $n } { incr i } {
# update selection for current frame
	$sel frame $i
	$sel update

# assign z coord
	set position [$sel get z]

# Write frame and position to file
	puts $file "$i $position"
}
	
close $file
```

> The selection must always be updated, if not it will go to the last frame.

## Class 6: Markov-Based Analysis of Biomolecular Systems

<!--Date: 26/04/2019-->

Toni Giorgino, from the National Research Council of Italy. [Website](giorginolab.it), [email](mailto:toni.giorgino@cnr.it),[GitHub](https://github.com/giorginolab).

### Introduction

States: when the system overgoes a transition (folded to folded, but also intermediate: ie open/closed but lso occupied or not). Normal sim has times going forward. We want to analyse where systems are passing to.

Kinetic info: the time it takes to transition.

Umbiased: we dont force anything, we use hysiological conditions and observe it forward in time. 

Due to tech (GPU) we can go to microseconds.

--------


Difference between d1 and d10: same motion, but observed at a time of (...). Kinetics are accelerated. Timescales for exchanges are cchanged. Same potential? same eigenvector.

----

State 0 is never observed.

Boltzmann inv: free energy profile. No kinetiks, no markov yet. O count states are not shown.

-----

Markov

Define tau first, i set to 10. Trick he uses to count stuff.

The TT matrix means how many times state has been seen followed bystate (colname). Ie state 31 followed by state 37 (col) is x times.

```{r fig.height=4}
ei <- eigen(t(P)) # all analysis in this list

# split the list
eva <- sort(abs(ei$values),decreasing = T)
eve <- ei$vectors

taus <- -tau/log(eva[-1]) # this means forget the first.

plot(eva,ylab="Eigenvalue magnitue",xlab="Eigenvalue index (sorted)")

plot(c(NA,taus),ylab="Implied timescale")
abline(h=tau,col="red")

```

eva: always one 1, the next is relaxation timescale. fast process that not so interest. we are interested in slower

eve: not needed an abs values.


if we see the taus, we see that the smallest?? relaxation is 3610. If the barrier is low (barrier for going left to right or whatever), the number is smaller. More barriers: we would found the numbers next to ours

We see in the implied... plot that only one main feature in the model, the others are faster and not so interestng.

-----

Eigenb¡vectrs: eval

eve: take the 1st column, take the real part, and then take -log (probabilities to free energies conversion), and whatever

then plot it: double well

# Part 2. Beyond Classical MD

## C8. Enhanced Sampling Techniques

<!--Date: 04/05/2019-->

### Introduction

Computer simulations of biomolecular systems have grown rapidly over the past few decades: from small molecules in vacuum (>20 atoms) to large protein complexes in a solvated lipid bilayer environment (>20000 atoms). 

However, despite its success, MD simulations are still limited in two regards:

* Inaccuracy of force fields: 
	* Force fields are simplifications based on collections of experimental data and ab initio calculations of how the system behaves.
	* Over the years have been refined, as this was needed to achieve longer simulation times with larger systems.
	* They describe the main energy function, but still can have some inaccuracy with certain kinetics.
	* Note that the stability of a system with time depends on its size: the smaller they are, the longer they can be simulated.
* High computational cost: it is needed half a year to simulate on the millisecond level, unless a supercomputer is used. 100 μs simulation of relatively small systems (approximately 25,000 atoms) running on state-of-the-art computing architecture requires a month of computation to complete.

Biological molecules are known to have rough energy landscapes, with many local minima frequently separated by high-energy barriers, as the following image of a protein folding process shows:

<img src="msi-notes.assets/8.1.png" alt=""
	title="" width="500"/>
	
It making it easy to fall into a non-functional state that is hard to jump out of in most conventional simulations. Replication can be used to escape those minima, but the difficulty increases when the event you want to observe is rare.

Also, we are not only interested in the global minimum, but rather in the ones that are biologically functional, as they can be relevant events. Examples of this are:

* More than one pathway to arrive to the same conformation: a main one and other biologically possible ones that are less used. Gathering information about those could be interesting.
* In case of the transport through membranes, channels and transporters have to undergo large conformational changes in the course of gating substrates.

But sampling those conformational states is difficult because of the high energy barriers, and limits our ability to analyse and reveal functional properties of the systems being examined.

The following image shows one of those states:

<img src="msi-notes.assets/8.2.png" alt=""
	title="" width="500"/>

In conclusion, escaping local minima is complicated and time consuming, and straightforward MD simulations cannot handle this. Algorithms that are able to sample conformational states are needed. They are divided in two categories:

<img src="msi-notes.assets/8.3.png" alt=""
	title="" width="500"/>

> CVs are collective variables, which are used by most sampling algorithms to describe the progress of a pathway. An example is "distance between the ends".

### CV-Dependent Methods

#### Metadynamics

This method was developed by Parrinello's group in order to improve sampling of systems where ergodicity is hindered by the form of the system's energy landscape by adding memory to the system. *A. Laio, M. Parrinello, Escaping free-energy minima, Proc. Natl. Acad. Sci. U. S. A. 99 (2002) 12562–12566*

**Ergodic Processes**

A process is ergodic when its statistical properties can be deduced from a single, sufficiently long, random sample of the process. For example, in a small system for which you can describe the whole phase space (its conformations) assuming you have enough time.

In the following image, the green rectangle is a phase space, and the yellow is another:

<img src="msi-notes.assets/8.4.png" alt=""
	title="" width="350"/>

* It is possible to sample the green it entirely because it is small. This represents how most molecular systems are.
* Sampling a part of the yellow doesn't allow to arrive to green space, as they are not connected. This applies in cases such as bond breaking. Quantum mechanics allows to define bond breaks, unlike with molecular mechanics.

Note that the big space could be subsetted into smaller spaces to solve the problem, but then it is not a single sample, so it cannot be considered ergodic.

**Procedure of Metadynamics**

We start by defining the CV and plotting the energy landscape of the system. In the example below we see a structure that was crystallised at a stage that is a local minimum (1st graph).

<img src="msi-notes.assets/8.5.png" alt=""
	title="" width="400"/>
	
After a certain time defined in the simulation protocol (for example, after 3000 integration steps of 2 or 4 fs, which are 6 ps).

When this time passes, a defined amount of energy is added. This energy is not just a value (which would be represented in the plot as a vertical line), but rather a gaussian distribution with a defined height and width. That is why the increases of energy are called hills.

As we add a hill, the potential energy function is increased. If we added just one, the structure would fall to the minimum again, then it would climb again and after 6 ps, fall again. To prevent this, more hills are added until the vale is filled and the structure cannot go to the bottom. The consequence of this is that the energy barrier is decreased. We keep the information about the number of hills added and their shape, so we know how much energy we are adding to the landscape. *I guess this means that we don't end up with a biased landscape*.

As the hills remain in the vale once it has been filled it prevents that it is resampled. That's why it is said that this method **introduces memory to the system**.

Defining the shape of the hills can be done with big energy values, allowing to get a rough idea of the landscape. To be more accurate, smaller values are used. The bigger the value, the faster a vale is filled.

Convergence in metadynamics is achieved when everything has flattened out, so any hill you add results in staying in the same place. In papers, convergence is demonstrated by showing the last plots of the profile: if they are the same, convergence is achieved. It is also possible that maybe one part has converged and others not. In that case, it would be acceptable to say that there is convergence if we are only interested in that part.

The procedure is computationally fast, as 6 ps is nothing, while with unbiased simulations we could never scape the minimum. At the end, we will have information about the whole energy landscape.

The idea of the process can be seen in the following [video](https://www.youtube.com/watch?v=IzEBpQ0c8TA).

The following example has two collective variables plotted. The entrance of ions is studied: `a`, `b` and `c` are different places where an ion can be found, with energetic barriers between each of them.

<img src="msi-notes.assets/8.6.png" alt=""
	title="" width="500"/>

What they did was to filled up the landscape for the Z coordinate to follow the entrance. This was done with respect to the torsion angle of one chain. They saw that between `a` and `b` there is 1 kcal of difference to overcome to pass from one state to the other.

The cause of the high barrier from `b` to `c` can be explained structurally with the following image. We see that there is a residue blocking the way. When energy is added, the ion can pass through because a residue moves.

<img src="msi-notes.assets/8.7.png" alt=""
	title="" width="350"/>

**Summary**

Metadynamics does depend on a low dimensionality of the system in order to produce an accurate description of the free energy surface, therefore using a small set of collective coordinates is essential. Ideally one CV is used, two if it is needed for a better description, but three makes it very difficult to converge.

Such characteristics allow this method to be quickly used to provide qualitative information about the overall topology of the free energy surface being examined.

#### Umbrella Sampling

To calculate the energetic landscape using umbrella sampling, the pathway needs to be known beforehand. The following example shows how a molecule moves closer to the other structure, until it binds. Each state corresponds to an energy, where the bound state is the global minimum and the unbound to a local minimum, with a barrier in between.

<img src="msi-notes.assets/8.8.png" alt=""
	title="" width="250"/>
	
By taking a sample the image below, then we would not get the binding. 

<img src="msi-notes.assets/8.9.png" alt=""
	title="" width="250"/>

Instead, the method starts from snapshots along the pathway (unbiased simulations), and energy or force is added in order to sample around an area. Simulations for the different conformations are then run in parallel. Note that the snapshots need to overlap.

<img src="msi-notes.assets/8.10.png" alt=""
	title="" width="250"/>
	
Once you know the landscape, you know the difference between the bound and unbound states. Thus is interesting for drug design, as the affinity of a drug is higher when the energetic barrier is big. This is usually solved by trying to lower the barrier (some cofactors do this). <!--There is other possibility instead of lowering the barrier, didn't listen to it.-->

<img src="msi-notes.assets/8.11.png" alt=""
	title="" width="250"/>
	
**Comparing Metadynamics and Umbrella Sampling**

While metadynamics and umbrella sampling are conceptually similar techniques to overcome free energy barriers, they have differences:

* Metadynamics is better suitable for finding reaction pathways. However, potential of mean force (PMF) calculations are highly dependent on input parameters (including the height and width of the Gaussian and !"). The appropriate choice for these parameters is crucial for accurate calculations.
* Umbrella sampling is useful to calculate accurate PMF. *It allows calculating very accurate energetic landscapes, as you indicate exactly what you want to sample in each simulation.*

A possible strategy would be to use metadynamics to know the pathway and have a general idea, and then pull from here the states we are interested in and do umbrella sampling. An alternative would be to generate artificially the pathway by pulling the ligand out. This can be done manually, but is not as good as you are defining the pathway, and also it doesn't always work, as some ligands first bind a recognition site, and then the actual binding site.

### CV-Free Methods: Replica Exchange

A set of non-interacting replicas runs at different values of an exchange variable, usually temperature (T-REMD), but there are other alternatives (pH, pressure, etc). At specific intervals, replicas at neighbouring values for the exchange variable are swapped randomly.

<img src="msi-notes.assets/8.12.png" alt=""
	title="" width="400"/>

In the case of temperature, the fact that it increases can make the simulation overcome the energetic barrier, while decreasing it means that the well is sampled better. Temperature can also be good in unbiased simulations to overcome the energetic barrier, as it adds energy and might make the system more flexible and with more movement. But in some cases, the system can misbehave, so lowering the temperature is a good idea.

In an efficient run, all trajectories will experience changing of the exchange variable value. At each value for the exchange variable, the trajectories will be discontinuous, but follow a proper Boltzmann distribution for the specific value being exchanged.

The main advantages of this method are a gain of flexibility, overcoming ergodicity, and also a good parallelisation.

Applications:

* Sampling of different conformational states (e.g protein folding). In protein folding it has been used as a benchmark by comparing how quick and accurately the folding is reproduced.
* Study of protein protonation states, which might have a role in the protein active conformation.
* Usage in non-ergodic systems such as ... <!--mdc-->

<img src="msi-notes.assets/8.13.png" alt=""
	title="" width="250"/>

### PLUMED

PLUMED is a consortium to provide code for different enhanced sampling approaches for different simulation software. Here is an [explanative video](https://www.youtube.com/watch?v=PxJP16qNCYs) and the [installation procedure](https://www.plumed.org/doc-v2.5/user-doc/html/_installation.html). To incorporate it into GROMACS, a patching procedure is needed.

![](msi-notes.assets/8.14.png)

## C10. Umbrella Sampling Simulations

<!--Date: 16/05/2019-->

Most trajectories are not ergodic. 

Steered MD is needed in most cases for umbrella sampling.

Sampling time issue

Unbiased collective variables don't touch the energies.

In umbrella, the energy is untouched. Allows to estimate free energy.

The initial trajectory is to separate the dimer and do the sampling along the separation.

gmx make_ndx init.gro

We create the groups of the protein and kill ourselves he keeps repeating something about dopamine receptors.

Use the gro file, the index file from the previous step, `A2A.ndx` is the output, and then use cat to merge both indexes.


At some point he does 14 | 15 | 16 to make whatever groups.

He does a pulling that makes sense.

Set up US windows

Use gmx dist to compute distance per configuration (doing like a for loop, he has a script somewhere in the folder).

End up deciding which windows you will use (slide 15). Not trivial, you have to play around

It doesnt matter if we have different energies in the windows, as when we finish we remove the bias.

S21

We are applying the bias in 

S24 Here we make an histogram of each window bc yes.

bsResult has a low energy, the dimer, as it is the most stable conformation. We can change the orientation of the dimer in order to see how this profile changes.

bsProfs is used for another stuff. we see that is has converged. It doesnt converge if the fat part is lower, but i don't get why. This is the attempt of the bootstrapping i think.

bsresult shows error, if the error is big, then we should sample that window for longer in a publication, the name is `profileX.xvg`.

the pmf is the real one


name BB "SC*:"
resname DSCP

both as points, make them bigger, etc

## C11. Interactive and Steered Molecular Dynamics

As replica exchange requires a recompilation of GROMACS, adding a layer difficulty, we will do steered molecular dynamics (SMD) instead.

SMD is applied to processes such as the unfolding pathways of proteins. Forces are applied to one side of the molecule, and its behaviour can be analysed. SMD allows the calculations of the **potential of mean force (PMF)** from the trajectories obtained.

Also, with **interactive molecular dynamics (IMD)** we are able to apply the forces to the molecule and see its reaction in real time.

For this class, we will perform the simulation in vacuum so it can be executed in our machines, but this is not how it is done actually. The system we will use is deca-alanine, a peptide composed by ten alanine residues that forms an alpha-helix in vacuum as a stable conformation (top figure), stabilised by six hydrogen bonds. The molecule will be stretched by applying external forces with IMD and SMD: this will make the molecule undergo a gradual conformational change from the alpha-helix to the random coil (bottom figure).

![](msi-notes.assets/11.1.png)

### IMD

The directory contains the following files:

* da.psf: protein structure
* imd ini.pdb: initial coordinates
* par all27 prot lipid cmap.prm: CHARMM parameters
* imd.namd: NAMD configuration
* imdfixedatoms.pdb: fixed atoms

In the `imd.namd` file we can see the simulation conditions. It defines the atom that will be fixed, by setting to 1 the beta-value of that atom on the pdb file. It also contains settings in order to allow the interactivity: it communicates with VMD through a predefined port.

> The beta value gives an idea on how flexible the side chain is. High values happen with low resolution molecules or areas.

```
# Fixed Atoms Constraint (set PDB beta-column to 1)
if {1} {
fixedAtoms          on
fixedAtomsFile      imdfixedatoms.pdb
fixedAtomsCol       B
}

# IMD Settings (can view sim in VMD)
if {1} {
IMDon           on
IMDport         3000    ;# port number (enter it in VMD)
IMDfreq         1       ;# send every 1 frame
IMDwait         yes     ;# wait for VMD to connect before running?
}
```

To start, run NAMD:

```
namd2 imd.namd > da_imd.log &
```

The simulation will not actually run until the connection between NAMD and VMD is established, so we must execute it:

```
vmd -e imd.vmd
```

This line loads the molecule into VMD and shows the ribbon and CPK representations. The alpha carbon of the first residue is colored in orange: that one is the fixed atom.

Next, we go to Extensions > Simulation > IMD Connect (NAMD) and we type localhost as hostname and 3000 as port. When clicking on Connect, the simulation will start.

We can apply forces to stretch the molecule with Mouse > Force > Atom and then click and drag on an atom in the opposite end of the one with the orange atom. A red arrow will appear, which represents the magnitude and direction of the force that is being applied.

We can stretch the molecule to half of its original length. If when it arrives to that point we remove the force (by middle-clicking on the atom to which the force is being applied), we will see how the molecule holds back.

### SMD

IMD is useful to explore the system, but SMD allows a more systematic way to apply forces and analyse the system.

The files needed are:

* da.psf – protein structure
* smd.namd – NAMD configuration
* smd.tcl – Tcl script
* par all27_prot_lipid_cmap.prm – CHARMM parameters
* smd_ini.pdb – initial coordinates

The `smd.tcl` script contains the external forces that will be applied. It is referenced in the extra parameters section of the `smd.namd` file, where we see a "Tcl interface":

```
########################
## EXTRA PARAMETERS   ##
########################

# Tcl interface
tclForces           on
tclForcesScript     smd.tcl
```

Those external forces design that:

* One end of the molecule (the N atom of the first residue) is constrained to the origin.
* The other end (the capping N atom at the C-terminus) is constrained to a point that moves along the z-axis from 13 Å to 33 Å with a constant speed of 1 Å/ps.
* An harmonic potential with a force constant of 7.2 kcal/ (molÅ^2 ) is used for the constraints. <!--Doubt: to a particular atom?-->

In other words, we want the molecule to move from `(0, 0, 13)` to `0, 0, 33` by one A per picosecond, meaning that the full extension will take 20 ps. So we see that what is defined is the end point and the speed to arrive there, not the force itself.

Next, the simulation is run with the following command:

```
namd2 smd.namd > da_smd.log
```

The following output files are generated:

* da_smd.log: standard output
* da_smd.dcd: trajectory
* da_smd tcl.out: Tcl script output

The next step is to analyse the SMD trajectories with VMD. We have a simulation time of 20 ps, as we wrote each picosecond.

To load the trajectory:

```
vmd da.psf da_smd.dcd
```

**Hydrogen Bonds**

The first thing we want to observe is how the hydrogen bonds are broken in the helix-coil transition. We can monitor them with VMD with the following procedure:

* Choose CPK from Drawing Methods
* Create Rep and select HBonds from Drawing Method
* Change the parameters to the following:
	* Distance Cutoff: 4.0
	* Angle Cutoff: 40
	* Line Thickness: 10

Then, the timeline tool (Extension > Analysis > Timeline) is used to follow the hydrogen bonding.

**Analysis of the SMD trajectory**

In the Tcl console, we open the Tcl output of the simulation, storing its content in the variable `mytraj`:

```
set infile [open da_smd_tcl.out]
set mytraj [read $infile]
close $infile
```

The file contains three columns: time (in ps), extension (end-to-end distance in Å) and the applied force ((kcal/ (mol Å)). The distance tends to increase, and the force is overall positive. 

The columns need to be asigned to variables to proceed the analysis, which can be done with the following script:

```
source read.tcl
```

Then, we define the parameters `v` (pulling velocity) and `dt` (time step of the data, which is stored each 0.1, as seen in the first column):

```
set v 1
set dt 0.1
```

Recall that one end of the molecule was constrained to the origin and the other end was constrained to a moving point. The distance between the two constraints was changed from 13 Å to 33 Å with the constant velocity `v`.

Then, an array `c` is created to represent the distance between the two constraint points at each time step:

```
set c {}
set i 13
while {$i < 33.1} {lappend c $i; set i [expr $i + $v * $dt]}
```

We want to plot two things:

* Time `t` vs extension `z`.
* Time `t` vs distance between the two constraint points `c`.

This is fone with the multiplot plugin of VMD, by typing in the console:

```
package require multiplot
set plot1 [multiplot -x $t -y $z -plot]
$plot1 add $t $c -plot
```

Here, we observe two lines: one is a straight line, which represents the distance between the two constraint points, and the other is the actual extension. This is the pulling velocity and how the molecule adapted to the force, and it usually lags behinf the straight line. This is due to the fact that before it reaches that distance, the next force is applied. By decreaing the speed, we would give more time to the molecule to reach that arget value (as what we define are speeds, not forces).

There are times where the force actually moves the molecule past the value of distance we wanted for that timestep. Then, the system tries to adapt by applying a negative force, so it goes back.

We can see the force that is being applied each time with the force-extension curve, where the x axis is the extension (from start 13 to end 33) and the y axis is the force:

```
set plot2 [multiplot -x $z -y $f -plot]
```

Both plots can be compared, as we see that when the extension/constraint curve is above the target, the force in the force curve is negative.

**Potential of Mean Force (PMF) calculation**

The work done on the system during the pulling simulation is:

Where `v` is the pulling velocity.

The work can be calculated with the followinf Tcl script by numerical integration, which will store the work in a list with variable name `w`:

```
source calcwork1.tcl
```

If a pulling simulation is performed very slowly, then the process is reversible. The work done during such a reversible pulling is equal to the free energy difference of the system between initial and final states:

Jarzynski’s equality (Wikipedia):

In thermodynamics, the free energy difference $\Delta F = F_B - F_A$ between two states A and B is connected to the work W done on the system through the inequality: $\Delta F \leq W$, with equality when one takes the system from A to B infinitely slowly (such that all intermediate states are in thermodynamic equilibrium).

The exact PMF obtained from a reversible pulling is stored in the file `Fexact.dat`. Read the content of the file, and store it in a variable called `Fexact`, by running the follwong Tcl script:

```tcl
source load-Fexact.tcl
```

This array contains two columns: the extension and the PMF, respectively. We can plot `W` vs `c` in blue, with `Fexact` in red:

```tcl
set plot3 [multiplot -x $c -y $w -linecolor blue -plot]
$plot3 add $Fexact(1) $Fexact(2) -linecolor red -plot
```

The trajectory used is not practical for the PMF calculation, as the pulling speed was too high. It could be reduced to 0.01 Å/ps.