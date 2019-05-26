# MSI Part 2. Beyond Classical MD

* [MSI Part 2\. Beyond Classical MD](#msi-part-2-beyond-classical-md)
  * [C8\. Enhanced Sampling Techniques](#c8-enhanced-sampling-techniques)
    * [Introduction](#introduction)
    * [CV\-Dependent Methods](#cv-dependent-methods)
      * [Metadynamics](#metadynamics)
      * [Umbrella Sampling](#umbrella-sampling)
    * [CV\-Free Methods: Replica Exchange](#cv-free-methods-replica-exchange)
    * [PLUMED](#plumed)
  * [C10\. Umbrella Sampling Simulations](#c10-umbrella-sampling-simulations)
  * [C11\. Interactive Steered Molecular Dynamics](#c11-interactive-steered-molecular-dynamics)
    * [IMD](#imd)
    * [SMD](#smd)

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

## C11. Interactive Steered Molecular Dynamics

As replica exchange requires a recompilation of GROMACS, adding a layer difficulty, we will do steered molecular dynamics (SMD) instead.

SMD is applied to processes such as the unfolding pathways of proteins. Forces are applied to one side of the molecule, and its behaviour can be analysed. SMD allows the calculations of the **potential of mean force (PMF)** from the trajectories obtained.

Also, with **interactive molecular dynamics (IMD)** we are able to apply the forces to the molecule and see its reaction in real time.

For this class, we will perform the simulation in vacuum so it can be executed in our machines, but this is not how it is done actually. The system we will use is deca-alanine, a peptide composed by ten alanine residues that forms an alpha-helix in vacuum as a stable conformation (top figure), stabilised by six hydrogen bonds. The molecule will be stretched by applying external forces with IMD and SMD: this will make the molecule undergo a gradual conformational change from the alpha-helix to the random coil (bottom figure).

![](msi-notes.assets/11-1.png)

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
set infile [open da_smd_tcl.out
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