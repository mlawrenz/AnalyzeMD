Suggested initial workflow:

1) Extract your MD from the desmond trajs. Use the protein only for
protein-ligand analysis. To analyze solvent you must run Extract MD again with the --solvent flag.
This will be faster to produce the dcd file but may contain jumps in the
protein. 

2) You will get some reference structures from the trajectory extract, these are
very important and require specific AMBER atom ordering and format, and will be needed for all analysis. 
Align this to a reference (-r) xtal.pdb using AnalyzeMD.py pdb_align. Do not just use PyMol,
the output PDB reorders hydrogens and will break the analysis.

You can load in trajectories into vmd to look at the movie directly from the
command line wiht this command:
vmd -pdb xtal.pdb -dcd traj.dcd

You the "Analysis --> RMSD Trajectory Tool" to align the trajectory.

3) Run rmsd_calc using rmsd-all, load the pml file in rmsd-analysis to visualize
residues that deviate the most from starting position (red) color. Using
rmsf-all will show you residues that are most *flexible* throughout the
simulation (most variance from the average). The two metrics tell you different
things. Useful to start an analysis.pse file by launching loadbfactor.pml.

4) Run solute hbond analysis and view the strong and medium (and weak) hydrogen bond
networks in a pymol. Load this into your analysis pse to compare with RMSD
results.

5) Select residues from RMSD analysis or Hbond analysis or from project knowledge for clustering.
You can pass in residues from different chains, just make sure your numbering is
according to your reference PDB file passed in.

6) Run clustering on these residues with a few different cutoffs. You will get a
report of the % of the the simulation represented by the clusters, and many
times there will be a "cliff" in RMSD cutoff criteria after which you get
"reasonable" number of clusters. If you have two different simulation
conditions, the same cutoff can give different results, indicating a change in
the sampling. Load up the centroid-pdbs/*pdb into pymol to examine the
structures.

7)Then the AnalyzeMD solvent_calc can look at occupancy of waters and
analyze solvent hydrogen bonding. I recommend loading the water0.8.pdb into the
pse file you have going. The pml file produced will show you which residues have
a bridged water, connected by stron, medium and weak lines in PyMol.

