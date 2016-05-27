Suggested initial workflow:

1) Run rmsd_calc using rmsd-all, load the pml file in rmsd-analysis to visualize
residues that deviate the most from starting position (red) color. Using
rmsf-all will show you residues that are most *flexible* throughout the
simulation (most variance from the average). The two metrics tell you different
things.
2) Run solute hbond analysis and view the strong and medium hydrogen bond
networks in a pymol. Can be useful to also load in the loadbfactor.pml file from
rmsd analysis to see what residues move. And save this as an initialized
analysis pse file.
3) Select residues from RMSD analysis or from project knowledge for clustering.
You can pass in residues from different chains, just make sure your numbering is
according to your reference PDB file passed in.
4) Run clusering on these residues with a few different cutoffs. You will get a
report of the % of the the simulation represented by the clusters, and many
times there will be a "cliff" in RMSD cutoff criteria after which you get
"reasonable" number of clusters. If you have two different simulation
conditions, the same cutoff can give different results, indicating a change in
the sampling. Load up the centroid-pdbs/*pdb into pymol to examine the
structures.
5) To analyze solvent you must run Extract MD again with the --solvent flag.
This will be faster to produce the dcd file but may contain jumps in the
protein. Then the AnalyzeMD solvent_calc can look at occupancy of waters and
analyze solvent hydrogen bonding. I recommend loading the water0.8.pdb into the
pse file you have going. Then you can also map out the solute atoms with solvent
hbonds, bc you will just get a general value for these.
