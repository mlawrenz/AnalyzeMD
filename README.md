Suggested initial workflow:

1) Run rmsd_calc using rmsd-all, load the pml file in rmsd-analysis to visualize
residues that deviate the most from starting position (red) color. Using
rmsf-all will show you residues that are most *flexible* throughout the
simulation (most variance from the average). The two metrics tell you different
things.
2) Select residues from RMSD analysis or from project knowledge for clustering.
You can pass in residues from different chains, just make sure your numbering is
according to your reference PDB file passed in.
3) Run clusering on these residues with a few different cutoffs. You will get a
report of the % of the the simulation represented by the clusters, and many
times there will be a "cliff" in RMSD cutoff criteria after which you get
"reasonable" number of clusters. If you have two different simulation
conditions, the same cutoff can give different results, indicating a change in
the sampling. Load up the centroid-pdbs/*pdb into pymol to examine the
structures.
