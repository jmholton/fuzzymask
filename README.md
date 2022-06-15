# fuzzymask
fuzzy-logic like bulk solvent mask generator

Normal bulk solvent masks treat the atoms in a protein as hard spheres.
Even atoms with occupancy of 0.01 still exclude just as much solvent as 
atoms that have an occupancy of 1.00. Atoms with big B factors and are
flapping around in a large volume still have the same "radius" as atoms
that are well-ordered.
The fuzzymask approach is to treat the provided PDB file as a probability 
distribution. By default, 500 new PDBs are generated from this distribution.
B factors become displacement distributions and atoms are randomly "kicked"
in random directions.  This breaks bonds, of course, so a follow-up step 
averages the shift vectors of neighboring bonded atoms so that local groups
of atoms tend to move in concert. The shift magnitudes are suppressed by 
this thru-bond averaging and so are then scaled up and the thru-bond averaging
repeated until the rms shifts are equal to those prescribed by the B factors.
Partial occupancy atoms are randomly deleted with a probability consistent
with their occupancy. By default, all atoms in the same residue blink on
and off together.
Finally, these 500 "jiggled" molecules are run through a standard bulk-solvent
calculation using refmac5. This will run in parallel on SGE, PQS or Slurm clusters.
These maps are then averaged together to form the "fuzzy" mask.  Its values will
range from 0 to 1 with 1 meaning 100% occupied by bulk solvent, and 0 at the center
of well-ordered protein atoms. The fuzzy mask is then converted into phased structure
factors and appended to the mtz file provided, and is ready to be input into a
standard refmac5 run with the original PDB file.  Refmac supports externally-supplied
solvent like this by using the FPART keywords on the LABIN instruction. Don't forget
to turn off the default bulk solvent calculation with "SOLVENT NO" in your refmac run.
Results with these "fuzzy" masks tend to be superior to the default bulk solvent mask,
and the fuzzy mask can be re-used for many subsequent runs. It is a good thing to try
and the end of your refinement to see how much further you can push your R factors down,
and the good news is: no new free parameters!  Fuzzy masks are derived from nothing but
the input coordinates, and so can't contribute to "over fitting" any more than the
default bulk solvent mask does.

