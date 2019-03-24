# QMMM_Utils
### 1. amber_to_mcce
- Convert AMBER parameters to MCCE.
### 2. prep_oniom
- Prepare oniom input file from AMBER parm7 and rst7 files.
- Extract coordinates from gaussian optimization log in AMBER rst7 format. 
- Convert AMBER force fields to GAUSSIAN
- ff14SB contains atom types 2C and 3C. Gaussian does not recognize types which start with a number. Workaround: `preponiom` renames 2C, 3C with Z2, Z3.
- ff14 has C\* and N\* entries in VDW table. Gaussian does not recognize wildcards, so all types like CA, CB .. etc must be manuall added to gaussian prm file.
