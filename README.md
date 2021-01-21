# t2diaphite-sfexpt-backup
- sequences.txt contains the code sequences for different diaphite structures

- scalefactors.txt contains the scale factors to sample

- sfexpt.sh accepts lcbop or tersoff as arguments and generates the input positions (crystal structure multiplied by 4 in x&y and 3 in z), runs a LAMMPS minimisation on each and returns a data table of scalefactors vs energy/atom and volume/atom

- read_data.sh also accepts lcbop or tersoff and compiles a table of minimum energy/atom and volume/atom as a function of structure
