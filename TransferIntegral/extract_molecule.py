import MDAnalysis as mda
import sys

catalyst_name = sys.argv[1]
radical_name  = sys.argv[2]
tpr = sys.argv[3]
gro = sys.argv[4]
output = sys.argv[5]

traj = mda.Universe(tpr, gro)
atoms = traj.select_atoms(f"moltype {catalyst_name} or moltype {radical_name}")

if __name__ == "__main__":
    atoms.write(output)
