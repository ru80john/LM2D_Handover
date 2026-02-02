import os
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import MDAnalysis as mda

# ---------------------------
# user‑configurable constants
# ---------------------------
calc_molecule = 'AA'      # molecule type to be (de‑)protonated
radical_name  = 'HA_rad'  # name for the radical species
calc_idx      = [num for num in range(136)]       # list of indices (starting from 0) to process
delete_atom_index = 11  # index of the atom to be deleted (starting from 1)

# Each mdrun will spawn this many OpenMP threads.
#   NTOMP  ×  MAX_WORKERS  ≦  available CPU cores (48 in your case)
NTOMP        = 2
MAX_WORKERS  = 4

# tables for dihedral interactions (leave '' if none)
table_files = (
    'system.ff/table_d0.xvg  system.ff/table_d1.xvg  system.ff/table_d2.xvg  '
    'system.ff/table_d3.xvg  system.ff/table_d4.xvg  system.ff/table_d5.xvg  '
    'system.ff/table_d6.xvg  system.ff/table_d7.xvg  system.ff/table_d8.xvg  '
    'system.ff/table_d9.xvg  system.ff/table_d10.xvg'
)
table_files = [string for string in table_files.split() if string]

# paths for static input files
MD_DIR      = Path('MD_FILES')
QC_DIR      = Path('QC_FILES')
MP_DIR      = Path('MP_FILES')
SYSTEM_TOP  = MD_DIR / 'system.top'
SYSTEM_TPR  = MD_DIR / 'system.tpr'
SYSTEM_GRO  = MD_DIR / 'system.gro'
MAP_MDP     = Path('map.mdp')
OPTION_XML  = Path('options.xml')

# --------------------------------------------------
# helper functions (mostly copied from original script)
# --------------------------------------------------

def extract_defaults_section(input_file: Path) -> str:
    defaults = []
    grab = False
    with input_file.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('[ defaults ]'):
                grab = True
            if line.startswith('[ system ]'):
                break
            if grab:
                defaults.append(line)
    return '\n'.join(defaults)

def extract_molecules_section(file_path: Path):
    molecules = []
    grab = False
    with file_path.open() as fh:
        for line in fh:
            line = line.strip()
            if '[ molecules ]' in line:
                grab = True
                continue
            if grab and line and not line.startswith(';'):
                mol, n = line.split()[:2]
                molecules.append((mol, int(n)))
    return molecules

def create_system_and_molecules(system_name, molecules_list, charge=None):
    """
    Create a GROMACS .top file with [ system ] and [ molecules ] sections.
    Inputs:
    - system_name: Name of the system (string)
    - molecules_list: List of tuples (molecule_name, number_of_molecules)
    """
    # Create [ system ] section
    system_section = f"[ system ]\n;   name\n  {system_name}\n"

    # Create [ molecules ] section
    molecules_section = "[ molecules ]\n; compound    n_mol\n"
    for molecule, count in molecules_list:
        if molecule.split('_')[-1] == 'charge':
            molecule = molecule.split('_')[0] + '_' + charge
        molecules_section += f"  {molecule:<10}{count}\n"

    # Combine sections
    result = system_section + "\n" + molecules_section
    return result

def write_file(output, *arg):
    line = ''
    for i in arg:
        line += i
    with open(output, 'w') as f:
        f.write(line)

def pbc_whole(tpr_file, gro_file, output_dir):
    import MDAnalysis as mda
    import MDAnalysis.transformations as trans
    traj = mda.Universe(tpr_file, gro_file)
    transforms = [
      trans.center_in_box(traj.atoms),
      trans.wrap(traj.atoms),
      trans.unwrap(traj.atoms)]
    traj.trajectory.add_transformations(*transforms)
    traj.atoms.write(f'{output_dir}/whole.gro')

def gen_gro_and_top_files(traj, molecule_list, calc_molgroup, calc_idx):
    universe_atom = traj.atoms[[]]  # Empty atom group
    new_molecule_list = []
    for mol_name, number in molecule_list:
        if mol_name == calc_molecule:
            # Create a copy of calc_molgroup to avoid modifying the original
            modified_molgroup = [mol for mol in calc_molgroup]
            
            # Deprotonate the molecule at calc_idx (remove specific H atom)
            deprotonated_mol = calc_molgroup[calc_idx] - calc_molgroup[calc_idx][delete_atom_index - 1]
            modified_molgroup[calc_idx] = deprotonated_mol
            
            # Add all molecules to the universe
            for mol in modified_molgroup:
                universe_atom += mol
                
            # Update molecule_list according to which molecule was deprotonated
            if calc_idx == 0:
                new_molecule_list.append((radical_name, 1))
                new_molecule_list.append((calc_molecule, number - 1))
            elif calc_idx == len(calc_molgroup) - 1:
                new_molecule_list.append((calc_molecule, number - 1))
                new_molecule_list.append((radical_name, 1))
            else:
                new_molecule_list.append((calc_molecule, calc_idx))
                new_molecule_list.append((radical_name, 1))
                new_molecule_list.append((calc_molecule, number - calc_idx - 1))
        else:
            # Add other molecule types directly
            universe_atom += traj.select_atoms(f"moltype {mol_name}")
            new_molecule_list.append((mol_name, number))
    
    return universe_atom, new_molecule_list

# --------------------------------------------------
# the heavy job: executed in a *separate* process for each idx
# --------------------------------------------------

def run_single_idx(idx: int, initial_dir: str):
    # Re‑parse topology & trajectory inside the child process to avoid pickling big objects
    traj = mda.Universe(str(SYSTEM_TPR), str(SYSTEM_GRO))
    calc_atoms     = traj.select_atoms(f'moltype {calc_molecule}')
    calc_molgroup  = calc_atoms.split('molecule')
    defaults_sec   = extract_defaults_section(SYSTEM_TOP)
    mol_list       = extract_molecules_section(SYSTEM_TOP)

    # ---------------- working directory ----------------
    dir_path = Path(initial_dir) / f'Results_{calc_molecule}' / f'{calc_molecule}_{idx}'
    shutil.rmtree(dir_path, ignore_errors=True)
    dir_path.mkdir(parents=True, exist_ok=True)

    # ------------ generate new .gro / .top -------------
    new_atoms, new_mol_list = gen_gro_and_top_files(traj, mol_list, calc_molgroup, idx)
    new_atoms.write(dir_path / f'{calc_molecule}_{idx}.gro')
    top_txt = defaults_sec + '\n' + create_system_and_molecules(f'{calc_molecule}_{idx}', new_mol_list)
    write_file(dir_path / f'{calc_molecule}_{idx}.top', top_txt)

    # -------------- copy static inputs -----------------
    shutil.copytree(MD_DIR / 'system.ff', dir_path / 'system.ff', dirs_exist_ok=True)
    shutil.copy(MAP_MDP, dir_path / 'map.mdp')

    # ------------------- grompp ------------------------
    subprocess.run([
        'gmx_mpi', 'grompp',
        '-f', 'map.mdp',
        '-c', f'{calc_molecule}_{idx}.gro',
        '-p', f'{calc_molecule}_{idx}.top',
        '-o', 'map.tpr',
        '-maxwarn', '5'
    ], cwd=dir_path, check=True)

    # ------------------- mdrun -------------------------
    mdrun_cmd = [
        'gmx_mpi', 'mdrun',
        '-ntomp', str(NTOMP),
        '-v', 'yes',
        '-s', 'map.tpr',
        '-deffnm', 'map'
    ]
    if table_files:
        mdrun_cmd.append('-tableb')
        mdrun_cmd.extend(table_files)
    subprocess.run(mdrun_cmd, cwd=dir_path, check=True)

    # -------------- post‑processing --------------------
    pbc_whole(dir_path / 'map.tpr', dir_path / 'map.gro', dir_path)

    # mapping directory
    mapping = dir_path / 'mapping'
    mapping.mkdir(exist_ok=True)
    os.system(f'mv {dir_path}/* {mapping}')
    shutil.rmtree(f'{dir_path}/sysyem.ff', ignore_errors=True)

    shutil.copytree(QC_DIR, mapping / 'QC_FILES', dirs_exist_ok=True)
    shutil.copy('system.xml', mapping / 'system.xml')

    subprocess.run(['ctp_map', '-t', 'map.tpr', '-c', 'whole.gro', '-s', 'system.xml', '-f', 'state.sql'], cwd=mapping, check=True)

    # back to parent dir for ctp_run & filtering
    #if os.path.exists(dir_path / 'state.sql'):
    #    os.remove(dir_path / 'state.sql')
    shutil.copy(mapping / 'state.sql', dir_path / 'state.sql')
    shutil.copytree(initial_dir/MP_DIR, dir_path / 'MP_FILES', dirs_exist_ok=True)
    subprocess.run(['ctp_run', '-v', '-e', 'jobwriter', '-o', f'{initial_dir}/{str(OPTION_XML)}', '-f', 'state.sql'], cwd=dir_path, check=True)

    shutil.copy(dir_path / 'jobwriter.mps.chrg.xml', dir_path / 'backup.xml')
    subprocess.run(['python', f'{initial_dir}/filter.py', radical_name], cwd=dir_path, check=True)
    # Run VOTCA
    shutil.copy(mapping / 'system.xml', dir_path / 'system.xml')
    subprocess.run(['ctp_parallel', '-o', f'{initial_dir}/{str(OPTION_XML)}', '-f', 'state.sql', '-e', 'xqmultipole', '-t', str(NTOMP)], cwd=dir_path, check=True)

    return f'{calc_molecule}_{idx} finished'

# --------------------------------------------------
# main driver
# --------------------------------------------------

def main():
    initial_dir = os.getcwd()
    (Path(initial_dir) / f'Results_{calc_molecule}').mkdir(exist_ok=True)

    # "spawn" 避免 MPI + fork 混用的潛在問題
    ctx = mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=MAX_WORKERS, mp_context=ctx) as pool:
        futures = [pool.submit(run_single_idx, idx, initial_dir) for idx in calc_idx]
        for fut in as_completed(futures):
            print(fut.result())

if __name__ == '__main__':
    main()
