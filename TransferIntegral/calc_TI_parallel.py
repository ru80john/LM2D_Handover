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
catalyst_name = 'CNP'      # name of the catalyst
calc_molecule = 'AA'      # molecule type to be (de‑)protonated
radical_name  = 'HA_rad'  # name for the radical species
calc_idx      = [num for num in range(136)]       # list of indices (starting from 0) to process
#calc_idx = [0,1,2,3]

# Each mdrun will spawn this many OpenMP threads.
#   NTOMP  ×  MAX_WORKERS  ≦  available CPU cores (48 in your case)
NTOMP        = 12
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
initial_dir = Path(os.path.dirname(os.path.abspath(__file__)))

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

# --------------------------------------------------
# the heavy job: executed in a *separate* process for each idx
# --------------------------------------------------

def run(Results_dir: int, idx: str):
    # ------------------ working directory ----------------
    work_dir = Path(os.path.join(Results_dir, f'{calc_molecule}_{idx}'))
    # move everything to the DOS directory
    DOS_dir = work_dir / 'DOS'
    DOS_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'mv {work_dir}/* {DOS_dir}')
    # Create TI
    TI_dir = work_dir / 'TI'
    shutil.rmtree(TI_dir, ignore_errors=True)
    TI_dir.mkdir(parents=True, exist_ok=True)
    # Extract molecules
    subprocess.run(['python', f'{initial_dir}/extract_molecule.py', catalyst_name, radical_name, DOS_dir / 'mapping' / 'map.tpr', DOS_dir / 'mapping' / 'map.gro', TI_dir / 'extract.gro'], cwd=initial_dir, check=True)
    # Extract top information
    SYSTEM_TOP = DOS_dir / 'mapping' / f'{calc_molecule}_{idx}.top'
    defaults_sec   = extract_defaults_section(SYSTEM_TOP)
    mol_list       = extract_molecules_section(SYSTEM_TOP)
    new_mol_list   = [molecule for molecule in mol_list if molecule[0] == radical_name or molecule[0] == catalyst_name]
    ## ------------ generate new .gro / .top -------------
    top_txt = defaults_sec + '\n' + create_system_and_molecules(f'{calc_molecule}_{idx}', new_mol_list)
    write_file(TI_dir / 'extract.top', top_txt)
    # Copy files
    shutil.copytree(DOS_dir / 'mapping' / 'system.ff', TI_dir / 'system.ff', dirs_exist_ok=True)
    shutil.copy(initial_dir / 'map.mdp', TI_dir / 'map.mdp')
    # ------------------- grompp ------------------------
    subprocess.run([
        'gmx_mpi', 'grompp',
        '-f', 'map.mdp',
        '-c', 'extract.gro',
        '-p', 'extract.top',
        '-o', 'map.tpr',
        '-maxwarn', '5'
    ], cwd=TI_dir, check=True)
#
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
    subprocess.run(mdrun_cmd, cwd=TI_dir, check=True)

    # mapping directory
    mapping = TI_dir / 'mapping'
    mapping.mkdir(exist_ok=True)
    os.system(f'mv {TI_dir}/* {mapping}')
    shutil.rmtree(TI_dir / 'sysyem.ff', ignore_errors=True)
    shutil.copytree(DOS_dir / 'mapping' / 'QC_FILES', mapping / 'QC_FILES', dirs_exist_ok=True)
    shutil.copy(DOS_dir / 'system.xml', mapping / 'system.xml')
    pbc_whole(mapping / 'map.tpr', mapping / 'map.gro', mapping)
#
    subprocess.run(['ctp_map', '-t', 'map.tpr', '-c', 'whole.gro', '-s', 'system.xml', '-f', 'state.sql'], cwd=mapping, check=True)
#
    shutil.copy(mapping / 'state.sql', TI_dir / 'state.sql')
    shutil.copy(mapping / 'system.xml', TI_dir / 'system.xml')
    #shutil.copy(initial_dir / 'options.xml', TI_dir / 'options.xml')
    os.system(f'cp {initial_dir}/package* {TI_dir}')
    subprocess.run(['ctp_run', '-o', initial_dir / 'options.xml', '-e', 'neighborlist','-f', 'state.sql'], cwd=TI_dir, check=True)
    subprocess.run(['ctp_parallel', '-o', initial_dir / f'option_{catalyst_name}.xml', '-e', 'edft', '-j', 'write', '-f', 'state.sql'], cwd=TI_dir, check=True)
    subprocess.run(['ctp_parallel', '-o', initial_dir / f'options.xml', '-e', 'idft', '-j', 'write', '-f', 'state.sql'], cwd=TI_dir, check=True)
    os.system(f'mv {TI_dir}/edft_{catalyst_name}.jobs {TI_dir}/backup.edft')
    os.system(f'mv {TI_dir}/idft.jobs {TI_dir}/backup.idft')
    subprocess.run(['python', f'{initial_dir}/filter.py', f'{radical_name}', f'{catalyst_name}'], cwd=TI_dir, check=True)
    # Submit the job
    subprocess.run(['ctp_parallel', '-o', initial_dir / f'option_{catalyst_name}.xml' ,'-f' ,'state.sql', '-e', 'edft', '-j', 'run', '-t', '1', '-c', '1', '-m', '-1'], cwd=TI_dir, check=True)
    subprocess.run(['ctp_parallel', '-o', initial_dir / f'option_{radical_name}.xml' ,'-f' ,'state.sql', '-e', 'edft', '-j', 'run', '-t', '1', '-c', '1', '-m', '-1'], cwd=TI_dir, check=True)
    subprocess.run(['ctp_parallel', '-o', initial_dir / 'options.xml' ,'-f' ,'state.sql', '-e', 'idft', '-j', 'run', '-t', '1', '-c', '1', '-m', '-1'], cwd=TI_dir, check=True)

# --------------------------------------------------
# main driver
# --------------------------------------------------

def main():
    current_directory = os.getcwd()
    Results_dir = os.path.join(current_directory, f'Results_{calc_molecule}')
    #for idx in calc_idx:
    #    run(Results_dir, idx)

    ctx = mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=MAX_WORKERS, mp_context=ctx) as pool:
        futures = [pool.submit(run, Results_dir, idx) for idx in calc_idx]
        for fut in as_completed(futures):
            print(fut.result())
if __name__ == '__main__':
    main()
