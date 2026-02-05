import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

current_directory = os.getcwd()
Results_dir = os.path.join(current_directory, 'Results_AA')
start = 3888

def get_energy_from_xml(xml_file):
    """
    Extract energy from the given XML file.
    """
    energy_info = {}
    try:
        tree = ET.parse(xml_file)
        for job in tree.findall('job'):
            job_tag = job.find('tag').text
            job_energy = job.find('output/summary/total').text
            energy_info[job_tag] = job_energy
        return energy_info
    except ET.ParseError as e:
        print(f"Error parsing {xml_file}: {e}")
        return None
# List all entries in the current directory
sort_dir = sorted(os.listdir(Results_dir), key = lambda x: int(x.split('_')[1]))
n = []
e = []
molecule = []
for entry in sort_dir:
    # Construct the full path
    full_path = os.path.join(Results_dir, entry)
    # Check if the entry is a directory
    if os.path.isdir(full_path):
        #print(f"Processing directory: {entry}")
        energy_info = get_energy_from_xml(os.path.join(full_path, 'jobwriter.mps.chrg.xml'))
        molecule.append(entry)
    for key, value in energy_info.items():
        if key[-1] == 'n':
            n.append(float(value))
        else:
            e.append(float(value))

for mol, neu, elec in zip(molecule, n, e):
    idx = int(mol.split('_')[1]) + start + 1
    molname = mol.split('_')[0]
    E_neu = neu - neu  # Neutral energy relative to itself
    E_elec = elec - neu  # Electron energy relative to neutral
    E_hole = 0  # Hole energy is set to zero as per the original code
    # Prepare the line to write
    line = f'{idx:4d}{molname:>11}{0:3d}{E_neu:>10.5f}{-1:4d}{E_elec:>10.5f}{1:3d}{E_hole:>10.5f}{1:3d}{32.0:>6.1f}{-1:4d}{32.0:>6.1f}{1:3d}{32.0:>6.1f}{"XXX 000":>9}\n'
    # Write the line to the output file
    with open(f'{molname}_E.txt', 'a') as f:
        f.write(line)