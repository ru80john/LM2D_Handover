import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

current_directory = os.getcwd()
Results_dir = os.path.join(current_directory, 'Results_AA')
mol_1, mol_2 = 'CNP', 'AA'

def get_TI_from_xml(xml_file):
    """
    Extract energy from the given XML file.
    """
    energy_info = []
    try:
        tree = ET.parse(xml_file)
        for job in tree.findall('job'):
            id = job.findall('input/segment')[0].text
            job_TI = job.findall('output/pair/overlap')
            jAB = float(job_TI[2].get('jAB'))
            energy_info.append((id, jAB**2))
        return energy_info
    except ET.ParseError as e:
        print(f"Error parsing {xml_file}: {e}")
        return None
    
sort_dir = sorted(os.listdir(Results_dir), key = lambda x: int(x.split('_')[1]))
info = {'Molecule A':[], 'Molecule B':[], 'TI (eV^2)':[]}
for entry in sort_dir:
    # Construct the full path
    full_path = os.path.join(Results_dir, entry)
    # Check if the entry is a directory
    if os.path.isdir(full_path):
        energy_info = get_TI_from_xml(os.path.join(full_path, 'TI', 'idft.jobs'))
        for id, TI in energy_info:
            info['Molecule A'].append(f'{mol_1}_{int(id) - 1}')
            info['Molecule B'].append(entry)
            info['TI (eV^2)'].append(TI)

# Create a DataFrame from the collected information
df = pd.DataFrame(info)
# Save the DataFrame to a CSV file
df.to_csv('TI_HA.csv', index=False)
