import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd


current_directory = os.getcwd()
Results_dir = os.path.join(current_directory, 'Results_AA')

# vertical
#gas_HA_rad = -684.379910046
#gas_HA_dep = -684.507459972
#gas_HA_rad = -684.392474607
#gas_HA_dep = -684.494207472

# adiabatic
gas_HA_rad = -684.1146498
gas_HA_dep = -684.2291761


hat_to_eV = 27.2107
EA_gas = (gas_HA_rad - gas_HA_dep) * hat_to_eV

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

n = np.array(n)
e = np.array(e)
EA_solution = n - e
EA_final = -(EA_solution + EA_gas)
summary = {'Molecule':molecule, 'EA': EA_final}
df = pd.DataFrame(summary)
df.to_csv('AA_EA.csv', index = False)



