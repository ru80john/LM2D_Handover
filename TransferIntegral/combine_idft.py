import os
import xml.etree.ElementTree as ET

current_directory = os.getcwd()
Results_dir = os.path.join(current_directory, 'Results_AA')
new_root = ET.Element('jobs')
start = 3888

sort_dir = sorted(os.listdir(Results_dir), key = lambda x: int(x.split('_')[1]))
n = 0
for entry in sort_dir:
    # Construct the full path
    full_path = os.path.join(Results_dir, entry)
    mol_idx = int(entry.split('_')[1])
    # Check if the entry is a directory
    if os.path.isdir(full_path):
        root = ET.parse(os.path.join(full_path, 'TI', 'idft.jobs'))
        n += len(root.findall('job'))
        for job in root.findall('job'):
            ###
            segment2 = job.findall('input/segment')[1]
            segment2.text = f'{mol_idx + 1 + start}'
            segment2.set('id', f'{mol_idx + 1 + start}')
            segment2.set('type', f'{entry.split("_")[0]}')
            ###
            pair = job.find('output/pair')
            pair.set('idB', f'{mol_idx + 1 + start}')
            pair.set('typeB', f'{entry.split("_")[0]}')
            #pair.set('homoA', f"{int(pair.get('homoA')) - 1}")
            pair.set('homoA', f"{int(pair.get('homoA'))}")
            pair.set('homoB', f"{int(pair.get('homoB')) - 1}")
            new_root.append(job)

print(f'Total number of jobs: {n}')
tree = ET.ElementTree(new_root)
tree.write('idft_HA.jobs', encoding='utf-8', xml_declaration=False)

