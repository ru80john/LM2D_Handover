import xml.etree.ElementTree as ET
import sys

calc_molecule = sys.argv[1]
catalyst_name = sys.argv[2]
edft_root = ET.parse('backup.edft')
idft_root = ET.parse('backup.idft')
new_root = ET.Element('jobs')
new_root_edft_CNP = ET.Element('jobs')
new_root_edft_calc_molecule = ET.Element('jobs')

edft_CNP = []
for job in idft_root.findall('job'):
    job_tag = job.find('tag').text
    segments = job.findall('input/segment')
    calc_molecule_inside = False
    for segment_element in segments:
        if segment_element.get('type') == calc_molecule:
            calc_molecule_inside = True
            break
    if calc_molecule_inside:
        new_root.append(job)
        for segment_element in segments:
            if segment_element.get('type') != calc_molecule:
                edft_CNP.append(segment_element.text)

for id in edft_CNP:
    for job in edft_root.findall('job'):
        segments = job.find('input/segment')
        if segments.get('id') == str(id):
            new_root_edft_CNP.append(job)

for job in edft_root.findall('job'):
    segments = job.find('input/segment')
    if segments.get('type') == calc_molecule:
        new_root_edft_calc_molecule.append(job)




tree = ET.ElementTree(new_root)
tree.write('idft.jobs', encoding='utf-8', xml_declaration=False)
tree = ET.ElementTree(new_root_edft_CNP)
tree.write(f'edft_{catalyst_name}.jobs', encoding='utf-8', xml_declaration=False)
tree = ET.ElementTree(new_root_edft_calc_molecule)
tree.write(f'edft_{calc_molecule}.jobs', encoding='utf-8', xml_declaration=False)

