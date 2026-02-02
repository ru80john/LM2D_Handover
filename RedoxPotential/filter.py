import xml.etree.ElementTree as ET
import sys

calc_molecule = sys.argv[1]
root = ET.parse('backup.xml')
new_root = ET.Element('jobs')

for job in root.findall('job'):
    job_tag = job.find('tag').text
    if (calc_molecule in job_tag) and (('n' in job_tag) or ('e' in job_tag)):
        new_root.append(job)
    else:
        pass

tree = ET.ElementTree(new_root)
tree.write('jobwriter.mps.chrg.xml', encoding='utf-8', xml_declaration=False)

