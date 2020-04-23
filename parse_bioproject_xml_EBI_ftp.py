import pandas as pd 
import xml.etree.ElementTree as etree
from Bio import Entrez

# Set URL for XML file with all BioSamples associated with project
from urllib.request import urlopen
url = 'https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=biosample&HistoryId=NCID_1_125617756_130.14.18.48_5555_1587647699_834710852_0MetA0_S_HStore&QueryKey=8&Sort=&Filter=all&CompleteResultCount=683&Mode=file&View=fullxml&p$l=Email&portalSnapshot=%2Fprojects%2FBioSample%2Fbiosample%401.33&BaseUrl=&PortName=live&RootTag=BioSampleSet&FileName=&ContentType=xml'
xml_file = urlopen(url)

tree = etree.parse(xml_file)
root = tree.getroot()

# Append list of each all attributes for each sample within XML
all_samples = []
for child in root:
  samples = []
  for child in child:
    attribute_element = child.findall('Attribute')
    if attribute_element is not None:
      elements = []
      for element in attribute_element:
        individual = [element.items()[0][1], element.text]
        elements.append(individual)
      samples.append(elements)
  all_samples.append(samples)

# Replace all spaces with underscores to remove duplicates
for i in range(len(all_samples)):
  for num in range(len(all_samples[i][5])):
    all_samples[i][5][num][0] = all_samples[i][5][num][0].replace(' ', '_')

# Full list of Attributes across all sample 
features = []
for i in range(len(all_samples)):
  for x, y in all_samples[i][5]:
    features.append(x)
features = set(features)
len(features)

#Left join each sample by feature to create dataframe with all samples

df = pd.DataFrame(features, columns = ['features'])
df.set_index(['features'], inplace = True)
              
for i in range(len(all_samples)):
    temp = pd.DataFrame(all_samples[i][5], columns = ['features', str('Sample' + '_' + str(i))])
    temp.set_index(['features'], inplace = True)
    df = df.join(temp, how = 'left')
    if i % 100 == 0:
      print(str(i))

# Transpose so each row is a subject and each column is sample
df = df.T

# Reset dataframe index
df.reset_index(inplace = True)
df.rename(columns = {'index' : 'samples'}, inplace= True)

# Build link to project on EBI and download sample data

base = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport'
project = 'PRJEB18471'
result = 'read_run'
fields = 'study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp'
download = 'txt'

url = base + '?accession=' + project + '&result=' + result + '&fields=' + fields + '&download=' + download

download_links = pd.read_csv(url, sep = '\t')

# Merge data parsed from XML with ftp links for sequencing data from EBI
df = df.merge(download_links, how = 'left', left_on = 'SRA_accession', right_on = 'secondary_sample_accession')

