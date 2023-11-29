#########################################################################################################################
# This program is to extract ChIP-seq data for HCT116 biosamples from ENCODE, and extract the accessions of raw         #
# sequencing files that were concated by '|'. Then indentify the matching control experiments for each ChIP-seq         #
# experiment, and extract the accessions of the mathcing control raw sequencing files, concated by '|'. And save the    # 
# results in two files:                                                                                                 #
#       HCT116 _ChIPSeq_ENCODE_results_v1.txt                                                                           #
#       HCT116 _ChIPSeq_ENCODE_results_v2.txt                                                                           #
# the v1 result file contains the raw sequencing file accessions for the experiment and the control, separated by a tab.#
# the v2 result file contains the experiment accession, the raw sequencing file accessions for the experiment, the      # 
# control accession, and the raw sequencing file accessions for the control, separated by a tab.                        #
# 06/23/2021 Meng                                                                                                       #
#########################################################################################################################

import requests, json  

# return from the server in JSON format
headers = {'Accept': 'application/json'}

# Query for ChIP-seq data for HCT116 biosamples 
URL = ("https://www.encodeproject.org/search/?searchTerm=HCT116&type=Experiment&biosample_ontology.term_name=HCT116&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&assay_title=Mint-ChIP-seq&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&assay_title=Control+Mint-ChIP-seq&limit=200&control_type!=*")

response = requests.get(URL, headers=headers)

search_results = response.json()['@graph']

output = open('HCT116 _ChIPSeq_ENCODE_results_v1.txt', 'w')
output2 = open('HCT116 _ChIPSeq_ENCODE_results_v2.txt', 'w')
output.write('Experiment_raw_sequencing_files\tControl_raw_sequencing_files\n')
output2.write('Experiment_accession\tExperiment_raw_sequencing_files\tControl_accession\tControl_raw_sequencing_files\n')

# process each experiment
for experiment in search_results:    
    # extract the experiment accession 
    experiment_accession = experiment['accession']

    # extract experiment details
    URL = ("https://www.encodeproject.org/experiments/" + experiment_accession + "/?format=json")
    response = requests.get(URL, headers=headers)
    experiment_details = response.json()

    # extract the raw sequencing file accessions for each experiment
    raw_file_accessions = [file['accession'] for file in experiment_details['files'] if file['output_type'] == 'reads']
    experiment_raw_sequencing_files = '|'.join(raw_file_accessions)
    
    # extract the matching control experiment accession 
    control_accession = experiment_details['possible_controls'][0]['accession']

    # extract raw sequencing file accessions for the control
    URL = ("https://www.encodeproject.org/experiments/" + control_accession + "/?format=json")
    response = requests.get(URL, headers=headers)
    control_details = response.json()
    control_raw_file_accessions = [file['accession'] for file in control_details['files'] if file['output_type'] == 'reads']
    control_raw_sequencing_files = '|'.join(control_raw_file_accessions)
    
    # write the results to the output v1 file
    output.write(experiment_raw_sequencing_files + '\t' + control_raw_sequencing_files + '\n')
    # write the results to the output v2 file
    output2.write(experiment_accession + '\t' + experiment_raw_sequencing_files + '\t' + control_accession + '\t' + control_raw_sequencing_files + '\n')

output.close()
output2.close()





    


