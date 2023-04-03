#!/usr/bin/env python
# coding: utf-8



import os, re
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from collections import Counter
import pandas as pd
import multiprocessing as mp
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",help = "Enter input file in space saperated format where first column represents fasta file name and second column contains comma separated mutation list")
    args = parser.parse_args()

# Open the FASTA file and read in the sequence
#filename = input("Enter the name of the input file: ")
filename = args.input_file
# filename = /home/fathima/Sample_code/input_file.txt"
inputfile = open(filename, "r").readlines()
filename1 = [i.strip().split(" ") for i in inputfile]


for lines in filename1:
    fastafile = lines[0]
    list_mut = lines[1].split(',')
    
motif_list = ['nM','Mc', 'n_M', 'M_c', 'n__M','M__c', 'tri']
netsurfp = pd.read_csv("./data/"+fastafile.split('.')[0]+"_netsurfp.csv")
physicochem_properties_actual = pd.read_csv("./data/49_properties_numerical_Values.csv")
aacon_header = ['mut_pos','KABAT', 'JORES', 'SCHNEIDER', 'SHENKIN', 'GERSTEIN',
                'TAYLOR_GAPS', 'TAYLOR_NO_GAPS', 'VELIBIL', 'KARLIN', 'ARMON',
                'THOMPSON', 'NOT_LANCET', 'MIRNY', 'WILLIAMSON', 'LANDGRAF', 
                'SANDER', 'VALDAR', 'SMERFS']



all_prop = pd.DataFrame()
# =============================================================================
# Sequence based properties calculated using Biopython module
record = SeqIO.read("./example/"+fastafile, "fasta")
sequence = str(record.seq)
# Calculate the sequence-based features
pp = ProtParam.ProteinAnalysis(sequence)
mol_weight = pp.molecular_weight()
aa_count = pp.count_amino_acids()
aa_percent = pp.get_amino_acids_percent()
charge = pp.charge_at_pH(7.4)
gravy = pp.gravy()
instability_index = pp.instability_index()
iso_point = pp.isoelectric_point()
secondary_struct = pp.secondary_structure_fraction()
list_prop_pp = [pp,mol_weight,aa_count, charge, gravy,
                instability_index, iso_point, secondary_struct]

# =============================================================================
# 
# =============================================================================

######## Parse the fasta file if not having BioPython  #################
def read_fasta(file, list_mut):
    """
    read fasta sequence
    :param file:
    :return:
    """
    msg = ''
    if not os.path.exists(file):
        msg = 'Error: file %s does not exist.' % file
        return [], None, msg
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].strip(), re.sub('[^ACDEFGHIKLMNPQRSTUVWY-]', '-', ''.join(array[1:]).upper())
        fasta_sequences.append([header, sequence, list_mut])
    return fasta_sequences
seq = read_fasta("./example/"+fastafile, list_mut)

######### Function to calcualate the amino acid composition ##########
def _AAC(fasta_list):
    try:
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        header = ['SampleName']
        encodings = []
        for i in AA:
            header.append('AAC_{0}'.format(i))
        encodings.append(header)
        for i in fasta_list:
            name, sequence = i[0], re.sub('-', '', i[1])
            count = Counter(sequence)
            for key in count:
                count[key] = count[key] / len(sequence)
            code = [name]
            for aa in AA:
                code.append(count[aa])
            encodings.append(code)
        encodings = np.array(encodings)
        encodings = pd.DataFrame(encodings[1:, 1:].astype(float),
                                 columns=encodings[0, 1:], index=encodings[1:, 0])
        return encodings
    except:
        error_msg = "There is an error with the input sequence"
        print(error_msg)
        return False
# =============================================================================
# Function for calculating mutation based properties
# =============================================================================

############  Preparation of a dataframe to obtain individual mutations in a sequence
list_names = ['name', 'sequence', 'mutations list']
df_input = pd.DataFrame(columns=list_names)
for i in seq:
    df_input.loc[len(df_input)] = i 
dataset = df_input.explode('mutations list').reset_index(drop=True)

# =============================================================================
# This set of snippet will calculate different physicochemical properties, 
# PSSM scores, and conservation scores calculated from different methods
# =============================================================================

def compute_features_mutation(i, all_prop = all_prop):
    uniprot_id = dataset['name'][i]
    site = dataset['mutations list'][i][:-1]
    sul_c = 0
    polar = 0
    aliphatic = 0
    neg_c = 0
    pos_c = 0
    arom = 0
    pos = int(site[1:])
    sequence = dataset['sequence'][i]
    if len(sequence) >= pos:
        dict_prop = dict()
        wild = site[0]
        mut = dataset['mutations list'][i][0]
        wild_mut = wild+mut
        dict_prop['pos'] = int(site[1:])
        dict_prop['Uniprot ID'] = uniprot_id
        dict_prop['Site'] = site  
        f4 = open("./data/"+ fastafile.split('.')[0]+".pssm", 'r').readlines()[2:]
        if pos<7:
            window_13 = sequence[:pos+6]
        elif pos+6>len(sequence):
            window_13=sequence[pos-7:]
        else:
            window_13=sequence[pos-7:pos+6]
        dict_prop['window_13'] = window_13
        
        sul_c = len(re.findall('[CM]', window_13))
        pos_c = len(re.findall('[KRH]', window_13))
        aliphatic = len(re.findall('[GALIV]', window_13))
        arom = len(re.findall('[YFW]', window_13))
        neg_c = len(re.findall('[DE]', window_13))
        polar = len(re.findall('[NQSTP]', window_13))
        if pos>=len(sequence):
            n_ter = '_'
            dipep_n = n_ter+wild[0]
            gap2_n = n_ter+wild[0]
        else:
            n_ter = sequence[pos-2]
            dipep_n = n_ter+wild[0]
            gap2_n = sequence[pos-4]+wild[0]
        
        if pos>=len(sequence):
            n_ter_gap = '_'
            dipep_gap_n = n_ter_gap+wild[0]
        else:
            n_ter_gap = sequence[pos-3]
            dipep_gap_n = n_ter_gap+wild[0]
        
        if pos+1 >= len(sequence):
                c_ter_gap = '-'
                dipep_gap_c = wild[0] + c_ter_gap
                gap2_c = wild[0] + c_ter_gap
        else:
            c_ter_gap = sequence[pos+1]
            dipep_gap_c = wild[0] + c_ter_gap
            try:
                gap2_c = wild[0] + sequence[pos+2]
            except IndexError:
                gap2_c = wild[0] +c_ter_gap
        if pos >= len(sequence):
            c_ter = '-'
            dipep_c = wild[0]+c_ter
        else:
            c_ter = sequence[pos]
            dipep_c = wild[0]+c_ter
        tripep = n_ter+ wild[0]+c_ter
        dict_prop['nM'] = dipep_n
        dict_prop['Mc']= dipep_c
        dict_prop['tri'] = tripep
        dict_prop['n_M']=dipep_gap_n
        dict_prop['M_c'] = dipep_gap_c
        dict_prop['n__M'] = gap2_n
        dict_prop['M__c'] = gap2_c
        b = sum(physicochem_properties_actual[window_13[k]] for k in range(len(window_13)))/sum(c.isalpha() for c in window_13)
        j = 0
        for property1 in physicochem_properties_actual['index']:
            dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] 
                                                                     == property1][wild].tolist()[0]  
            if property1 == "pK'":
                dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index']
                                                                         == "pK'"][wild].tolist()[0]
                dict_prop[property1] =  b[j]
            else:
                dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index']
                                                                         == property1][wild].tolist()[0]
                dict_prop[property1] =  b[j]
                dict_prop[property1+'_diff'] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] 
                                                                                 == property1][wild_mut[-1]].tolist()[0] - physicochem_properties_actual.loc[physicochem_properties_actual['index'] == property1][wild_mut[0]].tolist()[0]
            j += 1
        dict_prop['neg_charge'] = neg_c
        dict_prop['pos_charge'] = pos_c
        dict_prop['polar'] = polar
        dict_prop['aromatic'] = arom
        dict_prop['S_containing'] = sul_c
        dict_prop['aliphatic'] = aliphatic
        
        try:
            dict1 = {'A':f4[pos].strip().split()[2:][0],'R':f4[pos].strip().split()[2:][1],
                'N':f4[pos].strip().split()[2:][2], 'D':f4[pos].strip().split()[2:][3],
                      'C':f4[pos].strip().split()[2:][4], 'Q':f4[pos].strip().split()[2:][5], 
                  'E':f4[pos].strip().split()[2:][6], 'G':f4[pos].strip().split()[2:][7],
                      'H':f4[pos].strip().split()[2:][8], 'I':f4[pos].strip().split()[2:][9], 
                  'L':f4[pos].strip().split()[2:][10], 'K':f4[pos].strip().split()[2:][11],
                      'M':f4[pos].strip().split()[2:][12], 'F':f4[pos].strip().split()[2:][13], 
                  'P':f4[pos].strip().split()[2:][14], 'S':f4[pos].strip().split()[2:][15],
                      'T':f4[pos].strip().split()[2:][16], 'W':f4[pos].strip().split()[2:][17],
                  'Y':f4[pos].strip().split()[2:][18], 'V':f4[pos].strip().split()[2:][19]}
            dict_prop['pssm_score1'] =  int(dict1[site[0]])   
            dict_prop['pssm_score2'] = sum([int(i) for i in dict1.values()])/20
            dict_prop['pssm_score3'] = int(dict1[wild_mut[-1]])-int(dict1[wild_mut[0]])  
            cons= f4[pos].strip('\n')[90:].split()[21:22]
            if len(cons) ==0:
                dict_prop['conservation']= 0
            else:
                for element in cons:
                    dict_prop['conservation']= float(element)
        except (KeyError, IndexError) as error:
                dict_prop['pssm_score1'] =  0
                dict_prop['pssm_score2'] = 0
                dict_prop['conservation']= 0
        f_aacon = pd.read_csv("./data/"+ fastafile.split('.')[0] + ".features", sep = '\t', skiprows=1, 
                              names = aacon_header)
        try:
            mm = f_aacon.iloc[pos-1]
            for item in aacon_header:
                dict_prop[item] = mm[item]
        except IndexError:
            for item in aacon_header:
                dict_prop[item] = 0
        all_prop = all_prop.append(dict_prop, ignore_index =True)
        return all_prop

df11 = pd.DataFrame()

def mp_func(i):
    re1 = compute_features_mutation(i)
    return re1
    
pool = mp.Pool(mp.cpu_count()-2)
re11 = pool.map(mp_func,range(len(dataset)))
df11 = df11.append(re11, ignore_index=True)


#####################  Save output for mutation bsed properties ##########################

df11.to_csv('./output/output_mutation_based_features.csv', index = False)

#####################  Save output for mutation bsed properties ##########################

info = _AAC(seq)
dict_prop_pp = {'mol_weight': mol_weight,'aa_count': aa_count, 'charge' : charge, 
                'gravy': gravy, 'instability_index': instability_index, 'iso_point': iso_point, 
                'sec_struc': secondary_struct}
for item in dict_prop_pp.keys():
    if item == 'aa_count':
        for aa_res in dict_prop_pp[item].keys():
            info[aa_res] = dict_prop_pp[item][aa_res]
    elif item == 'sec_struc':
        info[item] = str(dict_prop_pp[item])
    else:
        info[item] = dict_prop_pp[item]
info.to_csv('./output/output_protein_based_features.csv', index = False)

