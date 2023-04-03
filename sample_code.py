# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:44:41 2023

@author: HP
"""
import os, re
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from collections import Counter
import pandas as pd
import multiprocessing as mp
# Open the FASTA file and read in the sequence
#filename = input("Enter the name of the FASTA file: ")
filename = "P04637.txt"
# mutation = input("List the mutations separated by comma: ")
mutation = 'P8V'
list_mut  = mutation.split()
all_prop = pd.DataFrame()
# =============================================================================
# Sequence based properties calculated using Biopython module
record = SeqIO.read(filename, "fasta")
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
    # for mut in list_mut:
        # fasta_sequences.append()
    return fasta_sequences
seq = read_fasta("P04637.txt", ['P8V'])
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
        encodings = pd.DataFrame(encodings[1:, 1:].astype(float), columns=encodings[0, 1:], index=encodings[1:, 0])
        return encodings
    except:
        error_msg = "There is an error with the input sequence"
        print(error_msg)
        return False

info = _AAC(seq)
# =============================================================================
# 
# =============================================================================

# Print out the results
print("Molecular Weight: ", mol_weight)
print("Amino Acid Count: ", aa_count)
print("Amino Acid Percentages: ", aa_percent)
print("Charge at pH 7.4: ", charge)
print("Gravy: ", gravy)
print("Instability Index: ", instability_index)
print("Isoelectric Point: ", iso_point)
print("Secondary Structure Fractions: ", secondary_struct)


list_names = ['name', 'sequence', 'mutations list']
df_input = pd.DataFrame(columns=list_names)
for i in seq:
    df_input.loc[len(df_input)] = i 
dataset = df_input.explode('mutations list').reset_index(drop=True)

#dataset = pd.read_csv('/home/fathima/Medha/Medh_MutBLESS_revision/generalized_dataset/dataset_other_cancers.csv')
Polar = ['N','Q', 'S', 'T', 'P']
aromatic = ['Y', 'F', 'W']
pos_charge = ['K', 'R', 'H']
Neg_charge = ['D', 'E']
Sul_containing = ['C','M']
Aliphatic = ['G','A', 'L','I', 'V']
motif_list = ['nM','Mc', 'n_M', 'M_c', 'n__M','M__c', 'tri']
netsurfp = pd.read_csv('P04637_netsurfp.csv')
#physicochem_properties_normalized = pd.read_csv('/home/fathima/Medha/Medh_MutBLESS_revision/July2021/gbm_work/49_properties_normalizedValues.csv')
physicochem_properties_actual = pd.read_csv('49_properties_numerical_Values.csv')
#new_463 = pd.read_csv('/home/fathima/Medha/Medh_MutBLESS_revision/July2021/gbm_work/Nithya_files/463_unique_numerical_properties.csv')
#physicochem_properties_actual1 = pd.concat([physicochem_properties_actual, new_463[physicochem_properties_actual.columns]]).reset_index(drop =True)
aacon_header = ['mut_pos','KABAT', 'JORES', 'SCHNEIDER', 'SHENKIN', 'GERSTEIN',
                'TAYLOR_GAPS', 'TAYLOR_NO_GAPS', 'VELIBIL', 'KARLIN', 'ARMON',
                'THOMPSON', 'NOT_LANCET', 'MIRNY', 'WILLIAMSON', 'LANDGRAF', 
                'SANDER', 'VALDAR', 'SMERFS']


# all_prop = open('/home/fathima/Medha/Medh_MutBLESS_revision/generalized_dataset/feature_file.csv')

# =============================================================================
# This set of snippet will calculate different physicochemical properties, 
# PSSM scores, and conservation scores calculated from different methods
# =============================================================================
# for i in range(5):

def compute_features_mutation(i, all_prop = all_prop):
    print("HEREEEEEEE",i)
    #filename = input("Enter the name of the FASTA file: ")
    # mutation = input("List the mutations separated by comma: ")
    # list_mut  = mutation.split()
    # all_prop = pd.DataFrame()
    # =============================================================================
    # Sequence based properties calculated using Biopython module
    # record = SeqIO.read(filename, "fasta")
    # sequence = str(record.seq)
    uniprot_id = dataset['name'][i]
    site = dataset['mutations list'][i][:-1]
    sul_c = 0
    polar = 0
    aliphatic = 0
    neg_c = 0
    pos_c = 0
    arom = 0
    pos = int(site[1:])
    #     sequence = next(iter(uniprot_seq.loc[uniprot_seq['uniprot_id'] == uniprot_id, 'Sequence']))
    sequence = dataset['sequence'][i]
    if len(sequence) >= pos:
        dict_prop = dict()
        wild = site[0]
        mut = dataset['mutations list'][i][0]
        wild_mut = wild+mut
        dict_prop['pos'] = int(site[1:])
        # dict_prop['class'] = dataset['label'][i]
        dict_prop['Uniprot ID'] = uniprot_id
        # dict_prop['Gene Name'] = dataset['gene'][i]
        dict_prop['Site'] = site  
        f4 = open('P04637.pssm', 'r').readlines()[2:]
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
            dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] == property1][wild].tolist()[0]  
            if property1 == "pK'":
                dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] == "pK'"][wild].tolist()[0]
                dict_prop[property1] =  b[j]
            else:
                dict_prop[property1] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] == property1][wild].tolist()[0]
                dict_prop[property1] =  b[j]
                dict_prop[property1+'_diff'] = physicochem_properties_actual.loc[physicochem_properties_actual['index'] == property1][wild_mut[-1]].tolist()[0] - physicochem_properties_actual.loc[physicochem_properties_actual['index'] == property1][wild_mut[0]].tolist()[0]
            j += 1
        dict_prop['neg_charge'] = neg_c
        dict_prop['pos_charge'] = pos_c
        dict_prop['polar'] = polar
        dict_prop['aromatic'] = arom
        dict_prop['S_containing'] = sul_c
        dict_prop['aliphatic'] = aliphatic
        
        try:
            dict1 = {'A':f4[pos].strip().split()[2:][0],'R':f4    [pos].strip().split()[2:][1],
                'N':f4[pos].strip().split()[2:][2], 'D':f4[pos    ].strip().split()[2:][3],
                      'C':f4[pos].strip().split()[2:][4], 'Q':f4[pos    ].strip().split()[2:][5], 
                  'E':f4[pos].strip().split()[2:][6], 'G':f4[pos    ].strip().split()[2:][7],
                      'H':f4[pos].strip().split()[2:][8], 'I':f4[pos    ].strip().split()[2:][9], 
                  'L':f4[pos].strip().split()[2:][10], 'K':f4[pos    ].strip().split()[2:][11],
                      'M':f4[pos].strip().split()[2:][12], 'F':f4[pos    ].strip().split()[2:][13], 
                  'P':f4[pos].strip().split()[2:][14], 'S':f4[pos    ].strip().split()[2:][15],
                      'T':f4[pos].strip().split()[2:][16], 'W':f4[pos    ].strip().split()[2:][17],
                  'Y':f4[pos].strip().split()[2:][18], 'V':f4[pos    ].strip().split()[2:][19]}
            dict_prop['pssm_score1'] =  int(dict1[site[0]])   
        #             dict_prop['pssm_score1'] =  dict1[site[0]] 
            dict_prop['pssm_score2'] = sum([int(i) for i in dict1.values()])/20
            dict_prop['pssm_score3'] = int(dict1[wild_mut[-1]])-int(dict1[    wild_mut[0]])  
            cons= f4[pos].strip('\n')[90:].split()[21:22]
            if len(cons) ==0:
                dict_prop['conservation']= 0
            else:
                for element in cons:
                    dict_prop['conservation']= float(element)
        except (KeyError, IndexError) as error:
                dict_prop['pssm_score1'] =  0
                dict_prop['pssm_score2'] = 0
            #         dict_prop['pssm_score3'] = 0
                dict_prop['conservation']= 0
        f_aacon = pd.read_csv('P04637.features', sep = '\t', skiprows=1, names = aacon_header)
        try:
            mm = f_aacon.iloc[pos-1]
            for item in aacon_header:
                dict_prop[item] = mm[item]
        except IndexError:
            for item in aacon_header:
                dict_prop[item] = 0
    #     print(uniprot_id)
        print(dict_prop)
        all_prop = all_prop.append(dict_prop, ignore_index =True)
        return all_prop
# # all_prop['pos'] = all_prop['pos'].astype(int)
# pool = mp.Pool(mp.cpu_count()-2, maxtasksperchild=1)
# re1 = pool.map(calls, range(1))
# # re = pool.map(calls, range(10))
# df11 = pd.concat(re)
df11 = pd.DataFrame()
# for i in range(len(dataset)):
def mp_fuc(i):
    re1 = compute_features_mutation(i)
    return re1
    # df11 = df11.append(re1,ignore_index = True)
    
pool = mp.Pool(mp.cpu_count()-2)
re11 = pool.map(mp_fuc,range(len(dataset)))

# =============================================================================
# Create a dataframe from the values obtained from input
# =============================================================================
