import re
# import psutil
import pandas as pd
from Bio import SeqIO
import argparse

#### FASTA PARSER FUNCTION
def seqParser(fasta_sequences):
        pp=dict()
        ss=dict()
#        mm=dict()

#### Assign place holder       
        eds='NAN'+' '+'99999'+' '+'NAN'
        pp['pp'] = [eds]
#        mm['mm'] = [eds]
        ss['ss'] = [eds]
        
#### sequence parsing        
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

            pre = re.compile(r'(?=(P[^GPAVLIMCFWHKRQNED]))')
#            med = re.compile(r'(?=(P.P))')
            suc = re.compile(r'(?=([^GPAVLIMCFWHKRQNED]P))')
            
            for m in pre.finditer(sequence):
                ids=name+' '+str(m.start(1)+1)+' '+m.group(1)
                if 'pp' in pp:
                    pp['pp'].append(ids)
                else:
                    pp['pp'] = [ids]

            # for m in med.finditer(sequence):
            #     ids=name+' '+str(m.start(1)+1)+' '+m.group(1)
            #     if 'mm' in mm:
            #         mm['mm'].append(ids)
            #     else:
            #         mm['mm'] = [ids]

            for m in suc.finditer(sequence):
                ids=name+' '+str(m.start(1)+1)+' '+m.group(1)
                if 'ss' in ss:
                    ss['ss'].append(ids)
                else:
                    ss['ss'] = [ids]

### CREATE DATAFRAME
        # print(len(mm))
        dfp=pd.DataFrame.from_dict(pp)
        dfp[['id','start_location','Pre_sequence(PX)']] = dfp.pp.str.split(' ',expand=True,)
        dfp.drop(['pp'], axis=1, inplace = True)
        # dfm=pd.DataFrame.from_dict(mm)
        # dfm[['id','start_location','Med_sequence(PXP)']] = dfm.mm.str.split(' ',expand=True,)
        # dfm.drop(['mm'], axis=1, inplace = True)
        dfs=pd.DataFrame.from_dict(ss)
        dfs[['id','start_location','Suc_sequence(XP)']] = dfs.ss.str.split(' ',expand=True,)
        dfs.drop(['ss'], axis=1, inplace = True)
        dfa=pd.merge(dfp,dfs,on=['id','start_location'],how='outer')

### drop placeholder
        dfa=dfa[dfa.id != 'NAN']

### assign datatpes
        convert_dict = {'id': object, 
                        'start_location': int,
                        'Pre_sequence(PX)': object,
                        'Suc_sequence(XP)': object,
                        
                    } 
        
        dfa = dfa.astype(convert_dict) 
        dfa.sort_values(by=['id','start_location'], inplace=True)
        dfa.reset_index(drop=True, inplace=True)
        return dfa

### FILE I/O

parser = argparse.ArgumentParser(description='FastaPatternFinder--> Scan and find tri peptide patterns from Fasta sequences')


parser.add_argument('-i', action='store',
                    dest='input_file',
                    help='Provide input file name --> ex: input.fa')


parser.add_argument('-of', action='store',
                    default='output.csv',
                    dest='output_file',
                    help='Specify output file name --> Default: output.csv')

parser.add_argument('--version', action='version',
                    version='%(prog)s 1.0')

results = parser.parse_args()
inp = results.input_file
opfile= results.output_file

### DATA PROCESSING

fasta_sequences = SeqIO.parse(open(inp),'fasta')
dfa = seqParser(fasta_sequences)

print(dfa)
print("File Processing Completed.")
dfa.to_csv(opfile, index=False)

# print(psutil.virtual_memory())
# 