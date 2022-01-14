#!/usr/bin/python -tt

"""
Script explanation
"""

__author__ = "Greg Linchangco"
__copyright__ = "Copyright 2021, Greg Linchangco"
__credits__ = ["Linchangco"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Linchangco"
__email__ = "gvl@lanl.gov"
__status__ = "Alpha-Testing"

import sys,os,getopt,collections,subprocess
import pandas as pd


def multi_parse(data_file, header_row=0):
    """
    Parses data file (accepts xlsx,tsv,csv) WITH header as first row.
    # read in chunks at a time(returns iterable if chunksize > 0)
    for df_chunk in pd.read_csv('data.csv', index_col=0, chunksize=8):
        print(df_chunk, end='\n\n')
    # find mean of duplicate IDs
    df.groupby('sample_id', as_index=False).mean()
    """
    df = pd.DataFrame()
    # reads csv, txt(tsv), and xlsx files
    if data_file.endswith('.csv'):
        df = pd.read_csv(data_file, header=0)
    elif data_file.endswith('.tsv') or data_file.endswith('.txt'):
        df = pd.read_csv(data_file, delimiter='\t', header=0)
    elif data_file.endswith('.xlsx'):
        df = pd.read_excel(data_file, header=header_row)
    else:
        print(data_file)
        print(f"\n\nUnsupported file format\n\nPlease reformat...{data_file}")
        sys.exit()
    return df


def question_existence(path,item_type):
    """
    Check if directory or file exists
    Inputs:
    path= string of path/file
    item_type= 'dir' or 'file'
    """
    check_result = False
    if item_type == 'dir':
        check_result = os.path.isdir(path)
    elif item_type == 'file':
        check_result = os.path.isfile(path)
    else:
        print("Path/File not recognized...")
        sys.exit()
    if check_result == False and item_type == 'dir':
        subprocess.check_call(["mkdir",path])
        #os.system('mkdir '+path)
        print("Creating directory\t{0}".format(path))


# def parse_multiline_fasta(fastaFile):
#     """
#     convert multiline fasta file to single sequence line fasta
#     """
#     out = fastaFile.split('.')[0]
#     with open(fastaFile) as inputFile, open(out+'.oneline.fasta', 'w') as f_output:
#         block = []
#         for line in inputFile:
#             if line.startswith('>'):
#                 if block:
#                     f_output.write(''.join(block) + '\n')
#                     block = []
#                 f_output.write(line)
#             else:
#                 block.append(line.strip())

#         if block:
#             f_output.write(''.join(block) + '\n')


# def make_fastaDict(fastaFile):
#     """
#     cache fasta to memory for lookup later
#     """
#     fastaDict = {}                  # {subtype:(orig_header,sequence)}
#     with open(fastaFile) as inputFile:
#         for line in inputFile:
#             line = line.strip()
#             if line.startswith('>'):
#                 seqHeader = line
#                 subtype = seqHeader.split('>CONSENSUS_')[1]
#                 #print(subtype,seqHeader)
#                 seq = next(inputFile).strip()
#                 #print(seq)
#                 if subtype not in fastaDict:
#                     fastaDict[subtype] =  (seqHeader,seq)
#                 else:
#                     print(f"DUPLICATE FOR {subtype}")
#     return fastaDict


# def create_fasta(subDF_tuple, outputfile, fasta_dict):
#     """
#     creates fasta files PER subtype
#     each sequence header follows format:
#     >Accession_Patient Id_HXB2/MAC239 start_HXB2/MAC239 stop
#     REFACTOR TO BE MODULAR TO ACCEPT ANY COLUMNS
#     """
#     nomatch = set()
#     subtype = subDF_tuple[0]
#     out = f"{outputfile}{subtype}.fasta"
#     fastaOut = open(out,'w')
#     if subtype in fasta_dict:
#         header,seq = fasta_dict[subtype]
#         fastaOut.write(f"{header}\n{seq}\n")
#     else:
#         print(f"{subtype} NOT found in consensus fasta file")
#         nomatch.add(subtype)
#     for index, row in subDF_tuple[1].iterrows():
#         #print(row["Accession"], row["Patient Id"], row["HXB2/MAC239 start"], row["HXB2/MAC239 stop"])
#         #header = f">{row['Accession']}_{row['Patient Id']}_{row['HXB2/MAC239 start']}_{row['HXB2/MAC239 stop']}\n{row['Sequence']}\n"
#         fastaOut.write(header)
#     fastaOut.flush()
#     fastaOut.close()
#     return nomatch


def readFasta(fastaFile):
    from Bio import SeqIO
    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    for fasta in fasta_sequences:
        yield fasta.id, str(fasta.seq)


def sortUniqueSubtypeSeqs(fastaFile:str,splitC:str,ind:list,coords:dict=None) -> collections.defaultdict:
    """
    Takes in fasta file and yields each header and sequence.
    Classifies unique subtypes by header parsing (based on format provided) into dictionary.
    Separate by:
    subtype (1)
    patient (7)
    Filter by:
    start/end coordinates (14,15)
    """
    fasta_dict = collections.defaultdict(dict)
    fcounter = 0
    for seq_id,seq in readFasta(fastaFile):
        header_list = seq_id.split(splitC)
        subtype = header_list[ind[0]]
        patient = header_list[ind[1]]
        unique_key = (subtype,patient)
        if unique_key not in fasta_dict:
            fasta_dict[unique_key][seq_id] = seq
        elif unique_key in fasta_dict and seq_id not in fasta_dict[unique_key]:
            fasta_dict[unique_key][seq_id] = seq
        else:
            print(f"DUPLICATE Sequence ID:\t{seq_id}")
        fcounter += 1
    print(f"Total seqs in {fastaFile}:\t{fcounter}")
    print(f"Unique subtypes in {fastaFile}:\t{len(fasta_dict)}")
    return fasta_dict


def findAndWriteSubtypes(fastaFile:str,outputfile:str,consensus_dict:collections.defaultdict):
    """
    Writes out separate fasta files of each subtype.
    If subtype found in consensus fasta file, it is appended as the last sequence of subtype file.
    """
    fasta_dict = sortUniqueSubtypeSeqs(fastaFile,'.',[0,6])
    for subt,fasta in fasta_dict.items():
        if subt in consensus_dict:                                      # ONLY those with matching consensus are written
            thisSub = open(f"{outputfile}{subt}.seqs.fasta","w")
            for con_header,con_seq in consensus_dict[subt].items():
                thisSub.write(f"{con_header}\n{con_seq}\n")            
            for header,seq in fasta.items():
                thisSub.write(f"{header}\n{seq}\n")
            thisSub.flush()
            thisSub.close()


def main(argv):
    """
    time python3 filter_subset.py -i /path/to/consensusFile.fasta -o /path/to/results/directory/ -f /path/to/fastaFile.fasta
    """
    inFile = ''
    outputfile = ''
    fastaFile = ''
    coordMap = ''

    try:
        opts, args = getopt.getopt(
            argv, "hi:o:f:m:", ["inputFile=", "ofile=","fastaFile=","coordMap="])
    except getopt.GetoptError:
        print('filter_subset.py -i <inputFile> -o <ofile> -f <fastaFile> -m <coordMap>\n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('filter_subset.py -i <inputFile> -o <ofile> -f <fastaFile> -m <coordMap>\n\n')
            sys.exit()
        elif opt in ("-i", "--inputFile"):
            inFile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
            # should  end with '/'
        elif opt in ("-f", "--fastaFile"):
            fastaFile = arg
        elif opt in ("-m", "--coordMap"):
            coordMap = arg
    print(f'Consensus file is {inFile}')
    print(f'Output directory is {outputfile}')
    print(f'Fasta file is {fastaFile}')
    print(f'Gene Coordinate file is {coordMap}')
    if not outputfile.endswith('/'):
        outputfile = f'{outputfile}/'
    question_existence(outputfile,'dir')
    # Read in gene coordinate map of tabular format: name    start    end
    coordDF = multi_parse(coordMap, header_row=0)
    coords = dict(zip((coordDF.start,coordDF.end), coordDF.name))
    consensus_dict = sortUniqueSubtypeSeqs(inFile,'.',[0,6])
    findAndWriteSubtypes(fastaFile,outputfile,consensus_dict,coords)
 

    # if '.oneline.fasta' not in fastaFile:
    #     parse_multiline_fasta(fastaFile)        # convert multiline fasta to single fasta
    #     fastaFile = fastaFile.split('.')[0]+'.oneline.fasta'
    # fasta_dict = make_fastaDict(fastaFile)      # cache consensus fasta file
    # df1 = multi_parse(inFile)                   # read in dataset
    # # list column headers to console
    # print(df1.columns.to_list())
    # # ask user for input()
    # x = input('\nSelect ONE(1) column to group sequences by from above:')
    # if x not in df1.columns.to_list():
    #     print("{} does not exist, please try again...")
    #     x = input()
    # subs = df1[x].value_counts()
    # print(f"\nCounts of Unique elements from {x}")
    # print(subs)                     # show unique counts of elements in that column
    # df_list = []
    # # group into separate dataframes per unique element of column
    # for sub in subs.index.to_list():
    #     df_list.append((sub,df1[df1[x] == sub]))
    # # iterate over all unique element dataframes and create fasta files
    # nomatch = set()
    # for subDF in df_list:
    #     nomatch.update(create_fasta(subDF, outputfile, fasta_dict))
    # matched = len(df_list) - len(nomatch)
    # print(f"{len(nomatch)} subtypes were NOT matched to consensus sequences")
    # print(f"{matched} subtypes WERE matched to consensus sequences")

if __name__ == "__main__":
    main(sys.argv[1:])
