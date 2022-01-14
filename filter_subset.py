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

from enum import unique
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


def readFasta(fastaFile):
    from Bio import SeqIO
    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    for fasta in fasta_sequences:
        yield fasta.id, str(fasta.seq)


def findWithinXPer(seq_ind:int,gene_ind:int,threshold:float=0.90):
    """
    Find if query position is within X% of consensus position
    """
    test_list = [seq_ind,gene_ind]
    if seq_ind != gene_ind:
        hi_val_ind = max(test_list)
        lo_val_ind = min(test_list)
        test_val = lo_val_ind / hi_val_ind
        if test_val >= threshold:
            return True
        else:
            return False
    else:
        return True


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
    coord_list = None
    if coords:
        coord_list = list(map(list, zip(*coords)))         # unpack dictionary keys, unpack list of tuples to list of lists
    fasta_dict = collections.defaultdict(dict)
    fcounter = 0
    kept_counter = 0
    for seq_id,seq in readFasta(fastaFile):
        header_list = seq_id.split(splitC)
        subtype = header_list[ind[0]]
        unique_key = subtype
        start = header_list[-2]
        end = header_list[-1]
        if len(ind) > 1:                                # differentiate between consensus and query
            patient = header_list[ind[1]]
            unique_key = (subtype,patient)
        if coords:                                      # filter by coordinate map
            #if coords.get((start,end)):                 # if seq coords in coord map, keep. Expects exact match...
            for i in len(coord_list[0]):                # checks if within x percent of start and end coordinates (takes longer)
                gene_start = coord_list[0][i]
                gene_end = coord_list[1][i]
                test_start = findWithinXPer(start,gene_start)
                test_end = findWithinXPer(end,gene_end)
                if test_start and test_end:
                    if unique_key not in fasta_dict:
                        fasta_dict[unique_key][seq_id] = seq
                    elif unique_key in fasta_dict and seq_id not in fasta_dict[unique_key]:
                        fasta_dict[unique_key][seq_id] = seq
                    else:
                        print(f"DUPLICATE Sequence ID:\t{seq_id}")
                    kept_counter += 1
        fcounter += 1
    print(f"Total seqs in {fastaFile}:\t{fcounter}")
    print(f"Total seqs kept {fastaFile}:\t{kept_counter}")
    print(f"Unique subtypes in {fastaFile}:\t{len(fasta_dict)}")
    return fasta_dict


def findAndWriteSubtypes(fastaFile:str,outputfile:str,consensus_dict:collections.defaultdict):
    """
    Writes out separate fasta files of each subtype.
    If subtype found in consensus fasta file, it is appended as the last sequence of subtype file.
    """
    fasta_dict = sortUniqueSubtypeSeqs(fastaFile,'.',[0,6])
    for subt,fasta in fasta_dict.items():
        if subt[0] in consensus_dict:                                      # Filter by matching consensus, only these are written
            thisSub = open(f"{outputfile}{subt[0]}.{subt[1]}.seqs.fasta","w")
            for con_header,con_seq in consensus_dict[subt[0]].items():
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
    consensus_dict = sortUniqueSubtypeSeqs(inFile,'_',[0])
    findAndWriteSubtypes(fastaFile,outputfile,consensus_dict,coords)

if __name__ == "__main__":
    main(sys.argv[1:])
