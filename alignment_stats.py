#!/usr/bin/python -tt

import os,csv,sys,getopt,collections
import pandas as pd


def globIt(crawlDir,extensions=[]):
    """
    stdlib version of a directory crawler searching for files with specific extensions
    """
    import glob
    
    list_of_files = []
    for extension in extensions:
        filePath = '{}/*{}'.format(crawlDir,extension)
        print(filePath)
        list_of_files += glob.glob(filePath)
    latest_file = max(list_of_files, key=os.path.getctime)
    return list_of_files, latest_file


def parse_mafft_file(mafft_file):
    """
    parses mafft pairwise alignment fasta file and yields fasta tuple
    """
    with open(mafft_file,'r') as mafft:
        seq_id = None
        seq = []
        for line in mafft:
            line = line.strip()
            if line.startswith('>'):
                if seq_id == None:
                    seq_id = line
                else:
                    yield seq_id, "".join(seq)
                    seq_id = line
                    seq = []                    
            else:
                seq_on_line = list(line)
                seq.extend(seq_on_line)


def compare_seqs(seq1,seq2):
    """
    seq1 will always be CONSENSUS to help determine insertion/deletion/substitution
    """
    insertion_count = 0
    deletion_count = 0
    substitution_count = 0
    length = len(seq1)
    if length != len(seq2):
        print("NOT SAME LENGTH")
        length = [len(seq1),len(seq2)]
    for a, b in zip(seq1, seq2):
        a = a.lower()
        b = b.lower()
        if a != b:                  # bases are different
            if a == '-':
                insertion_count += 1
            elif b == '-':
                deletion_count += 1
            else:
                substitution_count += 1
    return insertion_count,deletion_count,substitution_count,length


def main(argv):
    """
    time python3 alignment_stats.py -i /path/to/file.csv -m /path/to/topologies/
    """
    inputPath = ''
    outfile = ''

    try:
        opts, args = getopt.getopt(
            argv, "hi:o:", ["inputPath=", "outfile="])
    except getopt.GetoptError:
        print('alignment_stats.py -i <inputPath> -o <outfile> \n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('alignment_stats.py -i <inputPath> -o <outfile>\n\n')
            sys.exit()
        elif opt in ("-i", "--inputPath"):
            inputPath = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg
    
    print('Input path is ', inputPath)
    # add trailing directory separator if needed
    if not inputPath.endswith('/'):
        inputPath += '/'
    outPath = f"{inputPath}/stats/"
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    print('Output path is ', outPath)
    file_list, newest_file = globIt(inputPath,['.mafft'])
    out = open(outPath+'outfile.txt','w')
    out.write('Query_ID\tConsensus_ID\tInsertions\tDeletions\tSubstitutions\tSequence_Length\n')
    for mafft_file in file_list:                            # parallelize here
        this_fasta = []
        for seq_id,seq in parse_mafft_file(mafft_file):
            if 'consensus' in seq_id.lower():
               seq_id = f"0{seq_id[1:]}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()
        insertion_count,deletion_count,substitution_count,length = compare_seqs(this_fasta[0][1],this_fasta[1][1])
        outline = '\t'.join([this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length])+'\n'
        out.write(outline)
    out.flush()
    out.close()

        
            

if __name__ == "__main__":
    main(sys.argv[1:])
