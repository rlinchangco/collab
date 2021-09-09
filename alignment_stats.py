#!/usr/bin/python -tt

import os,csv,sys,getopt,collections
import pandas as pd
import openpyxl

def globIt(crawlDir,extensions=[]):
    """
    stdlib version of a directory crawler searching for files with specific extensions
    """
    import glob
    
    list_of_files = []
    for extension in extensions:
        filePath = '{}*{}'.format(crawlDir,extension)
        print(filePath)
        list_of_files += glob.glob(filePath)
    #latest_file = max(list_of_files, key=os.path.getctime)
    return list_of_files#, latest_file


def readFasta(fastaFile):
    from Bio import SeqIO
    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    for fasta in fasta_sequences:
        yield fasta.id, str(fasta.seq)


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


def plot_histo(x,outPath,col):
    """
    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(17,17))
    plt.hist(x)
    plt.title(f"{col} Histogram")
    plt.xlabel(f"{col}")
    plt.ylabel('Frequency')
    #plt.show() 
    figure_file = "%s%s.png" % (outPath,col)
    fig.savefig(figure_file,dpi=fig.dpi)
    plt.close('all')    


def plotly_plot(df,outPath,col):
    """
    """
    import plotly.express as px
    fig = px.histogram(df, x=col)
    fig.write_html(f"{outPath}{col}.html")


def main(argv):
    inputPath = ''

    try:
        opts, args = getopt.getopt(
            argv, "hi:", ["inputPath="])
    except getopt.GetoptError:
        print('alignment_stats.py -i <inputPath>\n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('alignment_stats.py -i <inputPath>\n\n')
            print('time python3 alignment_stats.py -i /path/to/mafft_files/')
            sys.exit()
        elif opt in ("-i", "--inputPath"):
            inputPath = arg
    
    print('Input path is ', inputPath)
    # add trailing directory separator if needed
    if not inputPath.endswith('/'):
        inputPath += '/'
    outPath = f"{inputPath}/stats/"
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    print('Output path is ', outPath)
    file_list = globIt(inputPath,['.mafft'])
    #print(file_list)
    out = open(outPath+'statsoutfile.txt','w')
    out.write('Query_ID\tConsensus_ID\tInsertions\tDeletions\tSubstitutions\tSequence_Length\n')
    df_row = 0
    df = pd.DataFrame(columns=['Query_ID','Consensus_ID','Insertions','Deletions','Substitutions','Sequence_Length'])    
    consensus_id = None
    for mafft_file in file_list:
        this_fasta = []
        for seq_id,seq in readFasta(mafft_file):
            if 'consensus' in seq_id.lower():       # modify to force order
               seq_id = f"0{seq_id}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()                           # sort for order
        insertion_count,deletion_count,substitution_count,length = compare_seqs(this_fasta[0][1],this_fasta[1][1])
        df.loc[df_row] = [this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length]
        df_row += 1
        consensus_id = this_fasta[0][0][1:]
        outline = '\t'.join([this_fasta[1][0],this_fasta[0][0],str(insertion_count),str(deletion_count),str(substitution_count),str(length)])+'\n'
        out.write(outline)
    for col in df.columns.to_list():
        if 'ID' not in col:
            plot_histo(df[col], outPath, col)
            plotly_plot(df, outPath, col)
    df_out = f"{outPath}{consensus_id}_statsoutfile.xlsx"
    df.to_excel(df_out,index=False)      # need excelwriter like openpyxl
    out.flush()
    out.close()

if __name__ == "__main__":
    main(sys.argv[1:])
