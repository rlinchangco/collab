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
        print(f"Files to parsed from:\n{filePath}")
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
    seq1 will always be CONSENSUS/reference to help determine insertion/deletion/substitution
    """
    insertion_count = 0
    deletion_count = 0
    substitution_count = 0
    length = len(seq1)
    if length != len(seq2):
        print("NOT SAME LENGTH")
        length = [len(seq1),len(seq2)]
    apreceding = 0
    atrailing = 0
    bpreceding = 0
    btrailing = 0    
    # account for preceding or trailing "gaps"    
    if seq1.startswith('-'):
        apreceding = len(seq1) - len(seq1.lstrip('-'))
    if seq1.endswith('-'):
        atrailing = len(seq1) - len(seq1.rstrip('-'))
    if seq2.startswith('-'):
        bpreceding = len(seq2) - len(seq2.lstrip('-'))
    if seq2.endswith('-'):
        btrailing = len(seq2) - len(seq2.rstrip('-'))
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
    insertion_count = insertion_count - apreceding - atrailing
    deletion_count = deletion_count - bpreceding - btrailing
    total_length = length - apreceding - atrailing - bpreceding - btrailing
    distance_ham = insertion_count+deletion_count+substitution_count
    distance_var = distance_ham*(total_length - distance_ham)/total_length
    return insertion_count,deletion_count,substitution_count,total_length,distance_ham,distance_var


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


def mafft_stats(file_list,desired_year,outPath):
    """
    """
    df_row = 0
    df = pd.DataFrame(columns=['Year','Query_ID','Consensus_ID','Insertions','Deletions','Substitutions','Sequence_Length','Hamming_Distance','Distance Variance'])    
    for mafft_file in file_list:
        this_fasta = []
        year = int(mafft_file.split('.')[2])
        for seq_id,seq in readFasta(mafft_file):
            if 'consensus' in seq_id.lower():       # modify to force order
                seq_id = f"!{seq_id}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()                           # sort for order
        insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var = compare_seqs(this_fasta[0][1],this_fasta[1][1])
        df.loc[df_row] = [year,this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var]
        df_row += 1
        consensus_id = this_fasta[0][0][1:]
        outline = '\t'.join([this_fasta[1][0],consensus_id,str(insertion_count),str(deletion_count),str(substitution_count),str(length)])+'\n'
    for col in df.columns.to_list():
        if 'ID' not in col:
            plot_histo(df[col], outPath, col)
            plotly_plot(df, outPath, col)
    print(df.shape)
    if desired_year:
        df = df[df.Year <= desired_year]
        print(f"Filtered by {desired_year}\nRows/Columns remaining:{df.shape}")
    df_csv = f"{outPath}{consensus_id}_statsoutfile.csv"
    df.to_csv(df_csv,index=False)
    # df_xlsx = f"{outPath}{consensus_id}_statsoutfile.xlsx"
    # df.to_excel(df_xlsx,index=False)      # need excelwriter like openpyxl    


def fasta_stats(file_list,outPath,qualifier=None):
    """
    """
    df_row = 0
    df = pd.DataFrame(columns=['Query_ID','Mapped_ID','Insertions','Deletions','Substitutions','Sequence_Length','Hamming_Distance','Distance Variance'])    
    for fasta_file in file_list:
        this_fasta = []
        for seq_id,seq in readFasta(fasta_file):
            if '2002' in seq_id.lower():       # modify to force order
                seq_id = f"!{seq_id}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()                           # sort for order
        insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var = compare_seqs(this_fasta[0][1],this_fasta[1][1])
        df.loc[df_row] = [this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var]
        df_row += 1
        consensus_id = this_fasta[0][0][1:]
        outline = '\t'.join([this_fasta[1][0],consensus_id,str(insertion_count),str(deletion_count),str(substitution_count),str(length)])+'\n'
    for col in df.columns.to_list():
        if 'ID' not in col:
            plot_histo(df[col], outPath, col)
            plotly_plot(df, outPath, col)
    print(df.shape)
    df_csv = f"{outPath}{consensus_id}_statsoutfile.csv"
    df.to_csv(df_csv,index=False)
    # df_xlsx = f"{outPath}{consensus_id}_statsoutfile.xlsx"
    # df.to_excel(df_xlsx,index=False)      # need excelwriter like openpyxl 


def main(argv):
    inputPath = ''
    desired_year = None
    file_end = 'mafft'

    try:
        opts, args = getopt.getopt(
            argv, "hi:y:f:", ["inputPath=", "desiredYear=", "fileEnding="])
    except getopt.GetoptError:
        print('alignment_stats.py -i <inputPath> -y <desiredYear>\n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('alignment_stats.py -i <inputPath>\n\n')
            print('time python3 alignment_stats.py -i /path/to/mafft_files/ -y year')
            sys.exit()
        elif opt in ("-i", "--inputPath"):
            inputPath = arg
        elif opt in ("-y", "--desiredYear"):
            desired_year = int(arg)
        elif opt in ("-f", "--fileEnding"):
            file_end = arg

    print('Input path is ', inputPath)
    # add trailing directory separator if needed
    if not inputPath.endswith('/'):
        inputPath += '/'
    outPath = f"{inputPath}stats/"
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    print('Output path is ', outPath)
    file_list = globIt(inputPath,[file_end])
    if file_end == 'mafft':
        mafft_stats(file_list,desired_year,outPath)
    elif file_end == 'fasta':
        fasta_stats(file_list,outPath)
    

if __name__ == "__main__":
    main(sys.argv[1:])
