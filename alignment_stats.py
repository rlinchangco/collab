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


def compare_seqs(seq1,seq2,distance_metric,genome_ref_seq=None,read_frame=0):
    """
    seq1 will always be CONSENSUS/reference to help determine insertion/deletion/substitution
    genome_ref_seq not null for ORF substitution location identification, should be a sequence
    read_frame not null for ORF index of +1,+2, or +3 (1,2, or 3 as argument)
    """
    genome_ref_loc = []                                             # genome reference location list to store substitutions positions
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
    ref_loc = 0
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if genome_ref_seq:
            if genome_ref_seq[i] != '-':
                a = a.lower()
                b = b.lower()
                if a != b:                                                  # bases are different
                    if a == '-':
                        insertion_count += 1
                    elif b == '-':
                        deletion_count += 1
                    else:
                        substitution_count += 1
                        genome_ref_loc.append(ref_loc+read_frame)
                ref_loc += 1                                                # increment ref_loc by a position
        else:
            a = a.lower()
            b = b.lower()
            if a != b:                                                  # bases are different
                if a == '-':
                    insertion_count += 1
                elif b == '-':
                    deletion_count += 1
                else:
                    substitution_count += 1
    insertion_count = insertion_count - apreceding - atrailing
    deletion_count = deletion_count - bpreceding - btrailing
    total_length = length - apreceding - atrailing - bpreceding - btrailing
    distance = None
    if distance_metric == 'percent':
        distance = insertion_count+deletion_count+substitution_count/total_length
    elif distance_metric == 'hamming':
        distance = insertion_count+deletion_count+substitution_count
    else:
        print("NO DISTANCE SPECIFIED, DEFAULTING TO HAMMING")
        distance = insertion_count+deletion_count+substitution_count
    distance_var = distance*(total_length - distance)/total_length
    if genome_ref_seq:
        return insertion_count,deletion_count,substitution_count,total_length,distance,distance_var,genome_ref_loc
    else:
        return insertion_count,deletion_count,substitution_count,total_length,distance,distance_var


def compare_seqs_by_codon(seq1,seq2,distance_metric,genome_ref_seq=None,read_frame=0):
    """
    """
    # print(len(genome_ref_seq),genome_ref_seq.count('-'))          # expectation of 9719 positions for genome_ref_seq
    genome_ref_loc = []                                             # genome reference location list to store substitutions positions
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
    ref_loc = 0
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if genome_ref_seq:
            if genome_ref_seq[i] != '-':
                a = a.lower()
                b = b.lower()
                if a != b:                                                  # bases are different
                    if a == '-':
                        insertion_count += 1
                    elif b == '-':
                        deletion_count += 1
                    else:
                        substitution_count += 1
                        genome_ref_loc.append(ref_loc+read_frame)
                ref_loc += 1                                                # increment ref_loc by a position
        else:
            a = a.lower()
            b = b.lower()
            if a != b:                                                  # bases are different
                if a == '-':
                    insertion_count += 1
                elif b == '-':
                    deletion_count += 1
                else:
                    substitution_count += 1
    print(seq1.count('-'),apreceding+atrailing,insertion_count)
    insertion_count = insertion_count - apreceding - atrailing
    # print(deletion_count,bpreceding,btrailing)
    deletion_count = deletion_count - bpreceding - btrailing
    total_length = length - apreceding - atrailing - bpreceding - btrailing
    distance = None
    if distance_metric == 'percent':
        distance = insertion_count+deletion_count+substitution_count/total_length
    elif distance_metric == 'hamming':
        distance = insertion_count+deletion_count+substitution_count
    else:
        print("NO DISTANCE SPECIFIED, DEFAULTING TO HAMMING")
        distance = insertion_count+deletion_count+substitution_count
    distance_var = distance*(total_length - distance)/total_length
    if genome_ref_seq:
        return insertion_count,deletion_count,substitution_count,total_length,distance,distance_var,genome_ref_loc
    else:
        return insertion_count,deletion_count,substitution_count,total_length,distance,distance_var


def codon_subs(fasta1,fasta2,genome_ref_seq=None):
    """
    finds substitutions only between old sequence and new
    returns dataframe of: alignment index, genome position, codon position
    """
    sub_list = []
    substitution_count = 0
    genome_ref_sans_gap = 0                                     # genome reference base pair position (ignores insertions)
    seq1 = fasta1[1]
    seq2 = fasta2[1]
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        sub = 'no'
        pos = 0
        if genome_ref_seq[i] != '-':
            a = a.lower()
            b = b.lower()
            genome_ref_sans_gap += 1
            if a != b:                                          # bases are different
                if a == '-':                                    # insertion
                    continue
                elif b == '-':                                  # deletion
                    continue
                else:                                           # substitution
                    substitution_count += 1
                    sub = 'yes'
                    if genome_ref_sans_gap % 3 == 0:            # no remainder means codon position 3
                        pos = 3
                    else:
                        div = genome_ref_sans_gap / 3
                        if str(div).endswith('7'):              # .66667 means 2/3, codon position 2
                            pos = 2
                        else:                                   # .33333 means 1/3, codon position 1
                            pos = 1
                    sub_list.append([i+1,genome_ref_sans_gap,pos])  # only add substitutions to list output
    sub_list.insert(0,[f"query:{fasta2[0]}",f"ref:{fasta1[0]}",f"TOTAL SUBS:{substitution_count}"])
    columns = ["alignment_index", "hxb2_genome_position", "codon_position"]
    return pd.DataFrame(sub_list,columns=columns)


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


def mafft_stats(file_list,desired_year,outPath,distance_metric):
    """
    """
    df_row = 0
    df = pd.DataFrame(columns=['Year','Query_ID','Consensus_ID','Insertions','Deletions','Substitutions','Sequence_Length','Hamming_Distance','Distance Variance'])    
    consensus_id = ''
    for mafft_file in file_list:
        this_fasta = []
        year = int(mafft_file.split('.')[2])
        for seq_id,seq in readFasta(mafft_file):
            if 'consensus' in seq_id.lower():       # modify to force order
                seq_id = f"!{seq_id}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()                           # sort for order
        insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var = compare_seqs(this_fasta[0][1],this_fasta[1][1],distance_metric)
        df.loc[df_row] = [year,this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var]
        df_row += 1
        consensus_id = this_fasta[0][0][1:]
    for col in df.columns.to_list():
        if 'ID' not in col:
            plot_histo(df[col], outPath, col)
            plotly_plot(df, outPath, col)
    if desired_year:
        df = df[df.Year <= desired_year]
        print(f"Filtered by {desired_year}\nRows/Columns remaining:{df.shape}")
    df_csv = f"{outPath}{consensus_id}_statsoutfile.csv"
    for col in df.columns[3:]:
        df[col] = pd.to_numeric(df[col])
    medians = df.median()
    medians = medians.to_frame().T
    print(medians)
    pd.concat([pd.Series(['Combined','Combined','Combined']), medians])
    stats = df.describe()
    stats.insert(0, df.columns[1], 'Combined')
    stats.insert(0, df.columns[0], 'Combined')
    stats.insert(0, df.columns[2], 'Combined')
    new_df = pd.concat([df,stats,medians])
    new_df.to_csv(df_csv)
    # df_xlsx = f"{outPath}{consensus_id}_statsoutfile.xlsx"
    # df.to_excel(df_xlsx,index=False)      # need excelwriter like openpyxl    


def fasta_stats(file_list,outPath,distance_metric,qualifier=None):
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
        insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var = compare_seqs(this_fasta[0][1],this_fasta[1][1],distance_metric)
        df.loc[df_row] = [this_fasta[1][0],this_fasta[0][0],insertion_count,deletion_count,substitution_count,length,distance_ham,distance_var]
        df_row += 1
        consensus_id = this_fasta[0][0][1:]
    for col in df.columns.to_list():
        if 'ID' not in col:
            plot_histo(df[col], outPath, col)
            plotly_plot(df, outPath, col)
    print(df.shape)
    df_csv = f"{outPath}{consensus_id}_statsoutfile.csv"
    for col in df.columns[2:]:
        df[col] = pd.to_numeric(df[col])
    stats = df.describe()
    stats.insert(0, df.columns[1], 'Combined')
    stats.insert(0, df.columns[0], 'Combined')
    new_df = pd.concat([df,stats])
    new_df.to_csv(df_csv)
    # df_xlsx = f"{outPath}{consensus_id}_statsoutfile.xlsx"
    # df.to_excel(df_xlsx,index=False)      # need excelwriter like openpyxl 


def codon_stats(file_list,outPath):
    """
    only finds substitutions
    """
    overall_dfs = []
    for fasta_file in file_list:
        this_fasta = []
        for seq_id,seq in readFasta(fasta_file):
            if 'consensus' not in seq_id.lower():       # modify to force order of REFERENCE
                seq_id = f"!REF_{seq_id}"
            this_fasta.append((seq_id,seq))
        this_fasta.sort()                           # sort for order (older sequence should be second index)
        sub_df = codon_subs(this_fasta[2],this_fasta[1],genome_ref_seq=this_fasta[0][1])
        print(f"There are {sub_df.shape[0]} substitutions between {this_fasta[1][0]} and {this_fasta[2][0]}")
        overall_dfs.append(sub_df)

    df_csv = f"{outPath}substitutions_positions_file.csv"
    overall = pd.concat(overall_dfs,sort=False)
    overall.to_csv(df_csv,index=False)


def main(argv):
    inputPath = ''
    desired_year = None                         # YYYY
    file_end = 'mafft'                          # or fasta
    distance_metric = 'hamming'                 # or percent
    reading_frames = None

    try:
        opts, args = getopt.getopt(
            argv, "hi:y:f:d:r:", ["inputPath=", "desiredYear=", "fileEnding=", "distanceMetric=", "readingFrames="])
    except getopt.GetoptError:
        print('alignment_stats.py -i <inputPath> -y <desiredYear> -f <fileEnding> -d <distanceMetric>\n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('alignment_stats.py -i <inputPath>\n\n')
            print('time python3 alignment_stats.py -i /path/to/mafft_files/ -y year -f mafft -d hamming')
            sys.exit()
        elif opt in ("-i", "--inputPath"):
            inputPath = arg
        elif opt in ("-y", "--desiredYear"):
            desired_year = int(arg)
        elif opt in ("-f", "--fileEnding"):
            file_end = arg
        elif opt in ("-d", "--distanceMetric"):
            distance_metric = arg
        elif opt in ("-r", "--readingFrames"):
            reading_frames = arg

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
        mafft_stats(file_list,desired_year,outPath,distance_metric)
    elif file_end == 'fasta' and reading_frames is None:
        fasta_stats(file_list,outPath,distance_metric)
    else:
        codon_stats(file_list,outPath)
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
