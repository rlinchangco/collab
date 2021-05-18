#!/usr/bin/python -tt

import os,csv,sys,getopt,collections
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
    # add special characters to remove
    #extras = ['/', '\\', '']
    # df.columns = [sanitize_string(col_name, extras) for col_name in list(
    #    df.columns)]    # sanitize header string

    return df


def find_match(match_map,matchFile,topo):
    out = open(matchFile+topo[:-4]+"added.txt",'w')
    header = ['Outgroup_monophyly','Ancestral_state','IA','IB','pwd.AB','pd.A','pwd.B','Topology','Number Donor Tips','Number Recipient Tips','Total number of tips']
    new_header = ['Outgroup_monophyly','Ancestral_state','IA','IB','pwd.AB','pd.A','pwd.B','Topology','Number Donor Tips','Number Recipient Tips','Total number of tips','Donor_Date','Recipient_Date','Estimated_Transmission_Time_to_donor','Estimated_Transmission_Time_to_recipient','Pair','donor_diff_date','recipient_diff_date']
    outheader = '\t'.join(new_header)+'\n'
    out.write(outheader)
    with open(matchFile+topo, newline='') as csvfile:
        data = csv.DictReader(csvfile, delimiter = '\t')
        for row in data:
            tot_tips = row['Total number of tips']
            don_tips = row['Number Donor Tips']
            rec_tips = row['Number Recipient Tips']
            mykey = (tot_tips,don_tips,rec_tips)
            if mykey in match_map:
                outlist = []
                for head in header:
                    outlist.append(row[head])
                outline = '\t'.join(outlist+match_map[mykey])+'\n'
                out.write(outline)
                out.flush()
    out.close()


def match_to_file(df1):
    match_map = {}
    for index, row in df1.iterrows():
        tot_tips = row['Total_number_of_tips']
        don_tips = row['Number_Donor_Tips']
        rec_tips = row['Number_Recipient_tips']
        add_cols = [str(row['Donor_Date']),str(row['Recipient_Date']),str(row['Estimated_Transmission_Time_to_donor']),str(row['Estimated_Transmission_Time_to_recipient']),str(row['Pair']),str(row['donor_diff_date']),str(row['recipient_diff_date'])]
        mykey = (str(tot_tips),str(don_tips),str(rec_tips))
        if mykey not in match_map:
            match_map[mykey] = add_cols
    return match_map


def main(argv):
    """
    time python3 add_cols.py -i /path/to/file.csv -m /path/to/topologies/
    """
    inFile = ''
    matchFile = ''

    try:
        opts, args = getopt.getopt(
            argv, "hi:m:o:c:j:a:", ["inputFile=", "matchFile="])
    except getopt.GetoptError:
        print('add_cols.py -i <inputFile> -m <matchFile> \n\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('add_cols.py -i <inputFile> -m <matchFile>\n\n')
            sys.exit()
        elif opt in ("-i", "--inputFile"):
            inFile = arg
        elif opt in ("-m", "--matchFile"):
            matchFile = arg
    
    print('Input file is ', inFile)
    print('Match directory is ', matchFile)
    df1 = multi_parse(inFile)
    topos = os.listdir(matchFile)
    print("Load matching")
    match_map = match_to_file(df1)
    for topo in topos:
        if not topo.startswith('.'):
            print(f"Running for {topo}")
            find_match(match_map,matchFile,topo)
        

if __name__ == "__main__":
    main(sys.argv[1:])
