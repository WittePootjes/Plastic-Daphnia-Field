import pandas as pd
import os, sys

"""
tested on: Mac/Linux

GOAL: TO RENAME SAMPLES FROM:
    
    GC121272_TGAGTACG-GCTCTAGT.220329_2022_03_29_Miseq_Run8_Nikola_Manon_Naina.220415.MiSeq.FCB.lane1.gcap_20_01.R1.fastq.gz
    
    TO:
    
    SAMPLE_NAME_{NUMBER}.R1.fastq ie: N23.R1.fastq

STEPS:
    1. Remove the unnecessary middle and front of the sample name = run tidyup.bash script
    2. Make sure following are ready:
        a) Excel_file_1 with:
            Index2 (i5) column with FORWARD primers
            Index1 (i7) column with REVERSE primers
        b) Excel_file_2 with:
            Name_1 column with FORWARD primer names ie. F1,F2,F3 etc.
            Index2 (i5) with their sequences ie. atcgtacg etc.
            Name_2 column with REVERSE primer names ie. R1,R2,R3 etc.
            Index1 (i7) column with their sequences ie. gacatagt etc
        c)Excel_file_3 with NAMES of the SAMPLE and their primer pairs in rows and columns:
            Columns: Reverse primer names -> R1,R2,R3 etc.
            Rows: Forward primer names -> F1,F2,F3 etc.
            So that each sample ie. N1, N2, N3 corresponds to different primer pairs ie. F1 + R3.
        d)path_to_reads -> a full path to the folder with your raw reads
        e)outputs -> a full path to the folder with your outputs

"""

def rename(Excel_file_1,Excel_file_2,Excel_file_3,path_to_reads,outputs): # each excel file needs to be provided as a full path
    
    #check if all files exist
    rules = [os.path.exists(Excel_file_1),
    os.path.exists(Excel_file_2),
    os.path.exists(Excel_file_3),
    os.path.exists(path_to_reads),
    os.path.exists(outputs)]
    
    if not all(rules):
        sys.exit('One of the files does not exit, make sure you supplied a full path ie. /Users/Anna/Projet/Excel_file_1')

    #read excel to pd
    barcode_list = pd.read_excel(Excel_file_1)
    names_list = pd.read_excel(Excel_file_2)

#Renaming raw barcodes to primer names
    #Forward
    new_column=[]
    #Forward primers
    for i in range(len(barcode_list)):
        for j in range(len(names_list)):
            if barcode_list['Index2 (i5)'][i] == names_list['Index2 (i5)'][j]:
                new_column.append(names_list['Name_1'][j])

    barcode_list['Index2 (i5)_Names']=new_column

    #Reverse primers
    new_column=[]
    #Forward primers
    for i in range(len(barcode_list)):
        for j in range(len(names_list)):
            if barcode_list['Index1 (i7)'][i] == names_list['Index1 (i7)'][j]:
                new_column.append(names_list['Name_2'][j])

    barcode_list['Index1 (i7)_Names']=new_column
    path_save = outputs + '/Barcodesheet_with_names.xlsx'
    barcode_list.to_excel(path_save)

#Adding sample's names for each primer name pair
    
    samples = pd.read_excel (Excel_file_3, header=0, index_col=0)

    columns=[] 
    for i in range(len(samples.index)):
        for j in range(len(samples.columns)):
            for k in range(len(barcode_list)):
                if samples.index[i]==barcode_list['Index2 (i5)_Names'][k] and samples.columns[j]==barcode_list['Index1 (i7)_Names'][k]:
                    columns.append(samples.iloc[i,j])

                
    barcode_list['Sample names']=columns

    #Join forward and reverse with "-" to construct pairs like in the actual sample names:

    barcode_list['primer_pair']=barcode_list['Index1 (i7)'] + '-' + barcode_list['Index2 (i5)']
    barcode_list['primer_pair'] = barcode_list['primer_pair'].apply(lambda row : row.upper())
    path_save = outputs + '/Barcodesheet_complete.xlsx'
    barcode_list.to_excel(path_save)

#Now, Renaming actual files in the folder:

    folder = path_to_reads #Path to the folder
    files = list(os.listdir(folder))
    # iterate all files from a directory
    # Construct old file name
    counter1 = 0
    counter2 = 0
    for i in range(len(barcode_list['primer_pair'])):
        for j in range(len(files)):
            if barcode_list['primer_pair'][i] + '.R1.fastq' in files[j]:
                counter1 = counter1 + 1
                old = folder + '/' + files[j]
                new = folder + '/' + barcode_list['Sample names'][i] + '.R1.fastq'
                os.rename(old, new)
                print('Iteration' + '\t' + str(counter1) + '\t' + 'from' + '\t' + old + '\t' +'to'+ '\t' + new)
            elif barcode_list['primer_pair'][i] + '.R2.fastq' in files[j]:
                counter2 = counter2 + 1
                old = folder + '/' + files[j]
                new = folder + '/' + barcode_list['Sample names'][i] + '.R2.fastq'
                os.rename(old, new)
                print('Iteration',str(counter2),'from',old,'to',new)
    return('All files changed!')

#Unit testing
if __name__ == "__main__":
    rename(Excel_file_1='/Users/u0145079/Desktop/Daphnia MiSeq/test/Excel_file_1.xlsx',
            Excel_file_2='/Users/u0145079/Desktop/Daphnia MiSeq/test/Excel_file_2.xlsx',
            Excel_file_3='/Users/u0145079/Desktop/Daphnia MiSeq/test/Excel_file_3.xlsx',
            path_to_reads='/Users/u0145079/Desktop/Daphnia MiSeq/test/seq',
            outputs='/Users/u0145079/Desktop/Daphnia MiSeq/test/output')
