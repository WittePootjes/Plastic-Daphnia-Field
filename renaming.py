import pandas as pd
import os
import time
"""
Adding barcode names to the sequence list
FROM: raw barcode 
TO: barcode names ie: F1 and R1
"""

"""
names_list = pd.read_excel('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet2-with primers names.xlsx')
barcode_list = pd.read_excel ('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet.xlsx', header=0)

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
barcode_list.to_excel('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet_with_names.xlsx') 
barcode_list = pd.read_excel ('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet_with_names.xlsx', header=0)
#Adding sample's names
samples = pd.read_excel ('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Sample_names.xlsx', header=0, index_col=0)

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
barcode_list.to_excel('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet_complete.xlsx')
"""
#Now, Renaming actual files in the folder:

barcode_list = pd.read_excel ('/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/Barcodesheet_complete.xlsx', header=0)
folder = r'/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/RAW' #Path to the folder
files = list(os.listdir(folder))
print(range(len(barcode_list['primer_pair'])))
# iterate all files from a directory
# Construct old file name
counter1 = 0
counter2 = 0
counter3 = 0
for i in range(len(barcode_list['primer_pair'])):
    for j in range(len(files)):
        if barcode_list['primer_pair'][i] + '.R1.fastq' in files[j]:
            counter1 = counter1 + 1
            old = folder + '/' + files[j]
            new = folder + '/' + barcode_list['Sample names'][i] + '_R1.fastq'
            os.rename(old, new)
            print('Iteration' + '\t' + str(counter1) + '\t' + 'from' + '\t' + old + '\t' +'to'+ '\t' + new)
        elif barcode_list['primer_pair'][i] + '.R2.fastq' in files[j]:
            counter2 = counter2 + 1
            old = folder + '/' + files[j]
            new = folder + '/' + barcode_list['Sample names'][i] + '_R2.fastq'
            os.rename(old, new)
            print('Iteration',str(counter2),'from',old,'to',new)

""" 
long_data = pd.read_csv ('/Users/u0145079/long_data3.csv', header=0)
print(long_data)

def progress_percent(i, total):
    return round((i/total)*100,2)

total = len(long_data['Type'])
start_time = time.time()
for i in range(len(long_data['Type'])):
    print(str(i) + '/' + str(total) + ' (' + str(progress_percent(i, total)) + '%)', end="\r", flush=True)
    if long_data['Type'][i] == 'Regular pond (accessible for everyone)':
        long_data['Type'][i] = 'Regular pond'
    elif long_data['Type'][i] == 'Drinking water Reservoir (closed for public)':
        long_data['Type'][i] = 'Drinking water Reservoir'

print("--- %s seconds ---" % (time.time() - start_time))


long_data.to_csv('/Users/u0145079/Desktop/Daphnia_MiSeq/longdata_final.csv', sep=',') """

