import os 
import pandas as pd


datesDF= pd.read_csv('../GISAID-crawler/datesDF.csv')
download_folder= "../GISAID-crawler/downloads"

gisaid_downloads= os.listdir(download_folder)


def load_fasta(pathTo):
    print("pathTo",pathTo)
    f=open(pathTo)
    x=f.readline()[:-1]
    seqs= []
    totalSeqs=0
    while x!='':    
        if ">" in x:
            influenza=x[1:].split("|")
            metadata={}
            metadata["ID1"]=influenza[0]
            metadata["Gene Symbol"]=influenza[1]
            metadata["Strain Name"]=influenza[2]
            metadata["isolate ID"]=influenza[3]
            metadata["INSDC"]=influenza[4]
            metadata["type"]=influenza[5]
            allLines = ''
            x=f.readline()[:-1]
            while x!='':
                if ">" in x:
                    break
                allLines+=x
                x=f.readline()[:-1]
            metadata['sequence']= allLines
            totalSeqs+=1
            seqs.append(metadata)
    return [ seqs,  totalSeqs]

def load_metadata(pathTo):
    metadata_file=pd.read_csv(pathTo,dtype=str)
    return metadata_file

def load_annotated_fasta_with_species(pathTo):
    print("pathTo",pathTo)
    f=open(pathTo)
    x=f.readline()[:-1]
    seqs= []
    totalSeqs=0
    while x!='':    
        if ">" in x:
            influenza=x[1:].split("|")
            metadata={}
            metadata["ID1"]=influenza[0]
            metadata["Gene Symbol"]=influenza[1]
            metadata["Strain Name"]=influenza[2]
            metadata["isolate ID"]=influenza[3]
            metadata["INSDC"]=influenza[4]
            metadata["type"]=influenza[5]
            metadata['Host']=influenza[6]
            allLines = ''
            x=f.readline()[:-1]
            while x!='':
                if ">" in x:
                    break
                allLines+=x
                x=f.readline()[:-1]
            metadata['sequence']= allLines
            totalSeqs+=1
            seqs.append(metadata)
    return [ seqs,  totalSeqs]

#returns dataset with only the gene symbol asked
def filter_by_gene (gene_symbol, seqs):
    filtered_by_gene=[]
    for item in seqs:
        if item['Gene Symbol']==gene_symbol:
            filtered_by_gene.append(item)
    return filtered_by_gene

def filter_by_type(types, seqs):
    filtered_by_type=[]
    for item in seqs:
        if item['type'][0]==types:
            filtered_by_type.append(item)
    return filtered_by_type

def filter_by_host(host, seqs):
    filtered_by_host=[]
    for item in seqs:
        if item['Host'][0]==host:
            filtered_by_host.append(item)
    return filtered_by_host

def write_host_annotated_fasta(seqs, output_location):
    save1= open(output_location, "w+")
    for item in seqs:
        header= (">"+
        item["ID1"] + "|" +
        item["Gene Symbol"]+ "|" +
        item["Strain Name"]+ "|" +
        item["isolate ID"]+ "|" +
        item["INSDC"]+ "|" +
        item["type"]+ "|" +
        item["Host"] + "\n")
        save1.write(header)
        save1.write(item['sequence']+ "\n")
    save1.close()

def write_fasta(seqs, output_location):
    save1= open(output_location, "w+")
    for item in seqs:
        header= (">"+
        item["ID1"] + "|" +
        item["Gene Symbol"]+ "|" +
        item["Strain Name"]+ "|" +
        item["isolate ID"]+ "|" +
        item["INSDC"]+ "|" +
        item["type"]+ "\n")
        save1.write(header)
        save1.write(item['sequence']+ "\n")
    save1.close()
#function to take all fasta and puts them into a single file
def merge_downloads():
    allData = ""
    for file in gisaid_downloads:
        # f=open("downloads/"+file, "r")
        # for count, line in enumerate(f):
        #     pass
        # f.close()
        # print(file)
        fp=open("downloads/"+file, "r")
        print(file)
        dt=fp.read()[:-2]
        allData += dt
        allData += "\n"    

    with open("output/merged_sequences.fasta", "w+") as f:
        f.write(allData)

def merge_metadata():
    appended_data = []
    total_rows=0
    for file in gisaid_downloads:
        if "gisaid" in file:
            print(file)
            df = pd.read_excel("downloads/"+file, sheet_name=None)
            total_rows=total_rows+(df['Tabelle1'].shape[0])
            print(df['Tabelle1'].shape[0])
            print(total_rows)
            appended_data.append(df['Tabelle1'])
    print(total_rows)
    appended_data = pd.concat(appended_data)
    
            # appended_data.reset_index(level=0, inplace=True)

    # with open("output/merged_metadata.csv", "w+") as f:
    #     f.write(str(appended_data))
    appended_data.to_csv('output/merged_metadata.csv')

def filter_df_by_value_eq(df, value):
    return df[df.eq(value).any(axis=1)]

def annotate_fasta_with_species(pathTo):
            response=load_fasta(pathTo)
            sequences=response[0]
            del response
            print("loaded fasta in memory")
        
            metadata= load_metadata("output/simple_metadata.csv")
            metadata = metadata.reset_index()

            print("loaded metadata in memory")
            annotated=[]
            i=0
            for sequence in sequences:
                query_res= filter_df_by_value_eq(metadata, sequence['isolate ID'])
                sequence['Host']=query_res['Host'].values[0]
                annotated.append(sequence)
                print(i)
                i=i+1
            return annotated

def separate_big_file_by_segments():
    response=load_fasta("output/merged_sequences.fasta")
    print("\n"+"\n"+"The total amount of initial sequences is " +str(response[1])+"\n")
    gene_symbols=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    for gene_symbol in gene_symbols:
        response_filter= filter_by_gene(gene_symbol, response[0])
        print("number of"+gene_symbol+ "sequences is: "+str(len(response_filter)))
        output_path= "output/"+gene_symbol+"_sequences.fasta"
        write_fasta(response_filter,output_path)

def separate_file_by_type(types, path):
    response=load_fasta(path)[0]
    response_filter= filter_by_type(types, response)
    path_split= path.split("sequences.fasta")
    output_path= path_split[0]+types+"sequences.fasta"
    write_fasta(response_filter,output_path)

def separate_sequencesObject_by_segment_type_host(types, host, gene, sequences):
    response= filter_by_type(types, sequences)
    response= filter_by_host(host, response)
    response= filter_by_gene(gene, response)
    return response
###ACTIONS
def separate_bigFile_by_segment_type_host():
    gene_symbols=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    hosts=["Human"]
    types=["A"]
    outputDir= "humanOutput/"
    response= load_annotated_fasta_with_species("output/sequences_annotated_with_species.fasta")
    print("\n"+"\n"+"The total amount of initial sequences is " +str(response[1])+"\n")
    for gene_symbol in gene_symbols:
        response_filter= separate_sequencesObject_by_segment_type_host(types[0], hosts[0], gene_symbol, response[0])
        print("number of"+gene_symbol+ "sequences is: "+str(len(response_filter)))
        write_host_annotated_fasta(response_filter,outputDir+gene_symbol+"_A_Human_sequences.fasta" )

def merge_metadata_to_a_single_csv():
    metadata= load_metadata("output/merged_metadata.csv")
    metadata = metadata.reset_index()
    metadata=metadata[['Isolate_Id','Host']]
    metadata.to_csv("output/simple_metadata.csv")
    print("done")

def cut_the_sequences():
    sequences=load_fasta("output/merged_sequences.fasta")[0]
    # sequences=load_fasta("downloads copy/1972-06-011972-06-30.fasta")[0]
    # 1900781
    end_i=10000
    range_of= range(0, end_i, 100)
    file_number=0
    for i in range(0,len(range_of)):
        print("")
        from_i=0
        to_i=0
        if i == 0:
             from_i=range_of[i]
             to_i=range_of[i+1]
            #  write_fasta(sequences[range_of[i]:range_of[i+1]], "output/split_sequences/"+str(file_number)+".fasta")
        elif not i==len(range(0,len(range_of)))-2 and not i==len(range(0,len(range_of)))-1:
            from_i= range_of[i]+1
            to_i=range_of[i+1]
        elif i==len(range(0,len(range_of)))-2:
            from_i= range_of[i]+1
            to_i=end_i


        if not to_i==0:
            print(from_i)
            print(to_i)
            print(file_number)
            write_fasta(sequences[from_i:to_i], "output/split_sequences/"+str(file_number)+".fasta")
        file_number=file_number+1


    # appended_data.to_csv('output/merged_metadata.csv')

