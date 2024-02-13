import os
import sqlite3
import pandas as pd
import sys

def merge_downloads(directory="downloads", output_path="output/merged_sequences.fasta"):
    # datesDF= pd.read_csv('datesDF.csv')
    download_folder= "downloads"
    gisaid_downloads= os.listdir(download_folder)
    allData = ""
    for file in gisaid_downloads:
        fp=open("downloads/"+file, "r")
        print(file)
        dt=fp.read()[:-2]
        allData += dt
        allData += "\n"
    with open(output_path, "w+") as f:
        f.write(allData)

def count_lines(filename):
    # Open the file in read mode
    with open(filename, 'r') as file:
        # Use a list comprehension to read all lines and count them
        line_count = len([line for line in file])
    return line_count
def reverse_complement_dna(sequences):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complements = []

    for sequence in sequences:
        complement = ''.join([complement_dict[nucleotide] for nucleotide in sequence])
        reverse_complement = complement[::-1]  # Reverse the complement sequence
        reverse_complements.append(reverse_complement)

    return reverse_complements


#general tasks
def flush_directories(*directory_paths):
    try:
        for directory_path in directory_paths:
            for filename in os.listdir(directory_path):
                file_path = os.path.join(directory_path, filename)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            print(f"All files in {directory_path} have been removed.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

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
            metadata['sequence']= allLines.replace("n", "")
            totalSeqs+=1
            seqs.append(metadata)
    return seqs

def load_cannonical_fasta(pathTo):
    f=open(pathTo)
    x=f.readline()[:-1]
    seqs= []
    totalSeqs=0
    while x!='':    
        if ">" in x:
            metadata=x[1:-1]
            allLines = ''
            x=f.readline()[:-1]
            while x!='':
                if ">" in x:
                    break
                allLines+=x
                x=f.readline()[:-1]
            sequence= allLines.replace("n", "")
            totalSeqs+=1
            seqs.append(metadata)
    return sequence, metadata

def load_aligned_fasta(pathTo):
    f=open(pathTo)
    x=f.readline()[:-1]
    seqs= []
    totalSeqs=0
    while x!='':    
        if ">" in x:
            influenza=x[1:].split("|")
            if len(influenza)<6:
                print("cannonical model skipped")
                x=f.readline()[:-1]
                continue
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
            metadata['sequence']= allLines.replace("n", "")
            totalSeqs+=1
            seqs.append(metadata)
        else:
            x=f.readline()[:-1]
    return seqs

def load_aligned_protein_fasta(pathTo):
    f=open(pathTo)
    x=f.readline()[:-1]
    seqs= []
    totalSeqs=0
    while x!='':    
        if ">" in x:
            influenza=x[1:].split("_")
            metadata={}
            metadata["ID1"]=influenza[0]
            allLines = ''
            x=f.readline()[:-1]
            while x!='':
                if ">" in x:
                    break
                allLines+=x
                x=f.readline()[:-1]
            metadata['sequence']= allLines.replace("n", "")
            totalSeqs+=1
            seqs.append(metadata)
        else:
            x=f.readline()[:-1]
    return seqs


def load_metadata(pathTo):
    metadata_file=pd.read_csv(pathTo,dtype=str)
    return metadata_file

def write_host_annotated_fasta(seqs, output_location):
    save1= open(output_location, "w+")
    for item in seqs:
        item2=list(item)
        header= (">"+
        item2[0] + "|" +
        item2[1]+ "|" +
        item2[2]+ "|" +
        item2[3]+ "|" +
        item2[4]+ "|" +
        item2[5]+ "|" +
        str(item2[8] or "none") + "\n")
        save1.write(header)
        save1.write(item2[6]+ "\n")
    save1.close()

def write_host_annotated_fasta_pos_7(seqs, output_location):
    save1= open(output_location, "w+")
    for item in seqs:
        item2=list(item)
        header= (">"+
        item2[0] + "_" +
        item2[1]+ "_" +
        item2[2]+ "_" +
        item2[3]+ "_" +
        item2[4] + "\n")
        if len(item2[5])==0:
            print("item2[5] is empty")
            print(item2[5])
            exit()
        save1.write(header)
        save1.write(item2[5]+ "\n")
    save1.close()


def write_host_annotated_fasta_prime_cds(seqs, output_location, max_length=10,set_illegal_DNA_characters=False):
    save1= open(output_location, "w+")
    i=0
    for item in seqs:
        item2=list(item)
        header= (">"+
        item2[0] + "_" +
        item2[1]+ "_" +
        item2[2]+ "_" +
        item2[3]+ "_" +
        item2[4] + "\n")
        string_size =len(item2[5])
        illegal_DNA_characters = set((item2[5]).upper()) - set("ATCG")
        #search if there are any characters that are not A, T, C, G
        do_write=True
        if illegal_DNA_characters and set_illegal_DNA_characters:
            do_write=False
            print(illegal_DNA_characters)
            print(item2[0],item2[5])
        if string_size!=0 and string_size<max_length and do_write: 
            save1.write(header)
            save1.write(item2[5]+ "\n")
    save1.close()

def write_minimap_fasta_for_mafft(seqs, output_location):
    save1= open(output_location, "w+")
    for item in seqs:
        item2=list(item)
        header= (">"+
        item2[0] + "|" +
        item2[1]+ "|" +
        item2[2]+ "|" +
        item2[3]+ "|" +
        item2[4]+ "|" +
        item2[5]+ "\n" )
        save1.write(header)
        #more insert a break every 60 characters
        sequence=item['sequence']
        sequence= sequence.replace(".", "-").replace("*", "-")
        for i in range(0, len(sequence), 60):
            save1.write(sequence[i:i+60]+"\n")
        # save1.write(item2[6]+ "\n")
    save1.close()

#database tasks
def create_gisaid_DB2():
    conn=sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE sequences (
                ID1 text,
                geneSymbol text,
                strainName text,
                isolateId text,
                INSDC text,
                type text,
                sequence text,
                subtype text,
                species text,
                year text
                )""")
    conn.commit()
    c.execute("""CREATE TABLE metadata (
                isolateId text,
                host text
                )""")
    c.execute("""CREATE TABLE processed (
              ID1 text,
              processedCategory text
                )""")
    conn.commit()

    conn.close()

def create_gisaid_DB3():
    conn=sqlite3.connect('db_management/gisaid3.db')
    c = conn.cursor()

    c.execute("""CREATE TABLE IF NOT EXISTS metadata (
                ID1 text,
                isolateId text,
                geneSymbol text,
                subtype text,
                species text,
                host text
                )""")
    commands = [
    """CREATE TABLE HA (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE NA (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE MP (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE NP (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE NS (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE PA (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE PB1 (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE PB2 (
        ID1 text,
        sequence text
    )""",
    """CREATE TABLE discarded_sequences (
        ID1 text,
        sequence text,
        geneSymbol text
    )""",
    ]
    for command in commands:
        c.execute(command)
    conn.commit()
    conn.close()

def populate_db3_with_sequences(cannonical_max_lengths):
    conn3=sqlite3.connect('db_management/gisaid3.db')
    c3 = conn3.cursor()
    conn2=sqlite3.connect('db_management/gisaid2.db')
    c2 = conn2.cursor()
    #drop metadata table if exists
    c3.execute("DROP TABLE IF EXISTS metadata")
    conn3.commit()
    #create metadata table if not exists
    c3.execute("""CREATE TABLE metadata (
                ID1 text,
                isolateId text,
                geneSymbol text,
                subtype text,
                species text,
                host text
                )""")
    conn3.commit()
    # populate metadata table
    query = """
                SELECT
                    s.ID1,
                    s.isolateId,
                    s.geneSymbol,
                    s.subtype,
                    s.species,
                    p.host
                FROM
                    metadata p
                JOIN
                    sequences s ON p.isolateId = s.isolateId
                """
    c2.execute(query)
    response = c2.fetchall()
    to_db = [(item[0], item[1], item[2], item[3], item[4], item[5]) for item in response]
    c3.executemany("INSERT INTO metadata (ID1, isolateId, geneSymbol, subtype, species, host) VALUES (?, ?, ?, ?, ?, ?);", to_db)
    conn3.commit()

    # # populate gene tables
    for gene in cannonical_max_lengths:
        gene_symbol = gene["segment"]
        max_length = gene["maxSize"]
        min_length = gene["minSize"]
        query = f"""
                SELECT
                    s.ID1,
                    s.sequence
                FROM
                    sequences s
                WHERE
                    s.geneSymbol = ?
                """
        c2.execute(query, (gene_symbol,))
        response = c2.fetchall()
        for item in response:
            illegal_DNA_characters = set((item[1]).upper()) - set("ATCG")
            if len(item[1]) > 0 and len(item[1]) <= max_length and len(item[1]) >=min_length and not illegal_DNA_characters:
                to_db = [(item[0], item[1])]
                c3.executemany(f"INSERT INTO {gene_symbol} (ID1, sequence) VALUES (?, ?);", to_db)
            else:
                to_db = [(item[0], item[1], gene_symbol)]
                c3.executemany("INSERT INTO discarded_sequences (ID1, sequence, geneSymbol) VALUES (?, ?, ?);", to_db)
            conn3.commit()
    # make a new sequence table called clean_sequences, where no illegl characters are included
    # #drop table if exists
    c3.execute("DROP TABLE IF EXISTS clean_metadata")
    #create table if not exists
    # c3.execute("CREATE TABLE clean_metadata AS SELECT * FROM metadata m WHERE m.ID1 NOT IN (SELECT ID1 FROM discarded_sequences );")
    c3.execute("""
    CREATE TABLE clean_metadata AS 
    SELECT m.*
    FROM metadata m
    JOIN (
        SELECT isolateId, geneSymbol, MIN(ID1) as minID1
        FROM metadata
        WHERE ID1 NOT IN (SELECT ID1 FROM discarded_sequences)
        GROUP BY isolateId, geneSymbol
    ) AS subquery ON m.isolateId = subquery.isolateId AND m.geneSymbol = subquery.geneSymbol AND m.ID1 = subquery.minID1;
    """)


    conn3.close()
    conn2.close()

def calculate_mean_sequence_length():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    query = """
            SELECT subtype, AVG(LENGTH(sequence))
            FROM sequences
            GROUP BY species AND geneSymbol AND subtype
            """
    c.execute(query)
    result = c.fetchall()

    conn.close()

    return result


def create_table_with_translated_sequences_DB2():
    conn=sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE translated_sequences (
                ID1 text,
                sequence text,
                translated_sequence text
                )""")
    conn.commit()
    conn.close()

def create_table_paired_data_DB2():
    conn=sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE paired_data (
                isolateId text,
                id1_pa text,
                id1_pb1 text
                )""")
    conn.commit()
    conn.close()

def write_metadata_to_db2():
    metadata= load_metadata("output/simple_metadata.csv")
    

    print("loaded metadata in memory")
    conn=sqlite3.connect('db_management/gisaid2.db')
    c= conn.cursor()

    to_db = [(metadata['Isolate_Id'].values[i], metadata['Host'].values[i]) for i in range(0,len(metadata)-1)]


    c.executemany("INSERT INTO metadata (isolateId, host) VALUES (?, ?);", to_db)
    conn.commit()
    print("written the metadata db")
    conn.close()

def write_merged_sequences_to_db2():
    conn=sqlite3.connect('db_management/gisaid2.db')
    c= conn.cursor()
    response=load_fasta("output/merged_sequences.fasta")
    # response=load_fasta("downloads copy/1972-06-011972-06-30.fasta")
    sequences=response
    del response
    print("loaded fasta in memory")

    for sequence in sequences:
        # print(sequence)

        toProcess=sequence['type']
        splitted_list=toProcess.split("_")
        subtype= splitted_list[-1]
        sequence['subtype']=subtype
        sequence['species']=splitted_list[0]
        sequence['year']=(sequence['Strain Name']).split("/")[-1]
      
        attrib_names = ", ".join(sequence.keys())
        attrib_values = ", ".join("?"*len(sequence.keys()))
        attrib_names= "ID1, geneSymbol, strainName, isolateId ,INSDC, type, sequence, subtype, species, year" 
        sql = f"INSERT INTO sequences ({attrib_names}) VALUES ({attrib_values})"

        c.execute(sql, list(sequence.values()))

        # c.execute(sql, list(sequence.values()))
    conn.commit()
    conn.close()

    print("written the sequence db")
        
def write_processed_data_to_db2(file_route, method):
    conn=sqlite3.connect('db_management/gisaid2.db')
    c= conn.cursor()
    sequences=load_aligned_fasta(file_route)
    print("loaded fasta in memory")
    to_db = [(sequence['ID1'], method) for sequence in sequences]
    c.executemany("INSERT INTO processed (ID1, processedCategory) VALUES (?, ?);", to_db)
    conn.commit()
    print("written the processed db")
    conn.close()

def write_csvs_of_paired_data_PA():
    species= ["A","B"]
    gene_symbols=["PA"]
    outputDir="paired_data/input_csvs_for_mafft/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    for input_gene in gene_symbols:
        for input_species in species:
            print(input_gene, input_species)
            query = """
            SELECT
                s.ID1,
                s.isolateId,
                s.geneSymbol,
                s.subtype,
                s.species,
                t.translated_sequence
            FROM
                paired_data p
            JOIN
                sequences s ON p.id1_pa = s.ID1
            JOIN
                translated_sequences t ON p.id1_pa = t.ID1
            WHERE
                s.geneSymbol = ?
                AND s.species = ?
            """
            c.execute(query, (input_gene, input_species))
            response = c.fetchall()
            
            # Assuming write_host_annotated_fasta is a function that writes to a file
            write_host_annotated_fasta_pos_7(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta")


    conn.close()

def write_csvs_of_paired_data_PB1():
    species= ["A","B"]
    gene_symbols=["PB1"]
    outputDir="paired_data/input_csvs_for_mafft/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    for input_gene in gene_symbols:
        for input_species in species:
            print(input_gene, input_species)
            query = """
            SELECT
                s.ID1,
                s.isolateId,
                s.geneSymbol,
                s.subtype,
                s.species,
                t.translated_sequence
            FROM
                paired_data p
            JOIN
                sequences s ON p.id1_pb1 = s.ID1
            JOIN
                translated_sequences t ON p.id1_pb1 = t.ID1
            WHERE
                s.geneSymbol = ?
                AND s.species = ?
            """
            c.execute(query, (input_gene, input_species))
            response = c.fetchall()
            
            # Assuming write_host_annotated_fasta is a function that writes to a file
            write_host_annotated_fasta_pos_7(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta")


    conn.close()

def write_csvs_of_prime5_data():
    species= ["A","B"]
    gene_symbols=["PB1", "PA"]
    #if directory does not exist, create it
    os.makedirs("paired_data/input_csvs_for_mafft_prime5/", exist_ok=True)
    outputDir="paired_data/input_csvs_for_mafft_prime5/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    for input_gene in gene_symbols:
        if input_gene=="PB1":
            for input_species in species:
                print(input_gene, input_species)
                query = """
                SELECT
                    s.ID1,
                    s.isolateId,
                    s.geneSymbol,
                    s.subtype,
                    s.species,
                    t.prime5
                FROM
                    paired_data p
                JOIN
                    sequences s ON p.id1_pb1 = s.ID1
                JOIN
                    prime5 t ON p.id1_pb1 = t.ID1
                WHERE
                    s.geneSymbol = ?
                    AND s.species = ?
                """
                c.execute(query, (input_gene, input_species))
                response = c.fetchall()
                
                # Assuming write_host_annotated_fasta is a function that writes to a file
                write_host_annotated_fasta_prime_cds(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta",set_illegal_DNA_characters=True)
        elif input_gene=="PA":
            for input_species in species:
                print(input_gene, input_species)
                query = """
                SELECT
                    s.ID1,
                    s.isolateId,
                    s.geneSymbol,
                    s.subtype,
                    s.species,
                    t.prime5
                FROM
                    paired_data p
                JOIN
                    sequences s ON p.id1_pa = s.ID1
                JOIN
                    prime5 t ON p.id1_pa = t.ID1
                WHERE
                    s.geneSymbol = ?
                    AND s.species = ?
                """
                c.execute(query, (input_gene, input_species))
                response = c.fetchall()
                
                # Assuming write_host_annotated_fasta is a function that writes to a file
                write_host_annotated_fasta_prime_cds(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta",set_illegal_DNA_characters=True)
        

    conn.close()

def write_csvs_of_prime3_data():
    species= ["A","B"]
    gene_symbols=["PB1", "PA"]
    os.makedirs("paired_data/input_csvs_for_mafft_prime3/", exist_ok=True)

    outputDir="paired_data/input_csvs_for_mafft_prime3/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    for input_gene in gene_symbols:
        if input_gene=="PB1":
            for input_species in species:
                print(input_gene, input_species)
                query = """
                SELECT
                    s.ID1,
                    s.isolateId,
                    s.geneSymbol,
                    s.subtype,
                    s.species,
                    t.prime3
                FROM
                    paired_data p
                JOIN
                    sequences s ON p.id1_pb1 = s.ID1
                JOIN
                    prime3 t ON p.id1_pb1 = t.ID1
                WHERE
                    s.geneSymbol = ?
                    AND s.species = ?
                """
                c.execute(query, (input_gene, input_species))
                response = c.fetchall()
                
                # Assuming write_host_annotated_fasta is a function that writes to a file
                write_host_annotated_fasta_prime_cds(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta", max_length=45, set_illegal_DNA_characters=True)
        elif input_gene=="PA":
            for input_species in species:
                print(input_gene, input_species)
                query = """
                SELECT
                    s.ID1,
                    s.isolateId,
                    s.geneSymbol,
                    s.subtype,
                    s.species,
                    t.prime3
                FROM
                    paired_data p
                JOIN
                    sequences s ON p.id1_pa = s.ID1
                JOIN
                    prime3 t ON p.id1_pa = t.ID1
                WHERE
                    s.geneSymbol = ?
                    AND s.species = ?
                """
                c.execute(query, (input_gene, input_species))
                response = c.fetchall()
                
                # Assuming write_host_annotated_fasta is a function that writes to a file
                write_host_annotated_fasta_prime_cds(response, outputDir + f"{input_gene}_{input_species}_sequences.fasta", max_length=45, set_illegal_DNA_characters=True)
        

    conn.close()


def add_aligned_files_to_db2(pathToFastas, processedCategoryTool):
    aligned_files=[file for file in os.listdir(pathToFastas) if "fasta"  in file]
    for file in aligned_files:
        write_processed_data_to_db2(pathToFastas+file, processedCategoryTool)

def add_aligned_protein_sequences_to_db2(pathToFastas, rewrite_db=False):
    conn=sqlite3.connect('db_management/gisaid2.db')
    c= conn.cursor()
    if rewrite_db:
        c.execute("DROP TABLE IF EXISTS aligned_protein_sequences")
        c.execute("""CREATE TABLE aligned_protein_sequences (
                    ID1 text,
                    alignedsequence text
                        )""")
    elif not rewrite_db:
        #check if table exists 
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='aligned_protein_sequences';")
        if c.fetchone() is None:
            c.execute("""CREATE TABLE aligned_protein_sequences (
                            ID1 text,
                            alignedsequence text,
                            )""")
        else:
            print("table already exists")
            exit()
    aligned_files=[file for file in os.listdir(pathToFastas) if "fasta"  in file]
    for file in aligned_files:
        sequences=load_aligned_protein_fasta(pathToFastas+file)
        print("loaded fasta in memory")
        conn.commit()
        for sequence in sequences:
            to_db = [(sequence['ID1'], sequence['sequence'])]
            c.executemany("INSERT INTO aligned_protein_sequences (ID1, alignedsequence) VALUES (?, ?);", to_db)
        conn.commit()
        print("written the translated sequence db")
    conn.close()

def generate_fasta_file_by_gene_host_species_and_subtype2():
    gene_symbols=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    input_host="Human"
    input_species="A"
    subtypes = ["H1N1", "H3N2", "H7N9", "H0N0"]
    conn=sqlite3.connect('db_management/gisaid2.db')
    c= conn.cursor()
    outputDir="humanSubtypeOutput/"
    query = """
    SELECT *
    FROM sequences
    LEFT JOIN metadata USING (isolateId)
    WHERE species = ? AND host = ? AND geneSymbol = ? and subtype = ?
    """
    for input_gene in gene_symbols:
        for input_subtype in subtypes:
            c.execute(query, (input_species,input_host,  input_gene, input_subtype))
            response = c.fetchall()
            write_host_annotated_fasta(response,outputDir+input_gene+"_"+input_subtype+"_"+input_species+"_"+input_host+"_sequences.fasta" )

    conn.close()

def generate_fasta_file_by_gene_host_species_and_subtype2_outer():
    gene_symbols = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    input_host = "Human"
    input_species = "A"
    excluded_subtypes = ["H1N1", "H3N2", "H7N9", "H0N0"]
    outer_symbol="HXNX"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    outputDir="humanSubtypeOutput/"
    query = """
        SELECT *
        FROM sequences
        LEFT JOIN metadata USING (isolateId)
        WHERE species = ? AND host = ? AND geneSymbol = ?
        AND subtype NOT IN (?, ?, ?, ?)
    """
    non_listed_subtypes = set()
    for input_gene in gene_symbols:
        c.execute(query, (input_species, input_host, input_gene, *excluded_subtypes))
        response = c.fetchall()
        write_host_annotated_fasta(response,outputDir+input_gene+"_"+outer_symbol+"_"+input_species+"_"+input_host+"_sequences.fasta" )


    conn.close()
    return non_listed_subtypes

def find_non_aligned_files_and_write_to_csv1():
    gene_symbols=[ "HA", "NA"]
    input_host="Human"
    input_species="A"
    subtypes = [ "H3N2"]
    outputDir="humanSubtypeOutput_non_aligned_1/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query="""
    SELECT *
    FROM sequences 
    LEFT JOIN metadata USING (isolateId) LEFT JOIN processed USING (ID1)
    WHERE processedCategory IS NULL and species = ? AND host = ? AND geneSymbol = ? and subtype = ?;
    """
    for input_gene in gene_symbols:
        for input_subtype in subtypes:
            c.execute(query, (input_species,input_host,  input_gene, input_subtype))
            response = c.fetchall()
            write_host_annotated_fasta(response,outputDir+input_gene+"_"+input_subtype+"_"+input_species+"_"+input_host+"_sequences.fasta" )

def find_non_aligned_files_and_write_to_csv2():
    gene_symbols = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    input_host="Human"
    input_species="A"
    subtypes = ["H1N1", "H3N2", "H7N9", "H0N0"]
    outputDir="humanSubtypeOutput_non_aligned_2/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query="""
    SELECT *
    FROM sequences 
    LEFT JOIN metadata USING (isolateId) LEFT JOIN processed USING (ID1)
    WHERE processedCategory IS NULL and species = ? AND host = ? AND geneSymbol = ? and subtype = ?;
    """
    for input_gene in gene_symbols:
        for input_subtype in subtypes:
            c.execute(query, (input_species,input_host,  input_gene, input_subtype))
            response = c.fetchall()
            write_host_annotated_fasta(response,outputDir+input_gene+"_"+input_subtype+"_"+input_species+"_"+input_host+"_sequences.fasta" )


def add_matching_isolateId_to_db3():
    conn=sqlite3.connect('db_management/gisaid3.db')
    c= conn.cursor()
    #drop table if exists
    c.execute("DROP TABLE IF EXISTS matched_genesymbols")
    conn.commit()
    query = """ CREATE TABLE IF NOT EXISTS matched_genesymbols (
                isolateId text
                )"""
    c.execute(query)
    conn.commit()
    query = """INSERT INTO matched_genesymbols (isolateId)
    SELECT subquery.isolateId
    FROM (
        SELECT m.isolateId
        FROM clean_metadata m
        WHERE m.geneSymbol IN ('HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2')
        GROUP BY m.isolateId
        HAVING COUNT(DISTINCT m.geneSymbol) = 8
    ) AS subquery;"""


    c.execute(query)
    conn.commit()


def print_species_and_host_in_paired_data():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    # Query to get total isolates per species
    total_isolates_query = """
    SELECT
        s.species,
        COUNT(DISTINCT p.isolateId) as total_isolates
    FROM
        paired_data p
    JOIN
        sequences s ON p.isolateId = s.isolateId
    GROUP BY
        s.species
    """

    c.execute(total_isolates_query)
    total_isolates_result = dict(c.fetchall())
    # Query to get detailed breakdown per host
    query = """
    SELECT
        s.species,
        m.host,
        COUNT(DISTINCT p.isolateId) as total_isolates
    FROM
        paired_data p
    JOIN
        sequences s ON p.isolateId = s.isolateId
    JOIN
        metadata m ON p.isolateId = m.isolateId
    GROUP BY
        s.species, m.host
    ORDER BY
        s.species, total_isolates DESC;
    """
    c.execute(query)
    result = c.fetchall()

    current_species = None

    for row in result:
        species, host, total_isolates = row
        if species != current_species:
            if current_species is not None:
                print(f"  Total Isolates for {current_species}: {total_isolates_result[current_species]}")
            print(f"Species: {species}")
            current_species = species
        total_isolates_species = total_isolates_result.get(species, 0)
        percentage = (total_isolates / total_isolates_species) * 100 if total_isolates_species > 0 else 0
        print(f"    Host: {host}, Percentage: {percentage:.2f}%")

    if current_species is not None:
        print(f"  Total Isolates for {current_species}: {total_isolates_result[current_species]}")

    conn.close()

def print_species_and_subtype_in_paired_data():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    # Query to get total isolates per species
    total_isolates_query = """
    SELECT
        s.species,
        COUNT(DISTINCT p.isolateId) as total_isolates
    FROM
        paired_data p
    JOIN
        sequences s ON p.isolateId = s.isolateId
    GROUP BY
        s.species
    """

    c.execute(total_isolates_query)
    total_isolates_result = dict(c.fetchall())

    # Query to get detailed breakdown per subtype
    query = """
    SELECT
        s.species,
        s.subtype,
        COUNT(DISTINCT p.isolateId) as total_isolates
    FROM
        paired_data p
    JOIN
        sequences s ON p.isolateId = s.isolateId
    GROUP BY
        s.species, s.subtype
    ORDER BY
        s.species, total_isolates DESC;
    """

    c.execute(query)
    result = c.fetchall()

    current_species = None

    for row in result:
        species, subtype, total_isolates = row
        if species != current_species:
            if current_species is not None:
                print(f"  Total Isolates for {current_species}: {total_isolates_result[current_species]}")
            print(f"Species: {species}")
            current_species = species
        total_isolates_species = total_isolates_result.get(species, 0)
        percentage = (total_isolates / total_isolates_species) * 100 if total_isolates_species > 0 else 0
        print(f"    Subtype: {subtype}, Percentage: {percentage:.2f}%")

    if current_species is not None:
        print(f"  Total Isolates for {current_species}: {total_isolates_result[current_species]}")

    conn.close()



def find_non_aligned_files_and_write_to_csv2_outer():
    gene_symbols=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    input_host="Human"
    input_species="A"
    excluded_subtypes = ["H1N1", "H3N2", "H7N9", "H0N0"]
    outer_symbol="HXNX"
    outputDir="humanSubtypeOutput_non_aligned_2/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query="""
    SELECT *
    FROM sequences 
    LEFT JOIN metadata USING (isolateId) LEFT JOIN processed USING (ID1)
    WHERE processedCategory IS NULL and species = ? AND host = ? AND geneSymbol = ? AND subtype NOT IN (?, ?, ?, ?);
    """
    for input_gene in gene_symbols:
        c.execute(query, (input_species, input_host, input_gene, *excluded_subtypes))
        response = c.fetchall()
        write_host_annotated_fasta(response,outputDir+input_gene+"_"+outer_symbol+"_"+input_species+"_"+input_host+"_sequences.fasta" )

def find_non_aligned_models_and_write_to_csv_aft_mafft():
    gene_symbols=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    input_host="Human"
    input_species="A"
    outputDir="leftover_after_mafft/"
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query="""
    SELECT *
    FROM sequences 
    LEFT JOIN metadata USING (isolateId) LEFT JOIN processed USING (ID1)
    WHERE processedCategory IS NULL and species = ? AND host = ? AND geneSymbol = ?;
    """
    for input_gene in gene_symbols:
        c.execute(query, (input_species, input_host, input_gene))
        response = c.fetchall()
        write_host_annotated_fasta(response,outputDir+input_gene+"_"+input_species+"_"+input_host+"_sequences.fasta" )
    
def print_processed_categories():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query = """
    SELECT
        processedCategory,
        COUNT(*) AS categoryCount
    FROM
        processed
    GROUP BY
        processedCategory
    UNION ALL
    SELECT
        'Null Processed Category',
        COUNT(*) AS categoryCount
    FROM
        processed
    WHERE
        processedCategory IS NULL;
    """
    c.execute(query)
    results = c.fetchall()
    
    for category, count in results:
        print(f"Processed Category: {category}, Count: {count}")

def print_length_of_all_elements_in_table():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    # Execute an SQL query to calculate the length of each sequence and get the total count
    query = """
    SELECT COUNT(*) as total_count, COUNT(DISTINCT geneSymbol) as distinct_geneSymbol
    FROM sequences
    LEFT JOIN metadata USING (isolateId)
    LEFT JOIN processed USING (ID1)
    WHERE species = "A" AND host = "Human";
    """
    c.execute(query)

    row = c.fetchone()
    print(row)

    total_count = row[0]
    distinct_geneSymbol = row[1]

    print(f"Total Elements: {total_count}")
    print(f"Distinct Gene Symbols: {distinct_geneSymbol}")

    # Query to get distinct gene symbols and their counts
    distinct_gene_symbols_query = """
    SELECT geneSymbol, COUNT(*) as symbol_count
    FROM sequences
    GROUP BY geneSymbol;
    """
    c.execute(distinct_gene_symbols_query)
    distinct_gene_symbols = c.fetchall()

    for symbol, count in distinct_gene_symbols:
        print(f"Gene Symbol: {symbol}, Count: {count}")
    # Query to get the total count of all elements
    all_query = """
    SELECT COUNT(*) as all_count
    FROM sequences LEFT JOIN metadata USING (isolateId) 
    """
    # where host=="Human" and species=="A";
    c.execute(all_query)
    all_symbols = c.fetchall()
    print(all_symbols)
    # print(f"Gene Symbol: ALL, Count: {all_symbols[0][0]}")


    conn.close()


def print_count_of_all_elements_in_table(table):
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    #sql query to get the total count of all elements
    query = f"SELECT COUNT(*) as all_count FROM {table}"
    # Execute the query without placeholders since table names cannot be parameterized
    c.execute(query)
    all_elements_count = c.fetchone()  # Use fetchone() since it's a single value
    print(all_elements_count[0])  # Access the count value from the result tuple

    # Close the database connection
    conn.close()

def search_cannonical_model_by_string(fastas_directory, search_string):
    #find all the fasta files in the directory
    fastas = os.listdir(fastas_directory)
    # Iterate over the FASTA files and print the matching sequences
    for fasta in fastas:
        if fasta.endswith(".fasta"):
            sequence, metadata=load_cannonical_fasta(fastas_directory+fasta)
            #for sequence, search how many times the strings appear
            total_counts = 0
            total_overlapping_counts = 0
            print(metadata)
            for motif in search_string:
                total_counts += sequence.count(motif)
                for i in range(len(sequence) - len(motif) + 1):
                    if sequence[i:i+len(motif)] == motif:
                        total_overlapping_counts += 1

            
            print("    total",total_counts)
            print("    total_overlapping",total_overlapping_counts)
    


def search_sequence_by_string(search_string):
    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    # Execute the SQL query to search for the string in the sequence column
    query = """
    SELECT geneSymbol, COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE sequence LIKE ? AND host = "Human" and species= "A" GROUP BY geneSymbol;
    """

    
    c.execute(query, ('%' + search_string + '%',))

    # Fetch and return the matching rows
    matching_count = c.fetchall()

    # Query for total counts
    total_count_query = """
    SELECT geneSymbol, COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = 'Human' and species= 'A' GROUP BY geneSymbol;
    """
    c.execute(total_count_query)
    total_count = c.fetchall()


    for (symbol, count) in matching_count:
        for (symbol2, total_counts) in total_count:
            if symbol==symbol2:
                percentage = (count / total_counts) * 100
                print(f"{symbol}, Count: {count}, ({percentage:.2f}%)")
    
    total_total = """
    SELECT COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = "Human" and species= "A" ;
    """
    c.execute(total_total)
    total_total = c.fetchall()
    total_with_string= """
    SELECT COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId)  WHERE sequence LIKE ? and host = "Human" and species= "A";
    """
    c.execute(total_with_string, ('%' + search_string + '%',))
    total_with_string = c.fetchall()
    percentage = (total_with_string[0][0] / total_total[0][0]) * 100
    print(f"Present in: {total_with_string[0][0]}, Total: {total_total[0][0]}, Percentage:  {percentage:.2f}%")


    # Close the database connection
    conn.close()

def search_sequence_by_multiple_strings(search_strings):
    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()

    # Format the strings for the SQL query
    formatted_strings = tuple('%' + string + '%' for string in search_strings)
    
    #make the query accept multiple search strings
    query_conditions = " OR ".join(["sequence LIKE ?" for _ in search_strings])
    # Execute the SQL query to search for the string in the sequence column
    matching_count_query = f"""
    SELECT geneSymbol, COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE ({query_conditions}) AND host = 'Human' and species= 'A' GROUP BY geneSymbol;
    """
    c.execute(matching_count_query, formatted_strings)
    matching_count = c.fetchall()

    # Query for total counts
    total_count_query = """
    SELECT geneSymbol, COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = 'Human' and species= 'A' GROUP BY geneSymbol;
    """
    c.execute(total_count_query)
    total_count = c.fetchall()


    # Calculate and print the percentages for matching counts
    for (symbol, count) in matching_count:
        for (symbol2, total_counts) in total_count:
            if symbol == symbol2:
                percentage = (count / total_counts) * 100
                print(f"{symbol}, Count: {count}, ({percentage:.2f}%)")
    # Query for total number of sequences
    total_total = """
    SELECT COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = "Human" and species= "A" ;
    """
    c.execute(total_total)
    total_total = c.fetchone()[0]


    # Query for total number of sequences with any of the search strings
    total_with_string_query = f"""
    SELECT COUNT(*) as symbol_count
    FROM sequences LEFT JOIN metadata USING (isolateId) WHERE ({query_conditions}) and host = 'Human' and species= 'A';
    """
    c.execute(total_with_string_query, formatted_strings)
    total_with_string = c.fetchone()[0]
    # print("total_count", total_count)
    # print("matching_count", matching_count)
    # print("total_total", total_total)
    # print("total_with_string",total_with_string)

    print("query conditions", query_conditions)


    percentage = (total_with_string / total_total) * 100
    print(f"Present in: {total_with_string}, Total: {total_total}, Percentage:  {percentage:.2f}%")


    # Close the database connection
    conn.close()

def search_sequence_by_multiple_strings_db3(search_strings, _gene_symbols="all"):
    if _gene_symbols=="all":
        gene_symbols= ["HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2"]
    else:
        gene_symbols=_gene_symbols
    


    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid3.db')
    c = conn.cursor()
    #drop table if exists
    c.execute("DROP TABLE IF EXISTS non_matched_sequences_human")
    conn.commit()
    #create a table for the non matched sequences
    c.execute("""CREATE TABLE IF NOT EXISTS non_matched_sequences_human (
                ID1 text,
                sequence text
                )""")
    conn.commit()

    #for each gene symbol, search for the strings
    #each gene symbol will query a different table
    #each table only has ID1 and sequence
    for gene_symbol in gene_symbols:
        # Format the strings for the SQL query
        formatted_strings = tuple('%' + string + '%' for string in search_strings)
        
        #make the query accept multiple search strings
        query_conditions = " OR ".join(["sequence LIKE ?" for _ in search_strings])
        # Execute the SQL query to search for the string in the sequence column
        matching_count_query = f"""
        SELECT COUNT(*) as symbol_count
        FROM {gene_symbol} g 
        JOIN clean_metadata m ON g.ID1 = m.ID1
        WHERE ({query_conditions}) AND m.host = 'Human';
        """
        c.execute(matching_count_query, formatted_strings)
        matching_count = c.fetchone()[0]
        # Query for total counts
        total_count_query = f"""
        SELECT COUNT(*) as symbol_count
        FROM {gene_symbol} g 
        JOIN clean_metadata m ON g.ID1 = m.ID1 
        WHERE m.host = 'Human';
        """
        c.execute(total_count_query)
        total_count = c.fetchone()[0]

        # Calculate and print the percentages for matching counts
        percentage = (matching_count / total_count) * 100
        print(f"{gene_symbol}, Count: {matching_count}, ({percentage:.2f}%), Total: {total_count}")

        matching_count_query = f"""
        SELECT g.ID1, g.sequence as symbol_count
        FROM {gene_symbol} g 
        JOIN clean_metadata m ON g.ID1 = m.ID1
        WHERE NOT ({query_conditions}) AND m.host = 'Human';
        """
        c.execute(matching_count_query, formatted_strings)
        non_matching_count = c.fetchall()
        to_db = [(sequence[0], sequence[1]) for sequence in non_matching_count]
        c.executemany("INSERT INTO non_matched_sequences_human (ID1, sequence) VALUES (?, ?);", to_db)
        conn.commit()

    conn.close()


        
        







    # # Format the strings for the SQL query
    # formatted_strings = tuple('%' + string + '%' for string in search_strings)
    
    # #make the query accept multiple search strings
    # query_conditions = " OR ".join(["sequence LIKE ?" for _ in search_strings])
    # # Execute the SQL query to search for the string in the sequence column
    # matching_count_query = f"""
    # SELECT geneSymbol, COUNT(*) as symbol_count
    # FROM sequences LEFT JOIN metadata USING (isolateId) WHERE ({query_conditions}) AND host = 'Human' and species= 'A' GROUP BY geneSymbol;
    # # """
    # c.execute(matching_count_query, formatted_strings)
    # matching_count = c.fetchall()

    # # Query for total counts
    # total_count_query = """
    # SELECT geneSymbol, COUNT(*) as symbol_count
    # FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = 'Human' and species= 'A' GROUP BY geneSymbol;
    # """
    # c.execute(total_count_query)
    # total_count = c.fetchall()


    # # Calculate and print the percentages for matching counts
    # for (symbol, count) in matching_count:
    #     for (symbol2, total_counts) in total_count:
    #         if symbol == symbol2:
    #             percentage = (count / total_counts) * 100
    #             print(f"{symbol}, Count: {count}, ({percentage:.2f}%)")
    # # Query for total number of sequences
    # # total_total = """
    # # SELECT COUNT(*) as symbol_count
    # # FROM sequences LEFT JOIN metadata USING (isolateId) WHERE host = "Human" and species= "A" ;
    # # """
    # # c.execute(total_total)
    # # total_total = c.fetchone()[0]


    # # Query for total number of sequences with any of the search strings
    # total_with_string_query = f"""
    # SELECT COUNT(*) as symbol_count
    # FROM sequences LEFT JOIN metadata USING (isolateId) WHERE ({query_conditions}) and host = 'Human' and species= 'A';
    # """
    # c.execute(total_with_string_query, formatted_strings)
    # total_with_string = c.fetchone()[0]
    # # print("total_count", total_count)
    # # print("matching_count", matching_count)
    # # print("total_total", total_total)
    # # print("total_with_string",total_with_string)

    # print("query conditions", query_conditions)


    # percentage = (total_with_string / total_total) * 100
    # print(f"Present in: {total_with_string}, Total: {total_total}, Percentage:  {percentage:.2f}%")


    # Close the database connection
    # conn.close()

def search_sequence_by_multiple_strings_db3_onlyMatched_genesymbols(search_strings, _gene_symbol="all"):
    if _gene_symbol=="all":
        gene_symbols= ["HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2"]
    else:
        gene_symbols=[_gene_symbol]

    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid3.db')
    c = conn.cursor()

    #if table exists, delete it and create a new one
    c.execute("DROP TABLE IF EXISTS non_matched_sequences_human_matched_geneSymbol")
    #create a table for the non matched sequences
    c.execute("""CREATE TABLE IF NOT EXISTS non_matched_sequences_human_matched_geneSymbol (
                ID1 text,
                sequence text
                )""")
    conn.commit()

    #for each gene symbol, search for the strings
    #each gene symbol will query a different table
    #each table only has ID1 and sequence
    for gene_symbol in gene_symbols:
        # Format the strings for the SQL query
        formatted_strings = tuple('%' + string + '%' for string in search_strings)
        
        #make the query accept multiple search strings
        query_conditions = " OR ".join(["sequence LIKE ?" for _ in search_strings])
        # Execute the SQL query to search for the string in the sequence column
        matching_count_query = f"""
            SELECT COUNT(DISTINCT m.ID1) as symbol_count
            FROM {gene_symbol} g 
            JOIN clean_metadata m ON g.ID1 = m.ID1
            JOIN matched_genesymbols mg ON mg.isolateId = m.isolateId
            WHERE ({query_conditions}) AND m.host = 'Human' and m.species= 'A';
        """
        c.execute(matching_count_query, formatted_strings)
        matching_count = c.fetchone()[0]
        # Query for total counts
        total_count_query = f"""
            SELECT COUNT(DISTINCT m.ID1) as symbol_count
            FROM {gene_symbol} g 
            JOIN clean_metadata m ON g.ID1 = m.ID1 
            JOIN matched_genesymbols mg ON mg.isolateId = m.isolateId
            WHERE m.host = 'Human' and m.species= 'A';
        """
        c.execute(total_count_query)
        total_count = c.fetchone()[0]

        # Calculate and print the percentages for matching counts
        #prevent division by zero
        if total_count==0:
            total_count=1
        percentage = (matching_count / total_count) * 100
        print(f"{gene_symbol}, Count: {matching_count},Total Count: {total_count} ,({percentage:.2f}%)")

        matching_count_query = f"""
        SELECT DISTINCT g.ID1, g.sequence as symbol_count
        FROM {gene_symbol} g 
        JOIN clean_metadata m ON g.ID1 = m.ID1
        JOIN matched_genesymbols mg ON mg.isolateId = m.isolateId
        WHERE NOT ({query_conditions}) AND m.host = 'Human' and m.species= 'A';
        """
        c.execute(matching_count_query, formatted_strings)
        non_matching_count = c.fetchall()
        to_db = [(sequence[0], sequence[1]) for sequence in non_matching_count]
        print("amount of non matching sequences", len(to_db))
        c.executemany("INSERT INTO non_matched_sequences_human_matched_geneSymbol (ID1, sequence) VALUES (?, ?);", to_db)
        conn.commit()

    conn.close()

def search_multiple_ID1_per_geneSymbol_matched():
    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid3.db')
    c = conn.cursor()
    query = """
    SELECT m.geneSymbol, m.isolateId, COUNT(DISTINCT m.ID1) as symbol_count
    FROM clean_metadata m
    GROUP BY m.geneSymbol, m.isolateId
    HAVING COUNT(DISTINCT m.ID1) > 1;
    """
    c.execute(query)
    results = c.fetchall()
    print(results[0:10])
    #find EPI_ISL_10058 and HA genesymbol
    query = """
    SELECT m.geneSymbol, m.isolateId, m.ID1, g.sequence
    FROM clean_metadata m
    JOIN HA g ON g.ID1 = m.ID1
    WHERE m.isolateId = "EPI_ISL_10058" AND m.geneSymbol = "HA";
    """
    c.execute(query)
    results = c.fetchall()
    print(results)

    query = """
    SELECT *
    FROM discarded_sequences
    WHERE ID1 = '275378';
    """
    c.execute(query)
    results = c.fetchall()
    print(results)

    conn.close()

def search_isolateId_with_at_least_one_match():
    # Connect to the SQLite database
    conn = sqlite3.connect('db_management/gisaid3.db')
    c = conn.cursor()

    # Find the total number of distinct isolateIds in matched_genesymbols
    total_count_query = f"""
    SELECT COUNT(DISTINCT isolateId) as total_count
    FROM matched_genesymbols;
    """
    c.execute(total_count_query)
    total_count = c.fetchone()[0]

    total_recalculated_query = f"""
    SELECT COUNT(DISTINCT isolateId)
    FROM (
        SELECT m.isolateId
        FROM clean_metadata m
        WHERE m.geneSymbol IN ('HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2') 
        AND m.host = 'Human' 
        AND m.species= 'A'
        GROUP BY m.isolateId
        HAVING COUNT(DISTINCT m.geneSymbol) = 8
    ) as subquery;
    """

    c.execute(total_recalculated_query)
    total_recalculated_count = c.fetchone()[0]

    query = """
    SELECT COUNT(*)
    FROM (
        SELECT m.isolateId
        FROM non_matched_sequences_human_matched_geneSymbol n
        JOIN clean_metadata m ON n.ID1 = m.ID1
        WHERE m.geneSymbol IN ('HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2')
        GROUP BY m.isolateId
        HAVING COUNT(DISTINCT m.geneSymbol) = 8
    ) AS subquery;
    """
    c.execute(query)
    total_non_matched_count = c.fetchone()[0]
    

    # the total number of distinct isolateIds in non_matched_sequences_human_matched_geneSymbol
    total_non_match_query = f"""
    SELECT COUNT(DISTINCT mg.isolateId) as total_count
    FROM non_matched_sequences_human_matched_geneSymbol g 
    JOIN clean_metadata m ON g.ID1 = m.ID1
    JOIN matched_genesymbols mg ON mg.isolateId = m.isolateId
    WHERE m.host = 'Human' and m.species= 'A';
    """
    c.execute(total_non_match_query)
    total_non_match_count = c.fetchone()[0]

    print(f"Total count: {total_count}")
    print("total_recalculated_count", total_recalculated_count)
    print(f"Total non match count: {total_non_match_count}")
    print(f"Total non matched count: {total_non_matched_count}")
    # Close the database connection
    

    #total number of distinct isolateIds in non_matched_sequences
    total_non_All_query = f"""
    SELECT COUNT(DISTINCT ID1) as total_count
    FROM non_matched_sequences_human;
    """
    c.execute(total_non_All_query)
    total_non_All_count = c.fetchone()[0]

    total_non_matched_ALL_query = f"""
    SELECT COUNT(DISTINCT mg.isolateId) as total_count
    FROM non_matched_sequences_human g 
    JOIN clean_metadata m ON g.ID1 = m.ID1
    JOIN matched_genesymbols mg ON mg.isolateId = m.isolateId
    WHERE m.host = 'Human' and m.species= 'A';
    """
    c.execute(total_non_matched_ALL_query)
    total_non_matched_ALL_count = c.fetchone()[0]

    query = """
    SELECT COUNT(*)
    FROM (
        SELECT m.isolateId
        FROM non_matched_sequences_human n
        JOIN clean_metadata m ON n.ID1 = m.ID1
        WHERE m.geneSymbol IN ('HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2')
        GROUP BY m.isolateId
        HAVING COUNT(DISTINCT m.geneSymbol) = 8
    ) AS subquery;
    """
    c.execute(query)
    total_non_matched_count_ALL = c.fetchone()[0]

    print(f"Total DISTINCT isolateId on on_matched_sequences_human: {total_non_All_count}")
    print(f"Total complete genomes on non_matched_sequences_human: {total_non_matched_ALL_count}")
    print(f"Total complete genomes on non_matched_sequences_human_matched_geneSymbol: {total_non_matched_count_ALL}")
    conn.close()
    
def translate(sequence):
    # Define the mapping from codons to amino acids
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }

    # Find the first "ATG" codon as the start codon
    start_index = sequence.upper().find("ATG")
    if start_index == -1:
        print(sequence)
        print(" No start codon found")
        return ""  # No start codon found


    # Initialize the protein sequence and start translating from the "ATG"
    protein_sequence = ""
    i = start_index
    while i < len(sequence):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")  # Use "X" for unknown codons
        protein_sequence += amino_acid

        if amino_acid == "*":  # Stop translation at the first "*"
            break

        i += 3

    return protein_sequence

def extract_prime3_sequence(sequence):
    # Define the mapping from codons to amino acids
    codon_table = {
        "TAA": "*", "TAG": "*","TGA": "*"
    }

    # Find the first "ATG" codon as the start codon
    start_index = sequence.upper().find("ATG")
    if start_index == -1:
        print(sequence)
        print(" No start codon found")
        return ""  # No start codon found

    # Initialize the protein sequence and start translating from the "ATG"
 
    i = start_index
    prime3_sequence = ""
    while i < len(sequence):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")  # Use "X" for unknown codons
        if amino_acid == "*":  # Stop translation at the first "*"
            prime3_sequence=sequence[i+3:]
            break
        i += 3

    return prime3_sequence

def extract_CDS_sequence(sequence):
    # Define the mapping from codons to amino acids
    codon_table = {
        "TAA": "*", "TAG": "*","TGA": "*"
    }

    # Find the first "ATG" codon as the start codon
    start_index = sequence.upper().find("ATG")
    if start_index == -1:
        print(sequence)
        print(" No start codon found")
        return ""  # No start codon found

    # Initialize the protein sequence and start translating from the "ATG"
 
    i = start_index
    cds_sequence = ""
    while i < len(sequence):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")  # Use "X" for unknown codons
        if amino_acid == "*":  # Stop translation at the first "*"
            cds_sequence=sequence[start_index:i+3]
            break
        i += 3

    return cds_sequence

def extract_prime5_sequence(sequence):
    # Define the mapping from codons to amino acids
    # Find the first "ATG" codon as the start codon
    start_index = sequence.upper().find("ATG")
    if start_index == -1:
        print(sequence)
        print(" No start codon found")
        return ""  # No start codon found
    # from the start to the first "ATG", without including it, take all the nucleotides 
    prime5_sequence = ""
    i = 0
    while i < start_index:
        if sequence[i:i+3].upper()=="ATG":
            break
        prime5_sequence += sequence[i]
        i += 1
    return prime5_sequence

def add_translated_sequences_to_db2():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query = """
    SELECT ID1, sequence
    FROM sequences """
    c.execute(query)
    sequences = c.fetchall()
    to_db = []
    for sequence in sequences:
        translated_sequence=translate(sequence[1])
        to_db.append((sequence[0], sequence[1], translated_sequence)) 
    c.executemany("INSERT INTO translated_sequences (ID1, sequence, translated_sequence) VALUES (?, ?, ?);", to_db)
    conn.commit()
    print("written the translated_sequences db")
    conn.close()

def add_cds_sequences_to_db2():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""DROP TABLE IF EXISTS cds""")
    c.execute("""CREATE TABLE cds (
                ID1 text,
                cds text
                )""")
    conn.commit()
    query = """
    SELECT ID1, sequence
    FROM sequences """
    c.execute(query)
    sequences = c.fetchall()
    to_db = []
    print("extracting cds sequences")
    i=0
    total_len=len(sequences)
    for sequence in sequences:
        i+=1
        #print progress i of total_len
        sys.stdout.write("\r %i of %i" % (i, total_len))
        sys.stdout.flush()
        cds_sequence=extract_CDS_sequence(sequence[1])
        to_db.append((sequence[0], cds_sequence)) 
    c.executemany("INSERT INTO cds (ID1, cds) VALUES (?, ?);", to_db)
    conn.commit()
    print("written the cds_sequence db")
    conn.close()

def add_prime3_sequences_to_db2():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""DROP TABLE IF EXISTS prime3""")
    c.execute("""CREATE TABLE prime3 (
                ID1 text,
                prime3 text
                )""")
    conn.commit()
    query = """
    SELECT ID1, sequence
    FROM sequences """
    c.execute(query)
    sequences = c.fetchall()
    to_db = []
    print("extracting prime3 sequences")
    i=0
    total_len=len(sequences)
    for sequence in sequences:
        i+=1
        #print progress i of total_len
        sys.stdout.write("\r %i of %i" % (i, total_len))
        sys.stdout.flush()
        prime3_sequence=extract_prime3_sequence(sequence[1])
        to_db.append((sequence[0], prime3_sequence)) 
    c.executemany("INSERT INTO prime3 (ID1, prime3) VALUES (?, ?);", to_db)
    conn.commit()
    print("written the prime3_sequence db")
    conn.close()

def add_prime5_sequences_to_db2():
    #if table prime5 doesnt exist, create it; otherwise, drop it and create it again
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    c.execute("""DROP TABLE IF EXISTS prime5""")
    c.execute("""CREATE TABLE prime5 (
                ID1 text,
                prime5 text
                )""")
    conn.commit()
    #select all RNA sequences
    query = """
    SELECT ID1, sequence
    FROM sequences """
    c.execute(query)
    sequences = c.fetchall()
    to_db = []
    for sequence in sequences:
        prime5=extract_prime5_sequence(sequence[1])
        to_db.append((sequence[0],  prime5)) 
    c.executemany("INSERT INTO prime5 (ID1, prime5) VALUES (?, ?);", to_db)
    conn.commit()
    print("written the prime5 sequences db")
    conn.close()

def add_paired_data_to_db2():
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    # First execute the SELECT part
    select_query = """
    SELECT s1.isolateId, s1.ID1 AS id1_pa, s2.ID1 AS id1_pb1
    FROM sequences s1
    JOIN sequences s2 ON s1.isolateId = s2.isolateId
    JOIN translated_sequences t1 ON s1.ID1 = t1.ID1 AND s1.geneSymbol = 'PA' AND LENGTH(t1.translated_sequence) >= 714 AND LENGTH(t1.translated_sequence) < 750
    JOIN translated_sequences t2 ON s2.ID1 = t2.ID1 AND s2.geneSymbol = 'PB1' AND LENGTH(t2.translated_sequence) >= 15 AND LENGTH(t2.translated_sequence) < 770;
    """
    c.execute(select_query)
    selected_data = c.fetchall()


    # Then execute the INSERT INTO part
    insert_query = "INSERT INTO paired_data (isolateId, id1_pa, id1_pb1) VALUES (?, ?, ?)"
    c.executemany(insert_query, selected_data)

    # Commit the changes and close the connection
    conn.commit()
    print("Written the paired_data to the database.")
    conn.close()

def print_specific_column_from_table(table, column, whichone):
    conn = sqlite3.connect('db_management/gisaid2.db')
    c = conn.cursor()
    query = f"""
    SELECT *
    FROM {table}
    WHERE {column} = '{whichone}'
    """
    c.execute(query)
    results = c.fetchall()
    print(results)
    conn.close()

