#!/usr/bin/env python
import sys
import datetime
import os
import csv
import argparse
import subprocess
from argparse import ArgumentParser

# Running CRISPRone
def runCRISPRone(crisprone, output_base, output_prefix, input_fna):
    # usage: 
    # CRISPRone argument
    crispr = [ 'python3', crisprone, output_base, output_prefix, input_fna ]
    output= subprocess.call(crispr)
    return output

# Extracting spacer sequence
def spacerSequence(infile, outfile):
    """
    The function spacerSequence converts a 
    .crt file from the crisprone output, into 
    a FASTA formatted file to perform BLAST.
    """
    # Creating a list 
    lis = [line.split() for line in infile]
    result = []  # Saving results of list in variable called result
    for lister in lis:
        if len(lister) ==7: # Selecting the length of the list
            result.append(lister[2]+":::"+lister[0]+":::"+lister[4]+":::"+lister[5]) 
        
        else:
            try:
               # print(lister[0])
                if 'SEQ'in lister[0]:
                    result.append(lister[1]+"###")
            except IndexError:
                pass
    with open(outfile, 'w') as outfile:
        count = 1

        print("Converting FASTA....")
        temp = ''
        for strLine in result:
         #   print(strLine)
           # temp = ''
            if strLine != "Range:" and strLine != "Average" and strLine != "1" and strLine != "SPACER" and "-" not in strLine and "CRISPR" not in strLine and "Bases" not in strLine:  # Removing unwanted strings

                strLine = strLine.rstrip("\n")
               # print(strLine)
            if "###" in strLine:
                temp = strLine.replace("###", "").strip()  #temp gives us the ID 
               #i print(temp)
            elif ":::" in strLine:
             #   print(strLine)
                strLine=strLine.replace(",","")
                #If negative positions in the file, not to consider and add start position[1] with 0.
                if int(strLine.split(":::")[1])+int(strLine.split(":::")[2]) < 0 or int(strLine.split(":::")[1])+int(strLine.split(":::")[2])+int(strLine.split(":::")[3]) < 0:
                   # print(strLine)
                    outfile.write(">" + temp+ "_" +str(1+int(strLine.split(":::")[2].replace(",",",")))+"_"+str(int(strLine.split(":::")[2].replace(",",""))+int(strLine.split(":::")[3]))+"\n")
                   # print(strLine)
                    outfile.write(strLine.split(":::")[0] + "\n")

               # strLine=strLine.replace(",","")
                #    print(strLine)
                else:
                    # Writing output to a fasta format
                    outfile.write(">" + temp+ "_" +str(0+int(strLine.split(":::")[1])+int(strLine.split(":::")[2].replace(",",",")))+"_"+str(-1+int(strLine.split(":::")[1])+int(strLine.split(":::")[2].replace(",",""))+int(strLine.split(":::")[3]))+"\n") 

                    outfile.write(strLine.split(":::")[0] + "\n")

    print("Done")


# Running Blast on the genome to create a database

def runblast(makeblastdb, fasta, dbtype, blastoutdir, flag):
    if flag == 0:
        blast_cmd = [makeblastdb, '-in', fasta, '-dbtype', 'nucl', '-out', blastoutdir]
        output = subprocess.call(blast_cmd)
    elif flag == 1:
        blast_cmd2 = [makeblastdb, '-in', fasta, '-dbtype','prot','-out', blastoutdir]
        output = subprocess.call(blast_cmd2)
    return output
    


# Running Blastn for the concerted crt file into fasta, against the genome database
def runblastn(blastn,  fasta, database_name, formating, dirr, flag):
    if flag == 0:
        cmd = [blastn, '-query', fasta, '-db', database_name, '-outfmt', formating, '-out', dirr]
        output = subprocess.call(cmd)
    elif flag == 1:
        cmd1 = [ blastn, '-query', fasta, '-db', database_name, '-outfmt', formating, '-out', dirr ]
        output = subprocess.call(cmd1)
    return output
    # print(cmd)

# Obtaining Self Targets

def selfTarget(inputfile1,inputfile2, outputfile):
    range_start = []
    range_end = []
    with open(inputfile1,'r') as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            try:
                if 'Range' in row[0]:
                    range_start.append(row[0].split(":")[1].split(row[0].split(":")[1][0])[1]) 
                    range_end.append(row[0].split(":")[1].split(row[0].split(":")[1][0])[3])
            except:
                pass
    l = len(range_start)
    file2write = open(outputfile, 'w')

    with open(inputfile2,'r') as f:
        range_list = []
        reader = csv.reader(f)
        for row in reader:
            temp_start = int(row[0].split("\t")[8])
            count = 0
            temp_end = int(row[0].split("\t")[9])
            for i in range(0, l):
                if int(range_start[i]) < temp_start < int(range_end[i]) or int(range_start[i]) < temp_end < int(range_end[i]):
                    count  = count + 1
            if count == 0:
                file2write.write(row[0])
                file2write.write("\n")
    file2write.close()

# Acranker

def ac_ranker(Acranker,input_filename, output_filename):
    ranker = ['python3', Acranker, input_filename, output_filename]
    output = subprocess.call(ranker)
    return output

# Phage regions

def phage_region(csv_file, tsv_file, result_file):
     ranges_0 = []
     ranges_1 = []
     with open(tsv_file, newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
         count = 0
         for row in spamreader:
             count = count + 1
             temp_value = row[0]
             if count != 1:
                 ranges_0.append(int(temp_value.split("\t")[3]))
                 ranges_1.append(int(temp_value.split("\t")[4]))

     file2write = open(result_file, "w")
     with open(csv_file, newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
         count = 0
         c = 0
         for row in spamreader:
             count = count + 1
             if count != 1:
                 temp_value = row[0]
                # print(row)
                 range_0 = int(temp_value.split("_")[1])
                 range_1 = int(temp_value.split("_")[2])
                 for i in range(0, len(ranges_0) - 1):
                     if range_0 >= ranges_0[i] and range_0 <= ranges_1[i] and range_1 >= ranges_0[i] and range_1 <= ranges_1[i]:
                         c = c + 1
                         #print(row)
                         res = str(row)[1:-1]
                         res = res[1:len(res)-1]
                         #print(res)
                         file2write.write(res)
                         file2write.write("\n")
         print(c)
         file2write.close()

# Protein fasta file 
def prot_results(acrFile, faa_file, prot_fasta):
    #Files_Path = "/home/sithata/phage/gen162/GCA_000006765.1_ASM676v1.txt" #.txt file ( Acranker and virsorter2)
    #FOLDER = "/home/team/ont/wastewater/siri_wd/scripts/FragGeneScan/file_165/FAA/GCA_000006765.1_ASM676v1-FSG.faa" # .faa file 
    final_result = []
    result = {}
    with open(acrFile) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        for row in csv_reader:
            temp_name = str(row)[1:len(str(row))-1].replace("'","").split(",")[0]
            with open(faa_file,encoding= "ISO-8859-1") as csv_file_1:
                csv_reader_1 = csv_file_1.readlines()
                count = 0
                for row1 in csv_reader_1:
                    if temp_name in row1 and count != 1:
                        count = 1
                    elif count == 1:
                        count = 0
                        result[temp_name] = row1
                        break

    f = open(prot_fasta,"w") # prot fasta file req to run blastp against acr db
    for key in result:
        f.write(">"+key)
        f.write("\n")
        f.write(result[key])
    f.close()

# Highest score from Blast result 
def high_score(params):
    result_csv = params[0]
    final_result_txt = params[1]
    values_dict = {}
    results_dict = {}
    with open(result_csv) as csvfile:
        reader = csvfile.readlines()
        for row in reader:
            #print(row.replace('\t',"#"))
            if len(row.split('\t'))>2:#row.split('\t')[0] not in values_dict and len(row.split('\t'))>=2:
               values_dict[row]= float(row.split('\t')[11])
               results_dict[row] = row
    #print(values_dict)
    dup_find =  []
    sorted_x = sorted(values_dict.values())
    #print(sorted_x)
    f = open(final_result_txt,"w")
    for key in values_dict:
        if values_dict[key] == sorted_x[len(sorted_x)-1] and results_dict[key] not in dup_find:
            dup_find.append(results_dict[key])
    final_res = []
    for item in dup_find:
        temp_ = item.split('\t')[0]
        if temp_ not in final_res:
            final_res.append(temp_)
            f.write("\n")
            f.write(item)
    f.close()

##Summary Table
def summary(params):
    hs_txt = params[0]
    final_viral = params[1]
    phage_region = params[2]
    self_target = params[3]
    summary_report = params[4]
    gene_name = params[5]

    #high_score
    f = open(hs_txt,"r")
    f1 = f.readlines()
    hs = {}
    for x in f1:
        #print(x)
        if len(x.strip())>0:
            ID = x.split('\t')[0]
            ACR_ID = x.split('\t')[1]
            Score_ID = x.split('\t')[2]
            Align_ID = x.split('\t')[3]
            e_value = x.split('\t')[10]
            bit_score = x.split('\t')[11]
            hs[ID] = str(ACR_ID) + "&&&" + str(Score_ID) + "&&&" + str(Align_ID)+ "&&&" + str(e_value) + "&&&" + str(bit_score)
    f.close()
    #final_viral 
    with open(final_viral, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter = ',',quotechar = '|')
        fv = {}
        count = 0
        for rows in reader:
            rowser = rows[0]            
            #print(rows)
            row = rowser.split('\t')
           # print(count)
            count =count + 1
            if count != 1:
                ID = row[0]
                bp_start = row[3]
                bp_end = row[4]
                fv[str(ID)+"&&&"+str(count)] = str(bp_start) + "&&&" + str(bp_end)
            
    #phage_region
    fp = open(phage_region,"r")
    f1 = fp.readlines()
    pr = {}
    for x in f1:
       if len(x.strip()) > 0: 
          #print(x)  
          ID = x.split(',')[0]
          #ACR_ID = x.split('\t')[1]
          Score_ID = x.split(',')[2]
          pr[ID] = str(Score_ID)
    fp.close()
    #Self_target
    f = open(self_target,"r")
    f1 = f.readlines()
    st = {}
    for x in f1:
      if len(x.strip())>0:  
         ID = x.split('\t')[0]
         #ACR_ID = x.split('\t')[1]
         #Score_ID = x.split('\t')[2]
         st[ID] = "True"#str(ACR_ID) + "&&&" + str(Score_ID)
    f.close()
    #write a summary report
    with open (summary_report, 'w', newline = "") as csvfile:
        new_header = ['Genome', 'Self_target', 'ACranker score','ACranker ID', 'Phage_start_position', 'Phage_end_position', 'Anti CRISPR ID', 'ID score', 'Alignment_length', 'E_Value', 'Bitscore']
        csvwriter = csv.DictWriter(csvfile,fieldnames=new_header)
        csvwriter.writeheader()
        bps = ''
        bpe = ''
 #       print(hs)
  #      print(fv)
        for item in hs:
            range_s =int(item.split('_')[1])
            range_e = int(item.split('_')[2])
           # print(str(range_s))
           # print(str(range_e))
            bool_st = 'Yes'
           # bps = ''
            #bpe = ''
            if len(st) > 0:
                bool_st = 'Yes'
            else:
                bool_st = 'No'
            score_AC = pr[item]
            for val in fv:
                temp_val = fv[val]
                temp_s = int(temp_val.split('&&&')[0])
                temp_e =int(temp_val.split('&&&')[1])
                #print(str(temp_s))
                #print(str(temp_e))
                #print(str(range_s)) 
                #print(str(range_e))
                if temp_s < range_s and temp_e > range_s and temp_s < range_e and temp_e > range_e:
                       bps = temp_s
                       bpe = temp_e
                       break   
            #bps = fv[item].split('&&&')[0]
            #bpe = fv[item].split('&&&')[1]
            acID = hs[item].split('&&&')[0]
            hs_ = hs[item].split('&&&')[1]
            a_l = hs[item].split('&&&')[2]
            e_val = hs[item].split('&&&')[3]
            b_s = hs[item].split('&&&')[4]
            csvwriter.writerow({"Genome":gene_name, 'Self_target':bool_st, 'ACranker score':score_AC,'ACranker ID':item, 'Phage_start_position':bps, 'Phage_end_position':bpe, 'Anti CRISPR ID':acID, 'ID score':hs_, 'Alignment_length':a_l, 'E_Value':e_val, 'Bitscore': b_s})
    return "Completed"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute a genome')
    #parser.add_argument('-i', '--indir', help = 'Path to input fasta file')
    parser.add_argument('-f', '--input_fna', help= 'name and path to fasta file', required=True)
    #parser.add_argument('-p', '--output_prefix', help = 'CRISPRone output prefix') 
    parser.add_argument('-b', '--output_base', help= 'directory of the output file', required=True) 
    parser.add_argument('--crisprone', help = 'CRISPRone path location', default = '/home/sithata/cripsrone/CRISPRone/crisprone-local.py')
    #parser.add_argument('-input', '--infile', help = 'input crt file', default = 'file.crt') 
    #parser.add_argument('-output', '--outfile', help = 'path to outfile', default = 'test_outfile.fasta')
    #parser.add_argument('-f', '--fasta', help = 'query sequence')
    #parser.add_argument('-in', '--indir', help = 'path to input directory')
    #parser.add_argument('-out', '--blastoutdir', help = 'Directory to Output file')
    #parser.add_argument('-dbtype', '--database_type', help= 'Type of database')
    parser.add_argument('--makeblastdb', help = 'BLAST path location', default = '/home/sithata/bin/ncbi-blast-2.10.1+/bin/makeblastdb')
    #parser.add_argument('-query', '--fasta', help = 'query sequence')
    #parser.add_argument('-outfmt', '--formating', help = 'Formating options')
    #parser.add_argument('-o', '--dirr', help = 'Directory to Output file')
    #parser.add_argument('-db', '--database_name', help= 'Database name')
    parser.add_argument('--blastn', help = 'BLAST path location', default = '/home/sithata/bin/ncbi-blast-2.10.1+/bin/blastn')
    parser.add_argument('--blastp', help = 'BLAST path location', default = '/home/sithata/bin/ncbi-blast-2.10.1+/bin/blastp')
    #parser.add_argument('-in1', '--inputfile1', help = 'input crt file')
    #parser.add_argument('-in2', '--inputfile2', help = 'BLAST output')
    #parser.add_argument('-selftar', '--outputfile', help = 'path to outputfile', default = 'test_outputfile.txt')
    #parser.add_argument('-in', '--input_filename', help = 'input file Acranker')
    #parser.add_argument('-O', '--output_filename', help = ' output file Acranker')
    parser.add_argument('--Acranker', help = 'AcRanker path location', default= '/home/sithata/AcRanker/acranker.py')
    #parser.add_argument('-csv', '--csv_file', help = 'input csv file from AcRanker')
    parser.add_argument('-tsv', '--tsv_file', help = 'input tsv file from Virsorter2')
    #parser.add_argument('-res', '--result_file', help = 'path and name of outfile', default = 'test_outfile.txt')
    args = parser.parse_args()

###############################################################################################################################################################

## CRISPRone Function

    Genome_input_file = args.input_fna
    Genome_output_file = "AcrPrediction_"+ str(datetime.datetime.now()).replace(" ","_").replace(":","_")+ "_file"
    Genome_output_file_path = args.output_base
    Crisprone = args.crisprone

## BLAST - makeblastdb
    
    #Genome_input_file = args.fasta
    #Genome_input_file_path = args.indir
    Genome_Blast_output_file = Genome_output_file_path +"/" + Genome_output_file + "/" + "genome_database" #args.output_base
    Genome_dbtype = "" #args.database_type
    Makeblastdb = args.makeblastdb

## Extracting Spacers

    Spacer_infile = Genome_output_file_path + "/" + Genome_output_file +"/"+ Genome_output_file +".crt" 
    Spacer_outfile = Genome_output_file_path + "/" + Genome_output_file + "/" +"sequence.fasta"    
    #print(Spacer_infile)
    #print(Spacer_outfile)


## BLAST - blastn

   #Query_sequece = args.outfile
   #Database_name = args.blastoutdir   # Genome_Blast_output_file
    
    Output_file = Genome_output_file_path +"/" + Genome_output_file + "/" + "Blast_result.m8"  #args.dirr
    Formatting = "6" #args.formating
    Blastn = args.blastn

## Obtaining Self-targets
   # inputfile1 = args.outdir + args.prefix + '/' + args.prefix + '.crt'
    #inputfile2 = args.dirr
    outputfile = Genome_output_file_path +"/" + Genome_output_file + "/" + "selfTarget.txt"#args.outputfile

#Acranker function
    Acranker_input = Genome_output_file_path + "/" + Genome_output_file +"/"+ Genome_output_file + '.faa'
    Acranker_output= Genome_output_file_path + "/" + Genome_output_file +"/"+ "AcRanker"
    AcRanker = args.Acranker

# Phage region
    phage_csv = Acranker_output + '.csv'
    #os.chdir(r'/home/sithata/virsorter2/results/GCA_014672795.1_ASM1467279v1.out')
    phage_tsv = args.tsv_file #os.getcwd()+'/final-viral-boundary.tsv' #args.tsv_file
    #os.chdir(Genome_output_file_path + "/" + Genome_output_file)
    result = Genome_output_file_path + "/" + Genome_output_file + "/" + "phage_region.txt" 
    #result = os.getcwd()+ '/final_result.txt'

#Acr Database - BLAST
    Acr_infile = '/home/sithata/Acr_Database/proteinSequence.faa'
    Acr_outfile = Genome_output_file_path +"/" + Genome_output_file + "/"+ "AcrDB"

#Protein fasta file
    #acrFile= result
    #faa_file= Acranker_input
    prot_fasta= Genome_output_file_path +"/" + Genome_output_file + "/" + "Protein.fasta"

#Blastp against acrdb
    Blastp = args.blastp
    Formating = "6"
    Outdir = Genome_output_file_path + "/" + Genome_output_file + "/" + "Protein_seq.m8"
    #DBname = ACRBD"/home/team/ont/wastewater/siri_wd/scripts/python/trial/AcrDB"
    #input_fasta = prot_fasta


### Calling the functions

## CRISPRone
    #IF CRISPR NOT FOUND,DONT RUN
    runCRISPRone(Crisprone, Genome_output_file_path, Genome_output_file, Genome_input_file)

## BLAST - makedb  runblast(makeblastdb,indir, fasta, database_type, blastoutdir)
    runblast(Makeblastdb,  Genome_input_file, Genome_dbtype, Genome_Blast_output_file,0)
    runblast(Makeblastdb,  Acr_infile, Genome_dbtype, Acr_outfile,1)
## Extracting Spacers
# open infile
    with open(Spacer_infile, 'r') as infile:
          # open outfile
#        with open(Spacer_outfile, 'w') as outfile:
              # call get spacer function
            spacerSequence(infile,Spacer_outfile)

## BLAST - blastn  runblastn(blastn,  fasta, database_name, formating, dirr)
    runblastn(Blastn, Spacer_outfile, Genome_Blast_output_file, Formatting, Output_file, 0)
    #runblastn(Blastp, Outdir, Acr_outfile, Formatting, prot_fasta, 1)  

## Selt-Target
# Open inputfile1
    #with open(inputfile1, 'r') as inputfile1:
        # Open inputfile2
       # with open(inputfile2, 'r') as inputfile2:
            # Call the function
    selfTarget(Spacer_infile, Output_file, outputfile)

#Acranker
    ac_ranker(AcRanker, Acranker_input, Acranker_output)

## Phage region
    #with open(phage_csv, 'r') as csv_file:
       # with open(args.tsv_file, 'r') as tsv_file:
    phage_region(phage_csv, phage_tsv, result)


#Acr Databas
    prot_results(result, Acranker_input, prot_fasta)

# Blastp
    runblastn(Blastp, prot_fasta, Acr_outfile, Formatting, Outdir, 1)

# high Score

    high_score([Outdir,Genome_output_file_path + "/" + Genome_output_file + "/"+ "Highscore_result.m8"])
    
#Summary Table
    
    #summary(highscore, final_viral, phage_Region, selftarget, outputfile, -f)
    name_gene = Genome_input_file.split("/")[len(Genome_input_file.split("/"))-1]
    summary([Genome_output_file_path + "/" + Genome_output_file + "/"+ "Highscore_result.m8", args.tsv_file, result, outputfile, Genome_output_file_path + "/" + Genome_output_file + "/"+ "Summary_report.csv", name_gene])
   
