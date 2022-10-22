import os
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord

header = "Strain, Synteny-Group, Genome size, Chromosome size, Plasmid size, 16S rRNA length, ITS-1 Sequence,\
    Length of ITS-1 Seq., Total No. of features, No. CDS, No. rRNA, No. tRNA, No. misc_RNA, GC %"
    
#"I:/G1A_gbk_files/"
#"Desktop/Masterarbeit/Programmieren/G1A_gbk_files/"
directory = r"I:/Dateien/all_g1A_gbk_new/all_without_vers_acc/" #give pathway "r"  
all_g1a = []
ITS_seq_records = []

for (dirpath, dirnames, filenames) in os.walk(directory):
    for filename in filenames:
            
        information = [] #create empty list
            
        if filename.endswith(".gbk"):

            for record in SeqIO.parse(os.sep.join([dirpath, filename]), "genbank"):
                print("Es wird folgender Strain verwendet: ", record.description)
                
                # # splits string by " " and takes last listentry = the strain number
                # strainnumber = int(record.description.split()[-1])
                
                information.append(record.description)
            
#                 strain2group = {118: 'A', 143: 'A', 158: 'A', 197: 'A', 198: 'A', 222: 'A', 260: 'A',
#                                   2: 'B',  33: 'B',  38: 'B', 342: 'B', 379: 'B',
#                                 151: 'C', 193: 'C', 217: 'C', 428: 'C', 405: 'C', 425: 'C', 505: 'C',
#                                 126: 'D', 324: 'D',
#                                 135: 'E', 221: 'E', 242: 'E',
#                                 509: 'F', 511: 'F', 525: 'F', 625: 'F',
#                                 567: 'EF'}
    
# #                strain2samplingloc = {}
# #Das gleiche möchte ich auch mit den Metadaten machen

                
#                 if strainnumber in strain2group:
#                     print("Strain gehört der Syntheniegruppe", strain2group[strainnumber], "an")
#                     information.append(strain2group[strainnumber])
                
#                 else:
#                     print("Dieser Strain gehört bisher keiner Syntheniegruppe an")
#                     information.append("No group assigned yet")   
            
            
            
                # calculates the total length of the sequences of the GenBank file
                temp_genome = sum(len(genome) for genome in SeqIO.parse(os.sep.join([dirpath, filename]), "genbank"))
                information.append(temp_genome) #append genome size
                
                #das funkt noch nicht wirklich
                temp_chromosome = max(len(genome) for genome in SeqIO.parse(os.sep.join([dirpath, filename]), "genbank"))
                information.append(temp_chromosome)    #max fragment , werte in array adden un dann max funktion 
                
                information.append(temp_genome - temp_chromosome) # calculate and append plasmid legth (should be 1-200 kb)
                
                total_features = 0
                cds = 0
                rrna = 0
                tRNA = 0
                misc_RNA = 0
                
                for gene in record.features: #calculate no. of cds/rrna/... features
                    
                    if gene.type == "gene": 
                        total_features += 1
                        
                    if gene.type == "CDS": 
                        cds += 1
                        
                    if gene.type == "rRNA": 
                        rrna += 1
                        
                    if gene.type == "tRNA": 
                        tRNA += 1
            
                    if gene.type == "misc_RNA": 
                        misc_RNA += 1      

                    if gene.type == "rRNA": #calculate 16S rRNA-length
                        if "product" in gene.qualifiers:
                            if "16S" in gene.qualifiers["product"][0]: 
                                
                                start_16s = gene.location.nofuzzy_start
                                end_16s = gene.location.nofuzzy_end
                                information.append(abs(start_16s - end_16s)) 

                            if "23S" in gene.qualifiers["product"][0]: 
                        
                                start_23s = gene.location.nofuzzy_start
                                end_23s = gene.location.nofuzzy_end 
                                
                print("Die 16S rRNA Sequenz startet ab %s und endet bei %s" % (start_16s + 1, end_16s))
                print("Die 23S rRNA Sequenz startet ab %s und endet bei %s" % (start_23s + 1, end_23s))
                
                if start_16s < start_23s: #für Berechnung der its-1 wenn 16S vor 23S
                    ITS_seq = record.seq[end_16s : start_23s]
                    print("Das heißt, dass die ITS-1 Region bei Position", end_16s + 1, "beginnt und bei", start_23s, "endet")
                    print("Die Sequenz lautet: ", repr(ITS_seq), "und hat eine länge von: ", len(ITS_seq), "\n")
                    information.append(ITS_seq)
                    information.append(len(ITS_seq))
                    
                if  start_16s > start_23s:   #für Berechnung der its-1 wenn 23S vor 16S 
                    ITS_seq = record.seq[end_23s : start_16s]
                    ITS_seq = ITS_seq.reverse_complement()
                    print("Das heißt, dass die ITS-1 Region bei Position", end_23s + 1, "beginnt und bei", start_16s, "endet")
                    print("Die Sequenz lautet: ", repr(ITS_seq), "und hat eine länge von: ", len(ITS_seq), "\n")
                    information.append(ITS_seq) #reversed complement, bc on seq. (-) strain 23S <- 16S
                    information.append(len(record.seq[end_23s : start_16s]))
                
                
                #daraus wird jetzt ein seqRecord gemacht, damit dass später in eine Multifasta datei umgewandelt werden kann
                # if strainnumber in strain2group:
                #     ITS_seq_records.append(SeqRecord(ITS_seq, id = str(strainnumber) + str(strain2group[strainnumber]), description = record.description))
                    
                # else:
                #     ITS_seq_records.append(SeqRecord(ITS_seq, id = str(strainnumber) + "NoGroup", description = record.description))  
                    
                #sig_peptide
                information.append(total_features) # calculates total number of features                              
                information.append(cds)
                information.append(rrna)
                information.append(tRNA)
                information.append(misc_RNA)
                information.append(round(GC(record.seq),)) #problem mit 55.5 in csv -> Excel 555, daher keine decimal.
                
                "Create separate fasta files for ITS-Sequence"
                # SeqIO.write(ITS_seq_records, "I:/G1A_gbk_files/ITS_Sequences_G1A_" + str(strainnumber) + ".fasta", "fasta",)

            all_g1a.append(information)  
            
    # print(all_g1a)
                
    
    
"Safe ITS-Sequences as multifasta-file"
SeqIO.write(ITS_seq_records, "I:/Dateien/all_g1A_gbk_new/all_without_vers_acc/ITS_Sequences_G1A_Synthgroup.fasta", "fasta",)


"Safe as csv-file"

output_handle = open("I:/Dateien/all_g1A_gbk_new/all_without_vers_acc/G1A.csv", "w")

output_handle.write(header + "\n") 

for line in all_g1a:
    output_line = "" #convert int into strings
    
    for n in line:
        output_line += str(n) + ","
    
    output_line = output_line.rsplit(",", 1)[0]
    output_line += "\n"    
    output_handle.write(output_line)

output_handle.close()

print("Extraktion abgeschlossen") 

   
"""
Noch Probleme:
    - wie finde ich die Plasmidgröße heraus ? (alle daten scheinen plasmidlos zu sein, außerdem fehler im code?)
    - nicht zirkularisierte gbk files machen komische sachen -> vermutlich ähnliches problem wie mit plasmiden
    - dict ldrnen, hash funktion
    - Metadaten integrieren!
    
"""