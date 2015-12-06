__author__ = 'MacUser'

#!/usr/bin/python

import sys

file_names = []
file_names.append("chrx_Full_table_tab_gene.txt newtable_final.txt chrx\ ")


# For each chromosome

for c in range(0, len(file_names)):

    file_name = file_names[c].split(' ')
    name1 = file_name[0].split('.')
    name2 = file_name[1].split('.')
    directory = file_name[2]

    # get chromosome number from file name

    chr_file_name = str(name1[0])
    chr_num = chr_file_name[chr_file_name.find("r")+1:chr_file_name.find("_")]


    ######################################### Step 1 #########################################
    ###################### Adding Introns information to the first file ######################


    chro_outfile = open(name1[0]+'_Analyzed.txt', 'w')
    chro_outfile.write('Gene Name\tEnsemble Gene ID\tEnsemble Transcript ID\tGene Start (bp)\tGene End (bp)\tTranscript Count\tDescription\tExon Rank in Transcript\tExon Chr Start (bp)\t Exon Chr End (bp)\tStrand\tIntron Rank in Transcript\tIntron Chr Start (bp)\tIntron Chr End (bp)\n')

    chro_infile = open(file_name[0], 'r')
    chro_indata = chro_infile.readlines()

    intron_rank = 0
    gene_dict = dict()

    # Going through all the lines in the file and find the intron between every two consecutive exons in the same transcript

    for line in range(1,len(chro_indata)):
        if line != len(chro_indata)-1:

            data = chro_indata[line].split('\t')
            next_data = chro_indata[line+1].split('\t')
            strand = data[10]
            # save unique genes with their start, end ,strand, and gene name
            # no need for the multiple Exon for same gene in this file.
            # key is the gene ensemble ID
            if data[1].strip() not in gene_dict.keys():
                gene_dict[data[1].strip()] = str(str(data[3].strip()) + ' ' + str(data[4].strip()) + ' ' + str(data[10].strip()) + ' ' + str(data[0].strip()))


            if int(next_data[7].strip()) == 1:
                intron_rank = 0
            else:
                if int(data[7].strip()) == 1:
                    intron_rank = 1
                else:
                    intron_rank += 1

                if int(strand) == 1:
                    intron_start = int(data[9].strip()) + 1
                    intron_end = int(next_data[8].strip()) - 1
                elif int(strand) == -1:
                    intron_start = int(next_data[9].strip()) + 1
                    intron_end = int(data[8].strip()) - 1

            if intron_rank != 0:
                chro_outfile.write(str(chro_indata[line].strip()) + '\t' + str(intron_rank) + '\t' + str(intron_start) + '\t' + str(intron_end) + '\n')
            else:
                chro_outfile.write(str(chro_indata[line].strip()) + '\t0\t0\t0\n')

        else:
            chro_outfile.write(str(chro_indata[line].strip()) + '\t0\t0\t0\n')
            if next_data[1].strip() not in gene_dict.keys():
                gene_dict[next_data[1].strip()] = str(str(next_data[3].strip()) + ' ' + str(next_data[4].strip()) + ' ' + str(next_data[10].strip()) + ' ' + str(next_data[0].strip()))




    chro_infile.close()
    chro_outfile.close()

    ######################################### Step 2 #########################################
    ############ Adding Exact gene location and patterns group to the second file ############


    extra_outfile = open(name2[0]+'_Analyzed.txt', 'w')
    extra_outfile.write('Ensemble Gene ID\tPattern Start (bp)\tPattern End (bp)\tPattern\tGene Start (bp)\tGene End (bp)\tStrand\tExact Pattern Location\tGroup\tCluster\tChromosome number\tGene Name\n')

    extra_infile = open(file_name[1], 'r')
    e_indata = extra_infile.readlines()
    extra_indata = sorted(set(e_indata))

    gene_pattern = dict() # Dictionary to link each gene with all the patterns associated with it (by pattern id)
    pattern_line = dict() # Dictionary to link the pattern id with the pattern details

    counter =0
    for line in extra_indata:
        pattern_info = line.split(' ')


        # Find the pattern group and cluster

        pattern_seq = pattern_info[3].strip()

        if len(pattern_seq) > 21:
            p_group = 1
            p_cluster = int((len(pattern_seq) - 1)/4)
        elif len(pattern_seq) >= 21:
            p_group = 1
            p_cluster = 5
        elif len(pattern_seq) == 17:
            p_group = 2
            p_cluster = 4
        elif len(pattern_seq) == 13:
            if pattern_seq[4] == 'A':
                mismatch = 0
                match = "ATTTATTTATTTA"
                for i in range(0,13):
                    if pattern_seq[i] != match[i]:
                        mismatch += 1
                    if mismatch > 1:
                        #Extra checking other cluster types
                        if pattern_seq[4:10] != 'ATTTAT':
                            # erro no cluster for this belong to  1, 2 ,3
                            counter+=1
                            print('error:'+str(counter) +' Lenght:13 [index 5]= A and not in cluster 1 A. SubSeq:' +str(pattern_seq[4:9]) +'Seq:' + pattern_seq + ' match:'+ match +'\n' )
                            break
                        #end of Extra

                        p_group = 5
                        p_cluster = 1
                        break
                if mismatch <= 1:
                    p_group = 3
                    p_cluster = 3

            elif pattern_seq[4] == 'T':
                mismatch = 0
                match = "ATTTATTTATTTA"
                for i in range(0,13):
                    if pattern_seq[i] != match[i]:
                        mismatch += 1

                    if mismatch > 1:
                        #Extra checking other cluster types
                        if pattern_seq[3:10] != 'TTTATTT':

                            # erro no cluster for this belong to  1, 2 ,3
                            counter+=1
                            print('error:'+str(counter) +' Lenght:13 [index 5]= T and not in cluster 2 T. SubSeq:' +str(pattern_seq[3:10]) +'Seq:' + pattern_seq + ' match:'+ match +'\n' )
                            break
                        #end of Extra

                        p_group = 4
                        p_cluster = 2
                        break
                if mismatch <= 1:
                    p_group = 3
                    p_cluster = 3
            else:
                #Extra checking other cluster types
                mismatch=0
                for i in range(0,13):
                    if pattern_seq[i] != match[i]:
                        mismatch += 1
                    if mismatch <= 1:
                        p_group = 3
                        p_cluster = 3
                    else:
                        # erro no cluster for this belong to  1, 2 ,3
                        counter+=1
                        print('error:'+str(counter) +' Lenght= 13 and not A/T and not in cluster 3. seq:' + pattern_seq + ' match:'+ match +'\n' )
                        break

                #end of Extra

                p_group = 3
                p_cluster = 3


        # Calculate the exact location of the pattern in the gene

        gene_loc = gene_dict[pattern_info[0].strip()]
        gene_location = gene_loc.split(' ')
        g_start = int(gene_location[0])
        g_end = int(gene_location[1])
        strand = int(gene_location[2])
        gene_name = str(gene_location[3])


        p_start = int(pattern_info[1]) - 15
        exact_location = 0
        if strand == 1:
            exact_location = int(g_start) + int(p_start)
        elif strand == -1:
            exact_location = int(g_end) - int(p_start)


        # Write the analyzed information in the output file

        extra_outfile.write(str(pattern_info[0]).strip() + '\t' + str(p_start).strip() + '\t' + str(pattern_info[2]).strip() + '\t' + str(pattern_info[3]).strip() + '\t' + str(g_start).strip() + '\t' + str(g_end).strip() + '\t' + str(strand) + '\t' + str(exact_location).strip() + '\t' + str(p_group) + '\t' + str(p_cluster) + '\t'+ str(chr_num) + '\t' + gene_name +'\n')


#               # Update the dictionaries to be used in step 3
#               if pattern_info[0] in gene_pattern.keys():
#                       gene_pattern[pattern_info[0].strip()].append(exact_location)
#               else:
#                       gene_pattern[pattern_info[0].strip()] = []
#                       gene_pattern[pattern_info[0].strip()].append(exact_location)

        # include only the largest one

        if pattern_info[0] in gene_pattern.keys():
            old_p_location = gene_pattern[pattern_info[0].strip()][0]
            old_p = pattern_line[old_p_location].split('\t')
#                        print(pattern_info[0])
#                        print(old_p[2].strip())
#                        print(pattern_seq)
#                        exit()
            if len(old_p[2].strip()) < len(pattern_seq):
                 gene_pattern[pattern_info[0].strip()][0] = exact_location
     #                   gene_pattern[pattern_info[0].strip()].append(exact_location)
            elif len(old_p[2].strip()) == len(pattern_seq):
                 gene_pattern[pattern_info[0].strip()].append(exact_location)
        else:
               # print(pattern_info[0])
               # print(pattern_seq)
               # exit()
            gene_pattern[pattern_info[0].strip()] = []
            gene_pattern[pattern_info[0].strip()].append(exact_location)

        pattern_l = str(str(p_start).strip() + '\t' + str(pattern_info[2]).strip() + '\t' + str(pattern_info[3]).strip() + '\t' + str(exact_location).strip() + '\t' + str(p_group) + '\t' + str(p_cluster) + '\t'+ str(chr_num))
        pattern_line[exact_location] = pattern_l

    extra_infile.close()
    extra_outfile.close()
    #print(gene_pattern)


    ######################################### Step 3 #########################################
    ######################### Combining the first and the seond files ########################


    # Writing the file with the full final informations

    full_data = open(name1[0]+'_Final.txt', 'w')
    full_data.write('Gene name\tV1\tEnsemble TranscriptID\tGene Start (bp)\tGene End (bp)\tTranscript Count\tDescription\tExon Rank in Transcript\tExon Chr Start (bp)\t Exon Chr End (bp)\tStrand\tIntron Rank in Transcript\tIntron Chr Start (bp)\tIntron Chr End (bp)\tPattern Start\tPattern End\tPattern\tExact Pattern Location\tGroup\tCluster\tChromosome number\tFound in\n')

    # Writing the file with only exon patterns final informations

    exons_data = open(name1[0]+'_Final_Exons.txt', 'w')
    exons_data.write('Gene name\tV1\tEnsemble TranscriptID\tGene Start (bp)\tGene End (bp)\tTranscript Count\tDescription\tExon Rank in Transcript\tExon Chr Start (bp)\t Exon Chr End (bp)\tStrand\tIntron Rank in Transcript\tIntron Chr Start (bp)\tIntron Chr End (bp)\tPattern Start\tPattern End\tPattern\tExact Pattern Location\tGroup\tCluster\tChromosome number\tFound in\n')

    # Writing the file with only intron patterns final informations

    introns_data = open(name1[0]+'_Final_Introns.txt', 'w')
    introns_data.write('Gene name\tV1\tEnsemble TranscriptID\tGene Start (bp)\tGene End (bp)\tTranscript Count\tDescription\tExon Rank in Transcript\tExon Chr Start (bp)\t Exon Chr End (bp)\tStrand\tIntron Rank in Transcript\tIntron Chr Start (bp)\tIntron Chr End (bp)\tPattern Start\tPattern End\tPattern\tExact Pattern Location\tGroup\tCluster\tChromosome number\tFound in\n')

    # Integrate the patterns file with the chromosome file to attach each pattern to its intron or exon

    chro_infile = open(name1[0]+'_Analyzed.txt', 'r')
    skip = chro_infile.readline()
    chro_indata = chro_infile.readlines()

    for lines in chro_indata:
        intron = 0
        exon = 0
        patterns = []
        line = lines.split('\t')

        if line[1].strip() in gene_pattern.keys():

            pattern_location = gene_pattern[line[1]]

            # For each line in the full data file, finding all the patterns associated with gene location using the pre-defined dictionary

            for pl in pattern_location:
                if int(pl) >= int(line[8]) and int(pl) <= int(line[9]):
                    p_info = pattern_line[pl]
                    details = lines.split('\t')
                    patterns.append(p_info  + '\t' + details[2] + '_Exon' + details[7])


                elif int(pl) >= int(line[12]) and int(pl) <= int(line[13]):
                    p_info = pattern_line[pl]
                    details = lines.split('\t')
                    patterns.append(p_info + '\t' + details[2] + '_Intron' + details[11])



            # If there is only one pattern for that line, wite it in the output files


            if len(patterns) == 1:
                info = patterns[0]
                full_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')
                in_or_ex = info.split('\t')

                if 'Exon' in in_or_ex[7]:
                    exons_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')
                elif 'Intron' in in_or_ex[7]:
                    introns_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')

            # If there is more than one pattern, prioritize them based on the cluster and write the pattern with the best cluster in the output files

            elif len(patterns) > 1:
                priority = 0
                index = 0
                for p in range(0, len(patterns)):
                    patt = patterns[p].split('\t')
                    if int(patt[5]) > priority:
                        priority = int(patt[5])
                        index = p

                info = patterns[index]
                full_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')

                in_or_ex = info.split('\t')
                if 'Exon' in in_or_ex[7]:
                    exons_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')
                elif 'Intron' in in_or_ex[7]:
                    introns_data.write(str(lines.strip()) + '\t' + str(info.strip()) + '\n')

    chro_infile.close()
    full_data.close()
    exons_data.close()
    introns_data.close()

    print ('The analysis of chromosome ' + str(chr_num) + ' is done')




