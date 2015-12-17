__author__ = 'f90089'
import os

rootdir=os.getcwd()
directory = rootdir+ '/OutPut2/'
######################################### Step 2 #########################################
############ Adding Exact gene location and patterns group to the second file ############



extra_infile = open("newtable_final.txt", 'r')
e_indata = extra_infile.readlines()
extra_indata = sorted(set(e_indata))

extra_outfile = open(directory +extra_infile.name.replace('.txt','')+'_Analyzed_onlyPatternRule2.txt', 'w')
extra_outfile.write('Ensemble Gene ID\tEnsemble Transcript ID\tPattern Start (bp)\tPattern End (bp)\tPattern\tGroup\tCluster\n')


# gene_pattern = dict() # Dictionary to link each gene with all the patterns associated with it (by pattern id)
# pattern_line = dict() # Dictionary to link the pattern id with the pattern details

counter =0
larger=0
for line in extra_indata:
    pattern_info = line.split()


    # # comment below code to make speerte in external file
    # # we won't write the group and cluster column in the generated file (extra_outfile) here.
    # #
    # # Find the pattern group and cluster
    #
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
                        print('error:'+str(counter) +' Lenght:13 [index 5]= A and not in cluster 1 A. SubSeq:' +str(pattern_seq[4:9]) +' Seq:' + pattern_seq + ' match:'+ match +'\n' )
                        p_group = 0
                        p_cluster = 0
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
                        print('error:'+str(counter) +' Lenght:13 [index 5]= T and not in cluster 2 T. SubSeq:' +str(pattern_seq[3:10]) +' Seq:' + pattern_seq + ' match:'+ match +'\n' )
                        p_group = 0
                        p_cluster = 0
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
                    p_group = 0
                    p_cluster = 0
                    break

            #end of Extra

            p_group = 3
            p_cluster = 3
    extra_outfile.write(line.replace('\n','') + '\t'+ str(p_group) + '\t'+ str(p_cluster) + '\n')
    #
    # End of clustering and grouping
