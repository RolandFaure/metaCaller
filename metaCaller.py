"""
SNPsnoop
Authors: Roland Faure, based on a previous program (strainminer) by Minh Tam Truong Khac
"""

__version__ = '0.3.2'

import pandas as pd 
import numpy as np
import math
import sys
import os
import pysam as ps

from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
from sklearn.metrics import pairwise_distances

from copy import deepcopy

import scipy.stats as stats


import time
from argparse import ArgumentParser

import threading
from concurrent.futures import ProcessPoolExecutor

fasta_lock = threading.Lock()
bam_lock = threading.Lock()
 
def log_comb(n, k):
    """Compute the logarithm of the combination n choose k. This is useful to avoid numerical errors when n is large"""
    if n < 1000:
        return np.log(float(math.comb(n, k)))
    else:
        return np.sum(np.log(np.arange(n, n - k, -1))) - np.sum(np.log(np.arange(1, k + 1)))

def statistical_test(a,n,b,m,error_rate):
    """Perform a statistical test to determine if a rectangle of 0s in a matrix is significant"""
    result_log = log_comb(n,a) + log_comb(m,b) + (a * b) * np.log(error_rate)
    if result_log > 0: #to avoid overflows
        return 1
    return np.exp(result_log)

def get_data(bamfile, originalAssembly, contig_name,start_pos,stop_pos, no_snp_threshold = 0.95):
    """
    Extracts a binary matrix of reads and columns from a BAM file, identifying suspicious positions in a specified contig region.
    Args:
        bamfile (str): Path to the sorted and indexed BAM file.
        originalAssembly (str): Path to the reference genome assembly in FASTA format.
        contig_name (str): Name of the contig to analyze.
        start_pos (int): Start position of the region to analyze.
        stop_pos (int): Stop position of the region to analyze.
        no_snp_threshold (float, optional): Threshold for considering a position as suspicious. Default is 0.95.
    Returns:
        tuple: A tuple containing:
            - matrix (dict): A dictionary where keys are positions and values are dictionaries mapping read names to binary values (1 for reference base, 0 for main alternative base, np.nan for others).
            - list_of_suspicious_positions (list): A list of positions identified as suspicious.
    """
    
    bamfile_ps = ps.AlignmentFile(bamfile,'rb')

    #first classify the reads in several groups, according to what they span in the window
    breakpoints_left = [stop_pos] #bp with reads left
    breakpoints_right = [start_pos] #bp with reads right
    past_aligned = -1
    for pileupcolumn in bamfile_ps.pileup(contig=contig_name, start=start_pos, stop=stop_pos, truncate=True, min_base_quality=0, multiple_iterators=True):
        if past_aligned == -1 :
            past_aligned = pileupcolumn.nsegments
        elif pileupcolumn.nsegments - past_aligned > max(5, 0.1*pileupcolumn.nsegments):
            breakpoints_right.append(pileupcolumn.reference_pos)
        elif past_aligned - pileupcolumn.nsegments > max(5, 0.1*pileupcolumn.nsegments):
            breakpoints_left.append(pileupcolumn.reference_pos)
        past_aligned = pileupcolumn.nsegments
    

    allreads_anywhere_on_the_window = set()
    reads_in_breakpoint_right = [set() for i in breakpoints_right]
    for b, bp_right in enumerate(breakpoints_right) :
        itercol_first_column = bamfile_ps.pileup(contig = contig_name,start = min(bp_right+10, stop_pos-1),stop = min(bp_right+11, stop_pos),truncate=True,min_base_quality=0, multiple_iterators=True)
        for pileupcolumn in itercol_first_column:
            allreads = pileupcolumn.get_query_names()
            reads_in_breakpoint_right[b] = set(allreads)
            allreads_anywhere_on_the_window = allreads_anywhere_on_the_window.union(set(allreads))
            break
        if bp_right != start_pos:
            itercol_first_column = bamfile_ps.pileup(contig = contig_name,start = max(bp_right-10, start_pos),stop = max(bp_right-9, start_pos+1),truncate=True,min_base_quality=0, multiple_iterators=True)
            for pileupcolumn in itercol_first_column:
                allreads = pileupcolumn.get_query_names()
                reads_in_breakpoint_right[b] = reads_in_breakpoint_right[b].difference(set(allreads))
                break

    reads_in_breakpoint_left = [set() for i in breakpoints_left]
    for b, bp_left in enumerate(breakpoints_left) :
        itercol_first_column = bamfile_ps.pileup(contig = contig_name,start = max(bp_left-11, start_pos),stop = max(bp_left-10, start_pos+1),truncate=True,min_base_quality=0, multiple_iterators=True)
        for pileupcolumn in itercol_first_column:
            allreads = pileupcolumn.get_query_names()
            reads_in_breakpoint_left[b] = set(allreads)
            allreads_anywhere_on_the_window = allreads_anywhere_on_the_window.union(set(allreads))
            break
        if bp_left != stop_pos:
            itercol_first_column = bamfile_ps.pileup(contig = contig_name,start = min(bp_left+9, stop_pos-1),stop = min(bp_left+10, stop_pos),truncate=True,min_base_quality=0, multiple_iterators=True)
            for pileupcolumn in itercol_first_column:
                allreads = pileupcolumn.get_query_names()
                reads_in_breakpoint_left[b] = reads_in_breakpoint_left[b].difference(set(allreads))
                break
    
    #list all combinations of possible breakpoints (2^len(bp_left)*2^len(bp_right)) and the corresponding intersection of reads
    combination_to_set_of_reads = {} #maps combination of breakpoints ([True, False],[True]) -> set of reads
    for read in allreads_anywhere_on_the_window:
        combination = ([False for i in breakpoints_left], [False for i in breakpoints_right])
        for i in range(len(breakpoints_left)):
            if read in reads_in_breakpoint_left[i]:
                combination[0][i] = True
        for i in range(len(breakpoints_right)):
            if read in reads_in_breakpoint_right[i]:
                combination[1][i] = True

        #make sure there is at least one breakpoint on the left and one on the right
        if True in combination[0] and True in combination[1]:
            combination_tuple = (tuple(combination[0]), tuple(combination[1]))
            if combination_tuple not in combination_to_set_of_reads:
                combination_to_set_of_reads[combination_tuple] = set()
            combination_to_set_of_reads[combination_tuple].add(read)

    
    itercol = bamfile_ps.pileup(contig = contig_name,start = start_pos,stop = stop_pos,truncate=True,min_base_quality=0, multiple_iterators=True)
    matrices = {} #map of pileup matrices, one per combination of breakpoints
    list_of_suspicious_positions = {} #map of list of suspicious positions, one per combination of breakpoints
    for i in combination_to_set_of_reads.keys():
        matrices[i] = {}
        list_of_suspicious_positions[i] = []
    
    ref = ps.FastaFile(originalAssembly)
    reference_sequence = ref.fetch(contig_name)
    ref.close()

    for pileupcolumn in itercol:
        if pileupcolumn.nsegments >=5:

            allreads = pileupcolumn.get_query_names()
            allbases_raw = [i.upper() for i in pileupcolumn.get_query_sequences(add_indels=True)]
            
            unique, counts = np.unique(allbases_raw, return_counts=True)
            idx_sort = np.argsort(counts)
            counts = counts[idx_sort]
            unique = unique[idx_sort]
            refbase = reference_sequence[pileupcolumn.reference_pos].upper()

            if len(counts) == 1 or (unique[-1]==refbase and max(4,(1-no_snp_threshold)*pileupcolumn.nsegments) > counts[-2]):
                continue

            #to obtain allbases2, convert allbases_raw: output only the A,C,G,T and *, forget the + and - and the following number
            allbases = []
            for a in allbases_raw:
                if len(a) == 1:
                    allbases.append(a)
                else:
                    allbases.append("".join([i for i in a if i in ['A','C','G','T','*']]))
            tmp_dict = {}
            allbases_2 = []
            for i in range(len(allbases)):
                if allreads[i] in allreads_anywhere_on_the_window :
                    tmp_dict[allreads[i]] = allbases[i]
                    allbases_2.append(allbases[i])

            #count and order the bases
            unique, counts = np.unique(allbases_2, return_counts=True)
            idx_sort = np.argsort(counts)
            unique = unique[idx_sort]
            counts = counts[idx_sort]

            #if there is a base different from the reference, add the position to the list of suspicious positions
            alternative_bases = set()
            for alt_b in range(len(unique)-1,-1,-1):
                if unique[alt_b] != refbase and counts[alt_b] > (1-no_snp_threshold)*len(tmp_dict):
                    alternative_bases.add(unique[alt_b])

            #in tmp_dict, replace the refbase by 1 and the alternative base by 0, and the rest by np.nan
            if len(alternative_bases) > 0 :
                for read in tmp_dict:
                    if tmp_dict[read] == refbase:
                        tmp_dict[read] = 1
                    elif tmp_dict[read] in alternative_bases:
                        tmp_dict[read] = 0
                    else:
                        tmp_dict[read] = np.nan

                for combination in combination_to_set_of_reads.keys():
                    tmp_dict_comb = {}
                    some_reads_here_in_this_combination = False
                    for read in combination_to_set_of_reads[combination]:
                        if read in tmp_dict:
                            tmp_dict_comb[read] = tmp_dict[read]
                            some_reads_here_in_this_combination = True
                    if some_reads_here_in_this_combination:
                        matrices[combination][pileupcolumn.reference_pos] = tmp_dict_comb
                        list_of_suspicious_positions[combination].append(pileupcolumn.reference_pos)

    bamfile_ps.close()

    matrices_list = []
    list_of_suspicious_positions_list = []
    for combination in matrices:
        matrices_list.append(matrices[combination])
        list_of_suspicious_positions_list.append(list_of_suspicious_positions[combination])

    return matrices_list, list_of_suspicious_positions_list

def find_SNPs_in_this_window(pileup, list_of_sus_pos, list_of_reads, max_error_on_a_column):
    """
    Identify Single Nucleotide Polymorphisms (SNPs) in a given window of genomic data.

    Parameters:
    pileup (numpy.ndarray): A matrix of reads (rows) and loci (columns) representing the genomic data.
    list_of_sus_pos (list): list_of_sus_pos[col] = (position of the column col in the original data).
    list_of_reads (list): list_of_reads[row] = (name of the read in the row row in the original data).
    max_error_on_a_column (float): The maximum error rate you can reasonably expect on a column. (sequencing+alignment errors)

    Returns:
    list: A map of validated SNP positions -> list_of_reads_with_alt_alleles_at_this_position

    The function performs the following steps:
    1. Fills missing values in the pileup matrix using K-Nearest Neighbors (KNN) imputation.
    2. Computes pairwise distances between columns
    3. Calculates chi-square statistics or proxies to determine the similarity between columns.
    4. Uses Agglomerative Clustering to group similar columns.
    5. Validates SNPs based on the statistical test
    6. Recovers additional SNPs that correlate well with validated SNPs.

    Note:
    - The function assumes that the input pileup matrix has more than one read and one locus.
    - The function uses a distance threshold of 0.05 for clustering and a p-value threshold of 0.001 for SNP validation.
    """

    ###Filling the missing values using KNN
    number_reads,number_loci = pileup.shape
    pileup_filled = deepcopy(pileup)

    if number_reads>1 and number_loci>1:
        imputer = KNNImputer(n_neighbors=5, copy=False, keep_empty_features=True)
        pileup_filled = imputer.fit_transform(pileup_filled)
        pileup_filled[(pileup_filled>=0.5)] = 1
        pileup_filled[(pileup_filled<0.5)] = 0
    else:
        pileup_filled = pileup  

    #now compute the pairwise number of 0-0, 0-1, 1-0, 1-1 between all the columns with matrix multiplication
    matrix_1_1 = np.dot(pileup_filled.T,pileup_filled)
    matrix_1_0 = np.dot(pileup_filled.T,1-pileup_filled)
    matrix_0_1 = np.dot(1-pileup_filled.T,pileup_filled)
    matrix_0_0 = np.dot(1-pileup_filled.T,1-pileup_filled)

    epsilon = 1e-10  #to avoid dividing 0 by 0
    matrix_0_0 += epsilon
    pairwise_distance_0 = matrix_0_0/(matrix_0_0+0.5*(matrix_0_1+matrix_1_0))
    matrix_1_1 += epsilon
    pairwise_distance_1 = matrix_1_1/(matrix_1_1+0.5*(matrix_0_1+matrix_1_0))
    matrix_0_1 += epsilon
    matrix_1_0 += epsilon

    #now compute the chisquare statistics between all the columns
    pvalues = np.zeros((number_loci,number_loci))
    best_pvalue = 0.0001
    for i in range(number_loci):
        pvalues[i,i] = best_pvalue
        for j in range(i+1,number_loci):
            if pairwise_distance_0[i,j] < 0.9: #if the two columns have not at least 90% similar 0s, dont link them
                pvalues[i,j] = 1
                pvalues[j,i] = 1
            elif pairwise_distance_1[i,j] > 0.9: #if they have highly similar 0s and highly similar 1s, then link them
                pvalues[i,j] = best_pvalue
                pvalues[j,i] = best_pvalue
            else: #this is a bit between the two
                pvalue = max(min(1,stats.chi2_contingency([[matrix_1_1[i,j],matrix_1_0[i,j]],[matrix_0_1[i,j],matrix_0_0[i,j]]])[1]), best_pvalue)
                pvalues[i,j] = pvalue
                pvalues[j,i] = pvalue
                # print("3 rrrrr", matrix_0_0[i,j],matrix_0_1[i,j],matrix_1_0[i,j],matrix_1_1[i,j], pvalue)

    #now we can perform the clustering using AgglomerativeClustering
    agglo_cl = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage = 'complete',distance_threshold=0.05)
    if len(pvalues) > 1:
        agglo_cl.fit(pvalues)
        labels = agglo_cl.labels_
    elif len(pvalues) == 1:
        labels = np.array([0])

    emblematic_snps = [] #keep one snp per cluster to compute correlations

    #extract each label of size at least two and find the biggest rectangle of 0s in the submatrix
    snps_res = {} #maps position -> list_of_reads_with_alt_alleles
    for label in np.unique(labels):
        idx = np.where(labels == label)[0] #idx is the list of positions where the group label is label
        if len(idx) >= 1:
            submatrix = pileup_filled[:,idx]
            #compute the average value of each row and store it in a vector
            row_means = submatrix.mean(axis = 1)
            #find the number of rows of 0s
            nb_rows_0s = len(np.where(row_means <= 0.1)[0])
            nb_rows_with_some_0s = len(np.where(row_means <= 0.9)[0])

            if nb_rows_with_some_0s <= max_error_on_a_column*number_reads and len(idx) == 1: #then the test will always fail anyway
                continue
            if number_reads==0:
                continue 
            error_rate = max(min(max_error_on_a_column, nb_rows_with_some_0s/number_reads),0.05)

            #now perform the statistical test if the obtained rectangle of 0s is significant
            p_value_of_the_column_cluster = statistical_test(nb_rows_0s, number_reads , len(idx), number_loci, error_rate)

            #validate the snps if the p-value is significant
            if p_value_of_the_column_cluster < 0.001:
                emblematic_snps.append(idx[0])
                for i in idx:
                    snps_res[list_of_sus_pos[i]] = set([list_of_reads[i] for i in list(np.where(row_means <= 0.1)[0])])
                
                # print("I found a rectangle of ", nb_rows_0s, " rows of 0s in a cluster of ", len(idx), " columns, among in total ", number_loci, " columns and ", number_reads, " rows with an error rate of ", error_rate, ", which gives a p-value of ", p_value_of_the_column_cluster)
                # for i in range(len(idx)):
                #     for j in range(number_reads):
                #         print(int(submatrix[j,i]), end = '')
                #     print(' : ', idx[i])
                # print('')

                # print("position ", list_of_sus_pos[i], " and reads ", set([list_of_reads[i] for i in list(np.where(row_means <= 0.1)[0])]))


    #as a last step, recover the SNPs which correlate well with validated SNPs, using the already computed pvalues

    for i in emblematic_snps:
        for j in range(number_loci):
            if list_of_sus_pos[j] not in snps_res:
                pvalue = stats.chi2_contingency([[matrix_1_1[i,j],matrix_1_0[i,j]],[matrix_0_1[i,j],matrix_0_0[i,j]]])[1]
                # if list_of_sus_pos[j] == 94 or list_of_sus_pos[j] == 95:
                #     print("I am looking at the correlation of the SNP ", list_of_sus_pos[i], " with the SNP ", list_of_sus_pos[j], " which has a p-value of ", pvalue)
                if pvalue < 0.0001:
                    snps_res[list_of_sus_pos[j]] = set([list_of_reads[alt_row] for alt_row in range(len(pileup)) if pileup_filled[alt_row,j]==0])

    #return the list of validated snps
    return snps_res

def output_VCF_of_this_window(bamfile, originalAssembly, variants, contig_name, vcf_file):
    """
    Generate a VCF file for a given window of a BAM file.
    Parameters:
    bamfile (pysam.AlignmentFile): The BAM file containing the read alignments.
    ref_file (pysam.FastaFile): The reference genome file.
    variants (map): A map of validated SNP positions -> set_of_reads_with_alt_alleles_at_this_position
    start_pos (int): The start position of the window.
    stop_pos (int): The stop position of the window.
    vcf_file (str): The name of the VCF file to write the output.
    Returns:
    None
    """

    #open the output file
    out = open(vcf_file,'w')

    list_of_variants = []

    #go through the variant positions, grouping positions if they are less than 5bp apart
    variantpositions = list(variants.keys())
    variantpositions.sort()

    index_pos_to_group = 0
    groups_of_variants = [] #tuple (start, end, (set) interesting_reads)
    size_of_groups = 5
    interesting_reads = set()
    while index_pos_to_group < len(variantpositions):
        start_pos_group = variantpositions[index_pos_to_group]
        end_pos_group = start_pos_group+1
        interesting_reads = variants[variantpositions[index_pos_to_group]]
        while index_pos_to_group < len(variantpositions)-1 and variantpositions[index_pos_to_group+1] - end_pos_group < size_of_groups:
            end_pos_group = variantpositions[index_pos_to_group+1]+1
            index_pos_to_group += 1
            interesting_reads.union(variants[variantpositions[index_pos_to_group]])
        index_pos_to_group += 1
        groups_of_variants.append((start_pos_group,end_pos_group, interesting_reads))

        interesting_reads = set()

    
    bamfile_ps = ps.AlignmentFile(bamfile,'rb')
    ref_file = ps.FastaFile(originalAssembly)
    index_of_variant_group = 0
    list_of_reads_there = {} #list of the sequence of all the reads in on these positions
    for pileupcolumn in bamfile_ps.pileup(contig = contig_name, start = variantpositions[0], stop = variantpositions[-1]+1,truncate=True, min_base_quality=0):

        if pileupcolumn.reference_pos < groups_of_variants[index_of_variant_group][0]:
            continue

        if pileupcolumn.reference_pos == groups_of_variants[index_of_variant_group][0]:
            #start of the group
            alleles = {}
            list_of_reads_there = {}

        #go through the reads and inventory the alleles
        for pileupread in pileupcolumn.pileups:
            name_of_read = pileupread.alignment.query_name
            if name_of_read in groups_of_variants[index_of_variant_group][2]:
                if pileupread.alignment.query_name not in list_of_reads_there:
                    list_of_reads_there[name_of_read] = ""
                if not pileupread.is_del:
                    insertion = max(0, pileupread.indel)
                    list_of_reads_there[name_of_read] += pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+1+insertion]

        if pileupcolumn.reference_pos == groups_of_variants[index_of_variant_group][1]-1: #exiting the group now, we can output the variants
            
            reference_sequence = ref_file.fetch(contig_name,groups_of_variants[index_of_variant_group][0], groups_of_variants[index_of_variant_group][1])
            if reference_sequence == '':
                reference_sequence = '.'

            alleles = {}
            alleles[reference_sequence] = 0
            for read in list_of_reads_there:
                if list_of_reads_there[read] not in alleles:
                    alleles[list_of_reads_there[read]] = 1
                else:
                    alleles[list_of_reads_there[read]] += 1
            if ''  in alleles:
                alleles['.'] = alleles['']
                alleles.pop('')

            #write the variants to the vcf file: the ref, the main alt and all the alts that are present in 5 reads or more
            out.write(contig_name+'\t'+str(groups_of_variants[index_of_variant_group][0]+1)+'\t.\t'+reference_sequence+'\t')
            #sort the alleles by decreasing order of frequency
            alleles = {k: v for k, v in sorted(alleles.items(), key=lambda item: item[1], reverse = True)}

            #delete the rarest alleles that are included in more frequent alleles: these are often artefacts of homopolymer errors
            to_pop = set()
            first_allele = True
            count_first_allele = 0

            for al in alleles :
                if first_allele:
                    count_first_allele = alleles[al]
                    first_allele = False
                elif alleles[al] >= 5 and alleles[al] >= 0.25*count_first_allele :
                    for al2 in alleles :
                        if alleles[al2] >= 5:

                            if al != reference_sequence and al != al2 and alleles[al] < alleles[al2] and (al == "." or al in al2):
                                to_pop.add(al)
                else:
                    to_pop.add(al)

            for i in to_pop:
                alleles.pop(i)


            first_allele = True
            counts = []
            count_first_allele = 0
            for allele in alleles:
                if allele != reference_sequence:
                    if not first_allele:
                        out.write(',')
                    else :
                        count_first_allele = alleles[allele]
                        first_allele = False
                    out.write(allele)
                    counts.append(alleles[allele])
            out.write('\t.\tAD=')
            for i in range(len(counts)):
                if i > 0:
                    out.write(',')
                out.write(str(counts[i]))
            out.write('\t.\n')
        
            index_of_variant_group += 1

    bamfile_ps.close()
    ref_file.close()

    #close the output file
    out.close()

def call_variants_on_this_window(contig_name, start_pos, end_pos, filtered_col_threshold, no_snp_threshold, bamfile, ref, vcf_file, max_error_on_column):
    """
    Calls variants on a specified window of a genomic contig. Encapsulates the entire process of variant calling.
    Parameters:
    contig_name (str): The name of the contig.
    start_pos (int): The starting position of the window.
    end_pos (int): The ending position of the window.
    filtered_col_threshold (float): The threshold for filtering reads with NaN values.
    no_snp_threshold (int): If a big proportion of reads in a column have the same base, it is considered as a no SNP position.
    bamfile (str): The path to the BAM file.
    ref (str): The reference genome sequence.
    vcf_file (str): The path to the VCF file where results will be written.
    max_error_on_column (float): The maximum error rate you can reasonably expect on a column. (sequencing+alignment errors)
    Returns:
    None
    """
     
    list_of_sus_pos = []

    #print(f'Parsing data on contig {contig_name} {start_pos}<->{end_pos}')
    #time get_data
    time_get_data = time.time()
    list_of_matrices, list_of_list_of_suspicious_positions = get_data(bamfile, ref, contig_name,start_pos,end_pos, no_snp_threshold) #for each window you can have different groups of reads e.g. if some of them align with very big gaps; Here, the lists represent the different groups of reads
    time_get_data2 = time.time() - time_get_data

    all_variants = {} #map of position -> set of reads with alt alleles
    time_call_variants = time.time()

    for dict_of_sus_pos, list_of_sus_pos in zip(list_of_matrices, list_of_list_of_suspicious_positions):

        ###create a matrix from the columns
        df = pd.DataFrame(dict_of_sus_pos) 
        df = df.dropna(axis = 0, thresh = filtered_col_threshold*(len(df.index))) #filter out the reads taht have NaN values in more than x% of the columns (usually corresponds to the reads that are not spanning the whole window)

        ###clustering
        if len(dict_of_sus_pos) > 0 :
            pileup = df.to_numpy(dtype=float)
            list_of_reads= df.index
            time_call_variants = time.time()
            variants_here = find_SNPs_in_this_window(pileup, list_of_sus_pos, list_of_reads, max_error_on_column)

            for pos in variants_here.keys():
                if pos not in all_variants:
                    all_variants[pos] = variants_here[pos]
                else:
                    all_variants[pos] = all_variants[pos].union(variants_here[pos])

    time_call_variants2 = time.time() - time_call_variants
    time_output_vcf2 = 0

    ###append to the vcffile
    if len(all_variants) > 0 :
        time_output_vcf = time.time()
        output_VCF_of_this_window(bamfile, ref, all_variants, contig_name, vcf_file)
        time_output_vcf2 = time.time() - time_output_vcf
    

    print("On contig ", contig_name, " from ", start_pos, " to ", end_pos, "Times: get_data ", time_get_data2, " call_variants ", time_call_variants2, " output_vcf ", time_output_vcf2)

def parse_arguments():
    """Parse the input arguments"""
    argparser = ArgumentParser()

    argparser.add_argument(
        '-r', '--reference', dest='reference', required=True, default='', type=str,
        help='Reference genome in FASTA format',
    )

    argparser.add_argument(
        '-b', '--bam', dest='bam', required=True, default='', type=str,
        help='Alignment file in BAM format',
    )

    argparser.add_argument(
        '-o', '--out-folder', dest='out', required=True, default='', type=str,
        help='Name of the output folder',
    )

    argparser.add_argument(
        '-t', '--threads', dest='threads', required=False, default=1, type=int,
        help='Number of threads to use [1]',
    )

    argparser.add_argument(
        '-e', '--error_rate', required=False, default=0.1, type=float,
        help='(Over)estimated error rate of the reads [0.1]',
    )

    argparser.add_argument(
        '--window', dest='window', required=False, default=5000, type=int,
        help='Size of window to perform read separation (must be at least twice shorter than average read length) [5000]',
    )

    arg = argparser.parse_args()


    return arg.bam, arg.out, arg.window, arg.reference, arg.threads, arg.error_rate


if __name__ == '__main__':

    print('metaCaller version ', __version__)
    print(" ".join(sys.argv))
    time_start = time.time()

    window = 5000
    no_snp_threshold = 0.95
    filtered_col_threshold = 0.9
    min_row_quality = 5
    min_col_quality = 3

    bamfile, out, window, originalAssembly, num_threads, error_rate = parse_arguments()
    #mkdir out
    if not os.path.exists(out):
        os.makedirs(out)

    max_error_rate_on_column = max(0.1,min(error_rate*5, 0.5))

    #check that the bam file is indexed, else index it
    if not os.path.exists(bamfile+'.bai'):
        print('Indexing the bam file')
        os.system('samtools index '+bamfile)
    
    #mkdir out/tmp
    if not os.path.exists(out+'/tmp'):
        os.makedirs(out+'/tmp')

    tmp_dir = out+'/tmp'

    #create the vcf file
    vcf_file = out+'/variants.vcf'
    if os.path.exists(vcf_file):
        os.remove(vcf_file)
    o = open(vcf_file,'w')
    o.write('##fileformat=VCFv4.2\n')
    o.write('##source=metaCaller\n')
    o.write('##reference='+originalAssembly+'\n')
    o.write('##INFO=<ID=AD,Number=A,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\tFILTER\n')
    o.close()

    ref = ps.FastaFile(originalAssembly)
    #if the reference is not indexed, index it
    if not os.path.exists(originalAssembly+'.fai'):
        print('Indexing the reference')
        os.system('samtools faidx '+originalAssembly)
    ref.close()


    #start calling variants
    start = time.time()
    bamfile_ps = ps.AlignmentFile(bamfile,'rb')
    contigs = (bamfile_ps.header.to_dict())['SQ']
    if len(contigs) == 0:
        print('ERROR: No contigs found when parsing the BAM file, check the bam file and the indexation of the bam file')
        sys.exit(1)

    windows_description = [] #list of the windows that will be processed, as tuples (contig_name, start_pos, stop_pos, tmp_vcf_file)
    for num in range(0,len(contigs)):
        contig_name = contigs[num]['SN']
        contig_length = contigs[num]['LN']
        for start_pos in range(0,contig_length,window):
            if start_pos+window <= contig_length:
                windows_description.append((contig_name,start_pos,start_pos+window, out+'/tmp/'+contig_name+'_'+str(start_pos)+'.vcf'))
            else:
                windows_description.append((contig_name,start_pos,contig_length, out+'/tmp/'+contig_name+'_'+str(start_pos)+'.vcf'))

    bamfile_ps.close()

    # windows_description = [i for i in windows_description if i[1]< 5000]
    #windows_description = [("CP075493.1_1", 12320000, 12325000, out+'/tmp/'+"CP075493.1"+'_'+str(70000)+'.vcf')] #DEBUGbc
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(call_variants_on_this_window, window[0], window[1], window[2], filtered_col_threshold, no_snp_threshold, bamfile, originalAssembly, window[3], max_error_rate_on_column) for window in windows_description]
        for future in futures:
            future.result()

    #concatenate all the vcf files in the order of the windows
    for window in windows_description:
        if os.path.exists(window[3]):
            os.system('cat '+window[3]+' >> '+vcf_file)
            os.remove(window[3])

    #print the time now, and the total time taken
    print("["+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())+"] Total time taken: ", time.time()-time_start)
    
