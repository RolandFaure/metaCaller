"""
SNPsnoop
Authors: Roland Faure, based on a previous program (strainminer) by Minh Tam Truong Khac
"""

__version__ = '0.4.0'

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
import shutil

fasta_lock = threading.Lock()
bam_lock = threading.Lock()
 
def log_comb(n, k):
    """Compute the logarithm of the combination n choose k. This is useful to avoid numerical errors when n is large"""
    if n == 0 or k == 0:
        return 0 #which is log 1
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

def get_data(bamfile, originalAssembly, contig_name,start_pos,stop_pos, no_snp_threshold = 0.95, max_error_rate_on_a_column = 0.1, variants_already_called = {}):
    """
    Extracts a binary matrix of reads and columns from a BAM file, identifying suspicious positions in a specified contig region.
    Args:
        bamfile (str): Path to the sorted and indexed BAM file.
        originalAssembly (str): Path to the reference genome assembly in FASTA format.
        contig_name (str): Name of the contig to analyze.
        start_pos (int): Start position of the region to analyze.
        stop_pos (int): Stop position of the region to analyze.
        no_snp_threshold (float, optional): Threshold for considering a position as suspicious. Default is 0.95.
        max_error_rate_on_a_column (float, optional): Maximum error rate allowed for a column to be considered. Default is 0.1.
        variants_already_called (dict, optional): Dictionary of variants already called: positions -> set of reads. Default is {}.
    Returns:
        tuple: A tuple containing:
            - matrix (dict): A dictionary where keys are positions and values are dictionaries mapping read names to binary values (1 for reference base, 0 for main alternative base, np.nan for others).
            - list_of_suspicious_positions (list): A list of positions identified as suspicious.
    """
    
    bamfile_ps = ps.AlignmentFile(bamfile,'rb')
    obvious_snps = {}

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

    breakpoints_left = sorted(breakpoints_left)
    breakpoints_right = sorted(breakpoints_right)
    
    allreads_anywhere_on_the_window = set()
    reads_in_breakpoint_right = [set() for i in breakpoints_right]
    reads_in_breakpoint_right = [set() for _ in breakpoints_right]
    reads_in_breakpoint_left = [set() for _ in breakpoints_left]

    for read in bamfile_ps.fetch(contig=contig_name, start=start_pos, stop=stop_pos):
        if read.is_unmapped:
            continue

        read_start = read.reference_start
        read_end = read.reference_end

        # Assign to the closest breakpoint to the right for breakpoints_left
        found_bp = False
        for b, bp_left in reversed(list(enumerate(breakpoints_left))):
            if read_end >= bp_left:
                reads_in_breakpoint_left[b].add(read.query_name)
                allreads_anywhere_on_the_window.add(read.query_name)
                found_bp = True
                break
        if not found_bp:
            reads_in_breakpoint_left[0].add(read.query_name)
            allreads_anywhere_on_the_window.add(read.query_name)

        # Assign to the closest breakpoint to the left for breakpoints_right
        found_bp = False
        for b, bp_right in enumerate(breakpoints_right):
            if read_start <= bp_right:
                reads_in_breakpoint_right[b].add(read.query_name)
                allreads_anywhere_on_the_window.add(read.query_name)
                found_bp = True
                break
        if not found_bp:
            reads_in_breakpoint_right[-1].add(read.query_name)
            allreads_anywhere_on_the_window.add(read.query_name)
    
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

    
    itercol = bamfile_ps.pileup(contig=contig_name, start=start_pos, stop=stop_pos, truncate=True, min_base_quality=0, multiple_iterators=True, stepper="all", flag_filter=4 | 256 | 512 | 1024)
    matrices = {} #map of pileup matrices, one per combination of breakpoints
    list_of_suspicious_positions = {} #map of list of suspicious positions, one per combination of breakpoints
    for i in combination_to_set_of_reads.keys():
        matrices[i] = {}
        list_of_suspicious_positions[i] = []
    
    ref = ps.FastaFile(originalAssembly)
    reference_sequence = ref.fetch(contig_name)
    ref.close()

    read_wise_error_rate = {}
    read_wise_nan_values = {}
    read_wise_number_of_bases = {}

    for pileupcolumn in itercol:
        if pileupcolumn.nsegments >=1:#filter on the coverage

            allreads = pileupcolumn.get_query_names()
            allbases_raw = [i.upper() for i in pileupcolumn.get_query_sequences(add_indels=True)]
            refbase = reference_sequence[pileupcolumn.reference_pos].upper()
            
            unique, counts = np.unique(allbases_raw, return_counts=True)
            idx_sort = np.argsort(counts)
            counts = counts[idx_sort]
            unique = unique[idx_sort]

            #compute the error rate of each read
            for read, base in zip(allreads, allbases_raw):
                if read not in read_wise_error_rate:
                    read_wise_error_rate[read] = 0
                    read_wise_number_of_bases[read] = 0
                if base.upper() != refbase and base.upper() in ['A', 'C', 'G', 'T']:
                    read_wise_error_rate[read] += 1
                if base.upper() in ['A', 'C', 'G', 'T']:
                    read_wise_number_of_bases[read] += 1    
            
            if (len(counts) == 1 and unique[-1] ==refbase) or (len(counts) > 1 and unique[-1]==refbase and max(2,(1-no_snp_threshold)*pileupcolumn.nsegments) > counts[-2]):
                continue

            #to obtain allbases2, convert allbases_raw: output only the A,C,G,T and *, forget the + and - and the following number
            allbases = []
            for a in allbases_raw:
                if len(a) == 1 and a in ['A','C','G','T','*']:
                    allbases.append(a)
                else:
                    allbases.append("".join([i for i in a if i in ['A','C','G','T','*']]))

            #remove already called bases and replace them by the ref
            if pileupcolumn.reference_pos in variants_already_called:
                # If the position is already called, replace the base by the reference for those reads
                for idx, read in enumerate(allreads):
                    if read in variants_already_called[pileupcolumn.reference_pos]:
                        allbases[idx] = refbase
                
            tmp_dict = {}
            allbases_2 = []
            for i in range(len(allbases)):
                if allreads[i] in allreads_anywhere_on_the_window :
                    tmp_dict[allreads[i]] = allbases[i]
                    allbases_2.append(allbases[i])

            # Remove empty strings from allbases_2
            allbases_2 = [base for base in allbases_2 if base != '']

            #count and order the bases
            unique, counts = np.unique(allbases_2, return_counts=True)
            idx_sort = np.argsort(counts)
            unique = unique[idx_sort]
            counts = counts[idx_sort]
        
            # print("couoqn", counts)

            #if there is a base different from the reference, add the position to the list of suspicious positions or to obvious snps if it's really too obvious
            alternative_bases = set()
            for alt_b in range(len(unique)-1,-1,-1):
                mean = max_error_rate_on_a_column * len(allbases_2)
                std_dev = math.sqrt(max_error_rate_on_a_column * (1 - max_error_rate_on_a_column) *len(allbases_2))
                if unique[alt_b] != refbase and unique[alt_b] != '' and counts[alt_b] >= 2 and 1 - stats.norm.cdf(counts[alt_b], loc=mean, scale=std_dev) < 0.001:  # Check if alternative bases are obviously a SNP
                    # print(f"selecting obvious SNP at position {pileupcolumn.reference_pos} with base {unique[alt_b]} (count: {counts[alt_b]}), error rate and number of reads : ", max_error_rate_on_a_column, " * ", len(allbases_2))
                    obvious_snps[pileupcolumn.reference_pos] = set([allreads[row] for row in range(len(allreads)) if allbases[row] == unique[alt_b]])
                    alternative_bases.add(unique[alt_b])
                elif unique[alt_b] != refbase and unique[alt_b] != '' and counts[alt_b] > (1-no_snp_threshold)*len(tmp_dict) and len(alternative_bases) == 0 :
                    alternative_bases.add(unique[alt_b])

            #in tmp_dict, replace the refbase by 1 and the alternative base by 0, and the rest by np.nan
            if len(alternative_bases) > 0 :
                for read in tmp_dict:
                    if tmp_dict[read] in alternative_bases:
                        tmp_dict[read] = 0
                    else:
                        tmp_dict[read] = 1

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
    read_error_rates = []
    errors_of_reads_list = []
    for combination in matrices:
        matrices_list.append(matrices[combination])
        list_of_suspicious_positions_list.append(list_of_suspicious_positions[combination])
        error_rates = []
        for read in combination_to_set_of_reads[combination]:
            if read in read_wise_error_rate :
                error_rates.append(read_wise_error_rate[read] / (read_wise_number_of_bases[read]+1))
        read_error_rates.append(error_rates)

    # # Debugging: Print the pileup matrices for positions between 13895 and 13900
    # for combination, matrix in matrices.items():
    #     print(f"Combination: {combination}")
    #     for pos in range(72660, 74670):
    #         if pos in matrix:
    #             print(f"Position {pos}: {''.join('*' if v is np.nan else str(int(v)) for k, v in sorted(matrix[pos].items()))}")

    return matrices_list, list_of_suspicious_positions_list, obvious_snps, read_error_rates

def find_SNPs_in_this_window(pileup, list_of_sus_pos, list_of_reads, max_error_on_a_column, error_rates_of_reads):
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

    if number_reads>1 and number_loci>1 :
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
            if matrix_1_1[i,j] > matrix_0_1[i,j] and pairwise_distance_0[i,j] < 0.7: #if the two columns have not at least 70% similar 0s, dont link them
                pvalues[i,j] = 1
                pvalues[j,i] = 1
            elif pairwise_distance_0[i,j] > 0.9 and pairwise_distance_1[i,j] > 0.9: #if they have highly similar 0s and highly similar 1s, then link them
                pvalues[i,j] = best_pvalue
                pvalues[j,i] = best_pvalue
            else: #this is a bit between the two
                # Ensure no zero values in the contingency table
                epsilon_chi = 1e-10
                contingency = [
                    [max(matrix_1_1[i,j], epsilon_chi), max(matrix_1_0[i,j], epsilon_chi)],
                    [max(matrix_0_1[i,j], epsilon_chi), max(matrix_0_0[i,j], epsilon_chi)]
                ]
                pvalue = max(min(1,stats.chi2_contingency(contingency)[1]), best_pvalue)
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
    # print([(int(labels[i]), list_of_sus_pos[i]) for i in range(len(labels))], "ooiuo")
    for label in np.unique(labels):
        idx = np.where(labels == label)[0] #idx is the list of positions where the group label is label
        if len(idx) >= 1:

            submatrix = pileup_filled[:,idx]

            # Ensure all columns in the group are almost identical (and not correlating like 1001 correlates with 0110)
            toggled_columns = [False for i in range(len(idx))]
            for j in range(1, len(idx)):
                col_i = submatrix[:, 0]
                col_j = submatrix[:, j]
                if np.sum(col_i != col_j) > 0.5 * number_reads:  # If more than 50% of rows differ
                    # Flip all bits in one of the columns to make them identical
                    submatrix[:, j] = 1 - submatrix[:, j]
                    toggled_columns[j] = True

            #compute the average value of each row and store it in a vector
            row_means = submatrix.mean(axis = 1)
            #find the number of rows of 0s
            nb_rows_0s = len(np.where(row_means <= 0.1)[0])
            nb_rows_with_some_0s = len(np.where(row_means <= 0.9)[0])

            if nb_rows_with_some_0s <= max_error_on_a_column*number_reads and len(idx) == 1: #then the test will always fail anyway
                continue
            if nb_rows_0s<=1:
                continue 

            # Extract error rates of reads having 0s
            error_rates_of_reads_with_zeros = []
            for row in range(number_reads):
                if row_means[row] <= 0.1:  # Same condition as used for nb_rows_0s
                    error_rates_of_reads_with_zeros.append(error_rates_of_reads[row])

            max_error_rate_of_reads = np.max(error_rates_of_reads_with_zeros)
            error_rate = max(min(max_error_on_a_column, nb_rows_with_some_0s/number_reads),0.05, max_error_rate_of_reads)

            # Do not count stretches of SNPs as distinct snps (can create false positives in deletions)
            number_of_separated_columns = 0
            pos_of_previous_snp = -2
            for i in range(len(idx)):
                if list_of_sus_pos[idx[i]] - pos_of_previous_snp > 1:
                    number_of_separated_columns += 1
                pos_of_previous_snp = list_of_sus_pos[idx[i]]

            #if there is only two rows, augment the number of columns (since if only one row is mutated the column was not selected)
            corrected_number_of_loci = number_loci - (len(idx)-number_of_separated_columns)
            if nb_rows_0s == 2:
                length_of_stretch = list_of_sus_pos[idx[-1]] - list_of_sus_pos[idx[0]]
                corrected_number_of_loci = max(corrected_number_of_loci, int( max_error_rate_of_reads * length_of_stretch ))

            #now perform the statistical test if the obtained rectangle of 0s is significant
            p_value_of_the_column_cluster = statistical_test(nb_rows_0s, number_reads , number_of_separated_columns, corrected_number_of_loci, error_rate)

            #validate the snps if the p-value is significant
            if p_value_of_the_column_cluster < 0.001:
                mean_snp_vector = np.mean(submatrix, axis=1)  # already computed above
                emblematic_snps.append((idx, mean_snp_vector))

                for idx_in_idx in range(len(idx)):
                    pos = idx[idx_in_idx]
                    # print("snp at ", list_of_sus_pos[i])
                    if toggled_columns[idx_in_idx]:
                        snps_res[list_of_sus_pos[pos]] = set([list_of_reads[i] for i in list(np.where(row_means > 0.5)[0])])
                    else:
                        snps_res[list_of_sus_pos[pos]] = set([list_of_reads[i] for i in list(np.where(row_means <= 0.1)[0])])

                # Example: If position 1612 is in the current rectangle, print a message
                if 12169 in [list_of_sus_pos[idx[i]] for i in range(len(idx))]:
                    print("Read names in rectangle:", [list_of_reads[j] for j in range(number_reads)])
                    print("I found a rectangle of ", nb_rows_0s, " rows of 0s in a cluster of ", len(idx), " columns, among in total ", corrected_number_of_loci, " columns and ", number_reads, " rows with an error rate of ", error_rate, ", which gives a p-value of ", p_value_of_the_column_cluster, " (label ", label, ")")
                    for i in range(len(idx)):
                        for j in range(number_reads):
                            print(int(submatrix[j,i]), end = '')
                        print(' : ', list_of_sus_pos[idx[i]])
                    print('')

                # print("position ", list_of_sus_pos[i], " and reads ", set([list_of_reads[i] for i in list(np.where(row_means <= 0.1)[0])]))


    #as a last step, recover the SNPs which correlate well with validated SNPs, using the mean SNP vector of each cluster

    for idxs, mean_snp_vector in emblematic_snps:
        for j in range(number_loci):
            if list_of_sus_pos[j] not in snps_res:
                # Compute correlation between mean_snp_vector and column j
                col_j = pileup_filled[:, j]
                # Only consider rows without NaN in mean_snp_vector or col_j
                valid_mask = ~np.isnan(mean_snp_vector) & ~np.isnan(col_j)
                if np.sum(valid_mask) == 0:
                    continue
                
                pvalue = stats.chi2_contingency([[np.sum((mean_snp_vector[valid_mask] < 0.5) & (col_j[valid_mask] < 0.5)),
                                                  np.sum((mean_snp_vector[valid_mask] < 0.5) & (col_j[valid_mask] >= 0.5))],
                                                 [np.sum((mean_snp_vector[valid_mask] >= 0.5) & (col_j[valid_mask] < 0.5)),
                                                  np.sum((mean_snp_vector[valid_mask] >= 0.5) & (col_j[valid_mask] >= 0.5))]])[1]
                # # Debug print for specific position
                # if list_of_sus_pos[j] == 12170:
                #     print("I am looking at the correlation of the mean SNP vector with the SNP ", list_of_sus_pos[j], " and p-value ", pvalue)
                # Rescue based on pvalue
                if (pvalue < 0.001 / max(1, len(emblematic_snps)) or pvalue<0.000001)and np.sum(col_j[valid_mask] < 0.5) > 0:
                    snps_res[list_of_sus_pos[j]] = set([list_of_reads[alt_row] for alt_row in range(len(pileup)) if pileup_filled[alt_row, j] == 0 and mean_snp_vector[alt_row] < 0.5])
    #return the list of validated snps
    return snps_res

def output_VCF_of_this_window(bamfile, originalAssembly, variants, contig_name, vcf_file, benchmark=False):
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
    if benchmark :
        size_of_groups = 0
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
    for pileupcolumn in bamfile_ps.pileup(contig = contig_name, start = variantpositions[0], stop = variantpositions[-1]+1,truncate=True, min_base_quality=0, stepper="all"):

        if pileupcolumn.reference_pos < groups_of_variants[index_of_variant_group][0]:
            continue

        if pileupcolumn.reference_pos == groups_of_variants[index_of_variant_group][0]:
            #start of the group
            alleles = {}
            list_of_reads_there = {}

        #go through the reads and inventory the alleles
        already_seen_name_in_the_column = set()
        for pileupread in pileupcolumn.pileups:

            name_of_read = pileupread.alignment.query_name
            while name_of_read in already_seen_name_in_the_column :
               name_of_read = name_of_read+"%"
            already_seen_name_in_the_column.add(name_of_read) 
            if name_of_read.rstrip('%') in groups_of_variants[index_of_variant_group][2]: #only take the read where there are mutations
                if name_of_read not in list_of_reads_there:
                    list_of_reads_there[name_of_read] = ""
                if not pileupread.is_del:
                    insertion = max(0, pileupread.indel)
                    seq = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+1+insertion]
                    # Skip the read if it contains a letter that is not A, C, G, T
                    if len(seq) == 1 and seq[0] not in "ACGT":
                        continue
                    list_of_reads_there[name_of_read] += seq

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

            #sort the alleles by decreasing order of frequency
            alleles = {k: v for k, v in sorted(alleles.items(), key=lambda item: item[1], reverse = True)}
            # if 1308 <= groups_of_variants[index_of_variant_group][0] <= 1314:
            #     print(f"Position: {groups_of_variants[index_of_variant_group][0]}, Alleles: {alleles}")

            #delete the rarest alleles that are included in more frequent alleles: these are often artefacts of homopolymer errors. Also delete if not frequent enough
            to_pop = set()
            first_allele = True
            count_first_allele = 0

            for al in alleles :
                if al == reference_sequence :
                    to_pop.add(al)
                    continue

                if first_allele:
                    count_first_allele = alleles[al]
                    first_allele = False
                elif alleles[al] >= 15 and alleles[al] >= 0.1*count_first_allele : #if allele is too secondary, don't output as variant
                    for al2 in alleles :
                        if alleles[al2] >= 5:
                            if al != reference_sequence and al != al2 and alleles[al] < alleles[al2] and (al == "." or al in al2): #if included in more numerous allele, dont output (homopolymer errors :/)
                                to_pop.add(al)
                else:
                    to_pop.add(al)

            for i in to_pop:
                alleles.pop(i)
            
            if len(alleles) == 0:
                index_of_variant_group += 1
                continue

            #write the variants to the vcf file: the ref, the main alt and all the alts that are present in 5 reads or more
            out.write(contig_name+'\t'+str(groups_of_variants[index_of_variant_group][0]+1)+'\t.\t'+reference_sequence+'\t')

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
            out.write('\t.\tDP=')
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

def call_variants_on_this_window(contig_name, start_pos, end_pos, filtered_col_threshold, no_snp_threshold, bamfile, ref, vcf_file, max_error_on_column, benchmark=False):
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
    bencmark (bool): if doing a benchmark against other assemblers, do not group close variants
    Returns:
    None
    """
     
    all_variants = {} #map of position -> set of reads with alt alleles
    number_of_new_variant_found = 1 #set to 1 to make it go through the loop
    number_of_iterations = 0
    time_get_data2 = 0
    time_output_vcf2 = 0
    time_call_variants2 = 0

    #do several iteration to unveil variant masked by other variants at the same position
    while number_of_new_variant_found > 0:
    
        #print(f'Parsing data on contig {contig_name} {start_pos}<->{end_pos}')
        number_of_new_variant_found = 0
        number_of_iterations += 1
        if number_of_iterations == 10:
            #this is a complex region, even a bit weird, but whatever let's break
            break

        #time get_data
        time_get_data = time.time()
        list_of_matrices, list_of_list_of_suspicious_positions, obvious_snps, error_rates\
                = get_data(bamfile, ref, contig_name,start_pos,end_pos, no_snp_threshold, max_error_on_column, all_variants) #for each window you can have different groups of reads e.g. if some of them align with very big gaps; Here, the lists represent the different groups of reads
        time_get_data2 += time.time() - time_get_data
        # # Print obvious SNPs for debugging
        # print("Obvious SNPs:")
        # for pos, reads in obvious_snps.items():
        #     print(f"Position: {pos}, Reads: {reads}")

        # print(list_of_list_of_suspicious_positions, "uisofuo")

        time_call_variants = time.time()

        for dict_of_sus_pos, list_of_sus_pos, error_rates_combination in zip(list_of_matrices, list_of_list_of_suspicious_positions, error_rates):

            ###create a matrix from the columns
            df = pd.DataFrame(dict_of_sus_pos) 
            
            # # Debugging: Print the pileup at position 68803
            # if 36 in dict_of_sus_pos:
            #     print(f"Pileup at position 68803: {''.join('*' if (v!=0 and v!=1) else str(int(v)) for v in dict_of_sus_pos[36].values())}")
            # df = df.dropna(axis=0, thresh=filtered_col_threshold * len(df.columns))  # Filter out rows (reads) with too many NaN values
            df = df.fillna(1) #na values are not defined, they should not help call variants

            ###clustering
            if len(dict_of_sus_pos) > 0 :
                pileup = df.to_numpy(dtype=float)
                list_of_reads= df.index
                time_call_variants = time.time()
                variants_here = find_SNPs_in_this_window(pileup, list_of_sus_pos, list_of_reads, max_error_on_column, error_rates_combination)
                number_of_new_variant_found = len(variants_here)
                # print("new variants found: ", sorted(variants_here.keys()))
                # sys.exit()

                for pos in variants_here.keys():
                    if pos not in all_variants:
                        all_variants[pos] = variants_here[pos]
                    else:
                        all_variants[pos] = all_variants[pos].union(variants_here[pos])

        # # Merge obvious SNPs with less obvious SNPs
        # for pos, reads in obvious_snps.items():
        #     if pos not in all_variants:
        #         all_variants[pos] = reads

        time_call_variants2 += time.time() - time_call_variants
        time_output_vcf2 = 0

    ###append to the vcffile
    if len(all_variants) > 0 :
        time_output_vcf = time.time()
        output_VCF_of_this_window(bamfile, ref, all_variants, contig_name, vcf_file, benchmark)
        time_output_vcf2 += time.time() - time_output_vcf 

    print("On contig ", contig_name, " from ", start_pos, " to ", end_pos, "Times: get_data ", int(time_get_data2), "s call_variants ", int(time_call_variants2), "s output_vcf ", int(time_output_vcf2), "s; in ", number_of_iterations, " iterations")

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
        '--window', dest='window', required=False, default=0, type=int,
        help='Size of window to perform read separation. Must be at least twice shorter than average read length. 0 for auto [0]',
    )

    argparser.add_argument(
        '--benchmark', dest='benchmark', action='store_true', required=False, default=False,
        help='Use if benchmarking to compare against other variant callers [False]',
    )

    arg = argparser.parse_args()


    return arg.bam, arg.out, arg.window, arg.reference, arg.threads, arg.benchmark


if __name__ == '__main__':

    print('metaCaller version ', __version__)
    print(" ".join(sys.argv))
    time_start = time.time()

    window = 0
    no_snp_threshold = 0.95
    filtered_col_threshold = 0.9
    min_row_quality = 5
    min_col_quality = 3

    bamfile, out, window, originalAssembly, num_threads, benchmark = parse_arguments()
    #mkdir out
    if not os.path.exists(out):
        os.makedirs(out)

    # Compute the error rate from the bam file
    total_bases = 0
    mismatched_bases = 0
    bamfile_ps = ps.AlignmentFile(bamfile, 'rb')
    for read in bamfile_ps.fetch():
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary and read.has_tag('NM'):  # NM tag indicates the number of mismatches
            mismatched_bases += read.get_tag('NM')
            total_bases += read.query_length
    bamfile_ps.close()
    if total_bases > 0:
        error_rate = mismatched_bases / total_bases
    else:
        error_rate = 0.1  # Fallback to default if no reads are found
    print(f"Computed error rate+divergence from the bam file: {error_rate}.")

    max_error_rate_on_column = max(0.1,min(error_rate*3, 0.5))

    #check that the bam file is indexed, else index it
    if not os.path.exists(bamfile+'.bai'):
        print('Indexing the bam file')
        os.system('samtools index '+bamfile)

    # Read all the contigs that exist in the BAM file
    bamfile_ps = ps.AlignmentFile(bamfile, 'rb')
    contigs = bamfile_ps.header['SQ']
    bamfile_ps.close()
    
    # Remove out/tmp if it exists, then recreate it
    if os.path.exists(out + '/tmp'):
        shutil.rmtree(out + '/tmp')
    os.makedirs(out + '/tmp')

    tmp_dir = out+'/tmp'

    #create the vcf file
    vcf_file = out+'/variants.vcf'
    if os.path.exists(vcf_file):
        os.remove(vcf_file)
    o = open(vcf_file,'w')
    o.write('##fileformat=VCFv4.2\n')
    o.write('##source=metaCaller\n')
    o.write('##reference='+originalAssembly+'\n')
    for contig in contigs:
        o.write(f"##contig=<ID={contig['SN']},length={contig['LN']}>\n")
    o.write('##INFO=<ID=DP,Number=.,Type=Integer,Description="Allelic depths for the variant in the order listed">\n')

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
        window_here = window
        if window == 0:
            average_read_length = 0
            read_count = 0
            for read in bamfile_ps.fetch(contig_name):
                average_read_length += read.query_length
                read_count += 1

            if read_count > 0:
                average_read_length //= read_count
            else:
                average_read_length = 200  # Default to 100 if no reads are found, we dont care anyway
            window_here = min(10000,max(200, average_read_length // 2))  # Ensure a window size between 200 and 10000
        for start_pos in range(0,contig_length,window_here):
            if start_pos+window_here <= contig_length:
                windows_description.append((contig_name,start_pos,start_pos+window_here, out+'/tmp/'+contig_name+'_'+str(start_pos)+'.vcf'))
            else:
                windows_description.append((contig_name,start_pos,contig_length, out+'/tmp/'+contig_name+'_'+str(start_pos)+'.vcf'))

    bamfile_ps.close()

    # Filter windows for contig_9843 encompassing position 35042
    windows_description = [window for window in windows_description if window[0] == 'edge_349' and window[1] <= 12173 < window[2]]
    print(windows_description)
    print("DEBUUUUUUG")

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(call_variants_on_this_window, window[0], window[1], window[2], filtered_col_threshold, no_snp_threshold, bamfile, originalAssembly, window[3], max_error_rate_on_column, benchmark) for window in windows_description]
        for future in futures:
            future.result()

    # Rescue easy SNPs using longshot
    longshot_vcf = tmp_dir + '/longshot.vcf'
    longshot_cmd = f'longshot --bam {bamfile} --ref {originalAssembly} --out {longshot_vcf}'
    print(f"Running longshot to rescue easy SNPs: {longshot_cmd}")
    os.system(longshot_cmd)

    #merge all VCFs
    tmp_vcf = tmp_dir + '/tmp_merge.vcf'
    # Format longshot VCF: keep first five columns, then add '.' and 'DP='
    with open(longshot_vcf, 'r') as infile, open(tmp_vcf, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    # Write first five columns, then '.', then 'DP='
                    outfile.write('\t'.join(fields[:5]) + '\t.\t'+ fields[7].split(';')[0] +'\n')
    for window in windows_description:
        if os.path.exists(window[3]):
            os.system('cat '+window[3]+' >> '+tmp_vcf )
            # os.remove(window[3])

    # Sort tmp_vcf by key (CHROM, POS) and remove duplicate keys
    sorted_vcf = out + '/variants.sorted.vcf'
    seen_keys = set()
    with open(tmp_vcf, 'r') as infile, open(sorted_vcf, 'w') as outfile:
        header_lines = []
        variant_lines = []
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) > 1:
                    key = (fields[0], fields[1])
                    if key not in seen_keys:
                        seen_keys.add(key)
                        variant_lines.append(line)
        # Write header
        for h in header_lines:
            outfile.write(h)
        # Sort variant lines by CHROM and POS
        variant_lines.sort(key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))
        for v in variant_lines:
            outfile.write(v)
    print(f"Sorted and deduplicated VCF written to {sorted_vcf}")

    # Concatenate sorted VCF to main VCF file using cat
    os.system(f"cat {sorted_vcf} >> {vcf_file}")

    #print the time now, and the total time taken
    print("["+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())+"] Total time taken: ", time.time()-time_start)
    
