"""
metaCaller
Authors: Roland Faure, based on a previous program by Minh Tam Truong Khac
"""

__version__ = '0.1.4'

import pandas as pd 
import numpy as np
import math
import sys
import os
import pysam

import gurobipy as grb
import pysam as ps

from sklearn.cluster import FeatureAgglomeration
from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
from sklearn.metrics import pairwise_distances

import scipy.stats as stats


import time
from argparse import ArgumentParser

def log_comb(n, k):
    """Compute the logarithm of the combination n choose k. This is useful to avoid numerical errors when n is large"""
    if n < 1000:
        return np.log(float(math.comb(n, k)))
    else:
        return np.sum(np.log(np.arange(n, n - k, -1))) - np.sum(np.log(np.arange(1, k + 1)))


def statistical_test(a,n,b,m,error_rate):
    return np.exp(log_comb(n,a) + log_comb(m,b) + (a * b) * np.log(error_rate))


def get_data(file, ref_file, contig_name,start_pos,stop_pos, no_snp_threshold = 0.95):
    #INPUT: a SORTED and INDEXED Bam file
    # Go through a window w on a contig and select suspicious positions
    #OUTPUT: a binary matrix of reads and columns, where the columns are the positions and the rows are the reads. In the matrix, 1 means the base is the same as the reference, 0 means the base is the main alternative base, and np.nan means the base is not the reference nor the main alternative base
    itercol = file.pileup(contig = contig_name,start = start_pos,stop = stop_pos,truncate=True,min_base_quality=10)
    matrix = {}
    list_of_suspicious_positions = []
    reference_sequence = ref_file.fetch(contig_name)

    for pileupcolumn in itercol:
        if pileupcolumn.nsegments >=5:

            allbases_raw = [i.upper() for i in pileupcolumn.get_query_sequences(add_indels=True)]
            allreads = pileupcolumn.get_query_names()
            unique, counts = np.unique(allbases_raw, return_counts=True)
            idx_sort = np.argsort(counts)
            counts = counts[idx_sort]

            if len(counts) == 1 or (1-no_snp_threshold)*pileupcolumn.nsegments > counts[-2]:
                continue

            #to obtain allbases2, convert allbases_raw: output only the A,C,G,T and *, forget the + and - and the following number
            allbases = []
            for a in allbases_raw:
                if len(a) == 1:
                    allbases.append(a)
                else:
                    allbases.append("".join([i for i in a if i in ['A','C','G','T','*']]))
            tmp_dict = {}
            for i in range(len(allbases)):
                tmp_dict[allreads[i]] = allbases[i]

            #count and order the bases
            unique, counts = np.unique(allbases, return_counts=True)
            idx_sort = np.argsort(counts)
            unique = unique[idx_sort]
            counts = counts[idx_sort]

            #if there is a base different from the reference, add the position to the list of suspicious positions
            refbase = reference_sequence[pileupcolumn.reference_pos]
            alternative_base = unique[-1]
            count_alternative_base = counts[-1]
            if alternative_base == refbase and len(unique) > 1:
                alternative_base = unique[-2]
                count_alternative_base = counts[-2]

            #in tmp_dict, replace the refbase by 1 and the alternative base by 0, and the rest by np.nan
            if count_alternative_base >= (1-no_snp_threshold)*pileupcolumn.nsegments:
                for read in tmp_dict:
                    if tmp_dict[read] == refbase:
                        tmp_dict[read] = 1
                    elif tmp_dict[read] == alternative_base:
                        tmp_dict[read] = 0
                    else:
                        tmp_dict[read] = np.nan

                matrix[pileupcolumn.reference_pos] = tmp_dict
                list_of_suspicious_positions.append(pileupcolumn.reference_pos)

    return matrix, list_of_suspicious_positions

def call_variants_in_this_window(pileup, list_of_sus_pos):
    '''
    INPUT: pileup: a matrix of reads and columns
           row_homogeneity: the minimum homogeneity of the rows to be considered as a cluster

    '''

    ###Filling the missing values using KNN
    number_reads,number_loci = pileup.shape
    if number_reads>1 and number_loci>1:
        imputer = KNNImputer(n_neighbors=5)
        matrix = imputer.fit_transform(pileup)
        matrix[(matrix>=0.5)] = 1
        matrix[(matrix<0.5)] = 0
    else:
        matrix = pileup  

    #now compute the pairwise number of 0-0, 0-1, 1-0, 1-1 between all the columns with matrix multiplication
    matrix_1_1 = np.dot(matrix.T,matrix)
    matrix_1_0 = np.dot(matrix.T,1-matrix)
    matrix_0_1 = np.dot(1-matrix.T,matrix)
    matrix_0_0 = np.dot(1-matrix.T,1-matrix)
    pairwise_distance_0 = matrix_0_0/(matrix_0_0+0.5*(matrix_0_1+matrix_1_0))
    pairwise_distance_1 = matrix_1_1/(matrix_1_1+0.5*(matrix_0_1+matrix_1_0))

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

        print("looking at column ", i, "out of ", number_loci, end = '\r')

    #now we can perform the clustering using AgglomerativeClustering
    agglo_cl = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage = 'complete',distance_threshold=0.05)
    agglo_cl.fit(pvalues)

    labels = agglo_cl.labels_
    emblematic_snps = [] #keep one snp per cluster to compute correlations

    #extract each label of size at least two and find the biggest rectangle of 0s in the submatrix
    validated_snps = [False for i in range(number_loci)]
    for label in np.unique(labels):
        idx = np.where(labels == label)[0]
        if len(idx) > 1:
            submatrix = matrix[:,idx]
            #compute the average value of each row and store it in a vector
            row_means = submatrix.mean(axis = 1)
            #find the number of rows of 0s
            nb_rows_0s = len(np.where(row_means <= 0.1)[0])
            nb_rows_with_some_0s = len(np.where(row_means <= 0.9)[0])
            error_rate = min(0.5, nb_rows_with_some_0s/number_reads)

            #now perform the statistical test if the obtained rectangle of 0s is significant
            p_value_of_the_column_cluster = statistical_test(nb_rows_0s, number_reads , len(idx), number_loci, error_rate)

            #validate the snps if the p-value is significant
            if p_value_of_the_column_cluster < 0.001:
                emblematic_snps.append(idx[0])
                for i in idx:
                    validated_snps[i] = True
                    # print("I found a rectangle of ", nb_rows_0s, " rows of 0s in a cluster of ", len(idx), " columns, among in total ", number_loci, " columns and ", number_reads, " rows with an error rate of ", error_rate, ", which gives a p-value of ", p_value_of_the_column_cluster)
                    # for i in range(len(idx)):
                    #     for j in range(number_reads):
                    #         print(int(submatrix[j,i]), end = '')
                    #     print(' : ', idx[i])
                    # print('')
    print("Found ", len(emblematic_snps), " clusters of columns")

    #as a last step, recover the SNPs which correlate well with validated SNPs, using the already computed pvalues
    for i in emblematic_snps:
        for j in range(number_loci):
            if not validated_snps[j]:
                pvalue = stats.chi2_contingency([[matrix_1_1[i,j],matrix_1_0[i,j]],[matrix_0_1[i,j],matrix_0_0[i,j]]])[1] 
                if pvalue < 0.0001:
                    validated_snps[j] = True

    #return the list of validated snps
    return [list_of_sus_pos[i] for i in range(number_loci) if validated_snps[i]]

def output_VCF_of_this_window(bamfile, ref_file, variantpositions, contig_name, start_pos, stop_pos, vcf_file):
    '''
    INPUT: bamfile: the bamfile
           variantpositions: a list of positions where variants were called
           contig_name: the name of the contig
           start_pos: the start position of the window
           stop_pos: the stop position of the window
           vcf_file: the name of the vcf file
    '''

    if len(variantpositions) == 0:
        return

    #open the output file
    out = open(vcf_file,'a')

    #go through the variant positions, grouping positions if they are less than 5bp apart
    index_pos_to_group = 0
    groups_of_variants = []
    while index_pos_to_group < len(variantpositions):
        start_pos_group = variantpositions[index_pos_to_group]
        end_pos_group = start_pos_group+1
        while index_pos_to_group < len(variantpositions)-1 and variantpositions[index_pos_to_group+1] - end_pos_group < 5:
            end_pos_group = variantpositions[index_pos_to_group+1]+1
            index_pos_to_group += 1
        index_pos_to_group += 1
        groups_of_variants.append((start_pos_group,end_pos_group))
    
    index_of_variant_group = 0
    list_of_reads_there = {} #list of the sequence of all the reads in on these positions
    for pileupcolumn in bamfile.pileup(contig = contig_name,start = variant_positions[0],stop = variant_positions[-1]+1,truncate=True):

        if pileupcolumn.reference_pos < groups_of_variants[index_of_variant_group][0]:
            continue

        if pileupcolumn.reference_pos == groups_of_variants[index_of_variant_group][0]:
            #start of the group
            alleles = {}
            list_of_reads_there = {}

        #go through the reads and inventory the alleles
        for pileupread in pileupcolumn.pileups:
            name_of_read = pileupread.alignment.query_name
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
            first_allele = True
            counts = [alleles[reference_sequence]]
            count_first_allele = 0
            for allele in alleles:
                if allele != reference_sequence:
                    if first_allele or (alleles[allele] >= 5 and alleles[allele] >= 0.25*count_first_allele):
                        if not first_allele:
                            out.write(',')
                        else :
                            count_first_allele = alleles[allele]
                            first_allele = False
                        out.write(allele)
                        counts.append(alleles[allele])
            out.write('\t.\t.\tAD\t')
            for i in range(len(counts)):
                if i > 0:
                    out.write(',')
                out.write(str(counts[i]))
            out.write('\n')

            index_of_variant_group += 1

    print("plisduf ")
    for pileupcolumn in bamfile.pileup(contig = contig_name, start=variant_positions[0], stop=variant_positions[-1], truncate=True):
        ...

    #close the output file
    out.close()


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
        '--window', dest='window', required=False, default=5000, type=int,
        help='Size of window to perform read separation (must be at least twice shorter than average read length) [5000]',
    )

    arg = argparser.parse_args()


    return arg.bam, arg.out, arg.window, arg.reference


if __name__ == '__main__':

    print('metaCaller version ', __version__)

    bamfile, out, window, originalAssembly = parse_arguments()
    #mkdir out
    if not os.path.exists(out):
        os.makedirs(out)

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
    o.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+bamfile+'\n')
    o.close()

    ref = ps.FastaFile(originalAssembly)
    #if the reference is not indexed, index it
    if not os.path.exists(originalAssembly+'.fai'):
        print('Indexing the reference')
        os.system('samtools faidx '+originalAssembly)


    #start calling variants
    start = time.time()
    bamfile = ps.AlignmentFile(bamfile,'rb')
    contigs = (bamfile.header.to_dict())['SQ']
    if len(contigs) == 0:
        print('ERROR: No contigs found when parsing the BAM file, check the bam file and the indexation of the bam file')
        sys.exit(1)
    for num in range(0,len(contigs)):
        contig_name = contigs[num]['SN']
        contig_length = contigs[num]['LN']

        print(contig_name, contig_length, ' length')
        window = 5000
        no_snp_threshold = 0.95
        filtered_col_threshold = 0.9
        min_row_quality = 5
        min_col_quality = 3
        list_of_reads = []
        index_of_reads = {}
        haplotypes = []
        list_of_sus_pos = []

        for start_pos in range(0,contig_length,window):

            haplotypes_here = {}

            if start_pos+window <= contig_length:
                # sol_file.write(f'CONTIG\t{contig_name} {start_pos}<->{start_pos+window} \n')

                print(f'Parsing data on contig {contig_name} {start_pos}<->{start_pos+window}')
                #time get_data
                time_get_data = time.time()
                dict_of_sus_pos, list_of_sus_pos = get_data(bamfile, ref, contig_name,start_pos,start_pos+window, no_snp_threshold)
                time_get_data2 = time.time() - time_get_data
            else : 
                # sol_file.write(f'CONTIG\t{contig_name} {start_pos}<->{contig_length} \n')

                print(f'Parsing data on contig  {contig_name} {start_pos}<->{contig_length}')
                dict_of_sus_pos, list_of_sus_pos = get_data(bamfile, ref, contig_name,start_pos,contig_length, no_snp_threshold)


            ###create a matrix from the columns
            df = pd.DataFrame(dict_of_sus_pos) 
            df = df.dropna(axis = 0, thresh = filtered_col_threshold*(len(df.index))) #filter out the reads taht have NaN values in more than 60% of the columns (usually corresponds to the reads that are not spanning the whole window)
            reads = list(df.index)

            ###clustering
            if len(dict_of_sus_pos) > 0 :
                pileup = df.to_numpy()
                time_call_variants = time.time()
                variant_positions = call_variants_in_this_window(pileup, list_of_sus_pos)
                time_call_variants2 = time.time() - time_call_variants

                ###append to the vcffile
                time_output_vcf = time.time()
                output_VCF_of_this_window(bamfile, ref, variant_positions, contig_name, start_pos, start_pos+window, vcf_file)
                time_output_vcf2 = time.time() - time_output_vcf

            print("Times: get_data ", time_get_data2, " call_variants ", time_call_variants2, " output_vcf ", time_output_vcf2)
