#The goal of this code is to count how many variants have been correctly and incorrectly identified in a VCF? given the original genomes

import sys
import argparse
import os
from copy import deepcopy
import pysam
import numpy as np
import pickle
import matplotlib.pyplot as plt


def reverse_complement(seq):
    rc = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return "".join([rc[i] for i in seq[::-1]])

def output_all_potential_kmers_represented_by_the_VCF(vcffile, reffile, outfile, l):
    """
    Generates all potential k-mers (of length 2*l) surrounding each variant in a VCF file and writes them to an output file.

    Args:
        vcffile (str): Path to the input VCF file containing variant information.
        reffile (str): Path to the reference genome file in FASTA format.
        outfile (str): Path to the output file where the k-mers will be written.
        l (int): Length of the k-mers to be generated on each side of the variant (total length will be 2*l).

    The function performs the following steps:
        1. Loads the reference genome sequences into memory.
        2. Parses the VCF file to extract variant positions and their corresponding reference and alternative alleles.
        3. For each variant, generates all possible k-mers (2*l-mers) that can be formed by combining sequences from the reference genome and the variants.
        4. Writes the generated k-mers to the output file in FASTA format.

    The output file will contain sequences in the following format:
        >contig$position$ref_allele$alt_allele$index
        sequence

    Where:
        - contig: The contig/chromosome name.
        - position: The 0-based position of the variant in the contig.
        - ref_allele: The reference allele at the variant position.
        - alt_allele: The alternative allele at the variant position.
        - index: A unique index for each k-mer sequence generated for the variant.
        - sequence: The k-mer sequence.
    """

    #load the reference in memory
    ref_seqs = {}
    with open(reffile) as f:
        name = ""
        for line in f :
            if line[0] == '>':
                name = line.strip().split()[0][1:]
            else:
                ref_seqs[name] = line.strip()

    #go through the vcf and inventoriate the variants
    variant_positions = {} #associates each contig with a list of variant positions
    variants = {} #associates (chrom,pos) to the list of variants there [(ref1, alt1), (ref2, alt2), ...]
    with open(vcffile) as f :
        for line in f:
            if line[0] == '#':
                continue

            else:
                ls = line.split()
                contig = ls[0]
                position = int(ls[1])-1 #-1 because VCF is 1-indexed !!
                ref_seq = ls[3]
                if ref_seq == '.':
                    ref_seq = ""
                alternative_alleles = ls[4].split(',')
                if contig not in variant_positions:
                    variant_positions[contig] = []
                variant_positions[contig].append(position)
                if (contig, position) not in variants :
                    variants[(contig, position)] = []
                for alt in alternative_alleles :
                    if alt == '.' :
                        alt = ""
                    variants[(contig, position)].append((ref_seq, alt))

    #mask regions where there are more than 10 varriants / 2*l+1-mer : there are too many possibilities, the file will be huge (2^10)
    masked_positions = {}
    for contig in variant_positions.keys():
        number_of_variants_around = [0 for i in range(len(ref_seqs[contig]))]
        for pos in variant_positions[contig]:
            for i in range(pos-l, pos+l+1):
                if i >= 0 and i < len(ref_seqs[contig]):
                    number_of_variants_around[i] += 1

        masked_positions[contig] = [(number_of_variants_around[i] > 10) for i in range(len(ref_seqs[contig]))]

    print("All variants indexed")


    #output for each variants all the possible 2*l mers surrounding the variant
    fo = open(outfile, "w")
    index_snp = 0
    for contig in variant_positions.keys():

        #list all the possible l-mers right of any base
        extensions_of_every_base_right = [[] for i in range(len(ref_seqs[contig]))] #associates to each position the list of possible l-mers starting right of this position
        extensions_of_every_base_left = [[] for i in range(len(ref_seqs[contig]))]
        for pos in range(len(ref_seqs[contig])):
            if masked_positions[contig][pos]: #this is to avoid to output too many sequences, even though we lose some information
                continue
            #now list all the extensions right
            print("computing ext for pos ", pos, " on contig ", contig, end='\r')
            extensions_right = [("", pos)] #extensions are a pair (seq, pos_on_ref)
            go_on = True
            while go_on:
                go_on = False
                new_extensions = []
                for past_ext in extensions_right :
                    if len(past_ext[0]) >= l:
                        new_extensions.append(past_ext)
                        continue
                    if past_ext[1] >= len(ref_seqs[contig]):
                        continue

                    go_on = True
                    new_extensions.append((past_ext[0]+ref_seqs[contig][past_ext[1]], past_ext[1]+1))
                    if (contig,past_ext[1]) in variants:
                        for variant in variants[(contig,past_ext[1])] :
                            new_extensions.append((past_ext[0]+variant[1], past_ext[1]+len(variant[0])))

                extensions_right = new_extensions

            for ext in extensions_right:
                extensions_of_every_base_right[pos].append(ext[0])
                if ext[1] < len(ref_seqs[contig]):
                    extensions_of_every_base_left[ext[1]].append(ext[0])

        #now, at each variant positions, combine all the extensions possible left and right to list all the sequences that can represent the 
        nb_variants_we_cannot_output = 0
        for variant_pos in variant_positions[contig] :
            if masked_positions[contig][variant_pos] or (variant_pos>l and masked_positions[contig][variant_pos-l]) or (variant_pos+l < len(ref_seqs[contig]) and masked_positions[contig][variant_pos+l]):
                nb_variants_we_cannot_output += 1
                # print("We cannot output variant at pos ", variant_pos, " because there are too many variants around")
                # fo.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$"+variant[1]+"$-1\n")
                # fo.write("NNNNN\n")
                continue
            for variant in variants[(contig, variant_pos)]:
                index_seq_for_this_snp = 0
                for left_ext in extensions_of_every_base_left[variant_pos]:
                    for right_ext in extensions_of_every_base_right[variant_pos+len(variant[0])]:
                        fo.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$"+variant[1]+"$" + str(index_seq_for_this_snp) +"\n")
                        fo.write(left_ext+variant[1]+right_ext+"\n")
                        index_seq_for_this_snp+=1
                        fo.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$REF\n")
                        fo.write(left_ext+ variant[0] +right_ext+"\n")
                index_snp +=1


    print("We could not output ", nb_variants_we_cannot_output, " variants because there were too many variants around")
    print("We could output ", index_snp, " variants")
    fo.close()
    return nb_variants_we_cannot_output

def map_the_potential_kmers_to_the_reference(potential_file, solution, mapped_file, tmp_dir):

    command = "bwa index " + solution + " -p " + tmp_dir+"/bwaidx"
    res_idx = os.system(command)
    if res_idx != 0:
        print("index failed, ", command)
        sys.exit(1)
    print("ref indexed !")

    command = "bwa fastmap " + tmp_dir + "/bwaidx " + potential_file + " > " + mapped_file
    res_fastmap = os.system(command)
    if res_fastmap != 0:
        print("fastmap failed, ", command)
        sys.exit(1)
    print("seqs aligned !")

    #delete the index
    os.system("rm " + tmp_dir + "/bwaidx*")

def classify_the_variants(mapped_file, output_file):

    variants_confirmed_or_not = {} #associates (contig, pos, ref, alt) to a int: 0 - not tested, 1 - not there, 2 - there
    contig = ""
    pos = 0
    ref = ""
    alt = ""
    length_seq = 0
    variant = ("", 0, "", "")
    untested_variant = False
    with open(mapped_file) as f :
        for line in f:

            ls = line.split()
            if len(ls) == 0 :
                continue
            if ls[0] == "SQ":
                lsls = ls[1].split("$")
                contig = lsls[0]
                pos = int(lsls[1])
                ref = lsls[2]
                alt = lsls[3]
                variant = (contig, pos, ref, alt)
                if lsls[4] == "-1":
                    untested_variant = True
                    variants_confirmed_or_not[variant] = 0
                elif variant not in variants_confirmed_or_not :
                    variants_confirmed_or_not[variant] = 1
                length_seq = int(ls[2].strip())

            elif ls[0] == "EM":
                if int(ls[2]) == length_seq and untested_variant == False:
                    
                    variants_confirmed_or_not[variant] = 2

                untested_variant = False

    nb_good = 0
    nb_bad = 0
    fo = open(output_file, "w")
    for variant in variants_confirmed_or_not.keys():
        print(variant, " : ", variants_confirmed_or_not[variant])
        if variants_confirmed_or_not[variant] == 2 :
            nb_good += 1
            fo.write(variant[0] + "\t" + str(variant[1]) + "\t" + variant[2] + "\t" + variant[3] + "\tCONFIRMED\n")
        elif variants_confirmed_or_not[variant] == 1 :
            nb_bad += 1
            fo.write(variant[0] + "\t" + str(variant[1]) + "\t" + variant[2] + "\t" + variant[3] + "\tNOT CONFIRMED\n")
        elif variants_confirmed_or_not[variant] == 0 :
            fo.write(variant[0] + "\t" + str(variant[1]) + "\t" + variant[2] + "\t" + variant[3] + "\tNOT TESTED\n")

    fo.close()

    print("In total, ", nb_good, " variants there and ", nb_bad, " not there")

def generate_solution_calls_with_minipileup(solution_genomes, ref_file, tmp_dir, solution_calls):

    solution_mapping = tmp_dir + "/solution_mapped.bam"
    command = "minimap2 -a -x asm10 " + ref_file + " " + solution_genomes + " | samtools view -bS - | samtools sort -o " + solution_mapping
    print("command minipileup: ", command)
    res = os.system(command)
    if res != 0:
        print("minimap failed, ", command)
        sys.exit(1)

    command = "samtools index " + solution_mapping
    res = os.system(command)

    #minipileup has a command that looks like this: /home/rfaure/Documents/software/minipileup/minipileup -f one_strain.fa -c refs_mapped.bam
    command = "minipileup -f " + ref_file + " -c " + solution_mapping + " > " + solution_calls
    res = os.system(command)
    if res != 0:
        print("minipileup failed, ", command)
        sys.exit(1)

    print("Generated solution calls")

def count_number_of_missing_variants(solution_calls, vcf_file):

    #inventory the positions of the solution variants on the ref
    variant_pos = []
    with open(solution_calls) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ls = line.split()
            pos = int(ls[1]) - 1  # Convert 1-based to 0-based
            if len(ls[2]) < 3 and len(ls[3]) < 3: #don't take into account big insertions or deletions, they are too hard to catch with just reads
                variant_pos.append((pos,len(ls[3])))

    #convert variant_pos in a bool list
    variant_presence = [False] * (max(pos + length for pos, length in variant_pos) + 1)
    for pos, length in variant_pos:
        for i in range(pos, pos + length):
            variant_presence[i] = True

    # inventory the positions of the VCF variants on the ref
    vcf_variant_pos = []
    vcf_variant_pos_conservative = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ls = line.split()
            pos = int(ls[1]) - 1  # Convert 1-based to 0-based
            vcf_variant_pos.append((pos, len(ls[3])))
            if len(ls[2]) < 3 and len(ls[3]) < 3:
                vcf_variant_pos_conservative.append((pos, len(ls[3])))

    # convert vcf_variant_pos in a bool list
    vcf_variant_presence = [False] * (max(pos + length for pos, length in vcf_variant_pos) + 1)
    for pos, length in vcf_variant_pos:
        for i in range(pos, pos + length):
            vcf_variant_presence[i] = True

    vcf_variant_presence_conservative = [False] * (max(pos + length for pos, length in vcf_variant_pos_conservative) + 1)
    for pos, length in vcf_variant_pos_conservative:
        for i in range(pos, pos + length):
            vcf_variant_presence_conservative[i] = True            

    # count the number of missing variants in the vcf
    missing_variants = []
    variant_there = 0
    for pos in range(len(variant_presence)):
        if variant_presence[pos]:
            if len(vcf_variant_presence) < pos or not vcf_variant_presence[pos] :
                # print("Missing variant at pos ", pos)
                missing_variants.append(pos)
            else :
                variant_there += 1

    #count the number of false positives in vcf_varaints_presence_conservative
    false_positives = []
    for pos in range(len(vcf_variant_presence_conservative)):
        if vcf_variant_presence_conservative[pos]:
            if len(variant_presence) < pos or not variant_presence[pos]:
                # print("False positive variant at pos ", pos)
                false_positives.append(pos)

    print("Precision: ", 1-len(false_positives)/(len(false_positives)+variant_there))

    print("Recall: ", 1-len(missing_variants)/(len(missing_variants)+variant_there))
    print("Number of true variants: ", variant_there)

    return missing_variants

def map_hifi_reads_on_the_kmers_to_check(kmers_to_check_fa, vcf_file, read_file, output_vcf, l):

    # Inventory variants in a dictionary
    variant_dict = {}  # Keys: (chrom, pos, ref), Values: {alt_allele: count}
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ls = line.strip().split()
            chrom = ls[0]
            pos = int(ls[1])
            ref = ls[3]
            alt_alleles = ls[4]
            for alt in alt_alleles:
                key = (chrom, pos, ref)
                if key not in variant_dict:
                    variant_dict[key] = {}
                    variant_dict[key][ref] = 0
                variant_dict[key][alt] = 0

    # Create another variant_dict for non-uniquely mapped reads
    ambiguous_variant_dict = deepcopy(variant_dict)

    # Index the k-mers using bwa
    # command = f"bwa index {kmers_to_check_fa}"
    # res_idx = os.system(command)
    # if res_idx != 0:
    #     print("Indexing k-mers failed:", command)
    #     sys.exit(1)
    # print("K-mers indexed!")

    # Map the reads to the k-mers using bwa fastmap
    mapped_reads_file = "reads_on_kmers_to_check.txt"
    # command = f"bwa fastmap {kmers_to_check_fa} {read_file} -l {2*l+1} | grep -B 1 --no-group-separator ^EM > {mapped_reads_file}"
    # res_map = os.system(command)
    # if res_map != 0:
    #     print("Mapping reads failed:", command)
    #     sys.exit(1)
    # print("Reads mapped to k-mers!")

    # Parse the mapped reads file and update variant_dict with counts
    with open(mapped_reads_file) as f:
        for line in f:
            #EM	11923	11955	1	CP075492.1_3$276$T$C$0:-1: that's the kind of output we have
            if line.startswith("EM") : 
                fields = line.strip().split(":")[0].split("\t")
                # print(fields)
                if fields[4] =="*":
                    continue
                kmer_info = fields[4].split("$")
                chrom = kmer_info[0]
                pos = int(kmer_info[1])
                ref = kmer_info[2]
                alt = kmer_info[3]
                if alt == "REF":
                    alt = ref

                key = (chrom, pos, ref)
                if key in variant_dict and alt in variant_dict[key]:
                    if int(fields[3]) == 1 :
                        variant_dict[key][alt] += 1
                    else:
                        ambiguous_variant_dict[key][alt] += 1

                # Increment and display line counter occasionally
                if 'line_counter' not in locals():
                    line_counter = 0
                line_counter += 1
                if line_counter % 1000 == 0:
                    print(f"Processed {line_counter} lines...", end="\r")
    
    print("finished countging; ")

    # Write the updated VCF with counts
    with open(output_vcf, "w") as fo:
        with open(vcf_file) as f:
            for line in f:
                if line.startswith("#"):  # Write header lines as is
                    fo.write(line)
                    continue
                ls = line.strip().split("\t")
                chrom = ls[0]
                pos = int(ls[1])
                ref = ls[3]
                alt = ls[4]
                key = (chrom, pos, ref)
                ref_count = 0
                alt_count = 0
                alt_ref_count=0
                alt_alt_count=0
                if key in variant_dict :
                    if ref in variant_dict[key]:
                        ref_count = variant_dict[key][ref]
                    if alt in variant_dict[key]:
                        alt_count = variant_dict[key][alt]
                    alt_ref_count = ambiguous_variant_dict[key][ref] if ref in ambiguous_variant_dict[key] else 0
                    alt_alt_count = ambiguous_variant_dict[key][alt] if alt in ambiguous_variant_dict[key] else 0
                fo.write("\t".join(ls) + f"\tCOUNTS:{ref_count},{alt_count}\tACOUNTS:{alt_ref_count},{alt_alt_count}\n")

def stats_on_variants(vcf_file_with_counts):

    all_variant_callers = ["bcfcalls_normalized.vcf","clair_normalized.vcf","longshotcalls_normalized.vcf","nanocaller_normalized.vcf","snpsnoop_normalized.vcf"]
    good_bad_unsure_missed_variants = {i:[0,0,0,0] for i in all_variant_callers}
    this_assembler_in_another_variant= {i:-1 for i in all_variant_callers}

    with open(vcf_file_with_counts) as f:

        for line in f :
            if line[0] == "#":
                continue

            ls = line.strip().split("\t")
            chrom = ls[0]
            position = int(ls[1])
            ref = ls[3]
            callers = ls[7].split(',')
            counts_sure = [int(i) for i in ls[8].lstrip("COUNTS:").split(',')]
            counts_unsure = [int(i) for i in ls[9].lstrip("ACOUNTS:").split(',')]
            
            variant = "good"
            if counts_sure[1] == 0 and counts_unsure[1] == 0 :
                variant = "bad"
            elif counts_sure[1] == 0 and counts_unsure[1] > 0:
                variant = "unsure"

            if variant == "bad" and "snpsnoop_normalized.vcf" in callers and counts_sure[0]+counts_sure[1]+counts_unsure[0]+counts_unsure[1] > 10 and chrom=="contig_10":
                print(f"Bad in snpsnoop: {chrom}, {position}, {ref}, ", line)

            for caller in all_variant_callers:
                if caller in callers:
                    if variant == "good":
                        good_bad_unsure_missed_variants[caller][0] += 1
                    elif variant == "unsure":
                        good_bad_unsure_missed_variants[caller][3] += 1
                    elif variant == "bad":
                        good_bad_unsure_missed_variants[caller][1] += 1
                    this_assembler_in_another_variant[caller] = position - 1 + len(ref)
                elif position > this_assembler_in_another_variant[caller] :
                    if variant == "good":
                        good_bad_unsure_missed_variants[caller][2] += 1 #missed variant
                    elif variant == "unsure":
                        good_bad_unsure_missed_variants[caller][3] += 1
    
    for caller in all_variant_callers:
        print(f"Variant caller: {caller}")
        print(f"Good variants: {good_bad_unsure_missed_variants[caller][0]}")
        print(f"Bad variants: {good_bad_unsure_missed_variants[caller][1]}")
        print(f"Missed variants: {good_bad_unsure_missed_variants[caller][2]}")
        print(f"Unsure variants: {good_bad_unsure_missed_variants[caller][3]}")
        print("-" * 40)

#input: vcf file containing all variants, bam file of hifi mapped on ref
def compare_calls_with_HiFi(all_variants, mapped_hifi):

    # # Index all positions (contig_name, position) from the input VCF
    # variants = {}
    # number_of_non_indexed_variants = 0
    # with open(all_variants) as vcf_file:
    #     for line in vcf_file:
    #         if line.startswith("#"):
    #             continue
    #         fields = line.strip().split("\t")
    #         if len(fields[3]) != len(fields[4]) or ',' in fields[4]: #insertion/deletion, we'll look at that later
    #             number_of_non_indexed_variants += 1
    #             continue
    #         fields[4] = fields[4].replace('.', '*')
    #         callers = fields[-2].split(',')
    #         depths = [int(i) if i != '*' else 0 for i in fields[-1].split(',')]
    #         contig_name = fields[0]
    #         position = int(fields[1]) - 1  # Convert 1-based to 0-based
    #         for i in range(len(fields[3])) :
    #             if fields[3][i] == fields[4][i]:
    #                 continue
    #             if (contig_name, position+i) not in variants :
    #                 variants[(contig_name, position+i)] = {}
    #                 variants[(contig_name, position+i)][("ref",0)] = fields[3][i]
    #             for c, caller in enumerate(callers) :
    #                 variants[(contig_name, position+i)][(caller, depths[c])] = fields[4][i]
    # print(f"Indexed {len(variants)} positions from the input VCF. Not indexed ", number_of_non_indexed_variants)

    # # Create a BED file describing all positions we need
    # bed_file_path = os.path.join(tmp_dir, "positions_to_check.bed")
    # with open(bed_file_path, "w") as bed_file:
    #     for (contig, pos), variant_info in variants.items():
    #         bed_file.write(f"{contig}\t{pos}\t{pos + 1}\n")
    # print(f"BED file created at {bed_file_path}")

    # # Step 1: Generate the pileup for all variant positions
    # pileup_file = os.path.join(tmp_dir, "variants.pileup")
    # variant_positions = ",".join([f"{contig}:{pos + 1}" for contig, pos in variants.keys()])
    # command = f"samtools mpileup -l {bed_file_path} {mapped_hifi} > {pileup_file}"
    # res = os.system(command)
    # if res != 0:
    #     print("Pileup generation failed:", command)
    #     sys.exit(1)
    # print("Pileup generation completed.")

    # # Step 2: Parse the pileup and count base occurrences for each variant position
    # base_counts = {}
    # with open(pileup_file) as pileup_in:
    #     for line in pileup_in:
    #         fields = line.strip().split("\t")
    #         if len(fields) >= 5:
    #             contig, pos = fields[0], int(fields[1]) - 1  # Convert 1-based to 0-based
    #             bases = fields[4]
    #             counts = {base: bases.upper().count(base) for base in "ATCG*"}
    #             base_counts[(contig, pos)] = counts
    #         else:
    #             contig, pos = fields[0], int(fields[1]) - 1
    #             base_counts[(contig, pos)] = {"A": 0, "T": 0, "C": 0, "G": 0, "*": 0}
    # print("Base counting completed.")

    # # Save base_counts and variants to pickle files to avoid recalculating
    # base_counts_pickle_file = os.path.join(tmp_dir, "base_counts.pkl")
    # variants_pickle_file = os.path.join(tmp_dir, "variants.pkl")

    # if not os.path.exists(base_counts_pickle_file) or not os.path.exists(variants_pickle_file):
    #     # Save base_counts
    #     with open(base_counts_pickle_file, "wb") as pkl_out:
    #         pickle.dump(base_counts, pkl_out)
    #     print(f"Base counts saved to {base_counts_pickle_file}")

    #     # Save variants
    #     with open(variants_pickle_file, "wb") as pkl_out:
    #         pickle.dump(variants, pkl_out)
    #     print(f"Variants saved to {variants_pickle_file}")
    # else:
    # Load base_counts

    base_counts_pickle_file = os.path.join(tmp_dir, "base_counts.pkl")
    variants_pickle_file = os.path.join(tmp_dir, "variants.pkl")
    # Load variants
    with open(variants_pickle_file, "rb") as pkl_in:
        variants = pickle.load(pkl_in)
    print(f"Variants loaded from {variants_pickle_file}")
    with open(base_counts_pickle_file, "rb") as pkl_in:
        base_counts = pickle.load(pkl_in)
    print(f"Base counts loaded from {base_counts_pickle_file}")

    # Initialize statistics for each caller
    caller_stats = {}
    total_number_of_true_variants = 0
    hist_cov_correctly_called_snpsnoop = [0 for i in range(100)]
    hist_cov_false_positive_snpsnoop = [0 for i in range(100)]
    hist_cov_false_negative_snpsnoop = [0 for i in range(100)]
    for variant_key, variant_info in variants.items():
        contig, pos = variant_key
        ref_base = variant_info[("ref", 0)]
        found_variant = False
        for caller_and_depth, alt_base in variant_info.items():
            caller, depth = caller_and_depth
            if caller == "ref":
                continue
            if variant_info[caller_and_depth] == "*":
                continue
            if alt_base == "*" or caller not in caller_stats:
                caller_stats[caller] = {"recall": 0, "false_positives": 0, "total": 0}
            caller_stats[caller]["total"] += 1

            # Check if the base is seen in the pileup
            if variant_key in base_counts:

                counts = base_counts[variant_key]
                if counts[alt_base] > 0:
                    caller_stats[caller]["recall"] += 1
                    found_variant = True 
                    if not any(caller == "snpsnoop_normalized.vcf" for caller, _ in variant_info.keys()):
                        if depth > 30:
                            print("Missing snp : ", variant_key, " ", variant_info, ", DP in other caller: ", depth)
                        hist_cov_correctly_called_snpsnoop[min(99,np.sum(list(counts.values())))] += 1
                    elif caller == "bcfcalls_normalized.vcf":
                        hist_cov_false_negative_snpsnoop[min(99,depth)] += 1

                else:
                    caller_stats[caller]["false_positives"] += 1
                    hist_cov_false_positive_snpsnoop[min(99,depth)] += 1
                    # if caller == "snpsnoop_normalized.vcf" :#and "bcfcalls_normalized" not in variant_info.keys():
                    #     if contig == "contig_3747" :
                    #         print("false positive snp ", variant_key, " ", variant_info)
                    # if caller == "snpsnoop_normalized.vcf" and not any(c == "bcfcalls_normalized" for c, _ in variant_info.keys()):
                caller_stats[caller]["false_positives"] += 1

        if found_variant:
            total_number_of_true_variants += 1

    # Print statistics for each caller
    for caller, stats in caller_stats.items():
        recall = stats["recall"] / total_number_of_true_variants
        false_positive_rate = stats["false_positives"] / (stats["recall"]+stats["false_positives"] )
        print(f"Caller: {caller}")
        print(f"  Recall: {recall:.2%}")
        print(f"  False Positive Rate: {false_positive_rate:.2%}")
        print(f"  Total Variants: {stats['total']}")
        print(f"  True Positives: {stats['recall']}")
        print(f"  False Positives: {stats['false_positives']}")
        print("-" * 40)

    # Plot histograms for SNP-Snoop
    x = np.arange(100)

    plt.figure(figsize=(12, 6))
    # plt.bar(x - 0.3, hist_cov_correctly_called_snpsnoop, width=0.3, label="Correctly Called", color="green")
    plt.bar(x, hist_cov_false_positive_snpsnoop, width=0.3, label="False Positives", color="red")
    plt.bar(x + 0.3, hist_cov_false_negative_snpsnoop, width=0.3, label="False Negatives", color="blue")

    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.title("Coverage Histogram for SNP-Snoop Variants")
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate VCF with solution")
    parser.add_argument("-v", "--vcf", required=True, help="Path to the input VCF file containing variant information")
    parser.add_argument("-r", "--ref", required=True, help="Path to the reference genome file in FASTA format")
    parser.add_argument("-f", "--reads", required=False, help="Path to the HiFi reads")
    parser.add_argument("-s", "--solution", required=False, help="Path to the solution file")
    parser.add_argument("-l", "--length", type=int, default=25, help="Length of the k-mers to be generated on each side of the variant (default: 15)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output VCF file with the counts")

    args = parser.parse_args()

    tmp_dir = os.path.join(os.path.dirname(args.output), "tmp")
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass

    reference_genome_available = False
    if not args.solution and not args.reads:
        print("Error: Either a solution genome file (--solution) or HiFi reads (--reads) must be provided.")
        sys.exit(1)
    elif args.solution:
        reference_genome_available = True

    #to normalize the variants: first use normalize_variants.py, then vt


    if reference_genome_available :
         # nb_of_variants_we_cannot_output = output_all_potential_kmers_represented_by_the_VCF(args.vcf, args.ref, "kmers_to_check.fa", args.length)
          # print("Outputted all variants as k-mers")

        map_the_potential_kmers_to_the_reference(os.path.join(tmp_dir, "kmers_to_check.fa"), args.solution, os.path.join(tmp_dir, "mapped.txt"), tmp_dir)
        classify_the_variants(os.path.join(tmp_dir, "mapped.txt"), args.output)
        print("Did not test ", nb_of_variants_we_cannot_output, " variants because there were too many variants around")

        # solution_calls = os.path.join(tmp_dir, "solution.vcf")
        # generate_solution_calls_with_minipileup(args.solution, args.ref, tmp_dir, solution_calls)
        # missing_variants = count_number_of_missing_variants(solution_calls, args.vcf)
    else:
        mapped_reads_bam = "hifi_on_ref.bam"

        # # Index the BAM file
        # command = f"samtools index {mapped_reads_bam}"
        # res = os.system(command)
        # if res != 0:
        #     print("Indexing BAM file failed:", command)
        #     sys.exit(1)

        # print("HiFi reads successfully mapped and indexed.")

        # map_hifi_reads_on_the_kmers_to_check(
        #     kmers_to_check_fa="kmers_to_check.fa",
        #     vcf_file=args.vcf,
        #     read_file=args.reads,  # Replace with the actual path to HiFi reads
        #     output_vcf=args.output,
        #     l = args.length
        # )
        # print("mapped the reads")

        # stats_on_variants(args.output)

        # # Map HiFi reads to the reference genome using minimap2
        # command = f"minimap2 -a -x map-hifi {args.ref} {args.reads} | samtools view -bS - | samtools sort -o {mapped_reads_bam}"
        # print("Mapping HiFi reads to the reference genome:", command)
        # res = os.system(command)
        # if res != 0:
        #     print("Mapping HiFi reads failed:", command)
        #     sys.exit(1)

        compare_calls_with_HiFi(args.vcf, mapped_reads_bam)






