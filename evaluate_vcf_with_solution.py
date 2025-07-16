#The goal of this code is to count how many variants have been correctly and incorrectly identified in a VCF? given the original genomes

import sys
import argparse
import os
from copy import deepcopy
# import pysam
import numpy as np
import pickle
import gzip
import shutil
import random
import pysam
import glob
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

    # Load variants called only by bcftools and only by snpsnoop if they exist
    only_bcftools_variants = []
    only_snpsnoop_variants = []
    only_bcftools_path = os.path.join(os.path.dirname(outfile), "tmp", "only_bcftools_variants.pkl")
    only_snpsnoop_path = os.path.join(os.path.dirname(outfile), "tmp", "only_snpsnoop_variants.pkl")
    if os.path.exists(only_bcftools_path):
        with open(only_bcftools_path, "rb") as f:
            only_bcftools_variants = pickle.load(f)
    if os.path.exists(only_snpsnoop_path):
        with open(only_snpsnoop_path, "rb") as f:
            only_snpsnoop_variants = pickle.load(f)

    outfile_snpsnoop = "kmers_to_check_only_snpsnoop.txt"
    outfile_bcftools = "kmers_to_check_only_bcftools.txt"
    # Create output files for variants only called by snpsnoop or bcftools
    outfile_snpsnoop = os.path.join(os.path.dirname(outfile), "tmp", "kmers_variants_only_snpsnoop.txt")
    outfile_bcftools = os.path.join(os.path.dirname(outfile), "tmp", "kmers_variants_only_bcftools.txt")

    #load the reference in memory
    ref_seqs = {}
    with open(reffile) as f:
        name = ""
        for line in f :
            if line[0] == '>':
                name = line.strip().split()[0][1:]
            else:
                ref_seqs[name] = line.strip()

    print("going through the vcf")
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

        masked_positions[contig] = [(number_of_variants_around[i] > 5) for i in range(len(ref_seqs[contig]))]

    print("All variants indexed")


    #output for each variants all the possible 2*l mers surrounding the variant
    fo = open(outfile, "w")
    fo1 = open(outfile_snpsnoop, "w")
    fo2 = open(outfile_bcftools, "w")
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
                        fo.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$REF\n")
                        fo.write(left_ext+ variant[0] +right_ext+"\n")

                        # Output if the variant is only in bcftools or only in snpsnoop
                        if (contig, variant_pos, variant[0], variant[1]) in only_bcftools_variants:
                            fo2.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$"+variant[1]+str(index_seq_for_this_snp) +"\n")
                            fo2.write(left_ext+variant[1]+right_ext+"\n")
                        if (contig, variant_pos, variant[0], variant[1]) in only_snpsnoop_variants:
                            fo1.write(">"+contig+"$"+str(variant_pos+1)+"$"+variant[0]+"$"+variant[1] + str(index_seq_for_this_snp) +"\n")
                            fo1.write(left_ext+variant[1]+right_ext+"\n")

                        index_seq_for_this_snp+=1

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

def generate_solution_calls_with_minipileup(solution_genomes, ref_file, tmp_dir, solution_calls_few_snps, solution_calls_many_snps):

    solution_mapping = tmp_dir + "/solution_mapped_few.bam"

    # Cut the reference genomes into overlapping 300bp chunks and write to a temporary FASTA
    chunk_size = 300
    overlap = 150  # 50% overlap
    chunked_solution_fasta = os.path.join(tmp_dir, "solution_chunks.fa")
    with open(solution_genomes) as f_in, open(chunked_solution_fasta, "w") as f_out:
        seq_name = ""
        seq = []
        for line in f_in:
            if line.startswith(">"):
                if seq_name and seq:
                    seq_str = "".join(seq)
                    for i in range(0, len(seq_str), chunk_size - overlap):
                        chunk = seq_str[i:i+chunk_size]
                        if len(chunk) < 50:
                            continue
                        f_out.write(f">{seq_name}_chunk_{i}\n{chunk}\n")
                seq_name = line.strip()[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if seq_name and seq:
            seq_str = "".join(seq)
            for i in range(0, len(seq_str), chunk_size - overlap):
                chunk = seq_str[i:i+chunk_size]
                if len(chunk) < 50:
                    continue
                f_out.write(f">{seq_name}_chunk_{i}\n{chunk}\n")
    command = (
        "minimap2 -a -x map-ont -E 10,9 "
        + ref_file + " " + chunked_solution_fasta +
        " 2> /tmp/trash.txt | samtools view -bS -F 256 -F 2048 - | samtools sort -o " + solution_mapping
    )
    #print("command minipileup: ", command)
    res = os.system(command)
    if res != 0:
        print("minimap failed, ", command)
        sys.exit(1)

    command = "samtools index " + solution_mapping
    res = os.system(command)

    #minipileup has a command that looks like this: /home/rfaure/Documents/software/minipileup/minipileup -f one_strain.fa -c refs_mapped.bam
    command = "minipileup -f " + ref_file + " -c " + solution_mapping + " > " + solution_calls_few_snps + " 2> /tmp/trash.txt"
    res = os.system(command)
    if res != 0:
        print("minipileup failed, ", command)
        sys.exit(1)

    solution_mapping = tmp_dir + "/solution_mapped_many.bam"
    # Generate solution calls with minipileup using both the full (uncut) solution genomes and cut genomes and without penalizing gaps too much
    # Concatenate the full solution genomes and the chunked solution fasta into a single file using cat
    concatenated_solution_fasta = os.path.join(tmp_dir, "solution_full_and_chunks.fa")
    os.system(f"cat {solution_genomes} {chunked_solution_fasta} > {concatenated_solution_fasta}")
    command = (
        "minimap2 -a -x map-ont " + ref_file + " " + concatenated_solution_fasta +
        " 2> /tmp/trash.txt | samtools view -bS - | samtools sort -o " + solution_mapping
    )
    #print("command minipileup (many snps): ", command)
    res = os.system(command)
    if res != 0:
        print("minimap failed, ", command)
        sys.exit(1)

    command = "samtools index " + solution_mapping
    res = os.system(command)

    command = "minipileup -f " + ref_file + " -c " + solution_mapping + " > " + solution_calls_many_snps + " 2> /tmp/trash.txt"
    res = os.system(command)
    if res != 0:
        print("minipileup failed, ", command)
        sys.exit(1)

    print("Generated solution calls")

def count_number_of_missing_variants(solution_calls_few_snps, solution_calls_many_snps, vcf_file):

    # Read VCF header to get contig lengths
    contig_lengths = {}
    with open(solution_calls_few_snps) as f:
        for line in f:
            if line.startswith("##contig="):
                # Example: ##contig=<ID=contig_12,length=15323>
                try:
                    id_part = line.split("ID=")[1]
                    contig_id = id_part.split(",")[0]
                    length_part = line.split("length=")[1]
                    contig_length = int(length_part.split(">")[0].replace(">", "").replace("]", "").replace("'", "").replace('"', '').strip())
                    contig_lengths[contig_id] = contig_length
                except Exception:
                    continue

    contigs_with_some_variants=set()

    #inventory the positions of the solution variants on the ref
    variant_pos_few = []
    with open(solution_calls_few_snps) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ls = line.split()
            pos = int(ls[1]) - 1  # Convert 1-based to 0-based
            contigs_with_some_variants.add(ls[0])
            if any([len(i)==len(ls[3]) and i != "." for i in ls[4].split(',')]) :  #don't take into account big insertions or deletions, they are too hard to catch with just reads
                variant_pos_few.append((ls[0], pos,len(ls[3])))

    #convert variant_pos in one bool list per reference
    variant_presence_few = {c : [False] * contig_lengths[c] for c in contig_lengths.keys()}
    for contig, pos, length in variant_pos_few:
        for i in range(pos, pos + length):
            variant_presence_few[contig][i] = True

    # inventory the positions of the solution variants (many snps) on the ref
    variant_pos_many = []
    with open(solution_calls_many_snps) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ls = line.split()
            pos = int(ls[1]) - 1  # Convert 1-based to 0-based
            if any([len(i) == len(ls[3]) and i != "." for i in ls[4].split(',')]):
                variant_pos_many.append((ls[0], pos, len(ls[3])))

    # convert variant_pos_many into one bool list per reference
    variant_presence_many = {c: [False] * contig_lengths[c] for c in contig_lengths.keys()}
    for contig, pos, length in variant_pos_many:
        for i in range(pos, pos + length):
            variant_presence_many[contig][i] = True

    # Accept a list of VCF files as input
    if isinstance(vcf_file, str):
        vcf_files = [vcf_file]
    else:
        vcf_files = vcf_file

    # Prepare output statistics for each VCF
    stats_per_vcf = []

    for vcf_path in vcf_files:
        # inventory the positions of the VCF variants on the ref
        vcf_variant_pos = []
        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                ls = line.split()
                pos = int(ls[1]) - 1  # Convert 1-based to 0-based
                if any([len(i)==len(ls[3]) and i != "." for i in ls[4].split(',')]) : 
                    vcf_variant_pos.append((ls[0], pos, len(ls[3])))

        # convert vcf_variant_pos into one bool list per reference
        vcf_variant_presence = {c: [False] * contig_lengths[c] for c in contig_lengths.keys()}
        for contig, pos, length in vcf_variant_pos:
            for i in range(pos, pos + length):
                vcf_variant_presence[contig][i] = True

        # Compare vcf_variant_presence and variant_presence to evaluate precision and recall
        missing_variants = []
        false_positives = []
        total_true_variants = 0
        total_variant_called = 0

        false_positive_vcf = f"{os.path.splitext(os.path.basename(vcf_path))[0]}_false_positive.vcf"
        missing_variant_vcf = f"{os.path.splitext(os.path.basename(vcf_path))[0]}_missing_variant.vcf"

        # # Clear previous output files if they exist
        open(false_positive_vcf, "w").close()
        # open(missing_variant_vcf, "w").close()

        for contig in contig_lengths:
            if contig not in contigs_with_some_variants:
                continue
            for i in range(contig_lengths[contig]):
                if variant_presence_few[contig][i]:
                    total_true_variants += 1
                    if not vcf_variant_presence[contig][i]:
                        missing_variants.append((contig, i))  
                        # with open(missing_variant_vcf, "a") as mvf:
                        #     mvf.write(f"{contig}\t{i+1}\t.\tN\tN\t.\t.\tMISSING\n")

                if vcf_variant_presence[contig][i]:
                    total_variant_called += 1
                    if not variant_presence_many[contig][i]:
                        false_positives.append((contig, i))
                        with open(false_positive_vcf, "a") as fpvf:
                            fpvf.write(f"{contig}\t{i+1}\t.\tN\tN\t.\t.\tFALSE_POSITIVE\n")

        precision = 1 - len(false_positives) / total_variant_called if total_variant_called > 0 else 0
        recall = 1 - len(missing_variants) / total_true_variants if total_true_variants > 0 else 0

        print(f"VCF: {vcf_path}")
        print("Precision: ", precision)
        print("Recall: ", recall)
        print("Number of true variants: ", total_true_variants, " and called ", total_variant_called)
        print("-" * 40)

        stats_per_vcf.append({
            "vcf": vcf_path,
            "precision": precision,
            "recall": recall,
            "true_variants": total_true_variants,
            "called_variants": total_variant_called,
            "missing_variants": missing_variants,
            "false_positives": false_positives,
            "missing_variant_vcf": missing_variant_vcf,
            "false_positive_vcf": false_positive_vcf
        })

    return stats_per_vcf


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
def compare_calls_with_HiFi(all_variants, mapped_hifi, false_positive_SNPs_file=None, solution_known = False):

    set_of_contig = {
        "contig_10182",
        "contig_10184",
        "contig_10482",
        "contig_12844",
        "contig_3090",
        "contig_5557",
        "contig_5905",
        "contig_8159",
        "contig_8582",
        "contig_8727"
    }
    set_of_contig={}
    tmp_dir = "tmp"

    # # Index all positions (contig_name, position) from the input VCF
    # variants = {}
    # number_of_non_indexed_variants = 0
    # with open(all_variants) as vcf_file:
    #     for line in vcf_file:
    #         if line.startswith("#"):
    #             continue
    #         fields = line.strip().split("\t")
    #         contig_name = fields[0]
    #         if len(set_of_contig) > 0 and contig_name not in set_of_contig:
    #             continue
    #         if len(fields[3]) != len(fields[4]) or ',' in fields[4]: #insertion/deletion, we'll look at that later
    #             number_of_non_indexed_variants += 1
    #             continue
    #         fields[4] = fields[4].replace('.', '*')
    #         callers = fields[-2].split(',')
    #         depths = [int(i) if i != '*' else 0 for i in fields[-1].split(',')]
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

    # variants_pickle_file = os.path.join(tmp_dir, "variants.pkl")
    # with open(variants_pickle_file, "wb") as pkl_out:
    #     pickle.dump(variants, pkl_out)
    # print(f"Variants saved to {variants_pickle_file}")

    # if not solution_known :
    #     # Create a BED file describing all positions we need
    #     bed_file_path = os.path.join(tmp_dir, "positions_to_check.bed")
    #     with open(bed_file_path, "w") as bed_file:
    #         for (contig, pos), variant_info in variants.items():
    #             bed_file.write(f"{contig}\t{pos}\t{pos + 1}\n")
    #     print(f"BED file created at {bed_file_path}")

    #     # Step 1: Generate the pileup for all variant positions
    #     pileup_file = os.path.join(tmp_dir, "variants.pileup")
    #     command = f"samtools mpileup -B -A -l {bed_file_path} {mapped_hifi} > {pileup_file}"
    #     res = os.system(command)
    #     if res != 0:
    #         print("Pileup generation failed:", command)
    #         sys.exit(1)
    #     print("Pileup generation completed.")
    # else:
    #     pileup_file = os.path.join(tmp_dir, "variants.pileup")
    #     command = f"samtools mpileup {mapped_hifi} > {pileup_file}"
    #     res = os.system(command)
    #     if res != 0:
    #         print("Pileup generation failed:", command)
    #         sys.exit(1)
    #     print("Pileup generation completed.")

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

    # # # Save base_counts
    # with open(base_counts_pickle_file, "wb") as pkl_out:
    #     pickle.dump(base_counts, pkl_out)
    # print(f"Base counts saved to {base_counts_pickle_file}")

    base_counts_pickle_file = os.path.join(tmp_dir, "base_counts.pkl")
    variants_pickle_file = os.path.join(tmp_dir, "variants.pkl")
    # Load variants
    with open(variants_pickle_file, "rb") as pkl_in:
        variants = pickle.load(pkl_in)
    print(f"Variants loaded from {variants_pickle_file}")
    with open(base_counts_pickle_file, "rb") as pkl_in:
        base_counts = pickle.load(pkl_in)
    print(f"Base counts loaded from {base_counts_pickle_file}")

    already_outputted_false_positives = set()

    # Initialize statistics for each caller
    caller_stats = {}
    total_number_of_true_variants = 0
    hist_cov_correctly_called_snpsnoop = [0 for i in range(100)]
    hist_cov_false_positive_snpsnoop = [0 for i in range(100)]
    hist_cov_false_negative_snpsnoop = [0 for i in range(100)]
    variants_called_by_both_tools = 0
    variant_called_by_only_snspnoop = set()
    variant_called_by_only_bcftools = set()
    for variant_key, variant_info in variants.items():
        contig, pos = variant_key
        if variant_key not in base_counts:
            # print("variant ker ", variant_key, variant_info)
            continue
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
            if alt_base == "N":
                continue
            caller_stats[caller]["total"] += 1

            # Check if the base is seen in the pileup
            if variant_key in base_counts:

                counts = base_counts[variant_key]
                if counts[alt_base] > 0:
                    caller_stats[caller]["recall"] += 1
                    found_variant = True 
                    # if not any(caller == "snpsnoop_normalized.vcf" or caller == "trash.vcf" for caller, _ in variant_info.keys()):
                    #     if depth > 30 or True:
                    #         print("Missing snp : ", variant_key, " ", variant_info, ", DP in other caller: ", depth)
                        # hist_cov_correctly_called_snpsnoop[min(99,np.sum(list(counts.values())))] += 1
                    # elif caller == "bcfcalls_normalized.vcf":
                    #     hist_cov_false_negative_snpsnoop[min(99,depth)] += 1

                else:
                    caller_stats[caller]["false_positives"] += 1
                    hist_cov_false_positive_snpsnoop[min(99,depth)] += 1
                    if caller == "metacaller_normalized.vcf" and "deepvariant_normalized.vcf" not in [i[0] for i in variant_info.keys()]:
                        print("not in deepvariant ", variant_key, " ", variant_info, " ", counts)
                        variant_called_by_only_snspnoop.add((variant_key, alt_base))
                    elif caller == "deepvariant_normalized.vcf" and "metacaller_normalized.vcf" not in [i[0] for i in variant_info.keys()]:
                        print("not in metacaller ", variant_key, " ", variant_info, " ", counts)
                        variant_called_by_only_bcftools.add((variant_key, alt_base))
                    elif caller == "deepvariant_normalized.vcf":
                        variants_called_by_both_tools+=1
                    
                    # Output false positive SNPs to file if not already outputted
                    false_positive_SNPs_file = "false_positives.vcf"
                    false_positive_key = (contig, pos, ref_base, alt_base, caller)
                    if false_positive_key not in already_outputted_false_positives:
                        callers_list = [c for c, _ in variant_info.keys() if c != "ref"]
                        callers_str = ",".join(callers_list)
                        with open(false_positive_SNPs_file, "a") as fp_out:
                            fp_out.write(f"{contig}\t{pos+1}\t.\t{ref_base}\t{alt_base}\t.\t.\tCALLER={callers_str}\n")
                        already_outputted_false_positives.add(false_positive_key)

        if found_variant:
            total_number_of_true_variants += 1

    if solution_known:
        # Identify true variants from HiFi pileup: record positions with no consensus
        true_variants = []
        for (contig, pos), counts in base_counts.items():
            if pos > 100000 : 
                continue
            total = sum(counts[base] for base in "ATCG")
            if total == 0:
                continue
            max_base = max("ATCG", key=lambda b: counts[b])
            if counts[max_base] < total:
                # No consensus: at least one other base is present
                true_variants.append((contig, pos))
        print(f"Number of true variant positions (no consensus): {len(true_variants)}")
        total_number_of_true_variants = len(true_variants)

        # Compute the number of false negatives for all callers and save in stats
        for caller in caller_stats.keys():
            false_negatives = 0
            for (contig, pos) in true_variants:
                found = False
                for (c, _), alt_base in variants.get((contig, pos), {}).items():
                    if c == caller and alt_base != "ref" and alt_base != "*":
                        found = True
                        break
                if not found:
                    false_negatives += 1
                    print(f"False negative for {caller} at {contig}:{pos}")
            caller_stats[caller]["recall"] = len(true_variants)-false_negatives
            
            print(f"Number of false negatives for {caller}: {false_negatives}")

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
    print("Variagnd  dsfkl : ", variants_called_by_both_tools, " ", len(variant_called_by_only_bcftools), " ", len(variant_called_by_only_snspnoop))

    # Output variants called only by bcftools and only by snpsnoop as pickles
    with open(os.path.join(tmp_dir, "only_bcftools_variants.pkl"), "wb") as f:
        pickle.dump(list(variant_called_by_only_bcftools), f)
    with open(os.path.join(tmp_dir, "only_snpsnoop_variants.pkl"), "wb") as f:
        pickle.dump(list(variant_called_by_only_snspnoop), f)

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

def create_Logan_queries_to_check_false_positives(false_positive_SNPs_file, mapped_nanopore):

    # Sample randomly from the false positive SNPs: 10 SNPs where snpsnoop_normalized.vcf and bcfcalls_normalized.vcf are in CALLER
    snps = []
    with open(false_positive_SNPs_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            info = fields[7] if len(fields) > 7 else ""
            if "CALLER=" in info or True:
                #callers = info.split("CALLER=")[1].split(",")
                #if ("snpsnoop_normalized.vcf" in callers and "bcfcalls_normalized.vcf" in callers) or True:
                snps.append(line.strip())
    sampled = random.sample(snps, min(10, len(snps)))
    # Instead of writing to a file, collect (contig, position, ref, alt) in a list
    sampled_variants = []
    for snp in sampled:
        fields = snp.split("\t")
        contig = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        sampled_variants.append((contig, pos, ref, alt))
    
    # Open the BAM file with HiFi reads mapped to the reference
    bamfile = pysam.AlignmentFile(mapped_nanopore, "rb")

    for contig, pos, ref, alt in sampled_variants:
        # 0-based position for pysam
        snp_pos = pos - 1

        # Fetch reads overlapping the SNP position
        reads_supporting_alt = []
        # Use pileup to iterate over all reads at the SNP position with specified flags
        itercol = bamfile.pileup(
            contig=contig,
            start=snp_pos,
            stop=snp_pos + 1,
            truncate=True,
            min_base_quality=0,
            multiple_iterators=True,
            stepper="all",
            flag_filter=4 | 256 | 512 | 1024
        )
        for pileupcolumn in itercol:
            allreads = pileupcolumn.get_query_names()
            allbases_raw = [i.upper() for i in pileupcolumn.get_query_sequences(add_indels=True)]
            # Determine the two most frequent bases at this position
            base_counts = {}
            for base in allbases_raw:
                if base not in base_counts:
                    base_counts[base] = 0
                base_counts[base] += 1
            # Sort bases by frequency
            sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_bases) >= 2:
                ref, alt = sorted_bases[0][0], sorted_bases[1][0]
            elif len(sorted_bases) == 1:
                ref, alt = sorted_bases[0][0], None
            else:
                ref, alt = None, None
            reads_supporting_alt = [allreads[i] for i in range(len(allreads)) if allbases_raw[i] == alt]


        # Build consensus from 50bp left and right of the SNP (total 100bp window)
        window_size = 50
        consensus_window = []
        for offset in range(-window_size, window_size):
            base_pos = snp_pos + offset
            if base_pos < 0:
                consensus_window.append("N")
                continue
            # Use pileup to get base calls at this position, but only from reads_supporting_alt
            base_counts = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
            for pileupcolumn in bamfile.pileup(
                contig=contig,
                start=base_pos,
                stop=base_pos + 1,
                truncate=True,
                min_base_quality=0,
                multiple_iterators=True,
                stepper="all",
                flag_filter=4 | 256 | 512 | 1024
                ):
                # Get query names and bases at this position
                query_names = pileupcolumn.get_query_names()
                bases = [b.upper() for b in pileupcolumn.get_query_sequences(add_indels=True)]
                for qname, b in zip(query_names, bases):
                    if qname in reads_supporting_alt:
                        if b in base_counts:
                            base_counts[b] += 1
                        else:
                            base_counts["N"] += 1
                # Choose the most frequent base, or N if no coverage
                if sum(base_counts.values()) == 0:
                    consensus_window.append("N")
                else:
                    consensus_base = max(base_counts, key=lambda k: base_counts[k])
                    consensus_window.append(consensus_base)

        consensus_seq = "".join(consensus_window)
        # if 'N' not in consensus_seq:
        print(f">{contig}:{pos}:{ref}>{alt}\n{consensus_seq}")

    bamfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate VCF with solution")
    parser.add_argument("-v", "--vcf", required=True, nargs='+', help="Path(s) to the input VCF file(s) (can use regex, e.g., '*.vcf' or 'dir/*.vcf.gz')")
    parser.add_argument("-r", "--ref", required=True, help="Path to the reference genome file in FASTA format")
    parser.add_argument("-b", "--mapped_reads", required=False, help="Path to the mapped HiFi reads (bam)")
    parser.add_argument("-s", "--solution", required=False, help="Path to the solution file")
    parser.add_argument("-l", "--length", type=int, default=25, help="Length of the k-mers to be generated on each side of the variant (default: 15)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output VCF file with the counts")

    args = parser.parse_args()

    tmp_dir = os.path.join(os.path.dirname(args.output), "tmp")
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass

    # Unzip gzipped VCFs if needed, and collect the list of uncompressed VCFs
    input_vcfs = []
    for vcf_path in args.vcf:
        if vcf_path.endswith(".gz"):
            unzipped_vcf = os.path.join(tmp_dir, os.path.basename(vcf_path).replace(".gz", ""))
            with gzip.open(vcf_path, "rt") as f_in, open(unzipped_vcf, "w") as f_out:
                shutil.copyfileobj(f_in, f_out)
            input_vcfs.append(unzipped_vcf)
        else:
            input_vcfs.append(vcf_path)
    args.vcf = input_vcfs

    reference_genome_available = False
    if not args.solution and not args.mapped_reads:
        print("Error: Either a solution genome file (--solution) or HiFi reads (--reads) must be provided.")
        sys.exit(1)
    elif args.solution:
        reference_genome_available = True

    #to normalize the variants: first use normalize_variants.py, then vt


    if reference_genome_available :
         # nb_of_variants_we_cannot_output = output_all_potential_kmers_represented_by_the_VCF(args.vcf, args.ref, "kmers_to_check.fa", args.length)
          # print("Outputted all variants as k-mers")

        # map_the_potential_kmers_to_the_reference(os.path.join(tmp_dir, "kmers_to_check.fa"), args.solution, os.path.join(tmp_dir, "mapped.txt"), tmp_dir)
        # classify_the_variants(os.path.join(tmp_dir, "mapped.txt"), args.output)
        # print("Did not test ", nb_of_variants_we_cannot_output, " variants because there were too many variants around")

        #create the solution calls for two categories of snps: the one that should be called and the one that can be called
        solution_calls_few = os.path.join(tmp_dir, "solution_few.vcf")
        solution_calls_many = os.path.join(tmp_dir, "solution_many.vcf")
        generate_solution_calls_with_minipileup(args.solution, args.ref, tmp_dir, solution_calls_few, solution_calls_many)
        missing_variants = count_number_of_missing_variants(solution_calls_few, solution_calls_many, args.vcf)
        # create_Logan_queries_to_check_false_positives("false_positive_vcf.vcf", "mapped.bam")
    else:

        if len(args.vcf) > 1 :
            print("ERROR: Several VCFs not supported yet")
        else:
            args.vcf = args.vcf[0]

        mapped_reads_bam = args.mapped_reads

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

        hifi_false_positives = os.path.join(tmp_dir, "hifi_false_positive.vcf")
        compare_calls_with_HiFi(args.vcf, mapped_reads_bam, hifi_false_positives, solution_known=False)
        # create_Logan_queries_to_check_false_positives(hifi_false_positives, "mapped.bam")
        # Output all potential k-mers represented by the VCF
        # output_all_potential_kmers_represented_by_the_VCF(args.vcf, args.ref, os.path.join(tmp_dir, "kmers_to_check.fa"), args.length)
        # print("Outputted all variants as k-mers")







