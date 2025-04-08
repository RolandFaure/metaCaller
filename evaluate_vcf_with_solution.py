#The goal of this code is to count how many variants have been correctly and incorrectly identified in a VCF? given the original genomes

import sys
import argparse
import os
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
            extensions_right = [("", pos+1)] #extensions are a pair (seq, pos_on_ref)
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
                fo.write(">"+contig+"$"+str(variant_pos)+"$"+variant[0]+"$"+variant[1]+"$-1\n")
                fo.write("NNNNN\n")
                continue
            for variant in variants[(contig, variant_pos)]:
                index_seq_for_this_snp = 0
                print("len ", extensions_of_every_base_left[variant_pos], " ", extensions_of_every_base_right[variant_pos+len(variant[0])-1])
                for left_ext in extensions_of_every_base_left[variant_pos]:
                    for right_ext in extensions_of_every_base_right[variant_pos+len(variant[0])-1]:
                        fo.write(">"+contig+"$"+str(variant_pos)+"$"+variant[0]+"$"+variant[1]+"$" + str(index_seq_for_this_snp) +"\n")
                        fo.write(left_ext+variant[1]+right_ext+"\n")
                        index_seq_for_this_snp+=1
                fo.write(">"+contig+"$"+str(variant_pos)+"$"+variant[0]+"$REF\n")
                fo.write(left_ext+ variant[0] +right_ext+"\n")
                index_snp +=1
                sys.exit()


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



# tmp_dir = "/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/out/tmp"
# try :
#     os.mkdir(tmp_dir)
# except: #if dir exist, good
#     ... 
# output_all_potential_kmers_represented_by_the_VCF("/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/out_trash/variants.vcf", "/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/one_strain.fa", "/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/out/all_seqs_to_test.fa", 15)
# map_the_potential_kmers_to_the_reference("/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/out/all_seqs_to_test.fa", "/home/rfaure/Documents/these/metaCaller/species/10_Ecoli/ref.fasta", tmp_dir+"/mapped.txt", tmp_dir)
# classify_the_variants(tmp_dir+"/mapped.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate VCF with solution")
    parser.add_argument("-v", "--vcf", required=True, help="Path to the input VCF file containing variant information")
    parser.add_argument("-r", "--ref", required=True, help="Path to the reference genome file in FASTA format")
    # parser.add_argument("-s", "--solution", required=True, help="Path to the solution file")
    parser.add_argument("-l", "--length", type=int, default=15, help="Length of the k-mers to be generated on each side of the variant (default: 15)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file where the output kmers will be written")

    args = parser.parse_args()

    tmp_dir = os.path.join(os.path.dirname(args.output), "tmp")
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass

    nb_of_variants_we_cannot_output = output_all_potential_kmers_represented_by_the_VCF(args.vcf, args.ref, args.output, args.length)
    # map_the_potential_kmers_to_the_reference(os.path.join(tmp_dir, "all_seqs_to_test.fa"), args.solution, os.path.join(tmp_dir, "mapped.txt"), tmp_dir)
    # classify_the_variants(os.path.join(tmp_dir, "mapped.txt"), args.output)
    # print("Did not test ", nb_of_variants_we_cannot_output, " variants because there were too many variants around")

    # solution_calls = os.path.join(tmp_dir, "solution.vcf")
    # generate_solution_calls_with_minipileup(args.solution, args.ref, tmp_dir, solution_calls)
    # missing_variants = count_number_of_missing_variants(solution_calls, args.vcf)




