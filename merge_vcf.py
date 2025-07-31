import argparse
from collections import defaultdict

def merge_vcfs(output_file, input_files):
    """
    Merge multiple VCF files into one, assuming they share the same reference.
    """
    
    # Dictionary to store variants indexed by contig and position
    variants = defaultdict(list)
    datasets = {}

    # Read each input VCF file
    for file in input_files:
        with open(file, 'r') as vcf:
            for line in vcf:
                # Skip header lines
                if line.startswith("#"):
                    continue
                
                # print(line, " iiiiuiuhuh")
                # Parse VCF line
                parts = line.strip().split("\t")
                contig = parts[0]
                position = int(parts[1])
                ref = parts[3].upper()
                alt = parts[4].upper()

                depth = "*"
                
                for info_field in parts[5].split(";"):
                    if info_field.startswith("DP=") or info_field.startswith("AD="):
                        depth = info_field.split("=")[1]
                        break
        
                # Construct a simplified variant line with file of origin
                simplified_variant_info = f"{contig}\t{position}\t.\t{ref}\t{alt}\t.\t"
                if simplified_variant_info not in datasets:
                    datasets[simplified_variant_info] = []
                    variants[(contig, position)].append(simplified_variant_info)
                datasets[simplified_variant_info].append((file,depth))
                
                

    # Write merged variants to the output file
    with open(output_file, 'w') as out_vcf:
        # Write a minimal VCF header
        out_vcf.write("##fileformat=VCFv4.2\n")
        out_vcf.write("#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFILES\n")
        
        # Write variants sorted by contig and position
        for (contig, position), entries in sorted(variants.items()):
            for entry in entries:
                out_vcf.write(entry+"\t")
                out_vcf.write(",".join(file for file, _ in datasets[entry]) + "\t")
                out_vcf.write(",".join(depth for _, depth in datasets[entry]) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Merge multiple VCF files sharing the same reference.")
    parser.add_argument("-o", "--output", required=True, help="Output merged VCF file.")
    parser.add_argument("input_files", nargs="+", help="Input VCF files to merge.")
    args = parser.parse_args()

    merge_vcfs(args.output, args.input_files)

if __name__ == "__main__":
    main()
