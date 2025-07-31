import argparse
import glob

def normalize_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fields = line.strip().split()
                if len(fields)>=8 :
                    chrom, pos, id_, ref, alt, qual, info, filter_ = fields[:8]
                    
                    for a, alt_allele in enumerate(alt.split(',')):
                        if input_vcf == "bcfcalls.vcf":
                            new_line = '\t'.join([chrom, pos, id_, ref, alt_allele, "AD="+fields[7].split(';')[0].split('=')[-1]]) + '\n'
                        elif input_vcf == "longshotcalls.vcf" or input_vcf == "nanocaller.vcf" or  input_vcf == "deepvariant.vcf" or input_vcf == "clair.vcf":
                            new_line = '\t'.join([chrom, pos, id_, ref, alt_allele, "AD="+fields[9].split(':')[3].split(',')[a+1]]) + '\n'
                        elif 'metacaller' in input_vcf:
                            new_line = '\t'.join([chrom, pos, id_, ref, alt_allele, "AD="+info.lstrip("DP=").split(',')[a]]) + '\n'
                        else:
                            print("ERROR: which caller is that?")
                        outfile.write(new_line)


def main():
    parser = argparse.ArgumentParser(description="Normalize one or more VCF files.")
    parser.add_argument("input_vcf", nargs='+', help="Path(s) to the input VCF file(s).")
    args = parser.parse_args()

    # Only include files that end with .vcf
    vcf_files = []
    for pattern in args.input_vcf:
        vcf_files.extend([f for f in glob.glob(pattern) if f.endswith('.vcf')])
    if not vcf_files:
        print(f"No files matched the pattern(s): {', '.join(args.input_vcf)}")
        return

    print("vcfiles: ", vcf_files)

    for vcf in vcf_files:
        print("normalizign ", vcf)
        if vcf.endswith('.vcf'):
            output_vcf = vcf[:-4] + '_normalized.vcf'
        else:
            output_vcf = vcf + '_normalized'
        normalize_vcf(vcf, output_vcf)

if __name__ == "__main__":
    main()