
# From a line of a VCF file, makes an object with the usefull information
def getChr (vcf_line):
    line = vcf_line.split('\t')
    return line[0]

def removeLeftContext(ref,alt):
    context = True
    while len(ref) and len(alt) and context:
        if ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
        else:
            context = False
    return ref, alt

def adjustAllelicInfo(allelic_info, alt_idx):
    info = allelic_info.split(':')

    first_allele = info[0][0]
    second_allele = info[0][2]
    
    if first_allele == str(alt_idx):
        first_allele = '1'
    elif first_allele == '.':
        first_allele == '.'
    else:
        first_allele = '0'

    if second_allele == str(alt_idx):
        second_allele = '1'
    elif second_allele == '.':
        second_allele == '.'
    else:
        second_allele = '0'

    new_allelic_info = info.copy()

    new_allelic_info[0] = first_allele + '/' + second_allele
    if len(new_allelic_info) == 1:
        new_allelic_info[0] += '\n' 
    new_info = ':'.join(new_allelic_info)
    return new_info

def makeBiallelic(vcf_line):
    line = vcf_line.split('\t')
    alts = line[4].split(',')
    lines = []
    alt_idx = 1
    for alt in alts:
        aux_line = line.copy()
        new_ref, new_alt = aux_line[3], alt # removeLeftContext(aux_line[3], alt)
        aux_line[3] = new_ref
        aux_line[4] = new_alt
        if len(alts) > 1:
            for i in range(9, len(aux_line)):
                aux_line[i] = adjustAllelicInfo(aux_line[i], alt_idx)
        lines.append('\t'.join(aux_line))
        alt_idx += 1
    return lines    


TEMP_PATH = './temp/'
VCF_FILE_PATH = './input/filtered_3A.vcf'
CHR_ORDER_PATH = './output/chr_order.txt'
OUTPUT_FOLDER = './output/'


common_header = []
chr_order = []


with open(VCF_FILE_PATH, 'r') as vcf:
    vcf_line = vcf.readline()
    while vcf_line[0] == "#":
        common_header.append(vcf_line)
        vcf_line = vcf.readline()

    while vcf_line:
        line_chr = getChr(vcf_line)
        print(line_chr)
        if line_chr not in chr_order:
            chr_order.append(line_chr)
        with open(f'{TEMP_PATH}output_{line_chr}.vcf', 'x') as chr_file:
            for h in common_header:
                chr_file.write(h)

            while line_chr in chr_order:
                vcf_lines = makeBiallelic(vcf_line)
                for line in vcf_lines:
                    chr_file.write(line)

                vcf_line = vcf.readline()
                line_chr = getChr(vcf_line)

print(chr_order)
