from email.policy import default
import os
import shutil
from sqlite3 import Time
import utilsDB as utils
import sys
import logging
import createDB
import subprocess
import argparse
import time 

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name)
        print('\n\nElapsed: %s \n\n' % (time.time() - self.tstart))


parser = argparse.ArgumentParser(
    description = "From a vcf file, create a generic database to be use as a database to annote vcf with annovar",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument(
    "-i","--in","--input", 
    type=str, 
    required=True, 
    dest='input',
    help="help" #TODO Make help text for input argument
)
parser.add_argument(
    "-o","--out","--output", 
    type=str, 
    required=True, 
    dest='output',
    help="help" #TODO Make help text for output argument
)

parser.add_argument(
    "--no-log",   
    default=True,
    action="store_false",
    dest='log',
    help="help" #TODO Make help text for log argument
)

parser.add_argument(
    "--annovar-folder",   
    type=str,
    dest='annovar_folder',
    help="help" #TODO Make help text for annovar-folder argument
)

parser.add_argument(
    "--plink-folder",   
    type=str,
    dest='plink_folder',
    help="help" #TODO Make help text for plink-folder argument
)

parser.add_argument(
    "--annoted-file",   
    type=str,
    dest='annoted',
    help="help" #TODO Make help text for annoted_file argument
)

parser.add_argument(
    "--freq-file",   
    type=str,
    dest='freq',
    help="help" #TODO Make help text for freq_file argument
)

parser.add_argument(
    "--keep-tmp-files",   
    default=False,
    action='store_true',
    dest='keep_tmp',
    help="help" #TODO Make help text for freq_file argument
)

args = parser.parse_args()
print(args)



TEST = 1
CHR_ORDER_PATH = './temp/chr_order.txt'
chr_order = [
        # 'output_chrM.vcf',
        'output_chr1.vcf',
        'output_chr2.vcf',
        'output_chr3.vcf', 
        'output_chr4.vcf', 
        'output_chr5.vcf', 
        'output_chr6.vcf', 
        'output_chr7.vcf', 
        'output_chr8.vcf', 
        'output_chr9.vcf', 
        'output_chr10.vcf', 
        'output_chr11.vcf', 
        'output_chr12.vcf',
        'output_chr13.vcf',
        'output_chr14.vcf',
        'output_chr15.vcf',
        'output_chr16.vcf',
        'output_chr17.vcf',
        'output_chr18.vcf',
        'output_chr19.vcf',
        'output_chr20.vcf',
        'output_chr21.vcf',
        'output_chr22.vcf',
        'output_chrX.vcf']
# por ahora lo escribo aca
with open(CHR_ORDER_PATH, 'x') as chr_order_file:
    for chr in chr_order:
        chr_order_file.write(chr+'\n')

with open(CHR_ORDER_PATH, 'r') as chr_order_file:
    chr_order = chr_order_file.readlines()


for chr_file in chr_order:
    if not utils.file_exits('./tmp'):
        os.mkdir('./tmp')

    with Timer():
        if TEST == 0:
            if args.annovar_folder:
                if not utils.file_exits("./humandb/hg19_clinvar_20180603.txt"):
                    script = f"{args.annovar_folder}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/ "
                    subprocess.run(script, shell=True)
                script = f"{args.annovar_folder}/table_annovar.pl {args.input} humandb/ -buildver hg19 -out ./tmp/annovar_file -remove -protocol clinvar_20180603 -operation f -nastring . -vcfinput"
                subprocess.run(script, shell=True)

            if args.plink_folder:
                script = f"{args.plink_folder}/plink --vcf {args.input} --freq --out ./tmp/freq_file --keep-allele-order"
                subprocess.run(script, shell=True)
        elif TEST == 1:
            if args.plink_folder:
                script = f"{args.plink_folder}/plink --vcf ./temp/{chr_file[:-1]} --make-bed --out ./tmp/bed_temp_file --keep-allele-order"
                subprocess.run(script, shell=True)
                script = f"{args.plink_folder}/plink --bfile ./tmp/bed_temp_file --recode vcf --out ./tmp/new_vcf --keep-allele-order"
                subprocess.run(script, shell=True)
                script = f"{args.plink_folder}/plink --vcf ./temp/{chr_file[:-1]} --freq --out ./tmp/freq_file --keep-allele-order"
                subprocess.run(script, shell=True)

            if args.annovar_folder:
                if not utils.file_exits("./humandb/hg19_clinvar_20180603.txt"):
                    script = f"{args.annovar_folder}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/ "
                    subprocess.run(script, shell=True)
                script = f"{args.annovar_folder}/table_annovar.pl ./tmp/new_vcf.vcf humandb/ -buildver hg19 -out ./tmp/annovar_file -remove -protocol clinvar_20180603 -operation f -nastring . -vcfinput"
                subprocess.run(script, shell=True)

        FREQ_FILE_PATH = args.freq if args.freq else './tmp/freq_file.frq'
        VCF_FILE_PATH = f'./temp/{chr_file[:-1]}'
        ANNOVAR_FILE_PATH = args.annoted if args.annoted else './tmp/annovar_file.hg19_multianno.txt'

        # createDB.makeTempFile(VCF_FILE_PATH, FREQ_FILE_PATH)
        # createDB.makeGenericDB(ANNOVAR_FILE_PATH)

        db_builder = createDB.DBBuilder()
        db_builder.set_vcf_file(VCF_FILE_PATH)
        db_builder.set_annoted_file(ANNOVAR_FILE_PATH)
        db_builder.set_freq_file(FREQ_FILE_PATH)
        db_builder.set_output_name(args.output)
        db_builder.set_append_mode(chr_file != chr_order[0])
        db_builder.makeGenericDB()
        shutil.rmtree('./tmp')

if not args.keep_tmp:
    shutil.rmtree('./tmp')

# # -in, -input, -i   Path to vcf file
# # -annovar          Path to annovar folder
# # -plink            Path to plink folder
# # -log              Log to a file, to console, both or none
# # -keep-tmp-files   Do not erase the temporary files
# # -out, -output, -o Output name of the DB           
# # -can-update       Keeps the number of observations in the DB file, otherwise create a diferent file with this information.
# # -snp              Path to snp db folder
# # -annovar_file     Path to Annovar multianno file
# # -plink-file       Path to freq plink file
# # folder_structure  Use folder structure.

# # Use Cases
# # 1 - The user provides the plink, annovar and vcf files. The program use them to proceed
# # 2 - The user provides the plink and vcf file. The program must create the annovar file to proceed
# # 3 - The user provides the annovar and vcf file. The program must create the plink freq file to proceed
# # 4 - the user provides the vcf file. The program must create the annovar and plink freq file to proceed
# # 5 - 


