import re
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



# From a line of a plink .frq file, makes an object with the usefull information
def makeFreqObject(freq_line):
    line = freq_line.split('\t')
    freq_object = {
        'chr': line[0],	
        'snp': line[1],	
        'alt': line[2] if line[2] != '*' else '0',	
        'ref': line[3],	
        'maf': line[4],	
        'nchrobs': line[5]
    } 
    return freq_object

# From a line of a VCF file, makes an object with the usefull information

def makeVcfObject (vcf_line):
    line = vcf_line.split('\t')
    vcf_object = {
        'chr': line[0],	
        'pos': line[1],
        'snp': line[2],	
        'ref': line[3] if line[3] != '.' and line[3] != '*' else '0',	
        'alt': line[4] if line[4] != '.' and line[4] != '*'  else '0',	
    }

    if vcf_object['chr'] == 'M' or  vcf_object['chr'] == 'chrM' or vcf_object['chr'] == 'm' or vcf_object['chr'] == 'ChrM':
        vcf_object['chr'] = '26'
    if vcf_object['chr'] == 'Y' or  vcf_object['chr'] == 'chrY' or vcf_object['chr'] == 'y' or vcf_object['chr'] == 'ChrY':
        vcf_object['chr'] = '24'
    if vcf_object['chr'] == 'X' or  vcf_object['chr'] == 'chrX' or vcf_object['chr'] == 'x' or vcf_object['chr'] == 'ChrX':
        vcf_object['chr'] = '23'
    for x in range(1, 23):
        if vcf_object['chr'] == f'Chr{x}' or  vcf_object['chr'] == f'chr{x}':
            vcf_object['chr'] = str(x)       
    return vcf_object

# From a line of a Merge file, makes an object with the usefull information
def makeMergeObject (merge_line):
    line = merge_line.split('\t')
    merge_object = {
        'chr': line[0],	
        'pos': line[1],
        'snp': line[2],	
        'ref': line[3] if line[3] != '.' else '0',	
        'alt': line[4] if line[4] != '.' else '0',	
        'frq': line[5],
        'nchrobs': line[6]
    }
    return merge_object

# From a line of a annovar annoted vcf .txt file, makes an object with the usefull information
def makeAnnovarObject(annovar_line):
    line = annovar_line.split('\t')
    annovar_object = {
        'chr': line[0],	
        'start': line[1],
        'end': line[2],	
        'ref': line[3] if line[3] != '-' else '',	
        'alt': line[4] if line[4] != '-' else '',	
    }
    if annovar_object['chr'] == 'M' or  annovar_object['chr'] == 'chrM' or annovar_object['chr'] == 'm' or annovar_object['chr'] == 'ChrM':
        annovar_object['chr'] = '26'
    if annovar_object['chr'] == 'Y' or  annovar_object['chr'] == 'chrY' or annovar_object['chr'] == 'y' or annovar_object['chr'] == 'ChrY':
        annovar_object['chr'] = '24'
    if annovar_object['chr'] == 'X' or  annovar_object['chr'] == 'chrX' or annovar_object['chr'] == 'x' or annovar_object['chr'] == 'ChrX':
        annovar_object['chr'] = '23'
    for x in range(1, 23):
        if annovar_object['chr'] == f'Chr{x}' or  annovar_object['chr'] == f'chr{x}':
            annovar_object['chr'] = str(x)   
    return annovar_object
# Compare a line of a VCF file and a line of a PLINK .frq file and return true if they correspond with the same variant. False otherwise
def compareVCFandFreq(vcf_object, freq_object):
    vcf_pseudo_id = vcf_object['snp'] + vcf_object['chr'] + vcf_object['ref'] + vcf_object['alt']
    freq_pseudo_id = freq_object['snp'] + freq_object['chr'] + freq_object['ref'] + freq_object['alt']
    return vcf_pseudo_id == freq_pseudo_id

# Given a line annovar annoted file and a line from a merge file, 
# verify if the difference between them corrispond with the modifications 
# made from annovar when annoting a deletion 
def isDeletion(annovar_object, merge_object):
    start = int(annovar_object['start']) == int(merge_object['pos']) + 1
    ref_len = len(annovar_object['ref']) == len(merge_object['ref']) - 1
    alt_len =  len(annovar_object['alt']) == len(merge_object['alt']) - 1

    return start and ref_len and alt_len

# Given a line annovar annoted file and a line from a merge file, 
# verify if the difference between them corrispond with the modifications 
# made from annovar when annoting a insertion 
def isInsertion(annovar_object, merge_object):
    start = int(annovar_object['start']) == int(merge_object['pos'])
    ref_len = len(annovar_object['ref']) == len(merge_object['ref']) - 1
    alt_len =  len(annovar_object['alt']) == len(merge_object['alt']) - 1

    return start and ref_len and alt_len

# Given a line annovar annoted file and a line from a merge file, 
# verify if the difference between them corrispond with the modifications 
# made from annovar when annoting a SNV 
def isSNV(annovar_object, merge_object):
    annovar_pseudo_id = annovar_object['chr'] + annovar_object['start'] + annovar_object['ref'] + annovar_object['alt']
    merge_pseudo_id = merge_object['chr'] + merge_object['pos'] + merge_object['ref'] + merge_object['alt']

    return annovar_pseudo_id == merge_pseudo_id

def deepCompare(annovar_object, merge_object):
    context_length = int(annovar_object['start']) - int(merge_object['pos'])
    if len(annovar_object['ref']) == 0:
        context_length += 1 
    txt = annovar_object['ref'] == merge_object['ref'][context_length:]
    if not txt:
    #     print(f'{bcolors.OKGREEN} True {bcolors.ENDC}')
    # else:
        print(context_length, annovar_object['ref'], merge_object['ref'][context_length:]) 

        # print(f'{bcolors.FAIL} False {bcolors.ENDC}')
        if len(annovar_object['ref']) > len(annovar_object['alt']):
            print(f'{bcolors.FAIL}>>>>> DELETION <<<<<{bcolors.ENDC}')
        if len(annovar_object['ref']) < len(annovar_object['alt']):
            print(f'{bcolors.OKGREEN}>>>>> INSERTION <<<<<{bcolors.ENDC}')
        if len(annovar_object['ref']) == len(annovar_object['alt']):
            print(f'{bcolors.OKBLUE}>>>>> SUBSTITUTION <<<<<{bcolors.ENDC}')
    return txt

# compare a line from a annoted annovar file and a line from a merge file. 
# If the lines correspond to the same variant it will return True. 
# False otherwise
def compareAnnovarAndMerge(annovar_object, merge_object):
    insertion = isInsertion(annovar_object, merge_object)
    deletion = isDeletion(annovar_object, merge_object)
    snv = isSNV(annovar_object, merge_object)
    if not (insertion or deletion or snv):
        
        compare = deepCompare(annovar_object, merge_object)
        if not compare:
            print('================================================')
            print(insertion, deletion, snv)
            print(annovar_object)
            print(merge_object)
            print('================================================\n\n')
    return insertion or deletion or snv

# Given a vcf_object and a freq_object, returns an 
# string containing all the data delimited with tabs
def makeMergeLine(vcf_object, freq_object):
    return vcf_object['chr'] + '\t' + vcf_object['pos'] + '\t' + vcf_object['snp'] + '\t' + vcf_object['ref'] + '\t' + vcf_object['alt'] + '\t' + freq_object['maf'] + '\t' + freq_object['nchrobs']


# Given a annovar_object and a merge_object, returns a 
# string containing all the data delimited with tabs
def makeGenericDBLine(annovar_object, merge_object):
    dbLine = annovar_object['chr'] + '\t' + annovar_object['start'] + '\t' + annovar_object['end'] + '\t' 
    dbLine += annovar_object['ref'] if len(annovar_object['ref']) > 0 else '-'
    dbLine += '\t'
    dbLine += annovar_object['alt'] if len(annovar_object['alt']) > 0 else '-' 
    dbLine += '\t' + merge_object['frq'] + '\t' + merge_object['nchrobs']
    return  dbLine
                    
class DBBuilder:
    def __init__(self, log = True, tmp='./tmp'):
        self.log = log
        self.tmp = tmp
        self.TEMP_TAB_FREQ_PATH = f'{tmp}/tabfrq.frq'
        self.TEMP_MERGE_PATH = f'{tmp}/mergefile.txt'
        self.TEMP_CHR_PATH = f'{tmp}/chr_order.txt'
        self.make_chr_order = False
        self.chr_order = []
        self.file_operation = 'x'
    
    def set_vcf_file(self, vcf_file):
        self.vcf_file = vcf_file
    
    def set_annoted_file(self, annoted_file):
        self.annoted_file = annoted_file

    def set_freq_file(self, freq_file):
        self.freq_file = freq_file
    
    def set_output_name(self, output):
        self.output_name = output

    def set_chr_ord(self, make_chr_order_file):
        self.make_chr_order = make_chr_order_file
    
    def set_append_mode(self, append):
        self.file_operation = 'a' if append else 'x'

    def makeTabFreqFile(self):
        with open(self.freq_file, 'r') as freq_file:
            with open(self.TEMP_TAB_FREQ_PATH, 'x') as tab_freq_file:
                original_line = freq_file.readline()
                while(original_line):
                    tab_freq_file.write(re.sub(" +", "\t", original_line.lstrip()))
                    original_line = freq_file.readline()

    def makeMergeFile(self):
        self.makeTabFreqFile()

        with open(self.vcf_file, 'r') as original_vcf_file:
            with open(self.TEMP_TAB_FREQ_PATH, 'r') as tab_freq_file:
                with open(self.TEMP_MERGE_PATH, 'x') as merge_file:

                    # Read through the header
                    vcf_line = original_vcf_file.readline()
                    while vcf_line[0:2] == '##':
                        vcf_line = original_vcf_file.readline()
                    vcf_line = original_vcf_file.readline()
                    freq_line = tab_freq_file.readline()
                    freq_line = tab_freq_file.readline()

                    while vcf_line and freq_line:
                        vcf_object = makeVcfObject(vcf_line)
                        freq_object = makeFreqObject(freq_line)
                        
                        if vcf_object['chr'] not in self.chr_order: 
                            self.chr_order.append(vcf_object['chr']) 

                        if compareVCFandFreq(vcf_object, freq_object):
                            merge_file.write(makeMergeLine(vcf_object, freq_object))
                        # else:
                        #     print("ERROR -> the lines does not correspond to the same variant")
                        #     print('>>> ' + str(vcf_object))
                        #     print()
                        #     print('>>> ' + str(freq_object))
                        #     print('===============================================\n\n')    
                        freq_line = tab_freq_file.readline()
                        vcf_line = original_vcf_file.readline()
        if self.make_chr_order:
            with open(f'{self.tmp}/chr_order.txt', 'x') as file:
                for chr in self.chr_order:
                    file.write(chr + '\n')

    def makeGenericDB(self):
        self.makeMergeFile()
        with open(self.annoted_file, 'r') as annovar_file:
            with open(f'{self.tmp}/mergefile.txt', 'r') as merge_file:
                with open(f'{self.output_name}.txt', self.file_operation) as db_file:
                    # Read 2 lines to discard the header
                    annovar_line = annovar_file.readline() 
                    annovar_line = annovar_file.readline()
                    merge_line = merge_file.readline()

                    while annovar_line and merge_line:    
                        annovar_object = makeAnnovarObject(annovar_line)
                        merge_object = makeMergeObject(merge_line)

                        if compareAnnovarAndMerge(annovar_object, merge_object):
                            db_file.write(makeGenericDBLine(annovar_object, merge_object))
                            
                        merge_line = merge_file.readline()
                        annovar_line = annovar_file.readline()