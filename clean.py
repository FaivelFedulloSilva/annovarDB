import shutil
import os
import utilsDB as utils

if utils.file_exits('./log.ans'):
    os.remove('./log.ans')

if utils.file_exits('./temp'):
    shutil.rmtree('./temp')

if utils.file_exits('./output/newDB.txt'):
    os.remove('./output/newDB.txt')

if utils.file_exits('./tmp'):
    shutil.rmtree('./tmp')

if not utils.file_exits('./temp'):
    os.mkdir('./temp')
#os.remove('./temp/chr_order.txt')