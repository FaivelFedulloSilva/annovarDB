import os
import logging


def file_exits(file_path):
    if os.path.exists(file_path):
        return True
    else:
        print("The file " + file_path.split('/')[-1] + " does not exist.")

    










