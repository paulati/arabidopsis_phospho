import gzip
import shutil
from os import path

def extractgz(input_file_path, output_file_path):
    with gzip.open(input_file_path, 'rb') as f_in:
        with open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
    
def compressgz(input_file_path, output_file_path):
    with open(input_file_path, 'rb') as f_in:
        with gzip.open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
            
def extractzip(input_file_path, output_folder_path):
    shutil.unpack_archive(input_file_path, output_folder_path)            
    
def compresszip(input_file_path, output_folder_path):    
    base_dir = path.basename(input_file_path)        
    root_dir = output_folder_path        
    shutil.make_archive(output_file_path, 'zip', root_dir, base_dir)
    
