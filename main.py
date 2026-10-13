from pathlib import Path
from helpers.fastq import (read_fastq, 
                           check_gc, 
                           check_length, 
                           check_quality,
                           write_filtered_fastq)

def filter_fastq(input_fastq,
                 output_fastq,
                 gc_bounds = (0,100),
                 length_bounds= (0, 2**32),
                 quality_threshold=0):
    output_dir=Path("filtered")
    output_dir.mkdir(exist_ok=True) # https://docs.python.org/3/library/pathlib.html    
    output_path=output_dir/Path(output_fastq).name # https://docs.python.org/3/library/pathlib.html    
    
    with open(input_fastq, 'r') as file:
        with open(output_path, "w") as output_file:
            while True:
                record=read_fastq(file)
                if record == None:
                    break
                
                header = record[0],
                seq = record[1], 
                separator = record[2], 
                quality = record[3]
                
                if check_gc(seq, gc_bounds) and check_length(seq, length_bounds) and check_quality(quality, quality_threshold):
                    write_filtered_fastq(output_file,record)
            
            
            