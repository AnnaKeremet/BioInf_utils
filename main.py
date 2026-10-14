from pathlib import Path
from helpers.fastq import (read_fastq, 
                           check_gc, 
                           check_length, 
                           check_quality,
                           write_filtered_fastq)
from helpers.bio_files import (write_fasta_record,
                               )

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


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    if output_fasta == None:
        input_path=Path(input_fasta)
        output_fasta=input_path.with_name(input_path.stem + "_oneline.fasta") #https://docs.python.org/3/library/pathlib.html
    
    with open(input_fasta, 'r') as input_file:
        with open(output_fasta, 'w') as output_file:
            header=None
            seq_line=[]        
            for line in input_file:
                line=line.strip() # https://ru.stackoverflow.com/questions/1143410/%D0%9A%D0%B0%D0%BA-%D1%83%D0%B1%D1%80%D0%B0%D1%82%D1%8C-%D0%BF%D0%B5%D1%80%D0%B5%D0%BD%D0%BE%D1%81%D1%8B-%D1%81%D1%82%D1%80%D0%BE%D0%BA%D0%B8-%D0%B8%D0%B7-%D1%81%D0%BF%D0%B8%D1%81%D0%BA%D0%B0-%D0%BF%D1%80%D0%B8-%D0%BE%D0%B1%D1%80%D0%B0%D0%B1%D0%BE%D1%82%D0%BA%D0%B5-%D0%B2-python3
                if line.startswith('>'):
                    if header != None:
                        write_fasta_record(output_file, header, seq_line)
                    header=line
                    seq_line=[]
                else:
                    seq_line.append(line)
            if header != None:
                write_fasta_record(output_file, header, seq_line) # записываем последнюю посл-ть в документе                            
                
                    
                
                
            
            