from pathlib import Path
from helpers.fastq import (read_fastq, 
                           check_gc, 
                           check_length, 
                           check_quality,
                           write_filtered_fastq)
from helpers.bio_files import (write_fasta_record, get_blast_description)

def filter_fastq(input_fastq:str,
                 output_fastq:str,
                 gc_bounds: tuple[int | float, int | float] | int | float = (0,100),
                 length_bounds: tuple[int | int] | int | float= (0, 2**32),
                 quality_threshold: int|float =0) -> None:
    """
    Filters .fastq reads by GC-content, length and mean quality

    Arguments:
        - input_fastq: Path to the input .fastq file;
        - output_fastq: Name of the output .fastq file;
        - gc_bounds: Filtration of reads by GC-content, allows reads with GC-range in percent (0-100).
        A single number - upper boundary;
        - length_bounds: Allowed read length range. A single number - upper boundary;
        - quality_threshold: Minimum mean Phred+33 quality score.

    Returns None, but saves filtered reads to the "filtered" folder
    """
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


def convert_multiline_fasta_to_oneline(input_fasta:str, 
                                       output_fasta: str | Path | None = None) -> None:
    """
    Converts multiline .fasta to single-line .fasta

    Argumentss:
        - input_fasta: Path to the input .fasta file;
        - output_fasta: Path to the output .fasta file. If None, the output filename is generated automatically.

    Returns None, but writes the converted sequences to .fasta file with same name+ "_oneline"
    """
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
                

def parse_blast_output(input_file:str, output_file:str) -> None:
    """
    Extracts the best BLAST hit description for each query in .txt file

    Arguments:
        - input_file: Path to the input BLAST text file;
        - output_file: Path to the output text file.

    Returns None, but saves protein descriptions sorted alphabetically in separate file
    """
    description_end=None
    description_part=False
    proteins=[]
    
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('Query #'):
                description_end=None
                description_part=False
            elif "Sequences producing significant alignments:" in line:
                description_part=True
            elif description_part and "Scientific" in line and "Common" in line:
                description_end= line.index("Scientific") # https://www.w3schools.com/python/ref_string_index.asp
            elif description_part and line.startswith('Description'):
                continue
            elif description_part and description_end != None and line.strip():
                protein=get_blast_description(line, description_end)
                proteins.append(protein)
                description_part=False
    proteins.sort() # изменяет исходный список
    
    with open(output_file, 'w') as file:
        for protein in proteins:
            file.write(protein + '\n')        
                
                
                           
                
                    
                
                
            
            