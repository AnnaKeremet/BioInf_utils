def write_fasta_record(output_file, header:str, seq_lines:list[str]) -> None:
    """
    Write a fasta record with the sequence on a single line

    Arguments:
        - output_file: Open output file for writing;
        - header: fasta sequence header;
        - seq_lines: List of sequence fragments.

    Returns None, but writes the fasta  record to the separate output file
    """
    output_file.write(header + '\n')
    output_file.write(''.join(seq_lines)+'\n') # https://skillbox.ru/media/code/kak-preobrazovat-spisok-v-stroku-v-python/
    

def get_blast_description(line:str, description_end:int) -> str:
    """
    Extract the protein description from a BLAST file from table with proteins names

    Arguments:
        - line: Text line containing a BLAST hit;
        - description_end: Index where the next table column begins.

    Returns str, which is a protein description without surrounding whitespace
    """
    description=line[:description_end].strip()
    return description