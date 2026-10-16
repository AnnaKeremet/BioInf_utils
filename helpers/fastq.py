from pathlib import Path

def read_fastq(file)-> tuple[str, str, str, str] | None:
    """
    Read one fastq record

    Arguments:
        - file: Open .fastq file for reading

    Returns tuple[str, str, str, str] | None: FASTQ header, sequence, separator, quality string.
    Returns None when the end of the file is reached
    """
    header=file.readline()
    if not header:
        return None
        
    seq=file.readline()
    separator=file.readline()
    quality=file.readline()
    
    #readline() сохраняет \n в конце каждой строки, если он есть в файле:
    header=header.rstrip('\r\n')
    seq=seq.rstrip('\r\n')
    separator=separator.rstrip('\r\n')
    quality=quality.rstrip('\r\n')
    
    if not all((header,seq,separator,quality)): # https://docs-python.ru/tutorial/vstroennye-funktsii-interpretatora-python/funktsija-all/
        raise ValueError('Одна из строк записи отсутствует в fastq-файле')
    if len(seq) != len(quality):
        raise ValueError('Длина сиквенса не соответствует длине поля качества рида')
    if not header.startswith('@') or not separator.startswith('+'):
        raise ValueError('Некорректны форматы полей заголовка и разделителя')
    return header,seq, separator, quality
        
        
def check_gc(seq:str, gc_bounds: tuple[int | float, int | float] | int | float) -> bool:
    """
    Check whether sequence GC content is within the given bounds

    Arguments:
        - seq: Nucleotide sequence;
        - gc_bounds: Allowed GC content range in percent. A single number specifies the upper boundary

    Returns bool, where True, if GC content is within the bounds, otherwise False
    """
    if isinstance(gc_bounds, (int,float)):
        gc_bounds=(0,gc_bounds)

    seq=str(seq).upper()
    if len(seq)==0:
        return False
    
    gc=seq.count('G') + seq.count('C')
    gc_percent= gc/len(seq) *100
    return gc_bounds[0] <= gc_percent <= gc_bounds[1]
        
        
def check_length(seq:str, length_bounds:tuple[int, int] | int | float)-> bool:
    """
    Check whether sequence length is within the given boundaries

    Arguments:
        - seq: Nucleotide sequence;
        - length_bounds: Allowed sequence length range. A single number specifies the upper boundary.

    Returns bool, where True if sequence length is within the bounds, otherwise False
    """
    if isinstance(length_bounds, (int,float)):
            length_bounds=(0,length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]
    
    
def check_quality(quality:str, quality_threshold:int | float)-> bool:
    """
    Check whether mean Phred+33 quality is filterd by threshold.

    Arguments:
        - quality: String containing ASCII-encoded quality scores;
        - quality_threshold: Minimum allowed mean Phred quality.

    Returns bool, where True if mean quality is at least the threshold, otherwise False.
    """
    if len(quality) == 0:
        return False
    
    sum_quality=0
    for symbol in quality:
        symbol=ord(symbol) # https://www.tutorialspoint.com/article/how-to-convert-a-single-character-to-its-integer-value-in-python
        phred= symbol - 33
        sum_quality += phred
    mean_quality=sum_quality/len(quality)
    return mean_quality >= quality_threshold
        
        
def write_filtered_fastq(file, record:tuple[str, str, str, str])-> None:
    """
    Write one .fastq record to an open output file.

    Arguments:
        - file: Open output file for writing;
        - record: Tuple containing the four fastq record lines.

    Returns None, but writes the fastq record to the output file
    """
    for line in record:
        file.write(line + '\n')
        
