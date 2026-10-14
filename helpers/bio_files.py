def write_fasta_record(output_file, header, seq_lines):
    output_file.write(header + '\n')
    output_file.write(''.join(seq_lines)+'\n') # https://skillbox.ru/media/code/kak-preobrazovat-spisok-v-stroku-v-python/
    
    