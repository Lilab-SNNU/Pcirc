def double(filename):
    path = filename.replace('aligned_seq', 'aligned_seq_double')
    file = open(filename, 'r')
    output_file = open(path, 'w')
    for line in file:
        if line[0] == '>':
            header = line
        elif line[0] != '>':
            seq = line.rstrip('\n') * 2 + '\n'
            output_file.write(header + seq)
    output_file.close()
    return path


