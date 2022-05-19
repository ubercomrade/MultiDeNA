import random


def read_fasta(path_in, path_out):
    fasta = list()
    append = fasta.append
    fasta_in = open(path_in, 'r')
    fasta_out = open(path_out, 'w')
    for index, line in enumerate(fasta_in):
        if not line.startswith('>'):
            line = line.strip().upper()
            line = clear_n(line)
            if line != '':
                fasta_out.write('>{0}\n'.format(int(index / 2)))
                fasta_out.write(line + '\n')       
    fasta_in.close()
    fasta_out.close()
    pass


def longest(ss):
    if len(ss[0]) > len(ss[1]):
        return(ss[0])
    else:
        return(ss[1])

    
def clear_n(string):
    while 1:
        position = string.find('N')
        if position == -1:
            break
        elif position == len(string) - 1:
            string = string[:position - 1]
            break
        elif string[position + 1] != 'N':
            string = string[:position] + random.choice('ACGT') + string[position + 1:]
        else:
            for index, n in enumerate(string[position:],position):
                if n != 'N':
                    string = longest([string[:position], string[index:]])
                    break
                elif index == len(string) - 1:
                    string = string[:position]
                    break
    return(string)
    

def clear_from_n(fasta_in, fasta_out):
    read_fasta(fasta_in, fasta_out)
    return(0)            