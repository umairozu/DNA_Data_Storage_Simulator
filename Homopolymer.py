# Error rates Sourced from MESA
# https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/error_sources/homopolymers.py

def error_func(homopolymer_length, base= None):
    if homopolymer_length < 3:
        return 0.0
    elif homopolymer_length < 4:
        return 0.3
    elif homopolymer_length < 5:
        return 0.6
    elif homopolymer_length < 6:
        return 0.9
    else:
        return 1.0

# take a sequence and error_func
# return list of tuple of (Base,homopolymer_error probability)
def homopolymer(sequence, include_chars=True):
    max_homopolymers = []
    sequence_length = len(sequence)
    start_pos = 0

    while start_pos < sequence_length:
        base = sequence[start_pos]
        end_pos = start_pos + 1
        while end_pos < sequence_length and sequence[end_pos] == base:
            end_pos += 1

        length = end_pos - start_pos
        error = error_func(length)
        if error > 0.0 and base != " ":
            entry = {'base': base}
            if include_chars:
                entry['chars'] = [base] * length
            entry.update({
                'start_pos': start_pos,
                'end_pos': end_pos - 1,
                'error': error
            })
            max_homopolymers.append(entry) #output currently: [{'base': 'A', 'chars': ['A', 'A', 'A', 'A'], 'start_pos': 0, 'end_pos': 3, 'error': 0.6}]
        start_pos = end_pos

    """
    for seq in max_homopolymers:
        length = len(seq)
        error = error_func(length)
        if error > 0.0:
            #result.append(("".join(seq),error))
            result.append(seq) # output currently: [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']]
    """
    return max_homopolymers

if __name__ == "__main__":
    print(homopolymer("AAG   TCA AAA"))


