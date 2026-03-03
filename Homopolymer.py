from itertools import groupby

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
def homopolymer(sequence):
    result = []

    # k (The Key): This is the character that is currently repeating (e.g., 'A')
    # g (The Group): This is an iterator containing all the items in that specific streak.
    # Because it is an iterator,you usually have to convert it to a list(g) to see the contents.

    "max_homopolymers = [list(g) for k, g in groupby(sequence)] # just collecting all homopolymers here"
    start_pos = 0
    max_homopolymers = []
    for k, g in groupby(sequence):
        group_list = list(g)
        length = len(group_list)
        end_pos = start_pos + length - 1

        #print(group_list)
        error = error_func(length)
        if error > 0.0 and group_list[0] != " ":
            #print(group_list)
            max_homopolymers.append({
                'base': k, 'chars': group_list, 'start_pos': start_pos, 'end_pos': end_pos, 'error': error
            }) #output currently: [{'base': 'A', 'chars': ['A', 'A', 'A', 'A'], 'start_pos': 0, 'end_pos': 3, 'error': 0.6}]
        start_pos += length

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




