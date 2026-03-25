
import random
import re
from typing import Final
import numpy as np
from ordered_set import OrderedSet
from Homopolymer import homopolymer


MUTATION_LIST: Final = ['insertion', 'deletion', 'substitution']

class Error_simulation:

    # user can provide their own mutation_attributes and error rates but should be of the same format
    def __init__(self, seq, process,attribute = None,error_rate = None, seed= None):
        self.bases = ['A', 'T', 'C', 'G']
        self.seq = seq
        self.process = process
        self.attributes = attribute
        self.error_rates = error_rate
        self.seed = seed if seed else np.random.seed()
        self.visited_bases = [{"base": char, "visited": False} for char in self.seq]

                                        # IMPORTANT IMPORTANT ⚠️ ⚠️ ⚠️

        # If intended to run rounds of mutation on the sequence, then reset the 'self.visited_bases' to False at each round
        # why ??
        # because for instance after some substitution in round 1, a homopolymer is introduced at the substitution position,
        # we need to allow homopolymer based mutation their for round 2.
        # keeping visited[pos] = True in this case prevents realistic behavior.



    def get_attributes(self, indels_type):
        # position is "Random" or "Homopolymer" location
        try:
            position = np.random.choice(list(indels_type["position"].keys()), p = list(indels_type["position"].values()))
        except KeyError:
            position = None

        try:
            pattern = indels_type["pattern"]
        except KeyError:
            pattern = None

        #position_range = [20,100] # this is just an example initialization
        position_range = None
        return position, pattern, position_range

    def insertion(self, ins_attr_dict = None):
        # optional insertion attribute dictionary
        if ins_attr_dict:
            position, pattern, position_range = self.get_attributes(ins_attr_dict)
        # default dictionary
        else:
            ins_dict = self.attributes['insertion']
            position, pattern, position_range = self.get_attributes(ins_dict)

        if not position or position == 'random':
            return self.indel(pattern,position_range,mode = 'insertion')

        if position == 'homopolymer':
            poly = homopolymer(self.seq)
            if poly:
                return self.indel_homopolymer(poly, pattern, mode = 'insertion')
        # if not position or position != random or position == homopolymer but poly is empty --> then random insertion
        return self.indel(pattern, position_range, mode='insertion')

    def deletion(self, del_attr_dict = None):
        # optional deletion attribute dictionary
        if del_attr_dict:
            position, pattern, position_range = self.get_attributes(del_attr_dict)
        # default dictionary
        else:
            del_dict = self.attributes['deletion']
            position, pattern, position_range = self.get_attributes(del_dict)

        if not position or position == 'random':
            return self.indel(pattern,position_range,mode = 'deletion')

        if position == 'homopolymer':
            poly = homopolymer(self.seq)
            if poly:
                return self.indel_homopolymer(poly, pattern, mode = 'deletion')
        # if not position or position != random or position == homopolymer but poly is empty --> then random insertion
        return self.indel(pattern, position_range, mode = 'deletion')


    def substitution(self, sub_attr_dict = None):
        # optional deletion attribute dictionary
        if sub_attr_dict:
            position, pattern, position_range = self.get_attributes(sub_attr_dict)
        # default dictionary
        else:
            sub_dict = self.attributes['substitution']
            position, pattern, position_range = self.get_attributes(sub_dict)

        if not position or position == 'random':
            return self.positional_sub(pattern, position_range=position_range)
        else:
            pass




    """
    If a pattern is provided e.g; {"pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}} :
    pick a target base type and delegate it to def random_indel()
    Else:
        If pattern not provided:
            check If positon range is provided:
                if so, generate a random position base ensuring a empty " " position is not retrieved,
            Otherwise: pick a random position from the entire sequence
            -delegate position to def indel_sub_base()
    """

    def indel(self, pattern, position_range, mode):
        if not pattern:
            if position_range:
                pos = random.randrange(position_range[0], position_range[1] + 1)
                while self.seq[pos] == " ": # if the chosen index is " ", chose a different index then
                    pos = random.randrange(position_range[0], position_range[1] + 1)
            else:
                pos = random.randrange(len(self.seq))
            return self.indel_sub_base(pos, mode)
        else:
            target_base = np.random.choice(list(pattern.keys()), p = list(pattern.values()))
            return self.random_indel(target_base,position_range,mode)

    """
    This method is trying to randomly find a index(with some rules ofc) that matches the target base
    - then pass that index/ pos to def indel_sub_base for mutation based on the mode ('insertion', deletion') specified
    """
    def random_indel(self,target_base, position_range, mode, count = 0):
        if position_range:
            sequence_indices = range(position_range[0], position_range[1] + 1)
        else:
            sequence_indices = range(len(self.seq))

        valid_indices = [i for i in sequence_indices if self.seq[i] != " "]

        if not valid_indices:
            return None ####

        indices = [i for i in valid_indices if self.seq[i] == target_base]

        pos = random.choice(valid_indices)

        if indices:
            pick_index = random.choice(indices)
            return self.indel_sub_base(pick_index, mode)
        else:
            if count < 12: # recursively try finding a random valid index using other bases (bases chosen randomly here also)
                target_base = random.choice(self.bases)
                return self.random_indel(target_base,position_range,mode, count = count + 1)

        # if unable to find the target base from amongst the valid sequence_indices, then mutate the
        # pos = random.choice(sequence_indices) anyways
        return self.indel_sub_base(pos,mode)



    """Insertion, deletion, substitution implementing methods"""
    def indel_sub_base(self, pos, mode):

        """2 improvements could be done later(if needed):
        - currently If position is already visited, nothing happens (silent skip),
            we can close this condition for increase mutation rate if needed
        - currently if base equals the original base, no substitution will happen, we can improve this
            by allowing any of the other 3 basses to go at that substitution place
        """


        #print(f"Original sequence: {self.seq}")


        assert mode in ("insertion", "deletion", "substitution")
        assert  0 <= pos <= len(self.seq)

        base = random.choice(self.bases)
        #print(base)
        new_mutation = {"base": base, "visited": True}

        if self.visited_bases[pos]["visited"] == False:
            print(fr"Mode: {mode.upper()}")
            if mode == 'insertion':  # pre-insertion to be specific!
                self.seq = self.seq[:pos] + base + self.seq[pos:]
                self.visited_bases.insert(pos, new_mutation)  # not touching original base (at pos) again
            elif mode == 'deletion':
                self.seq = self.seq[:pos] + " " + self.seq[pos + 1:]
                self.visited_bases[pos] = {"base": " ", "visited": True}
            else:  # substitution
                self.seq = self.seq[:pos] + base + self.seq[pos + 1:]
                self.visited_bases[pos] = new_mutation

        return self.seq # remove this return statement later (just for testing rn)

    def indel_homopolymer(self, poly, pattern, mode):
        print(f"provided polymer: {poly}")
        print(f"provided pattern: {pattern}")
        poly_bases = list(OrderedSet([item['base'] for item in poly]))

        # poly_bases = list(OrderedSet(poly['base']))
        print(f"Poly bases: {poly_bases}")
        new_pattern = []
        new_pattern_weights = {}

        if not poly:  # if no homopolymer available then go back to normal indel method
            return self.indel(pattern=None, position_range=None, mode=mode)

        for base in poly_bases:
            if base in pattern:
                new_pattern.append(base)  # this for loop is removing those bases from
                # pattern that are not in any of our homopolymers

            """        
            if not pattern: # if no pattern probabilities provided, then give equal weights to each base
            for base in new_pattern:
                new_pattern_weights[base] = 1 / len(new_pattern)
            else:
            """
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
            # instead i could make a validate pattern/ mutation_attributes method or class   #
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~⚠️IMPORTANT⚠️~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

        # sum of bases in new_pattern
        sum_bases = 0
        for base in new_pattern:
            sum_bases += pattern[base]

        """normalizing bases weight to 1"""
        # new_pattern_weights = {}
        for base in new_pattern:
            new_pattern_weights[base] = pattern[base] / sum_bases

        print(f"new pattern: {new_pattern}")
        print(f"new pattern weights: {new_pattern_weights}")

        chosen_base = np.random.choice(list(new_pattern_weights.keys()), p=list(new_pattern_weights.values()))
        print(f"chosen base from new pattern: {chosen_base}")

        # poly is like [{'base': 'A', 'chars': ['A', 'A', 'A', 'A'], 'start_pos': 0, 'end_pos': 3, 'error': 0.6}]
        # check homopolymer.py for details
        possible_mutables = [item['chars'] for item in poly if chosen_base in item['chars']]
        print(f"possible mutables: {possible_mutables}")

        chosen_mutable = random.choice(possible_mutables)  # a list here right now e.g ['T', 'T', 'T']
        chosen_mutable = "".join(chosen_mutable)  # converted to a string
        print(f"chosen mutable: {chosen_mutable}")

        # find position of the chosen mutable from the original poly list
        matching_entry = [entry for entry in poly if entry['chars'] == list(chosen_mutable)]
        if matching_entry:
            chosen_entry = random.choice(matching_entry)
            start_pos = chosen_entry['start_pos']
            print(f"chosen index:  {start_pos} ")

        """
        base = random.choice(self.bases)
        print(f"chosen base:  {base} ")

        # homopolymer is being mutated here based on mode
        if mode == 'insertion':
            chosen_mutable = chosen_mutable[:index] + base + chosen_mutable[index:]
        elif mode == 'deletion':
            chosen_mutable = chosen_mutable[:index] + " " + chosen_mutable[index + 1:]
        """

        """Output for insertion as an example:"""
        """
        provided polymer: [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']] 
        provided pattern: {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} 
        Poly bases: ['T', 'G', 'A']
        new pattern: ['T', 'G', 'A']
        new pattern weights: {'T': 0.3333333333333333, 'G': 0.3333333333333333, 'A': 0.3333333333333333}
        chosen base from new pattern: T
        possible mutables: ['TTTTT', 'TTT']
        chosen mutable: TTT
        chosen index:  0 
        chosen base:  G 
        mutated homopolymer: GTTT
        """
        return self.indel_sub_base(start_pos, mode)


    def positional_sub(self, pattern = None, position_range = None):
        if not pattern:
            return self.no_pattern_sub(position_range)
        else:
            return self.pattern_sub(pattern, position_range)

    def no_pattern_sub(self, position_range = None):
        if position_range:
            sequence_indices = range(position_range[0], position_range[1] + 1)
        else:
            sequence_indices = range(len(self.seq))

        valid_indices = [i for i in sequence_indices if self.seq[i] != " "]
        if not valid_indices:
            return None ####

        pos = random.choice(valid_indices)

        return self.indel_sub_base(pos, mode = 'substitution')




    def pattern_sub(self, pattern, position_range):
        if position_range:
            sequence_indices = range(position_range[0], position_range[1] + 1)
        else:
            sequence_indices = range(len(self.seq))

        print(sequence_indices)

        # original seq = TAGC
        # after deletion for instance it became "TA C"
        # then 'TAC' WOULD NOT match "TA C"
        # so doesn't need valid non-empty indices
        """valid_indices = [i for i in sequence_indices if self.seq[i] != " "]
        if valid_indices is None:
            return None"""

        motifs = pattern.keys()
        motifs = sorted(motifs, reverse = False)
        combined_motifs = "|".join(motifs)

        motif_matcher = re.compile(combined_motifs)  # e.g; re.compile('ATCA|ATCG|GG|T')
        print(motif_matcher)

        # Now searching the sequence region specified for any matching motifs
        # storing all matched motifs string if any with their positions
        # Then, randomly choosing one of those motifs and replacing it in the original sequence


        #matches = motif_matcher.finditer(self.seq, sequence_indices.start, sequence_indices.stop)
        matches = [{
            'base': m.group(),
            'start': m.start(),
            'end' : m.end() - 1
                    }
            for m in motif_matcher.finditer(self.seq, sequence_indices.start, sequence_indices.stop)
        ]

        if matches:
            #print(f"All matches =\n {matches}")

            chosen_match = random.choice(matches)
            chosen_base = chosen_match['base']
            chosen_span = [chosen_match['start'], chosen_match['end']]
            print(f" Chosen match = {chosen_match['base']}, Chosen Span = {chosen_span}")
            #chosen_match = chosen_match['base']

            if type(pattern[chosen_base]) == dict:
                #print("Dict as pattern")
                replacement = np.random.choice(list(pattern[chosen_base].keys()), p = list(pattern[chosen_base].values()))
            elif type(pattern[chosen_base]) == list:
                #print("List as pattern")
                replacement = np.random.choice(pattern[chosen_base])
            else:
                #print("String as pattern")
                replacement = pattern[chosen_base]

            #print(f"replacement: {replacement}")

            #Substituting here
            #print(f"Original Sequence: {self.seq}")
            self.seq = self.seq[:chosen_span[0]] + replacement + self.seq[chosen_span[1] + 1:]
            print("SUBSTITUTION")
            #logging mutated position to avoid re-mutation at the same position

            for pos, base in zip(range(chosen_span[0], chosen_span[1] + 1), replacement):
                new_mutation = {'base': base, 'visited': True}
                self.visited_bases[pos] = new_mutation

            #print(f"Substituted Sequence: {self.seq}")
            #print(self.visited_bases)

            """
            # EXAMPLE OUTPUT
            re.compile('TAC|TAG')
            All matches = [{'base': 'TAG', 'start': 11, 'end': 13}]
            Chosen match = TAG, Chosen Span = [11, 13]
            String as pattern
            replacement: TGG
            initial Sequence: ATCGAATCAGA TAG  ATAA
            after   Sequence: ATCGAATCAGA TGG  ATAA
            """

        else: # if no matched motif in the sequence, we are falling back to random single base substitution in def 'no_pattern_sub' function
            return self.no_pattern_sub(position_range)



    def run_mutations(self, mutation_list = MUTATION_LIST):
        #print(self.seq)
        for mutation_type in mutation_list:
            error_rate = self.error_rates["raw_rate"] * self.error_rates[str(mutation_type)]
            #attributes = self.get_attributes(mutation_type)
            #np.random.seed(self.seed)
            for _ in range(len(self.seq)):
                if np.random.random() <= error_rate:
                    eval('self.' + mutation_type)()
        """
        print(self.seq)
        """
        #return self.seed


    # To manually target a specific region in the sequence
    # Just a helper function
    # user input something like: [{'base': 'A', 'chars': ['A', 'A', 'A', 'A'], 'start_pos': 0, 'end_pos': 3, 'error': 0.6}]
    # and relevant things are retrieved from the input
    def manual_mutation(self, input_error):
        if np.random.random() <= input_error['error']:
            mutation_types = ['insertion', 'deletion', 'substitution']
            mutation_prob = [self.error_rates['insertion'], self.error_rates['deletion'], self.error_rates['substitution']]
            chosen_mutation = np.random.choice(mutation_types, p = mutation_prob)
            attributes = {'position_range': [input_error['start_pos'], input_error['end_pos']]}
            eval('self.' + chosen_mutation)(attributes)
        return None

    # read important notice in Constructor for reason of this method
    def reset_visited(self):
        self.visited_bases = [{"base": char, "visited": False} for char in self.seq]



if __name__ == "__main__":
    #my_poly = [['T', 'T', 'T', 'T', 'T'], ['T', 'T', 'T'], ['G', 'G', 'G'], ['A', 'A', 'A']]
    #my_pattern = {"ATCG": 0.25, "ATCA": 0.25, "GG": 0.25, "T": 0.25}

    #my_pattern = {"A": {"G": 0.5, "T": 0.25, "C": 0.25}}                     # dict
    #my_pattern = {"TAG": "TGG", "TAC": "TGC"}                                 # string
    #my_pattern = {"CG": ["CA", "TG"]}                                        # list
    my_pattern = {"CG": {"CA": 0.5, "TG": 0.5}}

    ins_mode = "insertion"
    del_mode = "deletion"
    sub_mode = "substitution"

    sequence = "ATCGAATCGGGATAGATAATCGAATCGGGATAGATA"
    #sequence = "ATCGAATCAGATAGATAA"
    my_poly2 = homopolymer(sequence)

    #sE = sequencingError(sequence, "sequencing", attribute = mutation_attributes["3"],error_rate = err_rates["3"])

    #indel_in_homopolymer = sE.indel_homopolymer(my_poly2,my_pattern,del_mode)
    #print(indel_in_homopolymer)

    #print(sE.pattern_sub(my_pattern,None))


