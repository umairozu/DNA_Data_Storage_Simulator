
import random
import re
from typing import Final
import numpy as np
from Homopolymer import homopolymer


MUTATION_LIST: Final = ['insertion', 'deletion', 'substitution']
BASES: Final = ('A', 'T', 'C', 'G')
BASE_ORDS: Final = tuple(ord(base) for base in BASES)
SPACE_ORD: Final = ord(" ")
_CHOICE_CACHE = {}
_REGEX_CACHE = {}

class Error_simulation:

    # user can provide their own mutation_attributes and error rates but should be of the same format
    def __init__(self, seq, process,attribute = None,error_rate = None, seed= None):
        self.bases = BASES
        self._base_ord = BASE_ORDS
        self.seq = seq
        self.process = process
        self.attributes = attribute
        self.error_rates = error_rate
        self.seed = seed

                                        # IMPORTANT IMPORTANT ⚠️ ⚠️ ⚠️

        # If intended to run rounds of mutation on the sequence, then reset the 'self.visited_bases' to False at each round
        # why ??
        # because for instance after some substitution in round 1, a homopolymer is introduced at the substitution position,
        # we need to allow homopolymer based mutation their for round 2.
        # keeping visited[pos] = True in this case prevents realistic behavior.


    @property
    def seq(self):
        if self._seq_cache is None:
            self._seq_cache = self._seq.decode("ascii")
        return self._seq_cache

    @seq.setter
    def seq(self, value):
        if isinstance(value, bytearray):
            self._seq = value
            self._seq_cache = value.decode("ascii")
        elif isinstance(value, bytes):
            self._seq = bytearray(value)
            self._seq_cache = value.decode("ascii")
        else:
            self._seq_cache = str(value)
            self._seq = bytearray(self._seq_cache, "ascii")
        self._visited = bytearray(len(self._seq))

    @property
    def visited_bases(self):
        return [
            {"base": chr(base), "visited": bool(visited)}
            for base, visited in zip(self._seq, self._visited)
        ]

    @visited_bases.setter
    def visited_bases(self, value):
        visited = bytearray(1 if item.get("visited") else 0 for item in value)
        if len(visited) < len(self._seq):
            visited.extend(b"\x00" * (len(self._seq) - len(visited)))
        self._visited = visited[:len(self._seq)]

    def _invalidate_seq_cache(self):
        self._seq_cache = None

    def _range_bounds(self, position_range):
        if position_range:
            return position_range[0], position_range[1] + 1
        return 0, len(self._seq)

    def _weighted_choice(self, pattern):
        pattern_id = id(pattern)
        cached = _CHOICE_CACHE.get(pattern_id)
        if cached is None or cached[0] is not pattern:
            cached = (pattern, tuple(pattern.keys()), tuple(pattern.values()))
            _CHOICE_CACHE[pattern_id] = cached
        _, keys, weights = cached
        return np.random.choice(keys, p=weights)

    def _base_counts(self, start, stop):
        counts = {base: 0 for base in self.bases}
        valid_count = 0
        for pos in range(start, stop):
            base = self._seq[pos]
            if base == SPACE_ORD:
                continue
            valid_count += 1
            base_char = chr(base)
            if base_char in counts:
                counts[base_char] += 1
        return counts, valid_count

    def _count_valid(self, start, stop):
        count = 0
        for pos in range(start, stop):
            if self._seq[pos] != SPACE_ORD:
                count += 1
        return count

    def _find_nth_valid(self, start, stop, target_offset):
        seen = 0
        for pos in range(start, stop):
            if self._seq[pos] == SPACE_ORD:
                continue
            if seen == target_offset:
                return pos
            seen += 1
        return None

    def _find_nth_base(self, start, stop, target_base, target_offset):
        target_ord = ord(target_base)
        seen = 0
        for pos in range(start, stop):
            if self._seq[pos] != target_ord:
                continue
            if seen == target_offset:
                return pos
            seen += 1
        return None

    def _choose_valid_index(self, position_range):
        start, stop = self._range_bounds(position_range)
        valid_count = self._count_valid(start, stop)
        if not valid_count:
            return None
        return self._find_nth_valid(start, stop, random.randrange(valid_count))


    def get_attributes(self, indels_type):
        # position is "Random" or "Homopolymer" location
        try:
            position = self._weighted_choice(indels_type["position"])
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
            poly = homopolymer(self.seq, include_chars=False)
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
            poly = homopolymer(self.seq, include_chars=False)
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
                pos = self._choose_valid_index(position_range)
                if pos is None:
                    return None
            else:
                pos = random.randrange(len(self._seq))
            return self.indel_sub_base(pos, mode)
        else:
            target_base = self._weighted_choice(pattern)
            return self.random_indel(target_base,position_range,mode)

    """
    This method is trying to randomly find a index(with some rules ofc) that matches the target base
    - then pass that index/ pos to def indel_sub_base for mutation based on the mode ('insertion', deletion') specified
    """
    def random_indel(self,target_base, position_range, mode, count = 0):
        start, stop = self._range_bounds(position_range)
        base_counts, valid_count = self._base_counts(start, stop)

        if not valid_count:
            return None ####

        while True:
            target_count = base_counts.get(target_base, 0)
            if target_count:
                pick_index = self._find_nth_base(start, stop, target_base, random.randrange(target_count))
                return self.indel_sub_base(pick_index, mode)
            if count < 12: # recursively try finding a random valid index using other bases (bases chosen randomly here also)
                target_base = random.choice(self.bases)
                count += 1
                continue
            break

        # if unable to find the target base from amongst the valid sequence_indices, then mutate the
        # pos = random.choice(sequence_indices) anyways
        pos = self._find_nth_valid(start, stop, random.randrange(valid_count))
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
        assert  0 <= pos <= len(self._seq)

        base = random.choice(self._base_ord)
        #print(base)

        if self._visited[pos] == 0:
            #print(fr"Mode: {mode.upper()}")
            if mode == 'insertion':  # pre-insertion to be specific!
                self._seq.insert(pos, base)
                self._visited.insert(pos, 1)  # not touching original base (at pos) again
            elif mode == 'deletion':
                self._seq[pos] = SPACE_ORD
                self._visited[pos] = 1
            else:  # substitution
                self._seq[pos] = base
                self._visited[pos] = 1
            self._invalidate_seq_cache()

        return self.seq # remove this return statement later (just for testing rn)

    def indel_homopolymer(self, poly, pattern, mode):
        seen_bases = set()
        poly_bases = []
        for item in poly:
            base = item['base']
            if base not in seen_bases:
                seen_bases.add(base)
                poly_bases.append(base)
        new_pattern = []
        new_pattern_weights = {}

        if not poly:  # if no homopolymer available then go back to normal indel method
            return self.indel(pattern=None, position_range=None, mode=mode)

        for base in poly_bases:
            if pattern and base in pattern:
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

        if not sum_bases:
            return self.indel(pattern=None, position_range=None, mode=mode)

        """normalizing bases weight to 1"""
        # new_pattern_weights = {}
        for base in new_pattern:
            new_pattern_weights[base] = pattern[base] / sum_bases

        chosen_base = np.random.choice(list(new_pattern_weights.keys()), p=list(new_pattern_weights.values()))

        # poly is like [{'base': 'A', 'chars': ['A', 'A', 'A', 'A'], 'start_pos': 0, 'end_pos': 3, 'error': 0.6}]
        # check homopolymer.py for details
        possible_mutables_count = 0
        for item in poly:
            if item['base'] == chosen_base:
                possible_mutables_count += 1

        # find position of the chosen mutable from the original poly list
        chosen_mutable_offset = random.randrange(possible_mutables_count)
        seen_mutables = 0
        start_pos = None
        for item in poly:
            if item['base'] != chosen_base:
                continue
            if seen_mutables == chosen_mutable_offset:
                start_pos = item['start_pos']
                break
            seen_mutables += 1

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
        start, stop = self._range_bounds(position_range)
        valid_count = self._count_valid(start, stop)
        if not valid_count:
            return None ####

        pos = self._find_nth_valid(start, stop, random.randrange(valid_count))

        return self.indel_sub_base(pos, mode = 'substitution')




    def pattern_sub(self, pattern, position_range):
        start, stop = self._range_bounds(position_range)

        #print(sequence_indices)

        # original seq = TAGC
        # after deletion for instance it became "TA C"
        # then 'TAC' WOULD NOT match "TA C"
        # so doesn't need valid non-empty indices
        """valid_indices = [i for i in sequence_indices if self.seq[i] != " "]
        if valid_indices is None:
            return None"""

        motifs = tuple(sorted(pattern.keys(), reverse = False))
        motif_matcher = _REGEX_CACHE.get(motifs)
        if motif_matcher is None:
            combined_motifs = "|".join(motifs)
            motif_matcher = re.compile(combined_motifs)  # e.g; re.compile('ATCA|ATCG|GG|T')
            _REGEX_CACHE[motifs] = motif_matcher
        #print(motif_matcher)

        # Now searching the sequence region specified for any matching motifs
        # storing all matched motifs string if any with their positions
        # Then, randomly choosing one of those motifs and replacing it in the original sequence


        seq = self.seq
        match_count = 0
        for _ in motif_matcher.finditer(seq, start, stop):
            match_count += 1

        chosen_match = None
        if match_count:
            chosen_match_offset = random.randrange(match_count)
            for current_offset, match in enumerate(motif_matcher.finditer(seq, start, stop)):
                if current_offset == chosen_match_offset:
                    chosen_match = match
                    break

        if chosen_match:
            #print(f"All matches =\n {matches}")

            chosen_base = chosen_match.group()
            chosen_span = [chosen_match.start(), chosen_match.end() - 1]
            #print(f" Chosen match = {chosen_match['base']}, Chosen Span = {chosen_span}")
            #chosen_match = chosen_match['base']

            if type(pattern[chosen_base]) == dict:
                #print("Dict as pattern")
                replacement = self._weighted_choice(pattern[chosen_base])
            elif type(pattern[chosen_base]) == list:
                #print("List as pattern")
                replacement = np.random.choice(pattern[chosen_base])
            else:
                #print("String as pattern")
                replacement = pattern[chosen_base]

            #print(f"replacement: {replacement}")

            #Substituting here
            #print(f"Original Sequence: {self.seq}")
            replacement_bytes = bytearray(str(replacement), "ascii")
            self._seq[chosen_span[0]:chosen_span[1] + 1] = replacement_bytes
            #print("SUBSTITUTION")
            #logging mutated position to avoid re-mutation at the same position

            self._visited[chosen_span[0]:chosen_span[1] + 1] = bytearray([1]) * len(replacement_bytes)
            self._invalidate_seq_cache()

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
        handlers = {
            'insertion': self.insertion,
            'deletion': self.deletion,
            'substitution': self.substitution,
        }
        for mutation_type in mutation_list:
            mutation_type = str(mutation_type)
            error_rate = self.error_rates["raw_rate"] * self.error_rates[mutation_type]
            if not (error_rate > 0):
                continue
            mutation_count = len(self._seq) if error_rate >= 1 else np.random.binomial(len(self._seq), error_rate)
            handler = handlers[mutation_type]
            #attributes = self.get_attributes(mutation_type)
            #np.random.seed(self.seed)
            for _ in range(mutation_count):
                handler()
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
            {
                'insertion': self.insertion,
                'deletion': self.deletion,
                'substitution': self.substitution,
            }[chosen_mutation](attributes)
        return None

    # read important notice in Constructor for reason of this method
    def reset_visited(self):
        self._visited = bytearray(len(self._seq))



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


