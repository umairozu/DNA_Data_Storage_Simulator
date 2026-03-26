
from dataclasses import dataclass
import random
from typing import List, Tuple
from collections import Counter
from Arrhenius_decay import Arrhenius_decay


xlsx = "RawData.xlsx"
sim = Arrhenius_decay.from_xlsx(xlsx)

@ dataclass
class Demo2:
    frac = 0.0
    prev_fracs = {}

    @staticmethod
    def copyFunc(sequence: str, num_copy: int, master_seed: int) -> List[Tuple[str, int]]:
        copyList = []
        for i in range(num_copy):
            seed = master_seed + i
            # res = ''.join(random.choices(string.digits,k=4))
            copyList.append((sequence, seed))
        return copyList

    @staticmethod
    def mutate(sub_list: List[Tuple[str, float]], eps: float = 1e-3) -> List[Tuple[str, float]]:
        # print(f"random = {rand}, (char, frac) = ({char},{frac})")

        for i, (char, frac) in enumerate(sub_list):
            if frac <= eps:
                sub_list[i] = ('X', 0.0)
            else:
                sub_list[i] = (char, frac)
            # print(sub_list[i])
        return sub_list

    @staticmethod
    def per_base_decay(seqs: str, f: float, copy_key: int, step_seed: int) -> List[Tuple[Tuple[str, float]]]:
        out: List[Tuple[Tuple[str, float]]] = []

        if copy_key not in Demo2.prev_fracs or len(Demo2.prev_fracs[copy_key]) != len(seqs):
            Demo2.prev_fracs[copy_key] = [1.0] * len(seqs)

        #
        rng = random.Random(step_seed)

        for i, ch in enumerate(seqs):
            index = rng.randrange(len(seqs))
            old_frac = Demo2.prev_fracs[copy_key][i]

            #### Making sure most of the bases are affected!!
            if index >= rng.randrange(len(seqs)):
                Demo2.prev_fracs[copy_key][i] = f
                out.append((ch, f))
            else:
                out.append((ch, old_frac))

        return out

    @staticmethod
    def get_output_base(data):
        seqs = [item for item in data]
        majority = []

        for chars in zip(*seqs):
            counts = Counter(chars)
            max_count = max(counts.values())

            winners = []
            for char, count in counts.items():
                winners.append(char)

            if len(winners) > 1 and ('X' or 'A' or 'G' or 'C' or 'T') in winners:
                winners.remove('X')
            majority.append(winners[0])
        return "".join(majority)


    @staticmethod
    def recovery_ratio(seq: str) -> float:
        seq_len = len(seq)
        num_mutate = 0
        for i, char in enumerate(seq):
            if char == 'X':
                num_mutate += 1
        ratio = ((seq_len - num_mutate) / seq_len) * 100
        return ratio


myList = []


def main():
    global myList

    RAID1 = 1
    RAID2 = 2
    RAID3 = 3
    RAID4 = 4
    RAID5 = 5

    with open("1K_sequence.txt", 'r') as file:
        sequences = file.read()

    copy_List = Demo2.copyFunc("ACGTAATTGGCCTCCAAGTC", 1 + RAID5  , 80)
    #copy_List = Demo2.copyFunc(sequences, 1 + RAID5, 2)
    print(f"Sequence/s, Seed: {copy_List}")
    temps_low = [4]
    temps_med = [25]
    temps_high = [55]

    for seq, master_seed in copy_List:
        print("-------------------------------------")
        copy_states = []
        for i, T in enumerate(temps_med):
            Demo2.frac = Arrhenius_decay.remaining_dna_frac(sim, T, True, 4000)
            state = Demo2.per_base_decay(seq, Demo2.frac, copy_key=master_seed,
                                         step_seed=master_seed + i)
            print(state, f"....degradation at Temp {T}")
            print("")
            state = Demo2.mutate(state)
            copy_states.append(state)
        myList.append(copy_states)


    majority_list = []
    print(f"Degradation @ Temp {temps_med[0]} for 4000 weeks")
    for i, copy_states in enumerate(myList):
        final_state = copy_states[-1]
        final_seq = "".join(ch for ch, _ in final_state)
        print(f"Copy {i} final: {final_seq}")
        majority_list.append(final_seq)
        # print(f"Majority seq: {majority_list}")

    result = Demo2.get_output_base(majority_list)
    print(f"final recovered sequence: {result}")
    print(f"Recovery Ratio: {Demo2.recovery_ratio(result)}")


if __name__ == "__main__":
    main()
