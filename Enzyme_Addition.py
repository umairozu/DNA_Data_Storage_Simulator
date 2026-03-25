import os

path = fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna'

lines_1 = []
with open(path) as file:
    #cutting_site = "GATC"
    primer_F = "TGGCTCATTT"
    primer_R = "ATAAATGACC"
    for line in file:
        if line[0] != '>':
            lines_1.append(primer_F + line.strip() + primer_R)

orig_length = len(lines_1[0])
pF_length = len(primer_F)
pR_length = len(primer_R)

new_path = fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt'

if not os.path.exists(new_path): # just creating empty file
        with open(new_path,"w+") as file:
            pass

with open(new_path, "a+") as file:  #appending text
    file.writelines(lines_1)

with open(new_path,"w") as file: #appending cutting site
    for line in lines_1:
        file.writelines(line + "\n")


#Users own primer's Input taken and added to the file
"""
pattern = re.compile(r'^[ACGTacgt]{20}$')

primer_F = input("Enter Forward primer: \n").upper()
while not pattern.fullmatch(primer_F):
    primer_F = input("Enter Forward primer: \n").upper()


primer_R = input("Enter Reverse primer: \n").upper()
while not pattern.fullmatch(primer_R):
    primer_R = input("Enter Reverse primer: \n").upper()


print("Valid Inputs --> Proceed forward")


with open(new_path) as file:
    lines_2 = file.readlines()

with open(new_path,"w+") as file:
    for line in lines_2:
        line = primer_F + line.strip() + primer_R
        file.writelines(line + "\n")
"""

# final FASTA file creation
# writing final Oligo sequence [Primer_F,RS,Payload,seed,cutting_site,binding_site] back with their
# identifiers as FASTA file for later use (Debugging/ decode etc.)
with open(path) as file:
    identifiers = []
    for line in file:
        if line[0] == '>':
            identifiers.append(line)

final_path = fr'{os.getcwd()}\dna-fountain\final_file.FASTA'
if not os.path.exists(final_path):
    with open(final_path, "w") as file:
        pass

with open(final_path,"a+") as file:
    for i, j in zip(identifiers,lines_1):
        file.writelines(i)
        file.writelines(j + "\n")


if __name__ == "__main__":
    print(length_sequence)