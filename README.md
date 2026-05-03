# DNA Data Storage Simulator 🧬
A simulation environment implementing concepts for long-term data storage using DNA sequences.

## Project Goal

The primary goal of this repository is to provide a flexible simulation environment for experimenting with DNA-based storage pipelines, translating theoretical ideas from literature into practical, testable code. 

**Specifically, this project focuses on:**

* **Exploring** how errors propagate across different stages of the pipeline.
* **Evaluating** how sequence constraints affect overall system reliability.
* **Experimenting** with various storage-decay assumptions and models.

---

## Conceptual Workflow


```text
Input file
   ↓
Compress to .tar.gz
   ↓
DNA Fountain encoding
   ↓
Primer / enzyme addition
   ↓
Synthesis simulation
   ↓
Storage simulation
   ↓
PCR simulation
   ↓
Sequencing simulation
   ↓
Cutadapt trimming
   ↓
PEAR or BWA-MEM2 processing
   ↓
DNA Fountain decoding
   ↓
Compare decoded file with original input
```
---

## 1. System Requirements

Recommended:

```text
Ubuntu Linux
Python 3.10+
Miniforge / Conda
gcc / build-essential
```

Install basic system packages:

```bash
sudo apt update
sudo apt upgrade -y

sudo apt install -y build-essential python3 python3-dev python3-venv python3-pip git wget curl
```
---

## 2. Install Miniforge and Required tools

If Conda/Miniforge is not installed:

```bash
cd ~
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O Miniforge3.sh
bash Miniforge3.sh
```
---


## 3. Create Bioinformatics Environment

This environment is used for:

```text
Cutadapt
PEAR
BWA-MEM2
Samtools
```

Create and activate the environment:

```bash
source ~/miniforge3/etc/profile.d/conda.sh

conda create -n dna_tools -c conda-forge -c bioconda python=3.11 cutadapt pear bwa-mem2 samtools -y

conda activate dna_tools
```
---

## 4. Clone and Set Up DNA Fountain

This project uses DNA Fountain for encoding and decoding. Clone the DNA Fountain repository inside the project root:

```bash
git clone https://github.com/jdbrody/dna-fountain.git dna-fountain
```

Now create and activate a Python virtual environment:

```bash
python3 -m venv venv_sim
source venv_sim/bin/activate
```

Install the required Python packages:

```bash
python -m pip install --upgrade pip setuptools wheel

python -m pip install numpy scipy cython reedsolo tqdm
```

Build the DNA Fountain Cython modules:

```bash
cd dna-fountain
python setup.py build_ext --inplace
```

Test that the decoder is working:

```bash
python decode.py -h
```

If the help message appears, DNA Fountain is ready.

## 5. Basic Simulation Workflow

### Step 1: Compress input file

Example:

```bash
tar -czvf test_file.tar.gz test_file.txt
```
---

### Step 2: Encode with DNA Fountain

```bash
cd dna-fountain

python encode.py --file_in test_file.tar.gz --size 16 -m 3 --gc 0.05 --rs 5 --delta 0.05 --out test_file.tar.gz.dna --stop 1000
```

Common parameters:

```text
--size 16     payload size in bytes
-m 3          maximum homopolymer length
--gc 0.05     GC content constraint, around 50 ± 5%
--rs 5        Reed-Solomon error correction bytes
--delta 0.05  DNA Fountain distribution parameter
--stop 1000   number of encoded oligos to generate
```

---

### Step 3: Add primers / enzyme region

Default:

```bash
python3 Enzyme_Addition.py
```

Example with custom primers:

```bash
python3 Enzyme_Addition.py --pf "FORWARD_PRIMER" --pr "REVERSE_PRIMER" --gc_min 0.45 --gc_max 0.55
```

This creates an orderable oligo file such as:

```text
test_file.tar.gz.dna_order.txt
```
---

### Step 4: Run synthesis simulation

```bash
python3 synthesis.py --mut 1 --in_file dna-fountain/test_file.tar.gz.dna_order.txt --out_file synthesis_LOW
```

Mutation modes:

```text
--mut 0    no mutation
--mut 1    low mutation
--mut 2    medium mutation
--mut 3    high mutation
--c        custom mutation
```

---

### Step 5: Run storage simulation

```bash
python3 storage.py --temp 30 --ph 7 --week 4 --encap 1 --mut 1 --in_file synthesis/synthesis_LOW.txt --out_file storage_LOW
```

---

### Step 6: Run PCR simulation

```bash
python3 pcr.py --s 0.5 --n 30 --mut 1 --in_file storage/storage_LOW.txt --out_file pcr_LOW
```

Parameters:

```text
--s      sampling fraction
--n      number of PCR cycles
--mut    PCR mutation intensity
```

---

### Step 7: Run sequencing simulation

Paired-end example:

```bash
python3 sequencing.py --type 1 --m 2 --s 1 --t 20000 --rl 100 --mut 1 --in_file pcr/pcr_LOW.txt --order_file dna-fountain/test_file.tar.gz.dna_order.txt
```

Parameters:

```text
--type 1    Illumina sequencing
--m 1       single-end sequencing
--m 2       paired-end sequencing
--s         sampling fraction
--t         total output reads
--rl        read length
--mut       sequencing mutation intensity
```

Sequencing output:

```text
sequencing/sequencing_R1.fastq
sequencing/sequencing_R2.fastq
```

---

## 6. Decoding Route 1: Cutadapt -> PEAR -> DNA Fountain Decoder

Use this route when paired-end reads can be merged into full-length sequences.

Run:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

This script performs:

```text
sequencing_R1.fastq + sequencing_R2.fastq
        ↓
Cutadapt primer/adapter trimming
        ↓
PEAR paired-end merging
        ↓
Filter full-length sequences
        ↓
Count and sort unique sequences
        ↓
Create decoder_input.seq
        ↓
DNA Fountain decode.py
```
---

## 7. Decoding Route 2: Cutadapt -> BWA-MEM2 -> DNA Fountain Decoder

Use this route for ## **'Reference-based Oligo reconstruction'**.

Run:

```bash
chmod +x run_pipeline_bwa.sh
./run_pipeline_bwa.sh
```

This script performs:

```text
sequencing_R1.fastq / sequencing_R2.fastq
        ↓
Cutadapt primer/adapter trimming
        ↓
Build BWA reference from designed oligos
        ↓
BWA-MEM2 alignment
        ↓
Convert aln.sam to decoder_input.seq
        ↓
DNA Fountain decode.py
```

BWA-MEM2 does not decode the file directly. It aligns sequencing reads to the expected oligo library. The pipeline then converts the SAM alignment output into a sequence-only file for DNA Fountain.

The DNA Fountain decoder expects:

```text
ACGT...
TGCA...
GATT...
```

one DNA sequence per line.

---

## Before Running the Pipeline

Before running `run_pipeline.sh` or `run_pipeline_bwa.sh`, replace the example values in the scripts with your own experiment-specific values.

Update these items:
- Cutadapt adapter/primer sequences:
  ```text
  -a YOUR_R1 3 PRIME_ADAPTERX
  -A YOUR_R2 3 PRIME_ADAPTERX
  ```

- Expected merged length in run_pipeline.sh:  
  ```text
   length($0)==100
  ```

- Primer/flank lengths in run_pipeline_bwa.sh:
  ```text
   PRIMER_F_LEN = 20
   PRIMER_R_LEN = 20
  ```
  

### 📌 Notes

> This repository is a **research-oriented simulator** and experimental codebase. 

It is designed as a platform for testing concepts derived from DNA storage literature rather than a finalized production framework. Use this environment to prototype, iterate, and validate theoretical models in a practical setting.

## Referenced Papers & Research

This simulator is built upon the methodologies and ideas described in the following research:

### 1. A Compact Cassette Tape for DNA-Based Data Storage

* **Published:** 2025
* **Journal:** *Science Advances*
* **Authors:** Jiankai Li, Cuiping Mao, Shuchen Wang, Xingjian Li, Xueqing Luo, Dou Wang, Shuo Zheng, Jialin Shao, Rui Wang, Chunhai Fan, Xingyu Jiang
* **Why it matters here:** Inspires the system-level view of DNA storage as a practical storage device with addressing, partitioning, recovery, redeposition, and long-term preservation.
* **Links:**
  * [Read the Full Paper](https://www.science.org/doi/10.1126/sciadv.ady3406)
  * [Original Implementation (GitHub)](https://github.com/JianKai-Lee/DNA-Cassette-Tape/tree/main/Editing_analysis)

---

### 2. DeSP: A Systematic DNA Storage Error Simulation Pipeline

* **Published:** 2022
* **Journal:** *BMC Bioinformatics*
* **Authors:** Lekang Yuan, Zhen Xie, Ye Wang, Xiaowo Wang
* **Why it matters here:** Provides the end-to-end simulation perspective across synthesis, decay, PCR, sampling, and sequencing, while emphasizing both sequence loss and within-sequence errors.
* **Links:**
  * [Read the Full Paper](https://doi.org/10.1186/s12859-022-04723-w)
  * [Official DeSP Repository](https://github.com/WangLabTHU/DeSP)

---

### 3. MESA: Automated Assessment of Synthetic DNA Fragments and Simulation of DNA Synthesis, Storage, Sequencing, and PCR Errors

* **Published:** 2020
* **Journal:** *Bioinformatics*
* **Authors:** Michael Schwarz, Marius Welzel, Tolganay Kabdullayeva, Anke Becker, Bernd Freisleben, Dominik Heider
* **Why it matters here:** Inspires sequence validation and constraint-aware simulation using GC content, homopolymers, repeats, motifs, and configurable error models.
* **Links:**
  * [Read the Full Paper](https://doi.org/10.1093/bioinformatics/btaa140)
  * [Official MESA Repository](https://github.com/umr-ds/mesa_dna_sim)

---

### 4. DNA Fountain Enables a Robust and Efficient Storage Architecture

* **Published:** 2017
* **Journal:** *Science*
* **Authors:** Yaniv Erlich, Dina Zielinski
* **Why it matters here:** Establishes a robust coding architecture for DNA storage under biochemical constraints such as GC balance, homopolymers, oligo dropout, and redundancy.
* **Links:**
  * [Read the Paper](https://www.science.org/doi/10.1126/science.aaj2038)

---

### 5. Dna-Storalator: A Computational Simulator for DNA Data Storage

* **Published:** 2025
* **Journal:** *BMC Bioinformatics*
* **Authors:** Gadi Chaykin, Omer Sabary, Nili Furman, Dvir Ben Shabat, Eitan Yaakobi
* **Why it matters here:** Extends the simulator viewpoint into clustering, reconstruction, error characterization, and benchmarking of retrieval pipelines.
* **Links:**
  * [Read the Full Paper](https://doi.org/10.1186/s12859-025-06222-0)
  * [DNA-Storalator Repository](https://github.com/gadihh/DNAStoralator)

---
