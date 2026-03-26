# DNA Data Storage Simulator 🧬
A simulation environment implementing concepts for long-term data storage using DNA sequences.

## Project Goal

The primary goal of this repository is to provide a flexible simulation environment for experimenting with DNA-based storage pipelines, translating theoretical ideas from literature into practical, testable code. 

**Specifically, this project focuses on:**

* **Exploring** how errors propagate across different stages of the pipeline.
* **Evaluating** how sequence constraints affect overall system reliability.
* **Experimenting** with various storage-decay assumptions and models.
* **Laying the groundwork** for future studies in encoding/decoding, and sequence reconstruction.

---

## Conceptual Workflow

```text
Digital Input
   ↓
Encoding / DNA Representation
   ↓
Synthesis Simulation
   ↓
Storage / Decay Modeling
   ↓
PCR Amplification
   ↓
Sequencing Simulation
   ↓
Constraint Analysis / Error Analysis
   ↓
Downstream Decoding or Benchmarking
```
---

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

## What This Repository Includes

This repository currently contains code for several parts of a DNA storage simulation workflow:

* **Synthesis simulation** — `synthesis.py`
* **Storage and decay modeling** — `storage.py`, `Arrhenius_decay.py`
* **PCR simulation** — `pcr.py`
* **Sequencing simulation** — `sequencing.py`
* **Sequence quality checks**
  * `GC_content.py`
  * `Homopolymer.py`
  * `K_mer.py`
  * `Undesired_sequences.py`
* **DNA Fountain-related assets** — `dna-fountain/`

Acknowledgment
#TODO
---
