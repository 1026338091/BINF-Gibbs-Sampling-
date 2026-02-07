# Introduction
Description of the project

# Pseudocode
What we know about the algorithm:
Initialization:
    choose a random position in each sequence and use that for the initial motif
    calculate the frequency of each nucleotide for each position within that motif = pfm
    calculate the frequency of every nucleotide in non-motif positions as background (for now we ignore it)
    convert frequency matrix to the probability weigth matrix

Iteration:
    choosing one of the sequences at random
    recalculate the pwm wihout the chosen sequence
    pwm gets used to score each motif in the removed sequence
    we choose the new motif that's going to get put back into the pool

Convergence:
    check the pfm from last iteration against current iteration, if they are similar enough then we've reached convergence
    return pfm
    


```python
bam_path = "assignment3/SRR9090854.subsampled_5pct.bam"
seqs = [read.seq for read in bs.AlignmentFile(bam_path)]

# consider the seeding because we might get stuck in local minima and never get out
def GibbsMotifFinder(seqs, k, seed=None): 
    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    for i in range(100000):
        
    
        check after 100 iterations to see if the comparison meets the threshold:
            if the dfiference between the scores is low enough (0.1)
            return the pfm




    return pfm



# we need a list for the motifs
def choose_motifs(seqs, k (length of the motif) - 10, rng):
    motif_list = []
    for seq in seqs:
        pick a random int between 0 and len(seq) - k
        choose motif from seq[int: int + k]
        append to motif_list
    choose a sequence at random from the motif_list
    return motif_list, motif

use this list to find the pfm and pwm
build_pfm(sequences: List[str], length: int)
build_pwm(pfm: np.ndarray)

next we need to pick the sequence to remove from the list

def random_seq(motifs, rng):
    choose a sequence at random - motif_list
    return motif

# get scores for all motifs in the sequence we removed

    
    

```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
