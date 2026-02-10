import math
import random


import numpy as np
import bamnostic as bs
import seqlogo as sl

from assignment3.seq_ops import reverse_complement
#import function for building sequence motif & idenfitying seqs matching to motif
from assignment3.data_readers import *
from assignment3.seq_ops import get_seq
from assignment3.motif_ops import *


# here is our bam file and our sequences loaded in using bamnostic
bam_path = "assignment3/SRR9090854.subsampled_5pct.bam"
seqs = [read.seq for read in bs.AlignmentFile(bam_path)]
clean_seqs = []
for seq in seqs:
    if 'N' not in seq:
        clean_seqs.append(seq)

# we're going to want to develop our program with a random subsample from our seqs, so let's do that first
subseqs = random.sample(clean_seqs, 100000)


# this is our function that will find the binding site motif for p53, but we'll make helper functions for it to so we can break our code up
def GibbsMotifFinder(seqs, k, seed=None):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler

    Args:
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find
        seed (int, default=None): seed for np.random

    Returns:
        pfm (numpy array): dimensions are 4xlength
    '''
    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    # first we'll want to choose our random motifs from the list for each seq (motifs = choose_motifs(seqs, k, rng))
    motifs = choose_motifs(seqs, k, rng)
    # calculate a pfm with pfm_builder (pfm = pfm_builder(seqs, length))
    pfm = build_pfm(motifs, k)
    print(pfm)

    # calculate a pwm with the pfm and pwm_builder (new_pwm = build_pwm(pmf))
    pwm = build_pwm(pfm)
    # where this needs to be an iterative process of unknown duration we could pick a stopping point like 10000 (i.e. if we're still going after 10,000 iterations let's stop)
    # could do a for i in range(x):
    for i in range(10001):
        moving_list = []
        # pick a random motif from the motifs (motifs = random_seq(motifs, rng))
        motif, random_int = random_seq(motifs, rng)
        # remove the sequence from both the motif list and the seqs list
        for j in range(len(motifs)):
            if j != random_int:
                moving_list.append(motifs[j])

        pfm = build_pfm(moving_list, 10)
        pwm - build_pwm(pfm)


        # every time we've iterated through say, 100 times:
        if i % 100 == 0 and i > 0:

            # compare the information content of the new pfm to the old pfm (old pfm should be one iteration behind):
            old_ifc = pfm_ic(pfm)
            pfm = build_pfm(moving_list, 10)
            pwm = build_pwm(pfm)
            new_ifc = pfm_ic(pfm)

            if math.isclose(new_ifc, old_ifc) == True:
                return pfm


        seq = seqs[random_int]
        # get the reverse_strand of the left out sequence with reverse_complement - use index created above to choose seq
        rev_seq = reverse_complement(seq)
        # get scores of all the motifs of the left out sequence with get_all_scores (use seqs[i] and the reverse strand, using the index from above)
        scores, seq_motifs = get_all_scores(seq, rev_seq, k, pwm)
        # choose the motif with select_motif
        new_motif = select_motif(scores, seq_motifs)
        motifs[random_int] = new_motif[0]
        # if the motif came from the reverse strand, replace it in the full seq pool by the random int value
        if seq_motifs.index(new_motif[0]) % 2 == 1:
            seqs[random_int] = rev_seq



    # if we reached the end of the loop, return the pfm
    return pfm


def choose_motifs(seqs, k, rng):
    '''
    Args:
        seqs (str list): a list of sequences
        k (int): the length of the motif to find
        rng: the random number generator to use to get motif positions
    :return: a list of the selected motifs from each position
    '''
    # initialize a motif list
    motif_list = []
    # for seq in seqs
    for seq in seqs:
        # pick a random integer between 0 and len(seq) - len(motif)
        random_int = rng.integers(0, len(seq) - k)
        # append the motif list with seq[int: int + k]
        motif = seq[random_int: random_int + k]
        motif_list.append(motif)
    # return the motif list
    return motif_list


def random_seq(motifs, rng):
    '''

    :param motifs: a list of motifs
    param rng: random number generator to select the sequence to remove
    :return: motif list without the selected sequence
    '''

    # choose a sequence at random (motif list and seq list should line up 1 to 1 so we can do this on the motif list)
    random_int = rng.integers(0, len(motifs) - 1)
    # return motif to remove so that it can be used to recalculate the pwm with calc pwm
    return motifs[random_int], random_int


def get_all_scores(seq, rev_seq, k, pwm):
    '''

    :param seq: sequence to score
    :param rev_seq: reverse compliment of the seq to score
    :param k: motif length
    :param pwm: probability weight matrix to score the sequence against
    :return: a list of scores for each motif in the forward and reverse sequence
    '''

    # initialize a score list for forward and reverse scores
    score_list = []
    seq_motifs = []

    # for each motif in the forward seq:
    for i in range(len(seq) - k + 1):
        f_motif = seq[i: i + k]
        r_motif = rev_seq[i: i + k]
        if len (r_motif) < 10:
            print(i)
            print(r_motif)


        seq_motifs.append(f_motif)
        seq_motifs.append(r_motif)
        # score the motif with the score_kmer function from motif_ops.py

        f_score = score_kmer(f_motif, pwm)
        r_score = score_kmer(r_motif, pwm)


        # add this score to the list
        score_list.append(f_score)
        score_list.append(r_score)

    # return the score list
    return score_list, seq_motifs


def select_motif(score_list, seq_motifs):
    '''

    :param score_list: a list of pwm based scores for each of the motifs in the forward and reverse seq
    :return: a new motif to be put back in the larger pool
    '''

    # first we need to turn the list of scores into a list of probabilities, so let's initialize that list
    log_scores = []
    # for score in score_list:
    for score in score_list:
        log_score = (2**score)
        log_scores.append(log_score)
    # make a random choice from the seq list with the corresponding prob list as an argument
    return random.choices(seq_motifs, log_scores)
    # return the new motif to put back into the pool


promoter_pfm = GibbsMotifFinder(subseqs, 10, seed=2076)
print(promoter_pfm)