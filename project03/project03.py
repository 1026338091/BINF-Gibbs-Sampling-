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

    '''

    # 2000 for test
    for i in range(2000):
        
        # random choose
        idx = rng.integers(0, len(seqs))
        
        # Motifs 
        # remove the motif
        # use pop() to grasp
        current_motifs = motifs.copy() 
        current_motifs.pop(idx)
        
        # caculate the new pfm pwm
        temp_pfm = build_pfm(current_motifs, k)
        temp_pwm = build_pwm(temp_pfm)
        
        # scan the choosen seq by the temp_pwm (both posti and rev)
        target_seq = seqs[idx]
        target_seq_rev = reverse_complement(target_seq)
        
        # get_all_scores
        scores, candidates = get_all_scores(target_seq, target_seq_rev, k, temp_pwm)
        
        # new Motif
        
        new_motif_list = select_motif(scores, candidates)
        new_motif = new_motif_list[0]
        
        # put back new motif
        motifs[idx] = new_motif
        
        # check each 100 times to find out if IC was developing?
        if i % 100 == 0:
            current_pfm = build_pfm(motifs, k)
            current_ic = pfm_ic(current_pfm)
            print(f"Iteration {i}, IC Score: {current_ic:.4f}")


    # use converge motifs to build PFM 
    final_pfm = build_pfm(motifs, k)
    return final_pfm


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
