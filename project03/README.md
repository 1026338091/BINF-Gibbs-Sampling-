# Introduction
This project implements a Gibbs Sampling algorithm to identify transcription factor binding sites (motifs) from ChIP-seq data, utilizing a simplified Peak Calling strategy to enhance signal-to-noise ratio. By iteratively updating Position Weight Matrices (PWMs) and employing probabilistic sampling on both strands, the model converges to identifying statistically significant sequence patterns. The final results are visualized via Sequence Logos, demonstrating the efficacy of stochastic methods in uncovering hidden biological signals within complex genomic backgrounds.

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

    check the pfm from a previous iteration against current iteration, if they are similar enough then we've reached convergence
    
    return pfm
    

How we can implement it:


```
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

    initialization:
    first we'll want to choose our random motifs for each seq in the list (motifs = choose_motifs(seqs, k, rng))
    next we'll need to build our frequency matrix and weight matrix (pfm = build_pfm(), pwm = build_pwm()
    
    iteration:
    
    where this needs to be an iterative process of unknown duration we could pick a stopping point like 10000 (i.e. if we're still going after 10,000 iterations let's stop)

    could do a for i in range(x):
        pick a random motif from the motifs and keep track of it's index (motifs = random_seq(motifs, rng))
        
        create a new moving list to add all of our motifs to, except for the one we just selected

        
        every time we've iterated through say, 100 times:
            find the information content of the old pfm
            build a new pfm and pwm with the newest list of motifs
            calculate the ifc of the new pfm

            compare the new ifc to the old ifc, if they are deemed to be similar enough:
                return the pfm

        using the index of the randomly selected and left out motif, find the original sequence
        get the reverse_strand of that sequence with reverse_complement()
        get scores of all the motifs of the left out sequence and it's compliment with get_all_scores ()
        choose the new motif with select_motif()
        determine which strand the motif came from, and place the new motif and proper strand back into the pool of sequences and motifs

    if we reached the end of the loop, return the pfm


# The GibbsMotifFinder function relies on helper functions which were mentioned above, they are documented in the rest of this pseudocode #
# The goal of these functions was to break up as many steps of each iteration as possible #


# we first need to choose a random position for each as the initial motif in each sequence
# the position can be anywhere on the sequence that is less than len(sequence) - len(motif)
# we could make a function choose_motifs
def choose_motifs(seqs, k, rng):
    '''
    Args:
        seqs (str list): a list of sequences
        k (int): the length of the motif to find
        rng: the random number generator to use to get motif positions
    :return: a list of the selected motifs from each position
    '''
    initialize a motif list
    for seq in seqs
        pick a random integer between 0 and len(seq)
        append the motif list with seq[int: int + k]
    return the motif list



# now that we have all of our initial motifs, we need to create a frequency matrix of all nucleotides, and convert it into a probability weight matrix
# this frequency matrix accounts for the proportion that each nucleotide accounts for across each position in the motif
# we can use the pfm and pwm functions he provides us in the utility scripts (build_pfm() and build pwm()


# with our newly created pwm we now need to start the iteration process
# for each iteration we need to:
# choose a sequence at random
# recalculate the pwm, omitting the chosen sequence
# this might be beneficial to put into its own function, random_seq

def random_seq(motifs, rng):
    '''

    :param motifs: a list of motifs
    :param rng: random number generator to select the sequence to remove
    :return: motif list without the selected sequence and the integer that was used for selection
    '''

    choose a sequence at random (motif list and seq list should line up 1 to 1 so we can do this on the motif list)
    return motif to remove and the integer used to select it
    


# after identifying out left out seq, we need to score  possible motif in the sequence that was omitted
# that is, we need something like a list of scores, for each k (motif length) seq in the whole sequence
# we need to be careful here because we also need to score the reverse complement where the binding sites aren't strand specific
# we can use his seq_ops function reverse_complement() to get the reverse compliment of our left out sequence



# now we can score all motifs in the sequence and reverse compliment, we should probably call the get_reverse() function accordingly in the main function
# that way we're not computing it twice if we need to place it back into the full seq list
def get_all_scores(seq, rev_seq, k, pwm):
    '''

    :param seq: sequence to score
    :param rev_seq: reverse compliment of the seq to score
    :param k: motif length
    :param pwm: probability weight matrix to score the sequence motifs against
    :return: a list of scores for each motif in the forward and reverse sequence
    '''

    initialize a score list for forward and reverse scores
    initialize a motif list to store all scored motifs

    for each motif in the forward seq and reverse seq:
        score the motif with the score_kmer function from motif_ops.py
        add this score to the score list
        add the motif to the motif list
   

    return the score list and the motif list
    

# to select the new motif, we need to create a list of probabilities for choosing all the motifs we've scored
# so this should be something like the score of the motif, divided by the sum of all scores - do this for every motif, then use it to select from probabilistically
# the selected motif is what will be placed back into the pool of motifs for continued iteration
# let's find this motif with a function called select_motif()

def select_motif(score_list):
    '''

    :param score_list: a list of pwm based scores for each of the motifs in the forward and reverse seq
    :return: a new motif to be put back in the larger pool
    '''

    # first, since our pwm scores are in log2 space, we need raise 2 to the power of the score to translate it to a value that can be used appropriately for probabilities

    initialize a log transformed list
    for score in score_list:
        transformed value = 2^score
        append transformed list with value

    # now we can make a list of probailities all adding to 1, and select the motif using these weights\

    initialize a probability list
    for val in transformed list:
        probability = val / sum of the vals
        append probability list with the val

    make a random choice from the all the scored motifs with the corresponding prob list as an argument
    return the new motif to put back into the pool

```

# Successes
Description of the team's learning points

# Struggles
The three biggest struggles we faced in this project with respect to algorithmic implementation were identifying when and when not to recalculate the matrices, how to account for the reverse compliment sequence (we're still not sure if we did this one correctly), and how to appropriately assign probabilities to motifs from the left out sequence based on their scores from the pwm.

The probability frequency matrices and probability weight matrices are the objects used to score all the motifs in the left out sequence each iteration through the algorithm, and because the composition of motifs changes each iteration, they need to be frequently updated to reflect this change. However, it's not exactly clear how often they should be updated. Initially we were updating the matrices after every iteration, but we ultimately decided this might not be appropriate - if the pwm matrix that is used to score the new motif that will replace the old one is updated every iteration, then it never has time to influence the motif selection. To account for this we chose to recalculate the pfm and pwm only when we wanted to make a comparison, which we did after every 100 iterations. The selection of 100 iterations was relatively random, and while it appeared to work for Shine-Dalgarno sequence selection from a list of sequences guranteed to have the Shine-Dalgarno sequence, it may not be appropriate for more ambiguous and complex data sets.

The motifs we were attempting to identify were DNA binding site motifs that can be located on the forward and reverse strand of DNA. Because of this we needed to account for the forward and reverse strand of the sequence we were scoring (if we didn't, we could miss sequences that contain the motif), while still selecting the sequence probabilistically and not just taking the motif that had the highest score between the two sequences. What we decided to do was combine the transformed scores from all motifs in the forward and reverse strand into one list (alternating between forward and reverse so strand could be identified), normalize them so that their probabilities all added up to 1, then use that list to weight the decision. This strategy does gurantee probabilistic selection, however, we are not completely sure including all scores in one list is appropriate. Given that the number of motifs in a given sequence can be variable and this system can still be used for motif selection, we assume that it is acceptable, but if for some reason it's not or if there is a better method, we'd love to hear what other groups did!

Our last struggle was determining how to assign probabilities based on the pwm scores of each motif in the left out sequence. pwms represent log2 transformed values based on the frequency of observing As, Ts, Cs, and Gs at each position in the motif, where the larger the value the greater the degree of agreement between the motif and the pwm. That means that largest negatively scored motifs should have had the lowest probability of being selected, however, we initially had the probability as the absolute value of the motif's score divided by the sum of the absolute value of each motif in the sequence. This means we were weighting the motifs with the largest disagreement with the pwm very favorably, and thus selecting bad motifs to put back in the motif list. You can imagine that if both really good and really bad motifs are weighted equally you aren't going to have much luck with convergence, which is a problem. Luckily we were able to quickly identify this problem, and take our pwm scores out of log space before using them to create a list of probabilities, but it certainly hindered our initial efforts.

# Personal Reflections
## Group Leader - Connor
My biggest takeaway from this project is how crucial it is to break down problems, ideas, and algorithmic implementation into their smallest and most essential parts. There was a lot going on in this algorithm/method, and if I tried to digest it as whole, even now that I understand it better, I would get lost very quickly. However, taking the problem step by step - first going through each step in gibbs sampling with pen and paper, and then doing conversion into pseudocode, and only then into actual code made it doable and understandable. In regards to actual implementation the most helpful thing was indentiying what steps had to happen every iteration (or every x iterations for the comparison of matrices) and seperating those steps into their own functions if possible or necessary. Instead of exposing all operations of each sub-component of the iteration it allowed us to call the sub-component in the form of a function, which I found helpful for readability and debugging purposes.

## Other member
Hongyuan Deng：Through this project, I gained a deep understanding of stochastic algorithms, specifically fixing critical errors in log-likelihood probability conversion and iterative model convergence. I also learned to bridge the gap between raw data and algorithm efficiency by applying a simplified peak-calling heuristic to reduce noise in the input dataset.

# Generative AI Appendix
No AI was used for the creation of this project.
