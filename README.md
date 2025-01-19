# 🧬 Viterbi Algorithm for gene prediction

In this bioinformatics assignment, I was tasked to use the Viterbi algorithm, which is a dynamic programming algorithm, to find the sequence of states (in this case: introns, exons, and splice sites) with the highest likelihood of occurring in a given DNA sequence.  The algorithm uses the initial probabilities of the states as well as the transition and emission probabilities in order to do this. 

## How does it work?

### Initialization:

For each state, it computes the probability of starting in that state and emitting the first observation. Store this in the P matrix and note the previous state as "Start" in the backpointer.

### Recursion:

For each time step and state, compute the maximum probability of transitioning from a previous state to the current state while emitting the current observation. Store the probabilities in P and the state leading to this maximum probability in backpointer.

### Finalization:

Compute the best final state by considering the transition to an end state. Use the backpointer to trace back the sequence of states leading to the highest probability path.

### Output:

Print the probability matrix P.
Retrieve the best path using the backpointer.
