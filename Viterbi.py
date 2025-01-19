# -*- coding: utf-8 -*-
"""

# Viterbi Algorithm to find the path with highest likelihood

Below you can find short explanations of the variables used :
- delta_i is the highest probability along a single path
- obs is observations
- T is length of observations
- states are all your states
- q is state q in all your states
- phi is start/initial probabilies
- trans_p is transition probabilities (A)
- emit_p is emission probabilities (B)
"""

# Helps visualize the steps of Viterbi.
def print_path_probability_matrix(P,states):
    for q in states: print(q+ ":\t" + "\t".join(("%.2e" % (p)) for p in P[q]))

"""Your assignment starts from here. You will fill out the missing lines in the implementation of viterbi algorithm below."""

def viterbi(obs, states, phi, trans_p, emit_p):
    # We will create two variables; P and backpointer.
    # Both variables (P and backpointer) are dictionaries.

    P = {}
    backpointer = {}
    T=len(obs)

    # Initialization step

    for q in states:
        P[q] = [phi[q] * emit_p[q][obs[0]]]
        backpointer[q] = ["Start"]


    # Recursion step
    for t in range(1,T):
      for s in states:
        probs = []
        for s0 in states:
          probs.append((P[s0][t-1]) * (trans_p[s0][s]) * (emit_p[s][obs[t]]))



        P[s].append(max(probs))
        max_idx = probs.index(max(probs))
        backpointer[s].append(states[max_idx])




    # Finalization step
    (delta_i, state)=max((P[q][-1] * trans_p[q]["end"], q) for q in states)
    P["End"]=delta_i
    backpointer["End"]=state

    print_path_probability_matrix(P,states)

    # To retrieve the best path
    path = []
    curr_state = state
    for idx in range(1,T+1):
      path.append(curr_state)
      curr_state = backpointer[curr_state][-(idx)]



    return ("probability of the best path to observe the sequence: %.2e \nthe best path: %s" % (delta_i, "\t".join(path[::-1])))

states = ('Exon', 'Splice', 'Intron' )
seq="CTTCATGTGAAAGCAGACGTAAGTCA"
observations=tuple(seq)
start_probability = {'Exon': 1.0, 'Splice': 0.0, 'Intron': 0.0}
transition_probability = {
   'Exon' : {'Exon': 0.9, 'Splice': 0.1, 'Intron': 0.0, 'end': 0.0},
   'Splice' : {'Exon': 0.0, 'Splice': 0.0, 'Intron': 1.0, 'end': 0.0},
   'Intron' : {'Exon': 0.0, 'Splice': 0.0, 'Intron': 0.9, 'end': 0.1}
   }
emission_probability = {
   'Exon' : {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
   'Splice' : {'A': 0.05, 'C': 0.0, 'G': 0.95, 'T': 0.0},
   'Intron' : {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T': 0.4},
   }


final=viterbi(observations, states, start_probability, transition_probability, emission_probability)
print(final)
