// Random DNA sequence in Go

package main 

import "fmt"
import gonum as np
func random_dna_sequence(length, gc_share=None, probas=None, seed=None) {
    if seed is not None {
    np.random.seed(seed))
    }
    if gc_share is not None {
    g_or_c = gc_share / 2.0
    not_g_or_c = (1 - gc_share) / 2.0
    type probas struct {
        "G" var
        "C" var
        "A" var
        "T" var
    }
    probas := probas{
        "G": g_or_c,
        "C": g_or_c,
        "A": not_g_or_c,
        "T": not_g_or_c,
    }
    if probas is None {
    sequence = np.random.choice(list("ATCG"), length)
    } else {
    bases, probas = zip(*probas.items())
    sequence = np.random.choice(bases, length, p=probas)
    }
    return "".join(sequence)
    }