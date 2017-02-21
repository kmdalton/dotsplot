from matplotlib import pyplot as plt
import numpy as np


def encode_sequence(sequence, k, **kwargs):
    """
    encode_sequence(sequence, k, **kwargs)

    Parameters
    ----------
    sequence : str
        Sequences you wish to encode
    k : int
        Word length
    Returns
    -------
    encoded_seq : numpy.ndarray
        array of integers describing the encoded sequence
    kmer_encoding : dict
        the kmers corresponding to each int in encoded_seq
    """
    kmer_encoding = kwargs.get('kmer_encoding', None)
    if kmer_encoding is not None:
        n = max(kmer_encoding.keys()) + 1
        kmer_encoding = {v:k for k,v in kmer_encoding.iteritems()}
    else:
        n = 0
        kmer_encoding={}

    encoded = []
    for i in range(len(sequence) - k):
        subsequence = sequence[i:i+k]
        if subsequence not in kmer_encoding:
            kmer_encoding[subsequence] = n
            n += 1
        encoded.append(kmer_encoding[subsequence])
    
    return np.array(encoded), {v:k for k,v in kmer_encoding.iteritems()}

def get_matches(encoded_seq):
    """
    get_matches(encoded_seq)

    Parameters
    ----------
    encoded_sequence : numpy.ndarray
        array of integers describing the sequence. generate with dotsplot.encode_sequence

    Returns
    -------
    X : numpy.ndarray
        array of matches, first coordinate
    Y : numpy.ndarray
        array of matches, second coordinate
    """
    X,Y = [],[]
    for x,kmer in enumerate(encoded_seq):
        #print np.where(encoded_seq == kmer)
        y = np.where(encoded_seq == kmer)[0]
        Y = np.concatenate((Y, y, [x]*len(y)))
        X = np.concatenate((X, [x]*len(y), y))
    return X,Y

def get_pairwise_matches(encoded_seqA, encoded_seqB):
    """
    get_matches(encoded_seqA, encoded_seqB)

    Parameters
    ----------
    encoded_seqA : numpy.ndarray
        array of integers describing the sequence. generate with dotsplot.encode_sequence
    encoded_seqB : numpy.ndarray
        array of integers describing the sequence. generate with dotsplot.encode_sequence

    Returns
    -------
    X : numpy.ndarray
        array of matches, first coordinate
    Y : numpy.ndarray
        array of matches, second coordinate
    """
    X,Y = [],[]
    for x,kmer in enumerate(encoded_seqA):
        #print np.where(encoded_seq == kmer)
        y = np.where(encoded_seqB == kmer)[0]
        Y = np.concatenate((Y, y))
        X = np.concatenate((X, [x]*len(y)))
    return X,Y



if __name__=="__main__":
    seqlenA = 1000 #shorter seq
    seqlenB = 1500 #longer seq
    k = 4
    seqA = ''.join(np.random.choice(['A', 'T', 'C', 'G'], seqlenA))

    encoded_seqA, kmer_encoding = encode_sequence(seqA, k)
    X,Y = get_matches(encoded_seqA)
    plt.scatter(X, Y)
    plt.title("{}-mer Self Similarity".format(k))
    plt.show()
    plt.title("{}-mer Self Similarity".format(k))

    #seqB will be seqB + a random string
    seqB = seqA + ''.join(np.random.choice(['A', 'T', 'C', 'G'], seqlenB - seqlenA))
    encoded_seqB, kmer_encoding = encode_sequence(seqB, k, kmer_encoding=kmer_encoding)

    X,Y = get_pairwise_matches(encoded_seqB, encoded_seqA)
    plt.scatter(X,Y)
    plt.title("{}-mer Pairwise Similarity".format(k))
    plt.show()
