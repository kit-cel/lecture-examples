import numpy as np
import galois
GF_2 = galois.GF(2)


def pnsequence(prim_poly, seed, sequence_length):
    '''
    binary representation of a primitive polynomial over F2 in galois is, eg
    x³+x+1 --> GF_2([1,0,1,1])
    '''
    # Thus:
    buffer = seed#GF_2.Zeros(order)
    output = np.zeros(sequence_length,dtype=int)
    coefficients = np.flip(prim_poly.coeffs)[1:]

    for i in range(sequence_length):
        output[i]=(buffer*coefficients).sum() # performs elementwise multiplication and then sum over the coefficient
        buffer=np.roll(buffer,1) #shift all elements in the buffer one to the right
        buffer[0]=output[i] # leftmost bufferstate corresponds to the current output
    return output


# def zdf_element(index,offset,M):

# vectorize_zdf_element=
def zdcsequence(sequence_length,M):
    offset=sequence_length%2
    increasing_index=np.arange(0,sequence_length,1)

    exponent=-1j*np.pi*M/sequence_length*increasing_index*(increasing_index+offset+2)

    return np.exp(exponent)




