import numpy as np
import math

# Q: What does this file do?
# A: It creates a binary file filled with random zeros and ones.

# Q: What is it good for?
# A: One intuition for the entropy of a source is that it is the amount of incompressible information produced by the source. The binary file created by this script contains data from a source with adjustable entropy. Compressing this binary file leads to a file size that is (usually) reduced w.r.t. the original, however, never below the entropy.

# Q: How to use this file?
# A: Run it (e.g., python mk_random_file.py) and it produces the file 'uncompressed.bin' of size 10 MB. Compress this file with a (lossless) compression tool of your choice (e.g., zip). Observe the resulting file size. Then, look for the line starting 'p = ...'. Change it, run again and compress again. Observe how the compressed file size changes.

def h(p):
    """ Binary entropy function """
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

def array_to_bytes(arr: np.ndarray):
    """ Convert an array of 0, 1 into bytes """
    num_bytes = math.ceil(len(arr) / 8)
    padded_arr = np.zeros(num_bytes * 8, dtype=int)
    padded_arr[:len(arr)] = arr[:]
    bits_arr = padded_arr.reshape(num_bytes, 8)
    values_arr = bits_arr * 2**np.arange(8)[::-1]
    bytes_list = [int(val) for val in np.sum(values_arr, axis=1)]
    return bytes(bytes_list)


rng = np.random.default_rng()

chunk_size = 8 * 2**20 # 1 MiB
num_chunks = 10

p = 0.2 # P(X = 1)

print(f'Entropy: {chunk_size*num_chunks} * {h(p)} bit == {chunk_size*num_chunks/8*h(p)/2**20:.2f} MB')

with open('uncompressed.bin', 'bw') as binary_file:
    for _ in range(num_chunks):
        data = rng.binomial(1, p, chunk_size)
        binary_file.write(array_to_bytes(data))
