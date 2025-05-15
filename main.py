import numpy as np
import math

def split_2real_fft(x):
    z = x[::2] + x[1::2] * 1j  # Splitting odd and even
    Z = np.fft.fft(z)
    Zminconj = np.roll(np.flip(Z), 1).conj()
    Zx =  0.5  * (Z + Zminconj)
    Zy = -0.5j * (Z - Zminconj)

    N = len(x)
    W = np.exp(-1j * 2 * math.pi * np.arange(N//2) / N)
    Zall = np.concatenate([Zx + W*Zy, Zx - W*Zy])
    print(W*Zy)
    return Zall

# Testing
x = np.array([1, 2, 3, 4, 5, 6, 7, 8])

X = split_2real_fft(x)
