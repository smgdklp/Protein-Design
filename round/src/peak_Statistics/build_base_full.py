"""Base = smoothed full empirical histogram (includes the peak regions).

Every context keeps its own base: circular Gaussian smoothing sigma=2.5 of the
empirical 72x72 histogram, normalised. Peaks are only markers in the figure.
"""
import json
import os

import numpy as np

NB = 72
BW = 5.0
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
ROUND = r"D:\CodingPrograms\Protein Design\round"
PROBS = os.path.join(ROUND, "result", "static", "triplet_probs.npy")
OUT = os.path.join(ROUND, "result", "ground", "base_full")
SIGMA = 2.5


def ksig(n, s):
    i = np.arange(n, dtype=float)
    d = np.minimum(i, n - i)
    k = np.exp(-0.5 * (d / s) ** 2)
    return k / k.sum()


def main():
    os.makedirs(OUT, exist_ok=True)
    probs = np.load(PROBS, mmap_mode="r")
    k = ksig(NB, SIGMA)
    K = np.fft.fft(k)[:, None] * np.fft.fft(k)[None, :]
    for b in range(20):
        ch = ALPHABET[b]
        arr = np.zeros((20, 20, NB, NB), dtype=np.float32)
        for a in range(20):
            for c in range(20):
                idx = a * 400 + b * 20 + c
                m = np.asarray(probs[idx], dtype=np.float64)
                m /= m.sum()
                m = np.fft.ifft2(np.fft.fft2(m) * K).real
                m = np.clip(m, 0.0, None)
                m /= m.sum()
                arr[a, c] = m.astype(np.float32)
        np.save(os.path.join(OUT, f"center_{ch}.npy"), arr)
        print("base_full", ch)


if __name__ == "__main__":
    main()
