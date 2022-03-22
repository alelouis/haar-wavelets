__author__ = "Alexis LOUIS"
__version__ = "1.0"
__description__ = "Haar Discrete Wavelet Transform and its Inverse pure Python implementation"

import random, math

def forward(data):
    # Forward Haar transform, 1-pass
    approx = list(map(lambda x: sum(x)/2, zip(data[::2], data[1::2])))
    detail = [first - mean for first, mean in zip(data[::2], approx)]
    return approx, detail

def dwt(data, depth):
    # Haar Discrete Wavelet Transform, depth passes
    depth = min(int(math.log2(len(data))), depth)
    details = []
    for d in range(depth):
        approx, detail = forward(data)
        details = detail + details
        data = approx
    return approx, details

def idwt(approx, details):
    # Haar Inverse Discrete Wavelet Transform, depth passes
    N = len(approx + details)
    sums = [sum([[k] * (N // len(approx)) for k in approx], [])]
    bs, i = N//2, 0
    while bs > 0 and i<len(details):
        values = sum([sum([[-a]*(N//(2*bs)), [a]*(N//(2*bs))], []) for a in details[::-1][i:i+bs]], [])
        sums.append(values[::-1])
        i += bs
        bs //= 2
    rec_data = [sum([s[j] for s in sums]) for j in range(N)]
    return rec_data

if __name__ == "__main__":
    # Generate data (size 2**n)
    data = [float(random.randint(-100, 100)) for _ in range(2**10)]
    # DWT & IDWT
    data_rec = idwt(*dwt(data, 3))
    # Check reconstruction
    assert data == data_rec