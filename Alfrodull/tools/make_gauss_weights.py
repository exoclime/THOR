#! /usr/bin/env python3
"""Print gauss-legendre polynomials weights """

from numpy.polynomial.legendre import leggauss as G

length = 100


outfile = open("gauss_legendre_weights.cpp", "w")

outfile.write("\n")
outfile.write(f"double gauss_legendre_weights[{length}][{length}] = {{\n")

for i in range(length):
    deg = i+1
    wg = G(deg)[1]
    outfile.write("{ ")

    l = len(wg)
    for w in range(0, l):
        outfile.write(f"{wg[w]:.16e}")
        if w != length - 1:
            outfile.write(", ")
    for w in range(l, length):
        outfile.write(f"0.0")
        if w != length - 1:
            outfile.write(", ")

    outfile.write("}")
    if i != length - 1:
        outfile.write(",")
    outfile.write("\n")


outfile.write("};\n")
