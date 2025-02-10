#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def main():
    lines = []
    with open("gaussian_contour_data.txt") as f:
        for line in f:
            if line.strip():
                lines.append(line.strip())
    
    data = np.array([list(map(float, ln.split())) for ln in lines])
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    plt.figure(figsize=(7,6))
    contour = plt.tricontour(x, y, z, levels=20, cmap='viridis')
    plt.clabel(contour, inline=True, fontsize=8)
    
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Contour Plot")
    plt.show()


if __name__ == '__main__':
    main()
