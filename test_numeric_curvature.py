#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import pandas as pd
from matplotlib.colors import TwoSlopeNorm, LogNorm

from curvature_solve import calc_curvature_LSQ, curvature_derivative

plt.rcParams['text.usetex'] = True

def case1():
    """Defines ellipse
    """
    sources = [(0,0)]
    
    # Load scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_ellipse_data.txt")

    return data, sources

def case2():
    """
    """
    sources = [(0.25,0.25), [-0.25, -0.25]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_2.txt")

    return data, sources


def case3():
    """
    """
    sources = [(1,0.5), [1, -0.5]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_3.txt")

    return data, sources


def case4():
    """
    """
    sources = [(1,1), [0,-1], [-1,1]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("data/gaussian_contours_data_4.txt")

    return data, sources


def euc_dist(x1,y1,x2,y2):
    return ((x2-x1)**2 + (y2-y1)**2)**0.5

def closest_source(x, y, sources):
    closest = sources[0]
    prev_dist = euc_dist(x,y,sources[0][0], sources[0][1])

    for source in sources[1:]:
        dist = euc_dist(x,y,source[0], source[1])
        if dist < prev_dist:
            closest = source

    return closest

def assumption2(curvature_info, sources):
    radius_of_curvs = []
    dists_from_source = []
    contour_num = []
    ratios = []
    for ctr in curvature_info.keys():
        curve_info = curvature_info[ctr]
        for x, y, k in curve_info:
            if not math.isnan(k):
                xs, ys = closest_source(x,y,sources)
                # print(x,y,k)
                # print(type(x), type(y), type(k))
                dist = euc_dist(x, y, xs, ys)
                print(dist)
                radius_of_curv = 1 / k
                radius_of_curvs.append(radius_of_curv)
                dists_from_source.append(dist)
                ratios.append(radius_of_curv/dist)
                contour_num.append(ctr)

    plt.figure()
    plt.scatter(contour_num, dists_from_source, edgecolors='black', label="DFNS")
    plt.scatter(contour_num, radius_of_curvs, edgecolors='black', label="Radius of Curvature")
    plt.xlabel("Contour Level f(x,y) = c")
    plt.ylabel("Distance from Nearest Source and Radius of Curv.")
    plt.legend()
    plt.show()

    plt.figure()
    plt.scatter(radius_of_curvs, dists_from_source, edgecolors='black',)
    min_val = min(min(radius_of_curvs), min(dists_from_source))  
    max_val = max(max(radius_of_curvs), max(dists_from_source))  
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--",)  
    plt.xlabel("Radius of Curvature")
    plt.ylabel("Distance from Nearest Source")
    plt.show()

    plt.figure()
    plt.plot(contour_num, ratios)
    plt.xlabel("Contour Level f(x,y) = c")
    plt.ylabel("Radius of Curvature / Distance from Nearest Source")
    plt.axhline(y = max(ratios), linestyle="--", color="black", 
                label = f"possible c = {np.round(max(ratios),3)}")
    plt.legend(loc = 'lower right')
    plt.show()

def plot_curvature_on_contours(data, curvature_vals):
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:,0], data[:,1], data[:,2], levels=10)
    
    curvature_vals = np.array(curvature_vals)
    positive_curvature = curvature_vals[curvature_vals[:, 2] > 0]
    negative_curvature = curvature_vals[curvature_vals[:, 2] < 0]

    ax.plot(positive_curvature[:, 0], positive_curvature[:, 1], 'go', label='Positive Curvature')
    ax.plot(negative_curvature[:, 0], negative_curvature[:, 1], 'ro', label='Negative Curvature')

    plt.title("Curvature Values on Contours")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend()
    plt.show()

def plot_curvature_heatmap_on_contours(data, curvature_vals):
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:,0], data[:,1], data[:,2], levels=10)
    
    curvature_vals = np.array(curvature_vals)
    x = curvature_vals[:, 0]
    y = curvature_vals[:, 1]
    z = curvature_vals[:, 2]

    # Create a scatter plot with a heatmap based on curvature values
    sc = ax.scatter(x, y, c=z, cmap='cividis_r', edgecolor='k', size=0.1)
    plt.colorbar(sc, label='Curvature')

    plt.title("Curvature Heatmap on Contours")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()

def plot_curvature_gradient_heatmap(data, curvature_vals, dk):
    """
    Plots a heatmap of the gradient of the curvature (dk) for each contour.
    """
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:, 0], data[:, 1], data[:, 2], levels=10,
                            colors='grey',
                            size = 5)

    curvature_vals = np.array(curvature_vals)
    x = curvature_vals[:, 0]
    y = curvature_vals[:, 1]
    dk = np.array(dk)
    
    norm = TwoSlopeNorm(vmin=np.min(dk), vcenter=0, vmax=np.max(dk)) 
    # norm = TwoSlopeNorm(vmin=np.min(dk), vcenter=0, vmax=np.max(dk)) 

    # Create a scatter plot with a heatmap based on the gradient of curvature
    sc = ax.scatter(x, y, c=dk, cmap='seismic', norm=norm, alpha=1, s=10)
    plt.colorbar(sc, label='Gradient of Curvature (dk)')

    plt.title("Gradient of Curvature Function")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()

def plot_abs_curvature_gradient_heatmap(data, curvature_vals, dk):
    """
    Plots a heatmap of the gradient of the curvature (dk) for each contour.
    """
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:, 0], data[:, 1], data[:, 2], levels=10,
                            colors='grey',
                            size = 5)

    curvature_vals = np.array(curvature_vals)
    x = curvature_vals[:, 0]
    y = curvature_vals[:, 1]
    dk = np.array(dk)
    dk_norm = dk / np.nanmean(dk)
    
    norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1) 

    # Create a scatter plot with a heatmap based on the gradient of curvature
    sc = ax.scatter(x, y, c=dk_norm, cmap='viridis_r', norm=norm, alpha=1, s=10)
    plt.colorbar(sc, label=r'$(\nabla \kappa)^2$')

    plt.title(r"Normalized $(\nabla \kappa)^2$")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()

def main():
    testCase = "case4"
    if testCase == "case1":
        data, sources = case1()
    elif testCase == "case2":
        data, sources = case2() 
    elif testCase == "case3":
        data, sources = case3()
    elif testCase == "case4":
        data, sources = case4()
    else:
        raise ValueError("Invalid test case")
    step = 0.01
    
    # Plot contours using `tricontour`
    fig, ax = plt.subplots()
    contour = ax.tricontour(data[:,0], 
                            data[:,1],
                            data[:,2], levels=10,
                            colors="grey",)

    # Compute curvature for extracted contour points
    cps_and_curves = {}
    cps_grads_curves = {}
    curvature_all = []
    max_curvs = []
    grad_zero = []
    dk_all = []
    dk2_all = []
    firstContour = True
    for c, collection in zip(contour.levels, contour.collections):
        curvature_vals = []
        gradient_vals = []
        
        already_visited = set()
        for path in collection.get_paths():
            v = path.vertices
            x_pts, y_pts = v[:, 0], v[:, 1]
            if len(x_pts) < 2:
                continue

            # Compute curvature at each extracted contour point
            for (x_pt, y_pt) in zip(x_pts, y_pts):
                x_pt = np.round(x_pt, 2)
                y_pt = np.round(y_pt, 2)

                if (x_pt, y_pt) in already_visited:
                    continue
                else:
                    already_visited.add((x_pt, y_pt))

                a = (x_pt - step, y_pt)
                b = (x_pt, y_pt + step)
                c = (x_pt + step, y_pt)
                d = (x_pt, y_pt - step)
                e = (x_pt, y_pt)
                f = (x_pt + step, y_pt + step)
                g = (x_pt - step, y_pt - step)
                h = (x_pt + step, y_pt - step)
                i = (x_pt - step, y_pt + step)

                gradB, curvature = calc_curvature_LSQ(a,
                                                      b,
                                                      c,
                                                      d,
                                                      e,
                                                      f,
                                                      g,
                                                      h,
                                                      i,
                                                      data)
                # curvature = np.round(curvature, 2)
                curvature_vals.append([x_pt, y_pt, curvature])
                # gradient_vals.append([x_pt, y_pt, gradB[0], gradB[1]])
                curvature_all.append([x_pt, y_pt, curvature])

            curvature_vals = np.array(curvature_vals)
            gradient_vals = np.array(gradient_vals).astype('float64')
            num_maxcurv_points = min(1, len(curvature_vals))

            dk, dk2 = curvature_derivative(curvature_vals)
            for d in dk:
                dk_all.append(d)
            for d2 in dk2:
                dk2_all.append(d2)
   
            if testCase != "case4":
                num_regions = 2
            else:
                num_regions = 3

            seg_length = len(curvature_vals) // num_regions

            if num_maxcurv_points > 0:
                for j in range(num_regions):
                    start_idx = j * seg_length
                    end_idx = (j + 1) * seg_length if j < num_regions - 1 else len(curvature_vals)  
                    
                    # get max curvature points and where kprime is zero
                    region = curvature_vals[start_idx:end_idx] 
                    maxcurv_local_indices = np.argsort(region[:, 2])[-num_maxcurv_points:][::-1]
                    maxcurv_global_indices = [start_idx + idx for idx in maxcurv_local_indices]

                    dk_region = dk2[start_idx:end_idx]
    
                    #eps = 1e-9
                    # zero_gradient_indices = np.where((np.abs(dkdx_region) < eps) 
                    #                                  & (np.abs(dkdy_region) < eps))

                    # zero_gradient_/indices = np.where(np.abs(dk_region) < eps)
                    min_gradient_indices = np.where(np.abs(dk_region) == np.min(np.abs(dk_region)))

                    
                    # kpzero_global_indices = [start_idx + idx for idx in zero_gradient_indices]
                    min_grad_global_indices = [start_idx + idx for idx in min_gradient_indices]


                    plt.plot(curvature_vals[maxcurv_global_indices, 0], 
                            curvature_vals[maxcurv_global_indices, 1], 
                            'ro', alpha=0.5,)  
                        

                    # plt.plot(curvature_vals[kpzero_global_indices, 0], 
                    #         curvature_vals[kpzero_global_indices, 1],
                    #         'go', alpha=0.5,) 
                    
                    plt.plot(curvature_vals[min_grad_global_indices, 0], 
                            curvature_vals[min_grad_global_indices, 1],
                            'mo', alpha=0.5,) 
                    
                    # print(f"Max Curvature Points: {curvature_vals[maxcurv_global_indices]}")

                    cps_and_curves[c] = curvature_vals[maxcurv_global_indices,:]

                    max_curvs.append(curvature_vals[maxcurv_global_indices, 2])
                    # grad_zero.append(curvature_vals[kpzero_global_indices, 2][0])


    plt.title("Contour Lines and Maximum Curvature Points")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend()
    plt.show()


    # Plot curvature gradient heatmap on contours
    plot_curvature_gradient_heatmap(data, curvature_all, dk_all)
    plot_abs_curvature_gradient_heatmap(data, curvature_all, dk2_all)


if __name__ == '__main__':
    main()
