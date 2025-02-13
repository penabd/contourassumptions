#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.interpolate import splprep, splev
import sympy as sp
import math


def symbolic_gaussian(mu1=0, sigma1=1, mu2=0, sigma2=1):
    """This functino uses SymPy to symbolically represent a
    gaussian function of form f(x,y) = exp(-((x-mu)**2/2sigma) + (y-mu2)**2/2sigma2)
    
    """
    x, y = sp.symbols('x y')
    f = sp.exp(- ((x-mu1)**2 / (2*sigma1)) -  ((y-mu2)**2 / (2*sigma2)))

    return f, x, y

def calc_sym_gradient(f,x,y):
    """Note that x and y are the symbols created with sympy.
    """
    fx = sp.diff(f,x)
    fy = sp.diff(f,y)

    return [fx, fy]
 
def calc_sym_hessian(fx, fy, x, y):
    """Note that x and y are the symbols created with sympy.
    """

    fxx = sp.diff(fx, x)
    fxy = sp.diff(fx, y)
    fyx = sp.diff(fy, x)
    fyy = sp.diff(fy, y)

    return [fxx, fxy, fyx, fyy]

def calc_curvature(x_val, y_val, x, y, hessian, gradient):
    """Calculate curvature using cylindrical coordinates. 
    


    """
    fx = gradient[0].subs({x:x_val, y:y_val})
    fy = gradient[1].subs({x:x_val, y:y_val})

    fxx = hessian[0].subs({x:x_val, y:y_val})
    fxy = hessian[1].subs({x:x_val, y:y_val})
    fyx = hessian[2].subs({x:x_val, y:y_val})
    fyy = hessian[3].subs({x:x_val, y:y_val})

    knum = (fy**2)*fxx - 2*fx*fy*fxy + (fx**2)*fyy
    kden = (fx**2) + (fy**2)**1.5

    k = (knum / kden) * -1
    kevaled = k.evalf()

    return k, float(kevaled)

def is_convex(x_val, y_val, x, y, hessian):
    hessian_eval = [f.subs({x:x_val, y:y_val}).evalf() for f in hessian]
    hessian_eval = np.array(hessian_eval).astype('float64').reshape(2,2)
    print(hessian_eval)
    print(hessian_eval.shape)

    return np.all(np.linalg.eigvals(hessian_eval) > 0)

def filter_local_maxima(curvature_vals):
    filtered_curv_vals = []
    for i, (x_pt, y_pt, k) in enumerate(curvature_vals):
        if not math.isnan(k) and k >= 0: 
            if curvature_vals[i-1][2] < curvature_vals[i][2]:
                if i+1 >= len(curvature_vals):
                    j = 0
                else:
                    j = i+1
                if curvature_vals[j][2] < curvature_vals[j][2]:
                    filtered_curv_vals.append([x_pt, y_pt, k])
    
    return filtered_curv_vals


def case1():
    """Defines ellipse
    """
    # create gaussians
    f, x, y = symbolic_gaussian(0,0.1,0,0.25)
    sources = [(0,0)]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("gaussian_ellipse_data.txt")

    # Extract x, y, z values
    x_vals, y_vals, z_vals = data[:, 0], data[:, 1], data[:, 2]

    return f, x, y, x_vals, y_vals, z_vals, sources

def case2():
    """
    """
        # create gaussians
    f, x, y = symbolic_gaussian(0.5,0.05,0.5,0.05)
    f2 = sp.exp(- ((x)**2 / (2*0.1)) -  ((y)**2 / (2*0.1)))
    f3 = sp.exp(- ((x-(-0.35))**2 / (2*0.05)) -  ((y-(-0.35))**2 / (2*0.05)))

    f = f + f2 + f3

    sources = [(0.5,0.5), [-0.175, -0.175]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("gaussian_contour_data.txt")

    # Extract x, y, z values
    x_vals, y_vals, z_vals = data[:, 0], data[:, 1], data[:, 2]

    return f, x, y, x_vals, y_vals, z_vals, sources

def case3():
    """
    """
    # create gaussians
    f, x, y = symbolic_gaussian(0.45,0.1,0.45,0.1)
    f2, _, _ = symbolic_gaussian(-0.45,0.1,-0.45,0.1)
    # f3 = sp.exp(- ((x-(-0.35))**2 / (2*0.05)) -  ((y-(-0.35))**2 / (2*0.05)))

    f = f + f2

    sources = [(0.45,0.45), [-0.45, -0.45]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("gaussian_contours_data_3.txt")

    # Extract x, y, z values
    x_vals, y_vals, z_vals = data[:, 0], data[:, 1], data[:, 2]

    return f, x, y, x_vals, y_vals, z_vals, sources

def case4():
    """
    """
    # create gaussians
    f, x, y = symbolic_gaussian(0.75,0.1,0.75,0.1)
    f2, _, _ = symbolic_gaussian(0,0.1,0,0.1)
    f3, _, _ = symbolic_gaussian(-0.75,0.1,-0.75,0.1)

    f = f + f2 + f3

    sources = [(0.75,0.75), [-0.75, -0.75], [0,0]]
    
    # Load your scattered (x, y, z) data from the Gaussian contour file
    data = np.loadtxt("gaussian_contours_data_4.txt")

    # Extract x, y, z values
    x_vals, y_vals, z_vals = data[:, 0], data[:, 1], data[:, 2]

    return f, x, y, x_vals, y_vals, z_vals, sources

def find_sources(x,y, gradient):
    sources = sp.solve(gradient, (x,y))
    return sources

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
    plt.scatter(radius_of_curvs, dists_from_source, edgecolors='black',)
    min_val = min(min(radius_of_curvs), min(dists_from_source))  
    max_val = max(max(radius_of_curvs), max(dists_from_source))  
    plt.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--",)  
    plt.xlabel("Radius of Curvature")
    plt.ylabel("Distance from Nearest Source")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(contour_num, ratios)
    plt.xlabel("Contour Level f(x,y) = c")
    plt.ylabel("Radius of Curvature / Distance from Nearest Source")
    plt.show()



def main():
    f, x, y, x_vals, y_vals, z_vals, sources = case4() 
    
    # Plot contours using `tricontour`
    fig, ax = plt.subplots()
    contour = ax.tricontour(x_vals, y_vals, z_vals, levels=10)

    # Compute curvature for extracted contour points
    cps_and_curves = {}
    cps_grads_curves = {}
    for c, collection in zip(contour.levels, contour.collections):
        curvature_vals = []
        gradient_vals = []

        # define func for level set
        g = f - c
        gradient = [gx, gy] = calc_sym_gradient(g, x, y)
        hessian = calc_sym_hessian(gx, gy, x, y)

        for path in collection.get_paths():
            v = path.vertices
            x_pts, y_pts = v[:, 0], v[:, 1]

            # Compute curvature at each extracted contour point
            for (x_pt, y_pt) in zip(x_pts, y_pts):
                k, keval = calc_curvature(x_pt, y_pt, x, y, hessian, gradient)
                geval = [gx.subs({x:x_pt, y:y_pt}).evalf(), gy.subs({x:x_pt, y:y_pt}).evalf()]
                curvature_vals.append([x_pt, y_pt, keval])
                gradient_vals.append([x_pt, y_pt, geval[0] + geval[1]])


            # Find local curvature maxima
            # change to include segemnt length??
            # Ensure at least `num_maxcurv_points` points exist before slicing
            curvature_vals = np.array(curvature_vals)
            gradient_vals = np.array(gradient_vals).astype('float64')
            num_maxcurv_points = min(1, len(curvature_vals))

            # grab max curvurtature in x different regions
            num_regions = 8
            seg_length = len(x_pts) // num_regions


            if num_maxcurv_points > 0:
                # Get indices of the highest curvature values
                for j in range(num_regions):
                    start_idx = j * seg_length
                    end_idx = (j + 1) * seg_length if j < num_regions - 1 else len(curvature_vals)  
                    
                    region = curvature_vals[start_idx:end_idx] 
                    maxcurv_local_indices = np.argsort(region[:, 2])[-num_maxcurv_points:][::-1]  
                    maxcurv_global_indices = [start_idx + idx for idx in maxcurv_local_indices]

                    # convex_indices = np.where([
                    #                         is_convex(x_val, y_val, x, y, hessian) 
                    #                         for x_val, y_val, _ in curvature_vals[maxcurv_global_indices]
                    #                     ])[0] 
                    
                    # maxcurv_global_indices = maxcurv_global_indices[convex_indices]
                    filtered_curv_vals = np.array(filter_local_maxima(curvature_vals))

                    if len(filtered_curv_vals) > 0:
                        plt.plot(filtered_curv_vals[:,0], 
                                filtered_curv_vals[:, 1], 
                                'ro')  
                    else:
                        plt.plot(curvature_vals[maxcurv_global_indices, 0], 
                                curvature_vals[maxcurv_global_indices, 1], 
                                'go')  

                    cps_and_curves[c] = curvature_vals[maxcurv_global_indices,:]


                # plot where gradient is zero
                gzero_indices = np.where(np.isclose(gradient_vals[:, 2], 
                                                           0, 
                                                           atol=1e-4),)[0]
                # try:
                #     convex_indices = np.where(curvature_vals[gzero_indices,2] >= 0)[0]
                # except:
                #     print(curvature_vals[gzero_indices,2])
                # gzero_indices = gzero_indices[convex_indices]
                
                if len(gzero_indices) > 0:
                    plt.plot(gradient_vals[gzero_indices, 0],
                             gradient_vals[gzero_indices, 1],
                             'bo', markersize=8, 
                             label="Gradient Zero")



    plt.title("Contour Lines and Sampled Curvature Points (tricontour)")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")


    plt.show()

    # Before we show plot, test assumptions
    assumption2(cps_and_curves, sources)
    print(cps_grads_curves)

    # deal with gradients
    # assumption2(cps_grads_curves, sources)



if __name__ == '__main__':
    main()
