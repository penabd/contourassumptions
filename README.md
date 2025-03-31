This repository contains a simulations used to evaluate an assumption made about the relationship between radius of curvature and distance from a source for the work being done for contour estimation with drones, as well as to compute curvature using LSQ. 

Files:
- gen_data.cpp : used to generate files in data folder
- util.hpp : contains functions used to create cpp gaussians and ellipses
- test_assumptions.py : used to test assumption about curvature and source
- test_numeric_curvature.py : computes curvature using LSQ and plots info about curvature and gradient of curvature


Test Cases:
- Case 1 elliptical contours
- Case 2 two gaussians with smaller distance between sources
- Case 3 two gaussians with larger distance (relative to case 2) between sources
- Case 4 three sources


The following files have been taken from the sketchalgorithm repository:

    util.hpp takes code from sketchalgorithautogaussian.cpp
