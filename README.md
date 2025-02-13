This repository contains a simulation used to evaluate an assumption made about the relationship between radius of curvature and distance from a source for the work being done for contour estimation with drones.

Test Cases:
[] circle (trivial)
[] ellipse
[] 2 Gaussians - vary variance but not mean, and vary distance between sources
[] 3 Gaussians - go for really large distance from source

Goal: Is the assumption valid?
Is radius of curvature $O(d)$?  If so, find asymptotic constant; if not, is there something we missed?


The following files have been taken from the sketchalgorithm repository:
- render.py
- util.hpp takes code from sketchalgorithautogaussian.cpp