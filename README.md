This repository contains a simulation used to evaluate an assumption made about the relationship between radius of curvature and distance from a source for the work being done for contour estimation with drones.

Contour data is generated in test_assumptions.cpp, and analysis and plotting done in test_assumptions.py.

Test Cases:
- Case 1 elliptical contours
- Case 2 two gaussians with smaller distance between sources
- Case 3 two gaussians with larger distance (relative to case 2) between sources
- Case 4 three sources

Goal: Is the assumption valid?
Is radius of curvature $O(d)$?  If so, find asymptotic constant; if not, is there something we missed?


The following files have been taken from the sketchalgorithm repository:
- render.py
- util.hpp takes code from sketchalgorithautogaussian.cpp
