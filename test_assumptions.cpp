#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "util.hpp"

using namespace std;

int main(){
    int n = 10;

    double xMin = -1.0, xMax = 1.0;
    double yMin = -1.0, yMax = 1.0;
    double step = 0.01;

    // Define gaussian vars and centers
    std::vector<utils::Point> gaussianCenters;
    gaussianCenters.push_back(utils::Point(0.5,0.5));
    gaussianCenters.push_back(utils::Point(-0.5,-0.5));

    std::vector<utils::Point> gaussianVars;
    gaussianVars.push_back(utils::Point(0.1,0.1));
    gaussianVars.push_back(utils::Point(0.5,0.5));

    // create file to save data and gen contours
    std::ofstream contourFile("gaussian_contour_data.txt");
    contourFile << std::fixed << std::setprecision(5);

    for (double x = xMin; x <= xMax; x += step) {
        for (double y = yMin; y <= yMax; y += step) {
            utils::Point p(x, y);
            double total = 0.0;
            for (int j = 0; j < gaussianCenters.size(); j++) {
                total += utils::gaussian(p, gaussianCenters[j], gaussianVars[j]); 
            }
            contourFile << x << " " << y << " " << total << "\n";
        }
        contourFile << "\n"; 
    }

    contourFile.close();
    
    return 0;
}