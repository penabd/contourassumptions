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

// std::vector<utils::Point> getCriticalPoints(){
//     return
// }


int main(){
    int n = 10;

    double xMin = -3.0, xMax = 3.0;
    double yMin = -3.0, yMax = 3.0;
    double step = 0.01;

    // Define gaussian vars and centers
    std::vector<utils::Point> gaussianCenters;
    gaussianCenters.push_back(utils::Point(1,1));
    gaussianCenters.push_back(utils::Point(0,-1));
    gaussianCenters.push_back(utils::Point(-1,1));

    std::vector<utils::Point> gaussianVars;
    gaussianVars.push_back(utils::Point(0.1,0.5));
    gaussianVars.push_back(utils::Point(0.5,0.1));
    gaussianVars.push_back(utils::Point(0.4,0.4));

    //define ellipse
    // std::vector<utils::Point> gaussianCenters;
    // gaussianCenters.push_back(utils::Point(0, 0));

    // std::vector<utils::Point> gaussianVars;
    // gaussianVars.push_back(utils::Point(0.05, 0.125));



    // create file to save data and gen contours
    std::ofstream contourFile("gaussian_contours_data_5.txt");
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