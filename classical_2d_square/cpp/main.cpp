#include "head.h"
#include "class_ising.hpp"

int main() {

        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_mag("mag.dat", std::ios_base::app | std::ios_base::binary);
        // std::ofstream file_mag_error("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);

        time_t start, end;
        start = time(NULL);

        int size = 32;
        int numSample = 100;
        int numMeasure = 100;
        double b = 0.01;
        for (int l = 0; l < numSample; ++l) {
            Ising ising(0.0, size, size, "P", "P", b);
            int nx = ising.NumSiteX();
            int ny = ising.NumSiteY();
            int numSite = nx*ny;
            std::vector<int> config(numSite);
            ising.Initialization(config);
            ising.MkvEvolve(config, 10000); // Initial thermalization. 

            std::vector<double> ene;
            std::vector<double> mag;
            for (int i = 0; i < numMeasure; ++i) {
                // ising.PrintConfig(config);
                ene.push_back(ising.EnergyDensity(config));
                mag.push_back(abs(ising.Magnetization(config)));
                int cs = ising.Wolff(config);
                int numInterval = 5*int(numSite / cs);
                // std::cout << numInterval << std::endl;
                for (int l = 0; l < numInterval; ++l) { ising.Wolff(config); }
                }
            double res = Mean(mag);
            double resErr = StdErr(mag);
            file_mag.write((char*)(&res), sizeof(double));
            file_mag.write((char*)(&resErr), sizeof(double));
            std::cout << b << " energy: " << Mean(ene) <<" mag: " << Mean(mag) << " " << StdErr(mag) << std::endl;
            b += 0.01;
            }

        end = time(NULL);
        file_log << "Time: " << double(end-start) / 60.0 << " min" << std::endl;
        file_log.close();
        file_mag.close();
        return 1;
        }

double Mean(std::vector<double> vec) { return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size(); }

double StdErr(std::vector<double> vec) {
        int len = vec.size();
        double m = Mean(vec);
        double acc = 0.0;
        for (int i = 0; i < len; ++i) { acc += (vec[i]-m)*(vec[i]-m); }
        return std::sqrt(acc / (len*(len-1)));
        }

