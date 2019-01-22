// The class Ising defines an object, namely a physical system, which consists of a Hamiltonian, a kind of lattice.
class Ising {
        private:
        // Parameters in the Hamiltonian.
        double mag;

        // Define the square lattice.
        int numSiteX;
        int numSiteY;
        std::string BCX; // Boundary conditions. Only valid for "P" or "O";
        std::string BCY;

        // Inversed temperature.
        double beta;

        public:
        int NumSiteX();
        int NumSiteY();

        int Coordinate(int x, int y);
        void Initialization(std::vector<int>& config);
        void SimpleUpdate(std::vector<int>& config);
        int Wolff(std::vector<int>& config);
        void MkvEvolve(std::vector<int>& config, int num);
        double EnergyDensity(std::vector<int> config);
        double Magnetization(std::vector<int> config);
        double MagnetizationQuad(std::vector<int> config);
        double MagnetizationBiquad(std::vector<int> config);
        void PrintConfig(std::vector<int> config);
        
        Ising(double h, int nx, int ny, std::string bcx, std::string bcy, double b);
        ~Ising(); 
        };

// Constructor.
Ising::Ising(double h, int nx, int ny, std::string bcx, std::string bcy, double b) {
        mag = h;
        numSiteX = nx;
        numSiteY = ny;
        BCX = bcx;
        BCY = bcy;
        beta = b;
        }

// Destructor.
Ising::~Ising() {}

int Ising::NumSiteX() { return numSiteX; }

int Ising::NumSiteY() { return numSiteY; }

// Define the mapping rule of coordinates on square lattice to one dimensional vector.
int Ising::Coordinate(int x, int y) { return y*numSiteX+x; }

void Ising::Initialization(std::vector<int>& config) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> Dis(0, 1);

        for (std::vector<int>::iterator it = config.begin(); it != config.end(); ++it) {
            config[std::distance(config.begin(), it)] = 2*Dis(gen)-1;
            }
        }

void Ising::SimpleUpdate(std::vector<int>& config) {
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_int_distribution<int> DisX(0, numSiteX-1);
        std::uniform_int_distribution<int> DisY(0, numSiteY-1);
        int x = DisX(gen);
        int y = DisY(gen);
        int s = config[Coordinate(x, y)];
        int nsx = config[Coordinate((x+1) % numSiteX, y)]+config[Coordinate((x-1+numSiteX) % numSiteX, y)];
        int nsy = config[Coordinate(x, (y+1) % numSiteY)]+config[Coordinate(x, (y-1+numSiteY) % numSiteY)];
        double cost = 2.0*s*(nsx+nsy);
        std::uniform_real_distribution<long double> Prob(0.0, 1.0);
        double r = Prob(gen);
        if (r < std::exp(-1.0*beta*cost)) { config[Coordinate(x, y)] = -1*s; }
        }

// Wolff cluster update scheme.
int Ising::Wolff(std::vector<int>& config) {
        double p = 1.0-std::exp(-2.0*beta); // Wolff's link activation probability. 
        std::random_device rd;
        std::mt19937 gen(rd());

        // Randomly choose the cluster seed.
        std::uniform_int_distribution<int> DisX(0, numSiteX-1);
        std::uniform_int_distribution<int> DisY(0, numSiteY-1);
        int x = DisX(gen);
        int y = DisY(gen);
        int seed = Coordinate(x, y);

        std::vector<int> cluster; // All the sites within the cluster. 
        std::vector<int> current; // Current sites to search ahead. 
        std::vector<int> forward; // Record the sites searhed to add to the cluster. 
        cluster.push_back(seed);
        current.push_back(seed);

        std::uniform_real_distribution<long double> Prob(0.0, 1.0);
        while (!current.empty()) { // If forward (now has been copied to current) is empty, then terminate the procedure.
            while (!current.empty()) { // Search around all the current loops.
                int cx = current.back() % numSiteX;
                int cy = current.back() / numSiteX;
                int cSpin = config[current.back()];

                // Search the nearest neighbor four sites.
                for (int i = 0; i < 2; ++i) { // Search along x-direction. 
                    int fx = (cx+(2*i-1)+numSiteX) % numSiteX;
                    int fy = cy;
                    int fSite = Coordinate(fx, fy);
                    double r = Prob(gen);
                    if (r < p && config[fSite] == cSpin) {
                        auto it = std::lower_bound(cluster.begin(), cluster.end(), fSite); // Return the first iterator which does not compare less than fSite.
                        if (*it != fSite) {
                            cluster.insert(it, fSite); // Insert before the specified position.
                            forward.push_back(fSite);
                            }
                        }
                    }
                for (int i = 0; i < 2; ++i) { // Search along y-direction. 
                    int fx = cx;
                    int fy = (cy+(2*i-1)+numSiteY) % numSiteY;
                    int fSite = Coordinate(fx, fy);
                    double r = Prob(gen);
                    if (r < p  && config[fSite] == cSpin) {
                        auto it = std::lower_bound(cluster.begin(), cluster.end(), fSite);
                        if (*it != fSite) {
                            cluster.insert(it, fSite);
                            forward.push_back(fSite);
                            }
                        }
                    }
                current.pop_back(); // Pop the last element until it is empty. 
                }
            current.insert(current.begin(), forward.begin(), forward.end());
            forward.clear();
            }

        // Perform simple update if cluster is not constructed, or flip the cluster with certainty.
        if (1 == cluster.size()) {
            int ns = config[Coordinate((x+1) % numSiteX, y)]+config[Coordinate((x-1+numSiteX) % numSiteX, y)]+config[Coordinate(x, (y+1) % numSiteY)]+config[Coordinate(x, (y-1+numSiteY) % numSiteY)];
            double cost = 2.0*config[ns]*ns;
            double rr = Prob(gen);
            if (rr < std::exp(-1.0*beta*cost)) { config[seed] = -1*config[seed]; }
            }
        else { for (auto it = cluster.begin(); it != cluster.end(); ++it) { config[*it] = -1*config[*it]; } }
        return cluster.size();
        }

// Evolve the config by num steps.
void Ising::MkvEvolve(std::vector<int>& config, int num) {
        for (int i = 0; i < num; ++i) { SimpleUpdate(config); }
        }

double Ising::EnergyDensity(std::vector<int> config){
        double e = 0.0;
        for (int x = 0; x < numSiteX; ++x) {
            for (int y = 0; y < numSiteY; ++y) {
                int site = Coordinate(x, y);
                int sitex = Coordinate((x+1) % numSiteX, y);
                int sitey = Coordinate(x, (y+1) % numSiteY);
                e += -1.0*(config[site]*config[sitex]+config[site]*config[sitey]);
                }
            }
        return e / config.size();
        }

double Ising::Magnetization(std::vector<int> config) { return abs(double(std::accumulate(config.begin(), config.end(), 0.0)) / config.size()); }

double Ising::MagnetizationQuad(std::vector<int> config) { return std::pow(Magnetization(config), 2.0); }

double Ising::MagnetizationBiquad(std::vector<int> config) { return std::pow(Magnetization(config), 4.0); }

// double Ising:MagSuseptibility()
void Ising::PrintConfig(std::vector<int> config) {
        for (int i = 0; i < numSiteX; ++i) {
            for (int j = 0; j < numSiteY; ++j) {
                int s = j*numSiteX+i;
                std::cout << config[s] << " ";
                }
            std::cout << std::endl;
            }
        std::cout << std::endl;
        }
