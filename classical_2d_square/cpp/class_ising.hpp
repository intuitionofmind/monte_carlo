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
        int site = Coordinate(x, y);

        std::vector<int> cluster; // All the sites within the cluster. 
        std::vector<int> current; // Current sites to search ahead. 
        std::vector<int> forward; // Record the sites searhed to add to the cluster. 
        cluster.push_back(site);
        current.push_back(site);
        // PrintConfig(config);
        std::uniform_real_distribution<long double> Prob(0.0, 1.0);
        while (!current.empty()) { // If forward is empty, terminate the loop. 
            while (!current.empty()) {
                int cx = current.back() % numSiteX;
                int cy = current.back() / numSiteX;
                int cspin = config[current.back()];

                // Search the nearest neighbor four sites.
                for (int i = 0; i < 2; ++i) { // Search along x-direction. 
                    int fx = (cx+(2*i-1)+numSiteX) % numSiteX;
                    int fy = cy;
                    int fsite = Coordinate(fx, fy);
                    int fspin = config[fsite];
                    double r = Prob(gen);
                    if (r < p && std::find(cluster.begin(), cluster.end(), fsite) == cluster.end() && fspin == cspin) {
                        cluster.push_back(fsite);
                        forward.push_back(fsite);
                        }
                    }
                for (int i = 0; i < 2; ++i) { // Search along y-direction. 
                    int fx = cx;
                    int fy = (cy+(2*i-1)+numSiteY) % numSiteY;
                    int fsite = Coordinate(fx, fy);
                    int fspin = config[fsite];
                    double r = Prob(gen);
                    if (r < p && std::find(cluster.begin(), cluster.end(), fsite) == cluster.end() && fspin == cspin) {
                        cluster.push_back(fsite);
                        forward.push_back(fsite);
                        }
                    }
                current.pop_back(); // Pop the last element until it is empty. 
                }
            current.insert(current.begin(), forward.begin(), forward.end());
            forward.clear();
            }

        // Perform simple update if cluster is not constructed, or flip the cluster with certainty.
        if (1 == cluster.size()) {
            int ns = config[Coordinate(x, (y+1) % numSiteY)]+config[Coordinate(x, (y-1+numSiteY) % numSiteY)]+config[Coordinate(x, (y+1) % numSiteY)]+config[Coordinate(x, (y-1+numSiteY) % numSiteY)];
            double cost = 2.0*config[site]*ns;
            double rr = Prob(gen);
            if (rr < std::exp(-1.0*beta*cost)) { config[site] = -1*config[site]; }
            }
        else { for (std::vector<int>::iterator it = cluster.begin(); it != cluster.end(); ++it) { config[*it] = -1*config[*it]; } }
        return cluster.size();
        }

// Evolve the config by num steps.
void Ising::MkvEvolve(std::vector<int>& config, int num) {
        for (int i = 0; i < num; ++i) {
            Wolff(config);
            }
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

double Ising::Magnetization(std::vector<int> config) { return double(std::accumulate(config.begin(), config.end(), 0.0)) / config.size(); }

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
