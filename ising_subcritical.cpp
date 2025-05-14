// ising_subcritical.cpp
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <chrono>

class IsingModel {
private:
    std::vector<int8_t> grid;  // Using int8_t to reduce memory usage
    std::vector<bool> cluster_map;  // Persistent cluster map to avoid reallocation
    int N;
    double beta;
    std::vector<std::mt19937> rngs;
    std::vector<std::uniform_real_distribution<double> > uniform_dists;

    inline size_t idx(int i, int j) const {
        return static_cast<size_t>(i) * N + j;
    }

    inline int8_t get_spin(int i, int j) const {
        i = (i + N) % N;
        j = (j + N) % N;
        return grid[idx(i, j)];
    }

public:
    IsingModel(int size, double temperature, int seed = 42) 
        : N(size), beta(1.0/temperature) {
        
        // Safely allocate grid
        try {
            std::cout << "Allocating " << (N * N * sizeof(int8_t) / (1024.0 * 1024.0)) << " MB for grid..." << std::endl;
            grid.resize(static_cast<size_t>(N) * N);
            cluster_map.resize(static_cast<size_t>(N) * N);
        } catch (const std::bad_alloc& e) {
            std::cerr << "Failed to allocate memory" << std::endl;
            throw;
        }

        // Initialize RNGs
        int num_threads = omp_get_max_threads();
        rngs.resize(num_threads);
        uniform_dists.resize(num_threads, std::uniform_real_distribution<double>(0.0, 1.0));
        
        for(int i = 0; i < num_threads; i++) {
            rngs[i].seed(seed + i);
        }
        
        // Initialize grid in chunks
        const int chunk_size = 1000;
        for(int chunk = 0; chunk < N; chunk += chunk_size) {
            int chunk_end = std::min(chunk + chunk_size, N);
            #pragma omp parallel for collapse(2)
            for(int i = chunk; i < chunk_end; i++) {
                for(int j = 0; j < N; j++) {
                    int thread_id = omp_get_thread_num();
                    grid[idx(i, j)] = (uniform_dists[thread_id](rngs[thread_id]) < 0.5) ? 1 : -1;
                }
            }
        }
    }

    void wolff_step() {
        int thread_id = 0;
        
        // Pick random starting point
        int i0 = static_cast<int>(uniform_dists[thread_id](rngs[thread_id]) * N);
        int j0 = static_cast<int>(uniform_dists[thread_id](rngs[thread_id]) * N);
        
        int8_t cluster_spin = grid[idx(i0, j0)];
        double p_add = 1 - exp(-2 * beta);
        
        // Reset cluster map
        std::fill(cluster_map.begin(), cluster_map.end(), false);
        
        // Use fixed-size vectors for the stack to prevent memory issues
        std::vector<std::pair<int, int> > stack;
        stack.reserve(1000000);  // Reserve reasonable size
        
        stack.push_back(std::make_pair(i0, j0));
        cluster_map[idx(i0, j0)] = true;
        
        size_t max_cluster_size = static_cast<size_t>(N) * N / 4;  // Limit cluster size
        size_t cluster_size = 0;
        
        while (!stack.empty() && cluster_size < max_cluster_size) {
            auto current = stack.back();
            stack.pop_back();
            int i = current.first;
            int j = current.second;
            
            size_t current_idx = idx(i, j);
            grid[current_idx] *= -1;  // Flip the spin
            cluster_size++;
            
            // Check neighbors
            const int di[] = {1, -1, 0, 0};
            const int dj[] = {0, 0, 1, -1};
            
            for(int d = 0; d < 4; d++) {
                int ni = (i + di[d] + N) % N;
                int nj = (j + dj[d] + N) % N;
                size_t nidx = idx(ni, nj);
                
                if (!cluster_map[nidx] && 
                    grid[nidx] == cluster_spin && 
                    uniform_dists[thread_id](rngs[thread_id]) < p_add) {
                    stack.push_back(std::make_pair(ni, nj));
                    cluster_map[nidx] = true;
                }
            }
        }
    }

    void metropolis_step() {
        for(int color = 0; color < 2; color++) {
            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                #pragma omp for collapse(2)
                for(int i = 0; i < N; i++) {
                    for(int j = 0; j < N; j++) {
                        if((i + j) % 2 == color) {
                            size_t pos = idx(i, j);
                            int8_t spin = grid[pos];
                            int sum = get_spin(i+1, j) + get_spin(i-1, j) +
                                    get_spin(i, j+1) + get_spin(i, j-1);
                            
                            double dE = 2.0 * spin * sum;
                            if(dE < 0 || uniform_dists[thread_id](rngs[thread_id]) < exp(-beta * dE)) {
                                grid[pos] *= -1;
                            }
                        }
                    }
                }
            }
        }
    }

    void hybrid_equilibrate(int steps) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Use a mix of Wolff and Metropolis steps based on temperature
        double T = 1.0/beta;
        double T_c = 2.27;
        double wolff_ratio;
        
        if (T < T_c) {
            wolff_ratio = 0.8;  // Mostly Wolff for subcritical
        } else if (T > T_c) {
            wolff_ratio = 0.2;  // Mostly Metropolis for supercritical
        } else {
            wolff_ratio = 0.5;  // Equal mix at critical point
        }
        
        for(int s = 0; s < steps; s++) {
            if (uniform_dists[0](rngs[0]) < wolff_ratio) {
                wolff_step();
            } else {
                metropolis_step();
            }
            
            if(s % 100 == 0) {
                auto current_time = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
                double steps_per_second = (s + 1) / static_cast<double>(elapsed);
                double estimated_total_time = steps / steps_per_second;
                double remaining_time = estimated_total_time - elapsed;
                
                std::cout << "Step " << s << " of " << steps 
                         << " (" << (s * 100.0 / steps) << "%) "
                         << "ETA: " << (remaining_time / 3600.0) << " hours"
                         << " (Wolff ratio: " << wolff_ratio << ")\r" << std::flush;
            }
        }
        std::cout << std::endl;
    }

    void save_grid(const std::string& filename) {
        std::cout << "Saving full grid to " << filename << "..." << std::endl;
        std::ofstream out(filename, std::ios::binary);
        out.write(reinterpret_cast<char*>(grid.data()), grid.size());
    }
};

int main() {
    const int N = 10000;
    const double T_c = 2.27;
    const double T = 0.7 * T_c;  // Subcritical temperature
    const int n_steps = 100000;
    
    try {
        std::cout << "Running with " << omp_get_max_threads() << " threads" << std::endl;
        std::cout << "Initializing " << N << "x" << N << " Ising model..." << std::endl;
        std::cout << "Temperature: " << T << " (T/T_c = " << T/T_c << ")" << std::endl;
        
        IsingModel model(N, T);
        std::cout << "Starting hybrid Wolff/Metropolis algorithm..." << std::endl;
        model.hybrid_equilibrate(n_steps);
        
        std::cout << "Saving grid..." << std::endl;
        model.save_grid("ising_10000_sub.bin");
        
        std::cout << "Done! Use the Python plotting script to visualize results." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}