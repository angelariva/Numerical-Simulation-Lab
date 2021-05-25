#ifndef _MOLECULAR_DYNAMICS_H_
#define _MOLECULAR_DYNAMICS_H_
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <map>
#include <iostream>
#include "random.h"
#include "library.h"

class MolecularDynamics {
  private:
    // thermodynamical state
    unsigned int npart; // # of particles in the Simulation
    double energy, temp, vol, rho, box, rcut; // physical parameters

    // simulation parameters:
    unsigned int nstep, iprint;
    double delta;

    Random rand;

    // time interval that separates each measurement:
    // compute block size (for the blocking averages) and setting block index
    // (L = M/N of the previous exercises, M is nstep/measure_time_interval, N is
    // nblocks, L is block_size)
    unsigned int measure_time_interval;
    // size of each block
    unsigned int block_size;
    // actual filling block (counterss)
    unsigned int imeasure, iblock;
    // # of blocks
    unsigned int n_blocks;
    //output streams
    std::ofstream Epot, Ekin, Etot, Temp, Press, Gave, Gerr;

    // configuration:
    // positions, old positions, velocities, forces acting on each particle
    std::vector<double> x, y, z, xold, yold, zold, vx, vy, vz, fx, fy, fz;
    // observable properties:
    double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
    // vectors to get the blocking average:
    std::vector<double> est_pot, est_kin, est_etot, est_temp, est_press;

    // to compute g(r) (you need histograms):
    // for blocking method averages:
    std::vector<std::vector <double> > histo;
    const unsigned int nbins = 100;
    double bin_size;
    //maps for energy and pressure
    std::vector<std::string> keys;
    std::map<std::string, double> walker;
    std::map<std::string, double> block_average;
    std::map<std::string, double> global_average, global_average2;


    // functions: I declare them as private because I will only call them
    // from methods inside the class itself
    void Input(std::string);
    // computes the forces on the particles (from the potential):
    void Force();
    // computes Periodic Boundary Conditions for box of length L=box
    double Pbc(double) const;
    void Move();  // Moves the particles with Verlet algorithm
    void ConfFinal(std::string filename) const;   // writes the final configuration
    void ConfXYZ(int) const;  // writes the configuration in .xyz format
    void Measure(); // gives values of properties (on file!!)
    // computes and prints averages and statistical uncertaintes with blocking method
    void BlockingResults();

  public:

    //CONSTRUCTOR #1: without old positions
    // - 'simParameters' is the filename for input parameters
    // - 'configFile' is the filename for initial molecular configurations
    MolecularDynamics(std::string simParameters, std::string configFile);
    //CONSTRUCTOR #2: with old positions
    // - 'simParameters' is the filename for input parameters
    // - 'configFile' is the filename for initial molecular configurations,
    //    usually is the last configuration of a previous simulation.
    // - 'oldConfigFile' is the filename for the molecular configuration
    //    prior to 'configFile': it is used to extrapolate velocities:
    //
    // (1) With both of the constructors we obtain  r(t), r(t-dt);
    // (2) We then compute r(t+dt) with one step of the Verlet algorithm
    //     with the function Move();
    // (3) rold=r(t), r=r(t+dt). We can infere v(t+dt/2) and obtain the corresponding
    //     temperature T(t+dt/2): actual T = <v^2>/3,
    //     that we will use to rescale the velocity.
    // (4) we now use the rescaled velocity to estimate a new OLD
    //     spatial config: r_new(t) =  r(t+dt) - dt*v_scaled
    //     To do that we use the PBC algorithm (See function above)

    MolecularDynamics(std::string inputFile, std::string configFile,
        std::string oldConfigFile);

    // The simulation uses VERLET ALGORITHM;
    void RunSimulation();

    ~MolecularDynamics();

};



#endif //_MOLECULAR_DYNAMICS_H_
