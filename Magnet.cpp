//
// This is the implementation file for magnet class
//  file: Magnet.h
//        Magnet.cpp
// 
//  Program defines the Matrix class(equalvalent with Matrix in Physics)
//  and the magnet class which used to store the non_doped system, or the base 
//  class for a doped system
//
//  Programmers:  Youwei Liu            liu.9639@osu.edu
//
//  Revision history:
//      30-Apr-2023  original version 
//
//  Notes:
//   * all the override of = (assign to) operation is used to obtain the same starting configuration
//
//******************************************************************
#include "Magnet.h"


// function prototypes
extern long int random_seed ();	// routine to generate a seed

//********************************************************************
// Matrix class body:

// Constructor
template<typename T>
Matrix<T>::Matrix() {
    rows = 5;
    cols = 5;
    data.resize(rows);
    for (int i = 0; i < rows; i++) {
        data[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            data[i][j] = 0;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(int size) {
    rows = size;
    cols = size;
    data.resize(rows);
    for (int i = 0; i < rows; i++) {
        data[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            data[i][j] = 0;
        }
    }
}

template<typename T>
Matrix<T>::~Matrix() {}

template<typename T>
T Matrix<T>::get(int r, int c) const {
    return data[r][c];
}

template<typename T>
void Matrix<T>::set(int r, int c, T val) {
    data[r][c] = val;
}

template<typename T>
void Matrix<T>::re_size(int r, int c) {
    rows = r;
    cols = c;
    data.resize(rows);
    for (int i = 0; i < rows; i++) {
        data[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            data[i][j] = 0;
        }
    }
}

template<typename T>
void Matrix<T>::print() const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    // Check for self-assignment
    if (this == &other) {
        return *this;
    }
    // Resize this matrix to match the size of the other matrix
    re_size(other.rows, other.cols);
    // Copy the values from the other matrix to this matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            data[i][j] = other.data[i][j];
        }
    }
    // Return a reference to this matrix
    return *this;
}

//********************************************************************
// Magnet class body:

// Constructor
Magnet::Magnet(int s)
{
        // Set the size of the magnet
    size = s;
    seed=random_seed();
    J_ising=1.;
    k_T=2.;
    configurations.re_size(s,s);
    // Create an instance of the GSL RNG
    rng = gsl_rng_alloc(gsl_rng_default);
    // Seed the RNG
    gsl_rng_set(rng, seed);
    // Set the configuration matrix elements randomly
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // Generate a random number from Bernoulli distribution with p=0.5
            bool random_val = gsl_ran_bernoulli(rng, 0.5);
            configurations.set(i, j, random_val);    //random set the value of the configuration at given site to be 0 or 1
        }
    }
    energy=calculate_energy();
    magnetization=calculate_magnetization();

}

Magnet::Magnet(int s,double J,double T)
{
        // Set the size of the magnet
    size = s;
    seed=random_seed();
    J_ising=J;
    k_T=T;
    configurations.re_size(s,s);
    // Create an instance of the GSL RNG
    rng = gsl_rng_alloc(gsl_rng_default);
    // Seed the RNG
    gsl_rng_set(rng, seed);
    // Set the configuration matrix elements randomly
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // Generate a random number from Bernoulli distribution with p=0.5
            bool random_val = gsl_ran_bernoulli(rng, 0.5);
            configurations.set(i, j, random_val);
        }
    }
    energy=calculate_energy();
    magnetization=calculate_magnetization();

}

Magnet::Magnet(int s,double J,double T, long int my_seed)
{
    // Set the size of the magnet
    size = s;
    seed=my_seed;
    J_ising=J;
    k_T=T;
    // Create an instance of the GSL RNG
    rng = gsl_rng_alloc(gsl_rng_default);
    // Seed the RNG
    gsl_rng_set(rng, my_seed);
    // Set the configuration matrix elements randomly
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // Generate a random number from Bernoulli distribution with p=0.5
            bool random_val = gsl_ran_bernoulli(rng, 0.5);
            configurations.set(i, j, random_val);
        }
    }
    energy=calculate_energy();
    magnetization=calculate_magnetization();

}; 

// Destructor
Magnet::~Magnet() {    
    // Free the memory associated with the RNG
    gsl_rng_free(rng);}

// Getter function
int Magnet::get_spin(int x, int y) const
{
    // value *2-1, so for 0 input gives -1, 1 input gives 1
    return(int((configurations.get(x,y)))*2.-1.);
}

// printer function
void Magnet::print() const
{
    configurations.print();
}

// Setter function
void Magnet::set_value(int x, int y,bool value) {
    configurations.set(x,y,value);
}


//Calculate the potential energy change at given sight to flip a spin(not including changing the nearest neighbor energy)
double Magnet::cal_delta_E(int x, int y) const
{
    double delta_E=0;
    //the 2 here came from |fliped spin - original spin| =2
    delta_E+=-2. * J_ising * double(-get_spin(x,y) * get_spin((x-1+size)%size,y));
    delta_E+=-2. * J_ising * double(-get_spin(x,y) * get_spin((x+1)%size,y));
    delta_E+=-2. * J_ising * double(-get_spin(x,y) * get_spin(x,(y+1)%size));
    delta_E+=-2. * J_ising * double(-get_spin(x,y) * get_spin(x,(y-1+size)%size));
    return(delta_E);
}

//special setter function
void Magnet::flip(int x, int y)
{
    bool value = get_value(x,y);
    //the 2 here exist due to flip a spin not only change the energy of the energy count from that specific spin, 
    //but also its nearest neighbors, (because in their energy sum, that spin will be flipped), and all four nearest neighbor,
    //the amount of energy change overall will equal to delta_E, so it pick up a 2 here
    energy += 2.*cal_delta_E(x,y);
    configurations.set(x,y,!value);
    //the magnetization change will be (spin flipped - orignal spin)/whole size of configuration
    magnetization+=2*get_spin(x,y)/(double(size)*double(size));
}

void Magnet::single_flip()
{
    int x = gsl_rng_uniform_int(rng, size); // generate random integer coordinate between 0 and N-1
    int y = gsl_rng_uniform_int(rng, size);

    double P_flip = min(1., exp(-cal_delta_E(x,y)/k_T));// notice if this lower the energy it will 100% flip, otherwise, it have exp() possiblilty to flip
    double rand_num = gsl_rng_uniform(rng); // generate random number between 0 and 1
    //cout<<"P_flip = "<<P_flip<<"      rand_num = "<< rand_num<<endl<<-cal_delta_E(x,y);
    if (rand_num <= P_flip) {
        //if the random num smaller than the probility then flip, otherwise, do nothing

        flip(x,y);
    }
    else {
        return;
    }
};

void Magnet::flip_N(int N)
{
    for (int i=0;i<N;i++)
    {
        single_flip();
    }
};

//************************* calculate_energy ******************************
//
//  Given the array of integers configuration[0...num_sites-1], which 
//   specifies the spin at each lattice point, find the energy of that 
//   configuration [eq.(13.6) in Session 13 notes].
//  Note that a freeboundary condition is specified here.
//
//*************************************************************************
// this is a total check only use to intialize the energy and check the update 
// energy is correct or not, because this calculation can waste a lot of time for 
// a large size matrix
double 
Magnet::calculate_energy () const
{
  double energy_sum = 0.;

  // go through the 2d lattice with free boundary conditions
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        //convert bool to actual spin
        //note the boundary condition have already been taken account by taking the mode of the size, 
        double right =double(get_spin(i,(j+1)%size));
        double center = double(get_spin(i,j));
        double down = double(get_spin((i+1)%size,j));

        // x direction
        energy_sum += - J_ising * (center*right);

        // y direction;
        energy_sum += - J_ising * (center*down);
      }
    }

  return (energy_sum);
} 


// this is a total check only use to intialize the magnetization and check the update 
// magnetization is correct or not, because this calculation can waste a lot of time for 
// a large size matrix and it is unneccessary
double Magnet::calculate_magnetization() const
{
    int spin_sum = 0;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        spin_sum+=get_spin(i,j);
      }
    }
    return(double(spin_sum)/(double(size)*double(size)));
}


Magnet& Magnet::operator=(const Magnet& other)
{
    //assign everything is the same, but the random pointer
    if (this != &other) {
        size = other.size;
        J_ising = other.J_ising;
        k_T = other.k_T;
        energy = other.energy;
        magnetization = other.magnetization;
        configurations = other.configurations;
    }
    return *this;
}