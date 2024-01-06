//
// This is the header file for magnet class
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
#ifndef MAGNET_H
#define MAGNET_H
#include <iostream>		// cout and cin
#include <vector>       // For matrix class
#include <algorithm>    // For function min() for determine the probability
#include <cmath>        // For expenential eqn

using namespace std;		// we need this when .h is omitted

#include <gsl/gsl_rng.h>     // GSL random number generators
#include <gsl/gsl_randist.h>    // for generating random numbers from distributions

template<typename T>
class Matrix {
private:
    int rows;
    int cols;
    vector<vector<T>> data;

public:
    //Matrix();
    Matrix();
    Matrix(int size);
    ~Matrix();

    // Getter functions
    T get(int r, int c) const;//get the value of matrix at position, [r,c]
    int get_r() const {return rows;};
    int get_c() const {return cols;};

    // Setter functions
    void set(int r, int c, T val);//set the value of matrix at [r,c] to T, notice, r,c range from 0~size-1
    void re_size(int r, int c); // resize the matrix to r rows, c column
    
    // operator= definition
    Matrix<T>& operator=(const Matrix<T>& other); //resize and assign all the element to its rhs matrix

    // printer
    void print() const;//print out the matrix
};


class Magnet
{
    public:
        Magnet(int s);
        Magnet(int s,double J, double T);                  // Constructor
        Magnet(int s,double J,double T, long int my_seed);                  // Constructor
        ~Magnet();                      // Destructor
        
        //getter:
        // note the x,y must be strictly smaller than size, it goes from 0 ~ (size-1) 
        // virtual because it will need to also return the info about the site has disorder or not
        virtual int get_spin(int x, int y) const;// get the spin of the configuration, this will print either -1 or 1
        double get_energy() const    {return energy;}; //
        double get_magnetization() const    {return magnetization;}; 
        //Matrix<int> get_spin_config() const;//for printing purpose

        //printer:
        //print out the configurations using 1,0,
        //1 for spin 1, 0 for spin -1
        //virtual because later we will need to print the disorder with special character
        virtual void print() const;

        //setter:
        void set_kT(double T) {k_T=T;};
        void set_J(double J) {J_ising=J;};


        virtual void flip_N(int N); //random flip one spin based on it's energy change N times

        //do a whole swipe to calculate energy or magnetization,(iterate through all configurations for calculation, 
        //designed for double check the energy calculation is right, and intialize the energy and magnetization count)
        // the virtual is unset due to we are only work with non_magnetic disorder, and the get spin already been updated
        // to give -1,1 for down up spin and 0 for disorder
        double calculate_energy() const;
        double calculate_magnetization() const;
        Magnet& operator=(const Magnet& other);

    protected:
        //all below class is designed for usage to generate user interface only, not avaible for others to use
        //Calculate the potential energy change at given sight to flip a spin(not including changing the nearest neighbor energy)
        //Delta_E is different if disorder is involved
        virtual double cal_delta_E(int x, int y) const; 
        void flip(int x, int y);       //flip the spin at the given position
        void set_value(int x, int y, bool value); //internal function, use to set the configuration to right bool#, should not be called by user
        //this is virtual becuase, in the case of high disorder, we want to only flip the spin that can be flipped
        virtual void single_flip(); //random flip one spin
        bool get_value(int x, int y) const    {return configurations.get(x,y);};     //return 0 for spin down, 1 for spin up
        void set_magnetization(double M){magnetization=M;};//setter for disorder class
        void set_energy(double E){energy=E;};//setter for disorder class
        //Matrix<bool> get_config() const{return configurations;};//for setting configuration purpose
        void set_config(Matrix<bool> my_config) {configurations=my_config;  size=my_config.get_c();};
        int size;//size of the grid
        gsl_rng* rng;
        long int seed;
        double k_T;//the temperature of current magnet happens

    private:
        double J_ising;//exchange interaction prameter
        double energy;//total energy of the system
        double magnetization;//average magnetic moment per spin
        Matrix<bool> configurations;// store all the spin information, 0 for down, 1 for up, using bool because save space
        //and it is easy to flip using !(bool)
        
};

#endif

