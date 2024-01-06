//
// This is the header file for derived disorder class
//  file: Disorder.h
//        Disorder.cpp
// 
//  Program defines the Disorder class with the information of where is the disorder, what types is it
//  and as well change some functions defination because the introduction of disorder into the system
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

//
#include <string>
#include <iostream>		// cout and cin

using namespace std;		// we need this when .h is omitted

#include "Magnet.h"

#include <gsl/gsl_rng.h>     // GSL random number generators
#include <gsl/gsl_randist.h>    // for generating random numbers from distributions


class Disorder: public Magnet
{
    public:
        Disorder(int my_size, string my_type, double concentrition); //constructor of disorder
        Disorder& operator=(const Magnet& other);//used to control intial condition

        //my_size is the size of the grid, my_type only have "nonmagnetic availble", concentrition is the disorder concentration
        void set_type(string my_type);
        void set_max_tries(int N){max_tries=N;};// set the maximum tries before gives a error because of no spin to flip

        int get_spin(int x, int y) const override;// return -1,1,0, -1 and 1 for down and up spin, 0 for disorder
        bool get_disorder(int x, int y) const {return position.get(x,y);};//getter function for disorder, 0 do not have disorder, 1 is disorder site

        void print_type();//what king of disorder is it?
        void print() const override;//"N"print for where every position gives 1, configurations value if position gives 0 at the site
        void print_disorder() const;//print position out with 0 print for 0, "N" print for value 1

        //double cal_delta_E(int x, int y) const;
        void single_flip() override; //random flip one spin at the site excluding the non_magnetic disorder site, or report error if do not find one for at max the max_tries number



    private:
        int max_tries;//the maximum try to flip an flippable spin before reporting an error
        string type; //so far only non-magnetic availible
        Matrix<bool> position; //matrix store postion of the disorder, 1 for disorder, 0 for normal magnetic material
};


