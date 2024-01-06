//
// This is the implementation file for derived disorder class
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


#include "Disorder.h"


Disorder::Disorder(int my_size, string my_type, double concentrition):Magnet(my_size)
{
    position.re_size(my_size,my_size);
    if (my_type == "nonmagnetic" ) {
        type = my_type;
        max_tries=50;//because for other type of disorder we do not have the problem of no spin to flip
    }
    else {
        cerr << "Invalid type! Please choose from: nonmagnetic" << endl;
    }

    // Create an instance of the GSL RNG
    rng = gsl_rng_alloc(gsl_rng_default);
    // Seed the RNG
    gsl_rng_set(rng, seed);
    // Set the configuration matrix elements randomly
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // Generate a random number from Bernoulli distribution with p=0.5
            bool random_val = gsl_ran_bernoulli(rng, concentrition);
            position.set(i, j, random_val);    //random set the value of the configuration at given site to be 0 or 1
        }
    }
    //the code initialize the magnet before it generate the disorder so it must be recalculate here
    set_magnetization(calculate_magnetization());
    set_energy(calculate_energy());

}

Disorder& Disorder::operator=(const Magnet& other) {
    if(this != &other) {
        // Copy the base class members
        Magnet::operator=(other);

        // Copy any additional members in the Magnet class that are not part of the base class
        // ...
    }
    return *this;
}


void Disorder::set_type(string my_type) {
    // if (my_type == "nonmagnetic" || my_type == "ferromagnetic" || my_type == "antiferromagnetic") {
    //     type = my_type;
    // }
    // else {
    //     cerr << "Invalid type! Please choose from: nonmagnetic, ferromagnetic, antiferromagnetic." << endl;
    // }
    if (my_type == "nonmagnetic" ) {
        type = my_type;
    }
    else {
        cerr << "Invalid type! Please choose from: nonmagnetic" << endl;
    }
}

void Disorder::print_type() {
    cout << "Type of disorder: " << type << endl;
}

//if the position at site gives 0, then output the original spin, if position at site gives 1, then output 0 for non magnetic disorder 
int Disorder::get_spin(int x, int y) const {
    int spin = Magnet::get_spin(x, y); // get the original spin value
    if (position.get(x, y)) { // if position at site gives 1
        if (type == "nonmagnetic") {
            spin = 0; // output 0 for non-magnetic disorder
        }
    }
    return spin;
}


void Disorder::print() const
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (position.get(i, j))
            {
                if (type == "nonmagnetic")
                {
                    std::cout << "N";
                }
            }
            else
            {
                std::cout << get_value(i, j);
            }
            std::cout << " ";
        }
        cout << endl;
    }
}

void Disorder::print_disorder() const
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
           if (position.get(i, j))
            {
                if (type == "nonmagnetic")
                {
                    std::cout << "N ";
                }
            }
            else
            {
                std::cout <<  "0 ";
            }
        }
        cout << endl;
    }
}

void Disorder::single_flip()
{
    int tries = 0;
    int x, y;
    double P_flip, rand_num;
    bool valid_flip = false;

    // Try flipping a random spin until a valid spin is found or max_tries is exceeded
    while (!valid_flip && tries < max_tries) {
        x = gsl_rng_uniform_int(rng, size);
        y = gsl_rng_uniform_int(rng, size);

        if (get_disorder(x, y) == 0) {
            P_flip = min(1., exp(-cal_delta_E(x,y)/k_T));
            rand_num = gsl_rng_uniform(rng);

            if (rand_num <= P_flip) {
                flip(x,y);
            }
            valid_flip = true;
        }
        tries++;
    }

    if (!valid_flip) {
        cerr << "Could not find a valid spin to flip after " << max_tries << " tries." << endl;
    }
}
