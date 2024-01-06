// This is the main file for 2D ising project
//  file: 2d_ising.cpp
// 
//  Program to explore aspects of Monte Carlo simulation with the
//   Metropolis algorithm using the two-dimensional Ising model.
//
//  Programmers:  Youwei Liu            liu.9639@osu.edu
//
//  Revision history:
//      30-Apr-2023  original version 
//
//  Notes:
//   * in order to control the varible, the same configuration is used over all
//     concentration of Disorder class
//
//******************************************************************

#include <iostream>		// cout and cin
#include <iomanip>		// manipulators like setprecision
#include <fstream>		// file input and output
#include <string>

using namespace std;

#include "Disorder.h"

// function prototypes
extern unsigned long int random_seed ();	// routine to generate a seed

string filename(double concentration)
{
  ostringstream my_filename_stream;
  my_filename_stream.str(""); // clear the string stream
  my_filename_stream << "Disorder_with_C=" << fixed << setprecision(3) << concentration
                     << "_test.dat";
  return (my_filename_stream.str());
}


int main(void)
{

    ofstream Mout("magnet.dat"); // open the output file
    
    ofstream energy_out("energy_output.dat");

    //A simple limit check to check if the magnet class work finely, at T->0 this gives about 1 or -1 magnetization
    // at T-> inf this gives magnetization about 0, the magnetization distribute randomly
    // Magnet my_magnet(10);
    // my_magnet.set_kT(0.01);
    // cout<<"the magnetization is :"<<endl;
    // my_magnet.print();
    // cout<<" E= " <<my_magnet.get_energy()<<"    M= "<<my_magnet.get_magnetization()<<endl;
    // my_magnet.flip_N(2000);
    // cout<<"the magnetization is :"<<endl;
    // my_magnet.print();
    // cout<<" E= " <<my_magnet.get_energy()<<"    M= "<<my_magnet.get_magnetization()<<endl;
    // my_magnet.set_kT(1000);
    // my_magnet.flip_N(2000);
    // cout<<"the magnetization is :"<<endl;
    // my_magnet.print();
    // cout<<" E= " <<my_magnet.get_energy()<<"    M= "<<my_magnet.get_magnetization()<<endl;
    // Magnet A=my_magnet;

    //size of the grid
    int my_size=50;
    int total_flip=1000000;
    //the step defines how often it read the data into file
    int flip_step=10000;
    string myfilename[6];
    //concentration of the code looping through
    double concentration[6]={0.005,0.01,0.03,0.1,0.2,0.4};
    ofstream Dout[6];
    //below code used because there is no default constructor defined for class disorder
    Disorder my_disorder[6] ={Disorder(my_size,"nonmagnetic",concentration[0]),Disorder(my_size,"nonmagnetic",concentration[1]),Disorder(my_size,"nonmagnetic",concentration[2]),Disorder(my_size,"nonmagnetic",concentration[3]),Disorder(my_size,"nonmagnetic",concentration[4]),Disorder(my_size,"nonmagnetic",concentration[5])};

    // initialize all the classes and print files
    Magnet my_mag(my_size);
    my_mag.set_kT(0.01);
    energy_out<<"energy of the configurations with given concentration of disorder"<<endl<<"          0"<<"          ";
    for(int i=0;i<6;i++)
    {
        my_disorder[i].set_kT(0.01);
        my_disorder[i]=my_mag;
        myfilename[i]=filename(concentration[i]);
        Dout[i].open(myfilename[i]); // open the output file
        Dout[i]<<"spin configurations of Magnet from high T to low T, Disordered Magnet"<<endl;
        energy_out <<concentration[i]<<"          ";
    }



    //print all the magnet's initial state
    energy_out<<endl<<fixed<<setw(9)<<"0" <<fixed<<setprecision(1)<<setw(9)<<my_mag.get_energy();
    Mout<<"spin configurations of Magnet from high T to low T, Magnet"<<endl;
    for (int k=0;k<6;k++)
    {
        for (int i=0;i<my_size;i++){
            for (int j=0;j<my_size;j++)
            {
                Mout<<my_mag.get_spin(i,j)<< " ";
                Dout[k]<<my_disorder[k].get_spin(i,j)<< " ";
            }
                Mout<<endl;
                Dout[k]<<endl;
        }
        energy_out<<fixed<<setw(9)<<my_disorder[k].get_energy();
        Mout<<endl;
        Dout[k]<<endl;

    }
    energy_out<<endl;

    //Monte Carol for evolving spins, and print out the energy and configuration in desinate files
    for (int flip=0;flip<total_flip;flip+=flip_step)
    {
        energy_out<<fixed<<setw(9)<<flip<<fixed<<setprecision(1)<<setw(9)<<my_mag.get_energy();
        for (int k=0; k<6;k++)
        {

            my_mag.flip_N(flip_step);
            my_disorder[k].flip_N(flip_step);
            for (int i=0;i<my_size;i++){
                for (int j=0;j<my_size;j++)
                {
                    Mout<<my_mag.get_spin(i,j)<< " ";
                    Dout[k]<<my_disorder[k].get_spin(i,j)<< " ";
                }
                

                Mout<<endl;
                Dout[k]<<endl;
                }
            Mout<<endl;
            Dout[k]<<endl;
            energy_out<<fixed<<setw(8)<<my_disorder[k].get_energy();
        }
        energy_out<<endl;
    }



    
    

    Mout.close();
    for (int i=0;i<6;i++)
    {
        Dout[i].close();
    }

    energy_out.close();
    


    return (0);
}