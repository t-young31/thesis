//
//  Eigenvalues of the radial part of a 3D exponential well
//
// compile with:
//     g++ -std=c++11 exp_well_3d.cpp -o exp_well_3d -I$HOME/miniconda3/include/ -O3
//
// using the a Miniconda install of boost for the headers

#include <cstdio>
#include <iostream>
#include <cmath>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< long double , 2 > vec_type;
typedef runge_kutta_cash_karp54< vec_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
typedef runge_kutta4< vec_type > stepper;


struct ode{

    long double m;
    long double energy;
    long double a;
    long double k;
    long double l;     // angluar quantum number

    void operator()( const vec_type &vec , vec_type &dvecdx , double r)
    {
        dvecdx[0] = vec[1];
        dvecdx[1] =  -2.0*vec[1]/r - (2*m*(energy - k*(exp(a*r) - 1.0)) - (l*(l+1)/(r*r)))*vec[0];
    }
};

struct nothing{

    void operator()( const vec_type &vec , const double r){

    }
};


void print_exp_well_psi_final(double k, double a, double m){
    controlled_stepper_type controlled_stepper;

    // Energy sweep range
    double min_energy = 0.000001;
    double max_energy = 10 * 0.00094;    // kBT = 0.00094 Ha

    int l_max = 50;                       // Maximum value of the l qunatum number

    double r_init = 0.001;               // Where to start the propogation from
    double r_final = 1.7;                // Where to end   '   '  as there is a divergence at r=0

    int n_steps = 20000;


    for (int l = 0; l <= l_max; l++ ){

        double prev_psi_final = 0.0;         // Final value of Psi i.e. Psi(r_final) for each eigenvalue


        for (int i = 0; i < n_steps; i++){

            double energy = min_energy + ((max_energy - min_energy) * i) / n_steps;

            vec_type vec = { 1.0 , 1E-20 }; // initial conditions
            integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-10 ),
                               ode {m , energy, a, k, static_cast<long double>(l)},
                               vec,
                               r_init,
                               r_final,
                               1E-10,
                               nothing {}
                               );

            // If the product of the previous Psi(r_final) and the current one is negative
            // there must be a zero between them, i.e. a normalaisable WF
            if (prev_psi_final * vec[0] < 0){
                // cout << energy << endl;
                cout << energy << '\t' << l << endl;

            }

            // Reset the previous Psi(r_final)
            prev_psi_final = vec[0];

        } // Energy
    }  // l
}

int main(int argc, char **argv){

    freopen("eigenvals_co2.txt","w", stdout);
    cout << "Eigenvalue / Ha    l number" << endl;

    double k_co2 = 0.00220873;                   // prefactor for the exponential well

    double k = 0.0001;
    double a = 2.832 * 0.529177;             // exponent in the expoential well
    double m = 48 * 1822.888486;             // mass in atomic units from amu

    print_exp_well_psi_final(k_co2, a, m);

    return 0;
}
