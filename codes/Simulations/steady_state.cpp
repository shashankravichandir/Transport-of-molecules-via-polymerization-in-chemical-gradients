// Program to generate the steady state profiles for active flexible chains with each monomer either being active or passive

#include <iostream>
#include <fstream>
#include <string>
#include "particle.h"
#include "random.h"

using namespace std;

int main(int argc, char** argv){
    seed_drand48();

    // Parameters
    double ks = 10.0, Dt = 1.0, Dr = 5.0, l0 = 0.0;
    double Lx = 100.0; double Ly = 100.0;

    //Intializing polymer chain
    int N = 2; 
    polymer p(N);
    vector<int> v{1};
    AssignActivity(p,v);
    for(int i=0; i<N; i++){
        p.chain[i].x = 50.0 - i*5.0;
        p.chain[i].y = 50.0 - i*5.0;
        p.chain[i].phi = 2*pi*uniform_rand();
    }
    polymer p_new = p;

    // Other variables
    double t = 0.0; double tmax = 10000000.0; double dt = 0.01; 
    double eta_x, eta_y, eta_phi;
    double l, Fx, Fy;

    // Averaging variables
    int timestep=0; double xc; int Nbins = 25; int count[Nbins]; double lbin = Lx/double(Nbins); int j; int timesteps_per_snapshot= 100;
    for(int k = 0; k<Nbins; k++){
        count[k] = 0;
    }
    
    // Opening files
    ofstream out1;
    out1.open("./Results/Densities/StateDiagramResults/N" + to_string(N) +"/alpha_1.dat");

    // Time loop
    
    while(t<tmax){
        timestep++;
        t = t+dt;

        for(int i = 0; i<N; i++){
            Fx = 0.0; Fy = 0.0;
            if(i!=0 && i!=N-1){
                l = length(p.chain[i], p.chain[i-1]);
                Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
                Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;

                l = length(p.chain[i], p.chain[i+1]);
                Fx += ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
                Fy += ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
            }
            else if(i==0){
                l = length(p.chain[i], p.chain[i+1]);
                Fx = ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
                Fy = ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
            }
            else if(i==N-1){
                l = length(p.chain[i-1], p.chain[i]);
                Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
                Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
            }
            
            eta_x = Gaussian();
            eta_y = Gaussian();
            eta_phi = Gaussian();

            p_new.chain[i].x += dt*(Fx + p.chain[i].alpha*fs(p.chain[i].x, Lx)*cos(p.chain[i].phi)) + sqrt(2*Dt*dt)*eta_x;
            p_new.chain[i].y += dt*(Fy + p.chain[i].alpha*fs(p.chain[i].x, Lx)*sin(p.chain[i].phi)) + sqrt(2*Dt*dt)*eta_y;
            p_new.chain[i].phi += sqrt(2*Dr*dt)*eta_phi;    
        }

        p.chain = p_new.chain;

        if(timestep%timesteps_per_snapshot == 0){
            xc = CentreOfMass(p);
            while(xc<0) {xc += Lx;}
            while(xc>Lx) {xc -= Lx;}
            j = int(xc/4);
            count[j]++;
        }
    }

    double b = tmax/dt/timesteps_per_snapshot;
    for(int k=0; k<25; k++){
        out1<<4*k + 2<<"\t"<<double((count[k]/b)*Nbins)<<endl;
    }
}