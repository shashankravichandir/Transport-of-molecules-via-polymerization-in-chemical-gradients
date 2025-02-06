// Program to generate the steady state profiles for active semiflexible chains with each monomer either being active or passive

#include <iostream>
#include <fstream>
#include <string>
#include "particle.h"
#include "random.h"

using namespace std;

int main(int argc, char** argv){
    seed_drand48();

    // Parameters
    double kb = 32.0, ks = 8.0, Dt = 1.0, Dr = 5.0, l0 = 0.0;
    double Lx = 100.0; double Ly = 100.0;

    //Intializing polymer chain
    int N = 6; vector<int> alpha{1,0,0,0,0,0};
    polymer p(N);
    for(int i=0; i<N; i++){
        p.chain[i].x = 50.0 - i*5.0;
        p.chain[i].y = 50.0 - i*5.0;
        p.chain[i].phi = 2*pi*uniform_rand();
    }
    polymer p_new = p;

    // Other variables
    double t = 0.0; double tmax = 5000000.0; double dt = 0.001; 
    double eta_x, eta_y, eta_phi;
    double l, Fx, Fy;

    // Averaging variables
    int timestep=0; double xc; const int Nbins = 25; int count[Nbins]; int j; int timesteps_per_snapshot = 1000;
    for(int k = 0; k<25; k++){
        count[k] = 0;
    }
    
    // Opening files
    ofstream out1, out2;
    out1.open("./Results/Semiflexible/l0_0/kb_" + to_string(int(kb)) +"/ks_" + to_string(int(ks)) +".dat");
    out2.open("Progress_kb" + to_string(int(kb)) + "_ks" + to_string(int(ks)) + ".dat");

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

                Fx += 2*kb*(p.chain[i-1].x - 2*p.chain[i].x + p.chain[i+1].x);
                Fy += 2*kb*(p.chain[i-1].y - 2*p.chain[i].y + p.chain[i+1].y);

                if(i!=1 && i!=N-2){
                    Fx += -kb*(p.chain[i].x - 2*p.chain[i+1].x + p.chain[i+2].x);
                    Fy += -kb*(p.chain[i].y - 2*p.chain[i+1].y + p.chain[i+2].y); 

                    Fx += -kb*(p.chain[i-2].x - 2*p.chain[i-1].x + p.chain[i].x);
                    Fy += -kb*(p.chain[i-2].y - 2*p.chain[i-1].y + p.chain[i].y);
                }
                else if(i==1){
                    Fx += -kb*(p.chain[i].x - 2*p.chain[i+1].x + p.chain[i+2].x);
                    Fy += -kb*(p.chain[i].y - 2*p.chain[i+1].y + p.chain[i+2].y); 
                }
                else if(i==N-2){
                    Fx += -kb*(p.chain[i-2].x - 2*p.chain[i-1].x + p.chain[i].x);
                    Fy += -kb*(p.chain[i-2].y - 2*p.chain[i-1].y + p.chain[i].y);
                }
            }
            else if(i==0){
                l = length(p.chain[i], p.chain[i+1]);
                Fx = ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
                Fy = ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;

                Fx += -kb*(p.chain[i].x - 2*p.chain[i+1].x + p.chain[i+2].x);
                Fy += -kb*(p.chain[i].y - 2*p.chain[i+1].y + p.chain[i+2].y); 
            }
            else if(i==N-1){
                l = length(p.chain[i-1], p.chain[i]);
                Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
                Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;

                Fx += -kb*(p.chain[i-2].x - 2*p.chain[i-1].x + p.chain[i].x);
                Fy += -kb*(p.chain[i-2].y - 2*p.chain[i-1].y + p.chain[i].y);
            }
            
            eta_x = Gaussian();
            eta_y = Gaussian();
            eta_phi = Gaussian();

            p_new.chain[i].x += dt*(Fx + alpha[i]*fs(p.chain[i].x, Lx)*cos(p.chain[i].phi)) + sqrt(2*Dt*dt)*eta_x;
            p_new.chain[i].y += dt*(Fy + alpha[i]*fs(p.chain[i].x, Lx)*sin(p.chain[i].phi)) + sqrt(2*Dt*dt)*eta_y;
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

        if(timestep%100000 == 0){
            out2<<timestep<<"\t"<<t<<endl;
        }
    }

    double b = tmax/dt/timesteps_per_snapshot;
    for(int k=0; k<25; k++){
        out1<<4*k + 2<<"\t"<<double(count[k]/b*Nbins)<<endl;
    }
}