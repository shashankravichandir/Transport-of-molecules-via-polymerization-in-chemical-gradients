// Program to calculate the mean first passage time for active flexible chains with each monomer either being active or passive

#include <iostream>
#include <fstream>
#include <string>
#include "particle.h"
#include "random.h"

using namespace std;

int main(int argc, char** argv){
    seed_drand48();

    // Parameters
    double ks = 8.0, Dt = 1.0, Dr = 5.0, l0 = 0.0;
    double Lx = 100.0; double Ly = 100.0;

    //Intializing polymer chain
    int N = 2; 
    polymer p(N), p_new(N);
    vector<int> v{1,2};
    AssignActivity(p,v);

    // Force variables;
    double Fx, Fy, l;

    // Time variables
    double t = 0.0; double dt = 0.01; int num_simulations = 100000;
    
    // Opening files
    ofstream out1;
    out1.open("./Results/MFPT/N" + to_string(N) +"/alpha_all.dat");

    // Time loop
    
    for(int k=1; k<=num_simulations; k++){
        // Initializing the polymer position
        for(int i=0; i<N; i++){
            p.chain[i].x = 50.0;
            p.chain[i].y = 50.0;
            p.chain[i].phi = 2*pi*uniform_rand();
        }
        p_new = p; 
        
        // Time loop
        t = 0;
        while(1){
            t = t+dt;

            for(int i = 0; i<N; i++){
                Fx = 0.0; Fy = 0.0;
                
                if(i!=0 && i!=p.N-1){
                    l = length(p.chain[i], p.chain[i-1]);
                    if(l>1e-6){
                        Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
                        Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
                    }
                    else{
                        Fx = 0.0;
                        Fy = 0.0;
                    }

                    l = length(p.chain[i], p.chain[i+1]);
                    if(l>1e-6){
                        Fx += ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
                        Fy += ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
                    }
                    else{
                        Fx += 0.0;
                        Fy += 0.0;
                    }
                }
                else if(i==0){
                    l = length(p.chain[i], p.chain[i+1]);
                    if(l>1e-6){
                        Fx = ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
                        Fy = ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
                    }
                    else{
                        Fx = 0.0;
                        Fy = 0.0;
                    }
                }
                else{
                    l = length(p.chain[i-1], p.chain[i]);
                    if(l>1e-6){
                        Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
                        Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
                    }
                    else{
                        Fx = 0.0;
                        Fy = 0.0;
                    }       
                }

                p_new.chain[i].x += dt*(Fx + p.chain[i].alpha*fs(p.chain[i].x, Lx)*cos(p.chain[i].phi)) + sqrt(2*Dt*dt)*Gaussian();
                p_new.chain[i].y += dt*(Fy + p.chain[i].alpha*fs(p.chain[i].x, Lx)*sin(p.chain[i].phi)) + sqrt(2*Dt*dt)*Gaussian();
                p_new.chain[i].phi += sqrt(2*Dr*dt)*Gaussian();    
            }

            p.chain = p_new.chain;

            double xc = CentreOfMass(p);
            while(xc>Lx){xc -= Lx;}
            while(xc<0){xc += Lx;}
            if(xc<Lx/4 || xc>5*Lx/4){
                out1<<k<<'\t'<<t<<endl;
                break;
            }
        }
    }
}