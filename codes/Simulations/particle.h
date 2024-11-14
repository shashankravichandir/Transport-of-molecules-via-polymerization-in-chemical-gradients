#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double pi = 3.1415926;

class particle{
    public:
    double x,y, phi;
    int alpha;
};

class polymer{
    public:
    unsigned int N;
    vector<particle> chain;

    polymer(int num){
        N = num;
        chain.resize(num);
        for(int i=0; i<N; i++){
            chain[i].alpha = 0;
        }
    }  
};

double CentreOfMass(polymer &p){
    double sum=0.0;
    for(int i=0; i<p.N; i++){
        sum+=p.chain[i].x;
    }
    return sum/p.N;
}

double length(particle p1, particle p2){
    double l;
    l = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return l;
}

double angle(particle p1, particle p2){
    double dx, dy, theta;
    dy = p2.y - p1.y;
    dx = p2.x - p1.x;
    theta = atan(dy/dx);

    if(dy<0 && dx<0){
        theta = theta + pi;
    }

    if(dy>0 && dx<0){
        theta = theta + pi;
    }

    return theta;
}

double fs(double x, double Lx){
    return 20*(1+ sin(2*pi*x/Lx));
}

void AssignActivity(polymer &p, vector<int> &v){
    for(int i=0; i<v.size(); i++){
        p.chain[v[i] - 1].alpha = 1;
    }
}

pair<double, double> SpringForces(polymer &p, int &i, double &ks, double &l0){
    double Fx, Fy, l;

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
    
        return {Fx, Fy};
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

        return {Fx, Fy};
        
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

        return {Fx, Fy};
    }
}