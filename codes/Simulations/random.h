#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

//const double pi = 3.1415926;

void seed_drand48(){
    time_t timer;
    time(&timer);
    srand48((long)timer);
}

inline double uniform_rand(){return drand48();}

double Gaussian(){
    return sqrt(-2*log(1-uniform_rand()))*cos(2*pi*uniform_rand());
}