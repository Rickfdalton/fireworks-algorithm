#include <iostream>
#include <algorithm>
#include "glm/glm.hpp"
using namespace std; 


/*
Fireworks algorithm (FWA) to find the minimum of sphere function d=2 ie f(x)=x1ˆ2 + x2ˆ2
x1,x2 in range -10 to 10
*/

const int NO_FIREWORKS = 5;
const int NO_ITERATIONS = 10;
const float m = 50;
const float A_hat= 5.0;
const int S_min = 5;
const int S_max = 20;
const double XI = 1e-12;   // ξ
const double ZETA = 1e-12; // ζ

double sphere(const glm::vec2& x){
    return (x.x)* (x.x)+ (x.y)*(x.y);
}

class Firework{
public:
    glm::vec2 x;
    double fx;
    double amplitude;
    int no_sparks;

    Firework(){};
    ~Firework(){};
    Firework(glm::vec2& x, double fx):x(x),fx(fx){}

};

void compute_explotion_params(Firework* fireworks){

    //find the max and min
    double y_max = fireworks[0].fx;
    double y_min = fireworks[0].fx;

    for (int i=0; i<NO_FIREWORKS ; i++){
        y_max = std::max(fireworks[i].fx,y_max);
        y_min = std::min(fireworks[i].fx,y_min);
    }

    //get the denominator ∑(𝑦𝑚𝑎𝑥− 𝑓(𝑥𝑖))+ 𝜁 for no_spark calc
    double denominator_spark =0;

    for (int i=0; i<NO_FIREWORKS ; i++){
        denominator_spark += y_max-fireworks[i].fx + ZETA;
    }

    //get the denominator ∑( 𝑓(𝑥𝑖))− 𝑦𝑚𝑖𝑛+ 𝜁 for amplitude calc
    double denominator_ampl =0;

    for (int i=0; i<NO_FIREWORKS ; i++){
        denominator_ampl += fireworks[i].fx-y_min+ ZETA;
    }

    //assign the no of sparks and the amplitude of each firework
    for (int i=0; i<NO_FIREWORKS ; i++){
        fireworks[i].no_sparks = int(m * (y_max-fireworks[i].fx + XI) / denominator_spark);
        fireworks[i].amplitude = A_hat * (fireworks[i].fx - y_min + XI)/ denominator_ampl;
        fireworks[i].no_sparks = std::clamp(fireworks[i].no_sparks,S_min, S_max);

    }

}

int main(){

    glm::vec2 firework_positions[NO_FIREWORKS];
    double y_values[NO_FIREWORKS];
    Firework fireworks[NO_FIREWORKS];

    //generate random initial fireworks
    for(int i=0;i<NO_FIREWORKS;i++){
        firework_positions[i]= glm::vec2((2*drand48()-1)*10, (2*drand48()-1)*10);
        y_values[i] = sphere(firework_positions[i]);
        fireworks[i] = Firework(firework_positions[i], y_values[i]);
    }
    
    compute_explotion_params(fireworks);

    return 0;
}