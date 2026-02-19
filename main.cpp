#include <iostream>
#include <algorithm>
#include "glm/glm.hpp"
using namespace std; 


/*
Fireworks algorithm (FWA) to find the minimum of sphere function d=2 ie f(x)=x1ˆ2 + x2ˆ2
x1,x2 in range -10 to 10
*/

const int NO_FIREWORKS = 5;
const int NO_ITERATIONS = 20;
const int M_GAUSSIAN = 2;
const float m = 50;
const float A_hat= 5.0;
const int S_min = 5;
const int S_max = 20;
const double XI = 1e-12;   // ξ
const double ZETA = 1e-12; // ζ

//function to minimize
double sphere(const glm::vec2& x){
    return (x.x)* (x.x)+ (x.y)*(x.y);
}
//get random gaussian number 
double rand_gaussian() {
    double u1 = drand48();
    double u2 = drand48();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
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

//find best firework
int best_firework(Firework* fireworks){
    int best = 0;
    for (int i =0; i<NO_FIREWORKS; i++){
        if (fireworks[i].fx < fireworks[best].fx){
            best=i;
        }
    }
    return best;
}

//we calculate the amplitude and the no of sparks per firework
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

//generate uniform sparks : NO_FIREWORKS
void generate_amplitude_sparks(Firework* fireworks, int& total_sparks, glm::vec2*& sparks){
    total_sparks=0;
    for (int i=0; i< NO_FIREWORKS;i++){
        total_sparks += fireworks[i].no_sparks;
    }
    sparks = new glm::vec2[total_sparks + M_GAUSSIAN*NO_FIREWORKS];
    int idx= 0;

    for(int i=0; i< NO_FIREWORKS ; i++){
        glm::vec2 x_i = fireworks[i].x;
        double A_i = fireworks[i].amplitude;
        //generate sparks per firework
        for(int s=0; s< fireworks[i].no_sparks ; s++){
            //randomly select a dimention
            int dim = rand() % 2;
            double delta = ((double)rand() / RAND_MAX) * 2.0 * A_i - A_i;
            sparks[idx] = x_i;
            if (dim==0){
                sparks[idx].x += delta;
            }else{
                sparks[idx].y += delta;
            }

            sparks[idx].x = std::clamp(sparks[idx].x, -10.0f, 10.0f);
            sparks[idx].y = std::clamp(sparks[idx].y, -10.0f, 10.0f);

            idx++;
        }
    }


}

//generate gasussian spark : m
void generate_gaussian_sparks(Firework* fireworks, int& total_sparks,glm::vec2*& sparks, int m_gaussian){
    int idx = total_sparks;
    //generate sparks for m_gaussian random fireworks
    for(int i=0;i<m_gaussian;i++){
        int fw_idx = rand() % NO_FIREWORKS;
        Firework best_fw = fireworks[best_firework(fireworks)];
        glm::vec2 x_i = fireworks[fw_idx].x;

        double epsilon = rand_gaussian();
        sparks[idx].x = x_i.x +  epsilon*(x_i.x - best_fw.x.x);
        sparks[idx].y = x_i.y +  epsilon*(x_i.y - best_fw.x.y);

        sparks[idx].x = std::clamp(sparks[idx].x ,-10.0f,10.0f);
        sparks[idx].y = std::clamp(sparks[idx].y ,-10.0f,10.0f);
        idx++;
       
    }
    total_sparks= idx;
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
    

    int total_sparks;
    //evaluate sparks
    for(int iter=0; iter< NO_ITERATIONS; iter++){
        cout << "Iteration " << iter << endl;
        compute_explotion_params(fireworks);
        total_sparks =0;
        glm::vec2* sparks = NULL;

        // uniform spark generation
        generate_amplitude_sparks(fireworks, total_sparks,sparks);

        //gaussian spark generation
        generate_gaussian_sparks(fireworks, total_sparks,sparks, M_GAUSSIAN);


        //pick next iteration candiates
        Firework* all_candidates= new Firework[total_sparks + NO_FIREWORKS];
        for (int j=0; j<NO_FIREWORKS; j++){
            all_candidates[j] = fireworks[j];
        }
        for (int j=0; j<total_sparks; j++){
            all_candidates[NO_FIREWORKS+j] = Firework(sparks[j],sphere(sparks[j]));
        }

        int best_idx= 0;
        for (int i = 0; i < NO_FIREWORKS; i++) {
            for (int k=0 ; k < total_sparks + NO_FIREWORKS ; k++){
                if (all_candidates[k].fx < all_candidates[best_idx].fx ){
                    best_idx = k;
                }
            }
            fireworks[i] = all_candidates[best_idx];
            all_candidates[best_idx].fx = 1e12;

        }

        cout << "Best firework: " << fireworks[0].x.x << ", " << fireworks[0].x.y  << " f=" << fireworks[0].fx << "\n\n";

        delete[] sparks;
        delete[] all_candidates;

    }

    

    return 0;
}