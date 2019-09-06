#ifndef BOX_H
#define BOX_H
#include <vector>


class Particle;
class Box
{
    public:
        unsigned int i,j,k;
        std::vector<Particle *>  particle;
        Box(void){}
        Box(Particle *inParticle, unsigned int inI, unsigned int inJ, unsigned int inK):i(inI),j(inJ),k(inK)
        {
            particle.push_back(inParticle);
        }

        unsigned int N(unsigned int n)
        {
            return n*n*i+n*j+k;
        }

};

#endif // BOX_H
