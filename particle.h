#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <ostream>


#include <Array3DBox.h>
#include <MathVector3d.h>

const double SPS_E=5e10;//110*100e9;   //0.1*100e9;   // Модуль упругости в паскалях

const double SPS_PARTICLE_L0=100e-6;  // 100 мкм равновесное расстояние между узлами
const double SPS_PARTICLE_K=SPS_PARTICLE_L0*SPS_E;
double SPS_PARTICLE_ALPHA=SPS_PARTICLE_K*SPS_PARTICLE_L0*0.03;

unsigned long long deadBondCount;
unsigned long long birthBondCount;

/*
 Трескается раньше времени
const double SPS_E=0.6*100e9;   // Модуль упругости в паскалях

const double SPS_PARTICLE_L0=0.8*100e-6;  // 100 мкм равновесное расстояние между узлами
const double SPS_PARTICLE_K=SPS_PARTICLE_L0*SPS_E;
double SPS_PARTICLE_ALPHA=SPS_PARTICLE_K*SPS_PARTICLE_L0*0.1/1;
*/

std::map<std::string, int>  compoundIndex;

void SPS_init()
{
        compoundIndex[""]=0;
        compoundIndex["cooper"]=1;
        compoundIndex["aluminium"]=2;
}

template <typename T> class eVector:public std::vector<T>
{
public:
    void erase(int i)
    {
        std::vector<T>::operator[] (i)=std::vector<T>::back();
        std::vector<T>::resize(std::vector<T>::size()-1);
    }
};


//class Bond;
class Particle
{
public:
    int index;
    int typeAtom;
    MathVector3D r,v,a;//r-radius-vector, v-velocity, a-acceleration
    eVector<Particle *> lp;//список указателей на ближайщие частицы(соседи).
    eVector<Particle *> lb;//список связей.
public:

    Particle(const int &inTypeAtom=0,
             const MathVector3D &inR=MathVector3D(0,0,0),
             const MathVector3D &inV=MathVector3D(0,0,0),
             const MathVector3D &inA=MathVector3D(0,0,0),
             const std::string & ch_name=std::string()) :index(compoundIndex[ch_name]), typeAtom(inTypeAtom),r(inR),v(inV),a(inA)  {}
};


class Array3DBox_ :public Array3DBox<Particle*>
{
public:
    Array3DBox_():Array3DBox<Particle*>(2*SPS_PARTICLE_L0){} // коробочки размером 2 lo
};

// строим списки соседей
Array3DBox_ birthNeighbors(std::vector<Particle> &vParticles)
{
    for(unsigned int i=0;i<vParticles.size();++i) // почистили списки соседей у каждой частицы
    {
            vParticles[i].lp.clear();
    }

    Array3DBox_ a;
    for(unsigned int i=0;i<vParticles.size();++i) // разместили указатели на частицы по коробочкам
    {
            const MathVector3D &r=vParticles[i].r;
            a.push_element(r.x, r.y, r.z, &vParticles[i]);
    }

    for(auto box=a.begin(); box!=a.end(); ++box)  // цикл по коробочкам в Array3dbox массиве
    {
        const Array3DPoint &p = box->first;  // координаты одной из частиц в коробочке
        for (int i=-1; i<=1; i++) // цикл по соседним коробочкам
        {
            for (int j=-1; j<=1; j++)
            {
                for (int k=-1; k<=1; k++)
                {
                    double x=p.x+p.l0*k;
                    double y=p.y+p.l0*j;
                    double z=p.z+p.l0*i;

                    auto box2=a.find(Array3DPoint(x,y,z,p.l0));
                    if (box2 == a.end() ){continue;}

                    std::vector<Particle*> &vpp1=box->second;  // vector указателей на частицы содержащийся в центральной коробочке
                    std::vector<Particle*> &vpp2=box2->second;  // vector указателей на частицы содержащийся в боковой коробочке

                    for(auto p1p=vpp1.begin(); p1p!=vpp1.end(); ++p1p)
                    {
                        Particle &p1=**p1p;
                        for(auto p2p=vpp2.begin(); p2p!=vpp2.end(); ++p2p)
                        {
                            if ((*p1p)==(*p2p)){continue;}
                            const Particle &p2=**p2p;
                            if(module2(p1.r-p2.r)<4*SPS_PARTICLE_L0*SPS_PARTICLE_L0)
                            {
                                p1.lp.push_back(*p2p);
                            }
                        }
                    }
                }
            }
        }
    }
    return a;
}

// функция вызывается после расстоновки частиц
void setBonds(std::vector<Particle> & vParticles)
{
    for(unsigned int i=0;i<vParticles.size();++i)
    {
        Particle & p1=vParticles[i];
        for(size_t j=0; j<p1.lp.size(); j++)
        {
            Particle & p2=*(p1.lp[j]);
            MathVector3D dr=p1.r-p2.r;
            if(module2(dr)<1.4*SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                p1.lb.push_back(&p2);
            }
        }
    }
    deadBondCount=birthBondCount=0.0;
}

// вклад сил через связи
void bondForce(std::vector<Particle> &vParticles)
{
    #pragma omp parallel for
    for(unsigned int i=0; i<vParticles.size(); ++i)
    {
        Particle &p1=vParticles[i];
        for(unsigned int j=0; j<p1.lb.size(); ++j)
        {
            const Particle &p2=*(p1.lb[j]);
            MathVector3D dr=p1.r-p2.r;
            double moduleR2=module2(dr);
            // связи работают только на притяжение
            if(moduleR2>SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                double moduleR=module(dr);
                double f=SPS_PARTICLE_K*(moduleR-SPS_PARTICLE_L0);
                p1.a+=(-f/moduleR)*dr;
            }
        }
    }
}

double cBondForce(std::vector<Particle> &vParticles)
{
    double returnF=0.0;
    #pragma omp parallel for
    for(unsigned int i=0; i<vParticles.size(); ++i)
    {
        Particle &p1=vParticles[i];
        for(unsigned int j=0; j<p1.lb.size(); ++j)
        {
            const Particle &p2=*(p1.lb[j]);
            MathVector3D dr=p1.r-p2.r;
            double moduleR2=module2(dr);
            // связи работают только на притяжение
            if(moduleR2>SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                double moduleR=module(dr);
                double f=SPS_PARTICLE_K*(moduleR-SPS_PARTICLE_L0);
                returnF+=f;
            }
        }
    }
    return returnF;
}

int NBond(std::vector<Particle> &vParticles)
{
    int nb=0;
    for(size_t i=0;i<vParticles.size();++i)
    {
        nb+=vParticles[i].lb.size();
    }
    return 0.5*nb;
}

void deadBond(std::vector<Particle> &vParticles)
{
    for(size_t i=0;i<vParticles.size();++i)
    {
        Particle & p1=vParticles[i];
        for(size_t j=0; j<p1.lb.size(); )
        {
            Particle & p2=*(p1.lb[j]);
            MathVector3D dr=p1.r-p2.r;
            if(module2(dr)>1.1*SPS_PARTICLE_L0*SPS_PARTICLE_L0)//if(module2(dr)>1.15*SPS_PARTICLE_L0*SPS_PARTICLE_L0)//1.2
            {
                deadBondCount++;
                p1.lb.erase(j);
                continue;
            }
            j++;
        }
    }
}

void birthBond(std::vector<Particle> &vParticles)
{
    for(unsigned int i=0;i<vParticles.size();++i)
    {
        Particle & p1=vParticles[i];
        for(unsigned int j=0; j<p1.lp.size(); j++)
        {
            Particle & p2=*(p1.lp[j]);
            MathVector3D dr=p1.r-p2.r;
            if(module2(dr)<0.95*SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                eVector<Particle*>::iterator p=find(p1.lb.begin(),p1.lb.end(), &p2);
                if (p==p1.lb.end())
                {
                    birthBondCount++;
                    p1.lb.push_back(&p2);
                }
            }
        }
    }

}


void repulsiveForce(std::vector<Particle> &vParticles)
{
    #pragma omp parallel for
    for(unsigned int i=0; i<vParticles.size(); ++i)
    {
        Particle &p1=vParticles[i];
        for(unsigned int j=0; j<p1.lp.size(); ++j)
        {
            Particle &p2=*(p1.lp[j]);
            MathVector3D dr=p1.r-p2.r;
            double moduleR2=module2(dr);
            // частицы только расталкиваются
            if(moduleR2<SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                double moduleR=sqrt(moduleR2);
                double f=SPS_PARTICLE_K*(moduleR-SPS_PARTICLE_L0);
                p1.a+=(-f/moduleR)*dr;
            }
        }
    }
}

void autoCorelation(std::vector<Particle> &vParticles,std::vector<double> &cor)
{
    cor.clear();
    const unsigned int sizeL=5*5;
    double C[2*sizeL];
    double dl=SPS_PARTICLE_L0/sizeL;
    for(unsigned int i=0;i<2*sizeL;++i)
    {
        C[i]=0.0;
    }
    for(unsigned int i=0;i<vParticles.size();++i)
    {
        Particle &p1=vParticles[i];
        for(unsigned int j=0;j<p1.lp.size();++j)
        {
            Particle &p2=*(p1.lp[j]);
            MathVector3D dr=p1.r-p2.r;
            unsigned int index=sqrt(module2(dr))/dl;
            if(index>2*sizeL-1)
            {
            }
            else
            {
                C[index]++;
            }
        }
    }

    for(unsigned int i=0;i<2*sizeL;++i)
    {
        cor.push_back(C[i]);
    }
}

double cRepulsiveForce(std::vector<Particle> &vParticles)
{
    double returnF;
    #pragma omp parallel for
    for(unsigned int i=0; i<vParticles.size(); ++i)
    {
        Particle &p1=vParticles[i];
        for(unsigned int j=0; j<p1.lp.size(); ++j)
        {
            Particle &p2=*(p1.lp[j]);
            MathVector3D dr=p1.r-p2.r;
            double moduleR2=module2(dr);
            // частицы только расталкиваются
            if(moduleR2<SPS_PARTICLE_L0*SPS_PARTICLE_L0)
            {
                double moduleR=sqrt(moduleR2);
                double f=SPS_PARTICLE_K*(moduleR-SPS_PARTICLE_L0);
                returnF+=f;
            }
        }
    }
    return returnF;
}

class Bond
{
public:
    int p1, p2;
    Bond(int inP1=-1, int inP2=-1):p1(inP1), p2(inP2) {}

    bool operator ==(Bond &a)
    {
        return (p1==a.p1&&p2==a.p2)||(p1==a.p2&&p2==a.p1);
    }

    bool operator !=(Bond &a)
    {
        return !((p1==a.p1&&p2==a.p2)||(p1==a.p2&&p2==a.p1));
    }
};


//void moveParticle(std::vector<Particle> vp, double dt)
//{
  
//}

// Добавить запись в траекторию частиц для WMD
void LAMMPSTRJ (const std::vector <Particle> & p, int step, const char * fname)
{
    std::ofstream os(fname, std::ios_base::app);
    double Rmax = 1;
    os << "ITEM: TIMESTEP\n";
    os << step << "\n";
    os << "ITEM: NUMBER OF ATOMS\n";
    os << p.size() << "\n";
    os << "ITEM: BOX BOUNDS pp pp pp\n";
    os << -Rmax << " " << Rmax << "\n";
    os << -Rmax << " " << Rmax << "\n";
    os << -Rmax << " " << Rmax << "\n";
    os << "ITEM: ATOMS id type xs ys zs\n";

    for (unsigned int i=0; i<p.size(); i++)
    {
        os << i << '\t' << p[i].typeAtom << '\t' << (p[i].r.x)/Rmax << '\t' << (p[i].r.y)/Rmax << '\t' << (p[i].r.z)/Rmax << "\n";
    }
}

//Запис в файл
//void Tinker (const std::vector <Particle> & p, int step, const char * fname)
//{
//    FILE *out;
//    out=fopen(fname,"w+");//открываем файл для записи.
//    if(!out)
//    {
//        printf("File no open inParticles, completion of the program");
//        return;
//    }
//    else
//    {
//        fprintf(out,"%d\n",p.size());
//        for(unsigned int i=0;i<p.size();++i)
//        {
//            fprintf(out,"%d Cu %lf %lf %lf %d\t",i+1,p[i].r.x,p[i].r.y,p[i].r.z,0);
//            for(unsigned int j=0;j<p[i].lb.size();++j)
//            {
//                fprintf(out,"%d ",p[i].lb[j]-&p[i]+i+1);
//            }
//            fprintf(out,"\n");
//        }
//    }
//    fclose(out);
//}

#endif // PARTICLE_H
