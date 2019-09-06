//18-10-2018
#include <omp.h>
#include <time.h>

#include <map>
#include <list>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>


#include <particle.h>
//#include <particle_utils.h>
#include <Press.h>

using namespace std;
const MathVector3D g(0,0,-0.05);


void Initial(double &inRd, double &inFinitP, double &inNparticles, double &inDparticle, double &inV);
bool saveFile(std::vector <Particle> &vParticles, int timeStep);

void placeGranul(std::vector<Particle> &vp, Press &press, double d, unsigned int count);
void markParticles(std::vector<Particle> &vp,std::vector<unsigned int> &saveIndexUpParticles,std::vector<unsigned int> &saveIndexDownParticles);

bool moveGranul(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, const double rStop);
bool createMaterial(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, const double finitP);
bool removeTopBounds(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v);
bool removeSideBounds(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v);
bool testCompression(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v);
bool testStreching(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, std::vector<unsigned int> &saveIndexUpParticles,std::vector<unsigned int> &saveIndexDownParticles);

int main(void)
{
    double v;//=0.0004;//0.0015
    double finitP;//=4e9;
    double rd;//=1.4*0.0014;
    double Nparticle;//=2*40000;
    double Dparticle;//=0.0006;
    vector<unsigned int> saveIndexUpParticles;
    vector<unsigned int> saveIndexDownParticles;
	{
        ofstream ofs1("sps.lammpstrj");
        ofstream ofs2("p.dat");
        ofstream ofs3("corel.dat");
	}

    Press press;// press.r=0.0012;  // 1 см
	vector<Particle> vp;

    Initial(rd,finitP,Nparticle,Dparticle,v);

    placeGranul(vp,press,Dparticle,Nparticle);//расположили гранулы

    LAMMPSTRJ(vp, 0, "sps.lammpstrj");
    saveFile(vp, 0);
	
	birthNeighbors(vp);  // построили списки соседей
    setBonds(vp);       // построили парные связи

	int timeStep=0;
	double timetime=0;
    double dt=0.02*2*M_PI*sqrt(1/SPS_PARTICLE_K);
    printf("\n Particles=%d",vp.size());
    printf("\n Bonds=%d",NBond(vp));
    printf("\n RealTime=%lf",dt);

    while(moveGranul(vp,press,timeStep,timetime,dt,v,rd)){}//0.00125)){}
    while(createMaterial(vp,press,timeStep,timetime,dt,v,finitP)){}
    while(removeTopBounds(vp,press,timeStep,timetime,dt,v)){}
    while(removeSideBounds(vp,press,timeStep,timetime,dt,v)){}//*0.00125)){}
    markParticles(vp,saveIndexUpParticles,saveIndexDownParticles);
    while(testCompression(vp,press,timeStep,timetime,dt,v)){}
    //while(testStreching(vp,press,timeStep,timetime,dt,v,saveIndexUpParticles,saveIndexDownParticles)){}
	

    printf("\nThe end;");
	return 0;
}

bool saveFile(std::vector <Particle> &vParticles, int timeStep)
{
    char fname[512];
    list<Bond> lBonds;
    FILE *outParticles, *outBonds;


    sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Data//outParticles_time=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
    outParticles=fopen(fname,"w+");//открываем файл для записи.
    if(!outParticles)
    {
        printf("File no open inParticles, completion of the program");
        return false;
    }
    else
    {
        fprintf(outParticles,"%d\n",vParticles.size());
        for(unsigned int count=0;count<vParticles.size();++count)
        {
            fprintf(outParticles,"%lf\n",vParticles[count].r.x);
            fprintf(outParticles,"%lf\n",vParticles[count].r.y);
            fprintf(outParticles,"%lf\n",vParticles[count].r.z);
        }
    }
    fclose(outParticles);

    lBonds.clear();
    for(unsigned int i=0;i<vParticles.size();++i)
    {
        for(unsigned int j=0;j<vParticles[i].lb.size();++j)
        {
            Bond JamesBond(i,vParticles[i].lb[j]-&vParticles[i]+i);
            lBonds.push_back(JamesBond);
        }
    }

    sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Data//outBonds_time=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
    outBonds=fopen(fname,"w+");//открываем файл для чтение.
    if(!outBonds)
    {
        //printf("File no open inParticles, completion of the program");
        return false;
    }
    else
    {
        fprintf(outBonds,"%d\n",lBonds.size());

        for(list<Bond>::iterator b=lBonds.begin();b!=lBonds.end();++b)
        {
            fprintf(outBonds,"%d\n",b->p1);
            fprintf(outBonds,"%d\n",b->p2);
        }
    }
    fclose(outBonds);

    return true;
}

void placeGranul(std::vector<Particle> &vp, Press &press, double d, unsigned int count)
{
    for(;vp.size()<count;)
    {
        //printf("\naddGranule");
        double r = rand() / double( RAND_MAX ) ;
        if(r>0.3)
        {
            if(r>0.6)
            {
                addGranule(vp, press, "cooper", d, 0);  // добавить медную круглую гранулу диаметром 0.0006 м
            }
            else
            {
                addGranule(vp, press, "cooper", d, 1);
            }
        }
        else
        {
            addGranule(vp, press, "cooper", d, 2);
        }
    }
    press.z=press.z*100;
}

bool moveGranul(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, const double rStop)
{
    double vr=-v;
    double vz=0.0;


    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForce(press,vp);

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");
        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;


    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    if(press.r<=rStop)
    {
        press.z=vp[0].r.z;
        for(unsigned int i=1;i<vp.size();++i)
        {
            if(vp[i].r.z>press.z)
            {
                press.z=vp[i].r.z;
            }
        }
        press.z=press.z+10*SPS_PARTICLE_L0;

        printf("\n-------------------------");
        printf("\n r=%lf",press.r);
        printf("\n-------------------------");
        return false;
    }
    return true;
}
bool createMaterial(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, const double finitP)
{
    double vr=0.0;
    double vz=-v;

    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForce(press,vp);

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");
        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;

    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    if(pressureUp>=finitP)
    {
        return false;
    }
    return true;
}
bool removeTopBounds(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v)
{
    double vr=0.0;
    double vz=v;

    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForce(press,vp);

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");

        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;


    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    if(pressureUp==0.0)
    {
        return false;
    }
    return true;
}
bool removeSideBounds(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v)
{
    double vr=v;
    double vz=0.0;

    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForce(press,vp);

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");
        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;

    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    if(pressPressureWall(press,vp)==0)
    {
        return false;
    }
//    if(press.r>=rStop)
//    {
//        return false;
//    }
    return true;
}
bool testCompression(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v)
{
    double vr=2*v;
    double vz=-v;
    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForce(press,vp);

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");
        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;

    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    if(press.z<0.001)
    {
        return false;
    }
    return true;
}
void markParticles(std::vector<Particle> &vp,std::vector<unsigned int> &saveIndexUpParticles,std::vector<unsigned int> &saveIndexDownParticles)
{
    unsigned int saveIndex=0;
    double maxZ=vp[saveIndex].r.z;
    for(unsigned int i=saveIndex;i<vp.size();++i)
    {
        //поиск самой верхней частицы.// псевдо очень псевдо уровень ферми XD
        if(vp[i].r.z>maxZ)
        {
            maxZ=vp[i].r.z;
            saveIndex=i;
        }
    }
    for(unsigned int i=0;i<vp.size();++i)
    {
        if(vp[i].r.z<3*SPS_PARTICLE_L0)
        {
            saveIndexDownParticles.push_back(i);
        }
        if(vp[i].r.z>vp[saveIndex].r.z-3*SPS_PARTICLE_L0)
        {
            saveIndexUpParticles.push_back(i);
        }
    }
}

bool testStreching(std::vector<Particle> &vp, Press &press, int &timeStep, double &timetime, const double dt, const double v, std::vector<unsigned int> &saveIndexUpParticles, std::vector<unsigned int> &saveIndexDownParticles)
{
    double vr=2*v;
    double vz=v;
    if(timeStep%100==0)
    {
        birthNeighbors(vp);
    }

    birthBond(vp);
    deadBond(vp);

    bondForce(vp);
    repulsiveForce(vp);
    pressForceStretching(press,vp,saveIndexUpParticles,saveIndexDownParticles);
//    for(unsigned int i=0;i<saveIndexUpParticles.size();++i)
//    {
//        vp[saveIndexUpParticles[i]].v=MathVector3D(0,0,vz);
//    }
//    for(unsigned int i=0;i<saveIndexDownParticles.size();++i)
//    {
//        vp[saveIndexDownParticles[i]].v=MathVector3D(0,0,0);
//    }

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        p.v+=(p.a-SPS_PARTICLE_ALPHA*p.v)*dt;
        p.r+=p.v*dt;
        p.a=vZero;
    }

    if(timeStep%500==0)
    {
        //saveFile(vp, timeStep);
        printf("\n timeStep=%d",timeStep);
        printf("\n Particles=%d",vp.size());
        printf("\n Bonds=%d",NBond(vp));
        LAMMPSTRJ(vp, timeStep, "sps.lammpstrj");
        //Tinker(vp,timeStep,"sps.arc");
    }

    press.z+=vz*dt;
    press.r+=vr*dt;

    timetime+=dt;
    timeStep++;

    double pressureUp=pressPressureUp(press, vp);
    double nb=NBond(vp)/(1.0*vp.size());
    double pressureDown=pressPressureDown(press,vp);
    if(timeStep%50==0)
    {
        ofstream ofs("p.dat", ios_base::app);
        //если записывать общию силу bondForce+ repulsiveForse быть может получиться померить на растежение
        //или только бонд
        ofs<<timeStep << '\t' << timetime << '\t' << pressureUp << '\t'<< pressureDown << '\t' << press.z << '\t' << nb<<'\t'<< deadBondCount <<'\t'<< birthBondCount << '\n';
        ofs.close();
        char fname[512];
        sprintf(fname,"D://QtProjects//SPS//SPS_Particle//Correlation//cor=%d.dat",timeStep);//в чаровую переменную fname помещаем имя файла
        ofstream ofs1(fname);
        vector<double> cor;
        autoCorelation(vp,cor);
        for(unsigned int i=0;i<cor.size();++i)
        {
            ofs1<<i<<'\t'<<cor[i]<<endl;
        }
        ofs1.close();
    }
    return true;

}

void Initial(double &inRd, double &inFinitP, double &inNparticles, double &inDparticle, double &inV)
{
    FILE * IN;
    IN=fopen("D://QtProjects//SPS//SPS_Particle//build-test_2-Desktop_Qt_5_6_2_MinGW_32bit3-Release//Configuration.dat","r");
    if(fscanf(IN,"rd=%lf;\n",&inRd))
    {
        printf("rd = %5.8lf;\n",inRd);
    }
    else
    {
        printf("data from the file have not been read");
    }

    if(fscanf(IN,"finitP=%lf;\n",&inFinitP))
    {
         printf("finitP = %5.8lf;\n",inFinitP);
    }
    else
    {
        printf("data from the file have not been read");
    }

    if(fscanf(IN,"Nparticles=%lf;\n",&inNparticles))
    {
        printf("Nparticles = %lf;\n",inNparticles);
    }
    else
    {
        printf("data from the file have not been read");
    }

    if(fscanf(IN,"Dparticle=%lf;\n",&inDparticle))
    {
        printf("Dparticle = %5.9lf;\n",inDparticle);
    }
    else
    {
        printf("data from the file have not been read");
    }

    if(fscanf(IN,"velocity=%lf;\n",&inV))
    {
        printf("velocity=%5.9lf;\n",inV);
    }
    else
    {
        printf("data from the file have not been read");
    }

    fclose(IN);
}
