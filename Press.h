#ifndef SPS_PRESS
#define SPS_PRESS

#include <particle.h>
#include <matrix.h>

class Press
{
public:
    double z;//крышка
	double r;
    double x_,y_,z_;// координаты коробочки частиы
    double D;
    Press()
    {
        r=3*0.0014;  // радиус 1 см
        D=0.0007+1.2*SPS_PARTICLE_L0;//3 мм
        x_=-r+D/2;
        y_=-r+D/2;
        z_=D/2;
        z=D;
	}
    MathVector3D next(void)
    {
        for(;;)
        {
            x_+=D;
            if(x_>r)
            {
                x_=-r+D/2;
                y_+=D;
            }
            if(y_>r)
            {
                x_=-r+D/2;
                y_=-r+D/2;
                z_+=D;
            }
            if(sqrt(x_*x_+y_*y_)+D/2<r)
            {
                z=z_+D/2;
                return MathVector3D(x_,y_,z_);
            }
//            printf("\nx=%lf\ty=%lf\tz=%lf",x_,y_,z_);
//            getchar();
        }
	}
};

void createIcosahedronParticle(std::vector<Particle> &vParticles, MathVector3D translate,const int typeAtom)
{
    const static double pi=3.14159265358979323846264338327950288419716939937, k=pi/180;
    std::vector<Particle> coorIco;

    //икосаэдр
    double r=SPS_PARTICLE_L0; // радиус сферы
    double x,y,z; //координаты

    // начальные значения для икосаэдра
    double a=4*r/sqrt(10+2*sqrt(5)); // сторона икосаэдра
    double alpha=acos((1-a*a/2/r/r)); // первый угол поворота по тэта
    // cos(alpha)=(1-a*a/2/R/R)

    // вычисляем точки икосаэдра
    //1 точка
    x=y=0.0;
    z=r;

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));
    //2 точка
    x=r*sin(alpha)*sin(0);
    y=r*sin(alpha)*cos(0);
    z=r*cos(alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));

    //3 точка
    x=r*sin(alpha)*sin(72*k);
    y=r*sin(alpha)*cos(72*k);
    z=r*cos(alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));

    //4 точка
    x=r*sin(alpha)*sin(2*72*k);
    y=r*sin(alpha)*cos(2*72*k);
    z=r*cos(alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));

    //5 точка
    x=r*sin(alpha)*sin(3*72*k);
    y=r*sin(alpha)*cos(3*72*k);
    z=r*cos(alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));

    //6 точка
    x=r*sin(alpha)*sin(4*72*k);
    y=r*sin(alpha)*cos(4*72*k);
    z=r*cos(alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));

    //7 точка
    x=r*sin(pi-alpha)*sin(-36*k);
    y=r*sin(pi-alpha)*cos(-36*k);
    z=r*cos(pi-alpha);


    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));


    //8 точка
    x=r*sin(pi-alpha)*sin(36*k);
    y=r*sin(pi-alpha)*cos(36*k);
    z=r*cos(pi-alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));



    //9 точка
    x=r*sin(pi-alpha)*sin((36+72)*k);
    y=r*sin(pi-alpha)*cos((36+72)*k);
    z=r*cos(pi-alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));


    //10 точка
    x=r*sin(pi-alpha)*sin((36+2*72)*k);
    y=r*sin(pi-alpha)*cos((36+2*72)*k);
    z=r*cos(pi-alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));


    //11 точка
    x=r*sin(pi-alpha)*sin((36+3*72)*k);
    y=r*sin(pi-alpha)*cos((36+3*72)*k);
    z=r*cos(pi-alpha);

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));


    //12 точка
    x=0;
    y=0;
    z=-r;

    coorIco.push_back(Particle(typeAtom,MathVector3D(x,y,z)));
//    coorIco.push_back(Particle(MathVector3D(x,y,z),
//                                  MathVector3D(0.01*(double(rand())/RAND_MAX-0.5),
//                                               0.01*(double(rand())/RAND_MAX-0.5),
//                                               0.01*(double(rand())/RAND_MAX-0.5))
//                               ));
    //перенос


    initializedMatrixRotatedX(pi*double(rand())/RAND_MAX);//90*(double(rand())/RAND_MAX)*k);
    initializedMatrixRotatedY(pi*double(rand())/RAND_MAX);//180*(double(rand())/RAND_MAX)*k);
    initializedMatrixTranslate(translate.x,translate.y,translate.z);
    double rezMatrix[4][4];
    productMatrix(MatrixRotatedY,MatrixRotatedX,rezMatrix);
    productMatrix(MatrixTranslete,rezMatrix,rezMatrix);
//    productMatrix(rezMatrix,MatrixRotatedZ,rezMatrix);
//    printf("\n");
//    for(unsigned int i=0;i<4;++i)
//    {
//        for(unsigned int j=0;j<4;++j)
//        {
//            printf("%lf\t",rezMatrix[i][j]);
//        }
//        printf("\n");
//    }
//    for(unsigned int i=0;i<4;++i)
//    {
//        for(unsigned int j=0;j<4;++j)
//        {
//            if(abs(rezMatrix[i][j])<1e-11)
//            {
//                rezMatrix[i][j]=0.0;
//            }
//        }
//    }


    for(unsigned int i=0;i<coorIco.size();++i)
    {
        MathVector3D rezVector;
        //double w=sqrt(coorIco[i].r.x*coorIco[i].r.x+coorIco[i].r.y*coorIco[i].r.y+coorIco[i].r.z*coorIco[i].r.z);
        //cout<<coorIco[i].r.x<<" "<<coorIco[i].r.y<<" "<<coorIco[i].r.z<<endl;
        rezVector.x=rezMatrix[0][0]*coorIco[i].r.x+rezMatrix[0][1]*coorIco[i].r.y+rezMatrix[0][2]*coorIco[i].r.z+rezMatrix[0][3];
        rezVector.y=rezMatrix[1][0]*coorIco[i].r.x+rezMatrix[1][1]*coorIco[i].r.y+rezMatrix[1][2]*coorIco[i].r.z+rezMatrix[1][3];
        rezVector.z=rezMatrix[2][0]*coorIco[i].r.x+rezMatrix[2][1]*coorIco[i].r.y+rezMatrix[2][2]*coorIco[i].r.z+rezMatrix[2][3];
        coorIco[i].r=rezVector;
    }


    for(unsigned int i=0;i<coorIco.size();++i)
    {
        vParticles.push_back(coorIco[i]);
    }

}

void createSphera(std::vector<Particle> &vParticles, MathVector3D translate,const double R,const int typeAtom)
{
    const static double pi=3.14159265358979323846264338327950288419716939937;
    std::vector<Particle> coorSphera;
    double x,y,z; //координаты
    double x0,y0,z0;//центр
    x0=y0=z0=0.0;
    for(int i=-50;i<50;++i)
    {
        for(int j=-50;j<50;++j)
        {
            for(int k=-50;k<50;++k)
            {
                x=x0+SPS_PARTICLE_L0*sqrt(2)*k;
                y=y0+SPS_PARTICLE_L0*sqrt(2)*j;
                z=z0+SPS_PARTICLE_L0*sqrt(2)*i;

                // Гранецентрированная упаковка с равновесныь расстоянием между центрами SPS_PARTICLE_L0
                if(x*x+y*y+z*z<R*R)
                {
                    MathVector3D r1(x,y,z);
                    MathVector3D r2(x+SPS_PARTICLE_L0/sqrt(2), y+SPS_PARTICLE_L0/sqrt(2), z);
                    MathVector3D r3(x, y+SPS_PARTICLE_L0/sqrt(2), z+SPS_PARTICLE_L0/sqrt(2));
                    MathVector3D r4(x+SPS_PARTICLE_L0/sqrt(2), y, z+SPS_PARTICLE_L0/sqrt(2));
                    coorSphera.push_back(Particle(typeAtom,r1));
                    coorSphera.push_back(Particle(typeAtom,r2));
                    coorSphera.push_back(Particle(typeAtom,r3));
                    coorSphera.push_back(Particle(typeAtom,r4));
                }

            }
        }
    }
    initializedMatrixRotatedX(2*pi*double(rand())/RAND_MAX);//90*(double(rand())/RAND_MAX)*k);
    initializedMatrixRotatedY(2*pi*double(rand())/RAND_MAX);//180*(double(rand())/RAND_MAX)*k);
    initializedMatrixTranslate(translate.x,translate.y,translate.z);
    double rezMatrix[4][4];
    productMatrix(MatrixRotatedY,MatrixRotatedX,rezMatrix);
    productMatrix(MatrixTranslete,rezMatrix,rezMatrix);


    for(unsigned int i=0;i<coorSphera.size();++i)
    {
        MathVector3D rezVector;
        //double w=sqrt(coorIco[i].r.x*coorIco[i].r.x+coorIco[i].r.y*coorIco[i].r.y+coorIco[i].r.z*coorIco[i].r.z);
        //cout<<coorIco[i].r.x<<" "<<coorIco[i].r.y<<" "<<coorIco[i].r.z<<endl;
        rezVector.x=rezMatrix[0][0]*coorSphera[i].r.x+rezMatrix[0][1]*coorSphera[i].r.y+rezMatrix[0][2]*coorSphera[i].r.z+rezMatrix[0][3];
        rezVector.y=rezMatrix[1][0]*coorSphera[i].r.x+rezMatrix[1][1]*coorSphera[i].r.y+rezMatrix[1][2]*coorSphera[i].r.z+rezMatrix[1][3];
        rezVector.z=rezMatrix[2][0]*coorSphera[i].r.x+rezMatrix[2][1]*coorSphera[i].r.y+rezMatrix[2][2]*coorSphera[i].r.z+rezMatrix[2][3];
        coorSphera[i].r=rezVector;
    }


    for(unsigned int i=0;i<coorSphera.size();++i)
    {
        vParticles.push_back(coorSphera[i]);
    }
}
void addGranule(std::vector<Particle> &vp, Press &press, const std::string & compound_name, double d, const int typeAtom)
{  // добавить медную круглую гранулу диаметром 0.001 м

    createSphera(vp,press.next(),d/2,typeAtom);
}

void pressForce(const Press & press, std::vector<Particle> &vp)
{
    #pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        const double &x=p.r.x;
        const double &y=p.r.y;
        const double rr=x*x+y*y;
        
        if (rr>(press.r*press.r))
        {
			const double r=sqrt(rr);
        	const double f=3*SPS_PARTICLE_K*(r-press.r);
     		p.a.x+=-f*p.r.x/r;
     		p.a.y+=-f*p.r.y/r;
        }

        if(p.r.z<0.0)
        {
            p.a.z+=-3*SPS_PARTICLE_K*p.r.z;
        }
        if(p.r.z>press.z)
        {
            p.a.z+=-3*SPS_PARTICLE_K*(p.r.z-press.z);
        }
    }

}

void pressForceStretching(const Press &press, std::vector<Particle> &vp, std::vector<unsigned int> &saveIndexUpParticles, std::vector<unsigned int> &saveIndexDownParticles)
{

#pragma omp parallel for
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        Particle &p=vp[i];
        const double &x=p.r.x;
        const double &y=p.r.y;
        const double rr=x*x+y*y;

        if (rr>(press.r*press.r))
        {
            const double r=sqrt(rr);
            const double f=3*SPS_PARTICLE_K*(r-press.r);
            p.a.x+=-f*p.r.x/r;
            p.a.y+=-f*p.r.y/r;
        }

//        if(p.r.z<0.0)
//        {
//            for(unsigned int j=0;j<saveIndexDownParticles.size();++j)
//            {
//                if(i^saveIndexDownParticles[j])
//                {
//                    p.a.z+=-3*SPS_PARTICLE_K*p.r.z;
//                }
//            }
//        }
//        if(p.r.z>press.z)
//        {
//            for(unsigned int j=0;j<saveIndexDownParticles.size();++j)
//            {
//                if(i^saveIndexDownParticles[j])
//                {
//                    p.a.z+=-3*SPS_PARTICLE_K*(p.r.z-press.z);
//                }
//            }
//        }

        if(p.r.z>0.0)
        {
            for(unsigned int j=0;j<saveIndexDownParticles.size();++j)
            {
                if(!(i^saveIndexDownParticles[j]))
                {
                    p.a.z+=-0.1*SPS_PARTICLE_K*p.r.z;
                }
            }

        }
        if(p.r.z<press.z)
        {
            for(unsigned int j=0;j<saveIndexUpParticles.size();++j)
            {
                if(!(i^saveIndexUpParticles[j]))
                {
                    p.a.z+=+0.1*SPS_PARTICLE_K*(p.r.z-press.z);
                }
            }
        }
    }
}

double pressPressureUp(const Press & press, const std::vector<Particle> &vp)
{
	double v=0;
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        const Particle &p=vp[i];
        if (p.r.z>press.z)
        {
            v+=3*SPS_PARTICLE_K*(p.r.z-press.z);
        }
    }

	return v/(M_PI*press.r*press.r);
}

double pressPressureDown(const Press & press, const std::vector<Particle> &vp)
{
    double v=0;
    for(unsigned int i=0; i<vp.size(); ++i)
    {
        const Particle &p=vp[i];
        if(p.r.z<0.0)
        {
            v+=-3*SPS_PARTICLE_K*p.r.z;
        }
    }
    return v/(M_PI*press.r*press.r);
}

double pressPressureWall(const Press & press, const std::vector<Particle> &vp)
{
    double v=0.0;
    for(unsigned int i=0;i<vp.size();++i)
    {
        const Particle &p=vp[i];
        const double &x=p.r.x;
        const double &y=p.r.y;
        const double rr=x*x+y*y;

        if (rr>(press.r*press.r))
        {
            const double r=sqrt(rr);
            v+=3*SPS_PARTICLE_K*(r-press.r);
        }
    }
    return v;
}

#endif
