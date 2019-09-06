#ifndef MATHVECTOR3D_H
#define MATHVECTOR3D_H
#include <cmath>
#include <iostream>


class MathVector3D
{
    public:
        double x;
        double y;
        double z;
    public:
        MathVector3D(const MathVector3D &obj)
        {
            x=obj.x;
            y=obj.y;
            z=obj.z;
        }

        explicit MathVector3D(const double x=0, const double y=0, const double z=0):x(x),y(y),z(z){}


        MathVector3D & operator +=(const MathVector3D & value)
        {
            x+=value.x;
            y+=value.y;
            z+=value.z;
            return *this;
        }
        inline MathVector3D& operator -=(const MathVector3D &value)
        {
            x-=value.x;
            y-=value.y;
            z-=value.z;
            return *this;
        }
        inline const MathVector3D & operator =(const MathVector3D &value)
        {
            x=value.x;
            y=value.y;
            z=value.z;
            return *this;
        }

        inline const MathVector3D &operator *=(const double &value)
        {
            x*=value;
            y*=value;
            z*=value;
            return *this;
        }

};

const MathVector3D vZero(0,0,0);

std::ostream& operator<<(std::ostream& os, const MathVector3D & v)
{
	return os << v.x << ' ' << v.y << ' ' << v.z;
}

inline MathVector3D operator -(const MathVector3D &left, const MathVector3D &right)
{
    return MathVector3D(left.x-right.x,left.y-right.y,left.z-right.z);
}

inline MathVector3D operator +(const MathVector3D &left, const MathVector3D &right)
{
    return MathVector3D(left.x+right.x,left.y+right.y,left.z+right.z);
}

inline MathVector3D operator *(const MathVector3D &left, const double &right)
{
    return MathVector3D(right*left.x,right*left.y,right*left.z);
}

inline MathVector3D operator *(const double &left, const MathVector3D &right)
{
    return MathVector3D(left*right.x,left*right.y,left*right.z);
}

inline double module2(const MathVector3D & v)
{
    return v.x*v.x+v.y*v.y+v.z*v.z;
}

       
inline double module(const MathVector3D & v)
{
    return sqrt(module2(v));
}

double scalarProduct(const MathVector3D &a, const MathVector3D &b) 
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

MathVector3D vectorProduct(const MathVector3D &a, const MathVector3D &b)
{
    double i,j,k;
    i=a.y*b.z-a.z*b.y;
    j=-(a.x*b.z-a.z*b.x);
    k=a.x*b.y-a.y*b.x;
    return MathVector3D(i,j,k);
}
        

inline void set_zero(MathVector3D & v)
{
	v.x=v.y=v.z=0.0;
}

#endif // MATHVECTOR3D_H
