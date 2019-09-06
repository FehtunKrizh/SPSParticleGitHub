#ifndef ARRAY3dBOX
#define ARRAY3dBOX

#include<map>
#include<vector>
#include<iostream>
#include<cmath>

class Array3DPoint
{
public:
	double x, y, z;
	double l0;
    Array3DPoint(double x, double y, double z, double l0):x(x), y(y), z(z), l0(l0){}
};

bool operator<(const Array3DPoint & a1, const Array3DPoint & a2) // a1<a2
{
	const double & l0=a1.l0;
    if ( int(floor(a1.z/l0)) < int(floor(a2.z/l0)) ) {return true;}
    if ( int(floor(a1.z/l0)) > int(floor(a2.z/l0)) ) {return false;}
    if ( int(floor(a1.y/l0)) < int(floor(a2.y/l0)) ) {return true;}
    if ( int(floor(a1.y/l0)) > int(floor(a2.y/l0)) ) {return false;}
    if ( int(floor(a1.x/l0)) < int(floor(a2.x/l0)) ) {return true;}
    if ( int(floor(a1.x/l0)) > int(floor(a2.x/l0)) ) {return false;}
	return false;
}

bool operator==(const Array3DPoint & a1, const Array3DPoint & a2)// a1==a2
{
	const double & l0=a1.l0;
    if ( (int(a1.z/l0) == int(a2.z/l0)) && (int(a1.y/l0) == int(a2.y/l0)) && (int(a1.x/l0) == int(a2.x/l0)) ){ return true;}
	return false;
}

std::ostream & operator<<(std::ostream & os, const Array3DPoint & p)
{
	return os << p.x << ' ' << p.y << ' ' << p.z << ' ' << p.l0 ;
}

template<typename T> class Array3DBox :public std::map<Array3DPoint, std::vector<T> >
{
public:
	double l0;
    Array3DBox(double l0):l0(l0){}
    void push_element(double x, double y, double z, T p)
    {
        std::map<Array3DPoint, std::vector<T> >::operator[](Array3DPoint(x,y,z, l0)).push_back(p);
	}
};

#endif
