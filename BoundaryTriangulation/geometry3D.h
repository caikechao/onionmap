#ifndef GEOMETRY_H_
#define GEOMETRY_H_


#include <vector>
#include <set>
#include "common.h"
#include "math.h"
#include <stdio.h>

#define PS_ONLINE 0
#define PS_TOLEFT 1
#define PS_TORIGHT 2

using namespace std;

//does not include boundary!

extern int pnpoly_C(int npol, long double *xp, long double *yp, long double x, long double y);

class Point
{
public:
	Point() {x = 0; y = 0; z = 0;boundaryType = 0;}
	Point(long double x_,long double y_):x(x_),y(y_),z(0),boundaryType(0){};
	Point(long double x_, long double y_, long double z_) : x(x_), y(y_), z(z_),boundaryType(0){};
	~Point(){};
	
	long double X() {return x;}
	long double Y() {return y;}
	long double Z() {return z;}

	long double x;
	long double y;
	long double z;

	int boundaryType;

	set<int> levelTypes;
	//
// 	void initialize()
// 	{
// 		x = 0;
// 		y = 0;
// 		z = 0;
// 		levelTypes.clear();
// 		boundaryType = 0;
// 	}
	// for vector...
	Point crossProduct(const Point& vecPoint)
	{
		
		long double cpx = y*vecPoint.z - z*vecPoint.y;
		long double cpy = z*vecPoint.x - x*vecPoint.z;
		long double cpz = x*vecPoint.y - y*vecPoint.x;

		return Point(cpx,cpy,cpz);

	}
	// for vector...
	Point operator - (const Point& subtractor_) const
	{
		long double lx = x-subtractor_.x;
		long double ly = y-subtractor_.y;
		long double lz = z-subtractor_.z;

		return Point(lx,ly,lz);
	}

	inline long double distance(Point P) { 
		long double xDist = pow(P.x - x, 2);
		long double yDist = pow(P.y - y, 2);
		long double zDist = pow(P.z - z, 2);
		return sqrt(xDist + yDist + zDist);
	};
	
	inline long double distance(double X, double Y, double Z) { 		
		return distance(Point(X, Y, Z));
	};


	inline Point & operator = (const Point &right)
	{
		x = right.x;
		y = right.y;
		z = right.z;

		return *this;
	}

	inline bool operator == (const Point &right) const
	{
		return (x == right.x && y == right.y && z == right.z);
	}

	inline bool operator != (const Point &right) const
	{
		return (x != right.x || y != right.y || z != right.z);
	}

	inline bool approxInequal(const Point &right)
	{
		return (!DOUBLE_EQUAL(x, right.x) || !DOUBLE_EQUAL(y, right.y) || !DOUBLE_EQUAL(z, right.z));
	}

	inline bool approxEqual(const Point &right)
	{
		return (DOUBLE_EQUAL(x, right.x) && DOUBLE_EQUAL(y, right.y) && DOUBLE_EQUAL(z, right.z));
	}

	inline bool isValid()
	{
		return (x != IMPOSSIBLE_VALUE && y != IMPOSSIBLE_VALUE && z != IMPOSSIBLE_VALUE);
	}

	static Point getMiddlePoint(const Point & p1, const Point & p2) 
	{
		return Point((p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0, (p1.z + p2.z)/2.0);
	}
};


#endif


