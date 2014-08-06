  
#ifndef COMMON_H_
#define COMMON_H_

#define SQUAREROOT_3 1.732
#define SQUAREROOT_2 1.414
#define PI 3.14159265358979323846

#define GG_PLANARIZATION 0
#define RNG_PLANARIZATION 1

#define myabs(a) ((a) > 0 ? (a) : (-(a)))
#define between(a, asmall, abig) ((a) >= (asmall) && (a) <= (abig))
#define EQUAL_DOUBLE_ERROR double(0.0001)
#define DOUBLE_EQUAL(a, b) (myabs((a)-(b)) < EQUAL_DOUBLE_ERROR)
#define DOUBLE_EQUAL_PT(a, b) (myabs((a.x)-(b.x)) < EQUAL_DOUBLE_ERROR && myabs((a.y)-(b.y)) < EQUAL_DOUBLE_ERROR)

#define CEILING(a) ((int)(a) == (a) ? (a) : ((a)+1))    //only for positive value!

#define mymax(a,b)    (((a) > (b)) ? (a) : (b))
#define mymin(a,b)    (((a) < (b)) ? (a) : (b))

#define INVALID_POINT (Point(IMPOSSIBLE_VALUE, IMPOSSIBLE_VALUE, IMPOSSIBLE_VALUE))
#define BIG_DISTANCE (double)0xffffffff
#define IMPOSSIBLE_VALUE double(0xffffff)
#define BIG_VALUE (double)0xffffffff

extern double COMM_RANGE;

extern int QUASI_UBG;

extern double QUASI_UBG_ALPHA;

#endif


