#ifndef _RANDOM_H
#define _RANDOM_H

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#ifndef MAXINT
#define	MAXINT	2147483647	// XX [for now]
#endif

/*
 * random number generator
 * Again these codes are stolen from ns-2
 */
class RNG {

public:
	RNG(){}
	RNG(int seed) { set_seed( seed); };
	void set_seed(int seed = 1) {srand(seed);}
	inline static RNG* defaultrng() { return (default_); }

	inline double uniform() {return (double)rand()/(double)RAND_MAX;}// range [0.0, 1.0)

    
    inline unsigned long uniform32() {return (unsigned long)(uniform() * (unsigned long)0xffffffff); }
    
	inline double uniform(double r) { return (r * uniform());}
	inline double uniform(double a, double b) { return (a + uniform(b - a)); }
	inline double exponential()	{ 
		double u = uniform();
		while (u == 0.0) {
			u = uniform();
		}
		return (-log(u)); 
	}
	inline double exponential(double r)	{ return (r * exponential());}

	inline double pareto(double scale, double shape) { 
		// When 1 < shape < 2, its mean is scale**shape, and its 
		// variance is infinite.
		double u = uniform();
		while (u == 0.0 || u == 1.0) {
			u = uniform();
		}
		return (scale * (1.0/pow(u, 1.0/shape)));
	}
	
	//a- shape, k - scale, p - bound
	inline double paretoBd(double a, double k, double p) {
		double z;
		double rv;

		do {
			z = uniform();
		} while(z == 0 || z == 1);

		rv = -(z * pow(p, a) - z * pow(k, a) - pow(p, a))/(pow(p, a) * pow(k, a));
		rv = pow(rv, (-1/a));

		return rv;
	}

    inline double paretoII(double scale, double shape) { 
		double z;
		do {
			z = uniform();
		} while(z == 0 );

		return (scale * ((1.0/pow(uniform(), 1.0/shape)) - 1));
	}
	double normal(double avg, double std) {
		static int parity = 0;
		static double nextresult;
		double sam1, sam2, rad;
		
		if (std == 0) return avg;
		if (parity == 0) {
			sam1 = 2*uniform() - 1;
			sam2 = 2*uniform() - 1;
			while ((rad = sam1*sam1 + sam2*sam2) >= 1) {
				sam1 = 2*uniform() - 1;
				sam2 = 2*uniform() - 1;
			}
			rad = sqrt((-2*log(rad))/rad);
			nextresult = sam2 * rad;
			parity = 1;
			return (sam1 * rad * std + avg);
		}
		else {
			parity = 0;
			return (nextresult * std + avg);
		}
	}
	inline double lognormal(double avg, double std) { return (exp (normal(avg, std))); }

	static RNG* default_;
}; 

class RandomVariable {
 public:
	virtual double value() = 0;
	virtual double avg() = 0;
	RandomVariable(){ rng_ = RNG::defaultrng();};

 protected:
	RNG* rng_;
};

class UniformRandVar : public RandomVariable {
 public:
	UniformRandVar(double min = 0, double max = 1){
		min_ = min;
		max_ = max;
	}

    virtual double value(){ return(rng_->uniform(min_, max_)); }
	virtual inline double avg() { return (max_-min_)/2; };

	double* minp()	{ return &min_; };
	double* maxp()	{ return &max_; };
	double min()	{ return min_; };
	double max()	{ return max_; };
	void setmin(double d)	{ min_ = d; };
	void setmax(double d)	{ max_ = d; };
 private:
	double min_;
	double max_;
};


class ExponentialRandVar : public RandomVariable {
 public:

	ExponentialRandVar(double avg = 1) {avg_ = avg; }

	virtual double value(){ return(rng_->exponential(avg_)); }
	double* avgp() { return &avg_; };
	virtual inline double avg() { return avg_; };
	void setavg(double d) { avg_ = d; };

 private:
	double avg_;
};

/*
 *	It has infinite mean if a <= 1 and mean ab/(a-1) if a > 1.

	The variance is ab**2/(a-1)**2(a-2) if a > 2, 
		infinite if 1 < a <= 2, and undefined if a <= 1. 

 */
class ParetoRandVar : public RandomVariable {
 public:

	ParetoRandVar(double avg = 1, double shape = 0){
		avg_ = avg;
		shape_ = shape;
		scale_ = avg_ * (shape_ -1)/shape_;
	};	

	ParetoRandVar(double scale, double shape, int dumb) {
		shape_ = shape;
		scale_= scale;
		avg_ = scale_ * shape_/(shape_ - 1);
	}

	virtual double value() { return(rng_->pareto(scale_, shape_));}
	double* avgp() { return &avg_; };
	double* shapep() { return &shape_; };
	virtual double avg()	{ return avg_; };
	double shape()	{ return shape_; };
	void setavg(double d)	{ avg_ = d; };
	void setshape(double d)	{ shape_ = d; };

 protected:
	double avg_;
	double shape_;
	double scale_;
};

class ShiftedParetoRandVar : public RandomVariable {
 public:

	ShiftedParetoRandVar(double scale, double shape) {
		shape_ = shape;
		scale_= scale;
		avg_ = scale_ * shape_/(shape_ - 1);
	}

	virtual double value() { return(rng_->pareto(scale_, shape_) - scale_);}
	virtual double avg()	{ return avg_ - scale_; };
	double shape()	{ return shape_; };
	void setavg(double d)	{ avg_ = d; };
	void setshape(double d)	{ shape_ = d; };

 protected:
	double avg_;
	double shape_;
	double scale_;
};

/*
 *	Bounded Pareto Distribution 
 */
class BdParetoRandVar : public ParetoRandVar {
 public:

	BdParetoRandVar(double shape, double scale, double bound) {
		shape_ = shape;
		scale_= scale;
		bound_ = bound;
		avg_ = pow(scale_, shape_)/(1-pow(scale_/bound_, shape_));
		avg_ *= (shape_/(shape_-1));
		avg_ *= (1/pow(scale_, shape_-1) - 1/pow(bound_, shape_-1));
	}

	virtual double value() { 
		double val = rng_->paretoBd(shape_, scale_,  bound_);
		return val;
	}

	double bound() { return bound_;	}
	double scale() { return scale_; }
	double shape() { return shape_;	}

	void setscale(double scale) { scale_ = scale; }
	
private:
	double bound_;
};

class NormalRandVar : public RandomVariable {
 public:
		NormalRandVar(double avg = 1, double std = 0) {avg_ = avg; std_ = std;}
	    virtual double value(){ return(rng_->normal(avg_, std_)); }        
        inline double* avgp() { return &avg_; };
        inline double* stdp() { return &std_; };
        virtual double avg()     { return avg_; };
        inline double std()     { return std_; };
        inline void setavg(double d)    { avg_ = d; };
        inline void setstd(double d)    { std_ = d; };
 private:
        double avg_;
        double std_;
};

class LogNormalRandVar : public RandomVariable {
public:	
        LogNormalRandVar(double logavg = 1, double logstd = 0) 
		{logavg_ = logavg; logstd_ = logstd;}

		virtual double value() { return(rng_->lognormal(logavg_, logstd_)); };
        //double* avgp() { return &avg_; };
        //double* stdp() { return &std_; };
        //virtual inline double avg()     { return avg_; };
        double std()     { return logstd_; };
        //void setavg(double d)    { avg_ = d; };
        //void setstd(double d)    { std_ = d; };
		//double mean() { return exp(avg_ + std_*std_/2); };
		double avg() { return exp(logavg_ + logstd_*logstd_/2); };
private:
        double logavg_;
        double logstd_;
};

double randouble(double low, double high);
int    randint(int low, int high);
unsigned int rand32();

#endif

