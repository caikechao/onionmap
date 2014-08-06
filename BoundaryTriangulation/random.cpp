      
#include "random.h"

RNG* RNG::default_ = new RNG(211); //21


double randouble(double low, double high)
{
    UniformRandVar urand(low, high);
    return urand.value();
}

int randint(int low, int high)
{
    int i;
    UniformRandVar urand(0, 1);
    
    if(low >= high) {
        i = low;
    }
    else {
        i = (int)(urand.value() * (high - low + 1) + low);
        if(i > high) i = high;
    }
    
    return i;
}

unsigned int rand32()
{
    return RNG::default_->uniform32();
}
