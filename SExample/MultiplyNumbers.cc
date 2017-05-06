


#include "MultiplyNumbers.h"

MultiplyNumbers::MultiplyNumbers ()
: _a(0), _b(0)
{
}

MultiplyNumbers::~MultiplyNumbers ()
{
}

void MultiplyNumbers::setA (int a)
{
        _a = a;
}

void MultiplyNumbers::setB (int b)
{
        _b = b;
}

int MultiplyNumbers::getA () const
{
        return _a;
}

int MultiplyNumbers::getB () const
{
        return _b;
}

int MultiplyNumbers::getProduct () const
{
        return _a * _b;
}
