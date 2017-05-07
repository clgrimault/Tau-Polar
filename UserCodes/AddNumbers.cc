


#include "AddNumbers.h"

AddNumbers::AddNumbers ()
: _a(0), _b(0)
{
}

AddNumbers::~AddNumbers ()
{
}

void AddNumbers::setA (int a)
{
        _a = a;
}

void AddNumbers::setB (int b)
{
        _b = b;
}

int AddNumbers::getA () const
{
        return _a;
}

int AddNumbers::getB () const
{
        return _b;
}

int AddNumbers::getSum () const
{
        return _a + _b;
}
