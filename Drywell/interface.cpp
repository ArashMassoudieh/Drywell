#include "interface.h"

Interface::Interface()
{

}

double Interface::getValue(const string quan) const
{
    if (quants.count(quan)>0)
        return quants.at(quan);
    else
        return 0;
}

Interface::~Interface()
{

}
Interface::Interface(const Interface &RHS)
{
    quants = RHS.quants;
}

Interface& Interface::operator=(const Interface &RHS)
{
    quants = RHS.quants;
    return *this;
}



