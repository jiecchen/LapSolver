#pragma once

class LinearOperator
{
public:
    virtual ~LinearOperator() {}

    virtual void apply(double *) = 0;
    virtual int getDimension() const = 0;
};