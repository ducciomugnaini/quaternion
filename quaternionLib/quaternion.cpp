
// shift + alt + cmd + (<-|->) fold/unfold

#include "iostream"
#include "math.h"
#include "quaternion.hpp"

typedef quaternion qtn;

// ----------------------------------------------------------------------------- costructor

quaternion::quaternion(double qa, double qx, double qy, double qz)
{
    a = qa;
    x = qx;
    y = qy;
    z = qz;
    
    V = new double[3];
    V[0] = x;
    V[1] = y;
    V[2] = z;
}

quaternion::quaternion(double qa, double* Va)
{
    a = qa;
    x = Va[0];
    y = Va[1];
    z = Va[2];
    
    V = new double[3];
    V[0] = x;
    V[1] = y;
    V[2] = z;
}

// default costructor for make array of quaternions
quaternion::quaternion(){}

// ----------------------------------------------------------------------------- public method

void quaternion::print()
{
    std::cout << "("<< a << ", [" << x <<"i, " << y << "j, " << z <<"k])\n";
}

// Overload + operator: QUATERNION + QUATERINON
quaternion quaternion::operator+(const quaternion& b)
{
    double a = this->a + b.a;
    double x = this->x + b.x;
    double y = this->y + b.y;
    double z = this->z + b.z;
    
    quaternion qRes (a,x,y,z);
    
    return qRes;
}

// Overload + operator: QUATERNION - QUATERINON
quaternion quaternion::operator-(const quaternion& b)
{
    double a = this->a - b.a;
    double x = this->x - b.x;
    double y = this->y - b.y;
    double z = this->z - b.z;
    
    quaternion qRes (a,x,y,z);
    
    return qRes;
}

// Overload * operator: QUATERNION * QUATERNION
quaternion quaternion::operator*(const quaternion& b)
{
    double a = (this->a * b.a) - quaternion::dProd(this->V, b.V);
    
    double* ba = quaternion::sProd(b.a, this->V);
    double* ab  = quaternion::sProd(this->a, b.V);
    double* cP = quaternion::cProd(this->V, b.V);
    
    double* Vres = quaternion::sum(quaternion::sum(ba, ab),cP);
    
    quaternion qRes (a,Vres[0],Vres[1],Vres[2]);
    
    return qRes;
}

// overloaded - operator:  conj(QUATERNION)
quaternion quaternion::operator!()
{
    return quaternion(a, -x, -y, -z);
}

// compute the magnitude of this quaternion
double quaternion::M(){
    double magn = pow(this->a,2);
    magn = magn + quaternion::normPow2(this->V);
    return magn;
}

// ----------------------------------------------------------------------------- public static method

void quaternion::sPrint(quaternion q)
{
    std::cout << "("<< q.a << ", [" << q.x <<"i, " << q.y << "j, " << q.z <<"k])\n";
}

// dot product between two arrays
double quaternion::dProd(double a[3], double b[3])
{
    double res = 0;
    for (int i=0; i<=2; i++) {
        res = res + (a[i]*b[i]);
    }
    return res;
}

// cross product between two arrays
double* quaternion::cProd(double a[], double b[])
{
    double *res = new double[3];
    
    res[0] = (a[1]*b[2]) - (a[2]*b[1]);
    res[1] = -((a[0]*b[2]) - (a[2]*b[0]));
    res[2] = (a[0]*b[1]) - (a[1]*b[0]);
    
    return res;
}

// scalar product between two arrays
double* quaternion::sProd(double a, double b[3])
{
    double* res = new double[3];
    
    res[0] = a*b[0];
    res[1] = a*b[1];
    res[2] = a*b[2];
    
    return res;
}

// sum between two arrays
double* quaternion::sum(double a[3], double b[3])
{
    double* c = new double[3];
    
    for (int i=0; i<3; i++) {
        c[i] = a[i] + b[i];
    }
    return c;
}

// diff between two arrays
double* quaternion::diff(double a[3], double b[3])
{
    double* c = new double[3];
    
    for (int i=0; i<3; i++) {
        c[i] = a[i] - b[i];
    }
    return c;
}

// |v|^2 of v array
double quaternion::normPow2(double a[3])
{
    double norm = 0;
    for (int i=0; i<3; i++) {
        norm = norm + pow(a[i], 2);
    }
    return norm;
}

double quaternion::norm(double a[3])
{
    double normA = quaternion::normPow2(a);
    return sqrt(normA);
}

double* quaternion::dirCos(double v[3])
{
    double* res = new double[3];
    
    double normV = qtn::norm(v);
    
    res[0] = v[0]/normV;
    res[1] = v[1]/normV;
    res[2] = v[2]/normV;
    
    return res;
}

// ----------------------------------------------------------------------------- private method

