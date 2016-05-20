
#ifndef quaternion_hpp
#define quaternion_hpp

#include <stdio.h>

class quaternion{

public:
    
    double a, x, y, z;  // quaternion all components
    double* V;          // quaternion vectorial part
    
    // ---- Costructors
    quaternion(double a, double x, double y, double z);
    quaternion(double a, double* V);
    quaternion();
    
    // ---- Public Methods
    quaternion operator*(const quaternion& b);
    quaternion operator+(const quaternion& b);
    quaternion operator-(const quaternion& b);
    quaternion operator!();
    
    double M();
    void print();
    
    // ---- friend not member function
    // (otherwise int*quaternion is not allowed)
    // (otherwise int+quaternion is not allowed)
    
    friend quaternion operator*(const quaternion &Q, const double &d)
    {
        return quaternion(Q.a * d, Q.V[0]*d, Q.V[1]*d, Q.V[2]*d);
    }
    
    friend quaternion operator*(const double &d, const quaternion &Q)
    {
        return Q*d;
    }
    
    friend quaternion operator+(const double &d, const quaternion &Q)
    {
        return quaternion(Q.a + d, Q.V[0],Q.V[1],Q.V[2]);
    }
    
    friend quaternion operator+(const quaternion &Q, const double &d)
    {
        return d + Q;
    }
    
    
//static:
    
    // ---- Public Static Methods
    static void    sPrint(quaternion q);
    static double* sProd(double a, double b[3]);
    static double  dProd(double a[3], double b[3]);
    static double* cProd (double a[3], double b[3]);
    static double* sum(double a[3], double b[3]);
    static double* diff(double a[3], double b[3]);
    static double  normPow2(double a[3]);
    static double  norm(double a[3]);
    static double* dirCos(double v[3]);

private:
    
protected:
    
};

#endif
