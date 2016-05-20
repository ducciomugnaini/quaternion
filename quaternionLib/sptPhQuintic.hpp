
#ifndef sptPhQuintic_hpp
#define sptPhQuintic_hpp

#include <stdio.h>
#include "quaternion.hpp"

typedef quaternion qtn;

typedef struct XYZCoordinates{
    
    double X[5000];
    double Y[5000];
    double Z[5000];
    
} XYZCoos;

typedef struct differentialMeasures{
    
    // ---- Container of differential measures at time t
    //      (pag. 178, 601)
    
    double T[3];    // Unit tangent
    double P[3];    // Unit principal normal
    double B[3];    // Unit binormal
    
    double K;       // Curvature
    double Tau;     // Torsion
    double RR;      // Rate of Rotation
    
} diffMeas;

class sptPhQuintic{
    
public:
    
    quaternion A[3];
    quaternion p[6];    // Bezier control points of PH quintic
    double sigma[5];    // parametric speed Bernstein coefficients
    double s[6];        // arc length Bernstein coefficients
    double w[5];        // weigths of tangent indicatrix
    quaternion wt[5];   // weigths and control points of tangent indicatrix
    std::string type;   // type of ph quintic
    double axis[3];     // axis for helical quintic ph
    
    // ---- Constructors
    sptPhQuintic(quaternion A0, quaternion A1, quaternion A2, double p0[3]);
    sptPhQuintic(qtn A0, qtn A2, double p0[3], double c0, double c2);
    sptPhQuintic(double pi[3], double pf[3], double di[3], double df[3], double tht0, double tht2);
    
    // ---- Public Methods
    XYZCoos evaluatePH(int nStep);
    XYZCoos evaluatePHonT(double t[], int tLen);
    quaternion evalAt(double t);
    quaternion evalDAt(double t);
    quaternion evalDDAt(double t);
    diffMeas evalDifferentialMeasures(double t);
    
    void printInfo();
    void plotOnFile(std::string fileName, int nStep, char* rwMode);
    void plotTanIndicatrixOnFile(std::string fileName, int nStep);
    void plotFrenetSerretFrame(std::string fileName, double t[], int tLen);
    void plot_K_T_RR(std::string fileName, double t[], int tLen);
    
//static:
    
    // ---- Public Static Methods
    // static double* sProd(double a, double b[3]);     //syntax example
    
private:
    
    void init(qtn A0, qtn A1, qtn A2, double p0[3]);
    
protected:
    
    double beval(double b[11], int n, double t);
};

#endif /* sptPhQuintic_hpp */
