
// shift + alt + cmd + (<-|->) fold/unfold
#include "iostream"
#include "math.h"

#include "sptPhQuintic.hpp"
#include "quaternion.hpp"
#include "readWriteLib.hpp"

#define MY_DEBUG 0

typedef quaternion qtn;

// ----------------------------------------------------------------------------- constructor

sptPhQuintic::sptPhQuintic(qtn A0, qtn A1, qtn A2, double p0[3])
{
    this->type = "General Ph Quintic";
 
    // TODO: CHANGE WITH INIT
    
    // init - building quaternion
    
    A[0] = A0;
    A[1] = A1;
    A[2] = A2;
    
    quaternion i = quaternion(0., 1., 0., 0.);
    
    // init - control points
    
    p[0] = quaternion(0, p0[0], p0[1], p0[2]);
    
    p[1] = p[0] + (1.0/5.0)  * (  A[0]*i*(!A[0]) );
    
    p[2] = p[1] + (1.0/10.0) * ( (A[0]*i*(!A[1])) + (A[1]*i*(!A[0])) );
    
    p[3] = p[2] + (1.0/30.0) * ( (A[0]*i*(!A[2])) + 4*(A[1]*i*(!A[1])) + (A[2]*i*(!A[0])) );
    
    p[4] = p[3] + (1.0/10.0) * ( (A[1]*i*(!A[2])) + (A[2]*i*(!A[1])) );
    
    p[5] = p[4] + (1.0/5.0)  * (  A[2]*i*(!A[2]) );
    
    
    // init - sigma coefficients
    
    quaternion s1Q = 1./2.*( (A[0]*(!A[1])) + (A[1]*(!A[0])));
    quaternion s2Q = 1./6.*( (A[0]*(!A[2])) + 4.*A[1].M() + A[2]*(!A[0]) );
    quaternion s3Q = 1./2.*(   A[1]*(!A[2]) + A[2]*(!A[1]) );
    
    sigma[0] = A[0].M();
    sigma[1] = s1Q.a;
    sigma[2] = s2Q.a;
    sigma[3] = s3Q.a;
    sigma[4] = A[2].M();
    
    // init - arc length coefficients
    
    s[0] = 0;
    double sigmaSum = 0;
    for (int i=1; i<=5; i++) {
        sigmaSum = sigmaSum + sigma[i-1];
        s[i] = 1./5. * sigmaSum;
    }
    /*
    s[1] = sigma[0];
    s[2] = 1./5. * (sigma[0] + sigma[1]);
    s[3] = 1./5. * (sigma[0] + sigma[1] + sigma[2]);
    s[4] = 1./5. * (sigma[0] + sigma[1] + sigma[2] + sigma[3]);
    s[5] = 1./5. * (sigma[0] + sigma[1] + sigma[2] + sigma[3] + sigma[4]);
    */
    
    // init - [weights] of tangent indicatrix (pag. 493)
    
    qtn w0 = A[0]*(!A[0]);
    qtn w1 = 1./2. * ( A[0]*(!A[1]) + A[1]*(!A[0]) );
    qtn w2 = 1./6. * ( A[0]*(!A[2]) + 4*A[1]*(!A[1]) + A[2]*(!A[0]) );
    qtn w3 = 1./2. * ( A[1]*(!A[2]) + A[2]*(!A[1]) );
    qtn w4 = A[2]*(!A[2]);
    
    this->w[0] = w0.a;
    this->w[1] = w1.a;
    this->w[2] = w2.a;
    this->w[3] = w3.a;
    this->w[4] = w4.a;
    
    // init - [weights * control points] of tangent indicatrix (pag. 493)
    
    this->wt[0] = A[0]*i*(!A[0]);
    this->wt[1] = 1./2. * ( A[0]*i*(!A[1]) + A[1]*i*(!A[0]) );
    this->wt[2] = 1./6. * ( A[0]*i*(!A[2]) + 4*A[1]*i*(!A[1]) + A[2]*i*(!A[0]) );
    this->wt[3] = 1./2. * ( A[1]*i*(!A[2]) + A[2]*i*(!A[1]) );
    this->wt[4] = A[2]*i*(!A[2]);
    
}

sptPhQuintic::sptPhQuintic(qtn A0, qtn A2, double p0[3], double c0, double c2) :
sptPhQuintic::sptPhQuintic(A0, c0*A0 + c2*A2, A2, p0)
{
    // TODO: CHANGE WITH INIT
    
    this->type = "General Helical Ph Quintic";
    
    // init - axis
    
    double* num = qtn::sProd(A[0].a, A[2].V);
    num = qtn::diff(num, qtn::sProd(A[2].a, A[0].V));
    num = qtn::sum(num, qtn::cProd(A[0].V, A[2].V));
    
    axis[0] = *(num);
    axis[1] = *(num+1);
    axis[2] = *(num+2);
    
    double denum = qtn::norm(num);
    
    for (int i=0; i<3; i++) {
        axis[i] = (*(num+i)) / denum;
    }
    
}

sptPhQuintic::sptPhQuintic(double pi[3], double pf[3], double di[3], double df[3], double tht0, double tht2)
{
    this->type = "Hermite Ph Quintic Interpolant";
    
    // ---- direction consines di[]
    double* di_dc = qtn::dirCos(di);
    double la_i = di_dc[0];
    double mu_i = di_dc[1];
    double nu_i = di_dc[2];
    
    // ---- direction consines df[]
    double* df_dc = qtn::dirCos(df);
    double la_f = df_dc[0];
    double mu_f = df_dc[1];
    double nu_f = df_dc[2];
    
    double di_norm = qtn::norm(di);
    double df_norm = qtn::norm(df);
    
    // ---- A0 components computation in accord to pag.597 (28.7)
    double A0_coeff = sqrt( 0.5 * ( 1.0 + la_i ) * di_norm );
    double A0_a     = -sin(tht0);
    double A0_i     = cos(tht0);
    double A0_j     = ( mu_i*cos(tht0) + nu_i*sin(tht0) ) / (1.0 + la_i);
    double A0_k     = ( nu_i*cos(tht0) - mu_i*sin(tht0) ) / (1.0 + la_i);
    
    A0_a = A0_coeff * A0_a;
    A0_i = A0_coeff * A0_i;
    A0_j = A0_coeff * A0_j;
    A0_k = A0_coeff * A0_k;
    
    qtn A0 (A0_a, A0_i, A0_j, A0_k);
    
    // ---- A2 components computation in accord to pag.597 (28.7)
    double A2_coeff = sqrt( 0.5 * (1.0 + la_f) * df_norm );
    double A2_a     = -sin(tht2);
    double A2_i     = cos(tht2);
    double A2_j     = ( mu_f*cos(tht2) + nu_f*sin(tht2) ) / (1.0 + la_f);
    double A2_k     = ( nu_f*cos(tht2) - mu_f*sin(tht2) ) / (1.0 + la_f);
    
    A2_a = A2_coeff * A2_a;
    A2_i = A2_coeff * A2_i;
    A2_j = A2_coeff * A2_j;
    A2_k = A2_coeff * A2_k;
    
    qtn A2 (A2_a, A2_i, A2_j, A2_k);
    
    // ---- C component computation (pag.597 - 28.9)
    quaternion i = quaternion(0, 1, 0, 0);
    double* Ca  = qtn::sProd(120.0, qtn::diff(pf, pi));
    double* Cb  = qtn::sProd(15.0, qtn::sum(di, df));
    qtn Cd      = (A0 * i*(!A2)) + (A2 * i*(!A0));
    Cd = 5.0 * Cd;
    
    // printf("\nca: %f %f %f\n",Ca[0], Ca[1], Ca[2]);
    // printf("cb: %f %f %f\n",  Cb[0], Cb[1], Cb[2]);
    // Cd.print();
    
    double C[3] = {0, 0, 0};
    C[0] = Ca[0] - Cb[0] + Cd.x;
    C[1] = Ca[1] - Cb[1] + Cd.y;
    C[2] = Ca[2] - Cb[2] + Cd.z;
    
    // ---- A1 components computation (pag.600 - 28.15)
    double* CDirCos = qtn::dirCos(C);
    double la = CDirCos[0];
    double mu = CDirCos[1];
    double nu = CDirCos[2];
    
    qtn     A1a = (-3.0/4.0) * (A0 + A2);
    double  A1b = sqrt(0.5 * (1.0+la) * qtn::norm(C) ) / 4.0;
    qtn     A1c (1., 0., -nu/(1+la), mu/(1+la));
    qtn     A1 = A1a + (A1b * A1c);
    
    // ---- call private init to construct the quintic ph object
    this->init(A0, A1, A2, pi);
    
}

// ----------------------------------------------------------------------------- private method

void sptPhQuintic::init(qtn A0, qtn A1, qtn A2, double p0[3])
{
    
    // init - building quaternion
    
    A[0] = A0;
    A[1] = A1;
    A[2] = A2;
    
    quaternion i = quaternion(0, 1, 0, 0);
    
    // init - control points
    
    p[0] = quaternion(0, p0[0], p0[1], p0[2]);
    
    p[1] = p[0] + 1.0/5.0 *(  A[0]*i*(!A[0]));
    
    p[2] = p[1] + 1.0/10.0*( (A[0]*i*(!A[1])) + (A[1]*i*(!A[0])) );
    
    p[3] = p[2] + 1.0/30.0*( (A[0]*i*(!A[2])) + 4*(A[1]*i*(!A[1])) + (A[2]*i*(!A[0])) );
    
    p[4] = p[3] + 1.0/10.0*( (A[1]*i*(!A[2])) + (A[2]*i*(!A[1])) );
    
    p[5] = p[4] + 1.0/5.0 *(  A[2]*i*(!A[2]));
    
    
    // init - sigma coefficients
    
    quaternion s1Q = 1./2.*( (A[0]*(!A[1])) + (A[1]*(!A[0])));
    quaternion s2Q = 1./6.*( (A[0]*(!A[2])) + 4.*A[1].M() + A[2]*(!A[0]) );
    quaternion s3Q = 1./2. * ( A[1]*(!A[2]) + A[2]*(!A[1]) );
    
    sigma[0] = A[0].M();
    sigma[1] = s1Q.a;
    sigma[2] = s2Q.a;
    sigma[3] = s3Q.a;
    sigma[4] = A[2].M();
    
    // init - arc length coefficients
    
    s[0] = 0;
    double sigmaSum = 0;
    for (int i=1; i<=5; i++) {
        sigmaSum = sigmaSum + sigma[i-1];
        s[i] = 1./5. * sigmaSum;
    }
    /*
     s[1] = sigma[0];
     s[2] = 1./5. * (sigma[0] + sigma[1]);
     s[3] = 1./5. * (sigma[0] + sigma[1] + sigma[2]);
     s[4] = 1./5. * (sigma[0] + sigma[1] + sigma[2] + sigma[3]);
     s[5] = 1./5. * (sigma[0] + sigma[1] + sigma[2] + sigma[3] + sigma[4]);
     */
    
    // init - [weights] of tangent indicatrix (pag. 493)
    
    qtn w0 = A[0]*(!A[0]);
    qtn w1 = 1./2. * ( A[0]*(!A[1]) + A[1]*(!A[0]) );
    qtn w2 = 1./6. * ( A[0]*(!A[2]) + 4*A[1]*(!A[1]) + A[2]*(!A[0]) );
    qtn w3 = 1./2. * ( A[1]*(!A[2]) + A[2]*(!A[1]) );
    qtn w4 = A[2]*(!A[2]);
    
    this->w[0] = w0.a;
    this->w[1] = w1.a;
    this->w[2] = w2.a;
    this->w[3] = w3.a;
    this->w[4] = w4.a;
    
    // init - [weights * control points] of tangent indicatrix (pag. 493)
    
    this->wt[0] = A[0]*i* (!A[0]);
    this->wt[1] = 1./2. * ( A[0]*i*(!A[1]) + A[1]*i*(!A[0]) );
    this->wt[2] = 1./6. * ( A[0]*i*(!A[2]) + 4*A[1]*i*(!A[1]) + A[2]*i*(!A[0]) );
    this->wt[3] = 1./2. * ( A[1]*i*(!A[2]) + A[2]*i*(!A[1]) );
    this->wt[4] = A[2]*i* (!A[2]);
    
}

// ----------------------------------------------------------------------------- public method

void sptPhQuintic::printInfo()
{
    printf("\n~~\n");
    
    std::cout << this->type << " :: " << "Spatial Quintic PH Infos:\n-\n";
    
    for (int i=0; i<3; i++) {
        printf("A[%d]: ",i);
        A[i].print();
    }
    
    printf("-\n");
    
    for (int i=0; i<6; i++) {
        printf("p[%d]: ",i);
        p[i].print();
    }
    
    printf("-\n");
    
    printf("sigma: ");
    for (int i=0; i<=4; i++) {
        printf("%3.3f, ",sigma[i]);
    }
    
    printf("\n-\n");
    
    printf("s: ");
    for (int i=0; i<=5; i++) {
        printf("%3.3f, ",s[i]);
    }
    
    printf("\n~~\n");
}

XYZCoos sptPhQuintic::evaluatePH(int nStep)
{
    XYZCoos xyzStruct;
    
    double nan = 0.0/0.0;
    double dt = 1./nStep;
    double t = 0;
    
    std::fill(xyzStruct.X, std::end(xyzStruct.X), nan);
    std::fill(xyzStruct.Y, std::end(xyzStruct.Y), nan);
    std::fill(xyzStruct.Z, std::end(xyzStruct.Z), nan);
    
    if(nStep <= 5000){
        
        double cpX[6] = {this->p[0].x, this->p[1].x, this->p[2].x, this->p[3].x, this->p[4].x, this->p[5].x};
        double cpY[6] = {this->p[0].y, this->p[1].y, this->p[2].y, this->p[3].y, this->p[4].y, this->p[5].y};
        double cpZ[6] = {this->p[0].z, this->p[1].z, this->p[2].z, this->p[3].z, this->p[4].z, this->p[5].z};
        
        for (int i=0; i<nStep; i++) {
            t = i*dt;
            xyzStruct.X[i] = beval(cpX, 5, t);
            xyzStruct.Y[i] = beval(cpY, 5, t);
            xyzStruct.Z[i] = beval(cpZ, 5, t);
        }
        
    }
    else
    {
        printf("(!) >> you are requiring the evaluation of too many points.");
    }
    
    return xyzStruct;
}

XYZCoos sptPhQuintic::evaluatePHonT(double t[], int tLen)
{
    XYZCoos xyzStruct;
    double nan = 0.0/0.0;
    std::fill(xyzStruct.X, std::end(xyzStruct.X), nan);
    std::fill(xyzStruct.Y, std::end(xyzStruct.Y), nan);
    std::fill(xyzStruct.Z, std::end(xyzStruct.Z), nan);
    
    if(tLen <= 5000){
        
        double cpX[6] = {this->p[0].x, this->p[1].x, this->p[2].x, this->p[3].x, this->p[4].x, this->p[5].x};
        double cpY[6] = {this->p[0].y, this->p[1].y, this->p[2].y, this->p[3].y, this->p[4].y, this->p[5].y};
        double cpZ[6] = {this->p[0].z, this->p[1].z, this->p[2].z, this->p[3].z, this->p[4].z, this->p[5].z};
        
        for (int i=0; i<tLen; i++) {
            xyzStruct.X[i] = beval(cpX, 5, t[i]);
            xyzStruct.Y[i] = beval(cpY, 5, t[i]);
            xyzStruct.Z[i] = beval(cpZ, 5, t[i]);
        }
        
    }
    else
    {
        printf("(!) >> you are requiring the evaluation of too many points.");
    }
    
    return xyzStruct;
}

quaternion sptPhQuintic::evalAt(double t)
{
     return (A[0]*pow(1.-t, 2.)) + (A[1]*2*(t*(1-t))) + (A[2]*pow(t, 2.));
}

quaternion sptPhQuintic::evalDAt(double t)
{
    return 2*( ((A[1]-A[0]) * (1.-t)) + ((A[2]-A[1]) * t) );
}

quaternion sptPhQuintic::evalDDAt(double t)
{
    return 2*( A[2] - (2*A[1]) + A[0] );
}

diffMeas sptPhQuintic::evalDifferentialMeasures(double t)
{
    if(MY_DEBUG)
    {
        printf(">> call::evalDifferentialMeasures( %3.5f )",t);
    }
    
    // ---- Implementation of differential measures of pag. 178, 601
    
    diffMeas diffM;
    
    // ---- computation of derivatives of curve
    qtn i      (0,1,0,0);
    qtn At   = this->evalAt(t);
    qtn DAt  = this->evalDAt(t);
    qtn DDAt = this->evalDDAt(t);
    
    qtn rDt  =  At*i*(!At);
    qtn rDDt = (DAt*i*(!At))  + (At*i*(!DAt));
    qtn rDDDt= (DDAt*i*(!At)) + (2*((!At)*i*(!DAt))) + (At*i*(!DDAt));
    
    
    // ---- computation of Frenet-Serret Frame
    qtn qT      = rDt * (1./qtn::norm(rDt.V));
    double T[3] = {qT.V[0], qT.V[1], qT.V[2]};
    
    double* rDt_x_rDDt  = qtn::cProd(rDt.V, rDDt.V);
    double* B           = qtn::dirCos(rDt_x_rDDt);
    double* P           = qtn::cProd(B, T);
    
    diffM.T[0] = T[0]; diffM.T[1] = T[1]; diffM.T[2] = T[2];
    diffM.P[0] = P[0]; diffM.P[1] = P[1]; diffM.P[2] = P[2];
    diffM.B[0] = B[0]; diffM.B[1] = B[1]; diffM.B[2] = B[2];
    
    // ---- computation of curvature
    double norm_rDt_x_rDDt  = qtn::norm(rDt_x_rDDt);
    double norm_rDt         = qtn::norm(rDt.V);
    diffM.K                 = norm_rDt_x_rDDt / (pow(norm_rDt,3.));
    
    // ---- computation of torsion
    double rD_x_rDD_d_rDDD = qtn::dProd(rDt_x_rDDt, rDDDt.V);
    double norm_rD_x_rDD   = qtn::norm(rDt_x_rDDt);
    diffM.Tau              = rD_x_rDD_d_rDDD / (pow(norm_rD_x_rDD, 2));
    
    // ---- computation of rate of rotation
    diffM.RR = sqrt( pow(diffM.K, 2.) + pow(diffM.Tau,2.) );
    
    return diffM;
}

void sptPhQuintic::plotFrenetSerretFrame(std::string fileName, double t[], int tLen)
{
    if(MY_DEBUG)
    {
        printf(">> call::plotFrenetSerretFrame");
    }
    if(tLen <= 5000){
        
        // ---- evaluation of PH on given points set t
        XYZCoos xyzPH = this->evaluatePHonT(t, tLen);
        
        // ---- init of all container of coordinates for T, P, B (Frenet-Serret)
        //      (note: alla values are of applayed T,P,B on ph)
        double nan = 0.0/0.0;
        
        XYZCoos xyzT;
        std::fill(xyzT.X, std::end(xyzT.X), nan);
        std::fill(xyzT.Y, std::end(xyzT.Y), nan);
        std::fill(xyzT.Z, std::end(xyzT.Z), nan);
        
        XYZCoos xyzP;
        std::fill(xyzP.X, std::end(xyzP.X), nan);
        std::fill(xyzP.Y, std::end(xyzP.Y), nan);
        std::fill(xyzP.Z, std::end(xyzP.Z), nan);
        
        XYZCoos xyzB;
        std::fill(xyzB.X, std::end(xyzB.X), nan);
        std::fill(xyzB.Y, std::end(xyzB.Y), nan);
        std::fill(xyzB.Z, std::end(xyzB.Z), nan);
        
        for (int i=0; i<tLen; i++) {
        
            if(MY_DEBUG)
            {
                printf("\nt: %f",t[i]);
            }
            
            // evaluation of frenet serret frame on one point t[i]
            diffMeas DFM = this->evalDifferentialMeasures(t[i]);
            
            xyzT.X[i] = xyzPH.X[i] + DFM.T[0];
            xyzT.Y[i] = xyzPH.Y[i] + DFM.T[1];
            xyzT.Z[i] = xyzPH.Z[i] + DFM.T[2];
            
            xyzP.X[i] = xyzPH.X[i] + DFM.P[0];
            xyzP.Y[i] = xyzPH.Y[i] + DFM.P[1];
            xyzP.Z[i] = xyzPH.Z[i] + DFM.P[2];
            
            xyzB.X[i] = xyzPH.X[i] + DFM.B[0];
            xyzB.Y[i] = xyzPH.Y[i] + DFM.B[1];
            xyzB.Z[i] = xyzPH.Z[i] + DFM.B[2];
            
        }
        
        char rwMode[] = "w";
        
        // ---- PHCoo application coos: print on file
        std::string PHAPPcooFileName = fileName;
        readWriteLib::printXYZCoo(PHAPPcooFileName.append("_FSF_PHAPPCoo.txt"), xyzPH.X, xyzPH.Y, xyzPH.Z, tLen, rwMode);
        
        // ---- T applyed coos: print on file
        std::string TcooFileName = fileName;
        readWriteLib::printXYZCoo(TcooFileName.append("_FSF_TCoo.txt"), xyzT.X, xyzT.Y, xyzT.Z, tLen, rwMode);
        
        // ---- P applyed coos: print on file
        std::string PcooFileName = fileName;
        readWriteLib::printXYZCoo(PcooFileName.append("_FSF_PCoo.txt"), xyzP.X, xyzP.Y, xyzP.Z, tLen, rwMode);
        
        // ---- B applyed coos: print on file
        std::string BcooFileName = fileName;
        readWriteLib::printXYZCoo(BcooFileName.append("_FSF_BCoo.txt"), xyzB.X, xyzB.Y, xyzB.Z, tLen, rwMode);
        
    }
    else
    {
        printf("(!) >> plotOnFile aborted: nStep is too big.");
    }
}

void sptPhQuintic::plot_K_T_RR(std::string fileName, double t[], int tLen)
{
    if(MY_DEBUG)
    {
        printf(">> call::plot_K_T_RR");
    }
    if(tLen <= 5000){
        
        // ---- init of all container of coordinates for T, P, B (Frenet-Serret)
        //      (note: alla values are of applayed T,P,B on ph)
        double nan = 0.0/0.0;
        
        XYZCoos xyK;
        std::fill(xyK.X, std::end(xyK.X), nan);
        std::fill(xyK.Y, std::end(xyK.Y), nan);
        std::fill(xyK.Z, std::end(xyK.Z), nan);
        
        XYZCoos xyTau;
        std::fill(xyTau.X, std::end(xyTau.X), nan);
        std::fill(xyTau.Y, std::end(xyTau.Y), nan);
        std::fill(xyTau.Z, std::end(xyTau.Z), nan);
        
        XYZCoos xyOmg;
        std::fill(xyOmg.X, std::end(xyOmg.X), nan);
        std::fill(xyOmg.Y, std::end(xyOmg.Y), nan);
        std::fill(xyOmg.Z, std::end(xyOmg.Z), nan);
        
        for (int i=0; i<tLen; i++) {
            
            if(MY_DEBUG)
            {
                printf("\nt: %f",t[i]);
            }
            
            // ---- evaluation of polynomial arc length S(t)
            double St = beval(this->s, 5, t[i]);
            
            double fracSt = St/this->s[5];
            
            // ---- evaluation of frenet serret frame on one point t[i]
            diffMeas DFM = this->evalDifferentialMeasures(fracSt);
            
            // ---- save Curvature value K(S(t)/S[5])
            xyK.X[i] = fracSt;
            xyK.Y[i] = DFM.K;
            
            // ---- save Torsion value Tau(S(t)/S[5])
            xyTau.X[i] = fracSt;
            xyTau.Y[i] = DFM.Tau;
            
            // ---- save Rate of rotation value Omega(S(t)/S[5])
            xyOmg.X[i] = fracSt;
            xyOmg.Y[i] = DFM.RR;
            
        }
        
        char rwMode[] = "w";
        
        // ---- Curvature K: print on file
        std::string KcooFileName = fileName;
        readWriteLib::printXYCoo(KcooFileName.append("_DM_Kpp.txt"), xyK.X, xyK.Y, tLen, rwMode);
        
        // ---- Torsion Tau: print on file
        std::string TaucooFileName = fileName;
        readWriteLib::printXYCoo(TaucooFileName.append("_DM_Tau.txt"), xyTau.X, xyTau.Y, tLen, rwMode);
        
        // ---- Rate of Rotation: print on file
        std::string RRcooFileName = fileName;
        readWriteLib::printXYCoo(RRcooFileName.append("_DM_Omg.txt"), xyOmg.X, xyOmg.Y, tLen, rwMode);
        
    }
    else
    {
        printf("(!) >> plotOnFile aborted: nStep is too big.");
    }
}

void sptPhQuintic::plotOnFile(std::string fileName, int nStep, char* rwMode)
{
    if(nStep <= 5000){
        
        // ---- print control points
        double cpX[6] = {this->p[0].x, this->p[1].x, this->p[2].x, this->p[3].x, this->p[4].x, this->p[5].x};
        double cpY[6] = {this->p[0].y, this->p[1].y, this->p[2].y, this->p[3].y, this->p[4].y, this->p[5].y};
        double cpZ[6] = {this->p[0].z, this->p[1].z, this->p[2].z, this->p[3].z, this->p[4].z, this->p[5].z};
        
        std::string cpFileName = fileName;
        readWriteLib::printXYZCoo(cpFileName.append("_sph5_cp.txt"), cpX, cpY, cpZ, 6, rwMode);
        
        // ---- print ph quintic
        XYZCoos xyzStruct = this->evaluatePH(nStep);
        
        std::string ph5FileName = fileName;
        readWriteLib::printXYZCoo(ph5FileName.append("_sph5.txt"), xyzStruct.X, xyzStruct.Y, xyzStruct.Z, nStep, rwMode);
        
    }
    else
    {
        printf("(!) >> plotOnFile aborted: nStep is too big.");
    }
}

void sptPhQuintic::plotTanIndicatrixOnFile(std::string fileName, int nStep)
{
    if(nStep <= 5000){
        
        double nan = 0.0/0.0;
        double dt = 1./nStep;
        double t = 0;
        
        double X[5000];
        double Y[5000];
        double Z[5000];
        double W[5000];
        
        std::fill(X, std::end(X), nan);
        std::fill(Y, std::end(Y), nan);
        std::fill(Z, std::end(Z), nan);
        
        // retrive the tangent indicatrix control points
        
        double cpX[5] = {this->wt[0].x, this->wt[1].x, this->wt[2].x, this->wt[3].x, this->wt[4].x};
        double cpY[5] = {this->wt[0].y, this->wt[1].y, this->wt[2].y, this->wt[3].y, this->wt[4].y};
        double cpZ[5] = {this->wt[0].z, this->wt[1].z, this->wt[2].z, this->wt[3].z, this->wt[4].z};
        
        double wk[5]  = {this->w[0], this->w[1], this->w[2], this->w[3], this->w[4]};
        
        for (int i=0; i<nStep; i++) {
            t = i*dt;
            
            // (pag.492 - 23.10) denumerator eval
            W[i] = beval(wk, 4, t);
            
            // (pag.492 - 23.10) numerator eval
            X[i] = beval(cpX, 4, t)/W[i];
            Y[i] = beval(cpY, 4, t)/W[i];
            Z[i] = beval(cpZ, 4, t)/W[i];
            
        }
        
        char wrMode[] = "w";
        
        // print control points
        std::string cpFileName = fileName;
        readWriteLib::printXYZCoo(cpFileName.append("_sph5_tiCp.txt"), cpX, cpY, cpZ, 5, wrMode);
        
        // print ph quintic
        std::string ph5FileName = fileName;
        readWriteLib::printXYZCoo(ph5FileName.append("_sph5_ti.txt"), X, Y, Z, nStep, wrMode);
        
        // print axis coos
        std::string axFileName = fileName;
        double ax0[1] = {this->axis[0]};
        double ax1[1] = {this->axis[1]};
        double ax2[1] = {this->axis[2]};
        readWriteLib::printXYZCoo(axFileName.append("_sph5_ax.txt"), ax0, ax1, ax2, 1, wrMode);
    }
    else
    {
        printf("(!) >> plotOnFile aborted: nStep is too big.");
    }
}


// ----------------------------------------------------------------------------- public static method


// ----------------------------------------------------------------------------- protected method

/*
 double coeffArray[11]  - array of coefficient of polynomial to evaluate
 int    iEnd            - last acceptable index on coeffArray 
                          (i.e. v[3] = {v[0] v[1] v[2]} (3 elements) => iEnd = 2)
 double t               - value of paramter t for evaluation
*/
double sptPhQuintic::beval(double coeffArray[11], int iEnd, double t)
{
    int k , r ; double blast[11] , bnext[11] , value;
    
    if ( iEnd < 0 ) return(0.0);
    if ( iEnd == 0 ) return(coeffArray[0]);
    
    for ( k = 0 ; k <= iEnd ; ++k ) {
        blast[k] = coeffArray[k];
    }
    
    for ( r = 1 ; r <= iEnd ; ++r )
    {
        for ( k = r ; k <= iEnd ; ++k )
        {
            bnext[k] = (1.-t)*blast[k-1] + t*blast[k];
        }
        
        for ( k = r ; k <= iEnd ; ++k )
        {
            blast[k] = bnext[k];
        }
    }
    value = bnext[iEnd];
    
    return(value) ;
}






