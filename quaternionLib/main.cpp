
#include <iostream>
#include "math.h"

#include "quaternion.hpp"
#include "sptPhQuintic.hpp"
#include "readWriteLib.hpp"

#define WRITE  "w "
#define APPEND 'a'

typedef quaternion qtn;

int main(int argc, const char * argv[]) {
    
    char write[]  = "w";
    char append[] = "a";
  
    // demoBranch asaddasd rovksjdskfjn
    
    // ------------------------------------------------------------------------- demo: quaternion
    /*
    qtn At (1,3,2,1);
    qtn it (0,1,0,0);
    
    qtn debug = it*(!At);
    
    debug.print();
    
    qtn debug2 = At*debug;
    debug2.print();
    
    qtn deb3 = At*it*(!At);
    deb3.print();
    */
    
    qtn A (0,2,3,4);
    qtn B (1,1,3,1);
    
    qtn C = A*B;
    
    C.print();
    
    qtn D = !B;
    D.print();
    
    printf("\n%f\n\n",B.M());
    
    
    // ---- rotation of pi along x-axis 180°
    
    XYZCoos xyzQtnRotation;
    
    double rotationAngle = (M_PI*1.95);
    int nStep = 100;
    
    double thet = rotationAngle/100;
    qtn V (0, 1, 0, 1);
    
    for (int k=0; k<nStep; k++) {
        
        qtn U (cos(thet/2.), sin(thet/2.) * 0.5774, sin(thet/2.)* 0.5774, sin(thet/2.)* 0.5774 );
        
        qtn R = U*V*(!U);
        
        xyzQtnRotation.X[k] = R.x;
        xyzQtnRotation.Y[k] = R.y;
        xyzQtnRotation.Z[k] = R.z;
        
        thet = thet + rotationAngle/nStep;
    }
    
    readWriteLib::printXYZCoo("xyzQtnRotation.txt", xyzQtnRotation.X, xyzQtnRotation.Y, xyzQtnRotation.Z, 100, write);
    
    /*
    qtn U (cos(thet2/2.), sin(thet2/2.)*1, 0, 0 );
    qtn cU = !U;
    U.print();
    cU.print();
    
    
    qtn R = ((U*V)*cU);
    
    printf("Rotated:\n");
    R.print();
    printf("\n");
    */
    
    
    // ----
    
    qtn UU (1,2,3,4);
    qtn i (0,1,0,0);
    
    qtn Ui = (UU*i);
    
    printf("\n~\n");
    UU.print();
    Ui.print();
    
    printf("\n~\n");
    
    qtn F (1,2,3,4);
    
    double a = 3.0;
    F = a*F;
    F.print();
    
    // ---- init t vector for evaluation of Frenet-Serret Frame
    
    int tLen = 100;
    double t[tLen];
    double dt = 1./tLen;
    for (int k=0; k<tLen-1; k++) {
        t[k] = dt * k;
    }
    t[tLen-1] = 1.0;
    
    // ------------------------------------------------------------------------- demo: spatial ph quintic
    
    printf("\n~\n");
    printf("Construction of spatial quintic ph by 3 building quaternions");
    
    qtn A0 (1,3,2,1);
    qtn A1 (1,2,7,5);
    qtn A2 (1,7,9,2);
    double v[3] = {8,6,9};
    sptPhQuintic ph1 (A0, A1, A2, v);
    
    ph1.printInfo();
    ph1.plotOnFile              ("ph1", 1000, write);
    ph1.plotTanIndicatrixOnFile ("ph1", 1000);
    ph1.plotFrenetSerretFrame   ("ph1", t, tLen);
    ph1.plot_K_T_RR             ("ph1", t, tLen);
    
    // ------------------------------------------------------------------------- demo: general helical ph quintic
    
    printf("\n~\n");
    printf("Construction of General helical quintic ph by 2 building quaternions and 2 coefficients");
    
    sptPhQuintic ph2Helical (A0, A2, v, 0.5, 0.7);
    
    ph2Helical.printInfo();
    ph2Helical.plotOnFile               ("pGHelical", 1000, write);
    ph2Helical.plotTanIndicatrixOnFile  ("pGHelical", 1000);
    ph2Helical.plotFrenetSerretFrame    ("pGHelical", t, tLen);
    ph2Helical.plot_K_T_RR              ("pGHelical", t, tLen);
    
    // ------------------------------------------------------------------------- demo: hermite quintic ph quintic
    
    printf("\n~\n");
    printf("Construction of Quintic ph by Hermite interpolation");
    
    double stretch = 10;
    
    double pi[3] = {0., 0., 0.};
    double pf[3] = {1.*stretch, 1.*stretch, 1.*stretch};
    
    double di[3] = {0.4, -1.5, -1.2};
    double df[3] = {-1.2, -0.6, -1.2};

    double tht0  =  -M_PI/2.0;
    double tht2  =  -M_PI/2.0;
    
    sptPhQuintic phHermite (pi, pf, di, df, tht0, tht2);
    
    phHermite.printInfo();
    phHermite.plotOnFile                ("phHermite", 1000, write);
    phHermite.plotTanIndicatrixOnFile   ("phHermite", 1000);
    phHermite.plotFrenetSerretFrame     ("phHermite", t, tLen);
    phHermite.plot_K_T_RR               ("phHermite", t, tLen);
    
    printf("\n\n");
    
    // ------------------------------------------------------------------------- demo: hermite - varying tht0 tht2
    
    // ---- varying tht0 tht2 multiwriting
    std::string inifileName = "phHermite_thtVar";
    std::string cpsFileName = inifileName;
    std::string xyzFileName = inifileName;
    remove( cpsFileName.append("_sph5_cp.txt").c_str() );
    remove( xyzFileName.append("_sph5.txt").c_str() );
    
    double thtStart = -M_PI/2.;
    double thtEnd   = +M_PI/2.;
    double incr     = M_PI/2.;
    int idPh        = 1;
    
    for (tht0 = thtStart; tht0 <= thtEnd; tht0 = tht0 + incr) {
        for (tht2 = thtStart; tht2 <= thtEnd; tht2 = tht2 + incr) {
            
            sptPhQuintic phHermTemp (pi, pf, di, df, tht0, tht2);
            
            printf("\nHERMITE PH [%d] >> theta0: %3.3f | theta2: %3.3f", idPh, tht0, tht2);
            phHermTemp.printInfo();
            
            phHermTemp.plotOnFile(inifileName, 100, append);
            
            idPh = idPh + 1;
            
        }
    }
    
    printf("\n~\n");
    return 0;
}



