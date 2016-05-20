//
//  readWriteLib.hpp
//  quaternionLib
//
//  Created by Duccio Mugnaini on 09/06/16.
//  Copyright © 2016 Duccio Mugnaini. All rights reserved.
//

#ifndef readWriteLib_hpp
#define readWriteLib_hpp

#include <stdio.h>
#include <string>

class readWriteLib{
    
    //method and variable protected
private:
    
protected:
    
public:
    
    // ---- Costructors
    
    // ---- Methods
    
    // ---- friend not member function (otherwise int*quaternion is not allowed)
    
    //static:
    static void printXYCoo(std::string fileName, double x[], double y[],  int XYLength, char* rwMode);
    static void printXYZCoo(std::string fileName, double x[], double y[], double z[], int XYLength,
                            char* rwMode);
    
};
#endif /* readWriteLib_hpp */
