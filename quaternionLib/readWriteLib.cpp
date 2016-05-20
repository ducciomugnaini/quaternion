
// shift + alt + cmd + (<-|->) fold/unfold

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "readWriteLib.hpp"

// ----------------------------------------------------------------------------- public static method

void readWriteLib::printXYCoo(std::string fileName, double x[], double y[],  int XYLength, char* rwMode)
{
    FILE * myFile;
    int i;
    
    const char * cFileName = fileName.c_str();

    myFile = fopen(cFileName, rwMode);
    for (i = 0 ; i < XYLength ; i++)
    {
        fprintf (myFile, "%3.5f %3.5f\n", x[i], y[i]);
    }
    fprintf (myFile, "\n");
    
    fclose (myFile);
}

void readWriteLib::printXYZCoo(std::string fileName, double x[], double y[], double z[],  int XYLength,
                               char* rwMode)
{
    FILE * myFile;
    int i;
    
    const char * cFileName = fileName.c_str();
    
    myFile = fopen(cFileName, rwMode);
    for (i = 0 ; i < XYLength ; i++)
    {
        fprintf (myFile, "%3.5f %3.5f %3.5f\n", x[i], y[i], z[i]);
    }
    fprintf (myFile, "\n");
    
    fclose (myFile);
}


