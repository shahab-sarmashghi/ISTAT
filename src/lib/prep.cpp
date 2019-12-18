#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>       // stoi, stoul     
#include <sstream>
#include <math.h>

using namespace std;

int scan_input(string path, string id)
{
    int n_int = 0;
    ifstream yourfile (path);
    if (yourfile.is_open())
    {
        string line;
        string token;
        n_int = 0;        
        while(getline(yourfile, line))
        {
            stringstream iss;
            iss << line;
            getline(iss, token, '\t');
            if(token==id)
                n_int++; 
        }
        yourfile.close();
    }
    else cout << "Unable to open file" << endl;
    return n_int;  
}

void read_input(unsigned int* I_start, unsigned int* I_end, string path, string id)
{
    I_start[0] = 0;
    I_end[0] = 0;
    ifstream yourfile (path);
    if (yourfile.is_open())
    {
        string line;
        string token;        
        int i = 1;
        while(getline(yourfile, line))
        {
            stringstream iss;
            iss << line;
            getline(iss, token, '\t');
            if(token==id)
            {
                getline(iss, token, '\t');
                I_start[i] = stoul(token);
                getline(iss, token);
                I_end[i] = stoul(token);
                i++; 
            }
        }
        yourfile.close();
    }
    else cout << "Unable to open file" << endl;
}

int merge_adj(unsigned int* I_start, unsigned int* I_end, int n_int)
{
    if (n_int<2)
        return n_int;
    int counter = 1;
    for (int i=2; i<=n_int; i++)
    {
        if (I_end[counter]>=I_start[i])
            I_end[counter] = max(I_end[i], I_end[counter]);
        else
        {
            counter++;
            I_start[counter] = I_start[i];
            I_end[counter] = I_end[i];
        }
    }
    return counter;
}

void compute_length(int* I_l, int* I_bar_l, unsigned int* I_start, unsigned int* I_end, int n_int)
{
    for(int i=1; i<=n_int; i++)
    {
        I_l[i] = I_end[i]-I_start[i];
        I_bar_l[i] = I_start[i]-I_end[i-1];
    }
}

void scale_interval(unsigned int* I_start, unsigned int* I_end, int* I_l, int* I_bar_l, int n_int, int scale)
{
    I_start[0] = 0;
    I_end[0] = 0;

    for(int i=1; i<=n_int; i++)
    {
        I_bar_l[i] = floor(double(I_bar_l[i]) / scale);
        I_l[i] = ceil(double(I_l[i]) / scale);
        I_start[i] = I_bar_l[i] + I_end[i-1];
        I_end[i] = I_start[i] + I_l[i];
    }
}