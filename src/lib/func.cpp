#include <iostream> //cout
#include <fstream>
#include <string>
#include <math.h> //log, exp
#include <algorithm> //min

using namespace std;

int sum(int *v, int len)
{
    int s=0;
    for(int i=1; i<=len; i++) s += v[i];
    return s;
}

double logsum(double a, double b)
{
    if(a==-1)
        return b;
    else if(b==-1)
        return a;
    else if(a>b)
        return a + log(1 + exp(b-a));
    else
        return b + log(1 + exp(a-b));
}

double logsum2(double a, double b)
{
    if(a==1000)
        return b;
    else if(b==1000)
        return a;
    else if(a>b)
        return a + log(1 + exp(b-a));
    else
        return b + log(1 + exp(a-b));
}

double lognchoosek(int nn, int k)
{
    double s=0;
    for (int i=1; i<=k; i++) s = s + log(nn-k+i)-log(i);
    return s;
}

int first_geq(unsigned int* I_start, int m, int h)
{
    int L = 0, R = m;
    int i;
    while(L<=R)
    {
        i = (L+R) / 2;
        if(I_start[i]<h)
            L = i + 1;
        else if(I_start[i]>h)
            R = i - 1;
        else
            return i;
    }
    return L;  
}

int first_g(unsigned int* I_end, int m, int h)
{
    int L = 0, R = m;
    int i;
    while(L<=R)
    {
        i = (L+R) / 2;
        if(I_end[i]<h)
            L = i + 1;
        else if(I_end[i]>h)
            R = i - 1;
        else
            return i+1;
    }
  return L;
}

int c(int li, int h, int m, unsigned int* I2_start, unsigned int* I2_end)
{//unsure about correctness of this as 10000scale gives a ratio slightly larger than one 
    //change it to use for preprocessing
    int i, j;
    i = first_geq(I2_start, m, h);
    j = first_g(I2_end, m, h-li);
    if (j>m || i==j) //first condition is necessary?
        return 0;
    int num = i-j;
    int z = ceil(double(I2_end[i-1]-I2_start[i-1])/2);
    z = 1;
    if(min(I2_end[i-1], unsigned(h))<z+max(I2_start[i-1], unsigned(h-li)))
        num--;
    z = ceil(double(I2_end[j]-I2_start[j])/2);
    z = 1;
    if(min(I2_end[j], unsigned(h))<z+max(I2_start[j], unsigned(h-li)))
        num--;
    if(num>0)
        return num;
    return 0;
}

int Indicator(int xj, unsigned int h, int n, unsigned int* I1_start, unsigned int* I1_end)
{
    int i, j;
    i = first_geq(I1_start, n, h);
    j = first_g(I1_end, n, h-xj);
    if (j>n || i==j) //first condition is necessary?
        return 0;
    int z = ceil(double(xj)/2);
    z = 1;
    for(int k=j; k<i; k++)
    {
        if(min(I1_end[k], h)>=z+max(I1_start[k], h-xj))
            return 1;
    }
    return 0;
}

int f(int h, int m, unsigned int* I2_start, unsigned int* I2_end)
{
    int i;
    i = first_geq(I2_start, m, h);
    if (i<1)
        return 0;
    int z = ceil(double(I2_end[i-1]-I2_start[i-1])/2);
    z = 1;
    if (I2_start[i-1]<=h-z && I2_end[i-1]>=h+z)
        return 1;
    return 0;
}

int overlap(int* I2_l, unsigned int* I2_end, unsigned int* I1_start, unsigned int* I1_end, int n, int m)
{
    int N_intersect = 0;
    for(int j=1; j<=m; j++)
    {
        N_intersect += Indicator(I2_l[j], I2_end[j], n, I1_start, I1_end);
    }
    return N_intersect;
}

double variance(int* I_l, int n_int)
{
    double l_mean = double(sum(I_l,n_int)) / n_int;
    double l_var = 0;

    for(int i=1; i<=n_int; i++)
    {
        l_var += pow(I_l[i]-l_mean, 2);
    }

    l_var /= n_int-1;
    return l_var;
}

void write2file(string output, double* p_value_dp, int m)
{
    ofstream myfile (output);
    if (myfile.is_open()) {
        myfile << "k" << "\t" << "DP" << endl;
        for(int k=0; k<=m; k++){
            myfile << k << "\t" << p_value_dp[k] << endl;
        }
        myfile.close();
        }
    else cout << "Unable to open file";
}

void write2file(string output, double* p_value_dp, double* p_value_poisbinom, int m)
{
    ofstream myfile (output);
    if (myfile.is_open()) {
        myfile << "k" << "\t" << "DP" << "\t" << "PB" << endl;
        for(int k=0; k<=m; k++){
            myfile << k << "\t" << p_value_dp[k] << "\t" << p_value_poisbinom[k] << endl;
        }
        myfile.close();
        }
    else cout << "Unable to open file";
}

void write2file(string output, double* p_value_dp, double* p_value_poisbinom, double* p_value_hypergeo, int m)
{
    ofstream myfile (output);
    if (myfile.is_open()) {
        myfile << "k" << "\t" << "DP" << "\t" << "PB" << "\t" << "HG" <<endl;
        for(int k=0; k<=m; k++){
            myfile << k << "\t" << p_value_dp[k] << "\t" << p_value_poisbinom[k] << "\t" << p_value_hypergeo[k] <<endl;
        }
        myfile.close();
        }
    else cout << "Unable to open file";
}