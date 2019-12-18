#ifndef FUNC_H
#define FUNC_H

int sum(int *v, int len);
// returns sum of v[1]:v[len]

double logsum(double a, double b);
// returns the result of sum of a and b in real scale

double logsum2(double a, double b);
// returns the result of sum of a and b in real scale

double lognchoosek(int nn, int k);
// returns nn choose k in logarithmic scale

int first_geq(unsigned int* I_start, int m, int h);
// finds the first element in I_start that is geq to h

int first_g(unsigned int* I_end, int m, int h);
// finds the first element in I_end that is greater than h

int c(int li, int h, int m, unsigned int* I2_start, unsigned int* I2_end);
// counts the number of intervals in I_2 intersected by interval with length
// li that ends at h

int Indicator(int xj, unsigned int h, int n, unsigned int* I1_start,
              unsigned int* I1_end);
// indicates whether interval with length xj ending at h is intersected by
// any interval in I_1

int f(int h, int m, unsigned int* I2_start, unsigned int* I2_end);
// indicates whether position h is covered by any interval in I_2

int overlap(int* I2_l, unsigned int* I2_end, unsigned int* I1_start, unsigned int* I1_end, int n, int m);
// finds the number of intervals in I2 intersected by intervals in I1

double variance (int* I_l, int n_int);
// compute the variance in I_l[1]:I_l[n_int]

void write2file(std::string output, double* p_value_dp, int m);
// writes tab-delimited p-values to output file

void write2file(std::string output, double* p_value_dp, double* p_value_poisbinom, int m);
// writes tab-delimited p-values to output file

void write2file(std::string output, double* p_value_dp, double* p_value_poisbinom, double* p_value_hypergeo, int m);
// writes tab-delimited p-values to output file
 
#endif