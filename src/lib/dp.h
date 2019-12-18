#ifndef DP_H
#define DP_H

void DP(double* output, int* I1_l, unsigned int* I2_start,
        int* I2_l, unsigned int* I2_end, int n, int m, unsigned int g,
        int mode=0);
// computes the p-value vector using dynamic programming algorithm

#endif