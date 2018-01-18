//
// Created by ggl on 2018/1/9.
//

#ifndef INC_5G_LDPC_LDPC_H
#define INC_5G_LDPC_LDPC_H

void LDPC_encode(int *ipInforBitsNoCRC,int* wTransBits,int iInfLength, int ulTransLength);

int LDPC_decode(double *dpLLR, int iCodeLength, int iInfLength, double *bypInforbits, double *OutputLLRWithoutChannel);
#endif //INC_5G_LDPC_LDPC_H
