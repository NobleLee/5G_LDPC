//
// Created by ggl on 2018/1/20.
//

#ifndef INC_5G_LDPC_CHANNAL_H
#define INC_5G_LDPC_CHANNAL_H

#include <stdlib.h>
#include<stdio.h>
#include <math.h>

const double _16QAM_real[16] = {0.3162, 0.3162, 0.9487, 0.9487, 0.3162, 0.3162, 0.9487, 0.9487, -0.3162, -0.3162, -0.9487, -0.9487, -0.3162, -0.3162, -0.9487, -0.9487};
const double _16QAM_imag[16] = {0.3162, 0.9487, 0.3162, 0.9487, -0.3162, -0.9487, -0.3162, -0.9487, 0.3162, 0.9487, 0.3162, 0.9487, -0.3162, -0.9487, -0.3162, -0.9487};
const int _16QAM[16][4] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1};
const int _16QAM_demodulation_index[8][8] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16, 1, 2, 5, 6, 9, 10, 13, 14,
                                             3, 4, 7, 8, 11, 12, 15, 16, 1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 10, 12, 14, 16};
const double _QPSK_real[4] = {0.7071, 0.7071, -0.7071, -0.7071};
const double _QPSK_imag[4] = {0.7071, -0.7071, 0.7071, -0.7071};
const double _64QAM_real[64] = {0.4629, 0.4629, 0.1543, 0.1543, 0.4629, 0.4629, 0.1543, 0.1543, 0.7715, 0.7715, 1.0801, 1.0801, 0.7715, 0.7715, 1.0801, 1.0801,
                                0.4629, 0.4629, 0.1543, 0.1543, 0.4629, 0.4629, 0.1543, 0.1543, 0.7715, 0.7715, 1.0801, 1.0801, 0.7715, 0.7715, 1.0801, 1.0801,
                                -0.4629, -0.4629, -0.1543, -0.1543, -0.4629, -0.4629, -0.1543, -0.1543, -0.7715, -0.7715, -1.0801, -1.0801, -0.7715, -0.7715,
                                -1.0801, -1.0801, -0.4629, -0.4629, -0.1543, -0.1543, -0.4629, -0.4629, -0.1543, -0.1543, -0.7715, -0.7715, -1.0801, -1.0801,
                                -0.7715, -0.7715, -1.0801, -1.0801};
const double _64QAM_imag[64] = {0.4629, 0.1543, 0.4629, 0.1543, 0.7715, 1.0801, 0.7715, 1.0801, 0.4629, 0.1543, 0.4629, 0.1543, 0.7715, 1.0801, 0.7715, 1.0801,
                                -0.4629, -0.1543, -0.4629, -0.1543, -0.7715, -1.0801, -0.7715, -1.0801, -0.4629, -0.1543, -0.4629, -0.1543, -0.7715, -1.0801,
                                -0.7715, -1.0801, 0.4629, 0.1543, 0.4629, 0.1543, 0.7715, 1.0801, 0.7715, 1.0801, 0.4629, 0.1543, 0.4629, 0.1543, 0.7715, 1.0801,
                                0.7715, 1.0801, -0.4629, -0.1543, -0.4629, -0.1543, -0.7715, -1.0801, -0.7715, -1.0801, -0.4629, -0.1543, -0.4629, -0.1543,
                                -0.7715, -1.0801, -0.7715, -1.0801};

const int _QPSK_demodulation_index[4][2] = {1, 2, 3, 4, 1, 3, 2, 4};
const int _64QAM_demodulation_index[12][32] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                               36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 33, 34, 35, 36,
                                               37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 1, 2, 3,
                                               4, 5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23, 24, 33, 34, 35, 36, 37, 38, 39, 40, 49, 50, 51, 52, 53, 54, 55, 56, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31, 32, 41, 42, 43,
                                               44, 45, 46, 47, 48, 57, 58, 59, 60, 61, 62, 63, 64, 1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20, 25, 26, 27, 28, 33, 34, 35, 36, 41, 42, 43, 44, 49, 50, 51, 52, 57, 58, 59, 60, 5, 6, 7, 8,
                                               13, 14, 15, 16, 21, 22, 23, 24, 29, 30, 31, 32, 37, 38, 39, 40, 45, 46, 47, 48, 53, 54, 55, 56, 61, 62, 63, 64, 1, 2, 5, 6, 9, 10, 13, 14, 17, 18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38,
                                               41, 42, 45, 46, 49, 50, 53, 54, 57, 58, 61, 62, 3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 36, 39, 40, 43, 44, 47, 48, 51, 52, 55, 56, 59, 60, 63, 64, 1, 3, 5, 7, 9,
                                               11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                                               42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64};

double uniform(double a, double b, long int *seed) {
    double t;
    *seed = 2045 * (*seed) + 1;
    *seed = *seed - (*seed / 1048576) * 1048576;
    t = (*seed) / 1048576.0;
    t = a + (b - a) * t;
    return (t);
}

bool InforBitsGen(int *t_bypInforBits, int t_ulInforLength) {
    /*srand(time(0));*/
    long seed = rand();
    double dRand;

    for (unsigned long ulTemp = 0; ulTemp < t_ulInforLength; ulTemp++) {
        dRand = uniform(0, 1, &seed);
        if (dRand > 0.5) {
            t_bypInforBits[ulTemp] = 0;
        } else {
            t_bypInforBits[ulTemp] = 1;
        }
    }

    return 1;
}

void Modulation(int *input, double *output, unsigned long inputlength, unsigned long modulation) {
    if ((inputlength % (int) modulation) != 0) {
        printf("data length is not match modulation type!\n");
        return;
    }
    int outputLength = inputlength / (int) modulation;
    int index;
    switch (modulation) {
        case 2:
            for (int i = 0; i < outputLength; i++) {
                index = input[2 * i] * 2 + input[2 * i + 1];
                output[2 * i] = _QPSK_real[index];
                output[2 * i + 1] = _QPSK_imag[index];
            }
            break;

        case 4:
            for (int i = 0; i < outputLength; i++) {
                index = input[4 * i] * 8 + input[4 * i + 1] * 4 + input[4 * i + 2] * 2 + input[4 * i + 3];
                output[2 * i] = _16QAM_real[index];
                output[2 * i + 1] = _16QAM_imag[index];
            }
            break;

        case 6:
            for (int i = 0; i < outputLength; i++) {
                index = input[6 * i] * 32 + input[6 * i + 1] * 16 + input[6 * i + 2] * 8 +
                        input[6 * i + 3] * 4 + input[6 * i + 4] * 2 + input[6 * i + 5];
                output[2 * i] = _64QAM_real[index];
                output[2 * i + 1] = _64QAM_imag[index];
            }
            break;

        default:
            printf("A wrong modulation type!\n");
    }


}

double LogForModulation(double *distance, int n, unsigned long modulation) {
    double temp = 0;
    switch (modulation) {
        case 2:
            for (int i = 0; i < 2; i++) {
                temp += exp(distance[_QPSK_demodulation_index[n][i] - 1]);
            }
            break;

        case 4:
            for (int i = 0; i < 8; i++) {
                temp += exp(distance[_16QAM_demodulation_index[n][i] - 1]);
            }
            break;

        case 6:
            for (int i = 0; i < 32; i++) {
                temp += exp(distance[_64QAM_demodulation_index[n][i] - 1]);
            }
            break;
        default:
            printf("A wrong modulation type!\n");
    }
    return temp;
}

void DeModulation1(double *input, double *output, unsigned long ulTransLength, unsigned long modulation, double dStan) {
    double sqr_dis[64];
    double max_sqr_dis1;
    double max_sqr_dis2;
    int inputlength = ulTransLength / (int) modulation;
    double variance = 2 * dStan * dStan;

    switch (modulation) {
        case 2:
            for (int i1 = 0; i1 < inputlength; i1++) {
                for (int i2 = 0; i2 < 4; i2++) {
                    sqr_dis[i2] = ((-1) / variance) * ((input[2 * i1] - _QPSK_real[i2]) * (input[2 * i1] - _QPSK_real[i2])
                                                       + (input[2 * i1 + 1] - _QPSK_imag[i2]) * (input[2 * i1 + 1] - _QPSK_imag[i2]));
                }
                max_sqr_dis1 = LogForModulation(sqr_dis, 0, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 1, modulation);
                output[(2 * i1)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 2, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 3, modulation);
                output[(2 * i1 + 1)] = log(max_sqr_dis1) - log(max_sqr_dis2);
            }
            break;

        case 6:
            for (int i1 = 0; i1 < inputlength; i1++) {
                for (int i2 = 0; i2 < 64; i2++) {
                    sqr_dis[i2] = ((-1) / variance) * ((input[2 * i1] - _64QAM_real[i2]) * (input[2 * i1] - _64QAM_real[i2])
                                                       + (input[2 * i1 + 1] - _64QAM_imag[i2]) * (input[2 * i1 + 1] - _64QAM_imag[i2]));
                }
                max_sqr_dis1 = LogForModulation(sqr_dis, 0, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 1, modulation);
                output[(6 * i1)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 2, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 3, modulation);
                output[(6 * i1 + 1)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 4, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 5, modulation);
                output[(6 * i1 + 2)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 6, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 7, modulation);
                output[(6 * i1 + 3)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 8, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 9, modulation);
                output[(6 * i1 + 4)] = log(max_sqr_dis1) - log(max_sqr_dis2);

                max_sqr_dis1 = LogForModulation(sqr_dis, 10, modulation);
                max_sqr_dis2 = LogForModulation(sqr_dis, 11, modulation);
                output[(6 * i1 + 5)] = log(max_sqr_dis1) - log(max_sqr_dis2);
            }
            break;

        default:
            printf("A wrong modulation type!\n");
    }

}

double MaxForModulation(double *distance, int n, unsigned long modulation) {
    double temp = -1e9;
    switch (modulation) {
        case 2:
            for (int i = 0; i < 2; i++) {
                if (distance[_QPSK_demodulation_index[n][i] - 1] > temp) {
                    temp = distance[_QPSK_demodulation_index[n][i] - 1];
                }
            }
            break;

        case 4:
            for (int i = 0; i < 8; i++) {
                if (distance[_16QAM_demodulation_index[n][i] - 1] > temp) {
                    temp = distance[_16QAM_demodulation_index[n][i] - 1];
                }
            }
            break;

        case 6:
            for (int i = 0; i < 32; i++) {
                if (distance[_64QAM_demodulation_index[n][i] - 1] > temp) {
                    temp = distance[_64QAM_demodulation_index[n][i] - 1];
                }
            }
            break;

        default:
            printf("A wrong modulation type!\n");
    }

    return temp;
}

/**
 * 输出似然比 log(p(x=0)/p(x=1))
 * @param input
 * @param output
 * @param ulTransLength
 * @param modulation
 * @param dStan
 */
void DeModulation(double *input, double *output, unsigned long ulTransLength, unsigned long modulation, double dStan) {
    double sqr_dis[64];
    double max_sqr_dis1;
    double max_sqr_dis2;
    int inputlength = ulTransLength / (int) modulation;
    double variance = 2 * dStan * dStan;

    switch (modulation) {
        case 2:
            for (int i1 = 0; i1 < inputlength; i1++) {
                for (int i2 = 0; i2 < 4; i2++) {
                    sqr_dis[i2] = ((-1) / variance) * ((input[2 * i1] - _QPSK_real[i2]) * (input[2 * i1] - _QPSK_real[i2])
                                                       + (input[2 * i1 + 1] - _QPSK_imag[i2]) * (input[2 * i1 + 1] - _QPSK_imag[i2]));
                }
                max_sqr_dis1 = MaxForModulation(sqr_dis, 0, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 1, modulation);
                output[(2 * i1)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 2, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 3, modulation);
                output[(2 * i1 + 1)] = max_sqr_dis1 - max_sqr_dis2;
            }
            break;


        case 6:
            for (int i1 = 0; i1 < inputlength; i1++) {
                for (int i2 = 0; i2 < 64; i2++) {
                    sqr_dis[i2] = ((-1) / variance) * ((input[2 * i1] - _64QAM_real[i2]) * (input[2 * i1] - _64QAM_real[i2])
                                                       + (input[2 * i1 + 1] - _64QAM_imag[i2]) * (input[2 * i1 + 1] - _64QAM_imag[i2]));
                }
                max_sqr_dis1 = MaxForModulation(sqr_dis, 0, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 1, modulation);
                output[(6 * i1)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 2, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 3, modulation);
                output[(6 * i1 + 1)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 4, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 5, modulation);
                output[(6 * i1 + 2)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 6, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 7, modulation);
                output[(6 * i1 + 3)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 8, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 9, modulation);
                output[(6 * i1 + 4)] = max_sqr_dis1 - max_sqr_dis2;

                max_sqr_dis1 = MaxForModulation(sqr_dis, 10, modulation);
                max_sqr_dis2 = MaxForModulation(sqr_dis, 11, modulation);
                output[(6 * i1 + 5)] = max_sqr_dis1 - max_sqr_dis2;
            }
            break;

        default:
            printf("A wrong modulation type!\n");
    }

}

double gauss(double mean, double sigma, long int *s) {
    int i;
    double x, y;
    for (x = 0, i = 0; i < 12; i++) {
        x += uniform(0.0, 1.0, s);
    }
    x = x - 6.0;
    y = mean + x * sigma;
    return (y);

}

void Channel_Gaussian(unsigned long t_ulLength, double *t_dpChannelIn, double *t_dpChannelOut, double t_dMean, double t_dStD) {
    double *dpNoise = new double[t_ulLength];
    long seed = rand();

    for (unsigned long ulTemp = 0; ulTemp < t_ulLength; ulTemp++) {
        dpNoise[ulTemp] = gauss(0, t_dStD, &seed);
        t_dpChannelOut[ulTemp] = t_dpChannelIn[ulTemp] + dpNoise[ulTemp];
    }
    delete[] dpNoise;
}


#endif //INC_5G_LDPC_CHANNAL_H
