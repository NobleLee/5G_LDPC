//
// Created by ggl on 2018/1/11.
//
#include "LDPC_5G.h"
#include<iostream>
#include"Channel.h"
#include"GF.h"

using namespace std;


int main() {
    getTime();
    int codelength = 4000;
    int inflength = 700;
    int type = 2;
    int rv = 0;
    int mod = 2;
    double dSNR = 20.0; // ϵͳSNR,ΪEs/N0

    double dLinearSNR = pow(10.0, dSNR / 10);
    double dStan;
    if (mod != 1) {
        dStan = sqrt(1 / (2 * dLinearSNR));
    } else {
        dStan = sqrt(1 / (dLinearSNR));
    }

    LDPC_5G *ldpc = new LDPC_5G(inflength, codelength, type, rv, mod);
    int *infbit = new int[inflength];
    int *code = new int[codelength];
    double *mods = new double[codelength];
    double *channelout = new double[codelength];
    double *channelLLr = new double[codelength];
    int *decodebit = new int[codelength];
    int count = 0;
    cout << "LDPC start: ";
    getTime();
    int errorNum = 0;
    while (count < 1000) {
        InforBitsGen(infbit, inflength);

        ldpc->encoder(infbit, code);

        Modulation(code, mods, codelength, mod);

        Channel_Gaussian(codelength, mods, channelout, 0, dStan);

        DeModulation(channelout, channelLLr, codelength, mod, dSNR);
        /*for (int i = 0; i < codelength; i++) {
        channelLLr[i] = code[i];
        }*/
        int iter = ldpc->decode(channelLLr, decodebit, 1, 50);


        for (int i = 0; i < inflength; i++) {
            if (infbit[i] != decodebit[i]) {
                cout << i << "  ";
                errorNum++;
            }
        }
        cout << "block: " << count++ << "error: " << errorNum << "\titer: " << iter << "\t\n";
        system("cls");

    }
    cout << "LDPC end:   ";
    getTime();

    cout << count << endl;
    system("pause");
    // main1();

    return 0;
}