//
// Created by ggl on 2018/1/11.
//
#include "LDPC_5G.h"
#include<iostream>
#include "LDPC_helper.h"
#include <windows.h>
#include"Channel.h"
#include"GF.h"

using namespace std;

int main() {
    getTime();
    int codelength = 2000;
    int inflength = 700;
    int type = 2;
    int rv = 3;
    int mod = 2;
    double dSNR = 10.0; //

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
    int *decodeInfBit = new int[inflength];
    cout << "LDPC start: ";
    getTime();

    const int limit_frame = 100;
    long total_Frame = 0;
    double dBER = 0;
    double dBLER = 0;
    long errorBitNum = 0;
    long errorFrameNum = 0;
    long frameErrorBit = 0;
    while (1) {
        total_Frame++;

        InforBitsGen(infbit, inflength);

        ldpc->encoder(infbit, code);

        Modulation(code, mods, codelength, mod);

        Channel_Gaussian(codelength, mods, channelout, 0, dStan);

        DeModulation(channelout, channelLLr, codelength, mod, dStan);


        int iter = ldpc->decode(channelLLr, decodeInfBit, 1, 50);

        frameErrorBit = Diffs(decodeInfBit, infbit, inflength);

        if (frameErrorBit > 0)
            errorFrameNum++;

        errorBitNum += frameErrorBit;

        dBER = (double) errorBitNum * 1.0 / (total_Frame * inflength);
        dBLER = (double) errorFrameNum * 1.0 / total_Frame;
        system("cls");
        cout << "Es/N0 = " << dSNR << "dB" << endl;
        cout << "Frame Number = " << total_Frame << "\titer = " << iter << endl;
        cout << "BLNumber = " << errorFrameNum << endl;
        cout << "BLER = " << dBLER << endl;
        cout << "BER = " << dBER << endl;
        cout << "total bit error =  " << errorBitNum << endl;
        cout << "Frame bit error = " << frameErrorBit << endl;
        cout << "---------------------------------------------------------------------" << endl;
        cout << "Mod type = QPSK   \texpend param = " << ldpc->zLength << endl;
        cout << "information bit number = " << inflength << "\tblock number = " << ldpc->blockNum << endl;
        cout << "code rate = " << ldpc->infLength * 1.0 / ldpc->codeLength << "\tencode block number = " << ldpc->codeLength << endl;

        if (errorFrameNum >= limit_frame) {
            break;
        }
    }
    cout << "LDPC end:   ";
    getTime();

    system("pause");

    return 0;
}