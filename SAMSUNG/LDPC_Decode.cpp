//
// Created by ggl on 2018/1/9.
//


#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;
#define wVarDeg 50
#define wCheDeg 75
#define GAP 0.0000000000001


/********************************文件读取****************************/
void lFileIn(long *t_lpInformation, string &t_strName)
{
    ifstream fFile(t_strName);
    unsigned long ulLength;
    fFile >> ulLength;
    for (unsigned long ulTemp = 0; ulTemp < ulLength; ulTemp++)
    {
        fFile >> t_lpInformation[ulTemp];
    }
    fFile.close();
    fFile.clear();
}

void ulFileIn(unsigned long *t_ulpInformation, string &t_strName)
{
    ifstream fFile(t_strName);
    unsigned long ulLength;
    fFile >> ulLength;
    for (unsigned long ulTemp = 0; ulTemp < ulLength; ulTemp++)
    {
        fFile >> t_ulpInformation[ulTemp];
    }
    fFile.close();
    fFile.clear();
}
/*************************************************************************/

/************************DerateMatch*****************************************/
void LDPC_DeRateMatch(double *input, int inLength, int outLength, double *output)
{
    if (inLength <= outLength) {
        for (int temp1 = 0; temp1<inLength; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = inLength; temp2<outLength; temp2++) {
            *(output + temp2) = 0.0;
        }
    }
    else {
        int numReTrans = inLength - outLength;
        int meanGap = outLength / numReTrans;
        for (int temp1 = 0; temp1<outLength; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = 0; temp2<numReTrans; temp2++) {
            *(output + temp2*meanGap) = *(output + temp2*meanGap) + input[temp2 + outLength];
        }
    }
}
/**********************************************************************************/

/******************************BP_DEC.cpp*********************************/
//ExtractInforBits(t_dpDecoding, ulpReArrangeCol, bypInforbits, t_ulCodeLength, t_ulCheLength);
bool ExtractInforBits(double* t_bypCodeWords, unsigned long* ulpReArrangeCol, double* t_bypInforbits, long t_ulCodeLength, long t_ulCheLength)
{
    unsigned long ulMiddle;
    for (long ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++)
    {
        if (ulpReArrangeCol[ulTemp1] != 0)
        {
            ulMiddle = t_bypCodeWords[ulTemp1];
            //  printf("%f,",t_bypCodeWords[ulTemp1]);
            t_bypCodeWords[ulTemp1] = t_bypCodeWords[ulpReArrangeCol[ulTemp1]];
            t_bypCodeWords[ulpReArrangeCol[ulTemp1]] = ulMiddle;
        }
    }

    for (unsigned long ulTemp1 = t_ulCheLength; ulTemp1 < t_ulCodeLength; ulTemp1++)
    {
        t_bypInforbits[ulTemp1 - t_ulCheLength] = t_bypCodeWords[ulTemp1];
    }
    return 1;
}

void LDPC_RateMatch(double *input, int numInBits, int numOutBits, double *output) {
    if (numInBits >= numOutBits) {
        for (int temp = 0; temp < numOutBits; temp++) {
            *(output + temp) = input[temp];
        }
    } else {
        int numReTrans = numOutBits - numInBits;
        int meanGap = numInBits / numReTrans;
        for (int temp1 = 0; temp1 < numInBits; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = 0; temp2 < numReTrans; temp2++) {
            *(output + numInBits + temp2) = input[temp2 * meanGap];
        }
    }
}


//    Layer_BPDEC(lpVarDis, lpCheDis, dpDECInputLLR, dpDecoding, dpDECOutputLLR, lCodeLength,lCheLength,ulIterMax, idecodeOut, ulpReArrangeCol);
bool Layer_BPDEC(long* t_lpVarDis, long* t_lpCheDis, double* t_dpChannelOut, double* t_dpDecoding, double* t_dpLLR, long t_ulCodeLength, long t_ulCheLength, unsigned long t_ulIterMax, double* bypInforbits, unsigned long* ulpReArrangeCol)
{
    const int layerLength[16] = { 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5 };
    unsigned long ulCount, ulNum, ulCheIndex, ulCheStart;
    unsigned long ulInfLength = t_ulCodeLength - t_ulCheLength;
    long lCheIndex, lVarIndex;
    unsigned int zLength = t_ulCodeLength / wCheDeg;
    double* dpLLRRec = (double *)malloc(t_ulCodeLength*sizeof(double));
    //  t_dpCheSta 校验节点状态，大小为校验节点长，作为先验信息输入，窗译码使用，译整个码直接全设为1
    double* t_dpCheSta = (double *)malloc(t_ulCheLength*sizeof(double));
    double* dpVar2Che = (double *)malloc(t_ulCheLength*wCheDeg*sizeof(double));
    double dMiddle = 0;
    int* bypSyn = (int *)malloc(t_ulCheLength*sizeof(int));


    unsigned char byMiddle;

    //stack<long> layerVar;
    double* dpProductForward = (double *)malloc(t_ulCheLength*wCheDeg*sizeof(double));
    double* dpProductBackward = (double *)malloc(t_ulCheLength*wCheDeg*sizeof(double));
    double* dpMiddle = (double *)malloc(t_ulCodeLength*wVarDeg*sizeof(double));
    double* dpSumForward = (double *)malloc(t_ulCodeLength*wVarDeg*sizeof(double));
    double* dpSumBackward = (double *)malloc(t_ulCodeLength*wVarDeg*sizeof(double));
    //------------------
    for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength*wVarDeg; ulTemp++)
    {
        dpMiddle[ulTemp] = 0.0;
    }
    //--------------------
    for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++)
    {
        dpLLRRec[ulTemp] = t_dpChannelOut[ulTemp];
    }
    for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++)
    {
        t_dpCheSta[ulTemp] = 1.0;
        for (unsigned short wDegTemp = 0; wDegTemp < wCheDeg; wDegTemp++)
        {
            lVarIndex = (long)t_lpCheDis[ulTemp*wCheDeg + wDegTemp];
            dpVar2Che[ulTemp*wCheDeg + wDegTemp] = 0.0;
            if (lVarIndex >= 0)
            {
                dpVar2Che[ulTemp*wCheDeg + wDegTemp] += dpLLRRec[lVarIndex];
            }
        }
    }



    ulCount = 0;
    while (1)
    {
        ulCount++;
        ulCheIndex = 0;
        ulCheStart = 0;
        for (int layer = 0; layer < 16; layer++) {
            unsigned long height = layerLength[layer] * zLength;
            for (unsigned long ulTemp = 0; ulTemp < height; ulTemp++)
            {
                ulCheIndex = ulCheStart + ulTemp;
                dpProductForward[ulCheIndex*wCheDeg] = t_dpCheSta[ulCheIndex];
                dpProductBackward[ulCheIndex*wCheDeg] = 1.0;
                for (unsigned short wDegTemp = 1; wDegTemp < wCheDeg; wDegTemp++)
                {
                    lVarIndex = (long)t_lpCheDis[ulCheIndex*wCheDeg + wDegTemp - 1];
                    //layerVar.push(lVarIndex);
                    if (lVarIndex >= 0)
                    {
                        dpProductForward[ulCheIndex*wCheDeg + wDegTemp] = dpProductForward[ulCheIndex*wCheDeg + wDegTemp - 1] * tanh(0.5 * dpVar2Che[ulCheIndex*wCheDeg + wDegTemp - 1]);
                    }
                    else
                    {
                        dpProductForward[ulCheIndex*wCheDeg + wDegTemp] = dpProductForward[ulCheIndex*wCheDeg + wDegTemp - 1];
                    }
                    lVarIndex = (long)t_lpCheDis[(ulCheIndex + 1)*wCheDeg - wDegTemp];
                    if (lVarIndex >= 0)
                    {
                        dpProductBackward[ulCheIndex*wCheDeg + wDegTemp] = dpProductBackward[ulCheIndex*wCheDeg + wDegTemp - 1] * tanh(0.5 * dpVar2Che[(ulCheIndex + 1)*wCheDeg - wDegTemp]);
                    }
                    else
                    {
                        dpProductBackward[ulCheIndex*wCheDeg + wDegTemp] = dpProductBackward[ulCheIndex*wCheDeg + wDegTemp - 1];
                    }
                }
            }

            for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++)
            {
                for (unsigned short wDegTemp = 0; wDegTemp < wVarDeg; wDegTemp++)
                {
                    lCheIndex = (long)t_lpVarDis[ulTemp*wVarDeg + wDegTemp];
                    if (lCheIndex > ulCheIndex || lCheIndex < ulCheStart) {
                        continue;
                    }
                    if (lCheIndex >= 0)
                    {
                        for (unsigned short wTemp = 0; wTemp < wCheDeg; wTemp++)
                        {
                            if ((long)t_lpCheDis[lCheIndex*wCheDeg + wTemp] == ulTemp)
                            {
                                dMiddle = dpProductForward[wCheDeg*lCheIndex + wTemp] * dpProductBackward[wCheDeg*(lCheIndex + 1) - wTemp - 1];
                                break;
                            }
                        }
                        if ((dMiddle + 1 < GAP) && (dMiddle + 1 >= 0))
                        {
                            dpMiddle[ulTemp*wVarDeg + wDegTemp] = -1e2;
                        }
                        else if ((1 - dMiddle < GAP) && (1 - dMiddle >= 0))
                        {
                            dpMiddle[ulTemp*wVarDeg + wDegTemp] = 1e2;
                        }
                        else
                        {
                            dpMiddle[ulTemp*wVarDeg + wDegTemp] = log((1 + dMiddle) / (1 - dMiddle));
                        }
                    }
                }
            }

            for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++)
            {
                dpSumForward[wVarDeg*ulTemp] = 0.0;

                dpSumBackward[wVarDeg*ulTemp] = 0.0;

                for (unsigned short wDegTemp = 1; wDegTemp < wVarDeg; wDegTemp++)
                {
                    if ((long)t_lpVarDis[ulTemp*wVarDeg + wDegTemp - 1] >= 0)
                    {
                        dpSumForward[ulTemp*wVarDeg + wDegTemp] = dpSumForward[ulTemp*wVarDeg + wDegTemp - 1] + dpMiddle[ulTemp*wVarDeg + wDegTemp - 1];
                    }
                    else
                    {
                        dpSumForward[ulTemp*wVarDeg + wDegTemp] = dpSumForward[ulTemp*wVarDeg + wDegTemp - 1];
                    }
                    if ((long)t_lpVarDis[(ulTemp + 1)*wVarDeg - wDegTemp] >= 0)
                    {
                        dpSumBackward[ulTemp*wVarDeg + wDegTemp] = dpSumBackward[ulTemp*wVarDeg + wDegTemp - 1] + dpMiddle[(ulTemp + 1)*wVarDeg - wDegTemp];
                    }
                    else
                    {
                        dpSumBackward[ulTemp*wVarDeg + wDegTemp] = dpSumBackward[ulTemp*wVarDeg + wDegTemp - 1];
                    }
                }
            }

            for (unsigned short wDegTemp = 0; wDegTemp < wCheDeg; wDegTemp++)
            {
                for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++)
                {
                    lVarIndex = (long)t_lpCheDis[ulTemp*wCheDeg + wDegTemp];
                    if (lVarIndex >= 0)
                    {
                        for (unsigned short wTemp = 0; wTemp < wVarDeg; wTemp++)
                        {
                            lCheIndex = (long)t_lpVarDis[lVarIndex*wVarDeg + wTemp];
                            if (lCheIndex == ulTemp)
                            {
                                dpVar2Che[ulTemp*wCheDeg + wDegTemp] = dpSumForward[lVarIndex*wVarDeg + wTemp] + dpSumBackward[(lVarIndex + 1)*wVarDeg - wTemp - 1] + dpLLRRec[lVarIndex];
                                break;
                            }
                        }
                    }
                }
            }

            for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++)
            {
                t_dpLLR[ulTemp] = dpSumForward[ulTemp*wVarDeg + 1] + dpSumBackward[(ulTemp + 1)*wVarDeg - 1] + dpLLRRec[ulTemp];
            }

            ulCheStart = ulCheStart + height;
        }

        for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++)
        {
            if (t_dpLLR[ulTemp] > 0)
            {
                t_dpDecoding[ulTemp] = 0;
            }
            else
            {
                t_dpDecoding[ulTemp] = 1;
            }
        }

        ulNum = 0;
        for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++)
        {
            if (t_dpCheSta[ulTemp] >= 0)
            {
                byMiddle = 0;
            }
            else
            {
                byMiddle = 1;
            }
            for (unsigned short wTemp = 0; wTemp < wCheDeg; wTemp++)
            {
                lVarIndex = (long)t_lpCheDis[ulTemp * wCheDeg + wTemp];
                if (lVarIndex >= 0)
                {
                    byMiddle += (unsigned char)t_dpDecoding[lVarIndex];
                }
            }
            bypSyn[ulTemp] = byMiddle % 2;
            if (bypSyn[ulTemp] != 0)
            {
                ulNum++;
            }
        }

        ExtractInforBits(t_dpLLR, ulpReArrangeCol, bypInforbits, t_ulCodeLength, t_ulCheLength);

        if (ulCount == t_ulIterMax)
        {
            break;
        }

    }

    free(bypSyn);
    free(dpLLRRec);
    free(t_dpCheSta);
    free(dpVar2Che);
    free(dpProductForward);
    free(dpProductBackward);
    free(dpSumForward);
    free(dpSumBackward);
    free(dpMiddle);

    if (ulNum == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
/**************************************************************************/
int LDPC_decode(double *dpLLR, int iCodeLength, int iInfLength, double *bypInforbits, double *OutputLLRWithoutChannel){
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*输入参数：
    *参数1：dpLLR----经过解调的编码序列
    *参数2：dpRate---编码速率
    */
    static double sRate = 0;
    static long llen = 0;
    /********************静态内存分配*********************/


    static double* dpDECInputLLR = NULL;
    static long *lpCheDis = NULL;
    static long *lpVarDis = NULL;
    static double* dpDecoding = NULL;
    static double* dpDECOutputLLR = NULL;
    static unsigned long* ulpReArrangeCol = NULL;
    /******************************************************/
//    /*****************检查输入参数********************************/
//    if (nrhs == 1) {
//        char *cpstr = mxArrayToString(prhs[0]);
//        if (strcmp(cpstr, "release") == 0) {
//            if (dpDECInputLLR != NULL) {
//                delete[] dpDECInputLLR;
//                dpDECInputLLR = NULL;
//            }
//            if (lpCheDis != NULL) {
//                delete[] lpCheDis;
//                lpCheDis = NULL;
//            }
//            if (lpVarDis != NULL) {
//                delete[] lpVarDis;
//                lpVarDis = NULL;
//            }
//            if (dpDecoding != NULL) {
//                delete[] dpDecoding;
//                dpDecoding = NULL;
//            }
//            if (dpDECOutputLLR != NULL) {
//                delete[] dpDECOutputLLR;
//                dpDECOutputLLR = NULL;
//            }
//            if (ulpReArrangeCol != NULL) {
//                delete[] ulpReArrangeCol;
//                ulpReArrangeCol = NULL;
//            }
//            sRate = 0;
//            llen = 0;
//        }
//        else if (strcmp(cpstr, "help") == 0) {
//            mexPrintf("    本模块输入为需要译码的软信息和码率，测试时，直接输入解调器输出的软信息\n    模块首先进行解速率适配，然后再用多层的BP进行迭代译码，迭代次数为15次，没有CRC校验\n    特别提醒：本模块需要和LDPC_EnCode.mexw64配合使用，当单独使用本模块时，要先用编码的模块生成矩阵，然后才能使用译码模块\n");
//            mexPrintf("//////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
//            mexPrintf("使用方法：\n   *方案1：\n      *  参数1：软解调信息（double）\n      *  参数2：编码速率(速率匹配后)和编码模块的速率一样（double）\n      *  输出：译码结果int32_t类型\n\n");
//            mexPrintf("   *方案2：\n      * 只输入'release',用于所有的编码完成后，内存的释放,在译码过程中使用\n\n   *注意事项：当码率改变或改变码长时，不需要做任何调整\n");
//            mexPrintf("//////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
//        }
//        else {
//            mexPrintf("输入的参数不正确！\n");
//        }
//        return;
//    }
//    else if (nrhs != 2) {
//        mexPrintf("输入的参数数目不正确！\n");
//        return;
//    }
//    /*************************************************************/



    double dCodeRate = 1.0 / 3; // 编码码率  没有凿孔之前
    /***************获得输入参数*************************/
    //double*dpLLR = mxGetPr(prhs[0]);
    //long iCodeLength = (long)mxGetN(prhs[0]);  //输入的编码序列的长度
    double dpRate =iInfLength*1.0/iCodeLength;   //编码的码率

    /***************************************************/
    /*******************参数初始化**********************************/
    unsigned long ulIterMax = 15;
    //long iInfLength = iCodeLength*(*dpRate);
    long lCodeLength = iInfLength / dCodeRate;
    long lCheLength = lCodeLength - iInfLength;
    /*******************************************************/





    /**********************矩阵读取***********************/
    if (sRate != dpRate || llen != iCodeLength) {
        sRate = dpRate;
        llen = iCodeLength;
        if (dpDECInputLLR == NULL) {
            dpDECInputLLR = new double[lCodeLength];
            lpCheDis = new long[wCheDeg*lCheLength];
            lpVarDis = new long[wVarDeg*lCodeLength];
            dpDecoding = new double[lCodeLength];
            dpDECOutputLLR = new double[lCodeLength];
            ulpReArrangeCol = new unsigned long[lCheLength];
        }
        else {
            delete[] dpDECInputLLR;
            delete[] lpCheDis;
            delete[] lpVarDis;
            delete[] dpDecoding;
            delete[] dpDECOutputLLR;
            delete[] ulpReArrangeCol;

            dpDECInputLLR = new double[lCodeLength];
            lpCheDis = new long[wCheDeg*lCheLength];
            lpVarDis = new long[wVarDeg*lCodeLength];
            dpDecoding = new double[lCodeLength];
            dpDECOutputLLR = new double[lCodeLength];
            ulpReArrangeCol = new unsigned long[lCheLength];
        }

        string Filename;
        Filename = ".\\data\\LDPC_CheDis.txt";
        lFileIn(lpCheDis, Filename);
        Filename = ".\\data\\LDPC_VarDis.txt";
        lFileIn(lpVarDis, Filename);
        Filename = ".\\data\\ReArrangeCol.txt";
        ulFileIn(ulpReArrangeCol, Filename);
    }
    /******************************************************/
    /***************************输出定义********************/
//    plhs[0] = mxCreateNumericMatrix(1, iInfLength, mxDOUBLE_CLASS, mxREAL);
//    double* idecodeOut = (double*)mxGetData(plhs[0]);

    /*******************************************************/
    LDPC_DeRateMatch(dpLLR, iCodeLength, lCodeLength, dpDECInputLLR);
    Layer_BPDEC(lpVarDis, lpCheDis, dpDECInputLLR, dpDecoding, dpDECOutputLLR, lCodeLength, lCheLength, ulIterMax, bypInforbits, ulpReArrangeCol);

    for (int i = 0; i < lCodeLength; i++) {
        dpDECOutputLLR[i] -= dpDECInputLLR[i];
    }
    LDPC_RateMatch(dpDECOutputLLR, lCodeLength, iCodeLength, OutputLLRWithoutChannel);
}