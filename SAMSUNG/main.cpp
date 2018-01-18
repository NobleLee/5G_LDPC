#include "LDPC.h"

void ulFileOut(unsigned long *t_ulpInformation, unsigned long t_ulLength, string &t_strName)
{
    ofstream fFile(t_strName);
    fFile << t_ulLength << " ";
    for (unsigned long ulTemp = 0; ulTemp < t_ulLength; ulTemp++)
    {
        fFile << t_ulpInformation[ulTemp] << " ";
    }
    fFile.close();
    fFile.clear();
}

void wFileOut(unsigned short *t_wpInformation, unsigned long t_ulLength, string &t_strName)
{
    ofstream fFile(t_strName);
    fFile << t_ulLength << " ";
    for (unsigned long ulTemp = 0; ulTemp < t_ulLength; ulTemp++)
    {
        fFile << t_wpInformation[ulTemp] << " ";
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

void lFileIn(long *t_lpInformation, string t_strName)
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

void wFileIn(unsigned short *t_wpInformation, string t_strName)
{
    ifstream fFile(t_strName);
    unsigned long ulLength;
    fFile >> ulLength;
    for (unsigned long ulTemp = 0; ulTemp < ulLength; ulTemp++)
    {
        fFile >> t_wpInformation[ulTemp];
    }
    fFile.close();
    fFile.clear();
}

const int dummy_LDPC = 10000;
const int perm[] = { 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31 };

void LDPC_RateMatch(unsigned short *input, int numInBits, int numOutBits, int *output)
{
    if (numInBits >= numOutBits) {
        for (int temp = 0; temp < numOutBits; temp++) {
            *(output + temp) = input[temp];
        }
    }
    else {
        int numReTrans = numOutBits - numInBits;
        int meanGap = numInBits / numReTrans;
        for (int temp1 = 0; temp1 < numInBits; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = 0; temp2 < numReTrans; temp2++) {
            *(output + numInBits + temp2) = input[temp2*meanGap];
        }
    }
}

void LDPC_RateMatch(double *input, int numInBits, int numOutBits, double *output) {
    if (numInBits >= numOutBits) {
        for (int temp = 0; temp < numOutBits; temp++) {
            *(output + temp) = input[temp];
        }
    }
    else {
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

//  根据LDPC码的矩阵，得到其按行顺序和列顺序排列的1的位置
//  输入：
//  t_bypMatrix  --->  LDPC码矩阵(该矩阵储存方式如MATLAB，按列读取)
//  t_ulRow,t_ulCol  --->  矩阵的行数和列数
//  t_VarDeg,t_CheDeg  --->  矩阵的变量节点度和校验节点度
//  得到的结果被保存至txt文件中

bool Position_Decetion(unsigned short* t_bypMatrix, unsigned long t_ulRow, unsigned long t_ulCol, unsigned short t_VarDeg, unsigned short t_CheDeg)
{
    unsigned long ulTemp1, ulTemp2, ulTemp3;
    long* ulpVarDis = new long[t_VarDeg*t_ulCol];
    long* ulpCheDis = new long[t_CheDeg*t_ulRow];

    for (ulTemp1 = 0; ulTemp1 < t_VarDeg*t_ulCol; ulTemp1++)
    {
        ulpVarDis[ulTemp1] = -1;
    }
    for (ulTemp1 = 0; ulTemp1 < t_CheDeg*t_ulRow; ulTemp1++)
    {
        ulpCheDis[ulTemp1] = -1;
    }

    //  按列查找到所有1的位置
    for (ulTemp1 = 0; ulTemp1 < t_ulCol; ulTemp1++)
    {
        ulTemp3 = 0;
        for (ulTemp2 = 0; ulTemp2 < t_ulRow; ulTemp2++)
        {
            if (t_bypMatrix[ulTemp1*t_ulRow + ulTemp2] == 1)
            {
                if (ulTemp3 < t_VarDeg)
                {
                    ulpVarDis[ulTemp1*t_VarDeg + ulTemp3] = ulTemp2;
                    ulTemp3++;
                }
                else
                {
                    return 0;
                }

            }
        }
    }

    //  按行查到所有1的位置
    for (ulTemp1 = 0; ulTemp1 < t_ulRow; ulTemp1++)
    {
        ulTemp3 = 0;
        for (ulTemp2 = 0; ulTemp2 < t_ulCol; ulTemp2++)
        {
            if (t_bypMatrix[ulTemp1 + ulTemp2*t_ulRow] == 1)
            {
                if (ulTemp3 < t_CheDeg)
                {
                    ulpCheDis[ulTemp1*t_CheDeg + ulTemp3] = ulTemp2;
                    ulTemp3++;
                }
                else
                {
                    return 0;
                }

            }
        }
    }

    //  将结果分别保存至txt文件中
    ofstream FileLDPCVarDis(".\\LDPC_matrix\\LDPC_VarDis.txt");
    ofstream FileLDPCCheDis(".\\LDPC_matrix\\LDPC_CheDis.txt");

    FileLDPCVarDis << t_VarDeg*t_ulCol << endl;
    FileLDPCCheDis << t_CheDeg*t_ulRow << endl;
    for (ulTemp1 = 0; ulTemp1 < t_VarDeg*t_ulCol; ulTemp1++)
    {
        FileLDPCVarDis << ulpVarDis[ulTemp1] << endl;
    }
    for (ulTemp1 = 0; ulTemp1 < t_CheDeg*t_ulRow; ulTemp1++)
    {
        FileLDPCCheDis << ulpCheDis[ulTemp1] << endl;
    }

    FileLDPCVarDis.close();
    FileLDPCVarDis.clear();
    FileLDPCCheDis.close();
    FileLDPCCheDis.clear();

    delete[] ulpVarDis;
    delete[] ulpCheDis;
    return 1;
}

bool Gaussian_Elimination(unsigned short* t_LDPC_Matrix, unsigned short* t_P_Matrix, unsigned long* t_ulpReArrangeCol, unsigned long t_ulCheLength, unsigned long t_ulCodeLength)
{
    unsigned long ulTemp1, ulTemp2, ulTemp3, ulMiddle1, ulMiddle2;

    for (ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++)
    {
        ulMiddle1 = ulTemp1;
        while (ulMiddle1 < t_ulCodeLength && t_LDPC_Matrix[ulTemp1 + ulMiddle1*t_ulCheLength] == 0)
        {
            for (ulTemp2 = ulTemp1 + 1; ulTemp2 < t_ulCheLength; ulTemp2++)
            {
                if (t_LDPC_Matrix[ulTemp2 + ulMiddle1*t_ulCheLength] != 0)
                {
                    break;
                }
            }
            ulMiddle1++;
        }

        if (ulMiddle1 >= t_ulCodeLength)
        {
            return 0;
        }

        if (ulMiddle1 != ulTemp1)
        {
            t_ulpReArrangeCol[ulTemp1] = ulMiddle1;
            for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++)
            {
                ulMiddle2 = t_LDPC_Matrix[ulTemp2 + ulTemp1*t_ulCheLength];
                t_LDPC_Matrix[ulTemp2 + ulTemp1*t_ulCheLength] = t_LDPC_Matrix[ulTemp2 + ulMiddle1*t_ulCheLength];
                t_LDPC_Matrix[ulTemp2 + ulMiddle1*t_ulCheLength] = ulMiddle2;
            }
        }

        if (t_LDPC_Matrix[ulTemp1*t_ulCheLength + ulTemp1] == 0)
        {
            for (ulTemp2 = ulTemp1 + 1; ulTemp2 < t_ulCheLength; ulTemp2++)
            {
                if (t_LDPC_Matrix[ulTemp2 + ulTemp1*t_ulCheLength] != 0)
                {
                    ulMiddle1 = ulTemp2;
                    for (ulTemp3 = ulTemp1; ulTemp3 < t_ulCodeLength; ulTemp3++)
                    {
                        t_LDPC_Matrix[ulTemp1 + ulTemp3*t_ulCheLength] = (t_LDPC_Matrix[ulMiddle1 + ulTemp3*t_ulCheLength] + t_LDPC_Matrix[ulTemp1 + ulTemp3*t_ulCheLength]) % 2;
                    }
                }
            }
        }

        for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++)
        {
            if (t_LDPC_Matrix[ulTemp1*t_ulCheLength + ulTemp2] != 0)
            {
                if (ulTemp2 != ulTemp1)
                {
                    for (ulTemp3 = ulTemp1; ulTemp3 < t_ulCodeLength; ulTemp3++)
                    {
                        t_LDPC_Matrix[ulTemp3*t_ulCheLength + ulTemp2] = (t_LDPC_Matrix[ulTemp3*t_ulCheLength + ulTemp2] + t_LDPC_Matrix[ulTemp3*t_ulCheLength + ulTemp1]) % 2;
                    }
                }
            }
        }
    }

    for (ulTemp1 = t_ulCheLength; ulTemp1 < t_ulCodeLength; ulTemp1++)
    {
        for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++)
        {
            t_P_Matrix[(ulTemp1 - t_ulCheLength)*t_ulCheLength + ulTemp2] = t_LDPC_Matrix[ulTemp1*t_ulCheLength + ulTemp2];
        }
    }

    return 1;
}

bool LDPC_Matrix_Generator(unsigned long ulInfLength, unsigned long ulCheLength)
{
    unsigned long ultemp1, ultemp2, ultemp3, temp;
    unsigned short wInfDeg = wCheDeg - wVarDeg;
    string Filename;
    long zLength = ulInfLength / wInfDeg;
    unsigned long ulCodeLength = ulInfLength + ulCheLength;

    long* BaseMatrix = new long[wCheDeg*wVarDeg];
    Filename = ".\\LDPC_matrix\\BaseMatrix.txt";
    lFileIn(BaseMatrix, Filename);

    for (ultemp1 = 0; ultemp1 < wVarDeg*wCheDeg; ultemp1++) {
        if (BaseMatrix[ultemp1] < zLength) {
            continue;
        }
        BaseMatrix[ultemp1] = BaseMatrix[ultemp1] % zLength;
    }

    unsigned short* LDPC_Matrix = new unsigned short[ulCheLength*ulCodeLength]();
    unsigned short* P_Matrix = new unsigned short[ulCheLength*ulInfLength];
    unsigned long* ulpReArrangeCol = new unsigned long[ulCheLength]();

    for (ultemp1 = 0; ultemp1 < wCheDeg; ultemp1++) {
        for (ultemp2 = 0; ultemp2 < wVarDeg; ultemp2++) {
            if (BaseMatrix[ultemp2*wCheDeg + ultemp1] == -1) {
                continue;
            }
            for (ultemp3 = 0; ultemp3 < zLength; ultemp3++) {
                temp = (BaseMatrix[ultemp2*wCheDeg + ultemp1] + ultemp3) % zLength;
                LDPC_Matrix[(ultemp1*zLength + ultemp3)*wVarDeg*zLength + ultemp2*zLength + temp] = 1;
            }
        }
    }

    string filename;
    filename = ".\\LDPC_matrix\\LDPC_Matrix.txt";
    wFileOut(LDPC_Matrix, ulCheLength*ulCodeLength, filename);

    Position_Decetion(LDPC_Matrix, ulCheLength, ulCodeLength, wVarDeg, wCheDeg);
    Gaussian_Elimination(LDPC_Matrix, P_Matrix, ulpReArrangeCol, ulCheLength, ulCodeLength);

    filename = ".\\LDPC_matrix\\PMatrix.txt";
    wFileOut(P_Matrix, ulCheLength*ulInfLength, filename);
    filename = ".\\LDPC_matrix\\ReArrangeCol.txt";
    ulFileOut(ulpReArrangeCol, ulCheLength, filename);

    delete[] BaseMatrix;
    delete[] LDPC_Matrix;
    delete[] P_Matrix;
    delete[] ulpReArrangeCol;
    return 1;
}

bool LDPC_Fast_Encoder(int * t_bypInforbits, unsigned short* t_bypCodeWords, unsigned short* P_Matrix, unsigned long* ulpReArrangeCol, unsigned long t_ulCodeLength, unsigned long t_ulCheLength)
{
    unsigned long ulInfBitLength = t_ulCodeLength - t_ulCheLength;
    unsigned long ulTemp1, ulTemp2, ulMiddle;

    unsigned short* ulCheckBits = new unsigned short[t_ulCheLength];

    for (ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++)
    {
        ulMiddle = 0;
        for (ulTemp2 = 0; ulTemp2 < ulInfBitLength; ulTemp2++)
        {
            ulMiddle += P_Matrix[ulTemp1 + ulTemp2*t_ulCheLength] * t_bypInforbits[ulTemp2];
        }
        ulCheckBits[ulTemp1] = ulMiddle % 2;
    }

    for (ulTemp1 = 0; ulTemp1 < t_ulCodeLength; ulTemp1++)
    {
        if (ulTemp1 < t_ulCheLength)
        {
            t_bypCodeWords[ulTemp1] = ulCheckBits[ulTemp1];
        }
        else
        {
            t_bypCodeWords[ulTemp1] = t_bypInforbits[ulTemp1 - t_ulCheLength];
        }
    }

    for (long ulTemp1 = t_ulCheLength - 1; ulTemp1 >= 0; ulTemp1--)
    {
        if (ulpReArrangeCol[ulTemp1] != 0)
        {
            ulMiddle = t_bypCodeWords[ulTemp1];
            t_bypCodeWords[ulTemp1] = t_bypCodeWords[ulpReArrangeCol[ulTemp1]];
            t_bypCodeWords[ulpReArrangeCol[ulTemp1]] = ulMiddle;
        }
    }

    delete[] ulCheckBits;

    return 1;
}

void LDPC_DeRateMatch(double *input, int inLength, int outLength, double *output)
{
    if (inLength <= outLength) {
        for (int temp1 = 0; temp1 < inLength; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = inLength; temp2 < outLength; temp2++) {
            *(output + temp2) = 0.0;
        }
    }
    else {
        int numReTrans = inLength - outLength;
        int meanGap = outLength / numReTrans;
        for (int temp1 = 0; temp1 < outLength; temp1++) {
            *(output + temp1) = input[temp1];
        }
        for (int temp2 = 0; temp2 < numReTrans; temp2++) {
            *(output + temp2*meanGap) = *(output + temp2*meanGap) + input[temp2 + outLength];
        }
    }
}


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

bool Layer_BPDEC(long* t_lpVarDis, long* t_lpCheDis, double* t_dpChannelOut, double* t_dpDecoding, double* t_dpLLR, long t_ulCodeLength, long t_ulCheLength, unsigned long t_ulIterMax, double* bypInforbits, unsigned long* ulpReArrangeCol)
{
    const int layerLength[16] = { 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5 };
    unsigned long ulCount, ulNum, ulCheIndex, ulCheStart;
    unsigned long ulInfLength = t_ulCodeLength - t_ulCheLength;
    long lCheIndex, lVarIndex;
    unsigned int zLength = t_ulCodeLength / wCheDeg;
    double* dpLLRRec = (double *)malloc(t_ulCodeLength * sizeof(double));

    //  t_dpCheSta 校验节点状态，大小为校验节点长，作为先验信息输入，窗译码使用，译整个码直接全设为1
    double* t_dpCheSta = (double *)malloc(t_ulCheLength * sizeof(double));
    double* dpVar2Che = (double *)malloc(t_ulCheLength*wCheDeg * sizeof(double));
    double dMiddle = 0;
    int* bypSyn = (int *)malloc(t_ulCheLength * sizeof(int));


    unsigned char byMiddle;

    //stack<long> layerVar;
    double* dpProductForward = (double *)malloc(t_ulCheLength*wCheDeg * sizeof(double));
    double* dpProductBackward = (double *)malloc(t_ulCheLength*wCheDeg * sizeof(double));
    double* dpMiddle = (double *)malloc(t_ulCodeLength*wVarDeg * sizeof(double));
    double* dpSumForward = (double *)malloc(t_ulCodeLength*wVarDeg * sizeof(double));
    double* dpSumBackward = (double *)malloc(t_ulCodeLength*wVarDeg * sizeof(double));
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


void LDPC_encode(int *ipInforBitsNoCRC, int* wTransBits, int iInfLength, int ulTransLength)
{

    /***********************内存不释放变量定义区******************/
    static unsigned short* wEncodedBits = NULL;
    //static unsigned short* wTransBits=NULL;
    static unsigned long* ulpReArrangeCol = NULL;
    static unsigned short* P_Matrix = NULL;

    static double sdRate = 0;
    static int ilen = 0;
    /***********************************************************/



    //检验输入参数
    static int iflag = 0;


    /**************************************************************/
    /*初始化：
    * dpInforBitsNoCRC：输入的信息bit
    * iInfLength：信息bit长度
    * dpRate：速率匹配后的码率
    * dcode_out:编码之后的bit
    */;
    double dpRate = iInfLength*1.0 / ulTransLength;                        //速率匹配后的码率

    /***********************************************************/
    double dCodeRate = 1.0 / 3; // 编码码率  没有凿孔之前
    unsigned long ulCodeLength = iInfLength / dCodeRate;
    unsigned long ulCheLength = ulCodeLength - iInfLength;


    /**********************生成矩阵******************************/
    if (sdRate != dpRate || ilen != iInfLength) {
        iflag = 0;
        ilen = iInfLength;
        sdRate = dpRate;
    }
    if (iflag == 0) {
        // mexPrintf("正在生成矩阵......\n");
        LDPC_Matrix_Generator(iInfLength, ulCheLength);
    }
    /******************************************************************/

    if (iflag == 0) {
        // mexPrintf("正在开辟内存空间.........\n");
        if (wEncodedBits != NULL) {
            delete[] wEncodedBits;
            wEncodedBits = NULL;
        }
        //  if(wTransBits!=NULL) delete[] wTransBits;
        if (ulpReArrangeCol != NULL) {
            delete[] ulpReArrangeCol;
            ulpReArrangeCol = NULL;
        }
        if (P_Matrix != NULL) {
            delete[] P_Matrix;
            P_Matrix = NULL;
        }

        wEncodedBits = new unsigned short[ulCodeLength];;
        ulpReArrangeCol = new unsigned long[ulCheLength];
        P_Matrix = new unsigned short[ulCheLength*iInfLength];
        // mexPrintf("正在读取矩阵.........\n");
        string Filename;
        Filename = ".\\LDPC_matrix\\PMatrix.txt";
        wFileIn(P_Matrix, Filename);
        Filename = ".\\LDPC_matrix\\ReArrangeCol.txt";
        ulFileIn(ulpReArrangeCol, Filename);
    }
    /**********************************************************/

    // mexPrintf("正在编码......\n");
    LDPC_Fast_Encoder(ipInforBitsNoCRC, wEncodedBits, P_Matrix, ulpReArrangeCol, ulCodeLength, ulCheLength);

    //RateMatch(wEncodedBits, ulCodeLength, ulTransLength, wTransBits);
    //mexPrintf("正在速率适配......\n");
    LDPC_RateMatch(wEncodedBits, ulCodeLength, ulTransLength, wTransBits);

    iflag = 1;
}

int LDPC_decode(double *dpLLR, int iCodeLength, int iInfLength, double *bypInforbits, double *OutputLLRWithoutChannel) {

    static double sRate = 0;
    static long llen = 0;

    /********************静态内存分配*********************/


    static double* dpDECInputLLR = NULL;
    static long *lpCheDis = NULL;
    static long *lpVarDis = NULL;
    static double* dpDecoding = NULL;
    static double* dpDECOutputLLR = NULL;
    static unsigned long* ulpReArrangeCol = NULL;

    double dCodeRate = 1.0 / 3; // 编码码率  没有凿孔之前

    double dpRate = iInfLength*1.0 / iCodeLength;   //编码的码率

    /***************************************************/
    /*******************参数初始化**********************************/
    unsigned long ulIterMax = 15;
    long lCodeLength = iInfLength / dCodeRate;
    long lCheLength = lCodeLength - iInfLength;

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

        LDPC_Matrix_Generator(iInfLength, lCheLength);

        string Filename;
        Filename = ".\\LDPC_matrix\\LDPC_CheDis.txt";
        lFileIn(lpCheDis, Filename);
        Filename = ".\\LDPC_matrix\\LDPC_VarDis.txt";
        lFileIn(lpVarDis, Filename);
        Filename = ".\\LDPC_matrix\\ReArrangeCol.txt";
        ulFileIn(ulpReArrangeCol, Filename);
    }

    LDPC_DeRateMatch(dpLLR, iCodeLength, lCodeLength, dpDECInputLLR);
    Layer_BPDEC(lpVarDis, lpCheDis, dpDECInputLLR, dpDecoding, dpDECOutputLLR, lCodeLength, lCheLength, ulIterMax, bypInforbits, ulpReArrangeCol);

    for (int i = 0; i < lCodeLength; i++) {
        dpDECOutputLLR[i] -= dpDECInputLLR[i];
    }

    LDPC_RateMatch(dpDECOutputLLR, lCodeLength, iCodeLength, OutputLLRWithoutChannel);

    return ulIterMax;
}




// LDPC MIMO
// Turbo MIMO使用
void ittc_LDPC_code_mimo_MQAM_transmit_PIC(int nt, ittc_vec *mapper, ittc_ivec *rate_match_bits, ittc_arr_cmat *tx_signal)
{
    /* Asumption: the average energy of the transmitted symbols is normalized to one.
    Parameter List:
    > nt     : number of transmitting antennas.
    > mapper : 1-dimensional PAM mapper with modulation order m/2.
    ...
    */

    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d);
    int m = p * 2;
    int N = rate_match_bits->size / nt / m; // here, N is the code length of the sub-codes

    static ittc_ivec *tmp_mod = NULL;
    if (tmp_mod == NULL || tmp_mod->size != m / 2)
    {
        if (tmp_mod != NULL)
            ittc_ivec_free(tmp_mod);
        tmp_mod = ittc_ivec_alloc(m / 2);
    }


    for (int i = 0; i < N; i++)  // 时隙
    {
        for (int j = 0; j < nt; j++)  // 天线序号
        {
            int sym_idx;
            ittc_complex c;
            for (int k = 0; k < p; k++)  // real part
            {
                tmp_mod->dat[k] = ittc_ivec_get(rate_match_bits, i*nt*m + j*m + k);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.real = ittc_vec_get(mapper, sym_idx);

            for (int k = 0; k < p; k++)  // imag part
            {
                tmp_mod->dat[k] = ittc_ivec_get(rate_match_bits, i*nt*m + j*m + k + p);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.imag = ittc_vec_get(mapper, sym_idx);
            ittc_cmat_set(tx_signal->arr[i], j, 0, c);
        }
    }
}

void ittc_LDPC_code_mimo_MQAM_detector_demodulator_LMMSE(ittc_arr_cmat *rx_signal, ittc_vec *mapper, ittc_arr_cmat *arrH, double nvar,
                                                         ittc_vec *bit_LLR)
{
    // 线性MMSE
    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d); // every symbol is transmitted over p parallel bit-channels.
    int m = p * 2;
    int M = M1d*M1d;

    int nr = arrH->arr[0]->rows; // number of receive antennas
    int nt = arrH->arr[0]->cols; // number of transmit antennas
    int N = arrH->size; // 时隙总量

    /* Step 1: ready the space*/
    //arrH_cha  H + eye(nt) -> H_cha (H_change)
    static ittc_cmat *H_cha = NULL;
    if (H_cha == NULL || H_cha->cols != nt || H_cha->rows != nr)
    {
        if (H_cha != NULL)
            ittc_cmat_free(H_cha);
        H_cha = ittc_cmat_alloc(nr, nt);
    }

    static ittc_cmat *H_mmse = NULL;
    if (H_mmse == NULL || H_mmse->cols != nr || H_mmse->rows != nt)
    {
        if (H_mmse != NULL)
            ittc_cmat_free(H_mmse);
        H_mmse = ittc_cmat_alloc(nt, nr);
    }

    static ittc_cmat *wh = NULL;
    if (wh == NULL || wh->cols != nt || wh->rows != nt)
    {
        if (wh != NULL)
            ittc_cmat_free(wh);
        wh = ittc_cmat_alloc(nt, nt);
    }

    static ittc_cmat *wy = NULL;
    if (wy == NULL || wy->rows != nt || wy->cols != 1)
    {
        if (wy != NULL)
            ittc_cmat_free(wy);
        wy = ittc_cmat_alloc(nt, 1);
    }
    ittc_cmat_set_zero(wy);

    static ittc_vec *signal_prob = NULL;
    if (signal_prob == NULL || signal_prob->size != M1d)
    {
        if (signal_prob != NULL)
            ittc_vec_free(signal_prob);
        signal_prob = ittc_vec_alloc(M1d);
    }
    ittc_vec_set_zero(signal_prob);  // Signal Prob


    /* Step 2: MIMO Detection and Demodulation */
    for (int timeslot = 0; timeslot < N; timeslot++)
    {
        ittc_cmat_memcpy(arrH->arr[timeslot], H_cha);
        // MMSE滤波矩阵计算
        ittc_cmat_calmmse(H_cha, nvar, 1.0, H_mmse);  // MMSE矩阵
        ittc_cmat_mult(H_mmse, H_cha, wh);
        ittc_cmat_mult(H_mmse, rx_signal->arr[timeslot], wy);

        // 逐天线计算
        for (int i = 0; i < nt; i++)
        {
            double tmp1 = 0;  // 干扰
            for (int k = 0; k < nt; k++)
            {
                if (k != i)
                {
                    tmp1 += ittc_complex_abs2(ittc_cmat_get(wh, i, k));
                }
            }

            double tmp2 = 0;  // 噪声
            for (int k = 0; k < nr; k++)
            {
                tmp2 += ittc_complex_abs2(ittc_cmat_get(H_mmse, i, k));
            }
            tmp2 = tmp2*nvar;

            double tmp3 = ittc_complex_abs2(ittc_cmat_get(wh, i, i));  // 有用信号
            double SINR = tmp3 / (tmp1 + tmp2);
            ittc_complex signal = ittc_div_cc(ittc_cmat_get(wy, i, 0), ittc_cmat_get(wh, i, i));
            double Eq_nvar = 1.0 / SINR;

            // Demodulation
            // Real Part
            for (int k = 0; k < M1d; k++)
            {
                signal_prob->dat[k] = -(pow((signal.real - mapper->dat[k]), 2)) / Eq_nvar;
            }

            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k] = p0 - p1;  // LLR Real Part
            }

            // Imag Part
            for (int k = 0; k < M1d; k++)
            {
                signal_prob->dat[k] = -(pow((signal.imag - mapper->dat[k]), 2)) / Eq_nvar;
            }

            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k + p] = p0 - p1;  // LLR Real Part
            }
        }
    }
}

void ittc_LDPC_code_mimo_MQAM_detector_demodulator_MMSE_PIC_Simple(ittc_arr_cmat *rx_signal, ittc_vec *decoded_LLR, ittc_vec *mapper, ittc_arr_cmat *arrH, double nvar,
                                                                   ittc_vec *bit_LLR)
{
    // MIMO检测不变 只变解调 性能很好
    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d); // every symbol is transmitted over p parallel bit-channels.
    int m = p * 2;
    int M = M1d*M1d;

    int nr = arrH->arr[0]->rows; // number of receive antennas
    int nt = arrH->arr[0]->cols; // number of transmit antennas
    int N = arrH->size; // 时隙总量

    /* Step 1: ready the space*/
    //arrH_cha  H + eye(nt) -> H_cha (H_change)
    static ittc_cmat *H_cha = NULL;
    if (H_cha == NULL || H_cha->cols != nt || H_cha->rows != nr)
    {
        if (H_cha != NULL)
            ittc_cmat_free(H_cha);
        H_cha = ittc_cmat_alloc(nr, nt);
    }

    static ittc_cmat *H_mmse = NULL;
    if (H_mmse == NULL || H_mmse->cols != nr || H_mmse->rows != nt)
    {
        if (H_mmse != NULL)
            ittc_cmat_free(H_mmse);
        H_mmse = ittc_cmat_alloc(nt, nr);
    }

    static ittc_cmat *wh = NULL;
    if (wh == NULL || wh->cols != nt || wh->rows != nt)
    {
        if (wh != NULL)
            ittc_cmat_free(wh);
        wh = ittc_cmat_alloc(nt, nt);
    }

    static ittc_cmat *wy = NULL;
    if (wy == NULL || wy->rows != nt || wy->cols != 1)
    {
        if (wy != NULL)
            ittc_cmat_free(wy);
        wy = ittc_cmat_alloc(nt, 1);
    }
    ittc_cmat_set_zero(wy);

    static ittc_vec *pre_signal_prob = NULL;
    if (pre_signal_prob == NULL || pre_signal_prob->size != M1d)
    {
        if (pre_signal_prob != NULL)
            ittc_vec_free(pre_signal_prob);
        pre_signal_prob = ittc_vec_alloc(M1d);
    }
    ittc_vec_set_zero(pre_signal_prob);  // Pre Signal Prob


    static ittc_vec *signal_prob = NULL;
    if (signal_prob == NULL || signal_prob->size != M1d)
    {
        if (signal_prob != NULL)
            ittc_vec_free(signal_prob);
        signal_prob = ittc_vec_alloc(M1d);
    }
    ittc_vec_set_zero(signal_prob);  // Signal Prob


    /* Step 2: MIMO Detection and Demodulation */
    for (int timeslot = 0; timeslot < N; timeslot++)
    {
        ittc_cmat_memcpy(arrH->arr[timeslot], H_cha);
        // MMSE滤波矩阵计算
        ittc_cmat_calmmse(H_cha, nvar, 1.0, H_mmse);  // MMSE矩阵
        ittc_cmat_mult(H_mmse, H_cha, wh);
        ittc_cmat_mult(H_mmse, rx_signal->arr[timeslot], wy);

        // 逐天线计算
        for (int i = 0; i < nt; i++)
        {
            double tmp1 = 0;  // 干扰
            for (int k = 0; k < nt; k++)
            {
                if (k != i)
                {
                    tmp1 += ittc_complex_abs2(ittc_cmat_get(wh, i, k));
                }
            }

            double tmp2 = 0;  // 噪声
            for (int k = 0; k < nr; k++)
            {
                tmp2 += ittc_complex_abs2(ittc_cmat_get(H_mmse, i, k));
            }
            tmp2 = tmp2*nvar;

            double tmp3 = ittc_complex_abs2(ittc_cmat_get(wh, i, i));  // 有用信号
            double SINR = tmp3 / (tmp1 + tmp2);
            ittc_complex signal = ittc_div_cc(ittc_cmat_get(wy, i, 0), ittc_cmat_get(wh, i, i));
            double Eq_nvar = 1.0 / SINR;

            // Demodulation
            // Real Part
            double p_sum = -polar_code_large_number;
            ittc_vec_set_zero(pre_signal_prob);
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算  vn_do->size == A->M
            {
                // 计算符号的先验概率
                for (int bit_idx = 0; bit_idx < p; bit_idx++)
                {
                    int bit_tmp = (k >> (p - bit_idx - 1)) & 1;
                    double lp0, lp1;
                    ittc_polar_code_llrs_to_logprobs(decoded_LLR->dat[timeslot*m*nt + i*m + bit_idx], &lp0, &lp1);
                    pre_signal_prob->dat[k] += (bit_tmp == 0 ? lp0 : lp1);  //先验概率
                }
                p_sum = ittc_polar_code_max_star(p_sum, pre_signal_prob->dat[k]);  // attention value range
            }

            //归一化
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算
                pre_signal_prob->dat[k] -= p_sum;

            for (int k = 0; k < M1d; k++)
                signal_prob->dat[k] = pre_signal_prob->dat[k] - ((pow((signal.real - mapper->dat[k]), 2)) / Eq_nvar);


            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k] = p0 - p1 - decoded_LLR->dat[timeslot*m*nt + i*m + k];  // LLR Real Part
            }

            // Imag Part
            p_sum = -polar_code_large_number;
            ittc_vec_set_zero(pre_signal_prob);
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算  vn_do->size == A->M
            {
                // 计算符号的先验概率
                for (int bit_idx = 0; bit_idx < p; bit_idx++)
                {
                    int bit_tmp = (k >> (p - bit_idx - 1)) & 1;
                    double lp0, lp1;
                    ittc_polar_code_llrs_to_logprobs(decoded_LLR->dat[timeslot*m*nt + i*m + bit_idx + p], &lp0, &lp1);
                    pre_signal_prob->dat[k] += (bit_tmp == 0 ? lp0 : lp1);  //先验概率
                }
                p_sum = ittc_polar_code_max_star(p_sum, pre_signal_prob->dat[k]);  // attention value range
            }

            //归一化
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算
                pre_signal_prob->dat[k] -= p_sum;

            for (int k = 0; k < M1d; k++)
                signal_prob->dat[k] = pre_signal_prob->dat[k] - ((pow((signal.imag - mapper->dat[k]), 2)) / Eq_nvar);


            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k + p] = p0 - p1 - decoded_LLR->dat[timeslot*m*nt + i*m + k + p];  // LLR Imag Part
            }
        }
    }
}

void ittc_LDPC_code_mimo_MQAM_detector_demodulator_MMSE_PIC_Complex(ittc_arr_cmat *rx_signal, ittc_vec *decoded_LLR, ittc_vec *mapper, ittc_arr_cmat *arrH, double nvar,
                                                                    ittc_vec *bit_LLR)
{
    // MIMO检测改变 软估计PIC 性能不好
    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d); // every symbol is transmitted over p parallel bit-channels.
    int m = p * 2;
    int M = M1d*M1d;

    int nr = arrH->arr[0]->rows; // number of receive antennas
    int nt = arrH->arr[0]->cols; // number of transmit antennas
    int N = arrH->size; // 时隙总量

    /* Step 1: ready the space*/
    //arrH_cha  H + eye(nt) -> H_cha (H_change)
    static ittc_cmat *H_cha = NULL;
    if (H_cha == NULL || H_cha->cols != 1 || H_cha->rows != nr)
    {
        if (H_cha != NULL)
            ittc_cmat_free(H_cha);
        H_cha = ittc_cmat_alloc(nr, 1);
    }

    static ittc_cmat *H_mmse = NULL;
    if (H_mmse == NULL || H_mmse->cols != nr || H_mmse->rows != 1)
    {
        if (H_mmse != NULL)
            ittc_cmat_free(H_mmse);
        H_mmse = ittc_cmat_alloc(1, nr);
    }

    static ittc_cmat *wh = NULL;
    if (wh == NULL || wh->cols != 1 || wh->rows != 1)
    {
        if (wh != NULL)
            ittc_cmat_free(wh);
        wh = ittc_cmat_alloc(1, 1);
    }

    static ittc_cmat *wy = NULL;
    if (wy == NULL || wy->rows != 1 || wy->cols != 1)
    {
        if (wy != NULL)
            ittc_cmat_free(wy);
        wy = ittc_cmat_alloc(1, 1);
    }
    ittc_cmat_set_zero(wy);

    static ittc_arr_vec *pre_signal_prob = NULL;
    if (pre_signal_prob == NULL || pre_signal_prob->size != N*nt * 2 || pre_signal_prob->arr[0]->size != M1d)
    {
        if (pre_signal_prob != NULL)
            ittc_arr_vec_free(pre_signal_prob);
        pre_signal_prob = ittc_arr_vec_alloc(N*nt * 2);
        for (int i = 0; i < pre_signal_prob->size; i++)
            pre_signal_prob->arr[i] = ittc_vec_alloc(M1d);
    }
    for (int i = 0; i < pre_signal_prob->size; i++)
        ittc_vec_set_zero(pre_signal_prob->arr[i]);  // Pre Signal Prob

    static ittc_vec *signal_prob = NULL;
    if (signal_prob == NULL || signal_prob->size != M1d)
    {
        if (signal_prob != NULL)
            ittc_vec_free(signal_prob);
        signal_prob = ittc_vec_alloc(M1d);
    }
    ittc_vec_set_zero(signal_prob);  // Signal Prob


    static ittc_cmat *xtmp = NULL;
    if (xtmp == NULL || xtmp->rows != nt || xtmp->cols != 1)
    {
        if (xtmp != NULL)
            ittc_cmat_free(xtmp);
        xtmp = ittc_cmat_alloc(nt, 1);
    }
    ittc_cmat_set_zero(xtmp);

    static ittc_cmat *xtmp_change = NULL;
    if (xtmp_change == NULL || xtmp_change->rows != nt || xtmp_change->cols != 1)
    {
        if (xtmp_change != NULL)
            ittc_cmat_free(xtmp_change);
        xtmp_change = ittc_cmat_alloc(nt, 1);
    }
    ittc_cmat_set_zero(xtmp_change);

    static ittc_cmat *ytmp = NULL;
    if (ytmp == NULL || ytmp->rows != nr || ytmp->cols != 1)
    {
        if (ytmp != NULL)
            ittc_cmat_free(ytmp);
        ytmp = ittc_cmat_alloc(nr, 1);
    }
    ittc_cmat_set_zero(ytmp);


    static ittc_cmat *ytmp_change = NULL;
    if (ytmp_change == NULL || ytmp_change->rows != nr || ytmp_change->cols != 1)
    {
        if (ytmp_change != NULL)
            ittc_cmat_free(ytmp_change);
        ytmp_change = ittc_cmat_alloc(nr, 1);
    }
    ittc_cmat_set_zero(ytmp_change);

    // 概率计算
    for (int timeslot = 0; timeslot < N; timeslot++)
    {
        for (int i = 0; i < nt; i++)
        {
            double p_sum = -polar_code_large_number;
            ittc_vec_set_zero(pre_signal_prob->arr[timeslot*nt * 2 + i * 2]);
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算  vn_do->size == A->M
            {
                // 计算符号的先验概率
                for (int bit_idx = 0; bit_idx < p; bit_idx++)
                {
                    int bit_tmp = (k >> (p - bit_idx - 1)) & 1;
                    double lp0, lp1;
                    ittc_polar_code_llrs_to_logprobs(decoded_LLR->dat[timeslot*m*nt + i*m + bit_idx], &lp0, &lp1);
                    pre_signal_prob->arr[timeslot*nt * 2 + i * 2]->dat[k] += (bit_tmp == 0 ? lp0 : lp1);  //先验概率
                }
                p_sum = ittc_polar_code_max_star(p_sum, pre_signal_prob->arr[timeslot*nt * 2 + i * 2]->dat[k]);  // attention value range
            }

            //归一化
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算
                pre_signal_prob->arr[timeslot*nt * 2 + i * 2]->dat[k] -= p_sum;


            p_sum = -polar_code_large_number;
            ittc_vec_set_zero(pre_signal_prob->arr[timeslot*nt * 2 + i * 2 + 1]);
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算  vn_do->size == A->M
            {
                // 计算符号的先验概率
                for (int bit_idx = 0; bit_idx < p; bit_idx++)
                {
                    int bit_tmp = (k >> (p - bit_idx - 1)) & 1;
                    double lp0, lp1;
                    ittc_polar_code_llrs_to_logprobs(decoded_LLR->dat[timeslot*m*nt + i*m + bit_idx + p], &lp0, &lp1);
                    pre_signal_prob->arr[timeslot*nt * 2 + i * 2 + 1]->dat[k] += (bit_tmp == 0 ? lp0 : lp1);  //先验概率
                }
                p_sum = ittc_polar_code_max_star(p_sum, pre_signal_prob->arr[timeslot*nt * 2 + i * 2 + 1]->dat[k]);  // attention value range
            }

            //归一化
            for (int k = 0; k < M1d; k++)   // 符号先验概率计算
                pre_signal_prob->arr[timeslot*nt * 2 + i * 2 + 1]->dat[k] -= p_sum;
        }
    }

    /* Step 2: MIMO Detection and Demodulation */
    for (int timeslot = 0; timeslot < N; timeslot++)
    {
        // 信号重构
        ittc_cmat_set_zero(xtmp);
        for (int j = 0; j < nt; j++)  // 天线序号
        {
            int sym_idx;
            ittc_complex c;
            c.real = 0;
            c.imag = 0;
            for (int k = 0; k < M1d; k++)  // real part 求均值
                c.real += exp(pre_signal_prob->arr[timeslot*nt * 2 + j * 2]->dat[k])*mapper->dat[k];

            for (int k = 0; k < M1d; k++)  // real part 求均值
                c.imag += exp(pre_signal_prob->arr[timeslot*nt * 2 + j * 2 + 1]->dat[k])*mapper->dat[k];

            ittc_cmat_set(xtmp, j, 0, c);
        }

        // PIC
        for (int i = 0; i < nt; i++)  // 天线序号
        {
            ittc_cmat_memcpy(xtmp, xtmp_change);
            ittc_complex tmp;
            tmp.real = 0;
            tmp.imag = 0;
            ittc_cmat_set(xtmp_change, i, 0, tmp);  // 待求位置置零

            ittc_cmat_mult(arrH->arr[timeslot], xtmp_change, ytmp);
            // 信号相减 信号配置
            for (int k = 0; k < nr; k++)  //sub
            {
                ittc_cmat_set(ytmp_change, k, 0, ittc_sub_cc(ittc_cmat_get(rx_signal->arr[timeslot], k, 0), ittc_cmat_get(ytmp, k, 0)));
                ittc_cmat_set(H_cha, k, 0, ittc_cmat_get(arrH->arr[timeslot], k, i));
            }

            ittc_cmat_H(H_cha, H_mmse);  //MRC合并准备
            // MRC
            ittc_cmat_mult(H_mmse, H_cha, wh);
            ittc_cmat_mult(H_mmse, ytmp_change, wy);

            ittc_complex w1 = ittc_cmat_get(wh, 0, 0);
            ittc_complex w2 = ittc_cmat_get(wy, 0, 0);

            ittc_complex rec_sig_processed = ittc_div_cd(w2, w1.real);  // 处理后信号
            // 等效噪声计算
            double Eqvar = nvar / w1.real;   // 等效噪声

            // 信号解调 注意先验概率有偏置
            // Real Part
            for (int k = 0; k < M1d; k++)
            {
                signal_prob->dat[k] = pre_signal_prob->arr[timeslot*nt * 2 + i * 2]->dat[k] - ((pow((rec_sig_processed.real - mapper->dat[k]), 2)) / Eqvar);
            }

            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k] = p0 - p1 - decoded_LLR->dat[timeslot*m*nt + i*m + k];  // LLR Real Part
            }

            // Imag Part
            for (int k = 0; k < M1d; k++)
            {
                signal_prob->dat[k] = pre_signal_prob->arr[timeslot*nt * 2 + i * 2 + 1]->dat[k] - ((pow((rec_sig_processed.imag - mapper->dat[k]), 2)) / Eqvar);
            }

            for (int k = 0; k < p; k++)
            {
                double p0 = -polar_code_large_number;
                double p1 = -polar_code_large_number;
                for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
                {
                    if (((s_cnt >> (p - k - 1)) & 1) == 0)
                    {
                        p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                    }
                    else
                    {
                        p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                    }
                }
                bit_LLR->dat[timeslot*m*nt + i*m + k + p] = p0 - p1 - decoded_LLR->dat[timeslot*m*nt + i*m + k + p];  // LLR Imag Part
            }
        }
    }
}

void ittc_LDPC_code_mimo_MQAM_detector_demodulator_LSD(ittc_arr_cmat *rx_signal, ittc_vec *mapper, ittc_arr_cmat *arrH, double nvar,
                                                       ittc_vec *bit_LLR)
{
    // channel_idx只控制一维调制星座p个比特的映射关系 其余的由于各天线均为一致 所以不用 因为码长较长 使用自适应模式
    // 此过程对应一个编码器 控制所有数据流
    // 自带速率适配功能 调制码长为2的幂次
    // 线性SD
    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d); // every symbol is transmitted over p parallel bit-channels.
    int m = p * 2;
    int M = M1d*M1d;

    int nr = arrH->arr[0]->rows; // number of receive antennas
    int nt = arrH->arr[0]->cols; // number of transmit antennas
    int N_sym = arrH->size;
    int N_final = N_sym*m*nt;  // 实际传输比特数

    /* Step 1: ready the space*/
    static ittc_vec *mapper_processed = NULL;
    if (mapper_processed == NULL || mapper_processed->size != mapper->size)
    {
        if (mapper_processed != NULL)
            ittc_vec_free(mapper_processed);
        mapper_processed = ittc_vec_alloc(mapper->size);
    }

    // 星座映射表预处理
    double coefficient = sqrt(2.0*(M1d*M1d - 1) / 3.0);
    for (int i = 0; i < mapper->size; i++)
    {
        int costellation = (int)(mapper->dat[i] * coefficient);
        mapper_processed->dat[i] = costellation*1.0;   // 放大为正常
    }

    static ittc_mat *H_processed = NULL;
    if (H_processed == NULL || H_processed->rows != nr * 2 || H_processed->cols != nt * 2)
    {
        if (H_processed != NULL)
            ittc_mat_free(H_processed);
        H_processed = ittc_mat_alloc(nr * 2, nt * 2);
    }

    static ittc_mat *rx_signal_processed = NULL;
    if (rx_signal_processed == NULL || rx_signal_processed->rows != nr * 2 || rx_signal_processed->cols != 1)
    {
        if (rx_signal_processed != NULL)
            ittc_mat_free(rx_signal_processed);
        rx_signal_processed = ittc_mat_alloc(nr * 2, 1);
    }

    static ittc_vec *bit_LLR_processed = NULL;
    if (bit_LLR_processed == NULL || bit_LLR_processed->size != nt*m)
    {
        if (bit_LLR_processed != NULL)
            ittc_vec_free(bit_LLR_processed);
        bit_LLR_processed = ittc_vec_alloc(nt*m);
    }

    /* Step 2: MIMO Detection and Demodulation */
    for (int timeslot = 0; timeslot < N_sym; timeslot++)
    {
        // 接收信号能量预处理
        for (int i = 0; i < nr; i++)
        {
            ittc_complex tmp = ittc_cmat_get(rx_signal->arr[timeslot], i, 0);
            ittc_mat_set(rx_signal_processed, i, 0, tmp.real*coefficient);
            ittc_mat_set(rx_signal_processed, i + nr, 0, tmp.imag*coefficient);
        }

        // 信道矩阵预处理
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nt; j++)
            {
                ittc_complex tmp = ittc_cmat_get(arrH->arr[timeslot], i, j);
                ittc_mat_set(H_processed, i, j, tmp.real);
                ittc_mat_set(H_processed, i, j + nt, -tmp.imag);
                ittc_mat_set(H_processed, i + nr, j, tmp.imag);
                ittc_mat_set(H_processed, i + nr, j + nt, tmp.real);
            }
        }

        double N0_processed = nvar*coefficient*coefficient;

        // 球译码检测
        SE_List_LLR(1e4, 10, nr * 2, nt * 2, N0_processed, mapper_processed, rx_signal_processed, H_processed, bit_LLR_processed);

        // 似然比分配
        for (int i = 0; i < nt; i++)
        {
            // Real Part
            for (int k = 0; k < p; k++)
                bit_LLR->dat[timeslot*m*nt + i*m + k] = bit_LLR_processed->dat[i*p + k];

            // Imag Part
            for (int k = 0; k < p; k++)
                bit_LLR->dat[timeslot*m*nt + i*m + k + p] = bit_LLR_processed->dat[i*p + k + p*nt];
        }
    }
}


// SIC模式
void ittc_LDPC_code_mimo_MQAM_transmit_SIC(int nt, ittc_vec *mapper, ittc_arr_ivec *Layer_rate_match_bits, ittc_arr_cmat *tx_signal)
{
    /* Asumption: the average energy of the transmitted symbols is normalized to one.
    Parameter List:
    > nt     : number of transmitting antennas.
    > mapper : 1-dimensional PAM mapper with modulation order m/2.
    ...
    */

    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d);
    int m = p * 2;
    int N = Layer_rate_match_bits->arr[0]->size / m; // here, N is the code length of the sub-codes

    static ittc_ivec *tmp_mod = NULL;
    if (tmp_mod == NULL || tmp_mod->size != m / 2)
    {
        if (tmp_mod != NULL)
            ittc_ivec_free(tmp_mod);
        tmp_mod = ittc_ivec_alloc(m / 2);
    }


    for (int i = 0; i < N; i++)  // 时隙
    {
        for (int j = 0; j < nt; j++)  // 天线序号
        {
            int sym_idx;
            ittc_complex c;
            for (int k = 0; k < p; k++)  // real part
            {
                tmp_mod->dat[k] = ittc_ivec_get(Layer_rate_match_bits->arr[j], i*m + k);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.real = ittc_vec_get(mapper, sym_idx);

            for (int k = 0; k < p; k++)  // imag part
            {
                tmp_mod->dat[k] = ittc_ivec_get(Layer_rate_match_bits->arr[j], i*m + k + p);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.imag = ittc_vec_get(mapper, sym_idx);
            ittc_cmat_set(tx_signal->arr[i], j, 0, c);
        }
    }
}

void ittc_LDPC_code_mimo_MQAM_detector_demodulator_MMSE_SIC(ittc_arr_cmat *rx_signal, ittc_arr_ivec *Layer_decoded_bits, ittc_vec *mapper, ittc_arr_cmat *arrH, double nvar, int now_detecting_idx,
                                                            ittc_vec *Layer_bit_LLR)
{
    // now_detecting_idx指示正在译码的数据流
    // 译码顺序 最后一根天线开始 第一根天线结束
    int M1d = mapper->size;
    int p = ittc_is_power_of_2(M1d); // every symbol is transmitted over p parallel bit-channels.
    int m = p * 2;
    int M = M1d*M1d;

    int nr = arrH->arr[0]->rows; // number of receive antennas
    int nt = arrH->arr[0]->cols; // number of transmit antennas
    int N = arrH->size; // 时隙总量

    /* Step 1: ready the space*/
    //arrH_cha  H + eye(nt) -> H_cha (H_change)
    static ittc_arr_cmat *arrH_cha = NULL;
    if (arrH_cha == NULL || arrH_cha->size != nt || arrH_cha->arr[0]->rows != nr || arrH_cha->arr[0]->cols != 1)
    {
        if (arrH_cha != NULL)
            ittc_arr_cmat_free(arrH_cha);
        arrH_cha = ittc_arr_cmat_alloc(nt);
        for (int i = 0; i < arrH_cha->size; i++)
            arrH_cha->arr[i] = ittc_cmat_alloc(nr, i + 1);  // for MMSE
    }
    for (int i = 0; i < arrH_cha->size; i++)
        ittc_cmat_set_zero(arrH_cha->arr[i]);


    static ittc_arr_cmat *arrH_mmse = NULL;
    if (arrH_mmse == NULL || arrH_mmse->size != nt || arrH_mmse->arr[0]->rows != 1 || arrH_cha->arr[0]->cols != nr)
    {
        if (arrH_mmse != NULL)
            ittc_arr_cmat_free(arrH_mmse);
        arrH_mmse = ittc_arr_cmat_alloc(nt);
        for (int i = 0; i < arrH_mmse->size; i++)
            arrH_mmse->arr[i] = ittc_cmat_alloc(i + 1, nr);  // for MMSE
    }
    for (int i = 0; i < arrH_mmse->size; i++)
        ittc_cmat_set_zero(arrH_mmse->arr[i]);

    static ittc_arr_cmat *arrwh = NULL;
    if (arrwh == NULL || arrwh->size != nt || arrwh->arr[0]->rows != 1 || arrwh->arr[0]->cols != 1)
    {
        if (arrwh != NULL)
            ittc_arr_cmat_free(arrwh);
        arrwh = ittc_arr_cmat_alloc(nt);
        for (int i = 0; i < arrwh->size; i++)
            arrwh->arr[i] = ittc_cmat_alloc(i + 1, i + 1);  // for wh
    }
    for (int i = 0; i < arrwh->size; i++)
        ittc_cmat_set_zero(arrwh->arr[i]);

    //arry
    static ittc_arr_cvec *arry = NULL;
    if (arry == NULL || arry->size != nt || arry->arr[0]->size != nr)
    {
        if (arry != NULL)
            ittc_arr_cvec_free(arry);
        arry = ittc_arr_cvec_alloc(nt);
        for (int i = 0; i < arry->size; i++)
            arry->arr[i] = ittc_cvec_alloc(nr);
    }
    for (int i = 0; i < arry->size; i++)
        ittc_cvec_set_zero(arry->arr[i]);

    static ittc_cmat *wy = NULL;
    if (wy == NULL || wy->rows != 1 || wy->cols != 1)
    {
        if (wy != NULL)
            ittc_cmat_free(wy);
        wy = ittc_cmat_alloc(1, 1);
    }
    ittc_cmat_set_zero(wy);


    static ittc_ivec *tmp_mod = NULL;
    if (tmp_mod == NULL || tmp_mod->size != m / 2)
    {
        if (tmp_mod != NULL)
            ittc_ivec_free(tmp_mod);
        tmp_mod = ittc_ivec_alloc(m / 2);
    }

    static ittc_vec *signal_prob = NULL;
    if (signal_prob == NULL || signal_prob->size != M1d)
    {
        if (signal_prob != NULL)
            ittc_vec_free(signal_prob);
        signal_prob = ittc_vec_alloc(M1d);
    }
    ittc_vec_set_zero(signal_prob);  // Signal Prob


    static ittc_cmat *xtmp = NULL;
    if (xtmp == NULL || xtmp->rows != nt || xtmp->cols != 1)
    {
        if (xtmp != NULL)
            ittc_cmat_free(xtmp);
        xtmp = ittc_cmat_alloc(nt, 1);
    }
    ittc_cmat_set_zero(xtmp);

    static ittc_cmat *ytmp = NULL;
    if (ytmp == NULL || ytmp->rows != nr || ytmp->cols != 1)
    {
        if (ytmp != NULL)
            ittc_cmat_free(ytmp);
        ytmp = ittc_cmat_alloc(nr, 1);
    }
    ittc_cmat_set_zero(ytmp);


    /* Step 2: MIMO Detection and Demodulation */
    for (int timeslot = 0; timeslot < N; timeslot++)
    {
        ittc_cmat_set_zero(arrH_cha->arr[now_detecting_idx]);
        for (int cols = 0; cols <= now_detecting_idx; cols++)  //copy data
        {
            for (int rows = 0; rows < nr; rows++)
            {
                ittc_cmat_set(arrH_cha->arr[now_detecting_idx], rows, cols, ittc_cmat_get(arrH->arr[timeslot], rows, cols));
            }
        }

        ittc_cmat_calmmse(arrH_cha->arr[now_detecting_idx], nvar, 1.0, arrH_mmse->arr[now_detecting_idx]);  // cal mmse matrix
        ittc_cmat_mult(arrH_mmse->arr[now_detecting_idx], arrH_cha->arr[now_detecting_idx], arrwh->arr[now_detecting_idx]);  // cal w*h

        double tmp1 = 0;   // interference   
        for (int s = 0; s < now_detecting_idx; s++)
        {
            ittc_complex r = ittc_cmat_get(arrwh->arr[now_detecting_idx], now_detecting_idx, s);
            tmp1 += r.real*r.real + r.imag*r.imag;
        }

        double tmp2 = 0;   // noise
        for (int s = 0; s < nr; s++)
        {
            ittc_complex r = ittc_cmat_get(arrH_mmse->arr[now_detecting_idx], now_detecting_idx, s);
            tmp2 += r.real*r.real + r.imag*r.imag;
        }

        ittc_complex tmp3 = ittc_cmat_get(arrwh->arr[now_detecting_idx], now_detecting_idx, now_detecting_idx);

        double nvar_1D = 0.5*((nvar*tmp2 + tmp1) / tmp3.real / tmp3.real);   // 等效噪声

        // 信号重构
        ittc_cmat_set_zero(xtmp);
        for (int j = nt - 1; j > now_detecting_idx; j--)  // 天线序号
        {
            int sym_idx;
            ittc_complex c;
            for (int k = 0; k < p; k++)  // real part
            {
                tmp_mod->dat[k] = ittc_ivec_get(Layer_decoded_bits->arr[j], timeslot*m + k);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.real = ittc_vec_get(mapper, sym_idx);

            for (int k = 0; k < p; k++)  // imag part
            {
                tmp_mod->dat[k] = ittc_ivec_get(Layer_decoded_bits->arr[j], timeslot*m + k + p);
            }
            sym_idx = ittc_ivec2int(tmp_mod);
            c.imag = ittc_vec_get(mapper, sym_idx);
            ittc_cmat_set(xtmp, j, 0, c);
        }

        ittc_cmat_mult(arrH->arr[timeslot], xtmp, ytmp);
        for (int k = 0; k < nr; k++)  //sub
        {
            ittc_cvec_set(arry->arr[now_detecting_idx], k, ittc_sub_cc(ittc_cmat_get(rx_signal->arr[timeslot], k, 0), ittc_cmat_get(ytmp, k, 0)));
        }

        ittc_complex w;
        w.real = 0;
        w.imag = 0;
        for (int k = 0; k < nr; k++)
        {
            w = ittc_add_cc(w, ittc_mult_cc(ittc_cmat_get(arrH_mmse->arr[now_detecting_idx], now_detecting_idx, k), ittc_cvec_get(arry->arr[now_detecting_idx], k)));
        }

        ittc_complex rec_sig_processed = ittc_div_cd(w, tmp3.real);  // 处理后信号

        // 信号解调
        // Real Part
        for (int k = 0; k < M1d; k++)
        {
            signal_prob->dat[k] = -(pow((rec_sig_processed.real - mapper->dat[k]), 2)) / nvar_1D / 2.0;
        }

        for (int k = 0; k < p; k++)
        {
            double p0 = -polar_code_large_number;
            double p1 = -polar_code_large_number;
            for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
            {
                if (((s_cnt >> (p - k - 1)) & 1) == 0)
                {
                    p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                }
                else
                {
                    p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                }
            }
            Layer_bit_LLR->dat[timeslot*m + k] = p0 - p1;  // LLR Real Part
        }

        // Imag Part
        for (int k = 0; k < M1d; k++)
        {
            signal_prob->dat[k] = -(pow((rec_sig_processed.imag - mapper->dat[k]), 2)) / nvar_1D / 2.0;
        }

        for (int k = 0; k < p; k++)
        {
            double p0 = -polar_code_large_number;
            double p1 = -polar_code_large_number;
            for (int s_cnt = 0; s_cnt < M1d; s_cnt++)
            {
                if (((s_cnt >> (p - k - 1)) & 1) == 0)
                {
                    p0 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p0);
                }
                else
                {
                    p1 = ittc_polar_code_max_star(signal_prob->dat[s_cnt], p1);
                }
            }
            Layer_bit_LLR->dat[timeslot*m + k + p] = p0 - p1;  // LLR Real Part
        }
    }

}