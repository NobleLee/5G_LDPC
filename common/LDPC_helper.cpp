//
// Created by ggl on 2018/1/14.
//
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdlib>
#include<fstream>
#include <cmath>
#include "LDPC_helper.h"
#include "LDPC_5G.h"


using namespace std;

/**
* 将两个向量做模2加操作，结果放到v1中
* @param v1 向量1
* @param v2 向量2
*/
inline void vec_mod2sum(int *v1, int *v2, const int size) {
    for (int i = 0; i < size; i++) {
        if (v1[i] == v2[i]) v1[i] = 0;
        else v1[i] = 1;
    }
}

void readMatrixFromFile(string &filename, vector<vector<int>> &parityMats) {
    ifstream in(filename);
    string s;
    vector<string> tmp;
    while (getline(in, s)) {
        SplitString(s, tmp, "\t");
        if (tmp.size() != 10) {
            std::cout << "read parity matrix error!!! \nexit" << endl;
            system("pause");
        }
        vector<int> tl;
        for (int i = 0; i < tmp.size(); i++)
            tl.push_back(atoi(tmp[i].c_str()));
        parityMats.push_back(tl);
        tmp.clear();
    }
    in.close();
}

/**
* 对矩阵进行高斯消元，矩阵右侧变换为单位阵
* @param matrix
* @param row
* @param column
*/
void Gaussian_Elimination(int **matrix, const int row, const int column) {

    const unsigned long matrix_start = column - row;

    for (unsigned long j = matrix_start, i = 0; j < column; j++, i++) {

        if (matrix[i][j] == 0) {
            int temp = 0;
            for (int i1 = i; i1 < row; i1++) {
                if (matrix[i1][j] == 0) continue;
                int *temp = matrix[i1];
                matrix[i1] = matrix[i];
                matrix[i] = temp;
                break;
            }
        }

        for (unsigned long k = 0; k < row; k++) {
            if (matrix[k][j] == 1 && k != i) {
                vec_mod2sum(matrix[k], matrix[i], column);
            }
        }
    }
    // 校验结果，查看校验矩阵，是否满秩
    for (int i = 0; i < row; i++) {
        if (matrix[i][matrix_start + i] == 0) {
            std::cout << "the parity matrix is non-singular matrix!!! exit" << std::endl;
            system("pause");
        }
    }
}


/**
* 将字符串按照某一字符进行分割
* @param s 要分割的字符串
* @param v 保存分割结果的容器
* @param c 分割字符
*/
void SplitString(const string &s, vector<string> &v, const string &c) {
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (string::npos != pos2) {
        v.push_back(s.substr(pos1, pos2 - pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

/**
* 获取变量节点与校验节点的关系，将结果存放在req中，结果的第j个元素，是通过校验节点与第j个变量节点连接的变量节点的索引
* 用途:BP译码迭代用
* @param P_mats 扩展之后的校验矩阵
* @param row 校验矩阵的行数
* @param column 校验矩阵的列数
* @param req 存放结果
*/
void getEdgeFrom_VNandCN(int **P_mats, const unsigned long row, const unsigned long column, vector<vector<vector<int>>> &req) {
    /// req 存放的数据结构解释：
    /// req的第i个元素为vector<vector<int>> 设为req_i，他的大小表示与这个变量节点相连接的校验节点的数目
    /// 那么req_i的第j个元素为vector<int> 设为req_ij,表示通过[与第i个变量节点相连接的第j个校验节点]来传递到第i个变量节点的元素的集合
    /// req_ij的第k个元素为int型数，表示通过某一校验节点与第i个变量节点相连接的变量节点的索引
    for (unsigned long col_index = 0; col_index < column; col_index++) {
        vector<vector<int>> C_V;
        for (unsigned long row_index = 0; row_index < row; row_index++) {
            if (P_mats[row_index][col_index] == 0) continue;
            vector<int> V;
            for (int edge_index = 0; edge_index < column; edge_index++) {
                if (P_mats[row_index][edge_index] == 1 && edge_index != col_index)
                    V.push_back(edge_index);
            }
            if (V.size() != 0)
                C_V.push_back(V);
        }
        req.push_back(C_V);
    }
}

/**
* 按列提取校验矩阵中不为0的点的索引
* @param P_mats 校验矩阵
* @param row
* @param column
* @param req
*/
void getParityMatrixPoint(int **P_mats, const unsigned long row, const unsigned long column, vector<vector<int>> &req) {
    for (unsigned long i = 0; i < row; i++) {
        vector<int> tmp;
        for (int j = 0; j < column; j++) {
            if (P_mats[i][j] == 1) tmp.push_back(j);
        }
        req.push_back(tmp);
    }
}

/**
* 将H_base中的矩阵按照扩展因子扩展，结果存放在P_Mats
* @param P_Mats 存放扩展结果
* @param H_base 存放基础矩阵
* @param row 基础矩阵的行
* @param columns 基础矩阵的列
* @param zLength 扩展因子
*/
void expendParityMatrix(int **P_Mats, int **H_base, const int row, const int columns, const int zLength) {
    int ly = 0;
    for (unsigned long i = 0; i < row; i++) {
        for (unsigned long j = 0; j < columns; j++) {
            if (H_base[i][j] == -1) continue;
            unsigned long start_x = i * zLength, start_y = j * zLength;
            for (int k = 0; k < zLength; k++) {
                ly = start_y + (k + H_base[i][j]) % zLength;
                P_Mats[start_x + k][ly] = 1;
            }
        }
    }
}

/***
 * 三星-对基础矩阵进行扩展
 * @param LDPC_Matrix
 * @param BaseMatrix
 * @param row
 * @param columns
 * @param zLength
 */
void SAM_expendParityMatrix(unsigned short *LDPC_Matrix, long *BaseMatrix, const int row, const int columns, const int zLength) {
    unsigned long ultemp1, ultemp2, ultemp3, temp;
    for (ultemp1 = 0; ultemp1 < columns; ultemp1++) {
        for (ultemp2 = 0; ultemp2 < row; ultemp2++) {
            if (BaseMatrix[ultemp2 * columns + ultemp1] == -1) {
                continue;
            }
            for (ultemp3 = 0; ultemp3 < zLength; ultemp3++) {
                temp = (BaseMatrix[ultemp2 * columns + ultemp1] + ultemp3) % zLength;
                LDPC_Matrix[(ultemp1 * zLength + ultemp3) * row * zLength + ultemp2 * zLength + temp] = 1;
            }
        }
    }
}


/***
 * 高斯消元
 * @param t_LDPC_Matrix
 * @param t_P_Matrix
 * @param t_ulpReArrangeCol
 * @param t_ulCheLength
 * @param t_ulCodeLength
 * @return
 */
bool SAM_Gaussian_Elimination(unsigned short *t_LDPC_Matrix, unsigned short *t_P_Matrix, unsigned long *t_ulpReArrangeCol, unsigned long t_ulCheLength, unsigned long t_ulCodeLength) {
    unsigned long ulTemp1, ulTemp2, ulTemp3, ulMiddle1, ulMiddle2;

    for (ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++) {
        cout << ulTemp1 << endl;
        ulMiddle1 = ulTemp1;
        while (ulMiddle1 < t_ulCodeLength && t_LDPC_Matrix[ulTemp1 + ulMiddle1 * t_ulCheLength] == 0) {
            for (ulTemp2 = ulTemp1 + 1; ulTemp2 < t_ulCheLength; ulTemp2++) {
                if (t_LDPC_Matrix[ulTemp2 + ulMiddle1 * t_ulCheLength] != 0) {
                    break;
                }
            }
            ulMiddle1++;
        }

        if (ulMiddle1 >= t_ulCodeLength) {
            return 0;
        }

        if (ulMiddle1 != ulTemp1) {
            t_ulpReArrangeCol[ulTemp1] = ulMiddle1;
            for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++) {
                ulMiddle2 = t_LDPC_Matrix[ulTemp2 + ulTemp1 * t_ulCheLength];
                t_LDPC_Matrix[ulTemp2 + ulTemp1 * t_ulCheLength] = t_LDPC_Matrix[ulTemp2 + ulMiddle1 * t_ulCheLength];
                t_LDPC_Matrix[ulTemp2 + ulMiddle1 * t_ulCheLength] = ulMiddle2;
            }
        }

        if (t_LDPC_Matrix[ulTemp1 * t_ulCheLength + ulTemp1] == 0) {
            for (ulTemp2 = ulTemp1 + 1; ulTemp2 < t_ulCheLength; ulTemp2++) {
                if (t_LDPC_Matrix[ulTemp2 + ulTemp1 * t_ulCheLength] != 0) {
                    ulMiddle1 = ulTemp2;
                    for (ulTemp3 = ulTemp1; ulTemp3 < t_ulCodeLength; ulTemp3++) {
                        t_LDPC_Matrix[ulTemp1 + ulTemp3 * t_ulCheLength] = (t_LDPC_Matrix[ulMiddle1 + ulTemp3 * t_ulCheLength] + t_LDPC_Matrix[ulTemp1 + ulTemp3 * t_ulCheLength]) % 2;
                    }
                }
            }
        }

        for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++) {
            if (t_LDPC_Matrix[ulTemp1 * t_ulCheLength + ulTemp2] != 0) {
                if (ulTemp2 != ulTemp1) {
                    for (ulTemp3 = ulTemp1; ulTemp3 < t_ulCodeLength; ulTemp3++) {
                        t_LDPC_Matrix[ulTemp3 * t_ulCheLength + ulTemp2] = (t_LDPC_Matrix[ulTemp3 * t_ulCheLength + ulTemp2] + t_LDPC_Matrix[ulTemp3 * t_ulCheLength + ulTemp1]) % 2;
                    }
                }
            }
        }
    }

    for (ulTemp1 = t_ulCheLength; ulTemp1 < t_ulCodeLength; ulTemp1++) {
        for (ulTemp2 = 0; ulTemp2 < t_ulCheLength; ulTemp2++) {
            t_P_Matrix[(ulTemp1 - t_ulCheLength) * t_ulCheLength + ulTemp2] = t_LDPC_Matrix[ulTemp1 * t_ulCheLength + ulTemp2];
        }
    }

    return 1;
}


/**
 * 信道快速编码
 * @param t_bypInforbits 待编码bit
 * @param t_bypCodeWords 编码之后的bit
 * @param P_Matrix 校验矩阵
 * @param ulpReArrangeCol 需要进行行列交换的位置
 * @param t_ulCodeLength 输出的码字长度
 * @param t_ulCheLength  校验bit的长度
 * @return
 */
bool SAM_LDPC_Fast_Encoder(int *t_bypInforbits, int *t_bypCodeWords, unsigned short *P_Matrix,
                           unsigned long *ulpReArrangeCol, unsigned long t_ulCodeLength, unsigned long t_ulCheLength) {
    unsigned long ulInfBitLength = t_ulCodeLength - t_ulCheLength;
    unsigned long ulTemp1, ulTemp2, ulMiddle;

    unsigned short *ulCheckBits = new unsigned short[t_ulCheLength];

    for (ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++) {
        ulMiddle = 0;
        for (ulTemp2 = 0; ulTemp2 < ulInfBitLength; ulTemp2++) {
            ulMiddle += P_Matrix[ulTemp1 + ulTemp2 * t_ulCheLength] * t_bypInforbits[ulTemp2];
        }
        ulCheckBits[ulTemp1] = ulMiddle % 2;
    }

    for (ulTemp1 = 0; ulTemp1 < t_ulCodeLength; ulTemp1++) {
        if (ulTemp1 < t_ulCheLength) {
            t_bypCodeWords[ulTemp1] = ulCheckBits[ulTemp1];
        } else {
            t_bypCodeWords[ulTemp1] = t_bypInforbits[ulTemp1 - t_ulCheLength];
        }
    }

    for (long ulTemp1 = t_ulCheLength - 1; ulTemp1 >= 0; ulTemp1--) {
        if (ulpReArrangeCol[ulTemp1] != 0) {
            ulMiddle = t_bypCodeWords[ulTemp1];
            t_bypCodeWords[ulTemp1] = t_bypCodeWords[ulpReArrangeCol[ulTemp1]];
            t_bypCodeWords[ulpReArrangeCol[ulTemp1]] = ulMiddle;
        }
    }

    delete[] ulCheckBits;

    return 1;
}

bool ExtractInforBits(double *t_bypCodeWords, unsigned long *ulpReArrangeCol, unsigned short *t_bypInforbits, unsigned long t_ulCodeLength, unsigned long t_ulCheLength) {
    unsigned long ulMiddle;
    for (long ulTemp1 = 0; ulTemp1 < t_ulCheLength; ulTemp1++) {
        if (ulpReArrangeCol[ulTemp1] != 0) {
            ulMiddle = t_bypCodeWords[ulTemp1];
            t_bypCodeWords[ulTemp1] = t_bypCodeWords[ulpReArrangeCol[ulTemp1]];
            t_bypCodeWords[ulpReArrangeCol[ulTemp1]] = ulMiddle;
        }
    }

    for (unsigned long ulTemp1 = t_ulCheLength; ulTemp1 < t_ulCodeLength; ulTemp1++) {
        t_bypInforbits[ulTemp1 - t_ulCheLength] = t_bypCodeWords[ulTemp1];
    }
    return 1;
}


/**
 * Decoding in AWGN Channel
 * @param t_lpVarDis 监督矩阵
 * @param t_lpCheDis 监督矩阵
 * @param t_dpChannelOut 信道输出，译码器输入（此处为LLR值）
 * @param t_dpDecoding 译码器输出
 * @param t_dpLLR  输出LLR值(外信息)
 * @param t_ulCodeLength 码长
 * @param t_ulCheLength 校验方程数
 * @param t_byVarDeg 变量节点度
 * @param t_byCheDeg 校验节点度
 * @param t_ulIterMax  最大迭代次数
 * @param bypInforbits
 * @param ulpReArrangeCol
 * @return 满足所有校验方程返回1，否则返回0
 */
bool Decoder_AWGN(long *t_lpVarDis, long *t_lpCheDis, double *t_dpChannelOut, double *t_dpDecoding, double *t_dpLLR, unsigned long t_ulCodeLength,
                  unsigned long t_ulCheLength, unsigned short t_byVarDeg, unsigned short t_byCheDeg, unsigned long t_ulIterMax, unsigned short *bypInforbits, unsigned long *ulpReArrangeCol) {

    unsigned short wVarDeg, wCheDeg;
    unsigned long ulCount, ulNum;
    unsigned long ulInfLength = t_ulCodeLength - t_ulCheLength;
    long lCheIndex, lVarIndex;

    wVarDeg = t_byVarDeg;
    wCheDeg = t_byCheDeg;
    double *dpLLRRec = (double *) malloc(t_ulCodeLength * sizeof(double));
//  t_dpCheSta 校验节点状态，大小为校验节点长，作为先验信息输入，窗译码使用，译整个码直接全设为1
    double *t_dpCheSta = (double *) malloc(t_ulCheLength * sizeof(double));
    double *dpVar2Che = (double *) malloc(t_ulCheLength * wCheDeg * sizeof(double));
    double dMiddle = 0;
    int *bypSyn = (int *) malloc(t_ulCheLength * sizeof(int));

    unsigned char byMiddle;

    double *dpProductForward = (double *) malloc(t_ulCheLength * wCheDeg * sizeof(double));
    double *dpProductBackward = (double *) malloc(t_ulCheLength * wCheDeg * sizeof(double));
    double *dpMiddle = (double *) malloc(t_ulCodeLength * wVarDeg * sizeof(double));
    double *dpSumForward = (double *) malloc(t_ulCodeLength * wVarDeg * sizeof(double));
    double *dpSumBackward = (double *) malloc(t_ulCodeLength * wVarDeg * sizeof(double));

//------------------
    for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength * wVarDeg; ulTemp++) {
        dpMiddle[ulTemp] = 0.0;
    }
//--------------------
    for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++) {
        dpLLRRec[ulTemp] = t_dpChannelOut[ulTemp];
    }
    for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++) {
        t_dpCheSta[ulTemp] = 1.0;
        for (unsigned short wDegTemp = 0; wDegTemp < wCheDeg; wDegTemp++) {
            lVarIndex = (long) t_lpCheDis[ulTemp * wCheDeg + wDegTemp];
            dpVar2Che[ulTemp * wCheDeg + wDegTemp] = 0.0;
            if (lVarIndex >= 0) {
                dpVar2Che[ulTemp * wCheDeg + wDegTemp] += dpLLRRec[lVarIndex];
            }
        }
    }


    ulCount = 0;
    while (1) {
        ulCount++;
        for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++) {
            dpProductForward[ulTemp * wCheDeg] = t_dpCheSta[ulTemp];
            dpProductBackward[ulTemp * wCheDeg] = 1.0;
            for (unsigned short wDegTemp = 1; wDegTemp < wCheDeg; wDegTemp++) {
                lVarIndex = (long) t_lpCheDis[ulTemp * wCheDeg + wDegTemp - 1];
                if (lVarIndex >= 0) {
                    dpProductForward[ulTemp * wCheDeg + wDegTemp] = dpProductForward[ulTemp * wCheDeg + wDegTemp - 1] * tanh(0.5 * dpVar2Che[ulTemp * wCheDeg + wDegTemp - 1]);
                } else {
                    dpProductForward[ulTemp * wCheDeg + wDegTemp] = dpProductForward[ulTemp * wCheDeg + wDegTemp - 1];
                }
                lVarIndex = (long) t_lpCheDis[(ulTemp + 1) * wCheDeg - wDegTemp];
                if (lVarIndex >= 0) {
                    dpProductBackward[ulTemp * wCheDeg + wDegTemp] = dpProductBackward[ulTemp * wCheDeg + wDegTemp - 1] * tanh(0.5 * dpVar2Che[(ulTemp + 1) * wCheDeg - wDegTemp]);
                } else {
                    dpProductBackward[ulTemp * wCheDeg + wDegTemp] = dpProductBackward[ulTemp * wCheDeg + wDegTemp - 1];
                }
            }

        }

        for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++) {
            for (unsigned short wDegTemp = 0; wDegTemp < wVarDeg; wDegTemp++) {
                lCheIndex = (long) t_lpVarDis[ulTemp * wVarDeg + wDegTemp];
                if (lCheIndex >= 0) {
                    for (unsigned short wTemp = 0; wTemp < wCheDeg; wTemp++) {
                        if ((long) t_lpCheDis[lCheIndex * wCheDeg + wTemp] == ulTemp) {
                            dMiddle = dpProductForward[wCheDeg * lCheIndex + wTemp] * dpProductBackward[wCheDeg * (lCheIndex + 1) - wTemp - 1];
                            break;
                        }
                    }
                    if ((dMiddle + 1 < GAP) && (dMiddle + 1 >= 0)) {
                        dpMiddle[ulTemp * wVarDeg + wDegTemp] = -1e3;
                    } else if ((1 - dMiddle < GAP) && (1 - dMiddle >= 0)) {
                        dpMiddle[ulTemp * wVarDeg + wDegTemp] = 1e3;
                    } else {
                        dpMiddle[ulTemp * wVarDeg + wDegTemp] = log((1 + dMiddle) / (1 - dMiddle));
                    }
                }
            }
        }

        for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++) {
            dpSumForward[wVarDeg * ulTemp] = 0.0;

            dpSumBackward[wVarDeg * ulTemp] = 0.0;

            for (unsigned short wDegTemp = 1; wDegTemp < wVarDeg; wDegTemp++) {
                if ((long) t_lpVarDis[ulTemp * wVarDeg + wDegTemp - 1] >= 0) {
                    dpSumForward[ulTemp * wVarDeg + wDegTemp] = dpSumForward[ulTemp * wVarDeg + wDegTemp - 1] + dpMiddle[ulTemp * wVarDeg + wDegTemp - 1];
                } else {
                    dpSumForward[ulTemp * wVarDeg + wDegTemp] = dpSumForward[ulTemp * wVarDeg + wDegTemp - 1];
                }
                if ((long) t_lpVarDis[(ulTemp + 1) * wVarDeg - wDegTemp] >= 0) {
                    dpSumBackward[ulTemp * wVarDeg + wDegTemp] = dpSumBackward[ulTemp * wVarDeg + wDegTemp - 1] + dpMiddle[(ulTemp + 1) * wVarDeg - wDegTemp];
                } else {
                    dpSumBackward[ulTemp * wVarDeg + wDegTemp] = dpSumBackward[ulTemp * wVarDeg + wDegTemp - 1];
                }
            }
        }

        for (unsigned short wDegTemp = 0; wDegTemp < wCheDeg; wDegTemp++) {
            for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++) {
                lVarIndex = (long) t_lpCheDis[ulTemp * wCheDeg + wDegTemp];
                if (lVarIndex >= 0) {
                    for (unsigned short wTemp = 0; wTemp < wVarDeg; wTemp++) {
                        lCheIndex = (long) t_lpVarDis[lVarIndex * wVarDeg + wTemp];
                        if (lCheIndex == ulTemp) {
                            dpVar2Che[ulTemp * wCheDeg + wDegTemp] = dpSumForward[lVarIndex * wVarDeg + wTemp] + dpSumBackward[(lVarIndex + 1) * wVarDeg - wTemp - 1] + dpLLRRec[lVarIndex];
                            break;
                        }
                    }
                }
            }
        }

        for (unsigned long ulTemp = 0; ulTemp < t_ulCodeLength; ulTemp++) {
            t_dpLLR[ulTemp] = dpSumForward[ulTemp * wVarDeg + 1] + dpSumBackward[(ulTemp + 1) * wVarDeg - 1] + dpLLRRec[ulTemp];
            if (t_dpLLR[ulTemp] > 0) {
                t_dpDecoding[ulTemp] = 0;
            } else {
                t_dpDecoding[ulTemp] = 1;
            }
        }


        ulNum = 0;
        for (unsigned long ulTemp = 0; ulTemp < t_ulCheLength; ulTemp++) {
            if (t_dpCheSta[ulTemp] >= 0) {
                byMiddle = 0;
            } else {
                byMiddle = 1;
            }
            for (unsigned short wTemp = 0; wTemp < wCheDeg; wTemp++) {
                lVarIndex = (long) t_lpCheDis[ulTemp * wCheDeg + wTemp];
                if (lVarIndex >= 0) {
                    byMiddle += (unsigned char) t_dpDecoding[lVarIndex];
                }
            }
            bypSyn[ulTemp] = byMiddle % 2;
            if (bypSyn[ulTemp] != 0) {
                ulNum++;
            }
        }

        ExtractInforBits(t_dpDecoding, ulpReArrangeCol, bypInforbits, t_ulCodeLength, t_ulCheLength);

        if ((ulNum == 0 && check_CRC(bypInforbits, ulInfLength)) || (ulCount == t_ulIterMax)) {
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

    if (ulNum == 0) {
        return 1;
    } else {
        return 0;
    }
}


/**
* 提取生成校验bit的信息位的索引
* @param res 索引存放容器
* @param P_Mats 经过高斯消元的校验矩阵
* @param row P_Mats的行数
* @param columns P_Mats的列数
*/
void getParityPoint(vector<vector<int>> &res, int **P_Mats, const int row, const int columns) {
    const int c_col = columns - row;
    for (int i = 0; i < row; i++) {
        vector<int> tmp;
        for (int j = 0; j < c_col; j++) {
            if (P_Mats[i][j] == 1)
                tmp.push_back(j);
        }
        res.push_back(tmp);
    }
};


/******************校验函数***************/
void getParityPoint(vector<vector<int>> &res, vector<vector<int>> &P_Mats) {
    const int row = P_Mats.size();
    const int columns = P_Mats[0].size();
    const int c_col = columns - row;
    for (int i = 0; i < row; i++) {
        vector<int> tmp;
        for (int j = 0; j < c_col; j++) {
            if (P_Mats[i][j] == 1)
                tmp.push_back(c_col - j - 1);
        }
        res.push_back(tmp);
    }
};


void vec_mod2sum(vector<int> &v1, vector<int> &v2) {
    for (int i = 0; i < v1.size(); i++) {
        if (v1[i] == v2[i]) v1[i] = 0;
        else v1[i] = 1;
    }
}

void Gaussian_Elimination(vector<vector<int>> &matrix) {

    const int row = matrix.size();
    const int column = matrix[0].size();
    const unsigned long matrix_start = column - row;

    for (unsigned long j = matrix_start, i = 0; j < column; j++, i++) {
        if (matrix[i][j] == 0) {
            for (int i1 = 0; i1 < row; i1++) {
                if (matrix[i1][j] == 0 || i1 == i) continue;
                vec_mod2sum(matrix[i], matrix[i1]);
                break;
            }
        }
        for (unsigned long k = 0; k < row; k++) {
            if (matrix[k][j] == 1 && k != i) {
                vec_mod2sum(matrix[k], matrix[i]);
            }
        }
    }
    // 校验结果，查看校验矩阵，是否满秩
    for (int i = 0; i < row; i++) {
        if (matrix[i][matrix_start + i] == 0) {
            std::cout << "the parity matrix is non-singular matrix!!! exit" << std::endl;
            system("pause");
        }
    }
}

void getTime() {
    SYSTEMTIME sys;
    GetLocalTime(&sys);
    printf("%4d/%02d/%02d %02d:%02d:%02d.%03d \n", sys.wYear, sys.wMonth, sys.wDay, sys.wHour, sys.wMinute, sys.wSecond, sys.wMilliseconds);

}

void coutmat(const vector<vector<int>> &matrix) {
    for (auto &i : matrix) {
        for (auto &j : i) {
            std::cout << j << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/*写入txt文件，供matlab读取验证*/
void writeToTxt(int *uncoded_bits, const char *filename, unsigned long length) {
    FILE *fid = fopen(filename, "w");
    for (int i = 0; i < length; i++) {
        fprintf(fid, "%d\n", uncoded_bits[i]);
    }
    fclose(fid);
}

void writeMatToText(int **uncoded_bits, const char *filename, unsigned long row, unsigned long col) {
    FILE *fid = fopen(filename, "w");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col - 1; j++) {
            fprintf(fid, "%d,", uncoded_bits[i][j]);
        }
        fprintf(fid, "%d\n", uncoded_bits[i][col - 1]);
    }
    fclose(fid);
}

void writeMatToText(double *uncoded_bits, const char *filename, unsigned long length) {
    FILE *fid = fopen(filename, "a+");
    for (int i = 0; i < length; i++) {
        fprintf(fid, "%.4f ", uncoded_bits[i]);
    }
    fprintf(fid, "\n");
    fclose(fid);
}