//
// Created by ggl on 2018/1/14.
//
#define _CRT_SECURE_NO_WARNINGS
#ifndef INC_5G_LDPC_LDPC_HELPER_H
#define INC_5G_LDPC_LDPC_HELPER_H

/**
* LDPC辅助操作函数
*/
#include<vector>
#include <string>
#include <windows.h>
#include <fstream>

using std::vector;
using std::string;

/**
* 对矩阵进行高斯消元，矩阵右侧变换为单位阵
* @param matrix
*/
void Gaussian_Elimination(int **matrix, const int row, const int column);

/**
* 将字符串按照某一字符进行分割
* @param s 要分割的字符串
* @param v 保存分割结果的容器
* @param c 分割字符
*/
void SplitString(const string &s, vector<string> &v, const string &c);

/**
* 获取变量节点与校验节点的关系，将结果存放在req中，结果的第j个元素，是通过校验节点与第j个变量节点连接的变量节点的索引
* 用途:BP译码迭代用
* @param P_mats 扩展之后的校验矩阵
* @param row 校验矩阵的行数
* @param column 校验矩阵的列数
* @param req 存放结果
*/
void getEdgeFrom_VNandCN(int **P_mats, const unsigned long row, const unsigned long column, vector<vector<vector<int>>> &req);

/**
* 读取文件名为filename的校验矩阵文件，并将结果放到parityMats中
* @param filename 文件名
* @param parityMats 校验矩阵
*/
void readMatrixFromFile(string &filename, vector<vector<int>> &parityMats);

/**
* 将H_base中的矩阵按照扩展因子扩展，结果存放在P_Mats
* @param P_Mats 存放扩展结果
* @param H_base 存放基础矩阵
* @param row 基础矩阵的行
* @param columns 基础矩阵的列
* @param zLength 扩展因子
*/
void expendParityMatrix(int **P_Mats, int **H_base, const int row, const int columns, const int zLength);

void SAM_expendParityMatrix(unsigned short *LDPC_Matrix, long *BaseMatrix, const int row, const int columns, const int zLength);


/***
 * 高斯消元
 * @param t_LDPC_Matrix
 * @param t_P_Matrix
 * @param t_ulpReArrangeCol
 * @param t_ulCheLength
 * @param t_ulCodeLength
 * @return
 */
bool SAM_Gaussian_Elimination(unsigned short *t_LDPC_Matrix, unsigned short *t_P_Matrix, unsigned long *t_ulpReArrangeCol, unsigned long t_ulCheLength, unsigned long t_ulCodeLength);

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
                           unsigned long *ulpReArrangeCol, unsigned long t_ulCodeLength, unsigned long t_ulCheLength);

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
                  unsigned long t_ulCheLength, unsigned short t_byVarDeg, unsigned short t_byCheDeg, unsigned long t_ulIterMax, unsigned short *bypInforbits,
                  unsigned long *ulpReArrangeCol);

/**
* 提取生成校验bit的信息位的索引
* @param res 索引存放容器
* @param P_Mats 经过高斯消元的校验矩阵
* @param row P_Mats的行数
* @param columns P_Mats的列数
*/
void getParityPoint(vector<vector<int>> &res, int **P_Mats, const int row, const int columns);

void getParityMatrixPoint(int **P_mats, const unsigned long row, const unsigned long column, vector<vector<int>> &req);

void getTime();

void getParityPoint(vector<vector<int>> &res, vector<vector<int>> &P_Mats);

void Gaussian_Elimination(vector<vector<int>> &matrix);

void coutmat(const vector<vector<int>> &matrix);

/*写入txt文件，供matlab读取验证*/
void writeToTxt(int *uncoded_bits, const char *filename, unsigned long length);

/* 将矩阵写到txt文件中 */
void writeMatToText(int **uncoded_bits, const char *filename, unsigned long row, unsigned long col);

void writeMatToText(double *uncoded_bits, const char *filename, unsigned long length);

#endif //INC_5G_LDPC_LDPC_HELPER_H