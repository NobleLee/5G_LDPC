//
// Created by ggl on 2018/1/14.
//

#ifndef INC_5G_LDPC_LDPC_HELPER_H
#define INC_5G_LDPC_LDPC_HELPER_H

/**
 * LDPC辅助操作函数
 */
#include<vector>
#include <string>
#include <windows.h>

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

#endif //INC_5G_LDPC_LDPC_HELPER_H