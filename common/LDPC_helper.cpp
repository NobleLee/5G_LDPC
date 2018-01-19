//
// Created by ggl on 2018/1/14.
//

#include <iostream>
#include <cstdlib>
#include<fstream>
#include "LDPC_helper.h"


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
            for (int i1 = 0; i1 < row; i1++) {
                if (matrix[i1][j] == 0 || i1 == i) continue;
                vec_mod2sum(matrix[i], matrix[i1], column);
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
                tmp.push_back(c_col - j - 1);
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

void getTime(){
    SYSTEMTIME sys;
    GetLocalTime(&sys);
    printf("%4d/%02d/%02d %02d:%02d:%02d.%03d \n", sys.wYear, sys.wMonth, sys.wDay, sys.wHour, sys.wMinute, sys.wSecond, sys.wMilliseconds);

}

void coutmat(const vector<vector<int>> &matrix) {
    for (auto &i:matrix) {
        for (auto &j:i) {
            std::cout << j << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}