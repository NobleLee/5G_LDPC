//
// Created by ggl on 2018/1/5.
//

#ifndef INC_5G_LDPC_LDPC_H
#define INC_5G_LDPC_LDPC_H


#include <iostream>
#include <cstdlib>
#include <string>

// 表示编码的bit的类型
class LDPC {
public:
    // 编码前长度
    const unsigned long infLength;
    // 编码之后输出长度
    const unsigned long codeLength;
    // 采用方案1 还是2
    const int type;


public:
    /**
    * @param coder 需要编码的bit指针
    * @param out   输出的编码bit的地址
    * @return     返回要编码结果的指针
    */
    int *encoder(int *infbit, int *out);

    /**
    * @param channelLLR   信道似然比
    * @param DECOutputLLR 经过信道译码之后的似然比-去掉信道信息
    * @param infoBitLLR   信息bit似然比
    * @return 迭代次数
    */
    int decode(double *channelLLR, double *DECOutputLLR, double *infoBitLLR, const int decodeType, const int maxIter);

    int decode(double *channelLLR, double *DECOutputLLR, const int decodeType, const int maxIter);


    // construct function
    LDPC(unsigned long infLength, unsigned codeLength, int type) : codeLength(codeLength), infLength(infLength), type(type) {

    }

};

#endif //INC_5G_LDPC_LDPC_H
