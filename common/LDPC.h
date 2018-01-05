//
// Created by ggl on 2018/1/5.
//

#ifndef INC_5G_LDPC_LDPC_H
#define INC_5G_LDPC_LDPC_H


#include <iostream>
#include <cstdlib>
#include <string>

// 表示编码的bit的类型
template<class EncodeType>
class LDPC {
public:
    // 编码前长度
    const unsigned long infLength;
    // 编码之后输出长度
    const unsigned long codeLength;
    // 采用方案1 还是2
    const int type = 1;

public:
    // 编码输出的内存指针
    EncodeType *coders = nullptr;
    // 信道译码之后的似然比--迭代用到
    double *DECOutputLLR = nullptr;
    // bit的似然比
    double *infoBitLLR = nullptr;
    // 译码的bit
    EncodeType *decodebit = nullptr;

public:

    EncodeType *encoder(EncodeType const *const encoder) {
        encoder(encoder, coders);
    }

    /**
     * @param coder 需要编码的bit指针
     * @param out   输出的编码bit的地址
     * @return     返回要编码结果的指针
     */
    virtual EncodeType *encoder(EncodeType const *const encoder, EncodeType *out);

    int softDecode(double *channelLLR) {
        if (DECOutputLLR == nullptr) {
            DECOutputLLR = new double[codeLength]();
        }
        if (infoBitLLR == nullptr) {
            infoBitLLR = new double[infLength]();
        }
        return decode(channelLLR, DECOutputLLR, infoBitLLR);
    };

    /**
     * @param channelLLR   信道似然比
     * @param DECOutputLLR 经过信道译码之后的似然比-去掉信道信息
     * @param infoBitLLR   信息bit似然比
     * @return 迭代次数
     */
    int decode(double *channelLLR, double *DECOutputLLR, double *infoBitLLR);


    // construct function
    LDPC(unsigned long infLength, unsigned codeLength, int type) : codeLength(codeLength), infLength(infLength) {
        LDPC(infLength, codeLength, type, true);
    }

    /**
    * @param infLength 信息bit长度
    * @param codeLength  编码后输出长度
    * @param type    采用哪个方案
    * @param ismalloc 输出编码时是传递空间，还是返回空间
    */
    LDPC(unsigned long infLength, unsigned codeLength, int type, bool ismalloc) : codeLength(codeLength), infLength(infLength) {
        if (ismalloc) {
            coders = new EncodeType[codeLength]();
        }
    }

    ~LDPC() {
        delete[] coders;
    }

    // getters
    EncodeType *getCoders() const {
        if (coders != 0)
            return coders;
        else {
            std::cout << "please init the encoders space!" << std::endl;
            exit(0);
        }
    }

    unsigned long getInfLength() const {
        return infLength;
    }

    unsigned long getCodeLength() const {
        return codeLength;
    }
};

#endif //INC_5G_LDPC_LDPC_H
