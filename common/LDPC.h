//
// Created by ggl on 2018/1/5.
//

#ifndef INC_5G_LDPC_LDPC_H
#define INC_5G_LDPC_LDPC_H


#include <iostream>
#include <cstdlib>
#include <string>

template<class EncodeType>
class LDPC {
private:
    // 编码前长度
    unsigned long infLength;
    // 编码之后输出长度
    unsigned long codeLength;
    // 编码输出的内存指针
    EncodeType *encoders = 0;
    // 编码方案的名字
    const std::string name;

public:

    int *encoder(EncodeType *encoder);

    double *decode(double *llr);

    // construct functions
    LDPC(unsigned long infLength, unsigned codeLength) {
        LDPC(infLength, codeLength, true);
    }

    LDPC(unsigned long const infLength, const unsigned codeLength, bool ismalloc) {
        this->codeLength = codeLength;
        this->infLength = infLength;
        if (ismalloc) {
            encoders = new EncodeType[infLength]();
        }
    }

    // abstract mothod
    virtual const std::string getName() = 0;

    // getters
    EncodeType *getEncoders() const {
        if (encoders != 0)
            return encoders;
        else {
            std::cout << "please init the encoders space!" << std::endl;
            exit(0);
        }
    }

    const std::string &getName() const {
        return name;
    }

    unsigned long getInfLength() const {
        return infLength;
    }

    unsigned long getCodeLength() const {
        return codeLength;
    }
};

#endif //INC_5G_LDPC_LDPC_H
