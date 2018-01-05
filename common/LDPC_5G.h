//
// Created by ggl on 2018/1/11.
//

#ifndef INC_5G_LDPC_LDPC_5G_H
#define INC_5G_LDPC_LDPC_5G_H

#include "LDPC.h"
#include <vector>

using namespace std;

template<class EncodeType>

class LDPC_5G : public LDPC<EncodeType> {
private:

    const string name = "R1-38.212";
    // CRC相关的
    int globalCRCType = 0;
    int blockCRCType = 1;
    EncodeType *codeCRC;
    int codeCRCLength;
    EncodeType *blockCRC;
    int blockCRCLength;

    int blockNum;
    int blockLength;

    //扩展因子
    int zLength;

    EncodeType *encoderCRC;

    vector<vector<int> *> *generateMatrix;

    void init();

    void crcInit();

    int getZlengthAndI_ls(const int Kb, const int K_);

public:
    LDPC_5G(unsigned long infLength, unsigned long codeLength, int type) : LDPC<EncodeType>(infLength, codeLength, type) {
        init();
    };


    EncodeType *encoder(EncodeType const *const encoder, EncodeType *out) {
        getCRC(encoder, globalCRCType);
        return nullptr;
    }

    void getCRC(EncodeType encoder, int crcType);

    /**
     * @param indexSet -i_ls
     * @param zlength  -矩阵扩展因子
     * @return
     */
    vector<vector<int> *> *getGenerateMatrix(int indexSet, int zlength) {
        vector<vector<int> *> *generateMatrix = new vector<vector<int> *>;

        return generateMatrix;

    }
};


#endif //INC_5G_LDPC_LDPC_5G_H

