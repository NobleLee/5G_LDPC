//
// Created by ggl on 2018/1/11.
//

#ifndef INC_5G_LDPC_LDPC_5G_H
#define INC_5G_LDPC_LDPC_5G_H

#include "LDPC.h"
#include <vector>

using namespace std;

class LDPC_5G : public LDPC {
private:
    const string name = "R1-38.212";
    // CRC相关的
    int globalCRCType = 0;
    int globalCRCLength;
    // 提案中固定为24
    int blockCRCType = 1;
    int blockCRCLength;

    int blockNum;
    int blockLength;
    /// 存储变量节点之间的关系，译码要用
    vector<vector<vector<int>>> edgeVNToVN;
    /// 存储信息bit的校验关系，编码使用
    vector<vector<int>> parityBit;
    /// 存放crc多项式数组
    static vector<int *> crcList;
    /// 临时空间
    int *CRCTemp;
    int *bitAddCRC;

private:
    void init();

    void crcInit();

    int getZlengthAndI_ls(const int Kb, const int K_, int &zLength);

    void tempSpaceInit();

public:
    LDPC_5G(unsigned long infLength, unsigned long codeLength, int type) : LDPC(infLength, codeLength, type) {
        init();
    };

    int *encoder(int *in, int *out);

    void getCRC(int *infbit, const int infbitLength, const int crcType);

    bool checkCRC(int *in, const int length, const int crcType);

    /**
     * @param indexSet -i_ls
     * @param zlength  -矩阵扩展因子
     * @return
     */
    void getGenerateMatrix(int indexSet, int zlength);

    ~LDPC_5G() {
        delete[] CRCTemp;
        delete[] bitAddCRC;
    }
};


#endif //INC_5G_LDPC_LDPC_5G_H

