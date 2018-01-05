//
// Created by ggl on 2018/1/11.
//

#include "LDPC_5G.h"
#include <cstdarg>
#include <unordered_map>
#include <math.h>

using std::unordered_map;

// 不同类型的crc长度列表
int crcLengthList[5] = {24, 24, 24, 16, 11};

// 存放不同码长的map
unordered_map<int, int> *initZLsMap() {
    static unordered_map<int, int> *map = new unordered_map<int, int>();
    static bool inited = false;
    if (!inited) {
        map->insert(std::make_pair<int, int>(2, 0));
        map->insert(std::make_pair<int, int>(4, 0));
        map->insert(std::make_pair<int, int>(8, 0));
        map->insert(std::make_pair<int, int>(16, 0));
        map->insert(std::make_pair<int, int>(32, 0));
        map->insert(std::make_pair<int, int>(64, 0));
        map->insert(std::make_pair<int, int>(128, 0));
        map->insert(std::make_pair<int, int>(256, 0));

        map->insert(std::make_pair<int, int>(3, 1));
        map->insert(std::make_pair<int, int>(6, 1));
        map->insert(std::make_pair<int, int>(12, 1));
        map->insert(std::make_pair<int, int>(24, 1));
        map->insert(std::make_pair<int, int>(48, 1));
        map->insert(std::make_pair<int, int>(96, 1));
        map->insert(std::make_pair<int, int>(192, 1));
        map->insert(std::make_pair<int, int>(384, 1));

        map->insert(std::make_pair<int, int>(5, 2));
        map->insert(std::make_pair<int, int>(10, 2));
        map->insert(std::make_pair<int, int>(20, 2));
        map->insert(std::make_pair<int, int>(40, 2));
        map->insert(std::make_pair<int, int>(80, 2));
        map->insert(std::make_pair<int, int>(160, 2));
        map->insert(std::make_pair<int, int>(320, 2));

        map->insert(std::make_pair<int, int>(7, 3));
        map->insert(std::make_pair<int, int>(14, 3));
        map->insert(std::make_pair<int, int>(28, 3));
        map->insert(std::make_pair<int, int>(56, 3));
        map->insert(std::make_pair<int, int>(112, 3));
        map->insert(std::make_pair<int, int>(224, 3));

        map->insert(std::make_pair<int, int>(9, 4));
        map->insert(std::make_pair<int, int>(18, 4));
        map->insert(std::make_pair<int, int>(36, 4));
        map->insert(std::make_pair<int, int>(72, 4));
        map->insert(std::make_pair<int, int>(144, 4));
        map->insert(std::make_pair<int, int>(288, 4));

        map->insert(std::make_pair<int, int>(11, 5));
        map->insert(std::make_pair<int, int>(22, 5));
        map->insert(std::make_pair<int, int>(44, 5));
        map->insert(std::make_pair<int, int>(88, 5));
        map->insert(std::make_pair<int, int>(176, 5));
        map->insert(std::make_pair<int, int>(352, 5));

        map->insert(std::make_pair<int, int>(13, 6));
        map->insert(std::make_pair<int, int>(26, 6));
        map->insert(std::make_pair<int, int>(52, 6));
        map->insert(std::make_pair<int, int>(104, 6));
        map->insert(std::make_pair<int, int>(208, 6));

        map->insert(std::make_pair<int, int>(15, 7));
        map->insert(std::make_pair<int, int>(30, 7));
        map->insert(std::make_pair<int, int>(60, 7));
        map->insert(std::make_pair<int, int>(120, 7));
        map->insert(std::make_pair<int, int>(240, 7));
    }

    return map;
}

/**
 * @tparam EncodeType
 * @param adress  需要填充1的起始地址
 * @param count   填充1的位置的数目
 * @param ...     填充1的位置
 * @return
 */
template<class EncodeType>
template<class T>
int fillAddress(EncodeType *addr, T fillThing, int count, ...) {
    va_list args;
    va_start(args, count);
    for (int i = 0; i < count; ++i) {
        int index = va_arg(args, int);
        addr[index] = fillThing;
    }
    va_end(args);
    return 0;
}

/**
 * @tparam EncodeType
 * @param type CRC校验选择的类型
 * @return CRC校验的数组
 */
template<class EncodeType>
EncodeType *generate_g(int type) {
    static EncodeType *crcList[5];
    static bool inited = false;
    if (!inited) {
        EncodeType *A_24 = new EncodeType[crcLengthList[0] + 1]();
        EncodeType *B_24 = new EncodeType[crcLengthList[1] + 1]();
        EncodeType *C_24 = new EncodeType[crcLengthList[2] + 1]();
        EncodeType *O_16 = new EncodeType[crcLengthList[3] + 1]();
        EncodeType *O_11 = new EncodeType[crcLengthList[4] + 1]();
        crcList[0] = A_24;
        crcList[1] = B_24;
        crcList[2] = C_24;
        crcList[3] = O_16;
        crcList[4] = O_11;
        //CRC_A初始化
        fillAddress(A_24, 1, 14, 0, 1, 3, 4, 5, 6, 7, 10, 11, 14, 17, 18, 23, 24);
        //CRC_B初始化
        fillAddress(B_24, 1, 6, 0, 1, 5, 6, 23, 24);
        //CRC_C初始化
        fillAddress(C_24, 1, 13, 0, 1, 2, 4, 8, 12, 13, 15, 17, 20, 21, 23, 24);
        //CRC_16
        fillAddress(O_16, 1, 4, 0, 5, 12, 16);
        //CRC_11
        fillAddress(O_11, 1, 5, 0, 5, 9, 10, 11);
    }
    return crcList[type];
}

template<class EncodeType>
void LDPC_5G::getCRC(EncodeType encoder, int crcType) {

    int g_length = crcLengthList[crcType];
    EncodeType *g = generate_g(crcType);
    int i, flag = infLength - g_length;

    EncodeType *temp = (EncodeType *) malloc(sizeof(EncodeType) * g_length);
    for (i = 0; i < g_length; i++) {
        temp[i] = encoder[flag + i]; //取出值
    }
    int j = 0;
    while (flag >= 0) {
        if (temp[g_length - 1] == 1)//最高位为1
        {
            if (flag != 0) {
                for (i = g_length - 2; i >= 0; i--) {
                    temp[i + 1] = (temp[i] + g[i]) % 2;//模2加再向高位移一位
                }
                temp[0] = encoder[flag - 1];
            } else {
                for (i = g_length - 1; i >= 0; i--) {
                    temp[i] = (temp[i] + g[i]) % 2;
                }
            }
        } else {
            if (flag != 0) {
                for (i = g_length - 2; i >= 0; i--) {
                    temp[i + 1] = temp[i];
                }
                temp[0] = encoder[flag - 1];
            }
        }
        flag--;
    }
    for (i = 0; i < (g_length - 1); i++) {
        encoder[i] = temp[i];//校验位在最左边
    }
    free(temp);
}

void LDPC_5G::crcInit() {
    // 初始化总码块的CRC矩阵
    codeCRC = generate_g(globalCRCType);
    codeCRCLength = crcLengthList[globalCRCType];
    // 初始化码块分割后的crc
    blockCRC = generate_g(blockCRCType);
    blockCRCLength = crcLengthList[blockCRCType];
}

inline int getKbGraph2(const int infLengthCRC) {
    if (infLengthCRC > 640)
        return 10;
    else if (infLengthCRC > 560)
        return 9;
    else if (infLengthCRC > 192)
        return 8;
    return 6;
}

/**
 * 根据码块长度和Kb来计算：扩展因子、基础矩阵的索引元素
 * @param Kb
 * @param K_ 码块长度
 * @return 返回索引
 */
int LDPC_5G::getZlengthAndI_ls(const int Kb, const int K_) {
    unordered_map<int, int> *map = initZLsMap();
    int temp = ceil(K_ * 1.0 / Kb);
    while (temp < 385) {
        auto zSet = map->find(temp);
        if (zSet != map->end()) {
            zLength = zSet->first;
            return zSet->second;
        }
        temp++;
    }
}

/**
 * 初始化数据，包括码块数目，扩展因子等
 */
void LDPC_5G::init() {
    crcInit();
    int infLengthCRC = infLength + codeCRCLength;

    //默认方案1配置
    int Kcb = type == 1 ? 8448 : 3840;
    int Kb = type == 1 ? 22 : getKbGraph2(infLengthCRC);

    if (infLengthCRC <= Kcb) {
        // 不进行码块分割时
        blockNum = 1;
        blockLength = infLengthCRC;
    } else {

    }
    const int I_ls = getZlengthAndI_ls(Kb, blockLength);


}


