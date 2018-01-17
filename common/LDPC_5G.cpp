//
// Created by ggl on 2018/1/11.
//

#include "LDPC_5G.h"
#include "LDPC_helper.h"
#include <unordered_map>
#include <math.h>

#include <cstdarg>

using namespace std;

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
 * @param adress  需要填充fillThing的起始地址
 * @param count   填充fillThing的位置的数目
 * @param ...     填充fillThing的位置
 * @return
 */
template<class EncodeType, class T>
int fillAddress(EncodeType *addr, T fillThing, T sub, int count, ...) {
    va_list args;
    va_start(args, count);
    for (int i = 0; i < count; ++i) {
        int index = va_arg(args, int);
        addr[sub - index] = fillThing;
    }
    va_end(args);
    return 0;
}

/**
 * 初始化CRC多项式
 * @param crcList
 */
void generate_g(vector<int *> &crcList) {
    int *A_24 = new int[crcLengthList[0] + 1]();
    int *B_24 = new int[crcLengthList[1] + 1]();
    int *C_24 = new int[crcLengthList[2] + 1]();
    int *O_16 = new int[crcLengthList[3] + 1]();
    int *O_11 = new int[crcLengthList[4] + 1]();
    //CRC_A初始化
    fillAddress(A_24, 1, 24, 14, 0, 1, 3, 4, 5, 6, 7, 10, 11, 14, 17, 18, 23, 24);
    //CRC_B初始化
    fillAddress(B_24, 1, 24, 6, 0, 1, 5, 6, 23, 24);
    //CRC_C初始化
    fillAddress(C_24, 1, 24, 13, 0, 1, 2, 4, 8, 12, 13, 15, 17, 20, 21, 23, 24);
    //CRC_16
    fillAddress(O_16, 1, 16, 4, 0, 5, 12, 16);
    //CRC_11
    fillAddress(O_11, 1, 11, 5, 0, 5, 9, 10, 11);
    crcList.push_back(A_24);
    crcList.push_back(B_24);
    crcList.push_back(C_24);
    crcList.push_back(O_16);
    crcList.push_back(O_11);
}

void LDPC_5G::tempSpaceInit() {
    CRCTemp = new int[infLength + globalCRCLength];
    bitAddCRC = new int[infLength + globalCRCLength];

    for (int i = 0; i < blockNum; i++) {
        blockBit.push_back(new int[blockLength]);
    }


    for (int i = 0; i < blockNum; i++) {
        afterEncode[i] = new int[blockCodeLength];
    }
}

void LDPC_5G::getCRC(int *infbit, const int infbitLength, const int crcType) {
    // 根据CRC的类型，取出CRC长度
    int g_length = crcLengthList[crcType];
    int *g = crcList[crcType];
    //int *temp = new int[infLength + g_length]();
    memcpy(CRCTemp, infbit, infbitLength * sizeof(int));
    memset(CRCTemp + infbitLength, 0, g_length * sizeof(int));
    int cur = 0, i = 0;//比较的游标，指示目前在哪个位置
    while (cur < infbitLength) {
        //i = cur;
        if (CRCTemp[cur] == 1) {
            for (i = cur; i < cur + g_length; i++) {
                CRCTemp[i] = (CRCTemp[i] + g[i - cur]) % 2;
            }
        }
        cur++;
    }
    memcpy(infbit + infbitLength, CRCTemp + infbitLength, g_length * sizeof(int));
}

bool LDPC_5G::checkCRC(int *in, const int length, const int crcType) {
    int g_length = crcLengthList[crcType];
    int *g = crcList[crcType];
    memcpy(CRCTemp, in, length * sizeof(int));
    int cur = 0, i = 0, sum = 0;//比较的游标，指示目前在哪个位置
    while (cur < length - g_length + 1) {
        //i = cur;
        if (CRCTemp[cur] == 1) {
            for (i = cur; i < cur + g_length; i++) {
                CRCTemp[i] = (CRCTemp[i] + g[i - cur]) % 2;
            }
        }
        cur++;

    }
    for (; cur < length; cur++) {
        sum += CRCTemp[cur];
    }
    return sum == 0 ? true : false;
}

void LDPC_5G::crcInit() {
    // 初始化总码块的CRC矩阵
    globalCRCLength = crcLengthList[globalCRCType];
    // 初始化码块分割后的crc
    blockCRCLength = crcLengthList[blockCRCType];
    if (crcList.size() == 0)
        generate_g(crcList);
}

/**
 * 对bit进行码块分割，同时在每个码块后加入CRC
 */
void LDPC_5G::blockSegmentation(int *bitAddCRC, vector<int *> &blockBit) {
    for (int i = 0; i < blockNum; i++) {
        if (i == 0) {
            memcpy(blockBit[i], bitAddCRC, blockInfBitLength[i] * sizeof(int));
            getCRC(blockBit[i], blockInfBitLength[i], blockCRCType);
            for (int j = blockInfBitLength[i] + blockCRCLength; j < blockLength; j++)
                blockBit[i][j] = 0;
        } else {
            memcpy(blockBit[i], bitAddCRC + i * blockInfBitLength[i - 1], blockInfBitLength[i] * sizeof(int));
            getCRC(blockBit[i], blockInfBitLength[i], blockCRCType);
            for (int j = blockInfBitLength[i] + blockCRCLength; j < blockLength; j++)
                blockBit[i][j] = 0;
        }
    }
}

/**
 * 利用经过高斯消元的校验矩阵的校验关系进行编码
 * @param blockBit 待编码码块
 * @param afterEncode 编码之后的码块
 */
void LDPC_5G::LDPC_Fast_Encode(vector<int *> &blockBit, vector<int *> &afterEncode) {
    int temp = 0, i = 0, j = 0, k = 0, l = 0;
    for (i = 0; i < blockNum; i++) {
        memcpy(afterEncode[i], blockBit[i] + 2 * zLength, (blockLength - 2 * zLength) * sizeof(int));
        for (k = 0, j = blockLength - 2 * zLength; j < blockCodeLength; j++, k++) {
            temp = 0;
            for (l = 0; l < parityBit[k].size(); l++)
                temp += blockBit[i][parityBit[k][l]];
            afterEncode[i][j] = temp % 2;
        }
    }
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


void LDPC_5G::getGenerateMatrix(const int I_ls, const int zLength) {
    // 存放提案中所有的矩阵
    vector<vector<int>> parityMats;
    string filename = type == 1 ? "LDPC_P1_38.212.txt" : "LDPC_P2_38.212.txt";
    readMatrixFromFile(filename, parityMats);

    const int row = type == 1 ? 46 : 42;
    const int columns = type == 1 ? 68 : 52;

    int **H_base = new int *[row]();
    for (int i = 0; i < row; i++) {
        H_base[i] = new int[columns];
        memset(H_base[i], -1, columns * sizeof(int));
    }
    // 提取对应的校验矩阵
    int loc_x = 0, loc_y = 0;
    for (int i = 0; i < parityMats.size(); i++) {
        loc_x = parityMats[i][0];
        loc_y = parityMats[i][1];
        H_base[loc_x][loc_y] = parityMats[i][I_ls + 2] % zLength;
    }
    //产生校验矩阵
    int **P_Mats = new int *[row * zLength]; //校验矩阵
    for (int i = 0; i < row * zLength; i++) {
        P_Mats[i] = new int[columns * zLength]();
    }
    expendParityMatrix(P_Mats, H_base, row, columns, zLength);
    /// 提取边连接关系-译码使用
    getEdgeFrom_VNandCN(P_Mats, row * zLength, columns * zLength, edgeVNToVN);

    /// 高斯消元，产生生成矩阵
    Gaussian_Elimination(P_Mats, row * zLength, columns * zLength);

    /// 产生生成矩阵校验关系-编码使用
    getParityPoint(parityBit, P_Mats, row * zLength, columns * zLength);

    for (int i = 0; i < row; i++) {
        free(H_base[i]);
    }
    for (int i = 0; i < row * zLength; i++) {
        free(P_Mats[i]);
    }
}

/**
 * 根据码块长度和Kb来计算：扩展因子、基础矩阵的索引元素
 * @param Kb
 * @param K_ 码块长度
 * @param zLength 扩展因子
 * @return 返回基础矩阵的索引
 */
int LDPC_5G::getZlengthAndI_ls(const int Kb, const int K_, int &zLength) {
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
 * 初始化数据，包括码块数目，边的连接关系，生成和校验矩阵，扩展因子等
 */
void LDPC_5G::init() {
    crcInit();
    int infLengthCRC = infLength + globalCRCLength;

    //默认方案1配置
    int Kcb = type == 1 ? 8448 : 3840;  //单个码块最大长度
    int Kb = type == 1 ? 22 : getKbGraph2(infLengthCRC);

    if (infLengthCRC <= Kcb) {
        // 不进行码块分割时
        blockNum = 1;
        blockLength = infLengthCRC;
    } else {
        blockNum = ceil(infLengthCRC * 1.0 / (Kb - blockCRCLength));

        int usualBlockInfLength = infLengthCRC / blockNum;
        for (int i = 0; i < blockNum - 1; i++)
            blockInfBitLength.push_back(usualBlockInfLength);
        int lastBlockLength = usualBlockInfLength + infLengthCRC % blockNum;

        blockInfBitLength.push_back(lastBlockLength);
    }

    const int I_ls = getZlengthAndI_ls(Kb, blockLength, zLength);
    blockLength = type == 1 ? 22 * zLength : 10 * zLength;
    blockCodeLength = type == 1 ? 66 * zLength : 50 * zLength;
    /*初始化编码矩阵**/
    getGenerateMatrix(I_ls, zLength);
    tempSpaceInit();
}


int *LDPC_5G::encoder(int *in, int *out) {
    memcpy(bitAddCRC, in, infLength * sizeof(int));
    /// 加入全局CRC
    getCRC(bitAddCRC, infLength, globalCRCType);
    /// 码块分割
    blockSegmentation(bitAddCRC, blockBit);
    /// 快速编码
    LDPC_Fast_Encode(blockBit, afterEncode);
    /// 速率匹配
}



