//
// Created by ggl on 2018/1/11.
//

#include "LDPC_5G.h"
#include <unordered_map>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include "GF.h"

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

/**
* 初始化速率匹配初始位置
*/
void LDPC_5G::initStartPosition() {
    if (type == 1) {
        switch (rvId) {
            case 0:
                startingPosition = 0;
                break;
            case 1:
                startingPosition = 17 * blockCodeLength / (66 * zLength) * zLength;
                break;
            case 2:
                startingPosition = 33 * blockCodeLength / (66 * zLength) * zLength;
                break;
            case 3:
                startingPosition = 56 * blockCodeLength / (66 * zLength) * zLength;
                break;
            default:
                cout << "rvId error" << endl;
                system("pause");
        }
    } else {
        switch (rvId) {
            case 0:
                startingPosition = 0;
                break;
            case 1:
                startingPosition = 13 * blockCodeLength / (50 * zLength) * zLength;
                break;
            case 2:
                startingPosition = 25 * blockCodeLength / (50 * zLength) * zLength;
                break;
            case 3:
                startingPosition = 50 * blockCodeLength / (50 * zLength) * zLength;
                break;
            default:
                cout << "rvId error" << endl;
                system("pause");
        }
    }

}

void LDPC_5G::tempSpaceInit() {
    CRCTemp = new int[infLength + globalCRCLength];
    bitAddCRC = new int[infLength + globalCRCLength];

    for (int i = 0; i < blockNum; i++) {
        blockBit.push_back(new int[blockLength]);
    }

    for (int i = 0; i < blockNum; i++) {
        afterEncode.push_back(new int[blockCodeLength]);
    }

    int usualLength = codeLength / blockNum;
    int lastLength = usualLength + codeLength % blockNum;
    for (int i = 0; i < blockNum - 1; i++) {
        rateMatchLength.push_back(usualLength);
        rateMatchPosition.push_back(new int[usualLength]);
    }
    rateMatchPosition.push_back(new int[lastLength]);
    rateMatchLength.push_back(lastLength);

    for (int i = 0; i < blockNum; i++) {
        deRateMatchLLR.push_back(new double[blockAfterEncodeLength]);
        bpDecodeLLR.push_back(new double[blockAfterEncodeLength]);
        bpIterLLR.push_back(new double[blockAfterEncodeLength]);
    }

    decodeLLRJudge = new int[blockAfterEncodeLength];
}

/**
* 根据码块长度和Kb来计算：扩展因子、基础矩阵的索引元素
* @param Kb
* @param K_ 码块长度
* @param zLength 扩展因子,会在这个函数中初始化
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
* CRC编码
* @param infbit 输入bit（包含CRC要填充的bit）
* @param infbitLength 有效信息长度
* @param crcType 加CRC的类型
*/
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

/**
* 进行CRC校验，返回校验结果
* @param in 要校验的输入
* @param length in指向地址长度
* @param crcType crc校验的类型
* @return CRC校验结果
*/
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

/**
* 有关CRC校验的参数
*/
void LDPC_5G::crcInit() {
    // 初始化总码块的CRC矩阵
    globalCRCLength = crcLengthList[globalCRCType];
    // 初始化码块分割后的crc
    blockCRCLength = crcLengthList[blockCRCType];
    if (crcList.size() == 0)
        generate_g(crcList);
}

/**
* 获得第二种编码方式的Kb
* @param infLengthCRC
* @return Kb
*/
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
* 初始化矩阵：包括编码要有的校验bit的位置，译码要用的校验关系
* @param I_ls 选择校验矩阵的索引
* @param zLength 扩展因子
*/
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
    /// 提取码字的校验关系,验证一个向量是否是一个码字使用
    getParityMatrixPoint(P_Mats, row * zLength, columns * zLength, checkH);
    cout << "before Gaussian Elimination: ";
    getTime();
    /// 高斯消元，产生生成矩阵
    Gaussian_Elimination(P_Mats, row * zLength, columns * zLength);
    cout << "after Gaussian Elimination: ";
    getTime();
    //coutmat(P_Mats, row * zLength, columns * zLength);
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
* 初始化速率匹配的索引，用于速率匹配和解速率匹配
*/
void LDPC_5G::rateMatchPositionInit() {
    const int blockNullRight = blockLength - 2 * zLength;
    // 分配最长的空间
    int *temp = new int[rateMatchLength[blockNum - 1]];

    for (int i = 0; i < blockNum; i++) {
        /// k 速率匹配结果索引
        /// j 待速率匹配索引
        // bit selection
        int k = 0, j = startingPosition;
        const int blockNullLeft = blockInfBitLength[i] - 2 * zLength;
        while (k < rateMatchLength[i]) {
            j = j % blockCodeLength;
            if (j < blockNullLeft || j >= blockNullRight) {
                temp[k++] = j;
            } else {
                j = blockNullRight - 1;
            }
            j++;
        }
        // bit interleaving
        k = 0;
        int tlength = ceil(1.0 * rateMatchLength[i] / modulationMod);
        for (int jj = 0; jj < tlength; jj++) {
            for (int ii = 0; ii < modulationMod; ii++) {
                if (ii * tlength + jj >= rateMatchLength[i]) continue;
                rateMatchPosition[i][k++] = temp[ii * tlength + jj];
            }
        }
    }

    delete[] temp;

}

/**
* 对bit进行码块分割，同时在每个码块后加入CRC
*/
void LDPC_5G::blockSegmentation(int *bitAddCRC, vector<int *> &blockBit) {
    if (blockNum > 1) {
        int shift = 0;
        for (int i = 0; i < blockNum; i++) {
            if (i == 0) {
                memcpy(blockBit[i], bitAddCRC, (blockInfBitLength[i] - blockCRCLength) * sizeof(int));
                getCRC(blockBit[i], blockInfBitLength[i] - blockCRCLength, blockCRCType);
                for (int j = blockInfBitLength[i]; j < blockLength; j++)
                    blockBit[i][j] = 0;
                shift += blockInfBitLength[i] - blockCRCLength;
            } else {
                memcpy(blockBit[i], bitAddCRC + shift, (blockInfBitLength[i] - blockCRCLength) * sizeof(int));
                getCRC(blockBit[i], blockInfBitLength[i] - blockCRCLength, blockCRCType);
                for (int j = blockInfBitLength[i]; j < blockLength; j++)
                    blockBit[i][j] = 0;
                shift += blockInfBitLength[i] - blockCRCLength;
            }
        }
    } else {
        memcpy(blockBit[0], bitAddCRC, blockInfBitLength[0] * sizeof(int));
        for (int j = blockInfBitLength[0]; j < blockLength; j++)
            blockBit[0][j] = 0;
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
        // 将信息部分信息bit复制到结果
        memcpy(afterEncode[i], blockBit[i] + 2 * zLength, (blockLength - 2 * zLength) * sizeof(int));
        for (k = 0, j = blockLength - 2 * zLength; j < blockCodeLength; j++, k++) {
            temp = 0;
            for (l = 0; l < parityBit[k].size(); l++)
                temp += blockBit[i][parityBit[k][l]];
            afterEncode[i][j] = temp % 2;
        }
    }
}


/**
* 进行解速率匹配
* @param channelInput
* @param deRateMatchLLR
*/
void LDPC_5G::deRateMatch(double *channelInput, vector<double *> &deRateMatchLLR) {

    int k = 0;
    double *rateMatchStartAdr = NULL;
    for (int i = 0; i < blockNum; i++) {
        memset(deRateMatchLLR[i], 0, sizeof(double) * (blockCodeLength + 2 * zLength));
        rateMatchStartAdr = deRateMatchLLR[i] + 2 * zLength;
        for (int j = 0; j < rateMatchLength[i]; j++) {
            //rateMatchStartAdr[rateMatchPosition[i][j]] += channelInput[k++];
            rateMatchStartAdr[rateMatchPosition[i][j]] = channelInput[k++];
        }
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
        blockInfBitLength.push_back(infLengthCRC);
        blockCRCType = globalCRCType;
        blockCodeLength = globalCRCLength;
    } else {
        blockNum = ceil(infLengthCRC * 1.0 / (Kcb - blockCRCLength));

        int usualBlockInfLength = (infLengthCRC + blockNum * blockCRCLength) / blockNum;
        for (int i = 0; i < blockNum - 1; i++)
            blockInfBitLength.push_back(usualBlockInfLength);
        int lastBlockLength = usualBlockInfLength + (infLengthCRC + blockNum * blockCRCLength) % blockNum;
        blockInfBitLength.push_back(lastBlockLength);
    }

    const int I_ls = getZlengthAndI_ls(Kb, blockInfBitLength[0], zLength);
    initStartPosition();
    blockLength = type == 1 ? 22 * zLength : 10 * zLength;
    blockCodeLength = type == 1 ? 66 * zLength : 50 * zLength;
    blockAfterEncodeLength = blockCodeLength + 2 * zLength;

    tempSpaceInit();
    rateMatchPositionInit();
    /*初始化编码矩阵**/
    getGenerateMatrix(I_ls, zLength);
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
    RateMatch(afterEncode, out, 0);
    return out;
}

/**
* 经典BP译码的单次迭代
* @param inputLLR 输入似然比
* @param inputOutLength 输入输出的长度
* @param outputLLR 输出地燃比
* @param edgeVNToVN 节点的校验关系
*/
void singleBPDecode(double *inputLLR, const int inputOutLength, double *outputLLR, const vector<vector<vector<int>>> &edgeVNToVN) {

    if (inputOutLength != edgeVNToVN.size()) {
        cout << "BP matrix error!!" << endl;
        system("pause");
        exit(1);
    }

    for (int i = 0; i < inputOutLength; i++) {
        outputLLR[i] = 0;
        for (int j = 0; j < edgeVNToVN[i].size(); j++) {
            double temp = 1;
            for (int k = 0; k < edgeVNToVN[i][j].size(); k++)
                temp *= tanh(0.5 * inputLLR[edgeVNToVN[i][j][k]]);
            outputLLR[i] += 2.0 / tanh(temp);
        }
    }
}

/**
* 简化的BP译码，MinSum算法
* @param inputLLR
* @param inputOutLength
* @param outputLLR
* @param edgeVNToVN
*/
void singleMinSumDecode(double *inputLLR, const int inputOutLength, double *outputLLR, const vector<vector<vector<int>>> &edgeVNToVN) {
    int count = 0;//计算有多少个-1
    double min = 99999;
    double temp = 0;
    if (inputOutLength != edgeVNToVN.size()) {
        cout << "BP matrix error!!" << endl;
        system("pause");
        exit(1);
    }

    for (int i = 0; i < inputOutLength; i++) {
        outputLLR[i] = 0;
        for (int j = 0; j < edgeVNToVN[i].size(); j++) {
            count = 0;//计算有多少个-1
            min = 9999999;//绝对值最小的数字
            for (int k = 0; k < edgeVNToVN[i][j].size(); k++) {
                temp = inputLLR[edgeVNToVN[i][j][k]];
                count += temp >= 0 ? 0 : 1;
                min = abs(temp) > min ? min : abs(temp);
            }
            outputLLR[i] += min * (count % 2 == 1 ? -1 : 1);
        }
    }
}


bool LDPC_5G::isVaildCode(int *decodeLLRJudge, vector<vector<int>> &checkH) {

    int count = 0;
    for (int i = 0; i < checkH.size(); i++) {
        count = 0;
        for (int &index : checkH[i]) {
            if (decodeLLRJudge[index] == 1) count++;
        }
        if (count % 2 == 1)
            return false;
    }
    return true;
}

/**
* BP译码迭代入口
* @param deRateMatchLLR
* @param bpDecodeLLR
* @param decodeType
* @param maxIter
* @return average iter times
*/
int LDPC_5G::BP_AWGNC(vector<double *> &deRateMatchLLR, vector<double *> &bpDecodeLLR, const int decodeType, const int maxIter) {

    ///choose the function address point according to decodetype
    void (*bp)(double *, const int, double *, const vector<vector<vector<int>>> &);
    bp = decodeType == 0 ? singleBPDecode : singleMinSumDecode;
    int i = 0, j = 0, k = 0;
    int iter = 0;  //compute mean iter times

    for (i = 0; i < blockNum; i++) {
        memcpy(bpDecodeLLR[i], deRateMatchLLR[i], blockAfterEncodeLength * sizeof(double));
        for (j = 0; j < maxIter; j++) {
            /// 将上次输出的结果放到输入
            double *temp = bpDecodeLLR[i];
            bpDecodeLLR[i] = bpIterLLR[i];
            bpIterLLR[i] = temp;
            // complete one iteration
            bp(bpIterLLR[i], blockAfterEncodeLength, bpDecodeLLR[i], edgeVNToVN);
            // get inforBit and crc check
            for (k = 0; k < blockAfterEncodeLength; k++) {
                bpDecodeLLR[i][k] += deRateMatchLLR[i][k];
                decodeLLRJudge[k] = bpDecodeLLR[i][k] >= 0 ? 0 : 1;
            }
            // 校验CRC和是否是一个码字
            if (checkCRC(decodeLLRJudge, blockInfBitLength[i], blockCRCType) && isVaildCode(decodeLLRJudge, checkH))
                break;
        }
        iter += j;
    }
    return ceil(iter * 1.0 / blockNum);
}

/**
* @param channelLLR   信道似然比
* @param DECOutputLLR 经过信道译码之后的似然比-去掉信道信息
* @param infoBitLLR   信息bit似然比
* @return 迭代次数
*/
int LDPC_5G::decode(double *channelLLR, double *DECOutputLLR, double *infoBitLLR, const int decodeType, const int maxIter) {
    deRateMatch(channelLLR, deRateMatchLLR);
    int iter = BP_AWGNC(deRateMatchLLR, bpDecodeLLR, decodeType, maxIter);

    int index = 0;
    for (int i = 0; i < blockNum; i++) {
        memcpy(infoBitLLR + index, bpDecodeLLR[i], blockInfBitLength[i] - blockCRCLength);
        index = blockInfBitLength[i];
    }

    for (int i = 0; i < blockNum; i++) {
        for (int j = 0; j < blockAfterEncodeLength; j++) {
            bpDecodeLLR[i][j] -= deRateMatchLLR[i][j];
        }
    }

    RateMatch(bpDecodeLLR, DECOutputLLR, 2 * zLength);

    return iter;
}

int LDPC_5G::decode(double *channelLLR, double *DECOutputLLR, const int decodeType, const int maxIter) {
    deRateMatch(channelLLR, deRateMatchLLR);
    int iter = BP_AWGNC(deRateMatchLLR, bpDecodeLLR, decodeType, maxIter);

    for (int i = 0; i < blockNum; i++) {
        for (int j = 0; j < blockAfterEncodeLength; j++) {
            bpDecodeLLR[i][j] -= deRateMatchLLR[i][j];
        }
    }

    RateMatch(bpDecodeLLR, DECOutputLLR, 2 * zLength);
    return iter;
}

int LDPC_5G::decode(double *channelLLR, int *outBit, const int decodeType, const int maxIter) {
    deRateMatch(channelLLR, deRateMatchLLR);




    /*for (int j = 0; j < blockNum; j++) {
    int count = 0;
    for (int i = 0; i < blockCodeLength; i++) {
    if (afterEncode[j][i] - deRateMatchLLR[j][512 + i] != 0)
    count++;
    }
    cout << "block " << j << " don't match point count " << count << endl;
    }*/

    int iter = BP_AWGNC(deRateMatchLLR, bpDecodeLLR, decodeType, maxIter);
    int index = 0;
    for (int i = 0; i < blockInfBitLength[0] - blockCRCLength; i++) {
        outBit[index++] = bpDecodeLLR[0][i] >= 0 ? 0 : 1;
    }

    /*int index = 0;
    for (int i = 0; i < blockInfBitLength[0] - blockCRCLength; i++) {
    outBit[index++] = bpDecodeLLR[0][i] >= 0 ? 0 : 1;
    }
    for (int i = 0; i < blockInfBitLength[0] - blockCRCLength-24; i++) {
    outBit[index++] = bpDecodeLLR[1][i] >= 0 ? 0 : 1;
    }*/
    return iter;
}
