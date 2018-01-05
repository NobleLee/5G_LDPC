//
// Created by ggl on 2018/1/11.
//
#include "LDPC_5G.h"
#include<iostream>

using namespace std;

int main() {
    LDPC_5G<int> *ldpc = new LDPC_5G<int>(512, 1024, 1);
    std::cout << ldpc->getCodeLength() << "   " << ldpc->getInfLength() << "  " << std::endl;
    system("pause");
}