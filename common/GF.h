//
// Created by ggl on 2018/1/21.
//

#ifndef INC_5G_LDPC_GF_H
#define INC_5G_LDPC_GF_H

#include <iostream>

template<class T, class U>
void Diffs(T *a, U *b, int l, char *str) {
    int errorNum = 0;
    for (int i = 0; i < l; i++) {
        if (a[i] - b[i] != 0) {
            errorNum++;
            std::cout << i << " ";
        }
    }
    std::cout << str << " errorNum: " << errorNum << std::endl;
};


template<class T>
void coutmat(T *matrix, unsigned long row, unsigned long column) {
    cout << "ROW " << row << "  column" << column << endl;
    for (unsigned long i = 0; i < row; i++) {
        std::cout << i << "\t";
        for (unsigned long j = 0; j < column; j++) {
            if (matrix[i][j] == 1)
                cout << j << ",";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#endif //INC_5G_LDPC_GF_H
