#ifndef ARR_H
#define ARR_H

#include "utils.h"

mat getArrMax(mat *lst, size_t nbEle);
size_t getArrMinEleIndex(mat* lst, size_t nbEle);
void addToArray(mat* arr, size_t nbMaxEle, mat *idx, mat ele);

#endif // ARR_H
