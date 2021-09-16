#include "utils.h"

mat getArrMax(mat *lst, size_t nbEle){
    NULL_PTR_CHK(lst);
    ABOVE_ZERO_CHK(nbEle);
    mat maxVal = 0;
    for (size_t i = 0; i < nbEle; ++i) {
        int val = lst[i];
        if(val > maxVal){
            maxVal = val;
        }
    }
    return maxVal;
}

size_t getArrMinEleIndex(mat* lst, size_t nbEle){
    size_t minIndex = 0;
    for (size_t i = 0; i < nbEle; ++i) {
        if(lst[minIndex] > lst[i]){
            minIndex = i;
        }
    }
    return minIndex;
}

void addToArray(mat* arr, size_t nbMaxEle, mat *idx, mat ele){
    NULL_PTR_CHK(arr);
    ABOVE_ZERO_CHK(nbMaxEle);
    NULL_PTR_CHK(idx);
    assert(*idx < (signed)nbMaxEle);
    arr[*idx] = ele;
    *idx += 1;
}
