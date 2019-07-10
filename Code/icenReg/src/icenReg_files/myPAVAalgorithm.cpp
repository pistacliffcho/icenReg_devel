//
//  myPAVAalgorithm.cpp
//  
//
//  Created by Cliff Anderson Bergman on 8/16/15.
//
//

void weighted_pool(double *y, double *w, int start, int stop){
    double sum_weights = 0, sum_weightedResponse = 0;
    for(int i = start; i <= stop; i++){
        sum_weightedResponse += y[i] * w[i];
        sum_weights += w[i];
    }
    
    sum_weightedResponse = sum_weightedResponse/sum_weights;
    for(int i = start; i <= stop; i++){y[i] = sum_weightedResponse;}
}

void weighted_pava(double *y, double *w, int *numberParameters){
    int numberPools, ind, k;
    int n = (*numberParameters);
    
    if(n <= 1) return;
    n--;
    do{
        ind = 0;
        numberPools = 0;
        while(ind < n){
            k = ind;
            while(k < n && y[k] >= y[k+1]) k++;
            if(y[ind] != y[k]){
                weighted_pool(y, w, ind, k);
                numberPools++;
            }
            ind = k + 1;
        }
    }while(numberPools > 0);
}


void pava(double *y, double *w, int *np){
    int n = (*np);
    double maxAbsW = 0;
    for(int i = 0; i < n; i++){
        maxAbsW = max(maxAbsW, abs(w[i]));
    }
    if(maxAbsW == 0){return;}
    weighted_pava(y, w, np);
}
