#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define SIZE 20 //种群大小
#define HSIZE 10 //每一代个体数量
#define MAXGEN 200 //种群繁殖次数
#define SLEN 200 //染色体数目
#define P_MUTATION 1.0/SLEN //基因突变概率

typedef struct {
    int x[SLEN]; //x:解的自变量，0-1串
    double y; //y=f(x),要优化问题的目标函数值
}Solution;

Solution pop[SIZE]; //解集，父代和子代都存储在这里，前HSIZE个是父代，后面HSIZE个是子代
Solution *current = pop; //当前代，也就是父代
Solution *offspring = pop + HSIZE; //子代解集

void PrintPop(Solution *P) {
    for (int i = 1; i <= HSIZE; i++) {
        printf("%d : %f\n", i, P[i].y);
    }
}

void Decode(int *x, double *xo){ //将0-1串x解码为实数*xo ,假定整数4bits，SLEN-4bits为小数部分长度
    *xo = 0.0;
    int i;
    for (i = 0; i < SLEN; i++) {
        double t = pow(2.0, 3.0-i);
        *xo += x[i] * t;
    }
}

double Function(int *x){
    double xo;
    Decode(x, &xo);  //将0-1串x解码成真正的解xo
    return xo * xo - 3 * xo + 2;     //计算目标函数值
}

//计算一个群体的所有解的目标函数值y ，给出了函数指针，支持个函数的优化
void Evaluate(Solution *P){
    int i;
    for (i = 0; i < HSIZE; i++) {
        P[i].y = Function(P[i].x);
    }
}

//算法初始化：分配两个解集所需的空间，随机生成currentPop中的解，并计算其y值
void Initialize(){
    for(int i = 0; i < HSIZE; i++) { //遍历currentPop.pop中每个解
        for(int j = 0; j < SLEN; j++) { //对每个解的0-1串，随机生成
            current[i].x[j] = rand()%2;
        }
    }
    Evaluate(current);
}

//从父代中选择两个解，通过杂交生成两个子代个体
//父代两个解通过PK选择出来（锦标选择）
void Crossover() { //交叉算子
    int k1, k2, k = 0;
    while (k < HSIZE){ //逐步生成子代，一次两个
        k1 = rand()%HSIZE;
        k2 = rand()%HSIZE;
        if (k1 != k2) {
            int point = rand()%SLEN, i;
            for (i = 0; i < point; i++) {
                offspring[k].x[i] = current[k1].x[i];
                offspring[k+1].x[i] = current[k2].x[i];
            }
            for (i = point; i < SLEN; i++) {
                offspring[k].x[i] = current[k2].x[i];
                offspring[k+1].x[i] = current[k1].x[i];
            }
            k = k+2;
        }
        else {
            continue;
        }
    }
}

//对offspring中的个体进行变异：变异概率为P_MUTATION
//所谓变异就是x[j]的取值 0-1互换： 0 <--> 1
void Mutate() { //变异算子
    for (int i = 0; i < HSIZE; i++) {
        for (int j = 0; j < SLEN; j++) {
            if ((rand()%10000)/10000.0 < P_MUTATION)
            {
                if (offspring[i].x[j] == 0) {
                    offspring[i].x[j] = 1;
                }
                else {
                    offspring[i].x[j] = 0;
                }
            }
        }
    }
}

//从current和offspring中选择下一代个体
//锦标选择，随机选择k个，相互pk，留下最好的放入下一代，依次选择HSIZE个
void Select(int k){ //选择算子 ：采用锦标选择
    int best;
    int besty[HSIZE];
    for (size_t i = 0; i < HSIZE; i++) {
        besty[i] = -1;
    }
    Solution tmp[HSIZE];
    for(int i = 0; i < HSIZE; i++){ //一个一个子代选择
        int choice[k];
        int j = 0;
        while (j < k) {
            choice[j] = rand()%SIZE;
            int flag = 1;
            for (size_t i = 0; i < HSIZE; i++) {
                if (choice[j] == besty[i]) {
                    flag = 0;
                    break;
                }
            }
            if (flag == 0) {
                continue;
            }
            else {
                j++;
            }
        }
        for (int j = 0; j < k; j++) {
            choice[j] = rand()%SIZE;
        }
        int times = 0;
        best = choice[times];
        while (1) {
            int flag = 0, j;
            for (j = 0; j < k; j++) {
                if (pop[best].y <= pop[choice[j]].y) {
                    flag++;
                }
                else {
                    break;
                }
            }
            if (flag == k) {
                break;
            }
            else {
                times++;
                best = choice[times];
            }
        }
        memcpy(&(tmp[i]), &(pop[best]), sizeof(Solution)); //选择出来的解，复制到临时解集中
        besty[i] = best;
    }
    memcpy(current, tmp, sizeof(Solution) * HSIZE);
}

int main(int argc, char const *argv[]) {
    int seed = 991;
    srand(seed);
    Initialize();
    
    for(int gen = 1; gen < MAXGEN; gen++){
        Crossover();
        Mutate();
        Evaluate(offspring);
        Select(2);
    }
    
    PrintPop(current);
    
    return 0;
}
