//another version of O(n^3)
#include<iostream>
#include<chrono>
using namespace std;

struct mat {
    int n, m;
    double** data;
};

mat* A = new mat;

int mul(mat* c, const mat* a, const mat* b) {
    int i, j, k;
    if (a->m != b->n)
        return 0;
    c->n = a->n;
    c->m = b->m;
    c->data = new double* [1414];
    for (int i = 0; i < 1414; i++) {
        c->data[i] = new double[1414];
    }
    for (i = 0; i < c->n; i++)
        for (j = 0; j < c->m; j++)
            for (c->data[i][j] = k = 0; k < a->m; k++)
                c->data[i][j] += a->data[i][k] * b->data[k][j];
    return 1;
}

mat* Build_Ones() {
    A->data = new double* [1414];
    for (int i = 0; i < 1414; i++) {
        A->data[i] = new double[1414];
    }
    A->m = A->n = 1414;
    for (int i = 0; i < 1414; i++) {
        for (int j = 0; j < 1414; j++) {
            A->data[i][j] = 1;
        }
    }
    cout << "build complete" << endl;
    return A;
}

int main() {
    Build_Ones();
    mat* c = new mat;
    using namespace literals;
    auto now = chrono::system_clock::now();
    auto t_c = chrono::system_clock::to_time_t(now - 24h);
    auto start = std::chrono::steady_clock::now();
    mul(c, A, A);
    auto end = std::chrono::steady_clock::now();
    std::cout << "ms ?" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << endl;
    return 0;
}
