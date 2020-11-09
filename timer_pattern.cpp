#include<iostream>
#include<chrono>
using namespace std;


int main() {
    const int len = 200000000;
    int* num = new int[len];
    for (int i = 0; i < len; i++) {
        num[i] = 1;
    }
    using namespace literals;
    auto now = chrono::system_clock::now();
    auto t_c = chrono::system_clock::to_time_t(now - 24h);
    auto start = std::chrono::steady_clock::now();
    long long int s = 0;
    for (int i = 0; i < 200000000; i++) {
        s += num[i] * num[i];
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "ms ˜ " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << endl;
    cout << s << endl;
    delete[] num;
    return 0;
}
