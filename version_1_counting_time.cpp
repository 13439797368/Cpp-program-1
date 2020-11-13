//单纯实现矩阵乘法，暴力模拟，生成全部变量为1的矩阵， 计时
#include<iostream>
#include<algorithm>
#include<chrono>
using namespace std;


struct Matrix {
	float* matrix;
	int column;
	int row;
};

Matrix* A = new Matrix;
//14143
Matrix* Build_Ones() {
	A->column = 1414;
	A->row = 1414;
	A->matrix = new float[A->column * A->row];//可以考虑使用malloc可能效率会提升 
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix[i] = 1;
	}
	cout << "build complete" << endl;
	return A;
}

Matrix* Build_Matrix(int row, int col) {
	Matrix* temp = new Matrix;
	temp->column = col;
	temp->row = row;
	temp->matrix = new float[col * row];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			cin >> temp->matrix[i * col + j];
		}
	}
	return temp;
}


float Dot_product(int len, float* vector_1, float* vector_2) {
	float sum = 0;
	for (int i = 0; i < len; i++) {
		sum += vector_1[i] * vector_2[i];
	}
	return sum;
}

void Print_Matrix(Matrix* a) {
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < a->column; j++) {
			cout << a->matrix[i * a->column + j] << " ";
		}
		cout << endl;
	}
	return;
}


Matrix* Matrix_Multiplication(Matrix* a,Matrix* b) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = b->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	//cout << sizeof(temp->matrix) << endl;
	//cout << 4 * sizeof(float) << endl;
	memset(temp->matrix, 0, temp->column * temp->row * sizeof(temp->matrix));
	//Print_Matrix(temp);
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			int flag = 0;
			while(flag < a->column) {
				temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag] * b->matrix[j + flag * b->column];
				flag++;
			}
		}
		cout << "row" << i << " complete" << endl;
	}
	return temp;
}



int main() {
	Build_Ones();
	//Print_Matrix(A);
	using namespace literals;
	auto now = chrono::system_clock::now();
	auto t_c = chrono::system_clock::to_time_t(now - 24h);
	auto start = std::chrono::steady_clock::now();
	Matrix* ans = Matrix_Multiplication(A, A);
	auto end = std::chrono::steady_clock::now();
	//Print_Matrix(ans);
	std::cout << "ms ≈ " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << endl;
	return 0;
}
