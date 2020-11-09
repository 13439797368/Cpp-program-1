//单纯实现矩阵乘法，暴力模拟，使用cin cout，没有错判 

#include<iostream>
#include<algorithm>
using namespace std;

struct Matrix {
	float* matrix;
	int column;
	int row;
};

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
	Print_Matrix(temp);
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			int flag = 0;
			while(flag < a->column) {
				temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag] * b->matrix[j + flag * b->column];
				flag++;
			}
		}
	}
	return temp;
}



int main() {
	int m, n;
	cin >> m >> n;
	Matrix* a = Build_Matrix(m, n);
//	Print_Matrix(a);
	cin >> m >> n;
	Matrix* b = Build_Matrix(m, n);
//	Print_Matrix(b);
	Print_Matrix(Matrix_Multiplication(a, b));
	return 0;
}
