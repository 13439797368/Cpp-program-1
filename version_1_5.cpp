//????????,????,??cin cout,??? 

#include<iostream>
#include<algorithm>
#include<cstring>
#include<chrono>
#include <stdlib.h>
#include<omp.h>
#include<stdio.h>
#include<opencv2/opencv.hpp>
#include"opencv2/core/simd_intrinsics.hpp"
using namespace std;

struct Matrix {
	float* matrix;
	float* matrix_c;
	int column;
	int row;
};

Matrix* Build_Ones(Matrix* A, float t,int col,int row) {
	A->column = col;
	A->row = row;
	A->matrix = new float[A->column * A->row];
	A->matrix_c = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix[i] = t;
		A->matrix_c[i / A->row + i % A->row * A->column] = A->matrix[i];
	}
	//printf_s("build complete\n");
	return A;
}

int check_int(string s) {
	int ans = 0;
	if (s[0] == '-' || s[0] == '0') {
		cout << "your input is illegal or nonpositive" << endl;
		exit(0);
	}
	for (int i = 0; i < s.length(); i++) {
		if (isdigit(s[i])) {
			ans = ans * 10 + s[i] - '0';
		}
		else {
			cout << "your input is illegal, this must be integer" << endl;
			exit(0);
		}
	}
	return ans;
}

float check_float(string s) {
	float ans = 0;
	int i = 0;
	int neg = 1;
	if (s[i] == '-') {
		neg = -1;
		i++;
	}
	if (s[i] == '0' && s[i + 1] != '.') {
		cout << "0 can not be initial number like you've printed" << endl;
		exit(0);
	}
	for (; i < s.length(); i++) {
		if (isdigit(s[i])) {
			ans = ans * 10 + s[i] - '0';
		}
		else if (s[i] == '.') {
			if (i == s.length() - 1) {
				cout << "numbers are missing after '.'" << endl;
				exit(0);
			}
			else {
				i++;
				break;
			}
		}
		else {
			cout << "this is not a float" << endl;
			exit(0);
		}
	}
	int j = i - 1;
	for (; i < s.length(); i++) {
		if (isdigit(s[i])) {
			ans += (s[i] - '0') / pow(10, j);
			j++;
		}
		else {
			cout << "this is not a float" << endl;
			exit(0);
		}
	}
	return neg * ans;
}

Matrix* Build_Matrix(int row, int col) {
	Matrix* temp = new Matrix;
	temp->column = col;
	temp->row = row;
	temp->matrix = new float[col * row];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			string s;
			cin >> s;
			temp->matrix[i * col + j] = check_float(s);
		}
	}
	return temp;
}

inline void Build_Matrix_c(Matrix* A) {
	A->matrix_c = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix_c[i / A->row + i % A->row * A->column] = A->matrix[i];
	}
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

Matrix* Matrix_Multiplication(Matrix* a, Matrix* b, Matrix* temp) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		return NULL;
	}
	Build_Matrix_c(b);
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			cv::v_float32 v_sum = cv::vx_setzero_f32();
			int p = a->column;
			if (a->column % 4 != 0) {
				p = a->column - a->column % 4;
			}
			for (int flag = 0; flag < p; flag += 4) {
				cv::v_float32x4 va = cv::v_load(a->matrix + i * a->column + flag);
				cv::v_float32x4 vb = cv::v_load(b->matrix_c + flag * b->row + flag);
				v_sum += va * vb;
				//temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag] * b->matrix[j + flag * b->column];
			}
			temp->matrix[i * temp->column + j] = cv::v_reduce_sum(v_sum);
			for (int flag = 0; flag < a->column % 4; flag++) {
				temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag + p] * b->matrix[j + (flag + p) * b->column];
			}
		}
	}
	return temp;
}



int main() {
	cout<<"please enter the column and row number of the first matrix, integers only:"<<endl;
	string s;
	cin >> s;
	int m = check_int(s);
	int M = m;
	cin >> s;
	int n = check_int(s);
	cout<<"please enter the elements of the first matrix, there are "<< m*n<<" in total, floats or integers only:"<<endl;
	Matrix* a = Build_Matrix(m, n);
	cout<<"please enter the column and row number of the first matrix, integers only:"<<endl;
	cin >> s;
	m = check_int(s);
	cin >> s;
	n = check_int(s);
	cout<<"please enter the elements of the first matrix, there are "<< m*n<<" in total, floats or integers only:"<<endl;
	Matrix* b = Build_Matrix(m, n);
	Matrix* c = new Matrix;
	Build_Ones(c, 0,M,n);
	cout<<"the result is:"<<endl;
	Print_Matrix(Matrix_Multiplication(a, b, c));
	return 0;
}

