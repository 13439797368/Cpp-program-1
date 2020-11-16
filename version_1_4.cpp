//this version change the recursive condition from a.row=a.column to 1 to 128
//use opencv to speed up
//opencv here is for debug, if i want to run release, reset the "opencv_world3412d.lib" to "opencv_world3412.lib"
//my computer is 128_bit
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

Matrix* A;

void Print_Matrix(Matrix* a) {
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < a->column; j++) {
			printf_s("%lf ", a->matrix[i * a->column + j]);
		}
		cout << endl;
	}
	return;
}

Matrix* Build_Ones(Matrix* A, float t) {
	A->column = 10;
	A->row = 10;
	A->matrix = new float[A->column * A->row];
	A->matrix_c = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix[i] = t;
		A->matrix_c[i / A->row + i % A->row * A->column] = A->matrix[i];
	}
	printf_s("build complete\n");
	return A;
}


Matrix* Matrix_Multiplication(Matrix* a, Matrix* b,Matrix*temp) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		return NULL;
	}
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			cv::v_float32 v_sum	= cv::vx_setzero_f32();
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
			for (int flag = 0; flag < a->column % 4; flag ++) {
				temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag+p] * b->matrix[j + (flag+p) * b->column];
			}
		}
	}
	return temp;
}

int main() {
	Matrix* a = new Matrix;
	Build_Ones(a,1);
	Matrix* ans = new Matrix;
	
	//memset(ans->matrix, 0, ans->column * ans->row * sizeof(ans->matrix));
	using namespace literals;
	auto now = chrono::system_clock::now();
	auto t_c = chrono::system_clock::to_time_t(now - 24h);
	auto start = std::chrono::steady_clock::now();
	Build_Ones(ans, 0);
	Matrix* Ans = Matrix_Multiplication(a, a,ans);

	auto end = std::chrono::steady_clock::now();
	Print_Matrix(ans);
	printf("ms = %d", (int)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
	return 0;
}
