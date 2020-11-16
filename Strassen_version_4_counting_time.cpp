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

Matrix* Build_Ones(Matrix* A, float t,int size) {
	A->column = size;
	A->row = size;
	A->matrix = new float[A->column * A->row];
	A->matrix_c = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix[i] = t;
		A->matrix_c[i / A->row + i % A->row * A->column] = A->matrix[i];
	}
	//printf_s("build complete\n");
	return A;
}

inline void Build_Matrix_c(Matrix* A) {
	A->matrix_c = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		A->matrix_c[i / A->row + i % A->row * A->column] = A->matrix[i];
	}
}

Matrix* Matrix_Multiplication(Matrix* a, Matrix* b, Matrix* temp) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		return NULL;
	}
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
//this function can pick out part of a matrix 
Matrix* Build_Matrix_s(Matrix* a, int row_begin, int row_end, int column_begin, int column_end) {
	Matrix* temp = new Matrix;
	temp->column = column_end - column_begin + 1;
	temp->row = row_end - row_begin + 1;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[(i + row_begin - 1) * a->column + j + column_begin - 1];
		}
	}
	return temp;
}
//this is no use in this code
float Dot_product(int len, float* vector_1, float* vector_2) {
	float sum = 0;
	for (int i = 0; i < len; i++) {
		sum += vector_1[i] * vector_2[i];
	}
	return sum;
}

//this function can make C11,C12,C21,C22 into C
Matrix* Rebuild_BlockMatrix(Matrix* c11, Matrix* c12, Matrix* c21, Matrix* c22) {
	Matrix* temp = new Matrix;
	temp->column = c11->column + c12->column;
	temp->row = c12->row + c21->row;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < c11->row; i++) {
		for (int j = 0; j < c11->column; j++) {
			temp->matrix[i * temp->column + j] = c11->matrix[i * c11->column + j];
		}
		for (int j = 0; j < c12->column; j++) {
			temp->matrix[i * temp->column + j + c11->column] = c12->matrix[i * c12->column + j];
		}
	}
	for (int i = 0; i < c21->row; i++) {
		for (int j = 0; j < c21->column; j++) {
			temp->matrix[(i + c11->row) * temp->column + j] = c21->matrix[i * c21->column + j];
		}
		for (int j = 0; j < c22->column; j++) {
			temp->matrix[(i + c11->row) * temp->column + j + c21->column] = c22->matrix[i * c22->column + j];
		}
	}
	return temp;
}

void Print_Matrix(Matrix* a) {
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < a->column; j++) {
			printf_s("%lf ", a->matrix[i * a->column + j]);
		}
		cout << endl;
	}
	return;
}


//----------------------------those plus or minus function below is designed especially for strassen-----------------------

//C11 = M5 + M4 - M2 + M6
Matrix* Matrix_Plus_s_for_C11(Matrix* a, Matrix* b, Matrix* c, Matrix* d) {
	if (a->row != b->row || a->column != b->column) {
		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[i * a->column + j] + b->matrix[i * b->column + j];

		}
	}
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] += d->matrix[i * a->column + j];

		}
	}
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] -= c->matrix[i * a->column + j];

		}
	}
	return temp;
}

//C12 = M1 + M2
//this is normal plus
Matrix* Matrix_Plus_s_forC12C21(Matrix* a, Matrix* b) {

	if (a->row != b->row || a->column != b->column) {
		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[i * a->column + j] + b->matrix[i * b->column + j];

		}
	}
	return temp;
}

//C22 = M5 + M1 - M3 - M7
Matrix* Matrix_Plus_s_for_C22(Matrix* a, Matrix* b, Matrix* c, Matrix* d) {
	if (a->row != b->row || a->column != b->column) {
		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		//#pragma omp parallel for
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[i * a->column + j] + b->matrix[i * b->column + j];

		}
	}
	for (int i = 0; i < temp->row; i++) {
		//#pragma omp parallel for
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] -= d->matrix[i * a->column + j];

		}
	}
	for (int i = 0; i < temp->row; i++) {
		//#pragma omp parallel for
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] -= c->matrix[i * a->column + j];

		}
	}
	return temp;
}


//this is normal minus, just for pattern
Matrix* Matrix_Minus(Matrix* a, Matrix* b) {

	if (a->row != b->row || a->column != b->column) {
		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		//#pragma omp parallel for
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[i * a->column + j] - b->matrix[i * b->column + j];

		}
	}
	return temp;
}



//realize matrix plus whitch can add part of two matrix to save the space and time creating new submatrix
//folding my code to avoid too long to read
Matrix* Matrix_Plus_s(Matrix* a, int a_row_begin, int a_row_end, int a_column_begin, \
	int a_column_end, Matrix* b, int b_row_begin, int b_row_end, int b_column_begin, int b_column_end) {

	if ((a_column_end - a_column_begin + 1) != (b_column_end - b_column_begin + 1) || \
		(a_row_end - a_row_begin + 1) != (b_row_end - b_row_begin + 1)) {

		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a_column_end - a_column_begin + 1;
	temp->row = a_row_end - a_row_begin + 1;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		//#pragma omp parallel for
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[(i + a_row_begin - 1) * a->column + j + a_column_begin - 1] \
				+ b->matrix[(i + b_row_begin - 1) * b->column + j + b_column_begin - 1];

		}
	}
	return temp;
}

Matrix* Matrix_Minus_s(Matrix* a, int a_row_begin, int a_row_end, int a_column_begin, int a_column_end, \
	Matrix* b, int b_row_begin, int b_row_end, int b_column_begin, int b_column_end) {

	if ((a_column_end - a_column_begin + 1) != (b_column_end - b_column_begin + 1) || \
		(a_row_end - a_row_begin + 1) != (b_row_end - b_row_begin + 1)) {

		cout << "These two matrix can not add together" << endl;
		return NULL;
	}
	Matrix* temp = new Matrix;
	temp->column = a_column_end - a_column_begin + 1;
	temp->row = a_row_end - a_row_begin + 1;
	temp->matrix = new float[temp->column * temp->row];
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			temp->matrix[i * temp->column + j] = a->matrix[(i + a_row_begin - 1) * a->column + j + a_column_begin - 1] \
				- b->matrix[(i + b_row_begin - 1) * b->column + j + b_column_begin - 1];

		}
	}
	return temp;
}

Matrix* c = new Matrix;

Matrix* Matrix_Strassen(Matrix* a, Matrix* b) {
	//divide matrix apart and build new matrix cost so much.

	/*-----------------------------
	C11 = M5 + M4 - M2 + M6

	C12 = M1 + M2

	C21 = M3 + M4

	C22 = M5 + M1 - M3 - M7
	------------------------------*/

	//Matrix* A_1_1 = Build_Matrix_s(a, 1, a->row / 2, 1, a->column / 2);
	//Matrix* A_1_2 = Build_Matrix_s(a, 1, a->row / 2, a->column / 2 + 1, a->column);
	//Matrix* A_2_1 = Build_Matrix_s(a, a->row / 2 + 1, a->row, 1, a->column / 2);
	//Matrix* A_2_2 = Build_Matrix_s(a, a->row / 2 + 1, a->row, a->column / 2 + 1, a->column);
	//Matrix* B_1_1 = Build_Matrix_s(b, 1, b->row / 2, 1, b->column / 2);
	//Matrix* B_1_2 = Build_Matrix_s(b, 1, b->row / 2, b->column / 2 + 1, b->column);
	//Matrix* B_2_1 = Build_Matrix_s(b, b->row / 2 + 1, b->row, 1, b->column / 2);
	//Matrix* B_2_2 = Build_Matrix_s(b, b->row / 2 + 1, b->row, b->column / 2 + 1, b->column);

	/*
	M1=A11(B12-B22)

	M2=(A11+A12)B22

	M3=(A21+A22)B11

	M4=A22(B21-B11)

	M5=(A11+A22)(B11+B22)

	M6=(A12-A22)(B21+B22)

	M7=(A11-A21)(B11+B12)
	*/

	if (a->row == 64 && b->column == 64) {
		Build_Ones(c, 0, 64);
		Build_Matrix_c(b);
		return Matrix_Multiplication(a, b, c);
	}
	Matrix* M1 = Matrix_Strassen(Build_Matrix_s(a, 1, a->row / 2, 1, a->column / 2), Matrix_Minus_s(b, 1, b->row / 2, b->column / 2 + 1, b->column, b, b->row / 2 + 1, b->row, b->column / 2 + 1, b->column));
	Matrix* M2 = Matrix_Strassen(Matrix_Plus_s(a, 1, a->row / 2, 1, a->column / 2, a, 1, a->row / 2, a->column / 2 + 1, a->column), Build_Matrix_s(b, b->row / 2 + 1, b->row, b->column / 2 + 1, b->column));
	Matrix* M3 = Matrix_Strassen(Matrix_Plus_s(a, a->row / 2 + 1, a->row, 1, a->column / 2, a, a->row / 2 + 1, a->row, a->column / 2 + 1, a->column), Build_Matrix_s(b, 1, b->row / 2, 1, b->column / 2));
	Matrix* M4 = Matrix_Strassen(Build_Matrix_s(a, a->row / 2 + 1, a->row, a->column / 2 + 1, a->column), Matrix_Minus_s(b, b->row / 2 + 1, b->row, 1, b->column / 2, b, 1, b->row / 2, 1, b->column / 2));
	Matrix* M5 = Matrix_Strassen(Matrix_Plus_s(a, 1, a->row / 2, 1, a->column / 2, a, a->row / 2 + 1, a->row, a->column / 2 + 1, a->column), Matrix_Plus_s(b, 1, b->row / 2, 1, b->column / 2, b, b->row / 2 + 1, b->row, b->column / 2 + 1, b->column));
	Matrix* M6 = Matrix_Strassen(Matrix_Minus_s(a, 1, a->row / 2, a->column / 2 + 1, a->column, a, a->row / 2 + 1, a->row, 1, a->column / 2), Matrix_Plus_s(b, b->row / 2 + 1, b->row, 1, b->column / 2, b, b->row / 2 + 1, b->row, b->column / 2 + 1, b->column));
	Matrix* M7 = Matrix_Strassen(Matrix_Minus_s(a, 1, a->row / 2, 1, a->column / 2, a, a->row / 2 + 1, a->row, 1, a->column / 2), Matrix_Plus_s(b, 1, b->row / 2, 1, b->column / 2, b, 1, b->row / 2, b->column / 2 + 1, b->column));

	return Rebuild_BlockMatrix(Matrix_Plus_s_for_C11(M5, M4, M2, M6), Matrix_Plus_s_forC12C21(M1, M2), Matrix_Plus_s_forC12C21(M3, M4), Matrix_Plus_s_for_C22(M5, M1, M3, M7));
}


int main() {
	Matrix* a = new Matrix;

	Build_Ones(a,1,1024);
	Matrix* c = new Matrix;

	using namespace literals;
	auto now = chrono::system_clock::now();
	auto t_c = chrono::system_clock::to_time_t(now - 24h);
	auto start = std::chrono::steady_clock::now();
	//Build_Ones(c, 0, 64);
	//Matrix* ans = Matrix_Multiplication(a, a, c);
	Matrix* ans = Matrix_Strassen(a, a);
	auto end = std::chrono::steady_clock::now();
	//Print_Matrix(ans);
	printf("ms = %d", (int)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
	return 0;
}
