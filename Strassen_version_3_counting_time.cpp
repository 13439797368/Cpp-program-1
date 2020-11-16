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
	int column;
	int row;
};

Matrix* Build_Ones(Matrix* A) {
	//128  9ms
	//256  60ms
	//512  451ms
	//1024 3100ms
	//2048 21263
	A->column = 128;
	A->row = 128;
	A->matrix = new float[A->column * A->row];
	for (int i = 0; i < A->column * A->row; i++) {
		//A->matrix[i] = (float)rand()/10000;
		A->matrix[i] = 1;
	}
	printf_s("build complete\n");
	return A;
}


Matrix* Matrix_Multiplication(Matrix* a, Matrix* b,Matrix* temp) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		return NULL;
	}
	temp->column = b->column;
	temp->row = a->row;
	//temp->matrix = new float[temp->column * temp->row];
	//cout << sizeof(temp->matrix) << endl;
	//cout << 4 * sizeof(float) << endl;
	float* c = (float*)malloc(4 * sizeof(float));
	//memset(temp->matrix, 0, temp->column*temp->row*sizeof(temp->matrix));
	//Print_Matrix(temp);
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
#pragma omp parallel for
			for (int flag = 0; flag < a->column; flag += 4) {
				cv::v_float32x4 va = cv::v_load(a->matrix + i * a->column + flag);
				c[0] = b->matrix[j + flag * b->column];
				c[1] = b->matrix[j + (flag + 1) * b->column];
				c[2] = b->matrix[j + (flag + 2) * b->column];
				c[3] = b->matrix[j + (flag + 3) * b->column];
				cv::v_float32x4 vb = cv::v_load(c);
				cv::v_float32x4 vc = va * vb;
				cv::v_store(c, vc);
				temp->matrix[i * temp->column + j] += c[0]+c[1]+c[2]+c[3];
				//temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag] * b->matrix[j + flag * b->column];
			}
		}
	}
	delete[] c;
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
		//	Print_Matrix(a);
		//	Print_Matrix(b);
		Matrix* c = new Matrix;
		c->matrix = new float[a->row * b->column];
		memset(c->matrix, 0, a->row * b->column * sizeof(float));
		return Matrix_Multiplication(a, b,c);
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

	Build_Ones(a);
	Matrix* c = new Matrix;
	c->matrix = new float[a->row * a->column];
	memset(c->matrix, 0, a->row * a->column * sizeof(float));
	using namespace literals;
	auto now = chrono::system_clock::now();
	auto t_c = chrono::system_clock::to_time_t(now - 24h);
	auto start = std::chrono::steady_clock::now();

	//Matrix*ans = Matrix_Multiplication(a, a,c);
	Matrix* ans = Matrix_Strassen(a, a);
	auto end = std::chrono::steady_clock::now();
	Print_Matrix(ans);
	printf("ms = %d", (int)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
	return 0;
}
