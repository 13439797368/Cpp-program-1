//????????,????,??cin cout,??? 

#include<iostream>
#include<algorithm>
#include<string.h>
using namespace std;

struct Matrix {
	float* matrix;
	int column;
	int row;
};

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

Matrix* Matrix_Multiplication(Matrix* a, Matrix* b) {
	if (a->column != b->row) {
		cout << "Their column and row do not match." << endl;
		exit(0);
	}
	Matrix* temp = new Matrix;
	temp->column = b->column;
	temp->row = a->row;
	temp->matrix = new float[temp->column * temp->row];
	//cout << sizeof(temp->matrix) << endl;
	//cout << 4 * sizeof(float) << endl;
	memset(temp->matrix, 0, temp->column * temp->row * sizeof(temp->matrix));
//	Print_Matrix(temp);
	for (int i = 0; i < temp->row; i++) {
		for (int j = 0; j < temp->column; j++) {
			int flag = 0;
			while (flag < a->column) {
				temp->matrix[i * temp->column + j] += a->matrix[i * a->column + flag] * b->matrix[j + flag * b->column];
				flag++;
			}
		}
	}
	return temp;
}



int main() {
	string s;
	cin >> s;
	int m = check_int(s);
	cin >> s;
	int n = check_int(s);
	Matrix* a = Build_Matrix(m, n);
	cin >> s;
	m = check_int(s);
	cin >> s;
	n = check_int(s);
	Matrix* b = Build_Matrix(m, n);
	Print_Matrix(Matrix_Multiplication(a, b));
	return 0;
}

