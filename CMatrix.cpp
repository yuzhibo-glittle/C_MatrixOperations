// CMatrix.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>

//报错提示
void edit_error(char* s1, char* s2);

//矩阵初始化
float** matrix_float(int nrh, int nch);

//释放内存
void free_matrix_float(float** m, int nrh);

//转置
void matrix_t(float** a_matrix, float** b_matrix, int krow, int kline);

//加减
void matrix_a(float** a_matrix, float** b_matrix, float** c_matrix,
	int krow, int kline, int ktrl);

//乘法
void matrix_m(float** a_matrix, float** b_matrix, float** c_matrix,
	int krow, int kline, int kmiddle, int ktrl);

//转置
void  matrix_inv(float** a_matrix, int ndimen);

int main()
{
	float** M_a;
	M_a = matrix_float(3, 3);
	M_a[0][0] = 1;
	M_a[0][1] = 0;
	M_a[0][2] = 1;
	M_a[1][0] = 1;
	M_a[1][1] = 1;
	M_a[1][2] = 1;
	M_a[2][0] = 0;
	M_a[2][1] = 0;
	M_a[2][2] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f ", M_a[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	float** M_b;
	M_b = matrix_float(3, 3);
	matrix_t(M_b, M_a, 3, 3);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f ", M_b[i][j]);
		}
		printf("\n");
	}

	free_matrix_float(M_a, 3);
	free_matrix_float(M_b, 3);
}

void edit_error(char* s1, char* s2)
{
	printf("\n A processing error occured ! \n %s%s\n", s1, s2);
	exit(1);
}

float** matrix_float(int nrh, int nch)
{
	int i, j;
	float** m;

	m = (float**)malloc((unsigned)(nrh) * sizeof(float*));
	if (!m)
		edit_error("allocation failure 1 in matrix()", "");

	for (i = 0; i < nrh; i++) {
		m[i] = (float*)malloc((unsigned)(nch) * sizeof(float));
		if (!m[i])
			edit_error("allocation failure 2 in matrix()", "");
	}
	for (i = 0; i < nrh; i++)
		for (j = 0; j < nch; j++)
			m[i][j] = 0.;
	return m;
}

void free_matrix_float(float** m, int nrh)
{
	int i;
	for (i = nrh - 1; i >= 0; i--)
		free((float*)(m[i]));
}

void matrix_t(float** a_matrix, float** b_matrix, int krow, int kline)

//	a_matrix:转置后的矩阵
//	b_matrix:转置前的矩阵
//	krow    :行数
//	kline   :列数

{
	int k, k2;

	for (k = 0; k < krow; k++)
	{
		for (k2 = 0; k2 < kline; k2++)
		{
			a_matrix[k2][k] = b_matrix[k][k2];
		}
	}
}

void matrix_a(float** a_matrix, float** b_matrix, float** c_matrix,
	int krow, int kline, int ktrl)

	//	a_matrix=b_matrix+c_matrix
	//	 krow   :行数
	//	 kline  :列数
	//	 ktrl   :大于0: 加法  不大于0:减法

{
	int k, k2;

	for (k = 0; k < krow; k++)
	{
		for (k2 = 0; k2 < kline; k2++)
		{
			a_matrix[k][k2] = b_matrix[k][k2]
				+ ((ktrl > 0) ? c_matrix[k][k2] : -c_matrix[k][k2]);
		}
	}
}

void matrix_m(float** a_matrix, float** b_matrix, float** c_matrix,
	int krow, int kline, int kmiddle, int ktrl)

	//	a_matrix=b_matrix*c_matrix
	//	krow  :b的行数
	//	kline :c的列数
	//  kmiddle:b的列数/c的行数
	//	ktrl  :	大于0:两个正数矩阵相乘 不大于0:正数矩阵乘以负数矩阵

{
	int k, k2, k4;
	double stmp;

	for (k = 0; k < krow; k++)
	{
		for (k2 = 0; k2 < kline; k2++)
		{
			stmp = 0.0;
			for (k4 = 0; k4 < kmiddle; k4++)
			{
				stmp += b_matrix[k][k4] * c_matrix[k4][k2];
			}
			a_matrix[k][k2] = stmp;
		}
	}
	if (ktrl <= 0)
	{
		for (k = 0; k < krow; k++)
		{
			for (k2 = 0; k2 < kline; k2++)
			{
				a_matrix[k][k2] = -a_matrix[k][k2];
			}
		}
	}
}

void  matrix_inv(float** a_matrix, int ndimen)

//	a_matrix:矩阵
//	ndimen :维数

{
	double tmp, tmp2, b_tmp[20], c_tmp[20];
	int k, k1, k2, k3, j, i, j2, i2, kme[20], kmf[20];
	i2 = j2 = 0;

	for (k = 0; k < ndimen; k++)
	{
		tmp2 = 0.0;
		for (i = k; i < ndimen; i++)
		{
			for (j = k; j < ndimen; j++)
			{
				if (fabs(a_matrix[i][j]) <= fabs(tmp2))
					continue;
				tmp2 = a_matrix[i][j];
				i2 = i;
				j2 = j;
			}
		}
		if (i2 != k)
		{
			for (j = 0; j < ndimen; j++)
			{
				tmp = a_matrix[i2][j];
				a_matrix[i2][j] = a_matrix[k][j];
				a_matrix[k][j] = tmp;
			}
		}
		if (j2 != k)
		{
			for (i = 0; i < ndimen; i++)
			{
				tmp = a_matrix[i][j2];
				a_matrix[i][j2] = a_matrix[i][k];
				a_matrix[i][k] = tmp;
			}
		}
		kme[k] = i2;
		kmf[k] = j2;
		for (j = 0; j < ndimen; j++)
		{
			if (j == k)
			{
				b_tmp[j] = 1.0 / tmp2;
				c_tmp[j] = 1.0;
			}
			else
			{
				b_tmp[j] = -a_matrix[k][j] / tmp2;
				c_tmp[j] = a_matrix[j][k];
			}
			a_matrix[k][j] = 0.0;
			a_matrix[j][k] = 0.0;
		}
		for (i = 0; i < ndimen; i++)
		{
			for (j = 0; j < ndimen; j++)
			{
				a_matrix[i][j] = a_matrix[i][j] + c_tmp[i] * b_tmp[j];
			}
		}
	}
	for (k3 = 0; k3 < ndimen; k3++)
	{
		k = ndimen - k3 - 1;
		k1 = kme[k];
		k2 = kmf[k];
		if (k1 != k)
		{
			for (i = 0; i < ndimen; i++)
			{
				tmp = a_matrix[i][k1];
				a_matrix[i][k1] = a_matrix[i][k];
				a_matrix[i][k] = tmp;
			}
		}
		if (k2 != k)
		{
			for (j = 0; j < ndimen; j++)
			{
				tmp = a_matrix[k2][j];
				a_matrix[k2][j] = a_matrix[k][j];
				a_matrix[k][j] = tmp;
			}
		}
	}
}
