#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// ���� ���:
//   1) ���� ����� �����
//       + ������ ���� ������ ����� ����� ����� �� ���� �� ����, ��� ������
//       + �������� �� ������������ ������������ �������, � ���������� �� �� ����������
//       - ���������, ��� �������� ��������� �������
// ���� ���:
//   2) �� ������ ���� ������ ���������� ��������� ��������, � �� ������ �� ������

#pragma region ��������� IR, XUV � �����
const double Up = 53.74;  // �������������� �������, ��
const double Ip = 21.55;  // ��������� ��������� (Ne) ��
const double Wxuv = 30;   // ������� XUV, ��
const double Wir = 1;     // ������� IR, ��
const double Txuv = 0.55; // ����� XUV, ��
const double Tir = 20;    // ����� IR, ��
const double fi = 0 * M_PI / 2; // ������������� ���� ���������
#pragma endregion

#pragma region ����� ������ (����������, ��������, IR, IR + XUV)
int monochromate = 1; // 1 - ���� ����������, 2 - ���� ��������
int xuv_ir = 1;       // 1 - ���� IR, 2 - ���� IR + XUV
#pragma endregion

#pragma region ��������� ���������
#define N 5     // ���������� ����������
#define M 100   // ���������� ����� �� ���� �����
#define D 10000 // ����������� �� ������ � ����  
#pragma endregion

#pragma region ��� ������ � ���� ����
double array_t[D];
double array_E[D];
int nextIndex = 0;
#pragma endregion

#pragma region ��������� ������� ���������
struct parameters
{
	double E; // ������� HHG
	double k; // �������� ������� (����������������)
	double C; // ���������
};
#pragma endregion

#pragma region �������, �������� IR-���� � ������� ���������
double F(double t)
{
	if (monochromate == 1)
		return cos(t);

	if (t < 0 || t > Tir)
		return 0;

	return cos(t - fi) * pow(sin(M_PI * t / Tir), 2);
}

double A(double t, double C)
{
	if (monochromate == 1)
		return -sin(t);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	double y = (sin(t - fi) / 2) - (sin(a1 * t - fi) / (4 * a1)) - (sin(a2 * t + fi) / (4 * a2)) + C;
	return -y;
}

double IntA(double t, double C)
{
	if (monochromate == 1)
		return cos(t);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	return (cos(t - fi) / 2) - (cos(a1 * t - fi) / (4 * a1 * a1)) - (cos(a2 * t + fi) / (4 * a2 * a2)) + (C * t);
}

double Q(double t1, double t2, double C)
{
	if (t1 == t2)
		return 0;

	return -(IntA(t2, C) - IntA(t1, C)) / (t2 - t1);
}

double P(double t, double t1, double t2, double C)
{
	return A(t, C) + Q(t1, t2, C);
}

double deltaE(double t1, double t2, double C)
{
	double sqrK = 2 * Ip;
	double a = sqrK / (t2 - t1);
	return a * P(t2, t1, t2, C) / (Wir * F(t1));
}
#pragma endregion

#pragma region ������������ ���������
double maxHHG(double t1, double t2, double C)
{
	double deltaHHG = xuv_ir == 1 ? 0 : F(t2) * (Wxuv - Ip) / F(t1);
	return (2 * Up * pow(A(t2, C) - A(t1, C), 2)) + deltaHHG;
}
#pragma endregion

#pragma region ������� ���������
int maxHHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	struct parameters * p = (struct parameters *)params;
	const double C = (p->C);

	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	const double e1 = P(t1, t1, t2, C);
	const double e2 = F(t2) + ((A(t2, C) - A(t1, C)) / (t2 - t1));

	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}

int HHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	struct parameters * p = (struct parameters *)params;
	const double E = (p->E);
	const double k = (p->k);
	const double C = (p->C);

	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	double e1, e2;
	if (k == 0)
	{
		e1 = P(t1, t1, t2, C);
		e2 = (2 * Up * pow(P(t2, t1, t2, C), 2)) - (E/* - Ip*/)/* + (deltaE(t1, t2, C) / 2)*/;
	}
	else
	{
		e1 = (2 * Up * pow(P(t1, t1, t2, C), 2)) - (Wxuv - Ip);
		e2 = (2 * Up * pow(P(t2, t1, t2, C), 2)) - (E/* - Ip*/);
	}
	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}
#pragma endregion

#pragma region ����� ����������
int max_t(double to1, double max_t1[N][2], double max_t2[N][2], double max_hhg[N], double C)
{
	if (to1 > Tir)
		return 0;

	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;

	int status;
	int count = 0;
	size_t iter = 0;

	const size_t n = 2;
	struct parameters p = { 0, 0, C };

	gsl_multiroot_function func =
	{
		&maxHHG_f,
		n,
		&p
	};

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, n);
	gsl_vector * x = gsl_vector_alloc(n);

	for (int i = 0; i < N; i++)
	{
		double to2 = to1 + ((2 + i) * M_PI) - (M_PI / 2);
		if (to2 > Tir)
			break;

		gsl_vector_set(x, 0, to1);
		gsl_vector_set(x, 1, to2);
		gsl_multiroot_fsolver_set(s, &func, x);

		do
		{
			status = gsl_multiroot_fsolver_iterate(s);
			if (status)
				break;

			status = gsl_multiroot_test_residual(s->f, 1e-7);
		} 
		while (status == GSL_CONTINUE && ++iter < 1000);

		if (status == GSL_SUCCESS)
		{
			double t1 = gsl_vector_get(s->x, 0);
			double t2 = gsl_vector_get(s->x, 1);
			if (t1 < 0 || t2 > Tir)
				continue;

			max_t1[i][0] = max_t1[i][1] = t1;
			max_t2[i][0] = max_t2[i][1] = t2;

			double hhg = maxHHG(t1, t2, C);
			max_hhg[i] = hhg/* + Ip*/;

			printf("t1 = %.3e \n", t1);
			printf("t2 = %.3e \n", t2);
			printf("t = %.3e \n", t2 - t1);
			if (xuv_ir == 1)
			{
				printf("maxHHG = %.3e Up = %.3e \n\n", 2 * pow(A(t2, C) - A(t1, C), 2), hhg);
			}
			else
			{
				printf("maxHHG = %.3e Up  +  %.3e (Wxuv - Ip) = %.3e \n\n", 2 * pow(A(t2, C) - A(t1, C), 2), (F(t2) / F(t1)), hhg);
			}
			count++;
		}
	}
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return count;
}
#pragma endregion

#pragma region ���������� ����������
int trajectories(double k, double C, double to1)
{
    #pragma region ���������� ����������
	double max_t1[N][2];
	double max_t2[N][2];
	double max_hhg[N];

	const int n = max_t(to1, max_t1, max_t2, max_hhg, C);
	if (n == 0)
		return 0;
    #pragma endregion

    #pragma region ����������� ��������
	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;

	const size_t size = 2;
	gsl_multiroot_function func =
	{
		&HHG_f,
		size,
		0
	};

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, size);
	gsl_vector * x = gsl_vector_alloc(size);
    #pragma endregion

    #pragma region ��������� ����������
	for (int i = 0; i < n; i++)
	{
		double maxE = max_hhg[i];
		double dE = maxE / M;

		for (int j = 0; j < 2; j++)
		{
			int sign = (2 * j) - 1;
			max_t2[i][j] += sign * 0.25;

			for (int l = M; l > 0; l--)
			{
				double E = l * dE;

				struct parameters p = { E, k, C };
				func.params = &p;

				gsl_vector_set(x, 0, max_t1[i][j]);
				gsl_vector_set(x, 1, max_t2[i][j]);
				gsl_multiroot_fsolver_set(s, &func, x);

				int status;
				size_t iter = 0;

				do
				{
					status = gsl_multiroot_fsolver_iterate(s);
					if (status)
						break;

					status = gsl_multiroot_test_residual(s->f, 1e-7);
				} 
				while (status == GSL_CONTINUE && ++iter < 1000);

				if (status == GSL_SUCCESS)
				{
					double t1 = gsl_vector_get(s->x, 0);
					double t2 = gsl_vector_get(s->x, 1);

					max_t1[i][j] = t1;
					max_t2[i][j] = t2;

					if (nextIndex == D)
						continue;

					array_t[nextIndex] = t2;
					array_E[nextIndex] = E;
					nextIndex++;
				}
			}
		}
	}
    #pragma endregion

    #pragma region ������������ ������ ��������
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return n;
    #pragma endregion
}
#pragma endregion

#pragma region ������ � ����
void writeFILE()
{
	FILE * fp;
	fp = fopen("HHG.txt", "w");
	for (int i = 0; i < nextIndex; i++)
	{
		fprintf(fp, "%e  %e\n", array_t[i], array_E[i]);
	}
	fclose(fp);
}
#pragma endregion

int main(void)
{
	do
	{
		monochromate = 0;
		xuv_ir = 0;
		nextIndex = 0;
		work();
	} 
	while (1 == 1);
	return 0;
}

int work()
{
    #pragma region �������� ��� ����������
	printf("1 - monochromate, 2 - impulses\n  ");
	do
	{
		scanf("%d", &monochromate);
		printf("  ");
	} 
	while (monochromate != 1 && monochromate != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region IR ��� IR + XUV
	printf("1 - IR, 2 - IR + XUV\n  ");
	do
	{
		scanf("%d", &xuv_ir);
		printf("  ");
	} 
	while (xuv_ir != 1 && xuv_ir != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region ����������� ��������� C
	double C = -1;
	for (int i = 0; i < M; i++)
	{
		double c = (i * 2.0 / M) - 1;
		if (abs(IntA(Tir, C) - IntA(0, C)) > abs(IntA(Tir, c) - IntA(0, c)))
		{
			C = c;
		}
	}
    #pragma endregion

    #pragma region �������� ������� (����������������, XUV + IR)
	const double k = xuv_ir == 1 ? 0 : sqrt((Wxuv - Ip) / (2 * Up)); 
    #pragma endregion

    #pragma region ���������� ����������
	for (int i = 0; ; i++)
	{
		double to1 = fi + (i * M_PI);
		if (to1 > Tir)
			break;

		printf("\n\nstart t1 = %d * PI\n", i);
		printf("------------------\n", i);
		if (trajectories(k, C, to1) == 0)
		{
			printf("nothing\n");
		}
	}
	writeFILE();
	printf("\n\n");
	return 0;
    #pragma endregion
}