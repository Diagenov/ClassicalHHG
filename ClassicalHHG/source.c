#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

const double Up = 53.74;  // �������������� �������, ��
const double Ip = 21.55;  // ��������� ��������� (Ne) ��
const double Wxuv = 30;   // ������� XUV, ��
const double Wir = 1;     // ������� IR, ��
const double Txuv = 0.55; // ����� XUV, ��
const double Tir = 20;    // ����� IR, ��
const double delay = 2 * M_PI; // �������� ����� XUV � IR, ��

#define monochromate 1 // 1 - ���� ����������, 0 - ���� ��������
#define N 5            // ���������� ����������
#define M 500          // ���������� ����� ������ (�����) ������� ���������
#define O 5 * 500 * 2  // ����� ������� ����� ������ � ����� ������� ���������

struct parameters
{
	double E; // ������� HHG
	double k; // �������� �������
	double C; // ���������
};

double F(double t)
{
	if (monochromate == 1)
		return cos(t);

	return cos(t) * pow(sin(M_PI * t / Tir), 2);
}

double A(double t, double C)
{
	if (monochromate == 1)
		return -sin(t);

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	double y = (sin(t) / 2) - (sin(a1 * t) / (4 * a1)) - (sin(a2 * t) / (4 * a2)) + C;
	return -y;
}

double IntA(double t, double C)
{
	if (monochromate == 1)
		return cos(t);

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	return (cos(t) / 2) - (cos(a1 * t) / (4 * a1 * a1)) - (cos(a2 * t) / (4 * a2 * a2)) + (C * t);
}

double Q(double t1, double t2, double C)
{
	if (t1 == t2)
		return 0;

	return -(IntA(t2, C) - IntA(t1, C)) / (t2 - t1);
}

double maxHHG(double t1, double t2, double C)
{
	return (2 * Up * pow(A(t2, C) - A(t1, C), 2)) + (F(t2) * (Wxuv - Ip) / F(t1));
}

int maxHHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	struct parameters* p = (struct parameters*)params;
	const double C = (p->C);

	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	const double e1 = A(t1, C) + Q(t1, t2, C);
	const double e2 = F(t2) + ((A(t2, C) - A(t1, C)) / (t2 - t1));

	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}

int HHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	struct parameters* p = (struct parameters*)params;
	const double E = (p->E);
	const double k = (p->k);
	const double C = (p->C);

	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	const double e1 = A(t1, C) + Q(t1, t2, C) - k;
	const double e2 = (2 * Up * pow(Q(t1, t2, C) + A(t2, C), 2)) - E/* + Ip*/;

	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}

int max_t(double max_t1[], double max_t2[], double max_hhg[], double C)
{
	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;

	int status;
	size_t iter = 0;

	const size_t n = 2;
	struct parameters p =
	{
		0,
		0,
		C
	};
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
		gsl_vector_set(x, 0, delay);       
		gsl_vector_set(x, 1, delay + ((2 + i) * M_PI) - (M_PI / 2));
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
			double t1 = max_t1[i] = gsl_vector_get(s->x, 0);
			double t2 = max_t2[i] = gsl_vector_get(s->x, 1);
			double hhg = max_hhg[i] = maxHHG(t1, t2, C);

			printf("t1 = %.3e \n", t1);
			printf("t2 = %.3e \n", t2);
			printf("t = %.3e \n", t2 - t1);
			printf("maxHHG = %.3e Up  +  %.3e (Wxuv - Ip) = %.3e \n\n", 2 * pow(A(t2, C) - A(t1, C), 2), (F(t2) / F(t1)), hhg);
		}
	}
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return 0;
}

void saveHHG(int n, double array_t2[], double array_E[], char* name)
{
	FILE* fp;
	fp = fopen(name, "w");
	for (int i = 0; i < n; i++)
	{
		fprintf(fp, "%e  %e\n", array_t2[i], array_E[i]);
	}
	fclose(fp);
}

int main(void)
{
	double C = -1;
	for (int i = 0; i < M; i++)
	{
		double c = (i * 2.0 / M) - 1;
		if (abs(IntA(Tir, C) - IntA(0, C)) > abs(IntA(Tir, c) - IntA(0, c)))
		{
			C = c;
		}
	}

	const double k = sqrt((Wxuv - Ip) / (2 * Up)); // �������� ������� (XUV + IR), ������������

	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;
	
	double max_t1[N];
	double max_t2[N];
	double max_hhg[N];

	max_t(max_t1, max_t2, max_hhg, C);

	int o = 0;
	double dE = max_hhg[0] / M;

	double array_E[O];
	double array_t2[O];

	const size_t n = 2;
	gsl_multiroot_function func = 
	{ 
		&HHG_f, 
		n, 
		0 
	};
	
	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, n);
	gsl_vector * x = gsl_vector_alloc(n);

	for (int i = M; i > 0; i--) 
	{
		double E = i * dE;

		for (int j = 0; j < N; j++)
		{
			if (E > max_hhg[j])
				continue;

			struct parameters p =
			{
				E,
				k,
				C
			};
			func.params = &p;

			for (int k = -1; k < 2; k += 2)
			{
				gsl_vector_set(x, 0, max_t1[j]);
				gsl_vector_set(x, 1, max_t2[j] + (k * 0.002 * (M - i)));
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
					array_t2[o] = gsl_vector_get(s->x, 1);
					array_E[o] = E;
					o++;
				}
			}
		}
	}

	saveHHG(o, array_t2, array_E, "HHG.txt");

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	system("PAUSE");
	return 0;
}













