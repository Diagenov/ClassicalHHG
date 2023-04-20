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
const double Up = 57.38;  // �������������� �������, ��
const double Ip = 21.55;  // ��������� ��������� (Ne) ��
const double Wxuv = 30;   // ������� XUV, ��
const double Wir = 1;     // ������� IR, ��
const double Txuv = 0.55; // ����� XUV, ��
double Tir = 20;    // ����� IR, ��
#pragma endregion

#pragma region ����� ������ (����������, ��������, IR, IR + XUV)
int monochromate = 1; // 1 - ���� ����������, 2 - ���� ��������
int xuv_ir = 1;       // 1 - ���� IR, 2 - ���� IR + XUV
#pragma endregion

#pragma region ��������� ���������
#define N 6     // ���������� ����������
#define M 100   // ���������� ����� �� ���� �����
#define D 10000 // ����������� �� ������ � ����  
#pragma endregion

#pragma region ��� ������ � ���� ����
double array_t[D];
double array_E[D];
int nextIndex = 0;
#pragma endregion

#pragma region ��������� ��� ������� ��������� ��� ��� ������� �������� S
struct parameters
{
	double E; // ������� HHG (��� ���������, ������� �� ������)
	double t1; // ����� ���������
	double t2; // ����� ������������ 
};
#pragma endregion

#pragma region �������, �������� IR-���� � ������� ���������
double C = 0;

double F(double t)
{
	if (monochromate == 1)
		return cos(t);

	if (t < 0 || t > Tir)
		return 0;

	return cos(t) * pow(sin(M_PI * t / Tir), 2);
}

double A(double t)
{
	if (monochromate == 1)
		return -sin(t);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;

	double y1 = -sin(t) / 2;
	double y2 = sin(a1 * t) / (4 * a1);
	double y3 = sin(a2 * t) / (4 * a2);
	return y1 + y2 + y3 + C;
}

double IntA(double t)
{
	if (monochromate == 1)
		return cos(t);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;

	double y1 = cos(t) / 2;
	double y2 = -cos(a1 * t) / (4 * a1 * a1);
	double y3 = -cos(a2 * t) / (4 * a2 * a2);
	return y1 + y2 + y3 + (C * t);
}

double Q(double t1, double t2)
{
	if (t1 == t2)
		return 0;

	return -(IntA(t2) - IntA(t1)) / (t2 - t1);
}

double P(double t, double t1, double t2)
{
	return A(t) + Q(t1, t2);
}

double deltaE(double t1, double t2)
{
	double a = Ip / (t2 - t1);
	return a * P(t2, t1, t2) / F(t1);
}
#pragma endregion

#pragma region ������������ ���������
double maxHHG(double t1, double t2)
{
	double deltaHHG = xuv_ir == 1 ? 0 : F(t2) * (Wxuv - Ip) / F(t1);
	return (2 * Up * pow(A(t2) - A(t1), 2)) + deltaHHG;
}
#pragma endregion

#pragma region ������� ���������
int maxHHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	const double e1 = P(t1, t1, t2);
	const double e2 = F(t2) + ((A(t2) - A(t1)) / (t2 - t1));

	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}

int HHG_f(const gsl_vector * x, void * params, gsl_vector * func)
{
	struct parameters * p = (struct parameters *)params;
	const double E = (p->E);

	const double t1 = gsl_vector_get(x, 0);
	const double t2 = gsl_vector_get(x, 1);

	double e1, e2;
	if (xuv_ir == 1)
	{
		e1 = P(t1, t1, t2);
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - E;
	}
	else
	{
		e1 = (2 * Up * pow(P(t1, t1, t2), 2)) - (Wxuv - Ip);
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - E;
	}
	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}
#pragma endregion

#pragma region ����� ����������
int max_t(double to1, double max_t1[N][2], double max_t2[N][2], double max_e[N])
{
	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;

	int count = 0;
	const size_t n = 2;

	gsl_multiroot_function func =
	{
		&maxHHG_f,
		n,
		0
	};

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, n);
	gsl_vector * x = gsl_vector_alloc(n);

	for (int i = 0; i < N; i++)
	{
		double to2 = to1 + (i * M_PI) + M_PI_2;
		if (to2 > Tir)
			break;

		gsl_vector_set(x, 0, to1);
		gsl_vector_set(x, 1, to2);
		gsl_multiroot_fsolver_set(s, &func, x);

		int status;
		int iter = 0;
		do
		{
			status = gsl_multiroot_fsolver_iterate(s);
			if (status)
				break;

			status = gsl_multiroot_test_residual(s->f, 1e-7);
		} 
		while (status == GSL_CONTINUE && ++iter < 1000);

		if (status != GSL_SUCCESS)
			continue;
		
		double t1 = gsl_vector_get(s->x, 0);
		double t2 = gsl_vector_get(s->x, 1);
		if (t1 < 0 || t2 > Tir || fabs(t2 - t1) < M_PI_4)
			continue;

		double hhg = maxHHG(t1, t2);
		if (hhg < 5)
			continue;

		max_t1[count][0] = max_t1[count][1] = t1;
		max_t2[count][0] = max_t2[count][1] = t2;
		max_e[count] = hhg;

		printf("t1 = %.3e \n", t1 / 1.6);
		printf("t2 = %.3e \n", t2 / 1.6);
		printf("t = %.3e \n", (t2 - t1) / 1.6);
		if (xuv_ir == 1)
		{
			printf("maxE = %.3e Up = %.3e \n\n", 2 * pow(A(t2) - A(t1), 2), hhg);
		}
		else
		{
			printf("maxE = %.3e Up  +  %.3e (Wxuv - Ip) = %.3e \n\n", 2 * pow(A(t2) - A(t1), 2), (F(t2) / F(t1)), hhg);
		}
		count++;
	}
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return count;
}
#pragma endregion

#pragma region ���������� ����������
int trajectories(double to1)
{
    #pragma region ���������� ����������
	double max_t1[N][2];
	double max_t2[N][2];
	double max_e[N];

	const int n = max_t(to1, max_t1, max_t2, max_e);
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
		array_t[nextIndex] = max_t2[i][0] / 1.6;
		array_E[nextIndex] = max_e[i];
		nextIndex++;

		for (int j = 0; j < 2; j++)
		{
			double E = max_e[i];
			int sign = (2 * j) - 1;
			max_t2[i][j] += (sign * 0.25);

			while (E > 0)
			{
				E -= 1.0;
				struct parameters p = { E, 0, 0 };
				func.params = &p;

				gsl_vector_set(x, 0, max_t1[i][j]);
				gsl_vector_set(x, 1, max_t2[i][j]);
				gsl_multiroot_fsolver_set(s, &func, x);

				int status;
				int iter = 0;
				do
				{
					status = gsl_multiroot_fsolver_iterate(s);
					if (status)
						break;

					status = gsl_multiroot_test_residual(s->f, 1e-7);
				} 
				while (status == GSL_CONTINUE && ++iter < 1000);

				if (status != GSL_SUCCESS)
					continue;
				
				double t1 = gsl_vector_get(s->x, 0);
				double t2 = gsl_vector_get(s->x, 1);

				max_t1[i][j] = t1;
				max_t2[i][j] = t2;

				array_t[nextIndex] = t2 / 1.6;
				array_E[nextIndex] = E;
				nextIndex++;
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
    #pragma region ����������
	FILE* fp;
	fp = fopen("HHG.txt", "w");
	for (int i = 0; i < nextIndex; i++)
	{
		fprintf(fp, "%e  %e\n", array_t[i], array_E[i] + Ip);
	}
	fclose(fp);
    #pragma endregion

	int K = 1000;
	double dt = Tir / K;

    #pragma region IR-�������������
	fp = fopen("F.txt", "w");
	for (int i = 0; i < K; i++)
	{
		double t = i * dt;
		fprintf(fp, "%e  %e\n", t, 35 * F(t));
	}
	fclose(fp);
    #pragma endregion

    #pragma region IR-���������
	fp = fopen("A.txt", "w");
	for (int i = 0; i < K; i++)
	{
		double t = i * dt;
		fprintf(fp, "%e  %e\n", t, 35 * A(t));
	}
	fclose(fp);
    #pragma endregion
}
#pragma endregion

int work()
{
    #pragma region �������� ��� ����������
	printf("1 - monochromate, 2 - impulses\n");
	do
	{
		printf("  ");
		scanf("%d", &monochromate);
	} 
	while (monochromate != 1 && monochromate != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region ����� IR
	if (monochromate == 2)
	{
		printf("Impulse length (fs)\n");
		do
		{
			printf("  ");
			scanf("%lf", &Tir);
		} 
		while (Tir < 3 * M_PI * 1.6 || Tir > 10 * M_PI * 1.6);
		printf("\n\n");
		Tir *= 1.6; // ��� �������� �� ����������� (fs) � ������������ ������� (W * t, �� ���� 1/s * s = 1, s - �������)
	}
    #pragma endregion

    #pragma region IR ��� IR + XUV
	printf("1 - IR, 2 - IR + XUV\n");
	do
	{
		printf("  ");
		scanf("%d", &xuv_ir);
	} 
	while (xuv_ir != 1 && xuv_ir != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region ������ ���������
	double t1 = -1;
	if (xuv_ir == 2)
	{
		printf("Ionization time (fs)\n");
		do
		{
			printf("  ");
			scanf("%lf", &t1);
		} 
		while (t1 < 0 || t1 * 1.6 > Tir);
		printf("\n\n");
		t1 *= 1.6;
	}
    #pragma endregion

    #pragma region ����������� ��������� C
	C = -(IntA(Tir) - IntA(0)); // 20: -0.1043234466544098, 7pi: 0
    #pragma endregion

    #pragma region ���������� ����������
	if (xuv_ir == 2) 
	{
		printf("\n\nstart t1 = %.3e fs\n", t1 / 1.6);
		printf("------------------\n");
		if (fabs(F(t1)) < 0.3 || trajectories(t1) == 0)
		{
			printf("nothing\n");
		}
	}

	for (int i = 0; xuv_ir == 1; i++)
	{
		double to1 = i * M_PI;
		if (to1 > Tir)
			break;

		printf("\n\nstart t1 = %.3e fs\n", i * M_PI / 1.6);
		printf("------------------\n");
		if (fabs(F(to1)) < 0.3 || trajectories(to1) == 0)
		{
			printf("nothing\n");
		}
	}

	writeFILE();
	printf("\n\n");
	return 0;
    #pragma endregion
}

int main(void)
{
	do
	{
		C = 0;
		monochromate = 0;
		xuv_ir = 0;
		nextIndex = 0;
		work();
	} 
	while (1 == 1);
	return 0;
}