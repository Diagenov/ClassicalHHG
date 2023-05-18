#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#pragma region параметры полей и атомной мишени 
double Up = 57.38; // пондермоторная энергия, эВ
double Ip = 21.55; // потенциал ионизации (Ne) эВ
double Wir = 1;   // частота IR, эВ
double Tir = 20;  // время IR, фс
double Wxuv = 30; // частота XUV, эВ
#pragma endregion

#pragma region вспомогательные параметры
int xuv_ir;           // 1 - поле IR, 2 - поле IR + XUV (+), 3 - поле IR + XUV (-)
double max_E;         // самый высокий максимум
double max_to1;       // время ионизации, при котором достигается самый высокий максимум
double units = 1.555; // переводная единица (для перевода фс в безразмерные и обратно) (1.5186)
#pragma endregion

#pragma region параметры программы
#define N 10     // количество максимумов
#define M 100   // количество точек на одну ветку
#define D 10000 // ограничение на запись в файл  
#define B 1000  // второе ограничение
#pragma endregion

#pragma region для записи в файл
double array_maxE[B];  
double array_Wxuv[B];
int nextIndex_max = 0;

double array_E[D];
double array_t[D];
int nextIndex = 0;
#pragma endregion

#pragma region параметры для системы уравнений
struct parameters
{
	double E; // энергия HHG (или электрона, зависит от записи)
};
#pragma endregion

#pragma region функции, задающие IR-поле и импульс электрона
double C = 0;

double F(double t)
{
	if (t < 0 || t > Tir)
		return 0;

	return cos(t) * pow(sin(M_PI * t / Tir), 2);
}

double A(double t)
{
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
#pragma endregion

#pragma region максимальная энергия электрона (от которой зависит HHG)
double deltaE(double t1, double t2)
{
	if (xuv_ir != 1)
		return 0;

	double a = Ip / (t2 - t1);
	return a * P(t2, t1, t2) / F(t1);
}

double maxHHG(double t1, double t2)
{
	double delta = xuv_ir == 1 ? 0 : F(t2) * (Wxuv - Ip) / F(t1);
	return (2 * Up * pow(A(t2) - A(t1), 2)) + delta;
}
#pragma endregion

#pragma region системы уравнений
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
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - (E/* + deltaE(t1, t2)*/);
	}
	else
	{
		if (xuv_ir == 2)
		{
			e1 = P(t1, t1, t2) - sqrt((Wxuv - Ip) / (2 * Up));
		}
		else
		{
			e1 = P(t1, t1, t2) + sqrt((Wxuv - Ip) / (2 * Up));
		}
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - E;
	}
	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}
#pragma endregion

#pragma region поиск максимумов
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
		if (t1 < 0 || t2 > Tir || t2 - t1 < M_PI_2)
			continue;

		double maxE = maxHHG(t1, t2);
		if (maxE < 5)
			continue;

		max_t1[count][0] = max_t1[count][1] = t1;
		max_t2[count][0] = max_t2[count][1] = t2;
		max_e[count] = maxE;
		count++;
	}
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return count;
}
#pragma endregion

#pragma region построение траекторий
int trajectories(double to1)
{
	if (fabs(F(to1)) < 0.3)
		return 0;

    #pragma region нахождение максимумов
	double max_t1[N][2];
	double max_t2[N][2];
	double max_e[N];

	const int n = max_t(to1, max_t1, max_t2, max_e);
	if (n == 0)
		return 0;
    #pragma endregion

    #pragma region определение решателя
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

    #pragma region отрисовка траекторий
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			int k = -1;
			double E = 5;

			int sign = (2 * j) - 1;
			max_t2[i][j] += (sign * M_PI / 5);

			while (E < max_e[i] + 10)
			{
				k++;
				E = 5 + (k * 0.1);
				struct parameters p = { E };
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
				{
					k++;
					E = 5 + (k * 0.1);
					continue;
				}
				
				if (E > max_E)
				{
					max_E = E;
					max_to1 = to1;
				}

				double t1 = gsl_vector_get(s->x, 0);
				double t2 = gsl_vector_get(s->x, 1);

				max_t1[i][j] = t1;
				max_t2[i][j] = t2;

				if (k % 10 == 0) 
				{
					array_t[nextIndex] = t2 / units;
					array_E[nextIndex] = E + deltaE(t1, t2);
					nextIndex++;
				}

				if (t2 - t1 <= 0.05)
				{
					break;
				}
			}
		}
	}
    #pragma endregion

    #pragma region освобождение памяти решателя
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return n;
    #pragma endregion
}
#pragma endregion

#pragma region запись в файл
void writeFILE_HHG(char* fileName)
{
	FILE* fp;
	fp = fopen(fileName, "w");
	for (int i = 0; i < nextIndex; i++)
	{
		fprintf(fp, "%e  %e\n", array_t[i], array_E[i] + Ip);
	}
	nextIndex = 0;
	fclose(fp);
}

void writeFILE_IR()
{
	int K = 2500;
	double dt = Tir / K;
	FILE* fp = fopen("IR.txt", "w");
	for (int i = 0; i < K; i++)
	{
		double t = i * dt;
		fprintf(fp, "%e  %e\n", t / units, F(t));
	}
	fclose(fp);
}

void writeFILE_maxHHG()
{
	FILE* fp;
	fp = fopen("maxHHG_XUV.txt", "w");
	for (int i = 0; i < nextIndex_max; i++)
	{
		fprintf(fp, "%e  %e\n", array_Wxuv[i], array_maxE[i] + Ip);
	}
	nextIndex_max = 0;
	fclose(fp);
}

void save_maxE()
{
	array_Wxuv[nextIndex_max] = Wxuv;
	array_maxE[nextIndex_max] = max_E;
	nextIndex_max++;
}
#pragma endregion

int work()
{
    #pragma region потенциал ионизации
	printf("Atom ionization energy (eV)\n");
	do
	{
		printf("  ");
		scanf("%lf", &Ip);
	} 
	while (Ip < 3 || Ip > 30);
	printf("\n\n");
    #pragma endregion

    #pragma region длина IR
	printf("Impulse length (fs)\n");
	do
	{
		printf("  ");
		scanf("%lf", &Tir);
	} 
	while (Tir * units < 3 * M_PI || Tir * units > 15 * M_PI);
	printf("\n\n");
	Tir *= units; // для перевода из фемтосекунд (fs) в безразмерные единицы (W * t, то есть 1/s * s = 1, s - секунда)
    #pragma endregion

    #pragma region определение константы C
	C = -(IntA(Tir) - IntA(0)); // 20: -0.1043234466544098, 7pi: 0
    #pragma endregion

    #pragma region построение траекторий
	max_E = 0;
	xuv_ir = 1;
	for (int i = 0;; i++)
	{
		double to1 = i * M_PI;
		if (to1 > Tir)
			break;

		printf("\n\nstart t1 = %.3e fs", i * M_PI / units);
		printf("\n------------------\n");
		printf("Count of maxs = %d", trajectories(to1));
	}
	writeFILE_IR();
	writeFILE_HHG("HHG_IR.txt");

	Wxuv = 0;
	for (int i = 0, count = 0; Wxuv < 80; i++)
	{
		Wxuv = rint(Ip) + i;

		printf("\n\nWxuv = %.3e eV", Wxuv);
		printf("\nstart t1 = %.3e fs", max_to1 / units);
		printf("\n------------------\n");

		xuv_ir = 2;
		max_E = 0;
		count = trajectories(max_to1);
		if (count > 0)
		{
			save_maxE();
		}
		printf("Count of maxs (+) = %d\n", count);
		

		xuv_ir = 3;
		max_E = 0;
		count = trajectories(max_to1);
		if (count > 0)
		{
			save_maxE();
		}
		printf("Count of maxs (-) = %d\n", count);

		if ((int)Wxuv % 5 == 0) 
		{
			char fileName[100];
			sprintf(fileName, "%s%d%s", "HHG_XUV_", (int)Wxuv, ".txt");
			writeFILE_HHG(fileName);
		}
		nextIndex = 0;
	}
	writeFILE_maxHHG();
    #pragma endregion
	
	return 0;
}

int main(void)
{
	do
	{
		work();
		printf("\n\n\n\n");
	} 
	while (1 == 1);
	return 0;
}