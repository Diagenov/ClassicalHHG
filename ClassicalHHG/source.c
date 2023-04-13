#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

// либо это:
//   1) идти снизу вверх
//       + сверху вниз вторая ветка может пойти по тому же пути, что первая
//       + позволит не рассчитывать максимальную энергию, а определять ее по построению
//       - непонятно, как выбирать начальные условия
// либо это:
//   2) на каждом шаге давать приращение начальным условиям, а не только на первом

#pragma region параметры IR, XUV и атома
const double Up = 53.74;  // пондермоторная энергия, эВ
const double Ip = 21.55;  // потенциал ионизации (Ne) эВ
const double Wxuv = 30;   // частота XUV, эВ
const double Wir = 1;     // частота IR, эВ
const double Txuv = 0.55; // время XUV, фс
double Tir = 20;    // время IR, фс
double fi = 0 * M_PI_2; // относительная фаза огибающей
#pragma endregion

#pragma region выбор случая (монохромат, импульсы, IR, IR + XUV)
int monochromate = 1; // 1 - если монохромат, 2 - если импульсы
int xuv_ir = 1;       // 1 - если IR, 2 - если IR + XUV
#pragma endregion

#pragma region параметры программы
#define N 5     // количество максимумов
#define M 100   // количество точек на одну ветку
#define D 10000 // ограничение на запись в файл  
#pragma endregion

#pragma region для записи в один файл
double array_t[D];
double array_E[D];
int nextIndex = 0;
#pragma endregion

#pragma region параметры для системы уравнений или для расчета действия S
struct parameters
{
	double E; // энергия HHG (или электрона, зависит от записи)
	double t1; // время ионизации
	double t2; // время рекомбинации 
};
#pragma endregion

#pragma region функции, задающие IR-поле и импульс электрона
double C = -1;

double F(double t)
{
	if (monochromate == 1)
		return cos(t - fi);

	if (t < 0 || t > Tir)
		return 0;

	return cos(t - fi) * pow(sin(M_PI * t / Tir), 2);
}

double A(double t)
{
	if (monochromate == 1)
		return -sin(t - fi);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	double y = (sin(t - fi) / 2) - (sin(a1 * t - fi) / (4 * a1)) - (sin(a2 * t + fi) / (4 * a2)) + C;
	return -y;
}

double IntA(double t, double c)
{
	if (monochromate == 1)
		return cos(t - fi);

	if (t < 0 || t > Tir)
		return 0;

	double a1 = (2 * M_PI / Tir) + 1;
	double a2 = (2 * M_PI / Tir) - 1;
	return (cos(t - fi) / 2) - (cos(a1 * t - fi) / (4 * a1 * a1)) - (cos(a2 * t + fi) / (4 * a2 * a2)) + (c * t);
}

double Q(double t1, double t2)
{
	if (t1 == t2)
		return 0;

	return -(IntA(t2, C) - IntA(t1, C)) / (t2 - t1);
}

double P(double t, double t1, double t2)
{
	return A(t) + Q(t1, t2);
}

double deltaE(double t1, double t2)
{
	double sqrK = 2 * Ip;
	double a = sqrK / (t2 - t1);
	return a * P(t2, t1, t2) / (2 * Wir * F(t1));
}
#pragma endregion

#pragma region для спектра
double sqrP(double t, void* params)
{
	struct parameters* p = (struct parameters*)params;
	const double t1 = (p->t1);
	const double t2 = (p->t2);

	return 2 * Up * pow(P(t, t1, t2), 2);
}

double S(double t1, double t2)
{
	gsl_integration_workspace* w
		= gsl_integration_workspace_alloc(1000);

	double result, error;
	struct parameters p = { 0, t1, t2 };

	gsl_function F;
	F.function = &sqrP;
	F.params = &p;

	gsl_integration_qags(&F, t1, t2, 0, 1e-7, 1000,
		w, &result, &error);

	gsl_integration_workspace_free(w);
	return (Ip * (t2 - t1)) + result;
}

double g_propagation(double OMEGA, double t1, double t2)
{
	double argument = (OMEGA * t2) - S(t1, t2);
	double ReZ = cos(argument);
	double ImZ = sin(argument);

	return ReZ / sqrt(pow(t2 - t1, 3));
}
#pragma endregion

#pragma region максимальная гармоника
double maxHHG(double t1, double t2)
{
	double deltaHHG = xuv_ir == 1 ? 0 : F(t2) * (Wxuv - Ip) / F(t1);
	return (2 * Up * pow(A(t2) - A(t1), 2)) + deltaHHG;
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
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - (E);
	}
	else
	{
		e1 = (2 * Up * pow(P(t1, t1, t2), 2)) - (Wxuv - Ip);
		e2 = (2 * Up * pow(P(t2, t1, t2), 2)) - (E);
	}
	gsl_vector_set(func, 0, e1);
	gsl_vector_set(func, 1, e2);

	return GSL_SUCCESS;
}
#pragma endregion

#pragma region поиск максимумов
int max_t(double to1, double max_t1[N][2], double max_t2[N][2], double max_e[N])
{
	if (to1 > Tir)
		return 0;

	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver * s;

	int status;
	int count = 0;
	size_t iter = 0;
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
		double to2 = to1 + ((2 + i) * M_PI) - M_PI_2;
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

			double hhg = maxHHG(t1, t2);
			if (hhg < 5)
				continue;

			max_t1[count][0] = max_t1[count][1] = t1;
			max_t2[count][0] = t2 - 0.3;
			max_t2[count][1] = t2 + 0.3;
			max_e[count] = hhg;

			printf("t1 = %.3e \n", t1);
			printf("t2 = %.3e \n", t2);
			printf("t = %.3e \n", t2 - t1);
			if (xuv_ir == 1)
			{
				printf("maxHHG = %.3e Up = %.3e \n\n", 2 * pow(A(t2) - A(t1), 2), hhg);
			}
			else
			{
				printf("maxHHG = %.3e Up  +  %.3e (Wxuv - Ip) = %.3e \n\n", 2 * pow(A(t2) - A(t1), 2), (F(t2) / F(t1)), hhg);
			}
			count++;
		}
	}
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return count;
}
#pragma endregion

#pragma region построение траекторий
int trajectories(double to1)
{
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
		array_t[nextIndex] = max_t2[i][0] + 0.3;
		array_E[nextIndex] = max_e[i];
		nextIndex++;

		for (int j = 0; j < 2; j++)
		{
			for (int E = (int)max_e[i]; E > 0; E--)
			{
				struct parameters p = { E, 0, 0 };
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

				if (status != GSL_SUCCESS)
					continue;
				
				double t1 = gsl_vector_get(s->x, 0);
				double t2 = gsl_vector_get(s->x, 1);

				if (abs(max_t2[i][j] - t2) > M_PI_4)
					break;

				max_t1[i][j] = t1;
				max_t2[i][j] = t2;

				array_t[nextIndex] = t2;
				array_E[nextIndex] = E + Ip;
				nextIndex++;
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
void writeFILE()
{
    #pragma region траектории
	FILE* fp;
	fp = fopen("HHG.txt", "w");
	for (int i = 0; i < nextIndex; i++)
	{
		fprintf(fp, "%e  %e\n", array_t[i], array_E[i]);
	}
	fclose(fp);
    #pragma endregion

	int K = 1000;
	double dt = Tir / K;

    #pragma region IR-напряженность
	fp = fopen("F.txt", "w");
	for (int i = 0; i < K; i++)
	{
		double t = i * dt;
		fprintf(fp, "%e  %e\n", t, 35 * F(t));
	}
	fclose(fp);
    #pragma endregion

    #pragma region IR-потенциал
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
    #pragma region импульсы или монохромат
	printf("1 - monochromate, 2 - impulses\n");
	do
	{
		printf("  ");
		scanf("%d", &monochromate);
	} 
	while (monochromate != 1 && monochromate != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region длина импульса
	if (monochromate == 2)
	{
		printf("Impulse length\n");
		do
		{
			printf("  ");
			scanf("%lf", &Tir);
		} 
		while (Tir < 0 || Tir > 30);
		printf("\n\n");
	}
    #pragma endregion

    #pragma region IR или IR + XUV
	printf("1 - IR, 2 - IR + XUV\n");
	do
	{
		printf("  ");
		scanf("%d", &xuv_ir);
	} 
	while (xuv_ir != 1 && xuv_ir != 2);
	printf("\n\n");
    #pragma endregion

    #pragma region момент ионизации
	int n = -1;
	if (xuv_ir == 2)
	{
		printf("Ionization time (N * PI)\n");
		do
		{
			printf("  N = ");
			scanf("%d", &n);
		} 
		while (n < 0 || n > 10);
		printf("\n\n");
	}
    #pragma endregion

    #pragma region относительная фаза огибающей
	printf("Carrier-envelope phase\n");
	do
	{
		printf("  ");
		scanf("%lf", &fi);
	} 
	while (fi < 0 || fi > 2 * M_PI);
	printf("\n\n");
    #pragma endregion

    #pragma region определение константы C
	for (int i = 0; i < 2000; i++)
	{
		double c = (i * 2.0 / 2000) - 1;
		if (abs(IntA(Tir, C) - IntA(0, C)) > abs(IntA(Tir, c) - IntA(0, c)))
		{
			C = c;
		}
	}
	//C = 0.60556;
    #pragma endregion

    #pragma region параметр Келдыша (модифицированный, XUV + IR)
	const double k = xuv_ir == 1 ? 0 : sqrt((Wxuv - Ip) / (2 * Up)); 
    #pragma endregion

    #pragma region построение траекторий
	for (int i = n < 0 ? 0 : n; n < 0 || i < n + 1; i++)
	{
		double to1 = fi + (i * M_PI);
		if (to1 > Tir)
			break;

		printf("\n\nstart t1 = %d * PI\n", i);
		printf("------------------\n", i);
		if (trajectories(to1) == 0)
		{
			printf("nothing\n");
		}
	}
	writeFILE();
	printf("\n\n");
	return 0;
    #pragma endregion
}