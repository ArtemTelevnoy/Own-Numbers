#include "return_codes.h"

#include <malloc.h>
#include <math.h>
#include <stdio.h>

void matrixMultiply(size_t n, size_t m, double *eps, double *a, double *b, double *d)
{
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			*(d + i * n + j) = 0;
			for (size_t k = 0; k < m; k++)
			{
				*(d + i * n + j) += *(a + i * n + k) * *(b + k * n + j);
			}
			if (fabs(*(d + i * n + j)) < *eps)
			{
				*(d + i * n + j) = 0;
			}
		}
	}
}

double len(size_t n, size_t high, double *a)
{
	double len = 0;
	for (size_t i = 0; i < n; i++)
	{
		len += pow(*(a + i * n + high), 2);
	}
	return sqrt(len);
}

void qrSolve(size_t n, double *eps, double *a, double *q, double *r, double *aj, double *qk)
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			*(q + i * n + j) = *(r + i * n + j) = 0;
		}
	}
	for (size_t j = 0; j < n; j++)
	{
		for (size_t i = 0; i < n; i++)
		{
			*(aj + i) = *(a + i * n + j);
		}

		for (size_t k = 0; k < j; k++)
		{
			for (size_t i = 0; i < n; i++)
			{
				*(r + k * n + j) += *(a + i * n + j) * *(q + i * n + k);
				*(qk + i) = *(q + i * n + k);
			}

			for (size_t i = 0; i < n; i++)
			{
				*(aj + i) -= *(qk + i) *= *(r + k * n + j);
			}
		}

		for (size_t k = 0; k < n; k++)
		{
			*(q + k * n + j) = *(aj + k);
		}

		double no = *(r + j * n + j) = len(n, j, q);

		for (size_t i = 0; i < n; i++)
		{
			*(q + i * n + j) = fabs(no) < *eps ? 0 : *(q + i * n + j) / no;
		}
	}
	matrixMultiply(n, n, eps, r, q, a);
}

void square(double *a, double *b)
{
	double be = *(a) + *(a + 3);
	double dis = pow(be, 2) - 4 * (*(a) * *(a + 3) - *(a + 1) * *(a + 2));
	*(b) = be / 2 + 0;
	*(b + 1) = sqrt(-dis) / 2;
	*(b + 2) = -*(b + 1) + 0;
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Not enough arguments");
		return ERROR_PARAMETER_INVALID;
	}
	FILE *fp;
	if (!(fp = fopen(argv[1], "r")))
	{
		fprintf(stderr, "File didn't open");
		//fclose(fp);
		return ERROR_CANNOT_OPEN_FILE;
	}
	size_t n;
	if (fscanf(fp, "%Iu", &n) != 1)
	{
		fprintf(stderr, "Didn't scan");
		fclose(fp);
		return ERROR_UNKNOWN;
	}

	double *a = (double *)malloc(pow(n, 2) * sizeof(double *));
	if (!a)
	{
		fprintf(stderr, "Not enough memory");
		fclose(fp);
		return ERROR_OUT_OF_MEMORY;
	}

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (fscanf(fp, "%lf", a + i * n + j) != 1)
			{
				fprintf(stderr, "Didn't scan");
				fclose(fp);
				free(a);
				return ERROR_UNKNOWN;
			}
		}
	}
	fclose(fp);

	double *aj = (double *)malloc(n * sizeof(double *));
	if (!aj)
	{
		free(a);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}
	double *qk = (double *)malloc(n * sizeof(double *));
	if (!qk)
	{
		free(a);
		free(aj);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}
	double *q = (double *)malloc(pow(n, 2) * sizeof(double *));
	if (!q)
	{
		free(a);
		free(aj);
		free(qk);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}
	double *r = (double *)malloc(pow(n, 2) * sizeof(double *));
	if (!r)
	{
		free(a);
		free(aj);
		free(qk);
		free(q);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}

	double *eps;
	double help21 = pow(0.1, pow(n, 2));
	eps = &help21;

	for (size_t i = 0; i < 10000; i++)
	{
		qrSolve(n, eps, a, q, r, aj, qk);
	}
	free(aj);
	free(qk);
	free(q);
	free(r);

	double *mac = (double *)malloc(4 * sizeof(double *));
	if (!mac)
	{
		free(a);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}
	double *mac2 = (double *)malloc(4 * sizeof(double *));
	if (!mac2)
	{
		free(a);
		free(mac);
		fprintf(stderr, "Not enough memory");
		return ERROR_OUT_OF_MEMORY;
	}

	FILE *wr;
	if (!(wr = fopen(argv[2], "w")))
	{
		fprintf(stderr, "File didn't create");
		//fclose(wr);
		free(a);
		return ERROR_CANNOT_OPEN_FILE;
	}

	for (size_t i = 0; i < n; i++)
	{
		if (i == n - 1 || *(a + (i + 1) * n + i) == 0)
		{
			if (fprintf(wr, "%g\n", *(a + i * n + i)) < 0)
			{
				fprintf(stderr, "Didn't write");
				free(a);
				free(mac);
				free(mac2);
				fclose(wr);
				return ERROR_OUT_OF_MEMORY;
			}
		}
		else
		{
			*(mac) = *(a + i * n + i);
			*(mac + 1) = *(a + i * n + i + 1);
			*(mac + 2) = *(a + (i + 1) * n + i);
			*(mac + 3) = *(a + (i + 1) * n + i + 1);
			square(mac, mac2);
			if (fprintf(wr, "%g%s%gi\n", *(mac2), *(mac2 + 1) > 0 ? " +" : " ", *(mac2 + 1)) < 0)
			{
				fprintf(stderr, "Didn't write");
				free(a);
				free(mac);
				free(mac2);
				fclose(wr);
				return ERROR_UNKNOWN;
			}
			if (fprintf(wr, "%g%s%gi\n", *(mac2), *(mac2 + 2) > 0 ? " +" : " ", *(mac2 + 2)) < 0)
			{
				fprintf(stderr, "Didn't write");
				free(a);
				free(mac);
				free(mac2);
				fclose(wr);
				return ERROR_OUT_OF_MEMORY;
			}
			i++;
		}
	}
	fclose(wr);
	free(a);
	free(mac);
	free(mac2);
	return SUCCESS;
}