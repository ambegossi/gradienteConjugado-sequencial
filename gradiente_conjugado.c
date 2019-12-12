#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iohb.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>

void imprime(int l, int colunas, double *x, double *r, int *linhas, double *valores, int *colptr, int n)
{
    int i, coluna;

    printf("\nvetor x\n");
    for (i = 0; i < l; i++)
    {
        printf("%lf\n", x[i]);
    }
    printf("\n");

    for (i = 0; i < n; i++)
        r[i] = 0;

    coluna = -1;
    i = 0;
    while ((coluna + 1) < n)
    {
        if (i + 1 == colptr[coluna + 1])
            coluna++;
        r[linhas[i] - 1] += valores[i] * x[coluna];
        if ((linhas[i] - 1) != coluna)
        {
            r[coluna] += valores[i] * x[linhas[i] - 1];
        }
        i++;
    }

    printf("Prova Real (Ax = b) \n");
    for (i = 0; i < n; i++)
    {
        printf("%.6f\n", r[i]);
    }
    printf("\n");
}

void subtrai(int linhas, double *b, double *vet_aux, double *r)
{
    int i;

    for (i = 0; i < linhas; i++)
    {
        r[i] = b[i] - vet_aux[i];
    }
}

void adicao(int linhas, double *b, double *matriz_aux, double *r)
{
    int i;

    for (i = 0; i < linhas; i++)
    {
        r[i] = b[i] + matriz_aux[i];
    }
}

void le_matriz(char *arquivo, int *M, int *N, int *naozeros, int **colptr, int **linhas, double **valores)
{
    int retorno, nrhs;
    char *tipo;
    retorno = readHB_info(arquivo, M, N, naozeros, &tipo, &nrhs);

    if (retorno == 0)
    {
        printf("Erro ao ler as informaçõess da matriz!\n");
        exit(-1);
    }

    printf("Linhas: %d \t Colunas: %d \t Não Zeros: %d \n\n", *M, *N, *naozeros);

    *valores = (double *)malloc(*naozeros * sizeof(double));
    *linhas = (int *)malloc(*naozeros * sizeof(int));
    *colptr = (int *)malloc((*N + 1) * sizeof(int));

    retorno = readHB_mat_double(arquivo, *colptr, *linhas, *valores);

    if (retorno == 0)
    {
        printf("Erro ao ler os valores da matriz!\n");
        exit(-1);
    }

    free(tipo);
}

void escreve_matriz(int M, int N, int naozeros, int *colptr, int *linhas, double *valores)
{
    int i;

    printf("VALORES:\n");
    for (i = 0; i < naozeros; i++)
    {
        printf("%f ", valores[i]);
    }
    printf("\n\n");

    printf("LINHAS:\n");
    for (i = 0; i < naozeros; i++)
    {
        printf("%d ", linhas[i]);
    }
    printf("\n\n");

    printf("PTR:\n");
    for (i = 0; i < M + 1; i++)
    {
        printf("%d ", colptr[i]);
    }
    printf("\n");
}

void multiplicacaoMatrizVetor(int M, int N, int *colptr, int *linhas, double *a, double *b, double *q)
{

    int i, j;

    for (j = 0; j < N; j++)
    {
        q[j] = 0;
    }

    for (i = 0; i < N; i++)
    {
        for (j = colptr[i] - 1; j < colptr[i + 1] - 1; j++)
        {
            if (j != colptr[i] - 1)
            {
                q[linhas[j] - 1] += a[j] * b[i];
            }

            q[i] += a[j] * b[linhas[j] - 1];
        }
    }
}

void gradiente_conjugado(int M, int N, int naozeros, int *colptr, int *linhas, double *valores, double *b)
{
    double *aux, *x, *r, *d, *q, *alphad, *betad;
    aux = (double *)malloc(M * sizeof(double));
    x = (double *)malloc(M * sizeof(double));
    r = (double *)malloc(M * sizeof(double));
    d = (double *)malloc(M * sizeof(double));
    q = (double *)malloc(M * sizeof(double));
    alphad = (double *)malloc(M * sizeof(double));
    betad = (double *)malloc(M * sizeof(double));

    double sigma_novo = 0, sigma0 = 0, sigma_velho = 0, alpha = 0, beta = 0;

    int j, i = 1, imax = 5000000;
    double erro = 0.00001;

    // inicializa o vetor aux e o vetor x
    for (j = 0; j < N; j++)
    {
        aux[j] = 0;
        x[j] = 0;
    }

    // aux = A * x
    multiplicacaoMatrizVetor(M, N, colptr, linhas, valores, x, aux);

    printf("\nvetor b\n");
    for (j = 0; j < M; j++)
    {
        printf("%lf\n", b[j]);
    }

    // r = b - aux
    subtrai(M, b, aux, r);

    // d = r
    for (j = 0; j < M; j++)
    {
        d[j] = r[j];
    }

    // sigma_novo = r' * r
    for (j = 0; j < M; j++)
    {
        sigma_novo += r[j] * r[j];
    }

    // sigma0 = sigma_novo
    sigma0 = sigma_novo;

    double dlinhaq;

    while (i < imax && sigma_novo > (erro * erro * sigma0))
    {
        // zera o vetor q
        for (j = 0; j < N; j++)
        {
            q[j] = 0;
        }

        // q = A * d
        multiplicacaoMatrizVetor(M, N, colptr, linhas, valores, d, q);

        dlinhaq = 0;
        // (d' * q)
        for (j = 0; j < M; j++)
        {
            dlinhaq += d[j] * q[j];
        }
        // alpha = sigma_novo/(d' * q);
        alpha = sigma_novo / dlinhaq;

        // alpha * d;
        for (j = 0; j < M; j++)
        {
            alphad[j] = alpha * d[j];
        }

        // x = x + alpha * d;
        adicao(M, x, alphad, x);

        // zera o vetor aux
        for (j = 0; j < N; j++)
        {
            aux[j] = 0;
        }

        if (i % 50 == 0)
        {
            // aux = A * x
            multiplicacaoMatrizVetor(M, N, colptr, linhas, valores, x, aux);
            // r = b - matriz_aux
            subtrai(M, b, aux, r);
        }
        else
        {
            // aux = alpha * q
            for (j = 0; j < M; j++)
            {
                aux[j] = alpha * q[j];
            }
            // r = r - matriz_aux
            subtrai(M, r, aux, r);
        }

        // sigma_velho = sigma_novo
        sigma_velho = sigma_novo;
        //printf("sigma = %lf\n", sigma_velho);
        sigma_novo = 0;

        // sigma_novo = r' * r
        for (j = 0; j < M; j++)
        {
            sigma_novo += r[j] * r[j];
        }

        // beta = sigma_novo / sigma_velho
        beta = sigma_novo / sigma_velho;

        // beta * d
        for (j = 0; j < M; j++)
        {
            betad[j] = beta * d[j];
        }
        // d = r + beta * d
        adicao(M, r, betad, d);

        // i = i + 1
        i = i + 1;
    }

    //imprime(M, N, x, r, linhas, valores, colptr, N);
    printf("i = %d\n", i);

    free(aux);
    free(x);
    free(r);
    free(d);
    free(q);
    free(alphad);
    free(betad);
}

int main(int argc, char **argv)
{
    double *valores = NULL;
    int *linhas = NULL, *colptr = NULL;
    int M, N, naozeros;

    if (argc != 2)
    {
        printf("%s < Arquivo HB >\n", argv[0]);
        exit(-1);
    }

    le_matriz(argv[1], &M, &N, &naozeros, &colptr, &linhas, &valores);

    double *b;
    b = (double *)malloc(M * sizeof(double));
    // inicializa o vetor b que seria o recebido por parametro no gradiente do octave
    int i;
    for (i = 0; i < M; i++)
    {
        b[i] = 7;
    }

    escreve_matriz(M, N, naozeros, colptr, linhas, valores);

    double ti = omp_get_wtime();

    gradiente_conjugado(M, N, naozeros, colptr, linhas, valores, b);

    double tf = omp_get_wtime();

    printf("tempo gradiente: %lf\n", tf - ti);

    free(valores);
    free(linhas);
    free(colptr);
    free(b);
}
