#ifndef MATRIZ
#define MATRIZ

# include <cassert>
# include <cstdlib>
# include <iostream>
using namespace std;

# include <cmath>
# include <iomanip>
# include <stdio.h>
# include <string>

// Desenvolvido por Daniel Gomes Ribeiro
// .
// .
// Classe de Matrizes
class Matriz
{
    public:
        // CONSTRUTORES DO OBJETO MATRIZ

        // Construtor Default
        Matriz();
        // Cria vetor
        Matriz (int);
        // Cria Matriz
        Matriz (int, int);
        // Cria Matriz a partir de uma lista estática
        Matriz (double *, int, int);
        // Cria um vetor a partir de uma lista estática
        Matriz (double *, int);


        // FUNÇÕES DE IMPRESSÃO, ACESSO E ALTERAÇÃO DE MATRIZES

        // Imprime uma matriz
        void imprime(int=4 , int=4);
        // acessa um elemento da matriz
        double saida(int, int);
        // Acessa um elemento do vetor
        double saida(int);
        // verifica a dimensão da matriz
        int dim(char);
        // altera um elemento específico da matriz
        void entrada(int, int, double);
        // altera um elemento específico do vetor
        void entrada(int, double);

        // Destruidor do objeto Matriz
        ~Matriz();

    private:
        double **elemento;
        int linha;
        int coluna;
};



#endif // MATRIZ
