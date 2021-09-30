#include "matriz.h"


// CLASSE DE MATRIZES
// Desenvolvida por Daniel Gomes Ribeiro
// ==================

// Construidor padrão
Matriz::Matriz()
{
    linha = 0;
    coluna = 0;
    elemento = NULL;
}

// Cria uma matriz quadrada de ordem definida pelo parâmetro tamanho
Matriz::Matriz(int tamanho)
{
    linha = 1;
    coluna = (tamanho>0 ? tamanho : 1);

    elemento = new double *[linha];

    if (elemento==NULL)
    {
        cout << "Problema de alocamento de memória. Abortando...";
        abort();
    }

    int i=0,j=0;

    while (i<linha)
    {
        elemento[i] = new double [coluna];
        if (elemento[i] == NULL)
        {
            cout << "Problema de alocamento de memória. Abortando...";
            abort();
        }
        i++;
    }

    i=0;
    while (i<linha)
    {
        j=0;
        while (j<coluna)
        {
            elemento[i][j] = 0.0;
            j++;
        }
        i++;
    }

}


Matriz::Matriz(int L, int C)
{
    linha  = (L>0 ? L : 1);
    coluna = (C>0 ? C : 1);

    elemento = new double *[linha];

    // Checa o alocamento de memória
    if (elemento==NULL)
    {
        cout << "Problema de alocamento de memória. Abortando...";
        abort();
    }

    // Parâmetros auxiliares
    int i=0,j=0;

    while (i<linha)
    {
        elemento[i] = new double [coluna];
        if (elemento[i] == NULL)
        {
            // Checa o alocamento de memória
            cout << "Problema de alocamento de memória. Abortando...";
            abort();
        }
        i++;
    }

    i=0;
    while (i<linha)
    {
        j=0;
        while (j<coluna)
        {
            elemento[i][j] = 0.0;
            j++;
        }
        i++;
    }
}


// Cria matriz a partir de uma lista 2d de valores
Matriz::Matriz(double *LISTA, int L, int C)
{
    linha  = (L>0 ? L : 1);
    coluna = (C>0 ? C : 1);

    elemento = new double *[linha];

    // Checa o alocamento de memória
    if (elemento==NULL)
    {
        cout << "Problema de alocamento de memória. Abortando...";
        abort();
    }

    // Parâmetros auxiliares
    int i=0, j=0;
    while (i<linha)
    {
        elemento[i] = new double [coluna];
        if (elemento[i] == NULL)
        {
            // Checa o alocamento de memória
            cout << "Problema de alocamento de memória. Abortando...";
            abort();
        }
        i++;
    }


    int contador = 0;
    i=0;
    while (i<linha)
    {
        j=0;
        while (j<coluna)
        {
            elemento[i][j] = static_cast <double> (LISTA[contador]);
            j++;
            contador+=1;
        }
        i++;
    }
}

// Cria matriz a partir de uma lista 2d de valores          
Matriz::Matriz(double *LISTA,int C)
{
    linha  = 1;
    coluna = (C>0 ? C : 1);

    elemento = new double *[linha];

    // Checa o alocamento de memória
    if (elemento==NULL)
    {
        cout << "Problema de alocamento de memória. Abortando...";
        abort();
    }

    // Par�metros auxiliares
    int i=0, j=0;
    while (i<linha)
    {
        elemento[i] = new double [coluna];
        if (elemento[i] == NULL)
        {
            // Checa o alocamento de memória
            cout << "Problema de alocamento de memória. Abortando...";
            abort();
        }
        i++;
    }


    int contador = 0;
    j=0;
    while (j<coluna)
    {
        elemento[0][j] = static_cast <double> (LISTA[contador]);
        j++;
        contador+=1;
    }
}



// Destruidor da Matriz
Matriz::~Matriz()
{
   for(int i=0; i<linha; i++)
   {
           delete[]elemento[i];
   }
   delete[]elemento;
};


void Matriz::imprime(int decimais, int space)
{
    int MaxMod = 0;
    int i = 0, j = 0, k=0, M=0;

    // Descobre a maior magnitude das componentes da matriz
    i = 0;
    while (i<linha)
    {
        j = 0;
        while (j<coluna)
        {

            M = int(  floor(log10( fabs(elemento[i][j]) ) )  );
            if ( M > MaxMod )
            {
                MaxMod = M;
            }

            j++;
        }
        i++;
    }


    // Imprime cada entrada do vetor
    i=0;
    while (i<linha)
    {
        j = 0;
        while (j < coluna)
        {
            // Magnitude do número atual
            M = int(  floor(log10( fabs(elemento[i][j]) ))   );

            // Imprime sinal de positivo, negativo ou espaço em branco
            if (elemento[i][j]>0.)
            {
                cout << "+";
            }
            else if (elemento[i][j]<0)
            {
                cout << "-";
            }
            else
            {
                cout << " ";
            }

            // Impressão dos números cujo módulo é maior do que 1
            if ( floor (fabs( elemento[i][j] )) > 1.)
            {
                k=0;
                while(k<(MaxMod-M))
                {
                    cout << "0";
                    k++;
                }
            }
            // impressão dos número menores que 1 (da forma 0.xxxx)
            else
            {
                k=0;
                while(k<MaxMod)
                {
                    cout << "0";
                    k++;
                }
            }

            // Imprime a quantidade de casas decimais pedidas, com precisão fixa
            if (elemento[i][j]>0)
            {
                cout << setprecision(decimais) << std::fixed << elemento[i][j];
            }
            else if (elemento[i][j]<0)
            {
                cout << setprecision(decimais) << std::fixed << -1.0*elemento[i][j];
            }
            else
            {
                cout << setprecision(decimais) << std::fixed << elemento[i][j];
            }


            k=0;
            while (k<space)
            {
                printf(" ");
                k++;
            }
            // Incremento de Coluna
            j++;
        }
        // Salta uma linha
        cout << endl;
        // Incremento de linnha
        i++;
    }

}

// Acessa dado interno à matriz
double Matriz::saida(int L, int C)
{
    if ( (L<linha) && (C < coluna) )
    {
        return elemento[L][C];
    }
    else
    {
        cout << "Erro de indice da matriz. Abortando...\n\n";
        abort();
    }
}


// Acessa dado interno à matriz - caso especial de quando a matriz é um vetor
double Matriz::saida(int POSITION)
{
    if (POSITION<coluna && linha == 1)
    {
        return elemento[0][POSITION];
    }
    else
    {
        cout << "Erro de indice da matriz. Abortando...\n\n";
        abort();
    }
}


// Altera dado interno à matriz
void Matriz::entrada(int L, int C, double IN)
{
    if (L<linha && C<coluna)
    {
        elemento[L][C] = IN;
    }
    else
    {
        cout << "Erro de indice da matriz. Abortando...\n\n";
        abort();
    }

}

// Altera dado interno ao vetor
void Matriz::entrada(int POSITION , double IN)
{
    if ( POSITION<coluna && linha == 1)
    {
        elemento[0][POSITION] = IN;
    }
    else
    {
        cout << "Erro de indice da matriz. Abortando...\n\n";
        abort();
    }
}


// Imprime a dimensão da Matriz
 int Matriz::dim( char LC)
 {
     if (LC=='l') { return linha;}
     else if (LC=='c') { return coluna; }
     else
     {
         cout << "Funcao dim deve receber 'l' (linha) ou 'c' (coluna) como parametro de entrada. Abortando...\n\n";
         abort();
     }
 }

 //


