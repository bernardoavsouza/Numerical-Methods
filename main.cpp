// Arquivo de programaçãoo contendo o main



// Inclusão de bibliotecas usadas em cada programa
#include <iostream>
using namespace std;

// Biblitoecas desenvolvidas
# include "matriz.h"
# include "minhabiblioteca.h"
# include "minhasfuncoes.h"


// Início da função main.


int main()
{

// Resolvendo uma regressão linear

int n=6;
double x[n]={2,3,4,5,6,7},r2;
double f[n]={1.4,1.75,2.6,3.25,3.7,4.4};

Matriz X(x,1,n), F(f,n,1), coeficientes(1,2);


MinimosQuad(X,F,1,coeficientes,r2);




coeficientes.imprime(7);

cout << "\n"<< r2 ;


}
