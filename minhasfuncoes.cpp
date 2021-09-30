# include "minhasfuncoes.h"

// .
// .
// .
// .
// .
// ========================================================================================

// R1 em R1
double funcaopadrao(double x)
{
   return 1/(x*x*x);
}


// Função Rn em R1
double funcaopadrao( Matriz &X )
{
    return 2*pow(X.saida(0,0),2)*X.saida((0,1));
}

// Função de Rn em Rn
void funcaopadrao(Matriz &X, Matriz &OUT)
{
    double temp, x,y,z,w;

    x= X.saida(0,0);
    y= X.saida(0,1);
    z= X.saida(0,2);
    w= X.saida(0,3);

    temp = 2*x*x+5*y-log(z*w)+exp(w);
    OUT.entrada(0, 0, temp);

    temp = x*y*z+cos(w);
    OUT.entrada(0, 1, temp);

    temp =  x+y+z-w;
    OUT.entrada(0, 2, temp);

    temp =  x*x+y*y+z*z+w*w;
    OUT.entrada(0, 3, temp);



}



// Função de elevar a 2
double xquadrado (double x)
{
    return x*x;
}

// Função de elevar a 3
double xcubo (double x)
{
    return x*x*x;
}

// Função que gera uma parábola
double parabol(Matriz X)
{
    return pow( X.saida(0,0) , 2) + pow(X.saida(0,1),2);
}

// Função que calcula fatorial
int fatorial (double n)
{
    if (n<0)
    {
        cout << "Numero negativo. Funcao fatorial. Abortando...";
        abort();
    }

    if(floor(n)!= ceil(n))

    {
        cout << "Numero nao inteiro. Funcao fatorial. Abortando...";
        abort();
    }

    if(n==0)
    {
        return 1;
    }

    int k, temp=1;
    for (k=1 ; k<=n ; k++)
    {
        temp = temp*k;

    }
    return temp;
}

// .
// .
// .
// .
// .