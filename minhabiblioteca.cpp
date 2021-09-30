
// inclusão do cabeçalho que aponta para esse arquivo
#include <iostream>
using namespace std;

# include "minhabiblioteca.h"
# include "matriz.h"
# include "minhasfuncoes.h"
# include <math.h>
// .
// .
// .
// .
// .

// LISTA DE FUNÇÕES A SEREM CRIADAS:

// Cria uma matriz identidade
void Identidade (Matriz &ID)

{

    int L, C, i;

        L = ID.dim('l');
        C = ID.dim('c');

   //Conferindo se a matriz é quadrada

   if(L!=C)
   {
       cout << "A matriz nao e quadrada. Abortando...";
       abort();
   }

   for(i=0;i<L;i++)
   {
        ID.entrada(i,i,1);

   }
}


// Copia a Matriz A para a matriz B.
void Copia (Matriz &A , Matriz &B)
{

    int LA, LB, CA, CB, i, j;

    //Descobre o número de linhas e colunas das matrizes

    LA=A.dim('l');
    LB=B.dim('l');
    CA=A.dim('c');
    CB=B.dim('c');

    if (LA!=LB || CA!=CB)
    {
        // Checa se linhas e colunas sao iguais
        cout << "As matrizes não tem dimensões iguais. Abortando...";
        abort();
    }
    for(i=0;i<LA;i++)
    {
        for(j=0;j<CA;j++)
        {
            B.entrada(i,j,A.saida(i,j));
        }
    }
}

// Soma duas matrizes
void SOMA(Matriz &A , Matriz &B, Matriz &C)
{
     int LA, LB, CA, CB, i, j, LC, CC;

    //Descobre o número de linhas e colunas das matrizes

    LA=A.dim('l');
    LB=B.dim('l');
    LC=C.dim('l');
    CA=A.dim('c');
    CB=B.dim('c');
    CC=C.dim('c');

    if (LA!=LB || CA!=CB)
    {
        cout << "As matrizes não tem dimensões iguais. Abortando...";
        abort();
    }

    for(i=0;i<LA;i++)
    {
        for(j=0;j<CA;j++)
        {
        C.entrada(i,j,A.saida(i,j)+ B.saida(i,j));
        }

    }
}

// Subtrai duas matrizes
void SUBTRAI(Matriz &A , Matriz &B, Matriz &C)
{
     int LA, LB, CA, CB, i, j, LC, CC;


    // Descobre o número de linhas e colunas das matrizes

    LA=A.dim('l');
    LB=B.dim('l');
    LC=C.dim('l');
    CA=A.dim('c');
    CB=B.dim('c');
    CC=C.dim('c');

    if (LA!=LB || CA!=CB)
    {
        cout << "As matrizes não tem dimensões iguais. Abortando...";
        abort();
    }

    for(i=0;i<LA;i++)
    {
        for(j=0;j<CA;j++)
        {
        C.entrada(i,j,A.saida(i,j)- B.saida(i,j));
        }

    }
}

// Multiplica duas matrizes
void MULTIPLICA(Matriz &A , Matriz &B , Matriz &C )
{
    // Checar se as dimensões SÃO COMPATÍVEIS
    int LA, CA, LB, CB, LC, CC;

    LA = A.dim('l');
    CA = A.dim('c');

    LB = B.dim('l');
    CB = B.dim('c');

    LC = C.dim('l');
    CC = C.dim('c');

    if ((CA!=LB) || (LC!=LA) || (CC!=CB))
    {
        cout << "Inconsistencia nas dimensoes das matrizes. Funcao, MULTIPLICA. Abortando...";
        abort();
    }

    int i,j,k;
    double somatorio;

    for (i=0;i<LA;i++)
    {

        for (j=0;j<CB;j++)
        {
            somatorio=0;
            for (k=0;k<CA;k++)
            {
                somatorio=somatorio+ A.saida(i,k)*B.saida(k,j);

            }

            C.entrada(i,j,somatorio);
        }
    }
}

// Transpõe a matriz
void TRANSPOSTA (Matriz &A, Matriz &B)
{
    if(A.dim('l')!= B.dim('c') || A.dim('c')!= B.dim('l'))
    {
        cout << "Erro de dimensão de matriz. Algoritmo TRANSPOSTA. Abortando...";
        abort();
    }
    int i,j;
    int n = A.dim('l');
    int m = A.dim('c');

    for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            B.entrada(j,i,A.saida(i,j));
        }
    }
}

// Série de Taylor da função e^(-2x)
double e2xTAYLOR (double x , double n)
{    

   if (n<=1 || floor(n)!= ceil(n))

   {
      cout << "Valor invalido de termos. Funcao e2xTAYLOR Abortando...";
      abort();
   }

   double i;
   double somat = 0;

   for (i=0; i<=n; i++)

   {
       somat = somat + ((pow(-2,i))/(fatorial(i)))*pow(x,i);
   }

   return somat;

}

// Série de Taylor da função e^(-ax)
double eaxTAYLOR (double x, double x0, double a, double n)
{
     if (n<=1 || floor(n)!= ceil(n))
   {
      cout << "Valor invalido de termos. Funcao eaxTAYLOR Abortando...";
      abort();
   }
    if(a>=0)
    {
        cout << "\nValor invalido de a. Funcao eaxTAYLOR. Abortando...";
        abort();
    }
    if(n>30)
    {
        cout << "\nValor de n excedeu o limite. Funcao eaxTAYLOR. Abortando...\n";
        abort();
    }

    int i;
    double somatorio = exp(a*x0);

    for(i=1; i<=n ; i++)
    {
        somatorio = somatorio+ (pow(a,i)*exp(a*x0)*pow(x-x0,i))/fatorial(i);

    }

    return somatorio;
}

double lnxTAYLOR (double x, double x0, double n)
{
    if (n<=1 || floor(n)!= ceil(n))
   {
      cout << "Valor invalido de termos. Funcao eaxTAYLOR Abortando...";
      abort();
   }
    if(x0<=0)
    {
        cout << "\nValor invalido de x0. Funcao lnxTAYLOR. Abortando...\n";
        abort();
    }
    if(n>30)
    {
        cout << "\nValor de n excedeu o limite. Funcao lnxTAYLOR. Abortando...\n";
        abort();
    }

    int i;
    double somatorio = log(x0);

    for(i=1; i<=n ; i++)
    {
        somatorio = somatorio+ (pow(-1,i+1))*pow(x-x0,i)/(i*pow(x0,i));
    }

    return somatorio;
}

double sinxTAYLOR (double x, double x0, double n)

{

    if (n<=1 || floor(n)!= ceil(n))
   {
      cout << "Valor invalido de termos. Funcao sinxTAYLOR Abortando...";
      abort();
   }

    if(n>30)
    {
        cout << "\nValor de n excedeu o limite. Funcao sinxTAYLOR. Abortando...\n";
        abort();
    }

    int i;
    double somatorio = sin(x0);

    for(i=1; i<=n; i++)
    {
        somatorio = somatorio + ( (pow(-1,i/2)*pow(sin(x0),((i+1)%2))*pow(cos(x0),i%2)*pow(x-x0,i) )/ fatorial(i) );
    }

    return somatorio;
}

double derivPROG (double x, double h)
{
    return (funcaopadrao(x+h) -funcaopadrao(x))/(h);
}

double derivREG (double x, double h)
{
    return (funcaopadrao(x) - funcaopadrao(x-h))/(h);
}
double derivCENT (double x, double h)
{
    return (funcaopadrao(x+h) - funcaopadrao(x-h))/(2*h);
}

double derivPARC (Matriz X0, int i, double h)
{
    int dim= X0.dim('c');
    Matriz H(1,dim), X0h(1,dim);
    H.entrada(0,i,h);
    double fh, f;
    SOMA(X0,H, X0h);
    fh = funcaopadrao(X0h);
    f = funcaopadrao(X0);

    return (fh-f)/h;
}

//Resolve um sistea triangular superior
void ResolTriangSup (Matriz A, Matriz B, Matriz &X)
{
    // Verificar se A é quadrada
    //Verificar se as dimensões de A e B são compatíveis
    //Verificar se as dimensões de B e X são compatíveis

    int i, j;
    int n=A.dim('c');
    double bn_nn, SOMAT, temp;

    bn_nn = B.saida(0,n-1)/A.saida(n-1,n-1);
    X.entrada(0,n-1, bn_nn);

    // for que percorre as linhas de baixo para cima
    for(i=n-2;i>=0;i--)
    {
        //Resolve o somatório
        SOMAT=0;
        for(j=i+1;j<=n-1;j++)
        {
            SOMAT = SOMAT + A.saida(i,j)*X.saida(0,j);
        }
        temp = (B.saida(0,i) - SOMAT)/A.saida(i,i);
        X.entrada(0,i,temp);
    }

}

//Realiza eleminação Gaussiana
void EliminGauss( Matriz A, Matriz B, Matriz &AA, Matriz &BB)
{
    //verificar se A é quadrada
    //Verificar se as dimensoes de a e b sao compativeis
    // verificar se as dimensoes  de A e  AA sao iguais
    //Verificar se as diumensoes de bb e b sao iguais

    int i,j,k,n;
    double m,temp_AA, temp_bb;

    Copia(A,AA);
    Copia(B,BB);

    n=B.dim('c');

    for(k=0;k<=n-2;k++)
    {
         for(i=k+1;i<=n-1;i++)
         {
             m = AA.saida(i,k)/ AA.saida(k,k);
             AA.entrada(i,k,0);
              for(j=k+1;j<=n-1;j++)
              {
                  temp_AA=AA.saida(i,j)- m*AA.saida(k,j);
                  AA.entrada(i,j, temp_AA);

                  temp_bb=BB.saida(0,i)- m*BB.saida(0,k);

              }
              BB.entrada(0,i,temp_bb);
         }
    }
}

// Eliminação gaussiana com pivoteamento parcial
void EliminGaussPivotParc( Matriz A, Matriz B, Matriz &AA, Matriz &BB, int DBG_VIEW)
{
    // Verifica se A é quadrada
    // Verifica se as dimensoes de A e B sao compativeis
    // Verifica se as dimensoes  de A e  AA sao iguais
    // Verifica se as diumensoes de BB e B sao iguais
    if((A.dim('c')!= A.dim('l')) || (A.dim('c')!= B.dim('c')) || (A.dim('c')!= AA.dim('c')) || (B.dim('c')!= BB.dim('c')) || (A.dim('l')!= AA.dim('l') || (B.dim('l')!= BB.dim('l'))))
    {
        cout << "Erro de dimensao das matrizes. Funcao EliminGaussPivotParc. Abortando...";
        abort();
    }

    int i,j,k,n,z,y;
    double m,temp_AA, temp_bb, pivo, dt;


    Copia(A,AA);
    Copia(B,BB);

    n=B.dim('c');
    Matriz VT (1,n);
    for(k=0;k<=n-2;k++)
    {

        pivo = AA.saida(k,k); // Valor do pivô padrão

        for(z=k+1;z<=n-1;z++) //
        {

            // Decorre a matriz
            if(abs(pivo) < abs(AA.saida(z,k)))
            {
                pivo = AA.saida(z,k); // Pivo novo escolhido
                for(y=0; y<=n-1;y++)
                {
                    VT.entrada(0,y,AA.saida(z,y)); // Vetor temporario contendo a linha do pivo copiada
                }

                dt = BB.saida(0,z); // Valor temporario de B na posição z copiado

                for(y=0; y<=n-1;y++)
                {
                        // Inserir a linha inicial na linha do futuro pivô
                    AA.entrada(z,y,AA.saida(k,y));
                }
                        // Inserir o valor de BB (posição z) na posição k
                BB.entrada(0,z,BB.saida(0,k));

                    for (y=0;y<=n-1;y++)
                    {
                        // Retorna o backup da linha z na linha k
                        AA.entrada(k,y,VT.saida(0,y));
                    }
                    // Retorna o valor do coeficiente de z em k
                BB.entrada(0,k,dt);
                if(DBG_VIEW ==1)
                {
                    cout << "\n";
                    AA.imprime();
                    cout << "\n";
                    BB.imprime();
                    cout << "\n";
                    cout << "" << pivo;
                    cout << "\n\n";

                }
            }
        }


         for(i=k+1;i<=n-1;i++)
         {
             m= AA.saida(i,k)/ AA.saida(k,k);
             AA.entrada(i,k,0);
              for(j=k+1;j<=n-1;j++)
              {
                  temp_AA=AA.saida(i,j)- m*AA.saida(k,j);
                  AA.entrada(i,j, temp_AA);

                  temp_bb=BB.saida(0,i)- m*BB.saida(0,k);

              }
              BB.entrada(0,i,temp_bb);
         }
    }
}

// Método de Jacobi para a resolução de equações
void Jacobi (Matriz A, Matriz b, Matriz X0, Matriz &X, double tol, int maximo)
{

    // Checa se o elemento da diagonal, em MODULO, é o maior elemento de sua linha (dominancia da diagonal).
    // Checa consistencia do problema (A tem que ser quadrada).
    // Checa se as dimensões de A, b, X0 e X são compatíveis.

    int n = A.dim('c'), k=0, i=0, j=0, m;
    double somatorio, temp, MAXI;

    Matriz Xk(1,n), Xk1(1,n);

    Copia(X0,Xk);

    // Contador de iterações
    for (k=0;k<maximo;k++)
    {

        // Incógnitas a serem calculadas
        for (i=0;i<n;i++) // for de linha
        {
            somatorio = 0;
            for (j=0;j<n;j++) //for de coluna
            {
                if (j!=i)
                {
                    somatorio = somatorio + A.saida (i,j)*Xk.saida (0,j); // Faz somatorio dos aij*xj (com exceção dos dos termos i=j)
                }
            }

            temp = (b.saida(0,i) - somatorio)/A.saida(i,i); //
            Xk1.entrada(0,i,temp);
        }
        // Criterio de parada
        // Descobrir o maior dos erros
        MAXI = abs(Xk1.saida(0,0) - Xk.saida(0,0)); // MAXI padr�o
        for (m=1 ; m<n ; m++)
        {
            if (abs(Xk1.saida(0,m) - Xk.saida(0,m)) > MAXI)
            {
                MAXI = abs(Xk1.saida(0,m) - Xk.saida(0,m));
            }
        }
        // Se o maior dos erros for menor que a tolerância, parar o loop de iterações
        if(MAXI < tol)
        {
            break;
        }

        Copia(Xk1,Xk);

    }
    if(k == maximo)
    {
        cout << "Nao houve convergencia. Algoritmo Jacobi\n";
    }

    Copia(Xk1,X);

}

// Método SOR de resolução de equações
void SOR (Matriz A, Matriz B, Matriz X0, Matriz &X, double w,double tol, double maximo, int DBG_VIEW)
{
    int n = X0.dim('c');
    int k, i, j,m;
    Matriz Xk(1,n), Xk1(1,n);
    Copia(X0, Xk);
    double soma1,soma2,temp,MAXI;

    for(k=1;k<maximo;k++)
    {
        for(i=1;i<n;i++)
        {
            soma1=0;
            soma2=0;
            for(j=0;j<i-1;j++)
            {
                soma1 = soma1 + A.saida(i,j)*Xk1.saida(0,j);

            }

            for(j=i+1;j<n;j++)
            {
                soma2 = soma2 + A.saida(i,j)*Xk.saida(0,j);

            }
            temp = (1-w)*Xk.saida(0,i) + (w*(B.saida(0,i)-soma1-soma2))/A.saida(i,i);
            Xk1.entrada(0,i, temp);
        }
        MAXI = abs(Xk1.saida(0,0) - Xk.saida(0,0)); // MAXI padrão
        for (m=1 ; m<n ; m++)
        {
            if (abs(Xk1.saida(0,m) - Xk.saida(0,m)) > MAXI)
            {
                MAXI = abs(Xk1.saida(0,m) - Xk.saida(0,m));
            }
        }
        // Se o maior dos erros for menor que a tolerância, parar o loop de iterações
        if(MAXI < tol)
        {
            break;
        }
    }
        if(k == maximo)
    {
        cout << "Nao houve convergencia. Algoritmo SOR";
    }
    Copia(Xk1,X);
}

// Método da bisseção para a resolução de equações
double Bissecao (double a_, double b_, int maximo, double tol, int DGB_VIEW)
{
    //maximo e tol> 0  e f(a)*f(b) <0
    int k;
    double a = a_, b = b_ , c;

    for(k=0;k<maximo;k++)
    {
        c=(a+b)/2;

        if(funcaopadrao(a)*funcaopadrao(c)<0)
        {
            b=c;
        }
        else
        {
            a=c;
        }
        if(abs(b-a)<tol)
        {
            return (a+b)/2;
        }
    }
        cout << "Nao houve convergencia. Algoritmo BISSECAO.\n";
        return (a+b)/2;
}

// Método de Newton para a resolução de equacões de uma variável
double MetNewton1 (double x0, int maximo, int DBG_VIEW, double e1, double e2)
{
    int k;
    double Xk, Xk1;
    Xk=x0;
    for(k=1;k<maximo;)
    {
        if(abs(derivCENT(Xk))<e1||abs(derivCENT(Xk))<e2)
        {
            cout << "Derivada proxima ou igual a zero. Funcao METNEWTON.\n";
            abort();
        }
        Xk1= Xk-funcaopadrao(Xk)/derivCENT(Xk);
        if(DBG_VIEW==1)
        {
            cout << Xk1;
            cout << "\n";
        }
        if(abs(funcaopadrao(Xk1))<e1)
        {
            return Xk1;
        }
        if(abs(Xk1-Xk)<e2)
        {
            return Xk1;
        }


        Xk=Xk1;
    }
    cout << "Nao convergiu. Algoritmo METNEWTON1.\n";
    return Xk1;
}

// Método da secante para a resolução de equações 
double MetSec (double x0,double x1,int maximo, int DBG_VIEW, double e1, double e2)
{

    if(abs(x0-x1)>e1 || abs(x0-x1)>e2)
    {
        cout <<"Valores iniciais muito proximos. Funcao METSEC.\n";
        abort();
    }
    //|x0-x1| > tol
    int k;
    double Xk0, Xk, Xk1;
    Xk=x1;
    Xk0=x0;
    for(k=1;k<maximo;)
    {
        Xk1= (Xk0*funcaopadrao(Xk)-Xk*funcaopadrao(Xk0))/(funcaopadrao(Xk)-funcaopadrao(Xk0));
        if(abs(funcaopadrao(Xk1))<e1)
        {
            return Xk1;
        }
        if(abs(Xk1-Xk)<e2)
        {
            return Xk1;
        }
        if(DBG_VIEW==1)
        {
            cout << "Xk-1: " <<Xk0;
            cout << "\n";
            cout << "Xk  : " <<Xk;
            cout << "\n";
            cout << "Xk+1: " <<Xk1;
            cout << "\n\n";
        }
        Xk0=Xk;
        Xk=Xk1;
    }
    cout << "Nao convergiu. Algoritmo METSEC.\n";
    return Xk1;
}

// Cálculo do vetor gradiente
void GradienteVet (Matriz &X, int indice, Matriz &OUT, double h)
{
    double temp;
    int n= X.dim('c') , i , j;
    Matriz H (1,n) , Xmh(1,n) , TEMPmh(1,n) , TEMP(1,n);

    for (i=0 ; i<n ; i++)
    {
        //zerar o VETOR  H
        for(j=0 ; j<n ; j++)
        {
            H.entrada(0,j, 0.0);
        }
        H.entrada(0,i,h);

        SOMA(X,H,Xmh);

        funcaopadrao(Xmh, TEMPmh);
        funcaopadrao(X, TEMP);

        temp = (TEMPmh.saida(0,indice-1) - TEMP.saida(0, indice-1))/h;

        OUT.entrada(0,i,temp);

    }


}

// Calculo da matriz Jacobiana
void Jacobiana (Matriz &X0, Matriz &MJ, int DBG_VIEW) 
{
    int i, j, n=X0.dim('c');
    Matriz Xd (1,n), MJt (MJ.dim('c'),MJ.dim('l'));
    for(i=1;i<=n;i++) //GradienteVet deriva em função da varivél i, começando de 1
    {
        GradienteVet(X0,i,Xd);

        for(j=0;j<n;j++)
        {
            MJ.entrada(i-1,j,Xd.saida(0,j));

            if(DBG_VIEW==1)
            {
                MJ.imprime();
                cout << "\n";
            }

        }

    }
}

// Regressão por mínimos quadrados
void MinimosQuad (Matriz &T, Matriz &Y, int n, Matriz &X, double &r2,int DBG_VIEW)
{
    int i, j, p = T.dim('c');
    Matriz X2 (X.dim('l'),X.dim('c'));
    if(n>=X.dim('c'))
    {
        cout << "Valor de 'n' menor que o numero de coeficientes. Funcao MINIMOSQUAD.";
        abort();
    }

   Matriz A(p,n+1), At (n+1,p), AtA(n+1,n+1), AtY (n+1,1);

    for(i=0;i<p;i++) // Colocando o valor 1 para a primeira coluna da matriz A
    {
        A.entrada(i,0,1);
    }

    for(i=0;i<p;i++) //Completando a Matriz A
    {
        for(j=1;j<=n;j++)
        {
            A.entrada(i,j,pow(T.saida(0,i),j));
        }
    }
    TRANSPOSTA(A,At);
    MULTIPLICA(At,A,AtA);
    MULTIPLICA(At,Y,AtY);
    Matriz AtYt (1,X.dim('c'));
    TRANSPOSTA(AtY,AtYt);
    Matriz AtA2(AtA.dim('l'),AtA.dim('c'));

    EliminGaussPivotParc(AtA,AtYt, AtA2, X);
    ResolTriangSup(AtA2,X,X);      // Resolvendo o sistema AtA*X = AtYt

    //r2
    double Somat1=0, Somat2=0, Somat3=0;

    if(n=1)
    {
        for (i=0;i<p;i++)
        {
            Somat1 = Somat1 + pow(Y.saida(i,0)-X.saida(0,0)-X.saida(0,1)*T.saida(0,i),2);
            Somat2 = Somat2 + Y.saida(i,0)*Y.saida(i,0);
            Somat3 = Somat3 + Y.saida(i,0);
        }

        double temp=(Somat1)/(Somat2-(1/p)*Somat3*Somat3);
    r2=1-temp;



    }

    if(DBG_VIEW==1)
    {
        A.imprime();
        cout << "\n";
        At.imprime();
        cout << "\n";
        AtA.imprime();
        cout << "\n";
        AtY.imprime();
    }

}

// Cálculo de integral por Quadratura Gaussiana
void IntegralQuadrGauss (double a, double b, int n ,double &S)
{
    int i;
    double v;

    Matriz X(1,n+1);
    Matriz W(1,n+1);
    int m=0;
    if (n==2)
    {

        double x[3] = {-0.77459666924148337704,  0,   0.77459666924148337704};
        double w[3] = {0.55555555,   0.88888888,0.55555555};
        for(i=0;i<n+1;i++)
        {
            X.entrada(0,i, x[i]);
            W.entrada(0,i, w[i]);
        }
        m=1;
    }
    if (n==4)
    {
        double x[5] = {-0.90617984593866399279,  -0.53846931010568309104,   0, 0.90617984593866399279, 0.53846931010568309104};
        double w[5] = {0.23692689,   0.47862867,   0.56888889,0.23692689,0.47862867};

        for(i=0;i<n+1;i++)
        {
            X.entrada(0,i, x[i]);
            W.entrada(0,i, w[i]);
        }
        m=1;
    }
        if (n==19)
    {
        double x[20] = {-0.0765265211334973, 0.0765265211334973, -0.2277858511416451, 0.2277858511416451, -0.3737060887154195, 0.3737060887154195, -0.5108670019508271, 0.5108670019508271, -0.6360536807265150, 0.6360536807265150, -0.7463319064601508, 0.7463319064601508, -0.8391169718222188, 0.8391169718222188, -0.9122344282513259, 0.9122344282513259 , -0.9639719272779138, 0.9639719272779138, -0.9931285991850949, 0.9931285991850949};
        double w[20] = {0.1527533871307258, 0.1527533871307258, 0.1491729864726037, 0.1491729864726037, 0.1420961093183820, 0.1420961093183820, 0.1316886384491766, 0.1316886384491766, 0.1181945319615184, 0.1181945319615184, 0.1019301198172404, 0.1019301198172404, 0.0832767415767048, 0.0832767415767048, 0.0626720483341091, 0.0626720483341091, 0.0406014298003869, 0.0406014298003869, 0.0176140071391521, 0.0176140071391521};

        for(i=0;i<n+1;i++)
        {
            X.entrada(0,i, x[i]);
            W.entrada(0,i, w[i]);
        }
        m=1;
    }

    if (m=0)
    {

        cout << "Insira um valor valido de n (2, 4 ou 19). Algoritmo INTEGRALQUADRGAUSS. Abortando...";
        abort();
    }


    double u= ((b-a)*X.saida(0,0)+a+b)/2;
    S = W.saida(0,0)*funcaopadrao(u);

    for (i=1;i<=n;i++)
    {
        u =((b-a)*X.saida(0,i)+a+b)/2;

        S = S + W.saida(0,i)*funcaopadrao(u);
    }
    S = (b-a)*S/2;


}






// Função de teste
void TesteDeFuncaoComCabecalho (void)
{
    cout << "Esta e uma funcao de teste." << endl;
}

