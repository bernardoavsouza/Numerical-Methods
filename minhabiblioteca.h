// Arquivo de cabeçalho minhabiblioteca.h. Contém os protótipos das funções
#ifndef MINHABIBLIOTECA
#define MINHABIBLIOTECA

# include <cassert>
# include <cstdlib>
# include "matriz.h"
// .
// .
// .
// .
// .
void Identidade (Matriz &);
void Copia (Matriz &, Matriz&);
void SOMA(Matriz & , Matriz &, Matriz &);
void SUBTRAI(Matriz & , Matriz &, Matriz &);
void MULTIPLICA(Matriz & , Matriz &, Matriz & );
double e2xTAYLOR (double, double);
double eaxTAYLOR (double, double, double, double);
double lnxTAYLOR (double, double, double);
double sinxTAYLOR (double, double, double);
double derivPROG (double, double = 0.01);
double derivREG (double , double = 0.01);
double derivCENT (double , double = 0.01);
double derivPARC (Matriz , int , double=0.01);
void ResolTriangSup (Matriz, Matriz, Matriz &);
void EliminGauss (Matriz, Matriz, Matriz &, Matriz &);
void EliminGaussPivotParc( Matriz , Matriz , Matriz &, Matriz &, int=0);
void Jacobi (Matriz, Matriz, Matriz, Matriz &, double=0.000001, int=1000);
void SOR (Matriz , Matriz , Matriz, Matriz&, double, double=0.000001, double=1000, int=0);
double Bissecao (double , double , int, double=0.000001, int=0);
double MetNewton1 (double ,int, int=0, double=0.000001, double=0.000001);
double MetSec (double ,double , int, int=0, double=0.000001, double=0.000001);
void GradienteVet (Matriz &, int , Matriz &, double=0.00001 );
void Jacobiana (Matriz &, Matriz &, int=0);
void TRANSPOSTA (Matriz &, Matriz &);
void MinimosQuad (Matriz &, Matriz &, int, Matriz &, double&,int=0);
void IntegralQuadrGauss (double, double, int, double&);



// .
// .
// .
// .
// .

// Função de teste.
void TesteDeFuncaoComCabecalho (void);

#endif // MINHABIBLIOTECA
