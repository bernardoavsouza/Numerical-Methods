# Numerical-Methods
These are some algorithms for troubleshooting mathematics.
The whole code is commented in portuguese and it was done in 2017 on CodeBlocks IDE.


This repository was coded based on a class called Matriz (that is Matrix) developed by Dr. Daniel Gomes Ribeiro.

The repository is separated by 6 main files: 
        BaseDeProgramacao.cbp   --> The project file that contains the following files
        minhabiblioteca.h       --> The header of the library created
        minhasfuncoes.h         --> The header of the used functions
        minhabiblioteca.cpp     --> Algorithms of the developed library
        minhasfuncoes.cpp       --> Definition of the used functions
        main.cpp                --> This is the main of the repository
        
The algorithms done here were:

        Identidade              --> It turns the input matrix into a indentity matrix
        Copia                   --> It creates a Copy of a matrix
        SOMA                    --> It sums two different matrices
        SUBTRAI                 --> It subtracts two matrices
        MULTIPLICA              --> It multiplies two matrices
        e2xTAYLOR               --> It generates Taylor series approximation of e^(-2x)
        eaxTAYLOR               --> It generates Taylor series approximation of e^(-ax) 
        lnxTAYLOR               --> It generates Taylor series approximation of ln(x)
        sinxTAYLOR              --> It generates Taylor series approximation of sin(x)
        derivPROG               --> It calculates the progressive dfferentiation
        derivREG                --> It calculates the regressive dfferentiation
        derivCENT               --> It calculates the centered dfferentiation
        derivPARC               --> It calculates the partial differentiation
        ResolTriangSup          --> It solves a upper triangular matrix system of equations
        EliminGauss             --> It does the Gaussian elimination to turn a system of equations into a upper triangular matrix
        EliminGaussPivotParc    --> It does the Gaussian elimination with pivoting to turn a system of equations into a upper triangular matrix
        Jacobi                  --> It troubleshoots a equation using the Jacobi method
        SOR                     --> It troubleshoots a equation using the SOR method
        Bissecao                --> It troubleshoots a equation using the bisection method
        MetNewton1              --> It troubleshoots a equation using the one variable Newton method
        GradienteVet            --> It calculates the gradient vector of a function
        Jacobiana               --> It calculates the jacobian matrix of a function
        TRANSPOSTA              --> It calculates the transpose of a given matrix
        MinimosQuad             --> It does a regression based on least squares method
        IntegralQuadrGauss      --> It calculates the integral of a function based on Gaussian quadrature method
        
