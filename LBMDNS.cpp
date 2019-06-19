#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Parameter for setup of the computational zone*/
const int Nx = 200;
const int Ny = 800;
const double g = 3.5e-7;	//gravity
int t;	//time step

/*Parameter for Bilinear Interpolation*/
double ESMatrix[121][121];	//solid concentration matrix, reading from SolidConcentration.txt
double ESMatrixInterval = 0.01;	//sampling interval for data in ESMatrix

/*Parameter for LBM*/
double w[9];	//weight
int e_x[9], e_y[9];	//discrete velocities
double f[Nx + 1][Ny + 1][9];	//distribution function in LBM
double ftemp[Nx + 1][Ny + 1][9];	//temp distribution function in LBM

/*Parameter for the gas phase*/
double rhog[Nx + 1][Ny + 1];	//gas density, which is corresponding to pressure in LBM
double ugx[Nx + 1][Ny + 1];	//gas velocity in x direction
double ugy[Nx + 1][Ny + 1];	//gas velocity in y direction
const double tao = 0.53;//1.796;	//gas viscosity

/*Parameter for the particle phase*/
const int PIS = 300;	//particle number in the system
const double Radius = 5.4;	//particle diameter
const double ms = 91608.84;	//particle mass, corrsponding to rhos=1000
double POS[PIS][2];	//particle position
double VEL[PIS][2];	//particle velocity

/*Parameter for the interaction between the gas and particle phase*/
int IB[Nx + 1][Ny + 1][2];	//[][][0]--->interaction check, check if interaction exists(=1) or only gas phase(=0)
							//[][][1]--->interaction index, store the particle number(interaction exists) or equal to -1(only gas phase)
double InterF[Nx + 1][Ny + 1][2];	//"interaction force" on one gas cell
double Drag[PIS][2];	//drag between the gas and particle phase

void EsDataRead()	//read solid concentration into ESMatrix
{
	FILE *fp;
	fp = fopen("SolidConcentration.txt", "r");
	int i, j;
	for (i = 0; i <= 120; i++)
	{
		for (j = 0; j <= 120; j++)
		{
			double temp1, temp2;
			fscanf(fp, "%lf,%lf,%lf", &temp1, &temp2, &ESMatrix[i][j]);
		}
	}
	fclose(fp);
}

double EsInterp(double x,double y)	//calculate solid concentration using bilinear interpolation
{
	double es;
	int I_x, I_y;
	I_x = int(x / ESMatrixInterval);
	I_y = int(y / ESMatrixInterval);

	int ii1, ii2, jj1, jj2;
	ii1 = I_x;
	ii2 = I_x + 1;
	jj1 = I_y;
	jj2 = I_y + 1;

	double temp1, temp2, temp3, temp4;
	temp1 = (ii2*0.01 - x)*(jj2*0.01 - y) / 0.01 / 0.01*ESMatrix[ii1][jj1];
	temp2 = (x - ii1 * 0.01)*(jj2*0.01 - y) / 0.01 / 0.01*ESMatrix[ii2][jj1];
	temp3 = (ii2*0.01 - x)*(y - jj1 * 0.01) / 0.01 / 0.01*ESMatrix[ii1][jj2];
	temp4 = (x - ii1 * 0.01)*(y - jj1 * 0.01) / 0.01 / 0.01*ESMatrix[ii2][jj2];
	es = temp1 + temp2 + temp3 + temp4;
	return es;
}

void SetLBMParameter()	//D2Q9
{
	w[0] = 4.0 / 9.0;
	w[1] = w[2] = w[3] = w[4] = 1.0 / 9.0;
	w[5] = w[6] = w[7] = w[8] = 1.0 / 36.0;

	e_x[0] =  0, e_y[0] =  0;
	e_x[1] =  1, e_y[1] =  0;
	e_x[2] =  0, e_y[2] =  1;
	e_x[3] = -1, e_y[3] =  0;
	e_x[4] =  0, e_y[4] = -1;
	e_x[5] =  1, e_y[5] =  1;
	e_x[6] = -1, e_y[6] =  1;
	e_x[7] = -1, e_y[7] = -1;
	e_x[8] =  1, e_y[8] = -1;
}

double f_equ(double ux, double uy, double rho, int b)  //calculate the value of the equilibrium distribution function
{
	double feq;
	double e_u = e_x[b] * ux + e_y[b] * uy;
	double usq = ux * ux + uy * uy;
	feq = w[b] * rho*(1 + 3 * e_u + 4.5*e_u*e_u - 1.5*usq);
	return feq;
}

void init()	//initialization
{
	int i, j;
	int b;
	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
		{
			rhog[i][j] = 1.0;
			ugx[i][j] = 0.0;
			ugy[i][j] = 0.1;//0.0;

			for (b = 0; b < 9; b++)
			{
				f[i][j][b] = f_equ(ugx[i][j], ugy[i][j], rhog[i][j], b);
			}
			IB[i][j][0] = 0;
			InterF[i][j][0] = 0.0;
			InterF[i][j][1] = 0.0;
			IB[i][j][1] = -1;
		}
	}

	int k;
	for (k = 0; k < PIS; k++)
	{
		POS[k][0] = (k % 20)*11.0 + 5.5 + ((k / 20) % 2)*5.0;
		POS[k][1] = (k / 20)*9.8 + 5.5;
		VEL[k][0] = 0.0;
		VEL[k][1] = 0.0;
		Drag[k][0] = 0.0;
		Drag[k][1] = 0.0;
	}
}

int main()
{
	EsDataRead();	//read solid concentration into ESMatrix
	SetLBMParameter();
	init();

	/*test code*/
	double x, y, es;
	x = 0.490;
	y = 0.840;
	es = EsInterp(x, y);
	printf("%lf\t%lf\t%lf\n", x, y, es);
	printf("w[1]=%lf\te_x[6]=%d\n", w[1], e_x[6]);
	printf("f[100][400][5]=%lf\n", f[100][400][5]);
	getchar();
	/*test code*/

	return 1;
}