#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Parameter for Bilinear Interpolation*/
double ESMatrix[121][121];	//solid concentration matrix, reading from SolidConcentration.txt
double ESMatrixInterval = 0.01;	//sampling interval for data in ESMatrix

/*Parameter for LBM*/
double w[9];	//weight
int e_x[9], e_y[9];	//discrete velocities

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

int main()
{
	EsDataRead();	//read solid concentration into ESMatrix
	SetLBMParameter();

	/*test code*/
	double x, y, es;
	x = 0.490;
	y = 0.840;
	es = EsInterp(x, y);
	printf("%lf\t%lf\t%lf\n", x, y, es);
	double temp = e_x[6];
	printf("w[1]=%lf\te_x[6]=%lf\n", w[1],temp);
	getchar();
	/*test code*/

	return 1;
}