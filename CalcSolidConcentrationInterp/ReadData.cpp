#include <STDIO.H>
#include <STDLIB.H>

const double r=1.0;
const double a=0.1;

float M[121][121];

void EsDataRead()
{
	FILE *fp;
	fp=fopen("dat.txt","r");
	int i,j;
	for (i=0;i<=120;i++)
	{
		for (j=0;j<=120;j++)
		{
			float temp1,temp2;
			fscanf(fp,"%f,%f,%f",&temp1,&temp2,&M[i][j]);
			//printf("%f\t%f\t%f\n",temp1,temp2,M[i][j]);
		}
	}
	fclose(fp);
}

int main()
{
	EsDataRead();
	double x,y,es;
	x=0.800;
	y=0.665;

	int I_x,I_y;
	I_x=x/0.01;
	I_y=y/0.01;

	int ii1,ii2,jj1,jj2;
	ii1=I_x;
	ii2=I_x+1;
	jj1=I_y;
	jj2=I_y+1;

	double temp1,temp2,temp3,temp4;
	temp1=(ii2*0.01-x)*(jj2*0.01-y)/0.01/0.01*M[ii1][jj1];
	temp2=(x-ii1*0.01)*(jj2*0.01-y)/0.01/0.01*M[ii2][jj1];
	temp3=(ii2*0.01-x)*(y-jj1*0.01)/0.01/0.01*M[ii1][jj2];
	temp4=(x-ii1*0.01)*(y-jj1*0.01)/0.01/0.01*M[ii2][jj2];
	es=temp1+temp2+temp3+temp4;

	printf("%f\t%f\t%f\n",x,y,es);
	return 1;
}