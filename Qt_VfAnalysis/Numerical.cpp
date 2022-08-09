#include <malloc.h>
#include "Numerical.h"

/****************************************************
Maybe better for 2X2 and 3X3 matrix inversion
****************************************************/
double *MatrixOpp(double A[],int m,int n) //inverse
{
    int i,j,x,y,k;
    double *SP=nullptr,*AB=nullptr,*B=nullptr,XX,*C;
    SP=(double *)malloc(m*n*sizeof(double));
    AB=(double *)malloc(m*n*sizeof(double));
    B=(double *)malloc(m*n*sizeof(double));

    XX=Surplus(A,m,n);
    XX=1/XX;

    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
        {
            for(k=0;k<m*n;k++)
                B[k]=A[k];
            {
                for(x=0;x<n;x++)
                    B[i*n+x]=0;
                for(y=0;y<m;y++)
                    B[m*y+j]=0;
                B[i*n+j]=1;
                SP[i*n+j]=Surplus(B,m,n);
                AB[i*n+j]=XX*SP[i*n+j];
            }
        }

    C=MatrixInver(AB,m,n);

    return C;
}


double * MatrixInver(double A[],int m,int n) //zhuanzhi
{
    int i,j;
    double *B=nullptr;
    B=(double *)malloc(m*n*sizeof(double));

    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            B[i*m+j]=A[j*n+i];

    return B;
}


double Surplus(double A[],int m,int n) //hanglieshi
{

    int i,j,k,p,r;
    double XX,temp=1,temp1=1,s=0,s1=0;

    if(n==2)
    {
        for(i=0;i<m;i++)
            for(j=0;j<n;j++)
                if((i+j)%2) temp1*=A[i*n+j];
                else temp*=A[i*n+j];
        XX=temp-temp1;
    }
    else{
        for(k=0;k<n;k++)
        {
            for(i=0,j=k;i<m,j<n;i++,j++)
                temp*=A[i*n+j];
            if(m-i)
            {
                for(p=m-i,r=m-1;p>0;p--,r--)
                    temp*=A[r*n+p-1];
            }
            s+=temp;
            temp=1;
        }

        for(k=n-1;k>=0;k--)
        {
            for(i=0,j=k;i<m,j>=0;i++,j--)
                temp1*=A[i*n+j];
            if(m-i)
            {
                for(p=m-1,r=i;r<m;p--,r++)
                    temp1*=A[r*n+p];
            }
            s1+=temp1;
            temp1=1;
        }

        XX=s-s1;
    }
    return XX;
}
