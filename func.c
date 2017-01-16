#include "func.h"
#include <omp.h>
#define PAD 8
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
	#pragma omp parallel
	{
	int i, j, k;
	#pragma omp for collapse (2) private (i,j,k)
	for (i=0; i<TSIZE; i++) {
   		for (j=0; j<TSIZE; j++) {
   			int r = a[i][j];
			for (k=0; k<TSIZE; k++) {
   				c[i][k]+=r*b[j][k];
			}
		}
	}
	}
}

void Func2(int d[][MSIZE], int c[][MSIZE])
{	
	int i, j, x, y;
	#pragma omp parallel for shared(d,c) private(i,j,x,y)
	for (i = 0; i < MSIZE; i+=16)
	{	for (j = 0; j < MSIZE; j+=16)
		{	for (x = i; x < MIN(i+16, MSIZE); x++)
			{
				for (y = j; y < MIN(j+16, MSIZE); y++)
					d[x][y] = c[y][x];
			}
		}
	}
}

void Func3(int z[][MSIZE], int d[][MSIZE])
{
        #pragma omp parallel
        {
        int y, x, temp;
        int near = 2;           // The weight of neighbor
        int itself = 84/2;        // The weight of center
        #pragma omp for collapse(2) private (y,x)
        for (y=0; y<MSIZE; y++) {
                for (x=0; x<MSIZE; x++) {
                        if (y==0) {
                                if (x==0) {
                                        temp = d[y][x] + 

                                                 d[y][x+1] +
                                                 d[y+1][x] + 
                                                 d[y+1][x+1] + 
                                                 d[y][x] +
                                                 d[y][x+1] +
                                                 d[y][x] + 
                                                 d[y+1][x] +
                                                itself * d[y][x];
                                } else if (x==MSIZE-1) {
                                        temp  =       d[y][x-1] +
                                                d[y][x] +
                                                d[y+1][x-1] +
                                                d[y+1][x] +
                                                d[y][x-1] +
                                                d[y][x] +
                                                d[y][x] +
                                                d[y+1][x] +
                                                itself * d[y][x];
                                } else {
                                        temp =       d[y][x-1] +
                                                d[y][x+1] +
                                                d[y+1][x-1] +
                                                d[y+1][x+1] +
                                                d[y][x-1] +
                                                d[y][x+1] +
                                                d[y][x] +
                                                d[y+1][x] +
                                                itself * d[y][x];
                                }
                        } else if (y==MSIZE-1) {
                                if (x==0) {
                                        temp = d[y-1][x] +
                                                d[y-1][x+1] +
                                                d[y][x] +
                                                d[y][x+1] +
                                                d[y][x] +
                                                d[y][x+1] +
                                                d[y-1][x] +
                                                d[y][x] +
                                                itself * d[y][x];
                                } else if (x==MSIZE-1) {
									temp =       d[y-1][x-1] +
                                                d[y-1][x] +
                                                d[y][x-1] +
                                                d[y][x] +
                                                d[y][x-1] +
                                                d[y][x] +
                                                d[y-1][x] +
                                                d[y][x] +
                                                itself * d[y][x];
                                } else {
                                        temp  =       d[y-1][x-1] +
                                                d[y-1][x+1] +
                                                d[y][x-1] +
                                                d[y][x+1] +
                                                d[y][x-1] +
                                                d[y][x+1] +
                                                d[y-1][x] +
                                                d[y][x] +
                                                itself * d[y][x];
                                }
                        } else {
                                if (x==0) {
                                        temp =      d[y-1][x] +
                                                d[y-1][x+1] +
                                                d[y+1][x] +
                                                d[y+1][x+1] +
                                                d[y][x] +
                                                d[y][x+1] +
                                                d[y-1][x] +
                                                d[y+1][x] +
                                                itself * d[y][x];
                                } else if (x==MSIZE-1) {
                                        temp =      d[y-1][x-1] +
                                                d[y-1][x] +
                                                d[y+1][x-1] +
                                                d[y+1][x] +
                                                d[y][x-1] +
                                                d[y][x] +
                                                d[y-1][x] +
                                                d[y+1][x] +
                                                itself * d[y][x];
                                } else {
                                        temp =       d[y-1][x-1] +
                                                d[y-1][x+1] +
                                                d[y+1][x-1] +
                                                d[y+1][x+1] +
                                                d[y][x-1] +
                                                d[y][x+1] +
                                                d[y-1][x] +
                                                d[y+1][x] +
                                                itself * d[y][x];
                                }
                        }
                        z[y][x]= temp/50;
                }
        }
        }
}



/*
void Func3(int z[][MSIZE], int d[][MSIZE])
{
	#pragma omp parallel
	{
	int y, x;
	int near = 2;  		// The weight of neighbor
	int itself = 84/2; 	// The weight of center
	#pragma omp for collapse(2) private (y,x)
	for (y=0; y<MSIZE; y++) {
		for (x=0; x<MSIZE; x++) {
			if (y==0) {
				if (x==0) {
					z[y][x] =  d[y][x] +

                                                 d[y][x+1] +
                                                 d[y+1][x] +
                                                 d[y+1][x+1] +
                                                 d[y][x] +
                                                 d[y][x+1] + 
                                                 d[y][x] + 
                                                 d[y+1][x] +
                                                itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = d[y][x-1] +
						d[y][x] +
						d[y+1][x-1] +
						d[y+1][x] +
						d[y][x-1] +
						d[y][x] +
						d[y][x] +
						d[y+1][x] +
						itself * d[y][x];
				} else {
					z[y][x] = d[y][x-1] +
						d[y][x+1] +
						d[y+1][x-1] +
						d[y+1][x+1] +
						d[y][x-1] +
						d[y][x+1] +
						d[y][x] +
						d[y+1][x] +
						itself * d[y][x];
				}
			} else if (y==MSIZE-1) {
				if (x==0) {
					z[y][x] = d[y-1][x] +
						d[y-1][x+1] +
						d[y][x] +
						d[y][x+1] +
						d[y][x] +
						d[y][x+1] +
						d[y-1][x] +
						d[y][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = d[y-1][x-1] +
						d[y-1][x] +
						d[y][x-1] +
						d[y][x] +
						d[y][x-1] +
						d[y][x] +
						d[y-1][x] +
						d[y][x] +
						itself * d[y][x];
				} else {
					z[y][x] = d[y-1][x-1] +
						d[y-1][x+1] +
						d[y][x-1] +
						d[y][x+1] +
						d[y][x-1] +
						d[y][x+1] +
						d[y-1][x] +
						d[y][x] +
						itself * d[y][x];
				}
			} else {
				if (x==0) {
					z[y][x] = d[y-1][x] +
						d[y-1][x+1] +
						d[y+1][x] +
						d[y+1][x+1] +
						d[y][x] +
						d[y][x+1] +
						d[y-1][x] +
						d[y+1][x] +
						d[y][x] + itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = d[y-1][x-1] +
						 d[y-1][x] +
						d[y+1][x-1] +
						d[y+1][x] +
						d[y][x-1] +
						d[y][x] +
						d[y-1][x] +
						d[y+1][x] +
						itself * d[y][x];
				} else {
					z[y][x] = d[y-1][x-1] +
						d[y-1][x+1] +
						d[y+1][x-1] +
						d[y+1][x+1] +
						d[y][x-1] +
						d[y][x+1] +
						d[y-1][x] +
						d[y+1][x] + itself * d[y][x];
				}
			}
			z[y][x]/=50;
		}
	}
	}
}*/
						
void Func4(int b[], int a[])
{
	int chuck_size=MSIZE;	 
	int array_size=VSIZE/chuck_size;
	int chuck[chuck_size];
   	int i, j;
	#pragma omp parallel for shared (a,b, chuck)  private (j,i)
	for(j=0; j<chuck_size; j++) {
		int jarr = j*array_size;
		b[jarr]=a[jarr];
		for (i=1; i<array_size; i++) {
			b[jarr+i]=b[jarr+i-1]+a[jarr+i];
		}
		chuck[j]=b[jarr + array_size-1];
	}
	
//	#pragma omp for private (j)	
	for(j=1; j<chuck_size; j++) {
		chuck[j]=chuck[j-1]+chuck[j];
	}
	#pragma omp parallel for shared (b) private (j,i)
	for(j=1; j<chuck_size; j++) {
		int jarr = j*array_size;
		int chuckj =chuck[j-1]/(j+1);
		for (i=0; i<array_size; i++) {
			b[jarr+i]+=chuckj;
		}
	}
}

void Func5(int b[], int a[])
{
//	#pragma omp parallel
	{
	int i=0, j,  stride, stride2, step;
    int temp;
	long log_N=log2(VSIZE);
	
	#pragma omp parallel for
	for(j=0; j<VSIZE; j+=2) {
		b[j]=a[j];
		b[j+1] = a[j] + a[j+1];
	}
	
//	#pragma omp parallel for	
	for(i=4; i<VSIZE; i*=2) {
		for(j=0; j<VSIZE; j+=i) {
				b[j+i-1] = b[j+i/2-1] + b[j+i-1];
		}
	}
	
	b[VSIZE-1]=0;
	for(i=(log_N-1); i>=0; i--) {
		//i = 1;
		stride2=(2<<i)-1;//1
		stride=(1<<i)-1; //0
		step=stride2+1; //2
		for(j=0; j<VSIZE; j+=step) {
                temp=b[j+stride];
			b[j+stride] = b[j+stride2];
			b[j+stride2] = temp+b[j+stride2];
		}
	}
	}
}
