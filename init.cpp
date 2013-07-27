#include "ran1.c"
void init_fcc(cestica *c)
{
     long int i, j, k, pa=0;
     double a02=2*a0;
     
     for (k=0; k<size ;k++)
     {
         for (j=0; j<size ;j++)
         {
	     for(i=0; i<size ;i++)
             {
              c[pa].r[0]=a0*i;
              c[pa].r[1]=a0*j;
              c[pa].r[2]=a0*k;
         	  pa++;
           	  if ((i<size) && (j<size))
                  {
               		c[pa].r[0]=a0*i;
               		c[pa].r[1]=a0*j;
               		c[pa].r[2]=a0*k;
               		c[pa].r[0]+=a02;
               		c[pa].r[1]+=a02;
             		pa++;
                  }
  	          if ((k<size) && (i<size))
                  {
               		c[pa].r[0]=a0*i;
               		c[pa].r[1]=a0*j;
               		c[pa].r[2]=a0*k;
               		c[pa].r[0]+=a02;
               		c[pa].r[2]+=a02;
             		pa++;
                  }
  	          if ((k<size) && (j<size))
		  {
               		c[pa].r[0]=a0*i;
               		c[pa].r[1]=a0*j;
               		c[pa].r[2]=a0*k;
               		c[pa].r[1]+=a02;
               		c[pa].r[2]+=a02;
             		pa++;
		  }
           }
         }
     }
}

void init_v(cestica *c)
{    
     int i, j, k;
     double sum[3]={0.};
     double x1, x2, x3, fs;
     
     long idum=-time(NULL);
     
     for (i=0; i<N; i++)
     {
         x1=(double)ran1(&idum);
         x2=(double)ran1(&idum);
         x3=(double)ran1(&idum);
         
         c[i].v[0]=sqrt(-log(x1))*cos(2*pi*x2)*2;
         c[i].v[1]=sqrt(-log(x1))*cos(2*pi*x3)*2;
         c[i].v[2]=sqrt(-log(x2))*cos(2*pi*x3)*2;
         
         for (j=0; j<3; j++) sum[j]+=c[i].v[j];
     }   //postavljanje nasumicnih brzina
     
	 for (i=0; i<3; i++) sum[i]/=N;  //izraèun prosjeène driftne brzine
	 
         for (i=0;i<N;i++)
         for (j=0;j<3;j++)
             c[i].v[j]-=sum[j];  //ponistavanje drifta
             
     termalizacija(c);           //temperaturno skaliranje brzina
}
