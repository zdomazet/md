inline void termalizacija(cestica *c)
{
     int i, j;
     double fs, Ek=0.;
     
     for(i=0; i<N; i++)
              for(j=0; j<3; j++)
                       Ek+=c[i].v[j]*c[i].v[j];     //kineticka energija
     Ek*=0.5*m;
     
     fs=sqrt(U/Ek);    //faktor skaliranja
     
     for(i=0; i<N; i++)
     for(j=0; j<3; j++) c[i].v[j]*=fs;     //temperaturno skaliranje
}

void provjere(cestica *c)
{
     int i, j;
     double temp, sum[3]={0.}, Ek=0.;
     
     printf("N=%d\n\n", N);
     
     for(i=0; i<N; i++)
     for (j=0; j<3; j++) sum[j]+=c[i].v[j];
     printf("<vx>=%lf\n<vy>=%lf\n<vz>=%lf\n\n", sum[0], sum[1], sum[2]); //provjera ukupne brzine
     
     for(i=0; i<N; i++)
     for(j=0; j<3; j++)
     Ek+=0.5*m*c[i].v[j]*c[i].v[j];
     
     temp=m*Ek/3/N/kb;
     
     printf("temperatura: %lfK (%lfK zadano)\n\n", temp, T); //provjera temperature
     
     printf("n_rdf=%d\n",n_rdf);
     
     printf("a0=%lf\nsize=%d\na0*size=%lf\na0*(size+1)=%lf\n\n", a0, size, a0*(size), a0*(size+1));
     
     printf("A=%.32lf\tB=%.32lf\n", A, B);
     
     printf("U_trunc = %lf\n", U_trunc);
}

void ispis(cestica *c)
{
     FILE *fid;
     int j;
     if((fid=fopen("brzine.txt","w"))!=NULL)
     {
          for (j=0; j<N; j++)
          fprintf(fid,"%+.8lf\t%+.8lf\t%+.8lf\t%+.8lf\n", c[j].v[0], c[j].v[1], c[j].v[2], sqrt(c[j].v[0]*c[j].v[0]+c[j].v[1]*c[j].v[1]+c[j].v[2]*c[j].v[2]));
          fclose(fid);
     }
     
     if((fid=fopen("koordinate.txt","w"))!=NULL)
     {
          for (j=0; j<N; j++)
          fprintf(fid,"%+.8lf\t%+.8lf\t%+.8lf\n", c[j].r[0],c[j].r[1],c[j].r[2]);
          fclose(fid);
     }
}

inline void F(cestica *c1, cestica *c2)
{
            double dx, dy, dz, d, Fu;
            
            dx=c1->r[0]-c2->r[0];         //
            dy=c1->r[1]-c2->r[1];         //
            dz=c1->r[2]-c2->r[2];         // racunanje medjusobne udaljenosti dvije cestice po koordinatama
            
            if(dx>a0*(size+1)/2)            dx-=a0*(size+1);   //
            else if(dx<-a0*(size+1)/2)      dx+=a0*(size+1);   //
            if(dy>a0*(size+1)/2)            dy-=a0*(size+1);   //
            else if(dy<-a0*(size+1)/2)      dy+=a0*(size+1);   //
            if(dz>a0*(size+1)/2)            dz-=a0*(size+1);   //
            else if(dz<-a0*(size+1)/2)      dz+=a0*(size+1);   // rubni uvjeti
            
            d=sqrt(dx*dx+dy*dy+dz*dz);
            
            if(d<cutoff)
            {
                      Fu=-12*A/pow(d,13)+6*B/pow(d,7);          // ukupna sila
                      c1->F[0]-=dx*Fu/d;                         //
                      c1->F[1]-=dy*Fu/d;                         //
                      c1->F[2]-=dz*Fu/d;                         // komponente sile cestice c2 na cesticu c1
            }
}

inline double P_current(cestica *c)
{
       int i, j;
       double P=0.;
       
       for(i=0; i<N; i++)
       for(j=0; j<3; j++)
       P+=c[i].r[j]*c[i].F[j];
       
       return N*kb*T/V - P/V/3/N;
       
}

inline double U_current(cestica *c)
{
       int i, j;
       double dx, dy, dz, d, Utot=0.;
       
       for(i=0; i<N-1; i++)
       for(j=i+1; j<N; j++)
       {
            dx=c[i].r[0]-c[j].r[0];         //
            dy=c[i].r[1]-c[j].r[1];         //
            dz=c[i].r[2]-c[j].r[2];         // racunanje medjusobne udaljenosti dvije cestice po koordinatama
            
            if(dx>a0*(size+1)/2)            dx-=a0*(size+1);   //
            else if(dx<-a0*(size+1)/2)      dx+=a0*(size+1);   //
            if(dy>a0*(size+1)/2)            dy-=a0*(size+1);   //
            else if(dy<-a0*(size+1)/2)      dy+=a0*(size+1);   //
            if(dz>a0*(size+1)/2)            dz-=a0*(size+1);   //
            else if(dz<-a0*(size+1)/2)      dz+=a0*(size+1);   // rubni uvjeti
            
            d=sqrt(dx*dx+dy*dy+dz*dz);
            
            if(d<cutoff) Utot += A/pow(d,12) - B/pow(d,6) - U_trunc;
       }
       
       return Utot;
}
