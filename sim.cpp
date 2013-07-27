#include "konstante.cpp"
#include "pomocne.cpp"
#include "init.cpp"
#include "rdf.cpp"

using namespace std;

main(void)
{
          int i, j, k;
          double a[3], a1[3];
          
          FILE *f_en=fopen("energija.txt", "w");
          
          cestica *c = (cestica*)malloc(N*sizeof(cestica));             //definicija èestica
          double *rdf = (double*)malloc(n_rdf*sizeof(double));          //definicija rdf polja
          for(i=0; i<n_rdf; i++) rdf[i]=0.;                             //inicijalizacija rdf polja
          
          init_fcc(c);
          init_v(c);        //inicijalizacija položaja i brzina èestica
          
          provjere(c);
          
          for(k=1; k<=n; k++)                   //petlja po koracima begin
          {
                   fprintf(stderr, "%d\n", k);
                   for(i=0; i<N; i++)          //petlja po cesticama begin
                   {
                            for(j=0; j<3; j++) c[i].F[j]=0.;        //postavljanje sve tri komponente sile na trenutnu cesticu na 0
                            
                            for(j=0; j<N; j++) if(j!=i) F(&c[i], &c[j]);                 // racunanje sile izmedju dvije cestice
                            
                            
                            /******** pomicanje cestice begin ********/
                            for(j=0; j<3; j++)
                            {
                                     a[j]=c[i].F[j]/m;
                                     c[i].r[j]+=c[i].v[j]*dt+0.5*a[j]*dt*dt;
                                     
                                     if(c[i].r[j]>a0*(size+1)) c[i].r[j]-=a0*(size+1);   //
                                     else if(c[i].r[j]<0) c[i].r[j]+=a0*(size+1);        //rubni uvjeti
                            } /******** pomicanje cestice end **********/
                            
                            
                            for(j=0; j<3; j++) c[i].F[j]=0.;           //postavljanje sve tri komponente sile na cesticu i na 0
                            
                            for(j=0; j<N; j++) if(j!=i) F(&c[i], &c[j]);                //racunanje sile izmedju dvije cestice
                            
                            /******** racunanje nove brzine cestice begin ********/
                            for(j=0; j<3; j++) 
                            {
                                     a1[j]=c[i].F[j]/m;
                                     c[i].v[j]+=0.5*(a1[j]+a[j])*dt;
                            } /******** racunanje nove brzine cestice end ***********/
                            
                   }        // petlja po cesticama end
                   
                   fprintf(f_en, "%d\t%lf\t%lf\n", k, U_current(c), P_current(c));
                   
                   termalizacija(c);
                   if(k>20000 && k%50==0) RDF(c, rdf, 0);
                   
          }        //petlja po koracima end
          
          ispis(c);
          
          fprintf(stderr, "gotovo ispisivanje\n");
          
          RDF(c, rdf, 1);
          
          fprintf(stderr, "gotov rdf\n");
          
          free(c);
          free(rdf);
          
          fclose(f_en);
          //getchar();
          fprintf(stderr, "END! ");
          return 0;
}
