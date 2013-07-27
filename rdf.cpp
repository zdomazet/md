inline double udaljenost(double r1[3], double r2[3])
{
       double dx, dy, dz, d;
       dx=fabs(r1[0]-r2[0]);    //
       dy=fabs(r1[1]-r2[1]);    //
       dz=fabs(r1[2]-r2[2]);    //raèunanje udaljenosti po svakoj koordinati izmeðu dvije èestice
       
       if(dx>(double)(size+1)*a0/2.) dx=fabs(dx-(double)(size+1)*a0);        //
       if(dy>(double)(size+1)*a0/2.) dy=fabs(dy-(double)(size+1)*a0);        //
       if(dz>(double)(size+1)*a0/2.) dz=fabs(dz-(double)(size+1)*a0);        //rubni uvjeti za udaljenosti u razmatranju rdf-a
       
       d = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
       
       return d;
}

void FFT(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    FFT(even);
    FFT(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

inline void RDF(cestica *c, double *rdf, int final)
{
     int i, j, nfft;
     static int nr=0;
     
     for(i=0; i<N-1; i++)
     for(j=i+1; j<N; j++)
                rdf[(int)ceil(udaljenost(c[i].r, c[j].r)/rdf_dr)]+=1.;    //skupljanje èestica u binove
     nr++;
     
     if(!final) return;         //ako nije zadnji poziv funkcije, vraæa se kontrola glavnom programu, a ako je zadnji poziv, funkcija raèuna statistiku i strukturni faktor.
     
     double ro=N/V;
     for(i=1; i<n_rdf; i++) rdf[i]/=(4*M_PI*i*i*ro*pow(rdf_dr,3)*nr);   //normiranje rdf-a
     
     FILE *fout=fopen("RDF.txt", "w");
     for(i=0; i<n_rdf; i++) fprintf(fout, "%lf\t%lf\n", i*rdf_dr/sigma, rdf[i]);     //ispis rdf-a
     fclose(fout);
     
     nfft = (int)pow(2, ceil(log(n_rdf)/log(2)));       //broj clanova polja za racunanje sf-a. najmanja visa potencija od 2
     
     Complex *fft=(Complex*)malloc(nfft*sizeof(Complex));      //definicija polja
     for(i=0; i<nfft; i++) fft[i]=0.;                          //popunjavanje nulama
     for(i=0; i<n_rdf; i++) fft[i]=rdf[i];                     //prebacivanje rdf podataka u novo polje
     
     CArray data(fft, nfft);                                   //oblikovanje podataka za lakšu obradu
     
     FFT(data);       //fft
     
     
     std::ofstream sfout;
     sfout.open("sf.txt");
     for (int i = 0; i < nfft; i++)
     {
        sfout << i << "\t" << abs(data[i]) << std::endl;      //ispis fft-a
     }
     sfout.close();
}
