/*  
    This program uses the original form of the Lorenz equations
    as stated in 

    Lorenz, Deterministic Nonperiodic Flow, Journal of the Atmospheric
    Sciences, Volume 20, 1963.

    X' = -sigma X + sigma Y
    Y' = -XZ + rX - Y
    Z' = XY - bZ

    The change of variables zeta = Z - r - sigma yields

    X' = -sigma X + sigma Y
    Y' = -sigma X - Y - X zeta
    zeta' = - b zeta + XY - b(r+sigma)

    We know theoretically that if (X,Y,zeta) is on the attractor, then

    X^2+Y^2+zeta^2 <= K = b^2(r+sigma)^2/4/(b-1) = 1540.27
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "xonorm.h"

#define BETA (8.0/3)
#define SIGMA 10.0
#define RHO 28.0

#ifndef M0
#define M0 512.0
#endif
#ifndef M3
#define M3 512.0
#endif
#ifndef N
#define N 2048
#endif
#ifndef EPSILON
#define EPSILON 0.0
#endif
#ifndef MU
#define MU 30
#endif
#ifndef M
#define M 50
#endif
#ifndef PMOD
#define PMOD (N/64)
#endif

double beta=BETA, sigma=SIGMA, rho=RHO, mu=MU, epsilon=EPSILON;
double x[6]={ -5.5751, -3.24345, 26.9272, 0,0,0 }, f[6];
double xmax[6] = { FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN };
double xmin[6] = { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX };
double vmax = 0;
double myeta = - SIGMA - RHO;
double Umax = 0;
double MM0=M0, MM3=M3;

// Cache of Brownian motion paths
#define L 5
static struct {
    int j,q;
    double b[N+1];
} p[L];
static double e[M+1];
static uint64_t s[M];

// Return normal N(mu,sigma) distributed random variable
static double normal(double mu,double sigma2){
    return sqrt(sigma2)*xonorm()+mu;
}

// Levy-Ciesielski Algorithm on the interval [l,l+1]
static void fill(int l){
    for(int delta=N;delta>1;delta/=2){
        int K=N/delta;
        for(int k=0;k<K;k++){
            double sigma=0.25*delta/N;
            p[l].b[delta/2+k*delta]=0.5*(p[l].b[k*delta]+p[l].b[(k+1)*delta])
                +normal(0,sigma);
        }
    }
}

// Initialize the endpoints of the intervals [l,l+1] on which
// we apply the Levy-Ciesieslki algorithm ahead of time
static void doinit(uint64_t seed){
    for(int l=0;l<L;l++) p[l].j=-1;
    xoseed(seed);
    e[0]=0;
    for(int i=0;i<M;i++){
        s[i]=xorand1024();
        e[i+1]=e[i]+xonorm();
    }
}

// Return standard Brownian motion W(t_i) where t_i=i/N
double getb(int i){
    static int q=0;
    int j=i/N;
    for(int l=0;l<L;l++) if(p[l].j==j){
        p[l].q=++q;
        return p[l].b[i%N];
    }
    int lm=0,qm=p[0].q;
    for(int l=1;l<L;l++) if(p[l].q<qm) {
        lm=l;
        qm=p[l].q;
    }
    p[lm].j=j;
    xoseed(s[j]);
    p[lm].b[0]=e[j];
    p[lm].b[N]=e[j+1];
    fill(lm);
    return p[lm].b[i%N];
}

// Forcing function for the Lorenz system
void force(double f[6],double x[6]){
    f[0]=sigma*(x[1]-x[0]);
    f[1]=rho*x[0]-x[1]-x[0]*x[2];
    f[2]=-beta*x[2]+x[0]*x[1];
    f[3]=sigma*(x[4]-x[3]);
    f[4]=rho*x[3]-x[4]-x[3]*x[5];
    f[5]=-beta*x[5]+x[3]*x[4];
}

// Norm of ||U-u|| where U=(x[0],x[1],x[2]) and u=(x[3],x[4],x[5])
double norm(double x[6]){
    double r=0;
    for(int i=0;i<3;i++){
        double z=x[i]-x[i+3];
        r+=z*z;
    }
    return sqrt(r);
}

// Compute time-averages observations and feedback
double avg(double mM,int m,double *xd,int n){
    double r=0;
    int k;
    for(k=n;k>n+1-m;k--) r+=xd[k%m];
    r+=xd[k%m]*(mM+1-m);
    return r/mM;
}

int main(int argc,char *argv[]){
    if(MM0<1) MM0=1;
    if(MM3<1) MM3=1;
    int m0=ceil(MM0)+0.5;
    int m3=ceil(MM3)+0.5;
    double x0a[m0],x3a[m3];
    int m=m0>m3?m0:m3;
    uint64_t seed=36819848450518068ULL;
    if(argc>1) seed=strtoull(argv[1],0,0);
    doinit(seed);
    double dt=1.0/N;
    int K=M*N;
    printf("#tn wn X Y Z x y z |U-u| avg(X) avg(x) noise\n");
    printf("# m0=%d\n",m0);
    printf("# m3=%d\n",m3);
    printf("# dt=%g\n",dt);
    printf("# delta_2=%g\n",MM0*dt);
    printf("# delta=%g\n",MM3*dt);
    printf("# epsilon=%g\n",epsilon);
    printf("# seed=%llu\n",seed);

    // Euler's method to solve for U and u
    for(int n=0;;n++){
        x0a[n%m0]=x[0];
        x3a[n%m3]=x[3];
        double tn=n*dt;
        for(int i=0;i<6;i++){
            if(xmin[i]>=x[i]) xmin[i]=x[i];
            if(xmax[i]<=x[i]) xmax[i]=x[i];
        }
        {
            double U=x[0]*x[0]+x[1]*x[1]+(x[2]+myeta)*(x[2]+myeta);
            if(Umax<U) Umax=U;
        }
        if(!(n%PMOD)) {
            if(n+1>=m){
                printf("%g %g %g %g %g %g %g %g %g %g %g %g\n",
                tn,epsilon*getb(n),x[0],x[1],x[2],x[3],x[4],x[5],norm(x),
                    avg(MM0,m0,x0a,n),avg(MM3,m3,x3a,n),
                    epsilon*(getb(n+1)-getb(n+1-m0))/m0/dt);
            } else {
                printf("%g %g %g %g %g %g %g %g %g 0 0 0\n",
                tn,epsilon*getb(n),x[0],x[1],x[2],x[3],x[4],x[5],norm(x));
            }
        }
        if(n>=K) break;
        force(f,x);
        if(n+1>=m) x[3]+=mu*(avg(MM0,m0,x0a,n)-avg(MM3,m3,x3a,n))*dt;
        for(int i=0;i<6;i++) x[i]+=f[i]*dt;
        if(n+1>=m) x[3]+=mu*epsilon*(getb(n+1)-getb(n+1-m0))/m0;
    }
    printf("# Max %g %g %g %g %g %g\n",
        xmax[0], xmax[1], xmax[2], xmax[3], xmax[4], xmax[5]);
    printf("# Min %g %g %g %g %g %g\n",
        xmin[0], xmin[1], xmin[2], xmin[3], xmin[4], xmin[5]);
    printf("# Umax %.8g\n", Umax);
    return 0;
}
