#include "xonorm.h"
#include <math.h>

/*  Written in 2014 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

Modified 2015 to produce normally distributed random numbers for use in
coupling simulations with noisy observations.  */

static uint64_t x=12455480610190296555ULL;
static uint64_t s[16] = {
    9899521803441930453ULL, 4037224574094532358ULL,
    13192621976843192702ULL, 10923948763034344365ULL,
    3888343346591572214ULL, 7485323648124204777ULL,
    8981086447562432722ULL, 1561717205308781184ULL,
    2780386954292984643ULL, 3863312353808779985ULL,
    12109031751610640964ULL, 11534438648781303991ULL,
    9340076221226405573ULL, 219462499049377941ULL,
    16106326203448552520ULL, 7245732521397473901ULL };
static int p,flag;

uint64_t xorand64(){
    x^=x>>12; x^=x<<25; x^=x>>27;
    return x*2685821657736338717ULL;
}

uint64_t xorand1024(){
    uint64_t s0=s[p];
    uint64_t s1=s[p=(p+1)&15];
    s1^=s1<<31; s1^=s1>>11; s0^=s0>>30;
    return (s[p]=s0^s1)*1181783497276652981ULL;
}

void xoseed(uint64_t seed){
    x=seed;
    for(int i=0;i<64;i++) xorand64();
    for(int i=0;i<16;i++) s[i]=xorand64();
    p=flag=0;
}

double xonorm(){
    static double x2;
    if(flag) {
        flag=0;
        return x2;
    }
    flag=1;
    double u1=(xorand1024()+1.0)/18446744073709551616.0;
    double u2=(xorand1024()+1.0)/18446744073709551616.0;
    double r=sqrt(-2*log(u1)), v=2*M_PI*u2;
    x2=r*sin(v);
    return r*cos(v);
}
