/* Independent brute force (same definition as brute.py), square table lookup. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { long q,a,b,c,d; } row;
static int cmp(const void*x,const void*y){const long*p=x,*r=y;for(int i=0;i<5;i++){if(p[i]<r[i])return -1;if(p[i]>r[i])return 1;}return 0;}
int main(int argc,char**argv){
  long N=atol(argv[1]),M=atol(argv[2]); long L=2*M*M+1;
  char*S=calloc(L,1); for(long r=0;r*r<L;r++) S[r*r]=1;
  size_t cap=1<<20,n=0; row*R=malloc(cap*sizeof(row)); long four=0;
  for(long Q=1;Q<=N;Q++) for(long A=Q-M;A<=M;A++){ long B=Q-A; if(labs(A)>M||labs(B)>M) continue; long a=labs(A),b=labs(B);
    for(long C=Q-M;C<=M;C++){ long D=Q-C; if(labs(C)>M||labs(D)>M) continue; long c=labs(C),d=labs(D);
      if(!S[a*a+c*c]||!S[b*b+c*c]) continue; int s3=S[a*a+d*d],s4=S[b*b+d*d];
      if(s3&&s4){four++; printf("FOUR Q=%ld A=%ld C=%ld\n",Q,A,C);}
      if(c>=1&&a>=1&&b>=1&&a!=b&&(s3||s4)){ if(n==cap){cap*=2;R=realloc(R,cap*sizeof(row));}
        R[n].q=Q;R[n].a=a<b?a:b;R[n].b=a<b?b:a;R[n].c=c;R[n].d=d;n++; } } }
  qsort(R,n,sizeof(row),cmp); FILE*f=fopen(argv[3],"w"); fprintf(f,"Q,a,b,c,d\n"); size_t u=0;
  for(size_t i=0;i<n;i++){ if(i&&!cmp(&R[i],&R[i-1])) continue; fprintf(f,"%ld,%ld,%ld,%ld,%ld\n",R[i].q,R[i].a,R[i].b,R[i].c,R[i].d); u++; }
  fclose(f); printf("%zu rows, %ld four-square points\n",u,four); return 0; }
