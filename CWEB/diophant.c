#define BLAS 1
#define USE_SSE 0
#define DEEPINSERT 1
#define DEEPINSERT_CONST 100
#define VERBOSE 1 \

#define GIVENS 1
#define LASTLINESFACTOR "1000000"
#define EPSILON 0.00001
#define LLLCONST_LOW 0.75
#define LLLCONST_HIGH 0.90
#define LLLCONST_HIGHER 0.999 \

#define SQRT sqrt
#define DOUBLE double
#define COEFF struct coe \

#define ROUND(r) ceil(r-0.5)  \

#define FINCKEPOHST 1
#define EIGENBOUND 0 \

/*2:*/
#line 45 "./diophant.w"

/*5:*/
#line 170 "./diophant.w"

#include <stdio.h> 
#include <time.h> 
#include <stdlib.h> 
#include <stdint.h> 
#include <string.h> 
#include <malloc.h> 
#include <math.h> 
#include <gmp.h> 
#include "diophant.h"
#if USE_SSE

#include<pmmintrin.h> 
#endif

#if BLAS
#include "OpenBLASsub/common.h"
#include "OpenBLASsub/cblas.h"
#endif

/*:5*/
#line 46 "./diophant.w"
;
/*6:*/
#line 191 "./diophant.w"

struct coe{
mpz_t c;
int p;
};

/*:6*/
#line 47 "./diophant.w"
;
/*7:*/
#line 201 "./diophant.w"

mpz_t matrix_factor;
mpz_t max_norm;
mpz_t max_norm_initial;
mpz_t max_up;
mpz_t dummy;
long nom,denom;
mpz_t lastlines_factor;

mpz_t snd_q,snd_r,snd_s;



/*:7*//*8:*/
#line 215 "./diophant.w"

int system_rows,system_columns;
int lattice_rows,lattice_columns;
COEFF**lattice;
int free_RHS;
int iszeroone;
mpz_t*upperbounds;
mpz_t upperbounds_max;
mpz_t upfac;

/*:8*//*9:*/
#line 226 "./diophant.w"

int*original_columns;
int no_original_columns;
int cut_after_coeff;
long stop_after_solutions;
long stop_after_loops;
long nosolutions;
int iterate;
int no_iterates;
int bkz_beta,bkz_p;
int SILENT;
int nboundvars;

/*:9*//*41:*/
#line 838 "./diophant.w"

mpz_t soltest_u;
mpz_t soltest_s;
mpz_t soltest_upfac;
/*:41*//*119:*/
#line 2900 "./diophant.w"

static FILE*fp;

/*:119*/
#line 48 "./diophant.w"
;
/*10:*/
#line 243 "./diophant.w"

#define put_to(i,j,val) mpz_set(lattice[i][j+1].c,val)
#define smult_lattice(i,j,factor) mpz_mul(lattice[i][j+1].c,lattice[i][j+1].c,factor)
#define get_entry(i,j) lattice[i][j+1].c

/*:10*/
#line 49 "./diophant.w"
;
/*28:*/
#line 563 "./diophant.w"

/*135:*/
#line 3651 "./diophant.w"

void stopProgram(){
#ifndef NO_OUTPUT
printf("Stopped after SIGALRM, number of solutions: %ld\n",nosolutions);
#endif
/*118:*/
#line 2896 "./diophant.w"

if(SILENT)fprintf(fp,"%ld solutions\n",nosolutions);
fflush(fp);

/*:118*/
#line 3656 "./diophant.w"
;
exit(11);
}

/*:135*/
#line 564 "./diophant.w"

/*29:*/
#line 575 "./diophant.w"

void debug_print(char*m,int l){
#ifndef NO_OUTPUT
if(VERBOSE>=l){
printf("debug>> %s\n",m);fflush(stdout);
}
#endif
return;
}

/*:29*/
#line 565 "./diophant.w"
;
/*30:*/
#line 587 "./diophant.w"

#if 1
void print_lattice(){
int i,j;
for(i= 0;i<lattice_columns;i++){
for(j= 0;j<lattice_rows;j++){
mpz_out_str(NULL,10,get_entry(i,j));
printf(" ");
}
printf("\n");
}
printf("\n");fflush(stdout);
return;
}
#else
void print_lattice(){
int i,j;
for(j= 0;j<lattice_rows;j++){
for(i= 0;i<lattice_columns;i++){
mpz_out_str(NULL,10,get_entry(i,j));
printf(" ");
}
printf("\n");
}
printf("\n");fflush(stdout);
return;
}
#endif

/*:30*/
#line 566 "./diophant.w"
;
/*37:*/
#line 769 "./diophant.w"

void shufflelattice(){
COEFF*tmp;
int i,j,r;
unsigned int s;

#if 1
s= (unsigned)(time(0))*getpid();
#else
s= 1300964772;
#endif
fprintf(stderr,"Seed=%u\n",s);
srand(s);

for(j= 0;j<100;j++){
for(i= lattice_columns-2;i> 0;i--){
r= rand()%i;
tmp= lattice[r];
lattice[r]= lattice[i];
lattice[i]= tmp;
}
}
return;
}

/*:37*/
#line 567 "./diophant.w"
;
/*31:*/
#line 618 "./diophant.w"

void print_NTL_lattice(){
int i,j;
#if 1
fprintf(stderr,"%d %d\n",lattice_columns,lattice_rows);
printf("%d\n",system_rows);
printf("\n[");
for(i= 0;i<lattice_columns-1;i++){
printf("[");
for(j= 0;j<lattice_rows;j++){
mpz_out_str(NULL,10,get_entry(i,j));
printf(" ");
}
printf("]");
printf("\n");
}
printf("]\n");fflush(stdout);
#if 1
printf("\n");
printf("%d ",lattice_columns-2);
mpz_out_str(NULL,10,upperbounds_max);
printf("\n\n[");
for(i= 0;i<lattice_columns-2;i++){
mpz_out_str(NULL,10,upperbounds[i]);
printf(" ");
}
printf("]\n");fflush(stdout);
#endif
#endif
return;
}

/*:31*/
#line 568 "./diophant.w"
;
/*34:*/
#line 709 "./diophant.w"

long gcd(long n1,long n2){
long a,b,c;

if(n1> n2){
a= n1;b= n2;
}else{
a= n2;b= n1;
}

while((c= a%b)> 0){
a= b;b= c;
}
return b;
}
/*:34*/
#line 569 "./diophant.w"
;
/*35:*/
#line 725 "./diophant.w"

void coeffinit(COEFF*v,int z)
{
short r= 0;
short i;
for(i= z;i>=0;i--){
v[i].p= r;
if(mpz_sgn(v[i].c)!=0)r= i;
}
return;
}

/*:35*/
#line 570 "./diophant.w"
;
/*32:*/
#line 659 "./diophant.w"

void read_NTL_lattice(){
int i,j,cols,rows,nbounds;

scanf("%d %d\n",&rows,&cols);

for(i= 0;i<rows;i++){
for(j= 0;j<cols;j++){
mpz_inp_str(lattice[i][j+1].c,NULL,10);
}
}
for(j= 0;j<cols;j++){
mpz_set_ui(lattice[rows][j+1].c,0);
}

scanf("%d",&nbounds);
mpz_inp_str(upperbounds_max,NULL,10);
for(j= 0;j<nbounds;j++){
mpz_inp_str(upperbounds[j],NULL,10);
}

lattice_columns= rows+1;
lattice_rows= cols;
for(i= 0;i<lattice_columns;i++)coeffinit(lattice[j],lattice_rows);

return;
}


/*:32*/
#line 571 "./diophant.w"
;
/*36:*/
#line 742 "./diophant.w"

int cutlattice(){
int j,i,flag;

/*38:*/
#line 795 "./diophant.w"

j= 0;
do{
if(lattice[j][0].p> system_rows)
j++;
else{
for(i= j+1;i<lattice_columns;i++)lattice[i-1]= lattice[i];
lattice_columns--;
}
}while(j<lattice_columns-1);

/*:38*/
#line 746 "./diophant.w"
;
/*39:*/
#line 809 "./diophant.w"

flag= 0;
for(i= 0;i<lattice_columns;i++)if(mpz_sgn(get_entry(i,lattice_rows-1))!=0){
flag= 1;
break;
}
if(flag==0){
#ifndef NO_OUTPUT
printf("Nonhomogenous solution not possible.\n");fflush(stdout);
#endif
exit(2);
return 0;
}
/*:39*/
#line 747 "./diophant.w"
;

for(j= 0;j<lattice_columns;j++){
if(nboundvars==0){
for(i= system_rows;i<lattice_rows;i++)
put_to(j,i-system_rows,get_entry(j,i));
}else{
for(i= system_rows;i<system_rows+nboundvars;i++)
put_to(j,i-system_rows,get_entry(j,i));
for(i= system_rows+system_columns;i<lattice_rows;i++)
put_to(j,i-system_rows-system_columns+nboundvars,get_entry(j,i));
}
}
lattice_rows-= system_rows;
lattice_rows-= (system_columns-nboundvars);

for(j= 0;j<lattice_columns;j++)coeffinit(lattice[j],lattice_rows);

return 1;
}

/*:36*/
#line 572 "./diophant.w"
;
/*40:*/
#line 823 "./diophant.w"

int solutiontest(int position){
int i,j;
int low,up;
int end;

/*43:*/
#line 848 "./diophant.w"

if(mpz_cmpabs(get_entry(position,lattice_rows-1),max_norm)!=0)return 0;
if(mpz_sgn(get_entry(position,lattice_rows-1-free_RHS))==0)return 0;
/*:43*/
#line 829 "./diophant.w"
;
/*44:*/
#line 858 "./diophant.w"

low= 0;
up= lattice_rows-1-free_RHS;
if(lattice_columns==system_columns+2+free_RHS){
for(i= 0;i<system_rows;i++)
if(mpz_sgn(get_entry(position,i))!=0)return 0;
low= system_rows;
}

if(iszeroone){
for(i= low;i<up;i++){
if(mpz_cmpabs(get_entry(position,i),max_norm)!=0)return 0;
}
}else{
for(i= low;i<up;i++){
if(mpz_cmpabs(get_entry(position,i),max_norm)> 0)return 0;
}
}
/*:44*/
#line 830 "./diophant.w"
;

mpz_set_si(upfac,1);
mpz_divexact(soltest_s,get_entry(position,lattice_rows-1),lastlines_factor);
/*45:*/
#line 879 "./diophant.w"

i= low;
if(cut_after_coeff==-1){
end= no_original_columns;
#if 0
if(nboundvars!=0){
end= nboundvars;
}
#endif
}else{
end= cut_after_coeff;
}

for(j= 0;j<end;j++){
if(original_columns[j]==0){
mpz_set_si(soltest_u,0);
}else{
if(!iszeroone){
if(mpz_cmp_si(upperbounds[i-low],0)!=0){
mpz_divexact(soltest_upfac,upperbounds_max,upperbounds[i-low]);
}else{
mpz_set(soltest_upfac,upperbounds_max);
}
}
mpz_set(soltest_u,get_entry(position,i));
mpz_sub(soltest_u,soltest_u,soltest_s);
mpz_divexact(soltest_u,soltest_u,max_norm_initial);
mpz_divexact(soltest_u,soltest_u,soltest_upfac);
mpz_divexact_ui(soltest_u,soltest_u,denom);
mpz_abs(soltest_u,soltest_u);
i++;
}
mpz_out_str(NULL,10,soltest_u);
printf(" ");
if(stop_after_solutions==1){
mpz_out_str(fp,10,soltest_u);
fprintf(fp," ");
}
}
if(free_RHS){
mpz_divexact(soltest_u,get_entry(position,up),max_up);
mpz_divexact(soltest_u,soltest_u,lastlines_factor);
mpz_abs(soltest_u,soltest_u);
printf(" L = ");
mpz_out_str(NULL,10,soltest_u);
}
printf("\n");fflush(stdout);
if(stop_after_solutions==1)fprintf(fp,"\n");

/*:45*/
#line 834 "./diophant.w"
;
/*46:*/
#line 929 "./diophant.w"

if(stop_after_solutions==1){
#ifndef NO_OUTPUT
printf("Stopped in phase 1 after finding a random solution\n");
#endif
exit(8);
}
/*:46*/
#line 835 "./diophant.w"
;
return 1;
}
/*:40*/
#line 573 "./diophant.w"
;

/*:28*/
#line 50 "./diophant.w"
;
/*47:*/
#line 940 "./diophant.w"

/*65:*/
#line 1275 "./diophant.w"

/*66:*/
#line 1289 "./diophant.w"

DOUBLE scalarproductlfp(COEFF*v,COEFF*w)
{
DOUBLE erg;
long t1,t2;
COEFF*vv,*ww;

erg= 0.0;
t1= v[0].p;
t2= w[0].p;
if((t1==0)||(t2==0))return 0;

do{
if(t2> t1){
t1= v[t2-1].p;
if(t2!=t1){
if(t1==0)break;
t2= w[t2].p;
if(t2==0)break;
}
else goto gleich;
}
else if(t2<t1){
t2= w[t1-1].p;
if(t2!=t1){
if(t2==0)break;
t1= v[t1].p;
if(t1==0)break;
}
else goto gleich;
}
else{
gleich:vv= &(v[t1]);
ww= &(w[t2]);
erg+= (DOUBLE)mpz_get_d(vv->c)*(DOUBLE)mpz_get_d(ww->c);
t1= vv->p;
if(t1==0)break;
t2= ww->p;
if(t2==0)break;
}
}
while(1);

return(erg);
}
/*:66*/
#line 1276 "./diophant.w"
;
/*67:*/
#line 1334 "./diophant.w"

DOUBLE scalarproductfp(DOUBLE*v,DOUBLE*w,int n)
{
#if BLAS
return cblas_ddot(n,v,1,w,1);
#else
DOUBLE r;
int i;
r= 0.0;
for(i= n-1;i>=0;i--)r+= v[i]*w[i];
return r;
#endif
}

/*:67*/
#line 1277 "./diophant.w"
;
/*68:*/
#line 1351 "./diophant.w"

int lllalloc(DOUBLE***mu,DOUBLE**c,DOUBLE**N,DOUBLE***bs,int s,int z)
{
int i,m;
#if USE_SSE
int zeven;
#endif

if((z<1)||(s<1))return 0;

(*c)= (DOUBLE*)calloc(s,sizeof(DOUBLE));
(*N)= (DOUBLE*)calloc(s,sizeof(DOUBLE));
(*mu)= (DOUBLE**)calloc(s,sizeof(DOUBLE*));
for(i= 0;i<s;i++)(*mu)[i]= (DOUBLE*)calloc(z,sizeof(DOUBLE));

m= (z> s)?z:s;
(*bs)= (DOUBLE**)calloc(m,sizeof(DOUBLE*));

#if USE_SSE
zeven= (m%8!=0)?(m/8+1)*8:m;
for(i= 0;i<m;i++)
(*bs)[i]= (DOUBLE*)calloc(zeven,sizeof(DOUBLE));
#else
for(i= 0;i<m;i++)(*bs)[i]= (DOUBLE*)calloc(z,sizeof(DOUBLE));
#endif

return 1;
}
/*:68*/
#line 1278 "./diophant.w"
;
/*69:*/
#line 1379 "./diophant.w"

int lllfree(DOUBLE**mu,DOUBLE*c,DOUBLE*N,DOUBLE**bs,int s)
{
int i;

for(i= 0;i<s;++i)free(bs[i]);
free(bs);
for(i= 0;i<s;++i)free(mu[i]);
free(mu);
free(N);
free(c);

return 1;
}

/*:69*/
#line 1279 "./diophant.w"
;

/*:65*/
#line 941 "./diophant.w"
;
/*48:*/
#line 954 "./diophant.w"

#define TWOTAUHALF 67108864.0 
int lllfp(COEFF**b,DOUBLE**mu,DOUBLE*c,DOUBLE*N,DOUBLE**bs,
int start,int s,int z,DOUBLE delta)
{
/*50:*/
#line 1014 "./diophant.w"

int i,j,k;
DOUBLE ss;
int counter;
/*:50*//*53:*/
#line 1068 "./diophant.w"

int Fc,Fr;
DOUBLE mus,cc;
mpz_t musvl;
mpz_t hv;
DOUBLE*swapd;

/*:53*//*55:*/
#line 1089 "./diophant.w"

int ii,iii;
COEFF*swapvl;
COEFF*bb;
/*:55*//*62:*/
#line 1206 "./diophant.w"

int Fi;

/*:62*/
#line 959 "./diophant.w"
;
mpz_init(musvl);
mpz_init(hv);

if((z<=1)||(s<=1)){
#ifndef NO_OUTPUT
printf("Wrong dimensions in lllfp\n");fflush(stdout);
#endif
return(0);
}

k= (start> 1)?start:1;

/*49:*/
#line 1004 "./diophant.w"

if(k<1)k= 1;
for(i= k-1;i<s;++i){
ss= 0.0;
for(j= 0;j<z;++j){
bs[i][j]= (DOUBLE)mpz_get_d(b[i][j+1].c);
ss+= bs[i][j]*bs[i][j];
}
N[i]= SQRT(ss);
}
/*:49*/
#line 972 "./diophant.w"
;
counter= 0;
while(k<s){
#if VERBOSE >  3
if((counter%500)==0){
#ifndef NO_OUTPUT
printf("LLL: %d k:%d\n",counter,k);fflush(stdout);
#endif
}
counter++;
#endif
/*51:*/
#line 1021 "./diophant.w"

if(k==1)c[0]= N[0]*N[0];
c[k]= N[k]*N[k];
for(j= 0;j<k;j++){
ss= scalarproductfp(bs[k],bs[j],z);
if(fabs(ss)<N[k]*N[j]/TWOTAUHALF){
ss= (DOUBLE)scalarproductlfp(b[k],b[j]);
}
for(i= 0;i<j;i++)ss-= mu[j][i]*mu[k][i]*c[i];
if(c[j]<EPSILON){
fprintf(stderr,"c[%d] is very small: %lf\n",j,c[j]);
}
mu[k][j]= ss/c[j];
c[k]-= ss*mu[k][j];
}
/*:51*/
#line 983 "./diophant.w"
;
/*52:*/
#line 1050 "./diophant.w"

Fc= Fr= 0;
for(j= k-1;j>=0;j--){
if(fabs(mu[k][j])> 0.5){
/*54:*/
#line 1078 "./diophant.w"

mus= ROUND(mu[k][j]);
mpz_set_d(musvl,mus);
if(fabs(mus)> TWOTAUHALF){
Fc= 1;
#if 0
printf("correct possible rounding errors\n");fflush(stdout);
#endif
}


/*:54*/
#line 1054 "./diophant.w"
;
Fr= 1;
/*56:*/
#line 1102 "./diophant.w"

switch(mpz_get_si(musvl)){
case 1:
/*57:*/
#line 1117 "./diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
mpz_sub(bb->c,bb->c,b[j][i].c);
iii= bb->p;
if((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
for(ii= i-1;(ii>=0)&&(b[k][ii].p==iii);ii--)b[k][ii].p= i;
else if(mpz_sgn(bb->c)==0){
for(ii= i-1;(ii>=0)&&(b[k][ii].p==i);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
for(i= 0;i<j;i++)mu[k][i]-= mu[j][i];
/*:57*/
#line 1105 "./diophant.w"
;
break;

case-1:
/*58:*/
#line 1133 "./diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
mpz_add(bb->c,bb->c,b[j][i].c);
iii= bb->p;
if((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
for(ii= i-1;(ii>=0)&&(b[k][ii].p==iii);ii--)b[k][ii].p= i;
else if(mpz_sgn(bb->c)==0){
for(ii= i-1;(ii>=0)&&(b[k][ii].p==i);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
for(i= 0;i<j;i++)mu[k][i]+= mu[j][i];
/*:58*/
#line 1109 "./diophant.w"
;
break;

default:
/*59:*/
#line 1149 "./diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
mpz_submul(bb->c,b[j][i].c,musvl);
iii= bb->p;
if((b[k][i-1].p!=i)&&(mpz_sgn(bb->c)!=0))
for(ii= i-1;(ii>=0)&&(b[k][ii].p==iii);ii--)b[k][ii].p= i;
else if(mpz_sgn(bb->c)==0){
for(ii= i-1;(ii>=0)&&(b[k][ii].p==i);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
#if 0
daxpy(j,-mus,mu[k],1,mu[j],1);
#endif
for(i= 0;i<j;i++)mu[k][i]-= mu[j][i]*mus;
/*:59*/
#line 1113 "./diophant.w"
;
}
/*:56*/
#line 1056 "./diophant.w"
;
mu[k][j]-= mus;
solutiontest(k);
}
}
/*60:*/
#line 1166 "./diophant.w"

{
N[k]= 0.0;
for(i= 0;i<z;i++){
bs[k][i]= (DOUBLE)mpz_get_d(b[k][i+1].c);
N[k]+= bs[k][i]*bs[k][i];
}
N[k]= SQRT(N[k]);
}
/*:60*/
#line 1061 "./diophant.w"
;
if(Fc==1){
k= (k-1> 1)?k-1:1;
}else{
/*61:*/
#line 1180 "./diophant.w"

if(N[k]<-EPSILON){
#ifndef NO_OUTPUT
fprintf(stderr,"Nk negativ! contact the author.\n");
fflush(stderr);
printf("Nk negativ! contact the author.\n");fflush(stdout);
#endif
exit(1);
}
if(N[k]<0.5){
swapvl= b[k];
ss= N[k];
swapd= bs[k];
for(i= k+1;i<s;i++){
b[i-1]= b[i];
N[i-1]= N[i];
bs[i-1]= bs[i];
}
b[s-1]= swapvl;
N[s-1]= ss;
bs[s-1]= swapd;
s= s-1;
k= 1;
continue;
}

/*:61*/
#line 1065 "./diophant.w"
;
}

/*:52*/
#line 984 "./diophant.w"
;
#if defined(DEEPINSERT)
/*63:*/
#line 1217 "./diophant.w"

cc= N[k]*N[k];
j= 0;
Fi= 0;
while(j<k){
#if 1
if((j> DEEPINSERT_CONST&&j<k-DEEPINSERT_CONST)||delta*c[j]<=cc){
#else
if(delta*c[j]<=cc){
#endif
cc-= mu[k][j]*mu[k][j]*c[j];
j++;
}else{
swapvl= b[k];
ss= N[k];
swapd= bs[k];
for(i= k-1;i>=j;i--){
b[i+1]= b[i];
N[i+1]= N[i];
bs[i+1]= bs[i];
}
b[j]= swapvl;
N[j]= ss;
bs[j]= swapd;

Fi= 1;
break;
}
}
if(Fi==1)k= (j-1> 1)?j-1:1;
else{
k++;
}

/*:63*/
#line 986 "./diophant.w"
;
#else
/*64:*/
#line 1258 "./diophant.w"

if(delta*c[k-1]> c[k]+mu[k][k-1]*mu[k][k-1]*c[k-1]){
swapvl= b[k];
b[k]= b[k-1];
b[k-1]= swapvl;
ss= N[k];
N[k]= N[k-1];
N[k-1]= ss;
swapd= bs[k];
bs[k]= bs[k-1];
bs[k-1]= swapd;

k= (k-1> 1)?k-1:1;
}else k++;

/*:64*/
#line 988 "./diophant.w"
;
#endif
}
mpz_clear(hv);
mpz_clear(musvl);
return(1);

}
/*:48*/
#line 942 "./diophant.w"
;
/*70:*/
#line 1395 "./diophant.w"

double logD(COEFF**lattice,DOUBLE*c,int s,int z){
double d= 0.0;
int i;
for(i= 0;i<s;i++){
d+= log(c[i])*(s-i);
}
d*= 0.5;
return d;
}

/*:70*/
#line 943 "./diophant.w"
;
/*71:*/
#line 1409 "./diophant.w"

double orthogonal_defect(COEFF**lattice,DOUBLE*c,int s,int z){
double defect= 0.0;

#if 0
int i;
for(i= 0;i<s;i++)defect+= log((double)normfp(lattice[i]))
-log((double)c[i]);
#endif
defect/= 2.0;
return defect;
}

/*:71*/
#line 944 "./diophant.w"
;
/*72:*/
#line 1423 "./diophant.w"

void lll(COEFF**b,int s,int z,DOUBLE quality)
{
DOUBLE**mu;
DOUBLE*c;
DOUBLE*N;
DOUBLE**bs;
int r;

lllalloc(&mu,&c,&N,&bs,s,z);
r= lllfp(b,mu,c,N,bs,1,s,z,quality);
lllfree(mu,c,N,bs,s);

return;
}

/*:72*/
#line 945 "./diophant.w"
;
/*73:*/
#line 1443 "./diophant.w"

DOUBLE iteratedlll(COEFF**b,int s,int z,int no_iterates,DOUBLE quality)
{
DOUBLE**mu;
DOUBLE*c;
DOUBLE*N;
DOUBLE**bs;
int r,l,i,j,runs;
COEFF*swapvl;
DOUBLE lD;

lllalloc(&mu,&c,&N,&bs,s,z);
r= lllfp(b,mu,c,N,bs,1,s,z,quality);

lD= logD(b,c,s,z);
#ifndef NO_OUTPUT
printf("   log(D)= %f\n",lD);
fflush(stdout);
#endif

for(runs= 1;runs<no_iterates;runs++){


for(j= s-1;j> 0;j--){
for(l= j-1;l>=0;l--){


if(N[l]> N[j]){
swapvl= b[l];
for(i= l+1;i<=j;i++)b[i-1]= b[i];
b[j]= swapvl;
}
}
}

r= lllfp(b,mu,c,N,bs,1,s,z,quality);
lD= logD(b,c,s,z);
#ifndef NO_OUTPUT
printf("%d: log(D)= %f\n",runs,lD);
fflush(stdout);
#endif
}

lllfree(mu,c,N,bs,s);

return lD;
}

/*:73*/
#line 946 "./diophant.w"
;
/*74:*/
#line 1492 "./diophant.w"

/*85:*/
#line 1912 "./diophant.w"

DOUBLE laurin(DOUBLE x){
static DOUBLE K1= 0.9181938533204672741780329736405620;
static DOUBLE K2= 0.0833333333333333333333333333333333333;
static DOUBLE K3= 0.0027777777777777777777777777777777777;
static DOUBLE K4= 0.000793650793650793650793650793650793650;
static DOUBLE K5= 0.0005952380952380952380952380952380952380;
static DOUBLE K6= 0.000841750841750841750841750841750841750;
static DOUBLE K7= 0.001917526917526917526917526917526917526;
static DOUBLE K8= 0.00641025641025641025641025641025641025;
static DOUBLE K9= 0.0295506529510021209716796875;
static DOUBLE K10= 0.17968122661113739013671875;
static DOUBLE K11= 1.39243221282958984375;
DOUBLE y;

y= 1.0/(x*x);
y= (x-0.5)*log(x)-x+K1+
(1.0/x)*(K2-y*(K3-y*(K4-y*(K5-y*(K6-y*(K7-y*(K8-y*(K9-y*(K10-y*K11)))))))));
return y;
}

DOUBLE log_gamma(DOUBLE x){
DOUBLE y;
int i,n;
static int MM= 13;

if(x<=0.0)return-1.0;

if(x> 100000000.0){
y= x*(log(x)-1.0);
}else{
if(x>=MM){
y= laurin(x);
}else{
n= MM-(int)(floor(x));
y= x-floor(x)+MM;
y= laurin(y);
for(i= 0;i<n;i++)y-= log(y+i);
}
}
return y;
}

/*:85*/
#line 1493 "./diophant.w"
;
/*82:*/
#line 1703 "./diophant.w"

DOUBLE enumerate(DOUBLE**mu,DOUBLE*c,long*u,int s,
int start_block,int end_block,int p)
{
DOUBLE cd,dum;
DOUBLE*y,*cs,*eta;

DOUBLE**sigma;
int*r;

long*us,*delta,*d,*v;
int t,i,t_up,len;
double alpha;
int tmax;
static DOUBLE pi= 3.141592653589793238462643383;
static DOUBLE e= 2.718281828459045235360287471352662497757247093;
int SCHNITT= 40;

if(c[start_block]<=EPSILON){
#ifndef NO_OUTPUT
fprintf(stderr,"Hier ist was faul! start_block=%d %f\n",start_block,
(double)c[start_block]);fflush(stderr);
printf("Hier ist was faul! start_block=%d %f\n",start_block,
(double)c[start_block]);fflush(stdout);
#endif
exit(1);
}

us= (long*)calloc(s+1,sizeof(long));
cs= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
y= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
delta= (long*)calloc(s+1,sizeof(long));
d= (long*)calloc(s+1,sizeof(long));
eta= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
v= (long*)calloc(s+1,sizeof(long));

sigma= (DOUBLE**)calloc(s,sizeof(DOUBLE*));
r= (int*)calloc(s+1,sizeof(int));
for(i= 0;i<s;i++){
sigma[i]= (DOUBLE*)calloc(s,sizeof(DOUBLE));
r[i]= i-1;
}


len= (end_block+1-start_block);
for(i= start_block;i<=end_block+1;i++){
cs[i]= y[i]= 0.0;
u[i]= us[i]= v[i]= delta[i]= 0;
d[i]= 1;
}
us[start_block]= u[start_block]= 1;
cd= c[start_block];

t= tmax= start_block;


/*84:*/
#line 1871 "./diophant.w"

eta[start_block]= 0.0;
if(end_block-start_block<=SCHNITT){
for(i= start_block+1;i<=end_block;i++)eta[i]= 0.0;
}
else{
dum= log(c[start_block]);



#if 0
for(i= start_block+1;i<=end_block;i++){
t_up= i-start_block;
eta[i]= exp((log_gamma(t_up/2.0+1)-p*log(2.0))*2.0/t_up+dum/t_up)/pi;
if(i<end_block)dum+= log(c[i]);
#if 0
if(i==start_block+1)printf("eta: \n");
printf("%0.2lf ",eta[i]);
if(i==end_block)printf("\n");
#endif
}
#else



dum= log(c[start_block]);
for(i= start_block+1;i<=end_block;i++){
t_up= i-start_block;
eta[i]= 0.5*t_up*exp((log(pi*t_up)-2.0*p*log(2.0)+dum)/t_up)/(pi*e);
if(i<end_block)dum+= log(c[i]);
#if 0
if(i==start_block+1)printf("eta: ");
printf("%0.2lf ",eta[i]);
if(i==end_block){
printf("\n");
fflush(stdout);
}
#endif
}
#endif
}
/*:84*/
#line 1759 "./diophant.w"
;
while(t<=end_block){
/*83:*/
#line 1779 "./diophant.w"

{
dum= us[t]+y[t];
cs[t]= cs[t+1]+dum*dum*c[t];

#if 0
if(cs[t]<cd-eta[t]){

#else
if(len<=SCHNITT){
alpha= 1.0;
}else{
alpha= sqrt(1.20*(end_block+1-t)/len);
if(alpha>=1.0)alpha= 1.0;
}
alpha*= cd;



if(cs[t]<alpha-EPSILON){
#endif
if(t> start_block){
t--;
if(r[t+1]> r[t])r[t]= r[t+1];

delta[t]= 0;
for(i= r[t+1];i> t;i--)sigma[i][t]= sigma[i+1][t]+us[i]*mu[i][t];
#if 0
dum= 0.0;
for(i= t+1;i<=tmax;i++)dum+= us[i]*mu[i][t];
if(fabs(dum-sigma[t+1][t])> 0.001){
printf("1diff: %0.6lf %0.6lf %0.8lf\n",dum,sigma[t+1][t],dum-sigma[t+1][t]);
}
#endif
y[t]= sigma[t+1][t];
us[t]= v[t]= (long)(ROUND(-y[t]));
d[t]= (v[t]> -y[t])?-1:1;
}else{
#if 0
printf("success %0.4lf %0.4lf %0.8lf\n",cd,cs[start_block],cd-cs[start_block]);
#endif
cd= cs[start_block];
for(i= start_block;i<=end_block;i++)u[i]= us[i];
goto nextstep;
}
}else{
t++;
r[t]= t;
nextstep:
if(tmax<t)tmax= t;
if(t<tmax)delta[t]= -delta[t];
if(delta[t]*d[t]>=0)delta[t]+= d[t];
us[t]= v[t]+delta[t];
}
}

/*:83*/
#line 1761 "./diophant.w"
;
}
free(us);
free(cs);
free(y);
free(delta);
free(d);
free(eta);
free(v);
for(i= s-1;i>=0;i--){
free(sigma[i]);
}
free(sigma);
free(r);
return(cd);
}

/*:82*/
#line 1494 "./diophant.w"
;
/*75:*/
#line 1513 "./diophant.w"

DOUBLE bkz(COEFF**b,int s,int z,DOUBLE delta,int beta,int p)
{
/*76:*/
#line 1574 "./diophant.w"

DOUBLE**mu,*c,*N;
DOUBLE**bs;
static mpz_t hv;
int zaehler;
int h,i,last;
int start_block,end_block;
long*u;
DOUBLE new_cj;
DOUBLE lD;

mpz_init(hv);

/*:76*//*78:*/
#line 1617 "./diophant.w"

int g,ui,q,j;
COEFF*swapvl;

/*:78*/
#line 1516 "./diophant.w"
;

last= s-2;
if(last<1){
#ifndef NO_OUTPUT
printf("BKZ: the number of basis vectors is too small.\n");
printf("Probably the number of rows is less or equal");
printf(" to number of columns in the original system\n");
printf("Maybe you have to increase c0 (the first parameter)!\n");
#endif
mpz_clear(hv);
return 0.0;
}

u= (long*)calloc(s,sizeof(long));
for(i= 0;i<s;i++)u[i]= 0;

lllalloc(&mu,&c,&N,&bs,s,z);
lllfp(b,mu,c,N,bs,1,s,z,delta);

start_block= zaehler= -1;
while(zaehler<last){

start_block++;
if(start_block==last)start_block= 0;
end_block= (start_block+beta-1<last)?start_block+beta-1:last;

#if 0
printf("start_block=%d, end_block=%d\n",start_block,end_block);
#endif
new_cj= enumerate(mu,c,u,s,start_block,end_block,p);

h= (end_block+1<last)?end_block+1:last;


if(delta*c[start_block]> new_cj){
/*77:*/
#line 1597 "./diophant.w"

#if defined(ORIGINAL_SCHNORR_EUCHNER)
/*80:*/
#line 1673 "./diophant.w"

for(l= 1;l<=z;l++)mpz_set_si(b[last+1][l].c,0);
for(i= start_block;i<=end_block;i++)
for(l= 1;l<=z;l++){
mpz_addmul_si(b[last+1][l].c,b[i][l].c,ui);
}
coeffinit(b[last+1],z);
solutiontest(last+1);
swapvl= b[last+1];
for(i= last;i>=start_block;i--)b[i+1]= b[i];
b[start_block]= swapvl;

/*:80*/
#line 1599 "./diophant.w"
;
#else
/*79:*/
#line 1628 "./diophant.w"

for(j= 1;j<=z;j++)mpz_set_si(b[last+1][j].c,0);
for(i= start_block;i<=end_block;i++){
if(u[i]!=0)for(j= 1;j<=z;j++){
if(u[i]> 0){
mpz_addmul_ui(b[last+1][j].c,b[i][j].c,u[i]);
}else{
mpz_submul_ui(b[last+1][j].c,b[i][j].c,-u[i]);
}
}
}
g= end_block;
while(u[g]==0)g--;

i= g-1;
while(labs(u[g])> 1){
while(u[i]==0)i--;
q= (int)ROUND((1.0*u[g])/u[i]);
ui= u[i];
u[i]= u[g]-q*u[i];
u[g]= ui;

for(j= 1;j<=z;j++){
mpz_set(hv,b[g][j].c);
mpz_mul_si(b[g][j].c,b[g][j].c,(long)q);
mpz_add(b[g][j].c,b[g][j].c,b[i][j].c);
mpz_set(b[i][j].c,hv);
}
coeffinit(b[g],z);
coeffinit(b[i],z);
}

swapvl= b[g];
for(i= g;i> start_block;i--)b[i]= b[i-1];
b[start_block]= b[last+1];
coeffinit(b[start_block],z);

b[last+1]= swapvl;
for(j= 1;j<=z;j++)mpz_set_si(b[last+1][j].c,0);
coeffinit(b[last+1],z);

/*:79*/
#line 1601 "./diophant.w"
;
#endif

lllfp(b,mu,c,N,bs,start_block-1,h+1,z,delta);

if(N[h]<-EPSILON){
#ifndef NO_OUTPUT
fprintf(stderr,"NN negativ\n");fflush(stderr);
printf("NN negativ\n");fflush(stdout);
#endif
exit(1);
}
#if defined(ORIGINAL_SCHNORR_EUCHNER)
/*81:*/
#line 1685 "./diophant.w"

if(N[h]<0.5){
swapvl= b[h];
for(i= h+1;i<=last+1;i++)b[i-1]= b[i];
b[last+1]= swapvl;
}else{
#ifndef NO_OUTPUT
printf("Not linear dependent; %f\n",(double)(N[h-1]));fflush(stdout);
#endif
exit(1);
}

/*:81*/
#line 1614 "./diophant.w"
;
#endif

/*:77*/
#line 1552 "./diophant.w"
;
zaehler= -1;
}
else{
if(h> 0){
lllfp(b,mu,c,N,bs,h-2,h+1,z,delta);

}
zaehler++;
}
}

lD= logD(b,c,s-1,z);
#ifndef NO_OUTPUT
printf("bkz: log(D)= %f\n",lD);
fflush(stdout);
#endif
lllfree(mu,c,N,bs,s);
free(u);
mpz_clear(hv);
return lD;
}
/*:75*/
#line 1495 "./diophant.w"
;

/*:74*/
#line 947 "./diophant.w"
;
/*86:*/
#line 1960 "./diophant.w"

/*87:*/
#line 2036 "./diophant.w"

#if 0
static FILE*fp;
#endif
/*:87*//*92:*/
#line 2157 "./diophant.w"

#if FINCKEPOHST
DOUBLE**muinv;
DOUBLE**fipo_UB,**fipo_LB;
#endif

long fipo_success;

/*:92*//*94:*/
#line 2192 "./diophant.w"

#if EIGENBOUND
DOUBLE**Rinv;
DOUBLE**R;
DOUBLE**eig_f;
DOUBLE**eig_RinvR;
DOUBLE eig_Rs;
DOUBLE*eig_min;
DOUBLE*eig_bound;
long eig_cut;

DOUBLE**Rinvi;
DOUBLE**bnd_up;
DOUBLE**bnd_lo;
#endif

/*:94*//*100:*/
#line 2283 "./diophant.w"

DOUBLE dum1,dum2;

/*:100*//*108:*/
#line 2571 "./diophant.w"

long only_zeros_no,only_zeros_success,hoelder_no,hoelder_success;
long cs_success;
long N_success;
long N2_success;
long N3_success;


/*:108*/
#line 1961 "./diophant.w"
;
/*120:*/
#line 2904 "./diophant.w"

DOUBLE compute_y(DOUBLE**mu_trans,DOUBLE*us,int level,int level_max){
#if BLAS
return cblas_ddot(level_max-level,&(us[level+1]),1,&(mu_trans[level][level+1]),1);
#else
int i;
DOUBLE dum;
i= level_max;
dum= 0.0;
while(i>=level+1){
dum+= mu_trans[level][i]*us[i];
i--;
}
return dum;
#endif
}

void compute_w2(DOUBLE**w,DOUBLE**bd,DOUBLE alpha,int level,int rows){
#if BLAS
cblas_daxpy(rows,alpha,bd[level],1,w[level],1);
#else
int i;
for(i= 0;i<rows;++i){
w[level][i]+= alpha*bd[level][i];
}
#endif
return;
}

void compute_w(DOUBLE**w,DOUBLE**bd,DOUBLE alpha,int level,int rows){
#if USE_SSE
__m128d a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
int l;

a= _mm_loaddup_pd(&alpha);
for(l= 0;l<rows;l+= 8){

x1= _mm_load_pd(&(bd[level][l]));
y1= _mm_load_pd(&(w[level+1][l]));

x2= _mm_load_pd(&(bd[level][l+2]));
y2= _mm_load_pd(&(w[level+1][l+2]));

x3= _mm_load_pd(&(bd[level][l+4]));
y3= _mm_load_pd(&(w[level+1][l+4]));

x4= _mm_load_pd(&(bd[level][l+6]));
y4= _mm_load_pd(&(w[level+1][l+6]));

x1= _mm_mul_pd(x1,a);
z1= _mm_add_pd(x1,y1);
_mm_storeu_pd((double*)&(w[level][l]),z1);

x2= _mm_mul_pd(x2,a);
z2= _mm_add_pd(x2,y2);
_mm_storeu_pd((double*)&(w[level][l+2]),z2);

x3= _mm_mul_pd(x3,a);
z3= _mm_add_pd(x3,y3);
_mm_storeu_pd((double*)&(w[level][l+4]),z3);

x4= _mm_mul_pd(x4,a);
z4= _mm_add_pd(x4,y4);
_mm_storeu_pd((double*)&(w[level][l+6]),z4);
}
return;
#else
#if BLAS
cblas_dcopy(rows,w[level+1],1,w[level],1);
cblas_daxpy(rows,alpha,bd[level],1,w[level],1);
#else
int l;

l= rows-1;
while(l>=0){
w[level][l]= w[level+1][l]+alpha*bd[level][l];
l--;
}
#endif
return;
#endif
}

/*:120*/
#line 1962 "./diophant.w"
;
/*121:*/
#line 2989 "./diophant.w"

/*122:*/
#line 2994 "./diophant.w"

void gramschmidt(COEFF**lattice,int columns,int rows,DOUBLE**mu,
DOUBLE**bd,DOUBLE*c,DOUBLE*N,DOUBLE Fq){
int i,l,j;
DOUBLE dum;

for(i= 0;i<columns;i++){
for(l= 0;l<rows;l++)bd[i][l]= (DOUBLE)mpz_get_d(get_entry(i,l));
N[i]= 0.0;
for(j= 0;j<i;j++){
dum= 0.0;
for(l= 0;l<rows;l++)dum+= (DOUBLE)mpz_get_d(get_entry(i,l))*bd[j][l];
mu[i][j]= dum/c[j];
for(l= 0;l<rows;l++)bd[i][l]-= mu[i][j]*bd[j][l];
}

c[i]= scalarproductfp(bd[i],bd[i],rows);
for(l= 0;l<rows;l++)N[i]+= fabs(bd[i][l]);
N[i]/= c[i];
N[i]*= Fq;
#ifndef NO_OUTPUT
#if VERBOSE >  0
printf("%lf ",(double)c[i]);
#endif
#endif

}
#ifndef NO_OUTPUT
#if VERBOSE >  0
printf("\n\n");fflush(stdout);
#endif
#endif
return;
}

/*:122*/
#line 2990 "./diophant.w"
;
/*123:*/
#line 3030 "./diophant.w"

void givens(COEFF**lattice,int columns,int rows,DOUBLE**mu,
DOUBLE**bd,DOUBLE*c,DOUBLE*N,DOUBLE Fq){
int i,l,j;
int mm;
DOUBLE d1,d2;
DOUBLE gc,gs;
DOUBLE t;





for(i= 0;i<columns;i++){
for(l= 0;l<rows;l++){
mu[i][l]= (DOUBLE)mpz_get_d(get_entry(i,l));
}
}

for(i= 0;i<rows;i++){
for(l= 0;l<rows;l++)bd[i][l]= 0.0;
bd[i][i]= 1.0;
}

for(j= 1;j<rows;j++){
mm= (j<columns)?j:columns;
for(i= 0;i<mm;i++){
if(mu[i][j]==0.0){

gc= 1.0;
gs= 0.0;
}else{



if(fabs(mu[i][j])>=fabs(mu[i][i])){
t= mu[i][i]/mu[i][j];
gs= 1.0/SQRT(1.0+t*t);
gc= gs*t;
}else{
t= mu[i][j]/mu[i][i];
gc= 1.0/SQRT(1.0+t*t);
gs= gc*t;
}

for(l= i;l<columns;l++){
d1= mu[l][i];
d2= mu[l][j];
mu[l][i]= gc*d1+gs*d2;
mu[l][j]= -gs*d1+gc*d2;
}

for(l= 0;l<rows;l++){
d1= bd[i][l];
d2= bd[j][l];
bd[i][l]= gc*d1+gs*d2;
bd[j][l]= -gs*d1+gc*d2;
}
}
}
}


for(i= 0;i<columns;i++){
c[i]= mu[i][i]*mu[i][i];
N[i]= 0.0;
for(j= 0;j<rows;j++){
bd[i][j]*= mu[i][i];
N[i]+= fabs(bd[i][j]);
}
N[i]/= c[i];
N[i]*= Fq;



for(j= columns-1;j>=i;j--)mu[j][i]/= mu[i][i];
#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("%6.3f ",(double)c[i]);
if(i> 0&&i%15==0)printf("\n");
#endif
#endif
}
#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("\n\n");fflush(stdout);
#endif
#endif
return;
}
/*:123*/
#line 2991 "./diophant.w"
;

/*:121*/
#line 1963 "./diophant.w"
;
/*124:*/
#line 3121 "./diophant.w"

#if FINCKEPOHST
void inverse(DOUBLE**mu,DOUBLE**muinv,int columns){
int i,j,k;
DOUBLE sum;
for(j= 0;j<columns;j++)
for(i= j;i>=0;i--){
sum= 0.0;
for(k= i+1;k<columns;k++)sum+= mu[k][i]*muinv[k][j];
if(i==j)
muinv[i][j]= 1.0-sum;
else
muinv[i][j]= -sum;
}
return;
}
#endif

/*:124*/
#line 1964 "./diophant.w"
;
/*125:*/
#line 3140 "./diophant.w"

/*126:*/
#line 3150 "./diophant.w"

int exacttest(DOUBLE*v,int rows,DOUBLE Fq){
int i;
i= rows-1;
do{
if(fabs(v[i])> Fq+EPSILON){
return 0;
}
i--;
}while(i>=0);

return 1;
}

/*:126*/
#line 3141 "./diophant.w"
;
/*127:*/
#line 3166 "./diophant.w"

int prune0(DOUBLE li,DOUBLE re){
if(li> re*(1+EPSILON)){
N_success++;
return 1;
}else{
return 0;
}
}

/*:127*/
#line 3142 "./diophant.w"
;
/*128:*/
#line 3177 "./diophant.w"

int prune(DOUBLE*w,DOUBLE cs,int rows,DOUBLE Fqeps){


#if BLAS
if(cs<Fqeps*cblas_dasum(rows,w,1))
return 0;
#else
DOUBLE reseite;
int i;
reseite= 0.0;
i= rows-1;
do{
reseite+= fabs(w[i]);
i--;
}while(i>=0);
if(cs<Fqeps*reseite)return 0;
#endif

return 1;
}

int pruneN(DOUBLE**w,DOUBLE*cs,int t,int rows,int cols,DOUBLE Fq){
int i,t_up;
DOUBLE sum;
DOUBLE r;
if(t>=cols-2)return 0;

if(t<cols/2+10)return 0;

t_up= t+1;
r= 0.4;
sum= 0.0;
for(i= 0;i<rows;i++){
sum+= fabs(w[t][i]-r*w[t_up][i]);
}
sum*= Fq*(1.000001);
sum-= fabs(cs[t]-r*cs[t_up]);
if(sum<0.0){
N2_success++;
#if 0
printf("PRUNEN %d %d %0.3lf\n",t,t_up,r);
#endif
return 1;
}

t_up= t+2;
r= 0.4;
sum= 0.0;
for(i= 0;i<rows;i++){
sum+= fabs(w[t][i]-r*w[t_up][i]);
}
sum*= Fq*(1.000001);
sum-= fabs(cs[t]-r*cs[t_up]);
if(sum<0.0){
N2_success++;
#if 0
printf("PRUNEN %d %d %0.3lf\n",t,t_up,r);
#endif
return 1;
}

t_up= cols-1;
r= 0.2;
sum= 0.0;
for(i= 0;i<rows;i++){
sum+= fabs(w[t][i]-r*w[t_up][i]);
}
sum*= Fq*(1.000001);
sum-= fabs(cs[t]-r*cs[t_up]);
if(sum<0.0){
N2_success++;
#if 0
printf("PRUNEN %d %d %0.3lf\n",t,t_up,r);
#endif
return 1;
}

if(0&&t> cols-20){
r= 0.2;
sum= 0.0;
for(i= 0;i<rows;i++){
sum+= fabs(w[t][i]-r*w[t_up][i]);
}
sum*= Fq*(1.000001);
sum-= fabs(cs[t]-r*cs[t_up]);
if(sum<0.0){
N3_success++;
#if 0
printf("PRUNEN %d %d %0.3lf\n",t,t_up,r);
#endif
return 1;
}
r= 0.3;
sum= 0.0;
for(i= 0;i<rows;i++){
sum+= fabs(w[t][i]-r*w[t_up][i]);
}
sum*= Fq*(1.000001);
sum-= fabs(cs[t]-r*cs[t_up]);
if(sum<0.0){
N3_success++;
#if 0
printf("PRUNEN %d %d %0.3lf\n",t,t_up,r);
#endif
return 1;
}
}
return 0;
}

/*:128*/
#line 3143 "./diophant.w"
;
/*129:*/
#line 3299 "./diophant.w"

int prune_only_zeros(DOUBLE**w,int level,int rows,DOUBLE Fq,
int*first_nonzero_in_column,int*firstp,
DOUBLE**bd,DOUBLE*y,DOUBLE*us,int columns){
int i;
int f;
DOUBLE u1,u2,swp;

only_zeros_no++;
for(i= 0;i<first_nonzero_in_column[firstp[level]];i++){
f= first_nonzero_in_column[firstp[level]+1+i];
u1= (Fq-w[level+1][f])/bd[level][f]-y[level];
u2= (-Fq-w[level+1][f])/bd[level][f]-y[level];

#if 0
if(u2<u1){
swp= u1;
u1= u2;
u2= swp;
}
fipo_LB[columns][level]= u1-EPSILON;
fipo_UB[columns][level]= u2+EPSILON;
#endif

if(iszeroone){
if(fabs(u1-round(u1))> EPSILON&&fabs(u2-round(u2))> EPSILON){
only_zeros_success++;
return-1;
}

if(fabs(fabs(w[level][f])-Fq)> EPSILON){
only_zeros_success++;
return 1;
}

#if 0
if(fabs(u1-us[level])> EPSILON&&fabs(u2-us[level])> EPSILON){
return 1;
}
#endif
}else{


if(u2-u1<=1.0+EPSILON&&
fabs(w[level][f])<UINT32_MAX&&
fabs(w[level][f]-round(w[level][f]))> 0.001){

only_zeros_success++;
return-1;
}

if(fabs(w[level][f])> Fq*(1+EPSILON)){
return 1;
}
#if 0
if(us[level]<u1-EPSILON||us[level]> u2+EPSILON){
return 1;
}
#endif
}
#if 0
if(iszeroone){
if(fabs(fabs(w[level][f])-Fq)> EPSILON){

return 1;
}
}else{
if(fabs(w[level][f])> Fq*(1+EPSILON)){
only_zeros_success++;
return 1;
}
}
#endif
}
return 0;
}

/*:129*/
#line 3144 "./diophant.w"
;
/*130:*/
#line 3377 "./diophant.w"

int prune_snd_nonzero(int columns,int rows,
int level,DOUBLE Fq,
int*first_nonzero,
int*snd_nonzero_in_column,int*sndp,
DOUBLE*us,
struct constraint*cons){
int i,k;
int ro;
int f1;

return 0;
for(i= 0;i<snd_nonzero_in_column[sndp[level]];i++){
ro= snd_nonzero_in_column[sndp[level]+1+i];
f1= first_nonzero[ro];
if(level-f1<=5){
continue;
}
mpz_set_si(soltest_s,0);
for(k= level;k<columns;k++){
if(ROUND(us[k])> 0){
mpz_addmul_ui(soltest_s,get_entry(k,ro),ROUND(us[k]));
}else{
mpz_submul_ui(soltest_s,get_entry(k,ro),-ROUND(us[k]));
}
}
mpz_sub(snd_s,max_norm,soltest_s);
mpz_fdiv_qr(snd_q,snd_r,snd_s,get_entry(f1,ro));
if(mpz_sgn(snd_r)!=0){
printf("Contradiction\n");
return 1;
}
#if 0
if(cons[f1].isSet==0){
cons[f1].isSet= 1;
cons[f1].val[0]= u1;
cons[f1].val[1]= u2;
}
#endif
}
return 0;
}

/*:130*/
#line 3145 "./diophant.w"
;
/*132:*/
#line 3523 "./diophant.w"

int print_solution(DOUBLE*w,int rows,DOUBLE Fq,DOUBLE*us,int columns){
int i,j,k;
int upper;
int end;


#if 1
if(fabs(fabs(w[rows-1])-Fq)> 0.5*Fq*EPSILON){
#else

if(fabs(w[rows-1])> Fq*(1+EPSILON)){
#endif
return 0;
}
upper= rows-1-free_RHS;
#if 0
if(free_RHS&&fabs(fabs(w[upper])-Fq)> 0.5*Fq*EPSILON){
#else
if(free_RHS&&fabs(w[upper])> Fq*(1+EPSILON)){
#endif
return 0;
}

if(!SILENT){
mpz_set_si(soltest_upfac,1);
#if 0
mpz_set_d(soltest_s,ROUND(w[rows-1]));
#else
mpz_set_si(soltest_s,0);
for(k= 0;k<columns;k++){
if(ROUND(us[k])> 0){
mpz_addmul_ui(soltest_s,get_entry(k,rows-1),ROUND(us[k]));
}else{
mpz_submul_ui(soltest_s,get_entry(k,rows-1),-ROUND(us[k]));
}
}
#endif

i= 0;
if(cut_after_coeff==-1){
end= no_original_columns;
#if 0
if(nboundvars!=0){
end= nboundvars;
}
#endif
}else{
end= cut_after_coeff;
}

for(j= 0;j<end;j++){
if(original_columns[j]==0){
mpz_set_si(soltest_u,0);
}else{
if(!iszeroone){
if(mpz_cmp_si(upperbounds[i],0)!=0){
mpz_divexact(soltest_upfac,upperbounds_max,upperbounds[i]);
}else{
mpz_set(soltest_upfac,upperbounds_max);
}
}
mpz_set_si(soltest_u,0);
for(k= 0;k<columns;k++){
if(ROUND(us[k])> 0){
mpz_addmul_ui(soltest_u,get_entry(k,i),ROUND(us[k]));
}else{
mpz_submul_ui(soltest_u,get_entry(k,i),-ROUND(us[k]));
}
}
mpz_sub(soltest_u,soltest_u,soltest_s);
mpz_divexact(soltest_u,soltest_u,max_norm_initial);
mpz_divexact(soltest_u,soltest_u,soltest_upfac);
mpz_divexact_ui(soltest_u,soltest_u,denom);
mpz_abs(soltest_u,soltest_u);
if(!iszeroone&&(mpz_cmp_si(soltest_u,0)<0||mpz_cmp(soltest_u,upperbounds[i])> 0)){
fprintf(stderr," rounding error -> this is not a solution!\n");
return 0;
}
i++;
}
mpz_out_str(NULL,10,soltest_u);
fflush(stdout);
mpz_out_str(fp,10,soltest_u);



printf(" ");
fprintf(fp," ");
}
if(free_RHS){
mpz_set_d(soltest_u,ROUND(w[i]));
mpz_divexact(soltest_u,soltest_u,max_up);
mpz_abs(soltest_u,soltest_u);
printf(" L = ");
mpz_out_str(NULL,10,soltest_u);
}
printf("\n");fflush(stdout);
fprintf(fp,"\n");fflush(fp);
}

nosolutions++;
#ifndef NO_OUTPUT
if(nosolutions%10000==0){
printf("%ld\n",nosolutions);fflush(stdout);
}
#endif

return 1;
}

/*:132*/
#line 3146 "./diophant.w"
;
/*131:*/
#line 3421 "./diophant.w"

DOUBLE Jacobi(DOUBLE**Ain,int n){
int i,j,k,nloops;
DOUBLE aa,si,co,tt,eps;
DOUBLE sum,ssum,amax,amin;
DOUBLE**V;
DOUBLE**A;

sum= 0.0;
eps= 0.000001;
nloops= 0;


V= (DOUBLE**)calloc(n,sizeof(DOUBLE*));
for(i= 0;i<n;++i)V[i]= (DOUBLE*)calloc(n,sizeof(DOUBLE));
A= (DOUBLE**)calloc(n,sizeof(DOUBLE*));
for(i= 0;i<n;++i)A[i]= (DOUBLE*)calloc(n,sizeof(DOUBLE));




for(i= 0;i<n;i++){
for(j= 0;j<n;j++){
V[i][j]= 0.0;
A[i][j]= 0.0;
for(k= 0;k<n;k++){
A[i][j]+= Ain[k][i]*Ain[k][j];
}

sum+= fabs(A[i][j]);
}
V[i][i]= 1.0;
}


if(n==1){
return A[0][0];
}

if(sum<=0.0){
return 0.0;
}

sum/= (n*n);


do{
ssum= 0.0;
amax= 0.0;
for(j= 1;j<n;j++){
for(i= 0;i<j;i++){

aa= fabs(A[i][j]);
if(aa> amax){
amax= aa;
}
ssum+= aa;
if(aa>=eps){

aa= atan2(2.0*A[i][j],A[i][i]-A[j][j])*0.5;
si= sin(aa);
co= cos(aa);

for(k= 0;k<n;k++){
tt= A[k][i];
A[k][i]= co*tt+si*A[k][j];
A[k][j]= -si*tt+co*A[k][j];
tt= V[k][i];
V[k][i]= co*tt+si*V[k][j];
V[k][j]= -si*tt+co*V[k][j];
}

A[i][i]= co*A[i][i]+si*A[j][i];
A[j][j]= -si*A[i][j]+co*A[j][j];
A[i][j]= 0.0;

for(k= 0;k<n;k++){
A[i][k]= A[k][i];
A[j][k]= A[k][j];
}

}
}
}
nloops++;
}
while(fabs(ssum)/sum> eps&&nloops<100000);
amin= A[0][0];
for(i= 1;i<n;i++)if(A[i][i]<amin)amin= A[i][i];

for(i= 0;i<n;i++)free(A[i]);
free(A);
for(i= 0;i<n;i++)free(V[i]);
free(V);

return amin;
}

/*:131*/
#line 3147 "./diophant.w"
;

/*:125*/
#line 1965 "./diophant.w"
;
/*133:*/
#line 3636 "./diophant.w"

void basis2poly(){
return;
}

/*:133*/
#line 1966 "./diophant.w"
;
/*134:*/
#line 3643 "./diophant.w"

void basis2LP(double*low,double*up){
return;
}

/*:134*/
#line 1967 "./diophant.w"
;

DOUBLE explicit_enumeration(COEFF**lattice,int columns,int rows)
{
/*89:*/
#line 2048 "./diophant.w"



int level,level_max;
int i,j,l;
long loops;

DOUBLE*y,*cs,*us;

long*delta,*d,*eta;
#if 0
mpz_t*v;
#else
long*v;
#endif
int*first_nonzero,*first_nonzero_in_column,*firstp;
int*snd_nonzero,*snd_nonzero_in_column,*sndp;
struct constraint*cons;

DOUBLE*N,**mu,*c,**w,**bd,**mu_trans;

DOUBLE Fd,Fq,Fqeps;
DOUBLE*dum;
DOUBLE tmp;
COEFF*swap_vec;
#if USE_SSE
int rowseven;
#endif

int isSideStep= 0;
DOUBLE stepWidth= 0.0;
DOUBLE olddum;

#if defined(FINCKEPOHST)
DOUBLE*fipo;
#endif

/*:89*//*104:*/
#line 2422 "./diophant.w"

#if EIGENBOUND
int k;
DOUBLE eig_s,eig_term2;
#endif

/*:104*/
#line 1971 "./diophant.w"
;


long nlow[1000];
for(i= 0;i<1000;i++)nlow[i]= 0;

/*90:*/
#line 2086 "./diophant.w"

#ifndef NO_OUTPUT
printf("Dimension of solution space (k): %d compared to s-z+2: %d\n",
columns,system_columns-system_rows+1+free_RHS);
fflush(stdout);
#endif

if(columns<system_columns-system_rows+1+free_RHS){
#ifndef NO_OUTPUT
fprintf(stderr,"LLL didn't succeed in computing a basis of the kernel.\n");
fprintf(stderr,"Please increase c0 (the first parameter)!\n");
printf("LLL didn't succeed in computing a basis of the kernel.\n");
printf("Please increase c0 (the first parameter)!\n");
#endif
return 0;
}


/*:90*/
#line 1977 "./diophant.w"
;
/*91:*/
#line 2106 "./diophant.w"

lllalloc(&mu,&c,&N,&bd,columns,rows);

us= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
cs= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
y= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
delta= (long*)calloc(columns+1,sizeof(long));
d= (long*)calloc(columns+1,sizeof(long));
first_nonzero= (int*)calloc(rows,sizeof(int));
first_nonzero_in_column= (int*)calloc(columns+rows+1,sizeof(int));
if(first_nonzero_in_column==NULL)return(0);
firstp= (int*)calloc(columns+1,sizeof(int));

snd_nonzero= (int*)calloc(rows,sizeof(int));
snd_nonzero_in_column= (int*)calloc(columns+rows+1,sizeof(int));
if(snd_nonzero_in_column==NULL)return(0);
sndp= (int*)calloc(columns+1,sizeof(int));

cons= (struct constraint*)calloc(columns,sizeof(struct constraint));
for(i= 0;i<columns;++i){
cons[i].isSet= 0;
}

eta= (long*)calloc(columns+1,sizeof(long));
v= (long*)calloc(columns+1,sizeof(long));
w= (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
#if USE_SSE
rowseven= (rows%8!=0)?(rows/8+1)*8:rows;
for(i= 0;i<=columns;i++)w[i]= (DOUBLE*)calloc(rowseven,sizeof(DOUBLE));
#else
for(i= 0;i<=columns;i++)w[i]= (DOUBLE*)calloc(rows,sizeof(DOUBLE));
#endif

mu_trans= (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
for(i= 0;i<=columns;i++)mu_trans[i]= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
dum= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));

#if FINCKEPOHST
fipo= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
muinv= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)muinv[i]= (DOUBLE*)calloc(rows,sizeof(DOUBLE));

fipo_LB= (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
fipo_UB= (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
for(i= 0;i<=columns;++i){
fipo_LB[i]= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
fipo_UB[i]= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
}
#endif

/*:91*/
#line 1978 "./diophant.w"
;
/*93:*/
#line 2166 "./diophant.w"

#if EIGENBOUND
Rinv= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)Rinv[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));

R= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)R[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));

eig_f= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)eig_f[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));

eig_RinvR= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)eig_RinvR[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));

eig_min= (DOUBLE*)calloc(columns,sizeof(DOUBLE));
eig_bound= (DOUBLE*)calloc(columns,sizeof(DOUBLE));

Rinvi= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)Rinvi[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));
bnd_up= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)bnd_up[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));
bnd_lo= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)bnd_lo[i]= (DOUBLE*)calloc(columns,sizeof(DOUBLE));
#endif

/*:93*/
#line 1979 "./diophant.w"
;
/*95:*/
#line 2209 "./diophant.w"

for(i= 0;i<=columns;i++){
cs[i]= y[i]= us[i]= 0.0;
delta[i]= 0;
#if 0
mpz_set_si(v[i],0);
#else
v[i]= 0;
#endif
eta[i]= d[i]= 1;
for(l= 0;l<rows;l++)w[i][l]= 0.0;
}

/*:95*/
#line 1980 "./diophant.w"
;
/*98:*/
#line 2250 "./diophant.w"

if(free_RHS){
i= 0;
for(j= columns-1;j>=0;j--)if(mpz_sgn(get_entry(j,rows-2))!=0)i++;
#ifndef NO_OUTPUT
printf("Number of nonzero entries in the second last row: %d\n",i);
fflush(stdout);
#endif
}

i= 0;
for(j= columns-1;j>=0;j--)if(mpz_sgn(get_entry(j,rows-1))!=0)i++;
#ifndef NO_OUTPUT
printf("Number of nonzero entries in the last row: %d\n",i);
fflush(stdout);
#endif

/*:98*/
#line 1981 "./diophant.w"
;
#if 0
/*96:*/
#line 2225 "./diophant.w"

for(j= columns-1;j> 0;j--){
for(l= j-1;l>=0;l--){
if(mpz_cmpabs(get_entry(l,rows-1),get_entry(j,rows-1))> 0){
swap_vec= lattice[l];
for(i= l+1;i<=j;i++)lattice[i-1]= lattice[i];
lattice[j]= swap_vec;
}
}
}

/*:96*/
#line 1983 "./diophant.w"
;
#endif

/*99:*/
#line 2272 "./diophant.w"

Fq= (DOUBLE)mpz_get_d(max_norm);
Fd= (rows*Fq*Fq)*(1.0+EPSILON);
Fqeps= (1.0+EPSILON)*Fq;
#ifndef NO_OUTPUT
#if VERBOSE >  0
printf("Fq: %f\n",(double)Fq);
printf("Fd: %f\n",(double)Fd);
#endif
#endif

/*:99*/
#line 1986 "./diophant.w"
;
/*101:*/
#line 2286 "./diophant.w"

#if GIVENS
givens(lattice,columns,rows,mu,bd,c,N,Fq);
#else
gramschmidt(lattice,columns,rows,mu,bd,c,N,Fq);
#endif


for(i= 0;i<columns;i++)
for(j= 0;j<columns;j++)
mu_trans[j][i]= mu[i][j];



/*:101*/
#line 1987 "./diophant.w"
;

#if 0
basis2poly();
#endif

#if FINCKEPOHST
/*102:*/
#line 2326 "./diophant.w"

fipo_success= 0;
inverse(mu,muinv,columns);

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("\nFincke-Pohst bounds:\n");fflush(stdout);
#endif
#endif

for(i= 0;i<columns;i++){
fipo[i]= 0.0;
dum1= 0.0;
for(j= 0;j<rows;j++){
tmp= 0.0;
for(l= i;l<columns;l++)tmp+= muinv[i][l]*bd[l][j]/c[l];
fipo[i]+= tmp*tmp;
dum1+= fabs(tmp);
}
fipo[i]= SQRT(fipo[i]*Fd);
dum1= fabs(dum1*Fq);
if(dum1<fipo[i])fipo[i]= dum1;
fipo[i]*= (1.0+EPSILON);

fipo_LB[columns][i]= -fipo[i];
fipo_UB[columns][i]= fipo[i];
#ifndef NO_OUTPUT

#endif

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("%0.3lf ",fipo[i]);
#endif
#endif
}

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("\n\n");fflush(stdout);
#endif
#endif

/*:102*/
#line 1994 "./diophant.w"
;
#endif






#if 1
for(i= columns-1;i>=0;i--){
if(fipo[i]<0.9){
printf("DEL\n");
columns--;
}else{
break;
}
}
#endif

#if 0
/*101:*/
#line 2286 "./diophant.w"

#if GIVENS
givens(lattice,columns,rows,mu,bd,c,N,Fq);
#else
gramschmidt(lattice,columns,rows,mu,bd,c,N,Fq);
#endif


for(i= 0;i<columns;i++)
for(j= 0;j<columns;j++)
mu_trans[j][i]= mu[i][j];



/*:101*/
#line 2014 "./diophant.w"
;
/*102:*/
#line 2326 "./diophant.w"

fipo_success= 0;
inverse(mu,muinv,columns);

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("\nFincke-Pohst bounds:\n");fflush(stdout);
#endif
#endif

for(i= 0;i<columns;i++){
fipo[i]= 0.0;
dum1= 0.0;
for(j= 0;j<rows;j++){
tmp= 0.0;
for(l= i;l<columns;l++)tmp+= muinv[i][l]*bd[l][j]/c[l];
fipo[i]+= tmp*tmp;
dum1+= fabs(tmp);
}
fipo[i]= SQRT(fipo[i]*Fd);
dum1= fabs(dum1*Fq);
if(dum1<fipo[i])fipo[i]= dum1;
fipo[i]*= (1.0+EPSILON);

fipo_LB[columns][i]= -fipo[i];
fipo_UB[columns][i]= fipo[i];
#ifndef NO_OUTPUT

#endif

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("%0.3lf ",fipo[i]);
#endif
#endif
}

#ifndef NO_OUTPUT
#if VERBOSE >  -1
printf("\n\n");fflush(stdout);
#endif
#endif

/*:102*/
#line 2015 "./diophant.w"
;
#endif

#if EIGENBOUND
/*103:*/
#line 2370 "./diophant.w"

#if EIGENBOUND
for(i= 0;i<columns;i++){
for(j= 0;j<columns;j++){
R[i][j]= mu[j][i]*SQRT(c[i]);
Rinv[i][j]= muinv[i][j]/SQRT(c[j]);
}
}
#if 0
printf("\nCheck\n");
for(i= 0;i<columns;i++){
for(j= 0;j<columns;j++){
for(l= 0,eig_s= 0.0;l<columns;l++)eig_s+= Rinv[i][l]*R[l][j];
printf("%0.3lf\t",eig_s);
}
printf("\n");
}
printf("\n");
#endif

for(k= 0;k<columns;k++){
for(i= 0;i<k;i++){
eig_s= 0.0;
for(j= 0;j<k;j++){
eig_s+= Rinv[i][j]*R[j][k];
}
eig_RinvR[k][i]= eig_s;
}
}

for(k= 0;k<columns;k++){
for(i= 0;i<k;i++){
Rinvi[k][i]= 0.0;
for(j= i;j<k;j++){
Rinvi[k][i]+= Rinv[i][j]*Rinv[i][j];
}
}
}

#if 0
printf("Jacobi: \n");
for(k= 0;k<columns;k++){
eig_s= Jacobi(R,k);
eig_min[k]= eig_s;
printf("%0.3lf ",eig_s);fflush(stdout);
}
printf("\n");
printf("End Jacobi \n");
eig_cut= 0;
#endif
#endif

/*:103*/
#line 2019 "./diophant.w"
;
#endif

#if 0
basis2LP(fipo_l,fipo_u);
#endif

/*106:*/
#line 2504 "./diophant.w"

for(l= 0;l<rows;l++){
for(i= 0;i<columns;i++)if(mpz_sgn(get_entry(i,l))!=0){
first_nonzero[l]= i;
break;
}
}

printf("First non-zero entries:\n");
j= 0;
for(l= 0;l<columns;l++){
firstp[l]= j;
first_nonzero_in_column[j]= 0;
j++;
for(i= 0;i<rows;i++){
if(first_nonzero[i]==l){
first_nonzero_in_column[j]= i;
first_nonzero_in_column[firstp[l]]++;
j++;
}
}
printf("%d ",first_nonzero_in_column[firstp[l]]);
}
printf(": %d\n",rows);
firstp[columns]= j;
first_nonzero_in_column[j]= 0;

/*:106*/
#line 2026 "./diophant.w"
;
/*107:*/
#line 2535 "./diophant.w"

for(l= 0;l<rows;l++){
for(i= first_nonzero[l]+1;i<columns;i++)if(mpz_sgn(get_entry(i,l))!=0){
snd_nonzero[l]= i;
break;
}
}

printf("Second non-zero entries:\n");
j= 0;
for(l= 0;l<columns;l++){
sndp[l]= j;
snd_nonzero_in_column[j]= 0;
j++;
for(i= 0;i<rows;i++){
if(snd_nonzero[i]==l){
snd_nonzero_in_column[j]= i;
snd_nonzero_in_column[sndp[l]]++;
j++;
}
}
printf("%d ",snd_nonzero_in_column[sndp[l]]);
}
printf(": %d\n\n",rows);
sndp[columns]= j;
snd_nonzero_in_column[j]= 0;


#if 0

for(l= 0;l<rows;l++){
printf("%d ",snd_nonzero[l]-first_nonzero[l]);
}
printf("\n");
#endif

/*:107*/
#line 2027 "./diophant.w"
;
/*109:*/
#line 2579 "./diophant.w"

level= first_nonzero[rows-1];
if(level<0)level= 0;
level_max= level;
us[level]= 1;
#if 0
mpz_set_si(v[level],1);
#else
v[level]= 1;
#endif
only_zeros_no= only_zeros_success= 0;
hoelder_no= hoelder_success= 0;
cs_success= nosolutions= loops= 0;
N_success= 0;
N2_success= 0;
N3_success= 0;

/*:109*/
#line 2028 "./diophant.w"
;
/*111:*/
#line 2634 "./diophant.w"

do{
/*110:*/
#line 2598 "./diophant.w"

loops++;
if((stop_after_loops> 0)&&(stop_after_loops<=loops))goto afterloop;

#ifndef NO_OUTPUT
#if VERBOSE >  -1
if(loops%100000000==0){
printf("%ld loops, solutions: %ld ",loops,nosolutions);
#if FINCKEPOHST
printf("fipo: %ld ",fipo_success);
#endif
#if EIGENBOUND
printf("eig_bound: %ld ",eig_cut);
#endif
printf("pruneN: %ld ",N2_success);
printf("pruneN2: %ld ",N3_success);
printf("\n");
#if 0

for(i= 0;i<level_max;i++){
printf("%03d: %ld\n",i,nlow[i]);
}
#endif







fflush(stdout);
}
#endif
#endif

/*:110*/
#line 2636 "./diophant.w"
;
/*112:*/
#line 2699 "./diophant.w"

olddum= dum[level];
dum[level]= us[level]+y[level];
cs[level]= cs[level+1]+dum[level]*dum[level]*c[level];

if(isSideStep){
stepWidth= dum[level]-olddum;
}

/*:112*/
#line 2637 "./diophant.w"
;
if((cs[level]<Fd)){
#if FINCKEPOHST
#if 1
if(fabs(us[level])> fipo[level]){
#else
if(level!=columns-1&&
(us[level]> fipo_UB[columns][level]||us[level]<fipo_LB[columns][level])
){
#endif
fipo_success++;
goto side_step;
}

#endif
#if EIGENBOUND
if(level<columns-1&&(us[level]<bnd_lo[level][level]||us[level]> bnd_up[level][level])){
printf("%d:\t%0.3lf %0.3lf %0.3lf -> cut\n",level,bnd_lo[level][level],us[level],bnd_up[level][level]);
goto side_step;
}
#endif
if(isSideStep){
compute_w2(w,bd,stepWidth,level,rows);
}else{
compute_w(w,bd,dum[level],level,rows);
}

if(level> 0){
/*113:*/
#line 2711 "./diophant.w"

i= prune_only_zeros(w,level,rows,Fq,first_nonzero_in_column,firstp,
bd,y,us,columns);

if(i<0){
goto step_back;
}else if(i> 0){
goto side_step;
}

i= prune_snd_nonzero(columns,rows,
level,Fq,
first_nonzero,
snd_nonzero_in_column,sndp,
us,
cons);
if(i> 0){
goto side_step;
}

#if 0
printf("%d ",level);
#endif
if(prune(w[level],cs[level],rows,Fqeps)){
if(eta[level]==1){
goto step_back;
}
eta[level]= 1;
delta[level]*= -1;
if(delta[level]*d[level]>=0)delta[level]+= d[level];
#if 0
us[level]= mpz_get_d(v[level])+delta[level];
#else
us[level]= v[level]+delta[level];
#endif
}else{
#if 0
if(pruneN(w,cs,level,rows,columns,Fq)){
goto side_step;
}
#endif

#if EIGENBOUND
if(2*level> 0.0*columns){
/*105:*/
#line 2428 "./diophant.w"

#if 1
if(level==columns-1){
printf("Oben\n");
eig_s= 0.0;
for(i= 0;i<level;i++){
eig_f[level][i]= eig_RinvR[level][i]*us[level];
dum1= fabs(eig_f[level][i]);
dum1-= round(dum1);
eig_s+= dum1*dum1;
}
eig_bound[level]= eig_s*eig_min[level];

for(i= 0;i<level;i++){
bnd_up[level-1][i]= fipo[i];
bnd_lo[level-1][i]= -fipo[i];
}
}else{
eig_term2= 0.0;
for(k= level+1;k<columns;k++){
eig_term2+= R[level][k]*us[k];
}

eig_s= 0.0;
for(i= 0;i<level;i++){
eig_f[level][i]= eig_f[level+1][i]-Rinv[i][level]*eig_term2+eig_RinvR[level][i]*us[level];
}

for(i= 0;i<level;i++){
dum1= SQRT(Rinvi[level][i]*(Fd-cs[level]));

bnd_up[level-1][i]= bnd_up[level][i];
if(bnd_up[level-1][i]> dum1-eig_f[level][i])bnd_up[level-1][i]= dum1-eig_f[level][i];

bnd_lo[level-1][i]= bnd_lo[level][i];
if(bnd_lo[level-1][i]<-dum1-eig_f[level][i])bnd_lo[level-1][i]= -dum1-eig_f[level][i];


if(bnd_lo[level-1][i]> bnd_up[level-1][i]){
printf("CUT\n");
}
}
}
#else
if(level==0){
eig_bound[level]= 0.0;
}else{
for(i= 0;i<level;i++){
eig_s= 0.0;
for(k= level;k<columns;k++){
eig_s+= R[i][k]*us[k];
}
eig_f[0][i]= eig_s;
}
eig_s= 0.0;
for(i= 0;i<level;i++){
eig_f[level][i]= 0.0;
for(k= 0;k<level;k++){
eig_f[level][i]+= Rinv[i][k]*eig_f[0][k];
}
dum1= fabs(eig_f[level][i]);
dum1-= round(dum1);
eig_s+= dum1*dum1;
}
if(fabs(eig_s*eig_min[level]-eig_bound[level])> 0.000001){
printf("WRONG: %0.5lf %0.5lf\n",eig_s*eig_min[level],eig_bound[level]);
}
eig_bound[level]= eig_s*eig_min[level];
}
#endif


/*:105*/
#line 2755 "./diophant.w"
;
#if 0
if(cs[level]+eig_bound[level]> Fd){
printf("%d:\t%0.3lf %0.3lf %0.3lf %0.1lf -> cut\n",level,cs[level],Fd,Fd-eig_bound[level],us[level]);
}else{
printf("%d:\t%0.3lf %0.3lf %0.3lf %0.1lf \n",level,cs[level],Fd,Fd-eig_bound[level],us[level]);
}
if(cs[level]+eig_bound[level]> Fd){
eig_cut++;
goto side_step;
}
#endif
}
#endif

level--;
#if 0
fipo_LB[columns][level]= -fipo[level];
fipo_UB[columns][level]= fipo[level];
#endif
delta[level]= eta[level]= 0;
y[level]= compute_y(mu_trans,us,level,level_max);
#if 0
mpz_set_d(v[level],ROUND(-dum));
us[level]= mpz_get_si(v[level]);
#else
us[level]= v[level]= ROUND(-y[level]);
#endif
d[level]= (v[level]> -y[level])?-1:1;
isSideStep= 0;
}

/*:113*/
#line 2665 "./diophant.w"
;
}else{
/*114:*/
#line 2788 "./diophant.w"

if(exacttest(w[0],rows,Fq)==1){
print_solution(w[level],rows,Fq,us,columns);
#if 0
printf("us: ");
for(i= columns-1;i>=0;i--){
printf("%d ",(int)us[i]);
}
printf("after %ld loops\n",loops);
#endif
if((stop_after_solutions> 0)&&(stop_after_solutions<=nosolutions))
goto afterloop;
}
goto side_step;

/*:114*/
#line 2667 "./diophant.w"
;
}
}else{
cs_success++;
step_back:

nlow[level]++;
level++;
if(level_max<level)level_max= level;
side_step:




if(eta[level]==0){
delta[level]*= -1;
if(delta[level]*d[level]>=0)delta[level]+= d[level];
}else{
delta[level]+= (delta[level]*d[level]>=0)?d[level]:-d[level];
}
#if 0
us[level]= mpz_get_d(v[level])+delta[level];
#else
us[level]= v[level]+delta[level];
#endif
isSideStep= 1;
}
}while(level<columns);
afterloop:

/*:111*/
#line 2029 "./diophant.w"
;

/*115:*/
#line 2804 "./diophant.w"

#ifndef NO_OUTPUT
printf("Prune_cs: %ld\n",cs_success);
printf("Prune_only_zeros: %ld of %ld\n",only_zeros_success,only_zeros_no);
printf("Prune_hoelder: %ld of %ld\n",hoelder_success,hoelder_no);
printf("Prune_N: %ld\n",N_success);
printf("Prune_N2: %ld\n",N2_success);
printf("Prune_N3: %ld\n",N3_success);
#if FINCKEPOHST
printf("Fincke-Pohst: %ld\n",fipo_success);
#endif
#if EIGENBOUND
printf("Eigen bound: %ld\n",eig_cut);
#endif
printf("Loops: %ld\n",loops);
#endif
if((stop_after_solutions<=nosolutions&&stop_after_solutions> 0)||
(stop_after_loops<=loops&&stop_after_loops> 0)){
#ifndef NO_OUTPUT
printf("Stopped after number of solutions: %ld\n",nosolutions);
#endif
/*118:*/
#line 2896 "./diophant.w"

if(SILENT)fprintf(fp,"%ld solutions\n",nosolutions);
fflush(fp);

/*:118*/
#line 2825 "./diophant.w"
;
if((stop_after_loops<=loops&&stop_after_loops> 0)){
exit(10);
}else{
exit(9);
}
}else{
#ifndef NO_OUTPUT
printf("Total number of solutions: %ld\n",nosolutions);
#endif
}
#ifndef NO_OUTPUT
printf("\n");fflush(stdout);
#endif

/*:115*/
#line 2031 "./diophant.w"
;
/*116:*/
#line 2840 "./diophant.w"

free(us);
free(cs);
free(y);
free(delta);
free(d);
free(first_nonzero);
free(first_nonzero_in_column);
free(firstp);

free(snd_nonzero);
free(snd_nonzero_in_column);
free(sndp);
free(cons);

free(eta);
free(v);
for(l= 0;l<=columns;l++)free(w[l]);
free(w);
free(original_columns);

#if FINCKEPOHST
free(fipo);
for(l= 0;l<columns;l++)free(muinv[l]);
free(muinv);
#endif

#if EIGENBOUND
for(l= 0;l<columns;l++)free(Rinv[l]);
free(Rinv);
for(l= 0;l<columns;l++)free(R[l]);
free(R);
for(l= 0;l<columns;l++)free(eig_f[l]);
free(eig_f);
for(l= 0;l<columns;l++)free(eig_RinvR[l]);
free(eig_RinvR);
free(eig_min);
free(eig_bound);

for(l= 0;l<columns;l++)free(Rinvi[l]);
free(Rinvi);
for(l= 0;l<columns;l++)free(bnd_up[l]);
free(bnd_up);
for(l= 0;l<columns;l++)free(bnd_lo[l]);
free(bnd_lo);
#endif

lllfree(mu,c,N,bd,columns);
for(l= 0;l<columns;l++)free(mu_trans[l]);
free(mu_trans);

/*:116*/
#line 2032 "./diophant.w"
;
return 1;
}

/*:86*/
#line 948 "./diophant.w"
;

/*:47*/
#line 51 "./diophant.w"
;

/*:2*//*3:*/
#line 54 "./diophant.w"

long diophant(mpz_t**a_input,mpz_t*b_input,mpz_t*upperbounds_input,
int no_columns,int no_rows,
mpz_t factor_input,mpz_t norm_input,mpz_t scalelastlinefactor,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
long stop_after_sol_input,long stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input,
int cut_after,int nboundedvars,FILE*solfile)
{
int i,j;
DOUBLE lD,lDnew;
COEFF*swap_vec;

/*11:*/
#line 250 "./diophant.w"

mpz_init_set(matrix_factor,factor_input);
mpz_init_set(max_norm,norm_input);
mpz_init(lastlines_factor);
mpz_init(upfac);

mpz_init(snd_q);
mpz_init(snd_r);
mpz_init(snd_s);

if(iterate){
no_iterates= iterate_no;
}else{
bkz_beta= bkz_beta_input;
bkz_p= bkz_p_input;
}
SILENT= silent;
stop_after_solutions= stop_after_sol_input;
stop_after_loops= stop_after_loops_input;
free_RHS= free_RHS_input;
nom= 1;
denom= 2;

system_rows= no_rows;
system_columns= no_columns;
nboundvars= nboundedvars;

/*:11*//*42:*/
#line 842 "./diophant.w"

mpz_init(soltest_u);
mpz_init(soltest_s);
mpz_init_set_ui(soltest_upfac,1);

/*:42*/
#line 68 "./diophant.w"
;
#if BLAS

#endif





/*13:*/
#line 303 "./diophant.w"

lattice_rows= system_rows+system_columns+1;
lattice_columns= system_columns+2;

if(free_RHS){
lattice_rows++;
lattice_columns++;
}else{
#ifndef NO_OUTPUT
fprintf(stderr,"The RHS is fixed !\n");fflush(stderr);
#endif
}
cut_after_coeff= cut_after;

/*:13*/
#line 77 "./diophant.w"
;
/*14:*/
#line 318 "./diophant.w"

lattice= (COEFF**)calloc(lattice_columns,sizeof(COEFF*));
for(j= 0;j<lattice_columns;j++){
lattice[j]= (COEFF*)calloc(lattice_rows+1,sizeof(COEFF));
for(i= 0;i<=lattice_rows;i++)mpz_init(lattice[j][i].c);
}
/*:14*/
#line 78 "./diophant.w"
;
/*15:*/
#line 326 "./diophant.w"

for(j= 0;j<system_rows;j++){
for(i= 0;i<system_columns;i++){
mpz_mul(lattice[i][j+1].c,a_input[j][i],matrix_factor);
}
mpz_mul(lattice[system_columns][j+1].c,b_input[j],matrix_factor);
}
/*:15*/
#line 79 "./diophant.w"
;
/*16:*/
#line 346 "./diophant.w"

mpz_init_set_si(upperbounds_max,1);
iszeroone= 1;
if(upperbounds_input==NULL){
#ifndef NO_OUTPUT
printf("No upper bounds: 0/1 variables are assumed \n");fflush(stdout);
#endif
}else{
upperbounds= (mpz_t*)calloc(system_columns,sizeof(mpz_t));
for(i= 0;i<system_columns;i++)mpz_init_set_si(upperbounds[i],1);
for(i= 0;i<nboundvars;i++){
mpz_set(upperbounds[i],upperbounds_input[i]);
if(mpz_sgn(upperbounds[i])!=0){
mpz_lcm(upperbounds_max,upperbounds_max,upperbounds[i]);
}
}
if(mpz_cmp_si(upperbounds_max,1)> 0)iszeroone= 0;
#ifndef NO_OUTPUT
fprintf(stderr,"upper bounds found. Max=");
mpz_out_str(stderr,10,upperbounds_max);
fprintf(stderr,"\n");
fflush(stderr);
#endif
}
/*:16*/
#line 80 "./diophant.w"
;
/*17:*/
#line 374 "./diophant.w"

if(org_col_input!=NULL)no_original_columns= no_org_col_input;
else no_original_columns= system_columns;

original_columns= (int*)calloc(no_original_columns,sizeof(int));

if(org_col_input!=NULL)
for(i= 0;i<no_original_columns;i++)original_columns[i]= org_col_input[i];
else{
for(i= 0;i<no_original_columns;i++)original_columns[i]= 1;
#ifndef NO_OUTPUT
printf("No preselected columns \n");fflush(stdout);
#endif
}

/*:17*/
#line 81 "./diophant.w"
;
/*18:*/
#line 397 "./diophant.w"

for(j= system_rows;j<lattice_rows;j++){
mpz_mul_si(lattice[j-system_rows][j+1].c,max_norm,denom);
mpz_mul_si(lattice[lattice_columns-2][j+1].c,max_norm,nom);
}
mpz_set(lattice[system_columns+free_RHS][lattice_rows].c,max_norm);

if(free_RHS){
mpz_set_si(lattice[system_columns][lattice_rows-1].c,1);
mpz_set_si(lattice[system_columns+1][lattice_rows-1].c,0);
}
mpz_set(lattice[system_columns+free_RHS][lattice_rows].c,max_norm);
for(i= 0;i<lattice_columns-1;i++)coeffinit(lattice[i],lattice_rows);

/*:18*/
#line 82 "./diophant.w"
;
/*117:*/
#line 2891 "./diophant.w"

fp= solfile;
if(SILENT)fprintf(fp,"SILENT\n");
fflush(fp);

/*:117*/
#line 83 "./diophant.w"
;

#if 0
printf("Before scaling\n");
print_lattice();
#endif
/*19:*/
#line 419 "./diophant.w"

mpz_init_set(max_norm_initial,max_norm);
mpz_init_set_si(max_up,1);
if(!iszeroone){
for(j= 0;j<nboundvars;j++){
if(mpz_sgn(upperbounds[j])!=0){
mpz_divexact(upfac,upperbounds_max,upperbounds[j]);
}else{
mpz_mul(upfac,upperbounds_max,upperbounds_max);
mpz_mul_si(upfac,upfac,10000);
}
smult_lattice(j,j+system_rows,upfac);
smult_lattice(system_columns+free_RHS,j+system_rows,upperbounds_max);
}
mpz_set(max_up,upperbounds_max);
mpz_mul(max_norm,max_norm,max_up);
if(free_RHS)
smult_lattice(system_columns,lattice_rows-2,max_up);

smult_lattice(system_columns+free_RHS,lattice_rows-1,max_up);
}

/*:19*/
#line 89 "./diophant.w"
;
#if 0
printf("After scaling\n");
print_lattice();
#endif

#if 1
#if 0
print_NTL_lattice();
return 0;
#endif

/*20:*/
#line 444 "./diophant.w"

swap_vec= lattice[lattice_columns-2];
for(i= lattice_columns-2;i> 0;i--)lattice[i]= lattice[i-1];
lattice[0]= swap_vec;

/*:20*/
#line 101 "./diophant.w"
;
#if 0
printf("After permute\n");
print_lattice();
#endif
shufflelattice();
/*21:*/
#line 451 "./diophant.w"

mpz_set_ui(lastlines_factor,1);
#ifndef NO_OUTPUT
printf("\n");fflush(stdout);
#endif
lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_LOW);

/*:21*/
#line 107 "./diophant.w"
;
#if 0
printf("After first reduction\n");
print_lattice();
#endif
/*23:*/
#line 475 "./diophant.w"

if(cutlattice()){
#ifndef NO_OUTPUT
printf("First reduction successful\n");fflush(stdout);
#endif
}else{
#ifndef NO_OUTPUT
printf("First reduction not successful\n");fflush(stdout);
#endif
return 0;
}

for(j= 0;j<lattice_columns-1;j++)solutiontest(j);
/*:23*/
#line 112 "./diophant.w"
;
#if 0
printf("After cutting\n");
print_lattice();
#endif

#if 1
shufflelattice();
/*22:*/
#line 460 "./diophant.w"

mpz_set_ui(lastlines_factor,1);
lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_HIGH);
#ifndef NO_OUTPUT
printf("Second reduction successful\n");fflush(stdout);
#endif

/*:22*/
#line 120 "./diophant.w"
;
#endif

#if 0
printf("After second reduction\n");
print_lattice();
#endif

#if 1
/*24:*/
#line 496 "./diophant.w"


mpz_set(lastlines_factor,scalelastlinefactor);
for(i= 0;i<lattice_columns;i++)
mpz_mul(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
if(free_RHS)
for(i= 0;i<lattice_columns;i++)
mpz_mul(lattice[i][lattice_rows-1].c,lattice[i][lattice_rows-1].c,lastlines_factor);

#if 0
for(i= 0;i<lattice_columns;i++){
for(j= 0;j<40;j++)
mpz_mul_ui(lattice[i][j+1].c,lattice[i][j+1].c,9);
}
#endif

/*:24*/
#line 129 "./diophant.w"
;
/*26:*/
#line 531 "./diophant.w"

#ifndef NO_OUTPUT
printf("\n");fflush(stdout);
#endif
if(iterate){
iteratedlll(lattice,lattice_columns-1,lattice_rows,no_iterates,LLLCONST_HIGH);
}else{
shufflelattice();
lDnew= bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGHER,40,bkz_p);

i= 0;
do{
lD= lDnew;
shufflelattice();
lDnew= bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGH,bkz_beta,bkz_p);
printf("%0.3lf %0.3lf %0.3lf\n",lD,lDnew,lD-lDnew);
i++;
}
while(i<1&&fabs(lDnew-lD)> 0.01);
}
#ifndef NO_OUTPUT
printf("Third reduction successful\n");fflush(stdout);
#endif

/*:26*/
#line 130 "./diophant.w"
;
/*25:*/
#line 514 "./diophant.w"

for(i= 0;i<lattice_columns;i++)
mpz_divexact(lattice[i][lattice_rows].c,lattice[i][lattice_rows].c,lastlines_factor);
if(free_RHS)
for(i= 0;i<lattice_columns;i++)
mpz_divexact(lattice[i][lattice_rows-1].c,
lattice[i][lattice_rows-1].c,lastlines_factor);

#if 0
for(i= 0;i<lattice_columns;i++){
for(j= 0;j<40;j++)
mpz_divexact_ui(lattice[i][j+1].c,lattice[i][j+1].c,9);
}
#endif

/*:25*/
#line 131 "./diophant.w"
;
#endif

#else
read_NTL_lattice();
#endif

#if 0
printf("Before enumeration\n");

print_lattice();
#endif

/*27:*/
#line 556 "./diophant.w"

#ifndef NO_OUTPUT
printf("\n");fflush(stdout);
#endif
nosolutions= explicit_enumeration(lattice,lattice_columns-1,lattice_rows);

/*:27*/
#line 144 "./diophant.w"
;
/*118:*/
#line 2896 "./diophant.w"

if(SILENT)fprintf(fp,"%ld solutions\n",nosolutions);
fflush(fp);

/*:118*/
#line 145 "./diophant.w"
;
/*12:*/
#line 277 "./diophant.w"

mpz_clear(matrix_factor);
mpz_clear(max_norm);
mpz_clear(lastlines_factor);
mpz_clear(upfac);
mpz_clear(max_norm_initial);
mpz_clear(max_up);
mpz_clear(soltest_u);
mpz_clear(soltest_s);
mpz_clear(soltest_upfac);
mpz_clear(upperbounds_max);

for(j= 0;j<lattice_columns;j++){
for(i= 0;i<=lattice_rows;i++)mpz_clear(lattice[j][i].c);
}
free(lattice);
if(upperbounds!=NULL){
for(i= 0;i<system_columns;i++)mpz_clear(upperbounds[i]);
free(upperbounds);
}

/*:12*/
#line 146 "./diophant.w"
;
return nosolutions;
}

/*:3*/
