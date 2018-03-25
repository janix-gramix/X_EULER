#include<stdio.h>
#include<string.h>
#define maxsum 3000
#define maxsize 8

#define startsum (1+2*1)
#define startsize 2

//#define startsum (51+6*11)
//#define startsize 6

//#define startsum (115+7*20)
//#define startsize 7

int min,size,sumsize;

long long int count=0;
char hit[maxsum];
int v[maxsize-1],vmin[maxsize-1];

char check(int *v, char n) {
    if(n<3) return 1;
    char ok=1; int sum=0;
    
    void test(int min, int todo){
        int x;todo--;
        for(x=min;ok&&x<n-todo;x++){
            sum+=v[x];
            if(todo)test(x+1,todo);
            if(!todo)if(hit[sum]++)ok=0;
            sum-=v[x];
        }
    }
    
    int sumsize=(n+1)/2,i;
    for(i=0;i<sumsize;i++)sum+=v[n-1-i];
    memset(hit,0,sum+1);
    
    sum=0;
    test(0,sumsize);
    if(ok)test(0,sumsize-1);
    return ok;
}

void loop(int *v, int i){
    int a,max,j,k,n,s,ds,dn;
    
    v[i]=a=(i?v[i-1]+1:1);
    
    for(j=i+1;j<size-1;j++)v[j]=v[j-1]+1;
    j=0;k=size-2;n=1;
    while(j<k) n+=v[k--]-v[j++];
    for(j=0,s=size*n;j<size-1;j++) s+=v[j];
    
    if(i<sumsize){dn=i;ds=(size+1)*sumsize-(size-1)*(sumsize-i);}
    else {dn+=size-i-1;ds=(size+1)*(size-1-i);}
    if(size%2==0)
        if(i<sumsize) ds++;
        else if(i==sumsize) {dn--;ds-=(size+1)-1;}
    
    for(;s<=min;v[i]=++a,s+=ds,n+=dn)if(check(v,i+1)){
        if(i==size-2) {
            for(j=0;j<size;j++)printf("%d ",vmin[j]=n+(j?v[j-1]:0));
            printf("(sum=%d)\n",s);
            min=s;
        }
        else
            loop(v,i+1);
    }
}

int main(){
    int i,j,k,n;
    for(i=0;i<1;i++){
        min=startsum;
        size=startsize;
        while(size<=maxsize&&min<maxsum){
            sumsize=(size-1)/2;
            printf("--\nsize=%d limit=%d\n",size,min);
            loop(v,0);
            for(j=0,k=++size-2,n=1;j<k;)n+=vmin[k--]-vmin[j++];
            min+=size*n;
        }
    }
}
