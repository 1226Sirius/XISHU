#include<stdio.h>
#include<stdlib.h>
#define MAXSIZE 20
#define MAXRC 20
typedef struct
{                //稀疏矩阵的三元组顺序表存储表示
int i,j;                  //该非零元的行下标和列下标
int e;
}Triple;

typedef struct
{
Triple data[MAXSIZE+1];      //非零元三元组表，data[0]未用
int rpos[MAXRC+1];        //各行第一个非零元的位置表
int mu,nu,tu;            //矩阵的行数列数和非零元的个数
}RLSMatrix;

void CreateSMatrix(RLSMatrix *T)          //输入创建稀疏矩阵
{
int k;
printf(" \n 请输入矩阵行数、列数及非零元个数: ");
scanf("%d%d%d",&T->mu,&T->nu,&T->tu);
printf("\n");
if(T->tu>MAXSIZE||T->mu>21)
{
printf(" 非零个数超出定义范围！出错！");
exit(0);
}
for(k=1;k<=T->tu;k++)
{
printf(" 请输入第%d个非零元素的行数,列数及其值: ",k);
scanf("%d%d%d",&T->data[k].i,&T->data[k].j,&T->data[k].e);
}
}

void AddRLSMatrix(RLSMatrix M,RLSMatrix N,RLSMatrix *Q) //稀疏矩阵相加
{
int p,q,k=1;
if(M.mu!=N.mu||M.nu!=N.nu)
{
printf(" 你的输入不满足矩阵相加的条件!\n");
exit(1);
}
Q->mu=M.mu;Q->nu=M.nu;
for(p=1,q=1;p<=M.tu&&q<=N.tu;)
{
if(M.data[p].i==N.data[q].i)
{
if(M.data[p].j==N.data[q].j)
{
Q->data[k].i=M.data[p].i;
Q->data[k].j=M.data[p].j;
Q->data[k].e=M.data[p].e+N.data[q].e;
p++;q++;k++;
}
else if(M.data[p].j<N.data[q].j)
{
Q->data[k].i=M.data[p].i;
Q->data[k].j=M.data[p].j;
Q->data[k].e=M.data[p].e;

k++;p++;

}

else if(M.data[p].j>N.data[q].j)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=N.data[q].e;

k++;p++;

}

}

else if(M.data[p].i<N.data[q].i)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e;

k++;p++;

}

else if(M.data[p].i>N.data[q].i)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=N.data[q].e;
k++;q++;

}

}

if(p!=M.tu+1)

for(;p<=M.tu;p++)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e;

k++;

}

if(q!=N.tu+1)

for(;q<=N.tu;q++)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=N.data[q].e;

k++;

}

}

void SubRLSMatrix(RLSMatrix M,RLSMatrix N,RLSMatrix *Q)  //稀疏矩阵相减

{

int p,q,k=1;

if(M.mu!=N.mu||M.nu!=N.nu)

{

printf(" 你的输入不满足矩阵相减的条件!\n");

exit(1);

}

Q->mu=M.mu;Q->nu=M.nu;

for(p=1,q=1;p<=M.tu&&q<=N.tu;)

{

if(M.data[p].i==N.data[q].i)

{

if(M.data[p].j==N.data[q].j)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e-N.data[q].e;

p++;q++;k++;

}

else if(M.data[p].j<N.data[q].j)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e;

k++;p++;

}

else if(M.data[p].j>N.data[q].j)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=-N.data[q].e;
k++;p++;

}

}

else if(M.data[p].i<N.data[q].i)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e;

k++;p++;

}

else if(M.data[p].i>N.data[q].i)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=-N.data[q].e;

k++;q++;

}

}

if(p!=M.tu+1)

for(;p<=M.tu;p++)

{

Q->data[k].i=M.data[p].i;

Q->data[k].j=M.data[p].j;

Q->data[k].e=M.data[p].e;

k++;

}

if(q!=N.tu+1)

for(;q<=N.tu;q++)

{

Q->data[k].i=N.data[q].i;

Q->data[k].j=N.data[q].j;

Q->data[k].e=-N.data[q].e;

k++;

}

}

int MulTSMatrix(RLSMatrix M,RLSMatrix N,RLSMatrix *Q)//稀疏矩阵相乘
{
int ccol=0,tp,brow,t,arow,p,q,i;
int ctemp[MAXSIZE+1];
if(M.nu!=N.mu)
{
printf(" 你的输入不满足矩阵相乘的条件!\n");
return 0;
}
Q->mu=M.mu;
Q->nu=N.nu;
Q->tu=0;
if(M.tu*N.tu!=0)
{
for(arow=1;arow<=M.mu;++arow)
{
for(i=1;i<=N.nu;i++)
ctemp[i]=0;
Q->rpos[arow]=Q->tu+1;
if(arow<M.mu)  tp=M.rpos[arow+1];
else  tp=M.tu+1;
for(p=M.rpos[arow];p<tp;++p)
{
brow=M.data[p].j;
if(brow<N.mu)  t=N.rpos[brow+1];
else  t=N.tu+1;

for(q=N.rpos[brow];q<t;++q)

{

ccol=N.data[q].j;

ctemp[ccol]+=M.data[p].e*N.data[q].e;

}

}

for(ccol=1;ccol<=Q->nu;++ccol)

{

if(ctemp[ccol])

{

if(++Q->tu>MAXSIZE)  return 0;

Q->data[Q->tu].i=arow;

Q->data[Q->tu].j=ccol;

Q->data[Q->tu].e=ctemp[ccol];

}

}

}

}

return 1;

}   

void PrintSMatrix(RLSMatrix Q)          //输出稀疏矩阵

{

int k=1,row,line;

printf("\n运算结果: ");

if(Q.tu==0)  printf("0");

else

{

for(row=1;row<=Q.mu;row++)

{

for(line=1;line<=Q.nu;line++)

{

if(Q.data[k].i==row&&Q.data[k].j==line)

	printf("%d ",Q.data[k++].e);

else 
	printf("0  ");

}

printf("\n\t  ");

}

}

}

int main()

{

RLSMatrix M,N,Q;

int i;

system("cls");

printf("        ┏━━━━━━━━━━━━━━━━━━━━━━━━━┓\n");

printf("        ┃㊣必做题，稀疏矩阵运算器㊣┃\n");

printf("        ┃姓名：刘*┃\n");

printf("        ┃学号：1934*******┃\n");
printf("        ┗━━━━━━━━━━━━━━━━━━━━━━━━━┛\n");

do

{

printf("╔════════════════════════════════╗\n");

printf("║欢迎进入稀疏矩阵运算器║\n");

printf("║-------------------------------------------    ║\n");

printf("║1.矩阵相加║\n");

printf("║2.矩阵相减║\n");

printf("║3.矩阵相乘║\n");

printf("║4.退出系统║\n");

printf("╚════════════════════════════════╝\n");

printf(">>>>>>请选择功能(1-4): ");

scanf("%d",&i);

if(i==4) goto end;

else

{

printf("\n 请输入第一个矩阵M:\n");

CreateSMatrix(&M);

printf("\n 请输入第二个矩阵N:\n");

CreateSMatrix(&N);

switch(i)

{

case 1:AddRLSMatrix(M,N,&Q);break;

case 2:SubRLSMatrix(M,N,&Q);break;

case 3:MulTSMatrix(M,N,&Q);break;

default:break;

}

}

PrintSMatrix(Q);

getchar();

getchar();

end: ;

}while(i!=4);
return 0;

}