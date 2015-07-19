#include <stdio.h>  
#include <conio.h>  
#include <stdlib.h>  
#include <time.h>  
#include <iostream>  

using namespace  std;

typedef struct Chrom // 结构体类型，为单个染色体的结构
{
	short int bit[6];//一共6bit来对染色体进行编码，其中1位为符号位。取值范围-64~+64  
	int fit ;//适应值
	double rfit;//相对的fit值，即所占的百分比  
	double cfit;//积累概率  
}chrom;

//定义将会用到的几个函数
void *evpop (chrom popcurrent[4]);//进行种群的初始化  
int x (chrom popcurrent);
int y (int x);
void *pickchroms (chrom popcurrent[4]);//选择操作,遗传
void *pickchroms_new (chrom popcurrent[4]); // 基于概率分布
void *crossover (chrom popnext[4]);//交叉操作
void *mutation (chrom popnext[4]);//突变，编译
double r8_uniform_ab ( double a, double b, int &seed );//生成a~b之间均匀分布的数字


chrom popcurrent [4];// 初始种群规模为；  
chrom popnext [4];// 更新后种群规模仍为；  


void main1()
{
	/*迭代次数*/
	int num ;
	int i ,j, l,Max ,k;
	/*函数最大值*/
	Max=0;

	printf("\nWelcome to the Genetic Algorithm！\n");  //   
	printf("The Algorithm is based on the function y = -x^2 + 5 to find the maximum value of the function.\n");  
	enter:printf ("\nPlease enter the no. of iterations\n请输入您要设定的迭代数 : ");  
	scanf_s("%d" ,&num);// 输入迭代次数，传送给参数 num

	if(num <1)                                    
		goto enter ;// 判断输入的迭代次数是否为负或零，是的话重新输入；  
	//不同的随机数可能结果不同？？那是当所设置的迭代次数过少时，染色体的基因型过早地陷入局部最优

	/*种下随机数种子*/
	srand(time(0));
	evpop(popcurrent );// 随机产生初始种群；  
	// 是否需要指定x的取值范围呢？6bit来表示数字，第一位为符号位，5bit表示数字大小。所以，取值范围为-32~+31
	Max = popcurrent[0].fit; // 对Max值进行初始化

	for(i =0;i< num;i ++)                          // 开始迭代；  
	{  

		printf("\ni = %d\n" ,i);// 输出当前迭代次数；  

		for(j =0;j<4; j++)  
		{  
			popnext[j ]=popcurrent[ j];// 更新种群；  
		}  

		pickchroms(popnext );// 挑选优秀个体；遗传
		crossover(popnext );// 交叉得到新个体；
		mutation(popnext );// 变异得到新个体；

		for(j =0;j<4; j++)   
		{  
			popcurrent[j]=popnext[j];// 种群更替；  
		}  

	}  // 等待迭代终止；

	//对于真正随机数是需要注意取较大的迭代次数  
	for(l =0;l<3; l++)  
	{  
		if(popcurrent [l]. fit > Max )  
		{  
			Max=popcurrent [l]. fit;  
			k=x(popcurrent [l]);//此时的value即为所求的x值  
		}  
	}  

	printf("\n当x等于%d时，函数得到最大值为： %d ",k ,Max);  
	printf("\nPress any key to end ! " );  

	flushall();// 清除所有缓冲区；  
	getche(); // 从控制台取字符，不以回车为结束；  

}                                               


/************************************************************************/
/* 随机生成种群                                                         */
/************************************************************************/
void *evpop (chrom popcurrent[4])
{  
	int i ,j, value1;  
	int random ;  
	double sum=0;  

	for(j =0;j<4; j++)                            // 从种群中的第1个染色体到第4个染色体  
	{  
		for(i =0;i<6; i++)                       // 从染色体的第1个基因位到第6个基因位  
		{  
			random=rand ();                     // 产生一个随机值  
			random=(random %2);                 // 随机产生0或者1  
			popcurrent[j ].bit[ i]=random ;       // 随机产生染色体上每一个基因位的值，或；  
		}    

		value1=x (popcurrent[ j]);// 将二进制换算为十进制，得到一个整数值；  
		popcurrent[j ].fit= y(value1); // 计算染色体的适应度值  
		sum = sum + popcurrent[j ].fit;  
		printf("\n popcurrent[%d]=%d%d%d%d%d%d  value=%d  fitness = %d",j, popcurrent[j ].bit[5], popcurrent[j ].bit[4], popcurrent[j ].bit[3], popcurrent[j ].bit[2], popcurrent[j ].bit[1], popcurrent[j ].bit[0], value1,popcurrent [j]. fit);   
		// 输出整条染色体的编码情况，  
	}  
	//计算适应值得百分比，该参数是在用轮盘赌选择法时需要用到的  
	for (j = 0; j < 4; j++)  
	{  
		popcurrent[j].rfit = popcurrent[j].fit/sum;  
		popcurrent[j].cfit = 0;//将其初始化为0  
	}  
	return(0);                  
}                                         


int x (chrom popcurrent)  // 函数：将二进制换算为十进制；  
{//此处的染色体长度为，其中个表示符号位  

	int z ;  
	z=(popcurrent .bit[0]*1)+( popcurrent.bit [1]*2)+(popcurrent. bit[2]*4)+(popcurrent .bit[3]*8)+( popcurrent.bit [4]*16);  

	if(popcurrent .bit[5]==1)  // 考虑到符号；  
	{  
		z=z *(-1);                               
	}  

	return(z );                             
}

//需要能能够从外部直接传输函数，加强鲁棒性  
int y (int x)// 函数：求个体的适应度；  
{  
	int y ;  
	y=-(x *x)+50;                                // 目标函数：y= - ( x^ 2 ) +5；  
	return(y );               
}   

//基于轮盘赌选择方法，进行基因型的选择  
void *pickchroms_new (chrom popnext[4])//计算概率  
{  
	int men;  
	int i;int j;  
	double p;  
	double sum=0.0;  
	//find the total fitness of the population  
	for (men = 0; men < 4; men++ )  
	{  
		sum = sum + popnext[men].fit;  
	}  
	//calculate the relative fitness of each member  
	for (men = 0; men < 4; men++ )  
	{  
		popnext[men].rfit = popnext[men].fit / sum;  
	}  
	//calculate the cumulative fitness,即计算积累概率  
	popcurrent[0].cfit = popcurrent[0].rfit;  
	for ( men = 1; men < 4; men++)  
	{  
		popnext[men].cfit = popnext[men-1].cfit + popnext[men].rfit;  
	}  

	for ( i = 0; i < 4; i++ )  
	{//产生0~1之间的随机数  
		//p = r8_uniform_ab ( 0, 1, seed );//通过函数生成0~1之间均匀分布的数字  
		p =rand()%10;//  
		p = p/10;  
		if ( p < popnext[0].cfit )  
		{  
			popcurrent[i] = popnext[0];        
		}  
		else  
		{  
			for ( j = 0; j < 4; j++ )  
			{   
				if ( popnext[j].cfit <= p && p < popnext[j+1].cfit )  
				{  
					popcurrent[i] = popcurrent[j+1];  
				}  
			}  
		}  
	}  
	//  Overwrite the old population with the new one.  
	//  
	for ( i = 0; i < 4; i++ )  
	{  
		popnext[i] = popcurrent[i];   
	}  
	return(0);  
}  

void *pickchroms (chrom popnext[4])          // 函数：选择个体；  
{  
	int i ,j;  
	chrom temp ;                                // 中间变量  
	//因此此处设计的是个个体，所以参数是  
	for(i =0;i<3; i++)                           // 根据个体适应度来排序；（冒泡法）  
	{  
		for(j =0;j<3-i; j++)  
		{  
			if(popnext [j+1]. fit>popnext [j]. fit)  
			{  
				temp=popnext [j+1];  
				popnext[j +1]=popnext[ j];  
				popnext[j ]=temp;  

			}    
		}                 
	}  
	for(i =0;i<4; i++)  
	{  
		printf("\nSorting:popnext[%d] fitness=%d" ,i, popnext[i ].fit);  
		//printf("\n" );                       
	}  
	flushall();/* 清除所有缓冲区 */                        
	return(0);  
}

double r8_uniform_ab( double a, double b, int &seed )  
{  
	{  
		int i4_huge = 2147483647;  
		int k;  
		double value;  

		if ( seed == 0 )  
		{  
			std::cerr << "\n";  
			std::cerr << "R8_UNIFORM_AB - Fatal error!\n";  
			std::cerr << "  Input value of SEED = 0.\n";  
			exit ( 1 );  
		}  

		k = seed / 127773;  

		seed = 16807 * ( seed - k * 127773 ) - k * 2836;  

		if ( seed < 0 )  
		{  
			seed = seed + i4_huge;  
		}  

		value = ( double ) ( seed ) * 4.656612875E-10;  

		value = a + ( b - a ) * value;  

		return value;  
	}  
}  


void *crossover (chrom popnext[4])              // 函数：交叉操作；  
{  

	int random ;  
	int i ;  
	//srand(time(0));   
	random=rand ();                             // 随机产生交叉点；  
	random=((random %5)+1);                     // 交叉点控制在0到5之间；  
	for(i =0;i< random;i ++)                     
	{  
		popnext[2].bit [i]= popnext[0].bit [i];   // child 1 cross over  
		popnext[3].bit [i]= popnext[1].bit [i];   // child 2 cross over  
	}  

	for(i =random; i<6;i ++)                      // crossing the bits beyond the cross point index  
	{  
		popnext[2].bit [i]= popnext[1].bit [i];    // child 1 cross over  
		popnext[3].bit [i]= popnext[0].bit [i];    // chlid 2 cross over  
	}    

	for(i =0;i<4; i++)  
	{  
		popnext[i ].fit= y(x (popnext[ i]));        // 为新个体计算适应度值；  
	}  

	for(i =0;i<4; i++)  
	{  
		printf("\nCrossOver popnext[%d]=%d%d%d%d%d%d    value=%d    fitness = %d",i, popnext[i ].bit[5], popnext[i ].bit[4], popnext[i ].bit[3], popnext[i ].bit[2], popnext[i ].bit[1], popnext[i ].bit[0], x(popnext [i]), popnext[i ].fit);   
		// 输出新个体；  
	}  
	return(0);  
}                                            


void *mutation (chrom popnext[4])               // 函数：变异操作；  
{  

	int random ;  
	int row ,col, value;  
	//srand(time(0));   
	random=rand ()%50;  // 随机产生到之间的数；  
	//变异操作也要遵从一定的概率来进行，一般设置为0到0.5之间  
	//  
	if(random ==25)                              // random==25的概率只有2%，即变异率为，所以是以小概率进行变异！！  
	{  
		col=rand ()%6;                            // 随机产生要变异的基因位号；  
		row=rand ()%4;                            // 随机产生要变异的染色体号；  

		if(popnext [row]. bit[col ]==0)             // 1变为；  
		{  
			popnext[row ].bit[ col]=1 ;  
		}  
		else if (popnext[ row].bit [col]==1)        // 0变为；  
		{  
			popnext[row ].bit[ col]=0;  
		}  
		popnext[row ].fit= y(x (popnext[ row]));     // 计算变异后的适应度值；  
		value=x (popnext[ row]);  
		printf("\nMutation occured in popnext[%d] bit[%d]:=%d%d%d%d%d%d    value=%d   fitness=%d", row,col ,popnext[ row].bit [5],popnext[ row].bit [4],popnext[ row].bit [3],popnext[ row].bit [2],popnext[ row].bit [1],popnext[ row].bit [0],value, popnext[row ].fit);  

		// 输出变异后的新个体；  
	}                                            

	return(0);  
}     