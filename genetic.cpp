//
// Created by thinkpad on 2021/1/4.
//遗传算法求解TSP问题

#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<ctime>
using namespace std;
typedef vector<int> VI;
typedef vector<VI> VVI;
#define PB push_back
#define MP make_pair

int next_int()
{
    return rand()*(RAND_MAX+1)+rand();
}
double next_double()
{
    return (double(rand()*(RAND_MAX+1)+rand()))/((RAND_MAX+1)*RAND_MAX+RAND_MAX);
}
void pmx(VI& a,VI& b,int pointcnt)//PMX交叉
{
    int sa=next_int()%pointcnt,sb=next_int()%pointcnt;//随机选择两交叉位
    int temp;
    if (sa>sb)
    {
        temp=sa;
        sa=sb;
        sb=temp;
    }//保证交叉位sa<=sb
    VI aa(pointcnt),bb(pointcnt);
    int i;
    for(i=0;i<pointcnt;i++)
    {
        aa[i]=a[i],bb[i]=b[i];
    }
    VI m1(pointcnt,-1);
    VI m2(pointcnt,-1);
    VI v1tov2(pointcnt,-1);
    for(i=sa;i<=sb;i++)
    {
        m1[aa[i]]=-2;	//m1存放aa非交叉段可代换的基因
        m2[bb[i]]=-2;	//m2存放aa非交叉段可代换的基因
    }
    for(i=0;i<pointcnt;i++)	{

        if(m2[i]==m1[i])//去掉m1和m2中重复代换的基因
        {
            m2[i]=-1;
            m1[i]=-1;
        }
    }
    int aaa=0;
    for(i=sa;i<=sb;i++)
    {
        if ((m1[aa[i]]==-2)&&(m2[bb[i]]==-2))
        {
            v1tov2[aa[i]]=bb[i];//v1tov2存放首先可以确定的互换基因
            m1[aa[i]]=-1;
            m2[bb[i]]=-1;
            aaa++;
        }
    }
    if(aaa!=(sb-sa+1))
    {
        for(i=0;i<pointcnt;i++)
            if (m1[i]==-2)
            {
                int aflag=0;
                for(int j=0;j<pointcnt;j++)
                    if( (m2[j]==-2) && (aflag==0))//寻找并确定可以互换的基因
                    {
                        v1tov2[i]=j;
                        aflag=1;
                        aaa++;
                        m2[j]=-1;
                        m1[i]=-1;
                    }
            }

    }
    for(i=sa;i<=sb;i++)	{

        swap(aa[i],bb[i]);	//交换sa到sb之间的基因串
    }

    for(i=0;i<pointcnt;i++)//查找染色体aa中sa之前和sb之后的基因是否有重复
    {
        if ((i<sa)||(i>sb))
            for (int j=sa;j<=sb;j++)
            {
                if(aa[i]==aa[j])  //有重复基因
                {
                    for(int k=0;k<pointcnt;k++)
                        if(aa[i]==v1tov2[k])
                            aa[i]=k;	//进行互换

                }
            }

    }

    for(i=0;i<pointcnt;i++)//查找染色体bb中sa之前和sb之后的基因是否有重复
    {
        if ((i<sa)||(i>sb))
            for (int j=sa;j<=sb;j++)
            {
                if(bb[i]==bb[j]) //有重复基因
                    bb[i]=v1tov2[bb[i]];	//进行互换
            }

    }
    a=aa;
    b=bb;

}

vector<double> x,y;
double fitness(const VI& v,int pointcnt)//计算适应度
{
    double r=0;
    for(int i=0;i<pointcnt;i++)
    {
        double dx=x[v[i]]-x[v[(i+1)%pointcnt]];
        double dy=y[v[i]]-y[v[(i+1)%pointcnt]];
        r+=sqrt(dx*dx+dy*dy);//个体的适应度为相邻两城市之间的距离平方的平方根和
    }
    return 1.0/r;
}
void change0(vector<int>& K,int N)//变异策略:两点互换
{
    int i=next_int()%N;
    int d=next_int()%(N-1);
    int j=(i+1+d)%N;
    swap(K[i],K[j]);
}

void mutate(VI& route,int mutate_type,int pointcnt)
{
    if(mutate_type==0)//两点互换
        change0(route,pointcnt);
}
bool pair_dec(const pair<double,VI*>& a,const pair<double,VI*>& b)
{
    return a>b;
}

class other_population
{
public:
    int popsize,pointcnt;//种群规模,染色体长度
    double pc,pm;//交叉概率,变异概率
    vector<pair<double,VI*> >pop;//种群
    pair<double,VI*> bestofpop;//最好个体
    int cross_type;//交叉类型
    int mutate_type;//变异类型
    int make_p;//个体概率分配策略类型
    int select_type;//个体选择类型
    int toursize;//竞赛规模
    double bestp;//最好个体选择概率
    other_population(int a,int b,int c,int f,int g,double d,double e,int h,double j,int m)
    {
        popsize=a,pointcnt=b,cross_type=c,mutate_type=f,make_p=g,pc=d,pm=e,toursize=h,bestp=j,select_type=m;
        for(int i=0;i<popsize;i++)//初始化种群
        {
            VI* v=new VI(pointcnt);
            for(int j=0;j<pointcnt;j++)
                (*v)[j]=j;
            random_shuffle(v->begin(),v->end());
            pop.PB(MP(fitness(*v,pointcnt),v));
        }
        sort(pop.begin(),pop.end(),pair_dec);
        bestofpop.first=pop[0].first;//初始时最好个体的适应度
        bestofpop.second=new VI(*pop[0].second);//初始时最好个体的染色体
    }
    ~other_population()
    {
        for(int i=0;i<pop.size();i++)
            delete pop[i].second;
        delete bestofpop.second;
    }
    void next()//产生下一代种群
    {
        vector<double> ps(popsize);
        if(make_p==0) //按适应度比例分配个体的选择概率
        {
            double sum=0;
            for(int i=0;i<popsize;i++)
                sum+=pop[i].first;
            for(int i=0;i<popsize;i++)
                ps[i]=pop[i].first/sum;
        }

        if(select_type==0)//轮盘赌选择个体
        {
            vector<pair<double,VI*> > select_res;
            vector<double> addsum(popsize);
            for(int i=0;i<popsize-1;i++)//计算个体的累计概率
            {
                if(i==0)
                    addsum[i]=ps[0];
                else
                    addsum[i]=addsum[i-1]+ps[i];
            }
            addsum[popsize-1]=1.5;
            for(int i=0;i<popsize;i++)
            {
                double rd=next_double();
                int r=lower_bound(addsum.begin(),addsum.end(),rd)-addsum.begin();
                VI* v=new VI(*pop[r].second);
                select_res.PB(MP(fitness(*v,pointcnt),v));
            }
            for(int i=0;i<popsize;i++)
                delete pop[i].second;
            pop=select_res;
        }
        for(int cc=0;cc<popsize/2;cc++)//随机选择两个个体,然后进行交叉
        {
            int a=next_int()%popsize;
            int b=(a+1+(next_int()%(popsize-1)))%popsize;
            if(next_double()<pc)//随机数小于交叉概率,进行交叉
            {
                if(cross_type==0)//pmx交叉
                    pmx(*pop[a].second,*pop[b].second,pointcnt);

                pop[a].first=fitness(*pop[a].second,pointcnt);//计算交叉后个体a的适应度
                if(bestofpop.first<pop[a].first)//更新最好个体
                {
                    bestofpop.first=pop[a].first;
                    delete bestofpop.second;
                    bestofpop.second=new VI(*pop[a].second);
                }
                pop[b].first=fitness(*pop[b].second,pointcnt);//计算交叉后个体b的适应度
                if(bestofpop.first<pop[b].first)//更新最好个体
                {
                    bestofpop.first=pop[b].first;
                    delete bestofpop.second;
                    bestofpop.second=new VI(*pop[b].second);
                }
            }
        }
        for(int i=pop.size()-1;i>=0;i--)//进行变异
            if(next_double()<pm)//随机数小于变异概率,进行变异
            {
                mutate(*pop[i].second,mutate_type,pointcnt);//变异
                pop[i].first=fitness(*pop[i].second,pointcnt);//计算变异后个体的适应度
            }
        sort(pop.begin(),pop.end(),pair_dec);//从大到小排序
        if(bestofpop.first<pop[0].first)//更新最好个体
        {
            delete bestofpop.second;
            bestofpop.first=pop[0].first;
            bestofpop.second=new VI(*pop[0].second);
        }
    }
};

int main()
{
    srand((unsigned)time(NULL));
    int CASNUM,POINTCNT,POPSIZE,GENERATIONS;
    //scanf("%d",&CASNUM);//输入实验次数
    CASNUM=10;//输入实验次数
    //scanf("%d%d%d",&POINTCNT,&POPSIZE,&GENERATIONS);//输入染色体长度（城市数），种群规模，最大迭代步数
    POINTCNT=10, POPSIZE=100,GENERATIONS=100;//输入染色体长度（城市数），种群规模，最大迭代步数
    x.resize(POINTCNT);
    y.resize(POINTCNT);
    x[0]=0, x[1]=1.1,x[2]=3.5,x[3]=3,x[4]=7,x[5]=8,x[6]=4,x[7]=4.5,x[8]=9,x[9]=2;
    y[0]=1.1,y[1]=3,y[2]=2,y[3]=4,y[4]=5.1,y[5]=8,y[6]=4,y[7]=4.5,y[8]=9,y[9]=2;
    cout<<"城市数="<<POINTCNT<<endl;
    cout<<"各城市坐标:"<<endl;
    for(int i=0;i<POINTCNT;i++)
    {
        //scanf("%lf%lf",&x[i],&y[i]);//输入各个城市的坐标
        cout<<"["<<x[i]<<", "<<y[i]<<"]"<<endl;//输出各个城市的坐标
    }
    int select_type,make_p_type,k,cross_type,mutate_type;
    double q,pc,pm;
    //scanf("%d%d%d",&select_type,&make_p_type,&k);//输入个体选择方法类型，个体选择概率分配类型，竞赛规模
    //scanf("%lf%lf%lf",&q,&pc,&pm);//输入最好个体选择概率，交叉概率，变异概率
    //scanf("%d%d",&cross_type,&mutate_type);//输入交叉类型，变异类型
    select_type=0,make_p_type=0,k=5;//输入个体选择方法类型，个体选择概率分配类型，竞赛规模
    q=0.5,pc=0.85,pm=0.15;//输入最好个体选择概率，交叉概率，变异概率
    cross_type=0,mutate_type=0;//输入交叉类型，变异类型

    double best=1e9,worst=0,sum=0;
    VI res;
    for(int cas=0;cas<CASNUM;cas++)//
    {
        other_population gen(POPSIZE,POINTCNT,cross_type,mutate_type,make_p_type,pc,pm,k,q,select_type);

        for(int g=0;g<GENERATIONS;g++)//进行迭代进化
            gen.next();
        if(best>1.0/gen.bestofpop.first)//更新历次最好适应度
        {
            best=1.0/gen.bestofpop.first;
            res=*gen.bestofpop.second;//存放最好个体的染色体
        }
        if(worst<1.0/gen.bestofpop.first)//更新历次最差适应度
            worst=1.0/gen.bestofpop.first;
        sum+=1.0/gen.bestofpop.first;//计算各次最好个体的适应度之和
    }
    sum/=CASNUM;//计算平均适应度
    cout<<endl;
    cout<<"历次最好适应度："<<best<<"\n"<<"历次最差适应度："<<worst<<"\n"<<"平均适应度："<<sum<<"\n";
    cout<<"输出最好解：";
    for(int i=0;i<POINTCNT;i++)//输出解
    {
        cout<<res[i];//输出各城市
        if (i<POINTCNT-1)
            cout<<"-";
    }
    cout<<endl;

    return 0;
}