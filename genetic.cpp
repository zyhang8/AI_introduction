//
// Created by thinkpad on 2021/1/4.
//�Ŵ��㷨���TSP����

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
void pmx(VI& a,VI& b,int pointcnt)//PMX����
{
    int sa=next_int()%pointcnt,sb=next_int()%pointcnt;//���ѡ��������λ
    int temp;
    if (sa>sb)
    {
        temp=sa;
        sa=sb;
        sb=temp;
    }//��֤����λsa<=sb
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
        m1[aa[i]]=-2;	//m1���aa�ǽ���οɴ����Ļ���
        m2[bb[i]]=-2;	//m2���aa�ǽ���οɴ����Ļ���
    }
    for(i=0;i<pointcnt;i++)	{

        if(m2[i]==m1[i])//ȥ��m1��m2���ظ������Ļ���
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
            v1tov2[aa[i]]=bb[i];//v1tov2������ȿ���ȷ���Ļ�������
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
                    if( (m2[j]==-2) && (aflag==0))//Ѱ�Ҳ�ȷ�����Ի����Ļ���
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

        swap(aa[i],bb[i]);	//����sa��sb֮��Ļ���
    }

    for(i=0;i<pointcnt;i++)//����Ⱦɫ��aa��sa֮ǰ��sb֮��Ļ����Ƿ����ظ�
    {
        if ((i<sa)||(i>sb))
            for (int j=sa;j<=sb;j++)
            {
                if(aa[i]==aa[j])  //���ظ�����
                {
                    for(int k=0;k<pointcnt;k++)
                        if(aa[i]==v1tov2[k])
                            aa[i]=k;	//���л���

                }
            }

    }

    for(i=0;i<pointcnt;i++)//����Ⱦɫ��bb��sa֮ǰ��sb֮��Ļ����Ƿ����ظ�
    {
        if ((i<sa)||(i>sb))
            for (int j=sa;j<=sb;j++)
            {
                if(bb[i]==bb[j]) //���ظ�����
                    bb[i]=v1tov2[bb[i]];	//���л���
            }

    }
    a=aa;
    b=bb;

}

vector<double> x,y;
double fitness(const VI& v,int pointcnt)//������Ӧ��
{
    double r=0;
    for(int i=0;i<pointcnt;i++)
    {
        double dx=x[v[i]]-x[v[(i+1)%pointcnt]];
        double dy=y[v[i]]-y[v[(i+1)%pointcnt]];
        r+=sqrt(dx*dx+dy*dy);//�������Ӧ��Ϊ����������֮��ľ���ƽ����ƽ������
    }
    return 1.0/r;
}
void change0(vector<int>& K,int N)//�������:���㻥��
{
    int i=next_int()%N;
    int d=next_int()%(N-1);
    int j=(i+1+d)%N;
    swap(K[i],K[j]);
}

void mutate(VI& route,int mutate_type,int pointcnt)
{
    if(mutate_type==0)//���㻥��
        change0(route,pointcnt);
}
bool pair_dec(const pair<double,VI*>& a,const pair<double,VI*>& b)
{
    return a>b;
}

class other_population
{
public:
    int popsize,pointcnt;//��Ⱥ��ģ,Ⱦɫ�峤��
    double pc,pm;//�������,�������
    vector<pair<double,VI*> >pop;//��Ⱥ
    pair<double,VI*> bestofpop;//��ø���
    int cross_type;//��������
    int mutate_type;//��������
    int make_p;//������ʷ����������
    int select_type;//����ѡ������
    int toursize;//������ģ
    double bestp;//��ø���ѡ�����
    other_population(int a,int b,int c,int f,int g,double d,double e,int h,double j,int m)
    {
        popsize=a,pointcnt=b,cross_type=c,mutate_type=f,make_p=g,pc=d,pm=e,toursize=h,bestp=j,select_type=m;
        for(int i=0;i<popsize;i++)//��ʼ����Ⱥ
        {
            VI* v=new VI(pointcnt);
            for(int j=0;j<pointcnt;j++)
                (*v)[j]=j;
            random_shuffle(v->begin(),v->end());
            pop.PB(MP(fitness(*v,pointcnt),v));
        }
        sort(pop.begin(),pop.end(),pair_dec);
        bestofpop.first=pop[0].first;//��ʼʱ��ø������Ӧ��
        bestofpop.second=new VI(*pop[0].second);//��ʼʱ��ø����Ⱦɫ��
    }
    ~other_population()
    {
        for(int i=0;i<pop.size();i++)
            delete pop[i].second;
        delete bestofpop.second;
    }
    void next()//������һ����Ⱥ
    {
        vector<double> ps(popsize);
        if(make_p==0) //����Ӧ�ȱ�����������ѡ�����
        {
            double sum=0;
            for(int i=0;i<popsize;i++)
                sum+=pop[i].first;
            for(int i=0;i<popsize;i++)
                ps[i]=pop[i].first/sum;
        }

        if(select_type==0)//���̶�ѡ�����
        {
            vector<pair<double,VI*> > select_res;
            vector<double> addsum(popsize);
            for(int i=0;i<popsize-1;i++)//���������ۼƸ���
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
        for(int cc=0;cc<popsize/2;cc++)//���ѡ����������,Ȼ����н���
        {
            int a=next_int()%popsize;
            int b=(a+1+(next_int()%(popsize-1)))%popsize;
            if(next_double()<pc)//�����С�ڽ������,���н���
            {
                if(cross_type==0)//pmx����
                    pmx(*pop[a].second,*pop[b].second,pointcnt);

                pop[a].first=fitness(*pop[a].second,pointcnt);//���㽻������a����Ӧ��
                if(bestofpop.first<pop[a].first)//������ø���
                {
                    bestofpop.first=pop[a].first;
                    delete bestofpop.second;
                    bestofpop.second=new VI(*pop[a].second);
                }
                pop[b].first=fitness(*pop[b].second,pointcnt);//���㽻������b����Ӧ��
                if(bestofpop.first<pop[b].first)//������ø���
                {
                    bestofpop.first=pop[b].first;
                    delete bestofpop.second;
                    bestofpop.second=new VI(*pop[b].second);
                }
            }
        }
        for(int i=pop.size()-1;i>=0;i--)//���б���
            if(next_double()<pm)//�����С�ڱ������,���б���
            {
                mutate(*pop[i].second,mutate_type,pointcnt);//����
                pop[i].first=fitness(*pop[i].second,pointcnt);//��������������Ӧ��
            }
        sort(pop.begin(),pop.end(),pair_dec);//�Ӵ�С����
        if(bestofpop.first<pop[0].first)//������ø���
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
    //scanf("%d",&CASNUM);//����ʵ�����
    CASNUM=10;//����ʵ�����
    //scanf("%d%d%d",&POINTCNT,&POPSIZE,&GENERATIONS);//����Ⱦɫ�峤�ȣ�������������Ⱥ��ģ������������
    POINTCNT=10, POPSIZE=100,GENERATIONS=100;//����Ⱦɫ�峤�ȣ�������������Ⱥ��ģ������������
    x.resize(POINTCNT);
    y.resize(POINTCNT);
    x[0]=0, x[1]=1.1,x[2]=3.5,x[3]=3,x[4]=7,x[5]=8,x[6]=4,x[7]=4.5,x[8]=9,x[9]=2;
    y[0]=1.1,y[1]=3,y[2]=2,y[3]=4,y[4]=5.1,y[5]=8,y[6]=4,y[7]=4.5,y[8]=9,y[9]=2;
    cout<<"������="<<POINTCNT<<endl;
    cout<<"����������:"<<endl;
    for(int i=0;i<POINTCNT;i++)
    {
        //scanf("%lf%lf",&x[i],&y[i]);//����������е�����
        cout<<"["<<x[i]<<", "<<y[i]<<"]"<<endl;//����������е�����
    }
    int select_type,make_p_type,k,cross_type,mutate_type;
    double q,pc,pm;
    //scanf("%d%d%d",&select_type,&make_p_type,&k);//�������ѡ�񷽷����ͣ�����ѡ����ʷ������ͣ�������ģ
    //scanf("%lf%lf%lf",&q,&pc,&pm);//������ø���ѡ����ʣ�������ʣ��������
    //scanf("%d%d",&cross_type,&mutate_type);//���뽻�����ͣ���������
    select_type=0,make_p_type=0,k=5;//�������ѡ�񷽷����ͣ�����ѡ����ʷ������ͣ�������ģ
    q=0.5,pc=0.85,pm=0.15;//������ø���ѡ����ʣ�������ʣ��������
    cross_type=0,mutate_type=0;//���뽻�����ͣ���������

    double best=1e9,worst=0,sum=0;
    VI res;
    for(int cas=0;cas<CASNUM;cas++)//
    {
        other_population gen(POPSIZE,POINTCNT,cross_type,mutate_type,make_p_type,pc,pm,k,q,select_type);

        for(int g=0;g<GENERATIONS;g++)//���е�������
            gen.next();
        if(best>1.0/gen.bestofpop.first)//�������������Ӧ��
        {
            best=1.0/gen.bestofpop.first;
            res=*gen.bestofpop.second;//�����ø����Ⱦɫ��
        }
        if(worst<1.0/gen.bestofpop.first)//�������������Ӧ��
            worst=1.0/gen.bestofpop.first;
        sum+=1.0/gen.bestofpop.first;//���������ø������Ӧ��֮��
    }
    sum/=CASNUM;//����ƽ����Ӧ��
    cout<<endl;
    cout<<"���������Ӧ�ȣ�"<<best<<"\n"<<"���������Ӧ�ȣ�"<<worst<<"\n"<<"ƽ����Ӧ�ȣ�"<<sum<<"\n";
    cout<<"�����ý⣺";
    for(int i=0;i<POINTCNT;i++)//�����
    {
        cout<<res[i];//���������
        if (i<POINTCNT-1)
            cout<<"-";
    }
    cout<<endl;

    return 0;
}