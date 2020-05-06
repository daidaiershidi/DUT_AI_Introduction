#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>

using std::sort;
using std::fabs;

const int MAX_DIMENSION = 2;//����ά��
const int MAX_SAMPLES = 4;//���ݸ��� 
double x[MAX_SAMPLES][MAX_DIMENSION];//���ݼ� 
double y[MAX_SAMPLES];//��ǩ 
double alpha[MAX_SAMPLES];//��������ϵ�� 
double w[MAX_DIMENSION];//ϵ��w 
double b;//ϵ��b 
double c;//�ɳ���(������k����ʽ���ÿ��ϵ������k��c����1) 
double eps = 1e-6;//��� 
struct _E{ 
    double val;//g - y,���ڸ���alpha 
    int index;//��� 
}E[MAX_SAMPLES]; 

bool cmp(const _E & a, const _E & b)//sort����ʱ�õıȽϺ��� 
{
    return a.val < b.val;
}

int num_dimension;
int num_samples;

double max(double a,double b)//ȡ���
{
    return a>b?a:b;
}

double min(double a,double b)//ȡ��С
{
    return a>b?b:a;
}

double kernal(double x1[], double x2[], double dimension)//���ڻ� 
{
    double ans = 0 ;
    for(int i = 0 ; i < dimension; i++)
    {
        ans += x1[i]*x2[i];
    }
    return ans;
}

double target_function()
{
	//���µ�alpha1,alpha2������W(a)��
	//W(a)=1/2*K11*a1*a1 + 1/2*K22*a2*a2 + g1 + g2 + ans + W(a3~an)
	//g1,g2��ָ���������g������ans���Ǳ�������ans 
    double ans = 0;
    for(int i = 0 ; i < num_samples; i++)
    {
        for(int j = 0 ; j < num_samples; j++)
        {
            ans += alpha[i]*alpha[j]*y[i]*y[j]*kernal(x[i],x[j],num_dimension);
        }
    }

    for(int i = 0 ; i < num_samples; i++)
    {
        ans -= alpha[i];
    }

    return ans;
}


double g(double _x[], int dimension) 
{
	//����W(a)��һ���֣��������target_function����,��ͬ����W(a) 
	//ͬʱgҲ��alpha����ʱ��һ���֣����������ﵥ���г� 
    double ans = b;

    for(int i = 0 ; i < num_samples; i++)
    {
        ans += alpha[i]*y[i]*kernal(x[i],_x,dimension);
    }

    return ans;
}

bool satisfy_constrains(int i, int dimension) 
{
	//����KKT�����Ͳ��ø��� 
	//�ж��Ƿ�Ϊ����KKT�����������������������ʱʵ���������ֵ� 
    if(alpha[i] == 0)//��һ���Ѿ��ֶ��˵ĵ� 
    {
        if(y[i]*g(x[i], dimension) >= 1)
        return true;
        else
        return false;
    }
    else if( alpha[i] > 0 && alpha[i] < c)//�߽��ϣ�֧������ 
    {
        if(y[i] * g(x[i], dimension) ==  1)
        return true;
        else
        return false;
    }
    else//alpha == C���߽�֮��ĵ� 
    {
        if(y[i] * g(x[i], dimension) <=  1)
        return true;
        else
        return false;
    }
}


double calE(int i, int dimension)//����E(����alphaʱ�Ƶ����ģ���alpha�Ƶ��ĺ���) 
{
    return g(x[i], dimension) - y[i];
}

void calW()//����W(a)����w (�����µ�alpha����)
{
    for(int i = 0 ; i < num_dimension; i++)
    {
        w[i] = 0;
        for(int j = 0 ; j < num_samples; j++)
        {
            w[i] += alpha[j] * y[j] * x[j][i];
        }
    }
    return ;
}

void calB()//�����ʼ�Ĺ�ʽ����b(�����µ�alpha����)
{
    double ans = y[0];
    for(int i = 0 ;  i < num_samples ; i++)
    {
        ans -= y[i]*alpha[i]*kernal(x[i], x[0], num_dimension);
    }
    b = ans;
    return;
}


void recalB(int alpha1index,int alpha2index, int dimension, double alpha1old, double alpha2old)
{
	//����Y(W^T * X + b) - 1 == 0�� ����֮ǰ�����õ�����alpha��������ʵ���Ͽɱ��滻���滻�������������ʽ  
	double alpha1new = alpha[alpha1index];
    double alpha2new = alpha[alpha2index];

    alpha[alpha1index] = alpha1old;
    alpha[alpha2index] = alpha2old;

    double e1 = calE(alpha1index, num_dimension);
    double e2 = calE(alpha2index, num_dimension);

    alpha[alpha1index] = alpha1new;
    alpha[alpha2index] = alpha2new;
	//b1 = �ɵ�b1 - E1 - y1(�µ�a1 - �ɵ�a1)*K(x1, x1) - y2(�µ�a2 - �ɵ�a2)*K(x1, x2) 
    //b2 = �ɵ�b2 - E2 - y1(�µ�a1 - �ɵ�a1)*K(x1, x2) - y2(�µ�a2 - �ɵ�a2)*K(x2, x2) 
    double b1new = -e1 - y[alpha1index]*kernal(x[alpha1index], x[alpha1index], dimension)*(alpha1new - alpha1old);
    b1new -= y[alpha2index]*kernal(x[alpha2index], x[alpha1index], dimension)*(alpha2new - alpha2old) + b;

    double b2new = -e2 - y[alpha1index]*kernal(x[alpha1index], x[alpha2index], dimension)*(alpha1new - alpha1old);
    b1new -= y[alpha2index]*kernal(x[alpha2index], x[alpha2index], dimension)*(alpha2new - alpha2old) + b;
	//���³�����b��ÿ��ȡƽ��ֵ 
    b = (b1new + b2new)/2;
}

bool optimizehelp(int alpha1index,int alpha2index)
{
	//����ֻ����������y==1��y==-1������Լ����alpha[i]*y[i]==alpha[j]*y[j]==eps 
    double alpha1new = alpha[alpha1index];
    double alpha2new = alpha[alpha2index];

    double alpha1old = alpha[alpha1index];
    double alpha2old = alpha[alpha2index];

    double H,L;
	//�µ�alpha�ķ�Χ 
    if(fabs(y[alpha1index] - y[alpha2index]) > eps)//ͬ�� 
    {
        L = max(0, alpha2old - alpha1old);
        H = min(c, c + alpha2old - alpha1old);
    }
    else//y[i]!=y[j]����ͬ�� 
    {
        L = max(0, alpha2old + alpha1old - c);
        H = min(c, alpha2old + alpha1old);
    }

    //���µ�alpha
	//��W(a)չ�����Ѱ���alpha1��alpha2��ʽ���ó�����������ɲ���
	//����alpha[i]*y[i]==alpha[j]*y[j]��ʽ����alpha2����alpha1������W(a)���ھ�ֻʣalpha1�ˣ���ʱ���� 
    double lena = kernal(x[alpha1index], x[alpha1index], num_dimension) + kernal(x[alpha2index], x[alpha2index], num_dimension) - 2*kernal(x[alpha1index], x[alpha2index], num_dimension);
    //�󵼺󣬿��Ƶ�����һ�еĹ�ʽ 
	alpha2new = alpha2old + y[alpha2index]*(calE(alpha1index, num_dimension) - calE(alpha2index, num_dimension))/lena;
    //�������L~H�͵���L��H 
	if(alpha2new > H)
    {
        alpha2new = H;
    }
    else if( alpha2new < L)
    {
        alpha2new = L;
    }
    alpha1new = alpha1old + y[alpha1index]*y[alpha2index]*(alpha2old - alpha2new);
    double energyold = target_function();

    alpha[alpha1index] = alpha1new;
    alpha[alpha2index] = alpha2new;

    double gap = 0.001;
	//����µ�b 
    recalB(alpha1index, alpha2index, num_dimension, alpha1old, alpha2old);
    return true;
}

bool optimize()//���� 
{
    int alpha1index = -1;
    int alpha2index = -1;
    double alpha2new = 0;
    double alpha1new = 0;

    //����E = g - y 
    for(int i = 0 ; i < num_samples; i++)
    {
        E[i].val = calE(i, num_dimension);
        E[i].index = i;
    }
	//ÿ�δӲ�����KKT�����Ķ��alpha�������������е���
    //����ÿһ��alpha1
    for(int i = 0 ; i < num_samples; i++)
    {
        alpha1new = alpha[i];
		//����Ǳ߽�֮��ĵ�
        if(alpha1new > 0 && alpha1new < c)
        {
			//�ж��Ƿ�����KKT 
            if(satisfy_constrains(i, num_dimension))
            continue;
			//�����Ա�����ѡalpha 
            sort(E, E+num_samples, cmp);
            //�������0�Ҳ��Ǳ���������Сֵ�͵�ǰalpha���е��� 
            if(alpha1new > 0)
            {
                if(E[0].index == i)
                {
                    ;
                }
                else
                {
                    alpha1index = i;
                    alpha2index = E[0].index;
                    if(optimizehelp(alpha1index, alpha2index))
                    {
                        return true;
                    }
                }
            }
            //���С��0�Ҳ��Ǳ����������ֵ�͵�ǰalpha���е���
            else
            {
                if(E[num_samples-1].index == i)
                {
                    ;
                }
                else
                {
                    alpha1index = i;
                    alpha2index = E[num_samples-1].index;
                    if(optimizehelp(alpha1index, alpha2index))
                    {
                        return true;
                    }
                }
            }
            //����ÿ��alpha2���ֱ��ڵ�ǰ��alpha1һ�����KKT���� 
            for(int j = 0 ; j < num_samples; j++)
            {
                alpha2new = alpha[j];

                if(alpha2new > 0 && alpha2new < c)
                {
                    alpha1index = i;
                    alpha2index = j;
                    if(optimizehelp(alpha1index , alpha2index))
                    {
                        return true;
                    }
                }
            }
            for(int j = 0 ; j < num_samples; j++)
            {
                alpha2new = alpha[j];

                if(!(alpha2new > 0 && alpha2new < c))
                {
                    alpha1index = i;
                    alpha2index = j;
                    if(optimizehelp(alpha1index , alpha2index))
                    {
                        return true;
                    }
                }
            }
        }
    }
    //����ÿһ��alpha1
    for(int i = 0 ; i < num_samples; i++)
    {
        alpha1new = alpha[i];
        //������Ǳ߽�֮��ĵ�
        if(!(alpha1new > 0 && alpha1new < c))
        {
        	//�ж��Ƿ�����KKT
            if(satisfy_constrains(i, num_dimension))
            continue;
            //�����Ա�����ѡalpha 
            sort(E, E+num_samples, cmp);
            //�������0�Ҳ��Ǳ���������Сֵ�͵�ǰalpha���е��� 
            if(alpha1new > 0)
            {
                if(E[0].index == i)
                {
                    ;
                }
                else
                {
                    alpha1index = i;
                    alpha2index = E[0].index;
                    if(optimizehelp(alpha1index, alpha2index))
                    {
                        return true;
                    }
                }
            }
            //���С��0�Ҳ��Ǳ����������ֵ�͵�ǰalpha���е���
			else
            {
                if(E[num_samples-1].index == i)
                {
                    ;
                }
                else
                {
                    alpha1index = i;
                    alpha2index = E[num_samples-1].index;
                    if(optimizehelp(alpha1index, alpha2index))
                    {
                        return true;
                    }
                }
            }
            //����ÿ��alpha2���ֱ��ڵ�ǰ��alpha1һ�����KKT����
            for(int j = 0 ; j < num_samples; j++)
            {
                alpha2new = alpha[j];

                if(alpha2new > 0 && alpha2new < c)
                {
                    alpha1index = i;
                    alpha2index = j;
                    if(optimizehelp(alpha1index , alpha2index))
                    {
                        return true;
                    }
                }
            }
            for(int j = 0 ; j < num_samples; j++)
            {
                alpha2new = alpha[j];

                if(!(alpha2new > 0 && alpha2new < c))
                {
                    alpha1index = i;
                    alpha2index = j;
                    if(optimizehelp(alpha1index , alpha2index))
                    {
                        return true;
                    }
                }
            }
        }
    }
  
    return false;
}

bool check()
{
	//���ÿһ���㶼�ǲ��Ǻϸ���  
    double sum = 0;
    for(int i = 0 ; i < num_samples; i++)
    {
        sum += alpha[i] * y[i];
        if(!(0 <= alpha[i] && alpha[i] <= c))//�жϸ����㶼�ǲ���֧������ 
        {
            printf("alpha[%d]: %lf wrong\n", i, alpha[i]);
            return false;
        }
        if(!satisfy_constrains(i, num_dimension))//�ж�ÿ�����Ƿ��Ѿ���ȷ���� 
        {
            printf("alpha[%d] not satisfy constrains\n", i);
            return false;
        }
    }
	//�ж��Ƿ������������Լ��alpha[i]*y[i] == 0 
    if(fabs(sum) > eps)
    {
        printf("Sum = %lf\n", sum);
        return false;
    }
    return true;
}
/*
min 1/2*||w||^2
s.t.  (w[i]*x[i] + b[i] - y[i]) >= 0;
*/
/*
step 1: cal alpha[]
step 2: cal w,b
*/

/*
min(para alpha) 1/2*sum(i)sum(j)(alpha[i]*alpha[j]*y[i]*y[j]*x[i]*x[j]) - sum(alpha[i])
s.t. sum(alpha[i] * y[i]) = 0
C>= alpha[i] >= 0
*/

int main()
{
    scanf("%d%d", &num_samples, &num_dimension);//�������ݸ���������ά��

	//��������
    for(int i = 0 ; i < num_samples; i++)
    {
        for(int j = 0; j < num_dimension; j++)
        {
            scanf("%lf",&x[i][j]);//��i�����ݵĸ���ά������ 
        }
        scanf("%lf",&y[i]);//��i�����ݵı�ǩ 
    }
    c = 1;

    //��������ϵ����ֵ��Ϊ0
    for(int i = 0 ; i < num_samples; i++)
    {
        alpha[i] = 0;
    }

    int count = 0;
    //���и��� 
    while(optimize()){
        calB();		//����b		
        count++;
    }
    printf("count = %d\n ",count);

    calW();//����W 
    calB();//����b 

    printf("y = ");
	//���ؽ�� 
    for(int i = 0 ; i < num_dimension; i++)
    {
        printf("%lf * x[%d] + ", w[i], i);
    }
    printf("%lf\n", b);
	//����� 
    if(!check())
    printf("Not satisfy KKT.\n");
    else
    printf("Satisfy KKT\n");
}
//3 2
//3 3 1
//4 3 1
//1 1 -1
/*
min 1/2*||w||^2
s.t.  (w[i]*x[i] + b[i] - y[i]) >= 0;
*/

/*
min(para alpha) 1/2*sum(i)sum(j)(alpha[i]*alpha[j]*y[i]*y[j]*x[i]*x[j]) - sum(alpha[i])
s.t. sum(alpha[i] * y[i]) = 0
C>= alpha[i] >= 0
*/
