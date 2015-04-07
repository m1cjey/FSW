#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[
#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
#include<omp.h>
#include<vector>
#include"function.h"


#define LOOP 0
#define UNLOOP 1

void set_BEM_model(mpsconfig *CON,vector<point2D> &static_NODE,vector<element2D> &static_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//�ÓI�ߓ_�y�їv�f�̏�񐶐�
	set_BEM_static_model(CON,static_NODE,static_ELEM,s_node_num,s_elem_num);
}

void set_BEM_static_model(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,int *s_node_num,int *s_elem_num)
{
	//���̈ȊO�̐ߓ_�E�v�f�����쐬

	int node_num=0;	//�ߓ_��
	int elemnum=0;	//�v�f��
	int ele_type=CON->get_BEM_elm_type();	//�v�f�^�C�v 0:��� 1:���`
	int B_flag=UNLOOP;		//B_flag��LOOP�Ȃ狫�E�͌q�����Ă���BUNLOOP�Ȃ�Ȃ����ĂȂ�
	int divN;	//1�ӂ̕�����
	double divL;	//������

	point2D NODE01;
	element2D ELEM01;

	double Xmin=CON->get_minX();//��͗̈�
	double Xmax=CON->get_maxX();
	double Ymin=CON->get_minY();
	double Ymax=CON->get_maxY();

	double V1=0;//CON->get_V();
	double V2=CON->get_V();		//�|�e���V�����l

	//���Ӎ쐬
	divN=50;
	divL=(Xmax-Xmin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmin+n*divL;
		NODE[node_num].r[A_Y]=Ymin;
		NODE[node_num].potential=V1;
		NODE[node_num].slop1=0;						//������
		NODE[node_num].boundary_condition=Diric;
		if(n==0) NODE[node_num].C=0.25;				//���p
		else NODE[node_num].C=0.5;
		NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		node_num++;
	}
	//�E�Ӎ쐬
	divL=(Ymax-Ymin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmax;
		NODE[node_num].r[A_Y]=n*divL+Ymin;
		NODE[node_num].potential=0;
		NODE[node_num].slop1=0;						//������
		NODE[node_num].boundary_condition=Neumn;
		NODE[node_num].C=0.5;
		if(n==0)
		{
			NODE[node_num].C=0.25;				//���p
			NODE[node_num].boundary_condition=Diric;
		}
		NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		node_num++;
	}
	//��Ӎ쐬
	divL=(Xmax-Xmin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmax-n*divL;
		NODE[node_num].r[A_Y]=Ymax;
		NODE[node_num].potential=V2;
		NODE[node_num].slop1=0;						//������
		NODE[node_num].boundary_condition=Diric;
		if(n==0) NODE[node_num].C=0.25;				//���p
		else NODE[node_num].C=0.5;
		NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		node_num++;
	}
	//���Ӎ쐬
	divL=(Ymax-Ymin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmin;
		NODE[node_num].r[A_Y]=Ymax-n*divL;
		NODE[node_num].potential=0;
		NODE[node_num].slop1=0;						//������
		NODE[node_num].boundary_condition=Neumn;
		NODE[node_num].C=0.5;
		if(n==0)
		{
			NODE[node_num].C=0.25;				//���p
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].potential=V2;
		}
		NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		node_num++;
	}
	elemnum=node_num;			//2D�ł͐ߓ_��=�v�f��
	for(int i=0;i<node_num;i++) NODE[i].particle=-1;	//�Ή�������̂��Ȃ��̂Ł|1���_�~�[�Ƃ��Ċi�[
	B_flag=LOOP;

	//if(B_flag==LOOP)
	{
		//�v�f���쐬
		for(int i=0;i<elemnum;i++)
		{
			ELEM.push_back(ELEM01);
			ELEM[i].node[0]=i;
			ELEM[i].node[1]=i+1;
			if(i+1==node_num) ELEM[i].node[1]=0;//�Ō�̗v�f�̏I�_�͐ߓ_0
			int n1=ELEM[i].node[0];
			int n2=ELEM[i].node[1];
			double X=NODE[n2].r[A_X]-NODE[n1].r[A_X];
			double Y=NODE[n2].r[A_Y]-NODE[n1].r[A_Y];
			double L=sqrt(X*X+Y*Y);			//�v�f����
			ELEM[i].L=L;
			ELEM[i].direct[A_X]=Y/L;
			ELEM[i].direct[A_Y]=-X/L;
			ELEM[i].r[A_X]=0.5*(NODE[n1].r[A_X]+NODE[n2].r[A_X]);
			ELEM[i].r[A_Y]=0.5*(NODE[n1].r[A_Y]+NODE[n2].r[A_Y]);
			ELEM[i].material=WALL;

			if(NODE[n1].boundary_condition==NODE[n2].boundary_condition) ELEM[i].boundary_condition=NODE[n1].boundary_condition;	//���[���������E�����Ȃ炻��ɏK��
			else
			{
				ELEM[i].boundary_condition=Neumn;//���[�ŋ��E�������قȂ�ꍇ�A������ިظڂ�ɲ�݂̍�ł���B���̂Ƃ���ɲ�݌^�Ƃ���B
				if(ele_type==CONSTANT)
				{
					int N=n1;
					if(NODE[n2].boundary_condition==Neumn) N=n2;
					ELEM[i].node[0]=N;		//���v�f�̏ꍇ�́A���ߓ_�Ƃ���ɲ�݌^�̕����i�[ �����ELEM[i]��ɲ�݌^�ŁA���̖@�������l�͐ߓ_N�̂���ɓ��������Ƃ��Ӗ�����
				}
			}
		}
	}

	cout<<"�Œ�ߓ_��="<<node_num<<" �Œ�v�f��="<<elemnum<<endl;
	
	*s_node_num=node_num;
	*s_elem_num=elemnum;

}

//���v�f�p�́A���I�ߓ_�E�v�f��񐶐��֐�
void set_BEM2D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<point2D> &NODE,vector<element2D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number)
{
	if(CON->get_dimention()==3) cout<<"set_BEM2D_dynaic_model_for_CONSTANT��3D�ł͎g�p�ł��܂���"<<endl;

	int node_num=0;
	int elem_num=0;
	double le=CON->get_distancebp();
	point2D NODE01;
	element2D ELEM01;
	
	int CLOSE_FLAG=ON;			//�̈悪������ƕ��Ă�����ON�@�����łȂ��Ȃ�OFF

	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//�������@���x�N�g��

	int *BEMID=new int[particle_number];		//BEM�ɐߓ_�Ƃ��ďo�͂���Ȃ炻�̐ߓ_�ԍ�����́@�o�͂��Ȃ��Ȃ�-1���i�[

	//�@���x�N�g������
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	}

	for(int i=0;i<particle_number;i++) BEMID[i]=-1;		//������

	//�ߓ_������
	for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)
		{
			NODE.push_back(NODE01);
			for(int D=0;D<2;D++) NODE[node_num].r[D]=PART[i].r[D];
			//NODE[node_num].boundary_condition=Diric;
			NODE[node_num].boundary_condition=BOTH;		//�Q�}���̂Ƃ��͂�����
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//������
			NODE[node_num].slop2=0; 				//������
			NODE[node_num].C=0.5;					//���v�f��0.5
			NODE[node_num].L=le;					//����ł����H
			NODE[node_num].particle=i;				//�Ή����闱�q�ԍ��L��
			BEMID[i]=node_num;							//BEM�ɏo�͂��邵�邵
			node_num++;
		}
	}///

	for(int n=0;n<node_num;n++)
	{
		int i=NODE[n].particle;//�Ή����闱�q
		double nx=direct[A_X][i];				//�������P�ʖ@���x�N�g���ł��邱�Ƃɒ���
		double ny=direct[A_Y][i];
		int J=i;								//�v�f�̑��[���\�����闱�q�ԍ��@i�ŏ�����
		double mindis=100;
		for(int k=0;k<PART[i].N3;k++)			//�ߗח��q�̂Ȃ��ŁA�@���x�N�g���̉E���ɑ��݂���ŒZ�����̗��q��T��
		{
			int j=PART[i].NEI3[k];
			if(BEMID[j]>=0)		//�Ή�����BEM�ߓ_�����݂���Ȃ�
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double dis=sqrt(X*X+Y*Y);
				double COS=(X*nx+Y*ny)/(dis*1.0);//���qi�̖@���x�N�g���ƁArij�̂Ȃ��p�x(cos��)
				NODE[n].C=acos(COS);				//���p
				//cout<<360*NODE[n].C/(2*PI)<<endl;
				double O_product=X*ny-Y*nx;			//�O�ρ@rij�~direct
				if(O_product>0)	//�O�ς����Ƃ������Ƃ́A���qj�͖@���x�N�g���̉E���ɑ��݂���
				{
					if(dis<mindis)
					{
						J=j;
						mindis=dis;
					}
				}
			}
		}
		if(J!=i)		//�ꏏ�ɗv�f���\�����闱�q�����������Ȃ�
		{
			double X=PART[J].r[A_X]-PART[i].r[A_X];
			double Y=PART[J].r[A_Y]-PART[i].r[A_Y];
			double dis=mindis;
			
			//�v�f�쐬
			ELEM.push_back(ELEM01);
			ELEM[elem_num].node[0]=n;
			ELEM[elem_num].node[1]=BEMID[J];
			for(int D=0;D<2;D++) ELEM[elem_num].r[D]=0.5*(PART[i].r[D]+PART[J].r[D]);//2�_�Ԃ̒��_
			ELEM[elem_num].L=dis;
			ELEM[elem_num].direct[A_X]=-Y/dis;		//�����̗̈�Ȃ̂ŁA��͗̈�[�Œ�`�����@���x�N�g���Ƃ͒�`�����΂ƂȂ�
			ELEM[elem_num].direct[A_Y]=X/dis;
			ELEM[elem_num].map=0;
			ELEM[elem_num].material=FLUID;
			ELEM[elem_num].boundary_condition=NODE[n].boundary_condition;
			
			elem_num++;
		}
		else
		{
			cout<<"�̈悪�Ƃ��Ă��Ȃ��H"<<endl;
			cout<<n<<" "<<i<<endl;
			CLOSE_FLAG=OFF;
		}
	}//*/

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	cout<<"���I�ߓ_��="<<node_num<<" ���I�v�f��="<<elem_num<<endl;
	*d_node_num=node_num;
	*d_elem_num=elem_num;

	delete [] BEMID;

	if(CLOSE_FLAG==OFF)
	{
		ofstream fl("checkclose1.dat");
		ofstream fl2("checkclose2.dat");
		for(int n=0;n<node_num;n++) fl<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<n<<endl;//�S�ߓ_�ʒu
		
		for(int n=0;n<elem_num;n++) fl2<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<n<<endl;
		fl.close();
		fl2.close();
	}
	
}

void couple_NODE_and_ELEM(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,vector<point2D> &s_NODE,vector<element2D> &s_ELEM,vector<point2D> &dy_NODE,vector<element2D> &dy_ELEM,vector<REGION> &region)
{
	//s_NODE��dy_NODE�����̂�����NODE�Ɋi�[����
	int s_node_num=int (s_NODE.size());
	int dy_node_num=int (dy_NODE.size());
	int s_elem_num=int (s_ELEM.size());
	int dy_elem_num=int (dy_ELEM.size());
	int countN=0;
	int countE=0;

	point2D NODE01;
	element2D ELEM01;
	REGION temp;
	region.push_back(temp);
	region[0].start=0;					//�ŏ��̌v�Z�̈�͕K���v�Z�_�O����n�܂�

	int *s_Nid=new int[s_node_num];		//�ړ���̐ߓ_�ԍ��i�[�Bs_NODE[i]�������ߓ_��NODE[s_Nid[i]]�Ɋi�[�����
	int *dy_Nid=new int[dy_node_num];


	//�ߓ_��񐶐�
	for(int i=0;i<s_node_num;i++)
	{
		NODE.push_back(NODE01);
		NODE[countN].r[A_X]=s_NODE[i].r[A_X];
		NODE[countN].r[A_Y]=s_NODE[i].r[A_Y];
		NODE[countN].boundary_condition=s_NODE[i].boundary_condition;
		NODE[countN].potential=s_NODE[i].potential;
		NODE[countN].slop1=s_NODE[i].slop1;
		NODE[countN].slop2=s_NODE[i].slop2;
		NODE[countN].C=s_NODE[i].C;
		NODE[countN].L=s_NODE[i].L;
		NODE[countN].particle=s_NODE[i].particle;
		s_Nid[i]=countN;			//i�Ԗڂ�s_NODE��countN�Ԗڂ̐ߓ_
		countN++;
	}
	region.push_back(temp);//��Q�̈�
	region[1].start=countN;
	for(int i=0;i<dy_node_num;i++)
	{
		NODE.push_back(NODE01);
		NODE[countN].r[A_X]=dy_NODE[i].r[A_X];
		NODE[countN].r[A_Y]=dy_NODE[i].r[A_Y];
		NODE[countN].boundary_condition=dy_NODE[i].boundary_condition;
		NODE[countN].potential=dy_NODE[i].potential;
		NODE[countN].slop1=dy_NODE[i].slop1;
		NODE[countN].slop2=dy_NODE[i].slop2;
		NODE[countN].C=dy_NODE[i].C;
		NODE[countN].L=dy_NODE[i].L;
		NODE[countN].particle=dy_NODE[i].particle;
		dy_Nid[i]=countN;			//i�Ԗڂ�dy_NODE��countN�Ԗڂ̐ߓ_
		countN++;
	}
	region[0].end=countN;
	region[1].end=countN;
	

	//�v�f���쐬
	for(int i=0;i<s_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		ELEM[countE].node[0]=s_ELEM[i].node[0];
		ELEM[countE].node[1]=s_ELEM[i].node[1];
		ELEM[countE].r[A_X]=s_ELEM[i].r[A_X];
		ELEM[countE].r[A_Y]=s_ELEM[i].r[A_Y];
		ELEM[countE].boundary_condition=s_ELEM[i].boundary_condition;
		ELEM[countE].L=s_ELEM[i].L;
		ELEM[countE].direct[A_X]=s_ELEM[i].direct[A_X];
		ELEM[countE].direct[A_Y]=s_ELEM[i].direct[A_Y];
		ELEM[countE].map=s_ELEM[i].map;
		ELEM[countE].material=s_ELEM[i].material;
		countE++;
	}
	for(int i=0;i<countE;i++)
	{
		int n1=ELEM[i].node[0];//���̒i�K�Ŋi�[����Ă���ߓ_�ԍ��́As_NODE�̂��̂ł���ANODE�̂��̂ł͂Ȃ�
		int n2=ELEM[i].node[1];
		n1=s_Nid[n1];			//�ߓ_�ԍ���NODE��̂��̂ɏ�������
		n2=s_Nid[n2];
		ELEM[i].node[0]=n1;
		ELEM[i].node[1]=n2;
	}
	for(int i=0;i<dy_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		ELEM[countE].node[0]=dy_ELEM[i].node[0];
		ELEM[countE].node[1]=dy_ELEM[i].node[1];
		ELEM[countE].r[A_X]=dy_ELEM[i].r[A_X];
		ELEM[countE].r[A_Y]=dy_ELEM[i].r[A_Y];
		ELEM[countE].boundary_condition=dy_ELEM[i].boundary_condition;
		ELEM[countE].L=dy_ELEM[i].L;
		ELEM[countE].direct[A_X]=dy_ELEM[i].direct[A_X];
		ELEM[countE].direct[A_Y]=dy_ELEM[i].direct[A_Y];
		ELEM[countE].map=dy_ELEM[i].map;
		ELEM[countE].material=dy_ELEM[i].material;
		countE++;
	}
	for(int i=s_elem_num;i<countE;i++)
	{
		int n1=ELEM[i].node[0];//���̒i�K�Ŋi�[����Ă���ߓ_�ԍ��́Ady_NODE�̂��̂ł���ANODE�̂��̂ł͂Ȃ�
		int n2=ELEM[i].node[1];
		n1=dy_Nid[n1];			//�ߓ_�ԍ���NODE��̂��̂ɏ�������
		n2=dy_Nid[n2];
		ELEM[i].node[0]=n1;
		ELEM[i].node[1]=n2;
		//cout<<n1<<" "<<n2<<endl;
	}
	
	cout<<"�S�ߓ_��="<<countN<<" �S�v�f��="<<countE<<endl;

	//�t�@�C���o�͂��č��W���m�F
	ofstream fp("BEM_input.dat");
	double times=0.0001;
	for(int i=0;i<countE;i++)
	{
		int n=ELEM[i].node[0];
	//	fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<endl;
		fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<ELEM[i].direct[A_X]*times<<"\t"<<ELEM[i].direct[A_Y]*times<<endl;
		//fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<NODE[n].potential<<endl;
	}//*/
	//for(int i=0;i<countN;i++) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Y]<<"\t"<<NODE[i].potential<<endl;
	fp.close();


	delete [] s_Nid;	
	delete [] dy_Nid;
}


