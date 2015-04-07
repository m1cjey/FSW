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

//�����̉�͗̈�쐬�֐�
void make_cube_region(mpsconfig *CON,vector<BEMpoint3D> &NODE,int *s_node_num);


void set_BEM3D_static_model(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *s_node_num,int *s_elem_num)
{
	//���̈ȊO�̐ߓ_�E�v�f�����쐬

	int node_num=0;								//�ߓ_��
	int elemnum=0;								//�v�f��
	int BEM_elm_type=CON->get_BEM_elm_type();	//�v�f�^�C�v 0:��� 1:���`
	int B_flag=UNLOOP;							//B_flag��LOOP�Ȃ狫�E�͌q�����Ă���BUNLOOP�Ȃ�Ȃ����ĂȂ�

	if(CON->get_region_shape()==0) make_cube_region(CON,NODE,s_node_num);//��͗̈�͒�����
	else cout<<"��͗̈�`�󂪖���"<<endl;

	int KTJ=*s_node_num;
	int KTE=12*KTJ;		//�ő�v�f���@3���������Ƃ߂�
	int nelm=0;			//���݂̗v�f��

	vector <point3D> NODE2;
	point3D NODE2_0;
	for(int i=0;i<KTJ+8+1;i++) NODE2.push_back(NODE2_0); //�ړ_�̍��W(+8���Ă���̂ͽ��߰�ޯ���̒��_��)
	for(int i=0;i<KTJ;i++)
	{
		for(int D=0;D<3;D++) NODE2[i+1].r[D]=NODE[i].r[D];//���W���R�s�[
		NODE2[i+1].material=AIR;
	}
	
	
	//�v�f�N���X�쐬
	vector <element3D> ELEM2;
	element3D ELEM2_0;
	for(int i=0;i<KTE;i++) ELEM2.push_back(ELEM2_0);

	/////////////�ߓ_���W�̐��K��
    double xmin=NODE2[1].r[A_X];
    double ymin=NODE2[1].r[A_Y];
    double zmin=NODE2[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    ///���W�̍ő�A�ŏ��l�����߂�
    for(int i=2;i<=*s_node_num;i++)
    {
        if(NODE2[i].r[A_X]<xmin) xmin=NODE2[i].r[A_X];
		else if(NODE2[i].r[A_X]>xmax) xmax=NODE2[i].r[A_X];
	
		if(NODE2[i].r[A_Y]<ymin) ymin=NODE2[i].r[A_Y];
		else if(NODE2[i].r[A_Y]>ymax) ymax=NODE2[i].r[A_Y];
	
		if(NODE2[i].r[A_Z]<zmin) zmin=NODE2[i].r[A_Z];
		else if(NODE2[i].r[A_Z]>zmax) zmax=NODE2[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X�������̐��@
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y�������̐��@
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z�������̐��@
    double rmax=rax;		///�ő吡�@
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //������else�ɂ�����_��

    ///���W�ϊ�
    double rrm=1.000000/rmax;///�������������������邱�ƂŁA���l�덷�����点��E�E�H
    for(int i=1;i<=*s_node_num;i++)
    {   //   A/B�Ƃ����v�Z�������Ƃ��A�`�̒l�ɂ���Ĕ�����1/B�Ƃ����{�����������Ă���̂ł͂Ȃ����ƍl���āA���̂悤�ȏ������ɂ��Ă���
        NODE2[i].r[A_X]=(NODE2[i].r[A_X]-xmin)*rrm;
		NODE2[i].r[A_Y]=(NODE2[i].r[A_Y]-ymin)*rrm;
		NODE2[i].r[A_Z]=(NODE2[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////

	int FINE_sw=OFF;
	delaun3D(CON,NODE2,ELEM2, KTJ, KTE, rax, ray, raz,s_node_num,&nelm, FINE_sw,rrm);
	/////ү����������

	///���b�V���������m�F
	double *val=new double[KTJ+1];
	for(int i=1;i<=*s_node_num;i++) val[i]=1;
	data_avs(*s_node_num,nelm,NODE2,ELEM2,KTJ,val,CON);
	delete [] val;

	//�\�ʗv�fELEM�𐶐�
	BEMelement3D ELEM0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=ELEM2[i].elm[j];
			if(jelm==0)					//�\�ʂȂ�
			{
				ELEM.push_back(ELEM0);
				int ia=ELEM2[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
				int ib=ELEM2[i].node[4-(j-1)/2*2];
				int ic=ELEM2[i].node[3-(j/2%2)*2];
				
				ia-=1;		//�f���[�j�����̃v���O�������ł́A�ߓ_�ԍ���1����n�܂邪�ABEM�ł�0����n�܂�̂ŁA������-1�����ߓ_�ԍ����A�Y������ߓ_�ԍ��ł���B
				ib-=1;
				ic-=1;

				ELEM[elemnum].node[0]=ia;
				ELEM[elemnum].node[1]=ib;
				ELEM[elemnum].node[2]=ic;
				double iaic[3];//ia��ic���޸�ِ����i�[
				double iaib[3];//ia��ib���޸�ِ����i�[
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
				}
				///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
				double S[3];//��L���޸�ِ����i�[
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				ELEM[elemnum].S=SS;
				////�ʐ�S�����Ƃ܂���
				for(int D=0;D<3;D++) ELEM[elemnum].direct[D]=S[D]/SS;//�O�����P�ʖ@���޸��
				for(int D=0;D<3;D++) ELEM[elemnum].r[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;	//�d�S
				ELEM[elemnum].map=0;			//������
				elemnum++;
			}
		}
	}////�v�f��񐶐�����

	//�v�f�̋��E��������эގ��쐬
	for(int i=0;i<elemnum;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];
		int BD=Diric;			//���E����
		for(int j=0;j<3;j++) if(NODE[N[j]].boundary_condition==Neumn) BD=Neumn;//�ЂƂł��m�C�}���^�̐ߓ_���܂�ł���΁A���̖ʂ̓m�C�}���^
		ELEM[i].boundary_condition=BD;
		ELEM[i].material=AIR;			//�����ō쐬�����\�ʗv�f�͂��ׂ�AIR   ???
	}

	*s_elem_num=elemnum;
}

//�����̉�͗̈�쐬�֐�
void make_cube_region(mpsconfig *CON,vector<BEMpoint3D> &NODE,int *s_node_num)
{
	int node_num=0;
	int divN[3];								//�e�ӂ̕�����
	double divL[3];								//������
	double V1=0;//CON->get_V();
	double V2=CON->get_V();						//�|�e���V�����l

	BEMpoint3D NODE01;

	double Xmin=CON->get_minX();				//��͗̈�
	double Xmax=CON->get_maxX();
	double Ymin=CON->get_minY();
	double Ymax=CON->get_maxY();
	double Zmin=CON->get_minZ();
	double Zmax=CON->get_maxZ();

	
	divN[A_X]=10;
	divN[A_Y]=10;
	divN[A_Z]=10;
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//��͗̈�̒��
			NODE[node_num].potential=V1;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}
	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//��͗̈�̏��
			NODE[node_num].potential=V2;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}

	
	//����Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}

	//����
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//������
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
			node_num++;
		}
	}
	*s_node_num=node_num;
}

//���v�f�p�́A���I�ߓ_�E�v�f��񐶐��֐�
void set_BEM3D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number)
{
	int node_num=0;
	//int elem_num=0;
	double le=CON->get_distancebp();
	BEMpoint3D NODE01;
	
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
			for(int D=0;D<3;D++) NODE[node_num].r[D]=PART[i].r[D];
		//	NODE[node_num].boundary_condition=Diric;
			NODE[node_num].boundary_condition=BOTH;		//�Q�}���̂Ƃ��͂�����
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//������
			NODE[node_num].slop2=0; 				//������
			NODE[node_num].C=0.5;					//���v�f��0.5
			NODE[node_num].particle=i;				//�Ή����闱�q�ԍ��L��
			BEMID[i]=node_num;							//BEM�ɏo�͂��邵�邵
			node_num++;
		}
	}///

	///�ꎞ�I�ɁA�����������ɂ��ߓ_���쐬���Ă���
	for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)
		{
			NODE.push_back(NODE01);
			for(int D=0;D<3;D++) NODE[node_num].r[D]=PART[i].r[D]+direct[D][i]*le;//���������̍��W�����
			NODE[node_num].boundary_condition=NOCALC;		//���Ƃŏ����Ƃ����Ӗ��ŁANOCALC
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//������
			NODE[node_num].slop2=0; 				//������
			NODE[node_num].C=0.5;					//���v�f��0.5
			NODE[node_num].particle=-1;				//�Ή����闱�q�͑��݂��Ȃ�
			node_num++;
		}
	}/////�ߓ_���쐬����*/

	//delaun3D�p�̗v�f�N���X�쐬
	int KTJ=node_num;
	int KTE=12*KTJ;		//�ő�v�f���@3���������Ƃ߂�
	int nelm=0;			//���݂̗v�f��
	vector <point3D> NODE2;
	point3D NODE2_0;
	for(int i=0;i<KTJ+8+1;i++) NODE2.push_back(NODE2_0); //�ړ_�̍��W(+8���Ă���̂ͽ��߰�ޯ���̒��_��)
	for(int i=0;i<KTJ;i++)
	{
		for(int D=0;D<3;D++) NODE2[i+1].r[D]=NODE[i].r[D];//���W���R�s�[
		NODE2[i+1].material=FLUID;
	}

	//�v�f�N���X�쐬
	vector <element3D> ELEM2;
	element3D ELEM2_0;
	for(int i=0;i<KTE;i++) ELEM2.push_back(ELEM2_0);

	/////////////�ߓ_���W�̐��K��
    double xmin=NODE2[1].r[A_X];
    double ymin=NODE2[1].r[A_Y];
    double zmin=NODE2[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    ///���W�̍ő�A�ŏ��l�����߂�
    for(int i=2;i<=node_num;i++)
    {
        if(NODE2[i].r[A_X]<xmin) xmin=NODE2[i].r[A_X];
		else if(NODE2[i].r[A_X]>xmax) xmax=NODE2[i].r[A_X];
	
		if(NODE2[i].r[A_Y]<ymin) ymin=NODE2[i].r[A_Y];
		else if(NODE2[i].r[A_Y]>ymax) ymax=NODE2[i].r[A_Y];
	
		if(NODE2[i].r[A_Z]<zmin) zmin=NODE2[i].r[A_Z];
		else if(NODE2[i].r[A_Z]>zmax) zmax=NODE2[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X�������̐��@
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y�������̐��@
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z�������̐��@
    double rmax=rax;		///�ő吡�@
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //������else�ɂ�����_��

    ///���W�ϊ�
    double rrm=1.000000/rmax;///�������������������邱�ƂŁA���l�덷�����点��E�E�H
    for(int i=1;i<=node_num;i++)
    {   //   A/B�Ƃ����v�Z�������Ƃ��A�`�̒l�ɂ���Ĕ�����1/B�Ƃ����{�����������Ă���̂ł͂Ȃ����ƍl���āA���̂悤�ȏ������ɂ��Ă���
        NODE2[i].r[A_X]=(NODE2[i].r[A_X]-xmin)*rrm;
		NODE2[i].r[A_Y]=(NODE2[i].r[A_Y]-ymin)*rrm;
		NODE2[i].r[A_Z]=(NODE2[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////

	int FINE_sw=OFF;
	delaun3D(CON,NODE2,ELEM2, KTJ, KTE, rax, ray, raz,&node_num,&nelm, FINE_sw,rrm);
	/////ү����������

	for(int i=1;i<=nelm;i++) ELEM2[i].material=FLUID;			//�����ō쐬�����v�f�͂��ׂ�FLUID

	//�����ŁA�\�ʂ����ȂƂ���ɕs�v�ȃ��b�V���������Ă���ꍇ�́ANOCALC���܂܂Ȃ��v�f�̂Ȃ��ŁA�ӂ̒����������������v�f�̍ގ�����C�ɂ���Ƃ���������ǉ����邱��

	///���b�V���������m�F
	double *val=new double[KTJ+1];
	for(int i=1;i<=node_num;i++) val[i]=1;
	data_avs(node_num,nelm,NODE2,ELEM2,KTJ,val,CON);
	delete [] val;

	//�\�ʗv�fELEM�𐶐�
	BEMelement3D ELEM0;
	int elemnum=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=ELEM2[i].elm[j];
			int flag=OFF;
			if(jelm==0) flag=ON;
			//else if(ELEM2[jelm].material==AIR) flag=ON;
			if(flag==ON)					//�\�ʂȂ�
			{
				ELEM.push_back(ELEM0);
				int ia=ELEM2[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
				int ib=ELEM2[i].node[4-(j-1)/2*2];
				int ic=ELEM2[i].node[3-(j/2%2)*2];
				
				ia-=1;		//�f���[�j�����̃v���O�������ł́A�ߓ_�ԍ���1����n�܂邪�ABEM�ł�0����n�܂�̂ŁA������-1�����ߓ_�ԍ����A�Y������ߓ_�ԍ��ł���B
				ib-=1;
				ic-=1;

				ELEM[elemnum].node[0]=ia;
				ELEM[elemnum].node[1]=ib;
				ELEM[elemnum].node[2]=ic;
				double iaic[3];//ia��ic���޸�ِ����i�[
				double iaib[3];//ia��ib���޸�ِ����i�[
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
				}
				///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
				double S[3];//��L���޸�ِ����i�[
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				ELEM[elemnum].S=SS;
				////�ʐ�S�����Ƃ܂���
				for(int D=0;D<3;D++) ELEM[elemnum].direct[D]=S[D]/SS;//�O�����P�ʖ@���޸��
				for(int D=0;D<3;D++) ELEM[elemnum].r[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;	//�d�S
				ELEM[elemnum].map=0;			//������
				elemnum++;
				//cout<<elemnum<<endl;
			}
		}
	}////�v�f��񐶐�����

	//NOCALC�ߓ_�̏���
	int erasenum=node_num/2;//�����̐ߓ_��NOCALC������A���̐��������������s
	for(int i=0;i<erasenum;i++) NODE.pop_back();
	node_num-=erasenum;

	//�v�f�̋��E��������эގ��쐬
	for(int i=0;i<elemnum;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];
		int BD=NODE[N[0]].boundary_condition;			//���܂̂Ƃ���A���E������N[0]��N[1]�������Ɖ��肵�āAN[0]�̂��g�p
		ELEM[i].boundary_condition=BD;
		ELEM[i].material=FLUID;			//�����ō쐬�����\�ʗv�f�͂��ׂ�FLUID  
	}

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	cout<<"���I�ߓ_��="<<node_num<<" ���I�v�f��="<<elemnum<<endl;
	*d_node_num=node_num;
	*d_elem_num=elemnum;

	delete [] BEMID;
	
}


void couple3D_NODE_and_ELEM(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,vector<BEMpoint3D> &dy_NODE,vector<BEMelement3D> &dy_ELEM,vector<REGION> &region)
{
	//s_NODE��dy_NODE�����̂�����NODE�Ɋi�[����
	int s_node_num=int (s_NODE.size());
	int dy_node_num=int (dy_NODE.size());
	int s_elem_num=int (s_ELEM.size());
	int dy_elem_num=int (dy_ELEM.size());
	int countN=0;
	int countE=0;

	BEMpoint3D NODE01;
	BEMelement3D ELEM01;
	REGION temp;
	region.push_back(temp);
	region[0].start=0;					//�ŏ��̌v�Z�̈�͕K���v�Z�_�O����n�܂�

	int *s_Nid=new int[s_node_num];		//�ړ���̐ߓ_�ԍ��i�[�Bs_NODE[i]�������ߓ_��NODE[s_Nid[i]]�Ɋi�[�����
	int *dy_Nid=new int[dy_node_num];


	//�ߓ_��񐶐�
	for(int i=0;i<s_node_num;i++)
	{
		NODE.push_back(NODE01);
		for(int D=0;D<3;D++) NODE[countN].r[D]=s_NODE[i].r[D];
		NODE[countN].boundary_condition=s_NODE[i].boundary_condition;
		NODE[countN].potential=s_NODE[i].potential;
		NODE[countN].slop1=s_NODE[i].slop1;
		NODE[countN].slop2=s_NODE[i].slop2;
		NODE[countN].C=s_NODE[i].C;
		NODE[countN].particle=s_NODE[i].particle;
		s_Nid[i]=countN;			//i�Ԗڂ�s_NODE��countN�Ԗڂ̐ߓ_
		countN++;
	}
	//region.push_back(temp);//��Q�̈�
	//region[1].start=countN;
	for(int i=0;i<dy_node_num;i++)
	{
		NODE.push_back(NODE01);
		for(int D=0;D<3;D++) NODE[countN].r[D]=dy_NODE[i].r[D];
		NODE[countN].boundary_condition=dy_NODE[i].boundary_condition;
		NODE[countN].potential=dy_NODE[i].potential;
		NODE[countN].slop1=dy_NODE[i].slop1;
		NODE[countN].slop2=dy_NODE[i].slop2;
		NODE[countN].C=dy_NODE[i].C;
		NODE[countN].particle=dy_NODE[i].particle;
		dy_Nid[i]=countN;			//i�Ԗڂ�dy_NODE��countN�Ԗڂ̐ߓ_
		countN++;
	}
	//region[0].end=countN;
	//region[1].end=countN;
	

	//�v�f���쐬
	for(int i=0;i<s_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		for(int j=0;j<3;j++) ELEM[countE].node[j]=s_ELEM[i].node[j];
		for(int D=0;D<3;D++) ELEM[countE].r[D]=s_ELEM[i].r[D];
		ELEM[countE].boundary_condition=s_ELEM[i].boundary_condition;
		ELEM[countE].S=s_ELEM[i].S;
		for(int D=0;D<3;D++) ELEM[countE].direct[D]=s_ELEM[i].direct[D];
		ELEM[countE].map=s_ELEM[i].map;
		ELEM[countE].material=s_ELEM[i].material;
		countE++;
	}
	for(int i=0;i<countE;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];	//���̒i�K�Ŋi�[����Ă���ߓ_�ԍ��́As_NODE�̂��̂ł���ANODE�̂��̂ł͂Ȃ�
		for(int j=0;j<3;j++) N[j]=s_Nid[N[j]];		//�ߓ_�ԍ���NODE��̂��̂ɏ�������
		for(int j=0;j<3;j++) ELEM[i].node[j]=N[j];
	}
	region.push_back(temp);//��Q�̈�
	region[1].start=countE;
	for(int i=0;i<dy_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		for(int j=0;j<3;j++) ELEM[countE].node[j]=dy_ELEM[i].node[j];
		for(int D=0;D<3;D++) ELEM[countE].r[D]=dy_ELEM[i].r[D];
		ELEM[countE].boundary_condition=dy_ELEM[i].boundary_condition;
		ELEM[countE].S=dy_ELEM[i].S;
		for(int D=0;D<3;D++) ELEM[countE].direct[D]=dy_ELEM[i].direct[D];
		ELEM[countE].map=dy_ELEM[i].map;
		ELEM[countE].material=dy_ELEM[i].material;
		countE++;
	}
	for(int i=s_elem_num;i<countE;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];	//���̒i�K�Ŋi�[����Ă���ߓ_�ԍ��́Ady_NODE�̂��̂ł���ANODE�̂��̂ł͂Ȃ�
		for(int j=0;j<3;j++) N[j]=dy_Nid[N[j]];		//�ߓ_�ԍ���NODE��̂��̂ɏ�������
		for(int j=0;j<3;j++) ELEM[i].node[j]=N[j];
	}
	region[0].end=countE;
	region[1].end=countE;
	cout<<"�S�ߓ_��="<<countN<<" �S�v�f��="<<countE<<endl;

	//�t�@�C���o�͂��č��W���m�F
	ofstream fp("BEM_input.dat");
	double times=0.0001;
	for(int i=0;i<countE;i++)
	{
		int n=ELEM[i].node[0];
		fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<ELEM[i].r[A_Z]<<endl;
	}//*/
	//for(int i=0;i<countN;i++) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Y]<<"\t"<<NODE[i].potential<<endl;
	fp.close();


	delete [] s_Nid;	
	delete [] dy_Nid;
}