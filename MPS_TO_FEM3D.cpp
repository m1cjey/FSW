#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[
#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
#include"FEM3Dclass.h"
#include<omp.h>
#include<vector>
#include"function.h"

//�����̉�͗̈�쐬�֐�
void make_cube_region(mpsconfig *CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2]);

//�t�H
void MPSTOFEM3D_droplet(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number);

//MPS_TO_FEM3Dmain�֐�
void MPS_TO_FEM3Dmain(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number)
{
	int model=CON->get_model_number();
	if(model==3) MPSTOFEM3D_droplet(CON,node_num,NODE,PART,  fluid_number,  particle_number);//�t�H

	ofstream fp("node_c.dat");
	for(int i=1;i<=*node_num;i++) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	fp.close();
	cout<<"MPSTOFEM3D���� �ߓ_��="<<*node_num<<endl;
}

//�t�H
void MPSTOFEM3D_droplet(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number)
{
	int num=0;						//�ߓ_��
	double le=CON->get_distancebp();
    double err=1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODE�͐ߓ_�ԍ�1����X�^�[�g���邩��A�����łЂƂm�ۂ��Ă���

	////���̗��q�o��
    for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//��͒��AFEM�ɓn�����q�ԍ��𓝈ꂵ����
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//�ߓ_i�ɑΉ����闱�q�ԍ���i
			NODE[num].remesh=ON;			//�����b�V��ON
		}
    }////////////*/

	if(CON->get_region_shape()==0)
	{
		int divN[3];
		divN[A_X]=10;
		divN[A_Y]=10;
		divN[A_Z]=10;
		double regionX[2]={CON->get_XL(),CON->get_XR()};
		double regionY[2]={CON->get_YD(),CON->get_YU()};
		double regionZ[2]={CON->get_ZD(),CON->get_ZU()};
		int count=num;//�����_�ł̐ߓ_��
		
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//��͗̈�͒�����

		//���ɁA���܍쐬������͗̈�̏��������ɐߓ_��ǉ�����B�����FINE3D�ŋ��E�����ʂɂ悯���Ȑߓ_���ǉ�����邱�Ƃ�h�����߂ł���B
		double dX=(regionX[1]-regionX[0])/divN[A_X];				//��قǍ쐬������͗̈�̃��b�V������dX
		double dY=(regionY[1]-regionY[0])/divN[A_Y];				//��قǍ쐬������͗̈�̃��b�V������dY
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//��قǍ쐬������͗̈�̃��b�V������dZ
		//�ȉ��̂悤�ɒ�`�������@�̔��̈���쐬����΁A��͗̈�ƐV�������Ƃ̌��Ԃɂ͗ǎ��Ȑ��l�ʑ̂ɋ߂����b�V�����쐬�����
		regionX[0]=CON->get_XL()+dX*sqrt(3.0)*0.5; regionX[1]=CON->get_XR()-dX*sqrt(3.0)*0.5;
		regionY[0]=CON->get_YD()+dY*sqrt(3.0)*0.5; regionY[1]=CON->get_YU()-dY*sqrt(3.0)*0.5;
		regionZ[0]=CON->get_ZD()+dZ*sqrt(3.0)*0.5; regionZ[1]=CON->get_ZU()-dZ*sqrt(3.0)*0.5;

		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//�����̔��쐬

		for(int i=count+1;i<=num;i++)				//���E����
		{
			double Z=NODE[i].r[A_Z];
			if(Z>CON->get_ZU()-err) NODE[i].boundary_condition=2;
			else if(Z<CON->get_ZD()+err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh�̈�쐬
		regionX[0]=CON->get_XL()*0.5; regionX[1]=CON->get_XR()*0.5;
		regionY[0]=CON->get_YD()*0.5; regionY[1]=CON->get_YU()*0.5;
		regionZ[0]=CON->get_ZD()*0.7; regionZ[1]=CON->get_ZU()*0.7;
		count=num;//�����_�ł̐ߓ_��
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//��͗̈�͒�����

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//���E����
			NODE[i].remesh=ON;
		}
		
	}

	

	*node_num=num;

}

//�����̉�͗̈�쐬�֐�
void make_cube_region(mpsconfig *CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2])
{
	int node_num=*node;
	//divN[3];								//�e�ӂ̕�����
	
	double Xmin=regionX[0];				//��͗̈�
	double Xmax=regionX[1];
	double Ymin=regionY[0];
	double Ymax=regionY[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];								//������
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//��͗̈�̒��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	//���
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//��͗̈�̏��
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	
	//����Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}

	//����
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//�Ή����闱�q�����݂��Ȃ�����-1���i�[
		}
	}
	*node=node_num;
}
