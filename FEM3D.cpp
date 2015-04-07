#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"	//#define �i�[
//#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
//#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
//#include"FEM3Dclass.h"
//#include<omp.h>
//#include<vector>
#include<complex>

#include"function.h"

#define FULL 3
#define REMESH 4
#define FULL_INPORT 5

int make_edge_element(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,vector<edge3D> &EDGE,int *branch_num,int **nei2,int KTE,vector<edge3D> &static_EDGE,int t,int node_sta);
void calc_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,int *jnb,int **nei,int *branch_num,double **current);
void denryu_side(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *T,double **current);
void VOLT3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei,double *RP);
void potential_calculation(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei);
void inport_J0_density(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0Je(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void Avector3D_node_eddy2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t);
void Avector3D_node_eddy2_jw(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_node_eddy2_jw2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_node_eddy2_jw3(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw);
void Avector3D_node_eddy2_jw5(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw);
void Avector3D_node_eddy2_jw5_ver2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta);
void Avector3D_node_eddy2_jw4(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP);
void Avector3D_edge_eddy(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_edge_eddy_jw(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_edge_eddy_jw2(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,int node_sta,vector<edge3D> &static_EDG,double *Am,double *phi);
void Avector3D_edge_eddy_jw_with_parabolic_node_element(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void calc_transitional_EM_field(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,int *jnb,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int **nei,int particle_node,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta,vector<edge3D> &static_EDGE);
void Bflux3D_node(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B,int t,int flag);
void Bflux3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *A,double **B,int t,int flag);
///���E�����K�p�֐�(�R�c�ӗv�f�p)
void set_boundary_condition3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A);
void NODE_F3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpsparticle> &PART,double **F,int fluid_number,int t);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);
void calc_eddy_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **old_A,double dt,double *V,double **Je,int t,double *sigma);
void calc_node_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **AR,double **AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t, double TIME,double *sigma, double omega);
void calc_edge_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,double *AR,double *AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t, double TIME,double *sigma, double omega);
int poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void laplacian_smoothing(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int node0,int nelm0);
void laplacian_smoothing2(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int *jnb, int **nei,int node0);
void node_sorting(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei);
void node_sorting2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei);
void arrange_matrix_complex(int pn,int *NUM,int **ROW,complex<double>**G);
void check_matrix_symmetry_complex(int pn,int *NUM,int **ROW,complex<double> **G);
void diagonal_distribution(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int pn,int *NUM,int **ROW,double **G);
void COCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void parallel_ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void DS_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void parallel_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);

void modify_node_info(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int**nei,int *newID);
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void calc_jw_field_node(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm);
void calc_jw_field_edge(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm,int nedge,double *Am, double *phi);



int FEM3D_calculation(mpsconfig *CON,int *static_node,int *static_nelm,int *static_nedge,vector<point3D> &static_NODE, vector<element3D> &static_ELEM,vector<edge3D> &static_EDGE,int particle_node,double **F,int t,double TIME,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,mpsconfig& CONF)
{
	/*//////////////////

	//COCG���܂Ƃ����ǂ������m�F���邽�߁A�ȒP�ȕ��f���̕�����������������
	 
	cout<<"test.cocg�@�X�^�[�g"<<endl;

	int pn3=2;
	int number3=4;
	complex<double> Re=(2,5.8);
	complex<double> Im;
	complex<double> z;
	double x=1; double y=2;
	complex<double> *val3 = new complex<double> [number3];
    int *ind3 = new int [number3];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr3 = new int [pn3+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	complex<double> *B3=new complex<double> [pn3];//���s��
	complex<double> *XX3=new complex<double> [pn3];//�s��̓����i�[

	double aaa=0;
	aaa=cos(2*PI*1000000);
	cout<<"aaa="<<aaa<<endl;
	Im=complex<double> (0,1);
	cout<<"Re="<<Re<<endl;
	Im*=2;
	cout<<"Im="<<Im<<endl;
	y*=4;
	z=complex<double> (1,1);
	complex<double> xx;
	xx=complex<double> (2,0);
	z*=xx;
	//z+=complex<double> (0,y);
	//////
	
	cout<<"z="<<z<<endl;

	//////
	//3x+(1+i)y=-i;  (1+i)x+(2-i)y=4;  =>   x=-i; y=1+i; 
	val3[0]=complex<double> (3.0,0.0);
	val3[1]=complex<double> (1.0,1.0);
	val3[2]=complex<double> (1.0,1.0);
	val3[3]=complex<double> (2.0,-1.0);

	ind3[0]=0;
	ind3[1]=1;
	ind3[2]=0;
	ind3[3]=1;

	ptr3[0]=0;
	ptr3[1]=2;
	ptr3[2]=4;

	B3[0]=complex<double> (0.0,-1.0);
	B3[1]=complex<double> (4.0,0.0);
	////////

	//COCG(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	//cs_MRTR(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	//ICCOCG(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	cs_ICMRTR(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);

	for(int i=0;i<pn3;i++) cout<<"xx="<<XX3[i]<<endl;

	delete [] val3;
    delete [] ind3;
    delete [] ptr3;
	delete [] B3;
	delete [] XX3;

	///////*/	

	//j�֖@�Ōv�Z������̃X�e�b�v�̊Ԃ́A���g�������𗘗p���ĕK�v�ȕ����l�����߁A�֐��𔲂��� //
	
	if(CON->get_m_A()==1)
	{
		if(CON->get_jw_Faverage()==OFF)//���̌v�Z�̎��ԍ��ݕ����ƂɁA���g�����������Ƃɓd���͂��v�Z����//���ݔp�~
		{   
			/*
			int dh=(int) (1.0/(CON->get_dt()*CON->get_Hz()));//�P����������������Ă��邩�Bdt��1/f��dh�����ł���悤�ɐݒ肷�邱��
			if(t%(CON->get_jw_interval()*dh)!=1)
			{
				cout<<"���g����������A���v�Z"<<endl;
				calc_jw_field(CON,NODE_jw,ELEM_jw,dt,TIME,PART,fluid_number,F,t,node,nelm);
				return 0;
			}
			*/
		}
		else if(CON->get_jw_Faverage()==ON)//�����⋭���d����1�����ɓ����d���͂̕��ϒl�𗘗p����
		{/*
			if(t%(CON->get_jw_interval())!=1)
			{
				//���ϒl�����߂鏈����FEM���s�����X�e�b�v�ɂĂ��łɏI���Ă���̂ŁA�O�̃X�e�b�v�̓d���͂�ǂݍ���
				cout<<"�d���͂̕��ϒl�̓ǂݍ���"<<endl;
				ifstream f("F_FEM.dat");
				if(!f) cout<<"cannot open F_FEM.dat"<<endl;
				f.unsetf(ifstream::dec);
				f.setf(ifstream::skipws);
				for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) f>>F[D][i];			
				f.close();
				return 0;
			}
		*/
		}
	}

	vector<point3D> NODE;
	vector<element3D> ELEM;
	vector<edge3D> EDGE;

	int node=0;					//�S�ߓ_��
	int nelm=0;					//�S�v�f��
	int nedge=0;				//�S�Ӑ�
	int KTJ;					//�ő�ߓ_��
	int KTE;					//�ő�v�f���@3���������Ƃ߂�
	double err=1.0e-14;			//�덷����̂������l
	

	
	if(CON->get_mesher()==0) //magnet�ō�������b�V�������Ƃɂ�ړ����̃��b�V������
	{
		int delaun_flag;			//�f���[�j�������s�����A�s��Ȃ���
		int node0=0;

		if(CON->get_mesh_input()==0)							//MPSTOFEM�ɂ��ߓ_�A�v�f����
		{
			if(CON->get_remesh_sw()==OFF) delaun_flag=FULL;		//remesh�̈��z�肹���A��ɂ��ׂĂ��f���[�j����
			else if(CON->get_remesh_sw()==ON)
			{
				if(t==1) delaun_flag=FULL;	//�S���f�����f���[�j����
				else delaun_flag=REMESH;	//remesh�̈�̂݃f���[�j����
			}
		}
		else if(CON->get_mesh_input()==1)						//Magnet���ǂݍ���
		{
			if(CON->get_remesh_sw()==OFF) delaun_flag=FULL_INPORT;		//���Magnet�̗v�f��ǂݍ��݉�� �Ǘ��җp�H
			else if(CON->get_remesh_sw()==ON)
			{
				if(t==1) delaun_flag=FULL_INPORT;	//�S���f����Magnet�t�@�C�����ǂݍ���
				else delaun_flag=REMESH;	//remesh�̈�̂݃f���[�j����
			}
		}

		if(delaun_flag==FULL) cout<<"FULL �f���[�j�������s"<<endl;
		else if(delaun_flag==REMESH) cout<<"remesh�̈�̂݃f���[�j�������s"<<endl;
		else if(delaun_flag==FULL_INPORT) cout<<"Magnet�����t�@�C�����v�f��񓙓ǂݍ���"<<endl;

		if(delaun_flag==FULL)
		{
			MPS_TO_FEM3Dmain(CON,&node,NODE,PART,  fluid_number,  particle_number);//���q�z�u���ߓ_�z�u�����
			KTJ=node;	
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;
		}
		else if(delaun_flag==REMESH)
		{
			node=(int) static_NODE.size()-1;	//�ÓI�ߓ_��
			nelm=(int) static_ELEM.size()-1;
			KTJ=node+particle_number;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;
		}
		else if(delaun_flag==FULL_INPORT)
		{
			ifstream fin("input_from_MAGNET.dat");
			if(!fin) cout<<"cannot open input_from_MAGNET.dat"<<endl;
			fin.unsetf(ifstream::dec);
			fin.setf(ifstream::skipws);

			fin>>node;			//�ߓ_���ǂݍ���
			fin>>nelm;			//�v�f���ǂݍ���
			KTJ=node+particle_number;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;

			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			int ID;

			for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//�z����m��
			for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);	//�z����m��

			//�ߓ_���ǂݍ���
			for(int i=1;i<=node;i++)
			{
				fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
				for(int D=0;D<3;D++) fin>>NODE[i].r[D];
				fin>>NODE[i].material;
				if(NODE[i].material==21) NODE[i].material=CRUCIBLE;
				NODE[i].boundary_condition=0;			//�Ƃ肠�����[�����i�[
				NODE[i].particleID=-1;					//�Ή����闱�q�͑��݂��Ȃ�
				NODE[i].remesh=OFF;						//non-remesh�̈�̐ߓ_�ł���Bremesh�̈�Ƃ̋��E�Ɉʒu����ߓ_�Ɋւ��Ă͌�ɏ������{��
				NODE[i].BD_node=OFF;	
			}
			//�v�f-�ߓ_���ǂݍ���
			for(int i=1;i<=nelm;i++)
			{
				fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
				for(int j=1;j<=4;j++) fin>>ELEM[i].node[j];
			}
			//�v�f-�v�f���ǂݍ���
			for(int i=1;i<=nelm;i++)
			{
				fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
				for(int j=1;j<=4;j++) fin>>ELEM[i].elm[j];
			}
			//�v�f�ގ����ǂݍ���
			for(int i=1;i<=nelm;i++)
			{
			fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
				fin>>ELEM[i].material;
				//////
				if(ELEM[i].material==21) ELEM[i].material=CRUCIBLE;
				/////
				ELEM[i].map=0;			//������
				for(int D=0;D<3;D++)ELEM[i].r[D]=0;
				ELEM[i].RR=0;
				ELEM[i].volume=0;		
			}

			fin.close();

			//���E�����ݒ�
			int *BD_flag=new int[node+1];			//BD_flag=ON�Ȃ狫�E�ߓ_(��͋��E��������Ȃ����Aremesh���E��������Ȃ�)
			for(int i=0;i<=node;i++) BD_flag[i]=OFF;//������
			for(int i=1;i<=nelm;i++)
			{
				for(int j=1;j<=4;j++)
				{
					if(ELEM[i].elm[j]==0)
					{
						int ia=ELEM[i].node[j%4+1];
						int ib=ELEM[i].node[4-(j-1)/2*2];
						int ic=ELEM[i].node[3-(j/2%2)*2];
						BD_flag[ia]=ON;						//���E�ߓ_�Ƃ�����
						BD_flag[ib]=ON;
						BD_flag[ic]=ON;
					}
				}
			}
			vector <int> BD_NODE_ID;						//���E�ߓ_�ԍ��i�[
			for(int i=1;i<=node;i++)
			{
				if(BD_flag[i]==ON) BD_NODE_ID.push_back(i);
			}
			int BD_num=(int) BD_NODE_ID.size();
			ofstream fs("remesh.dat");
			if(CON->get_region_shape()==0)				//��͗̈悪�����̂Ȃ�
			{
				double Xmax=0;
				double Ymax=0;
				double Zmax=0;

				

				for(int i=0;i<BD_num;i++)
				{
					int n=BD_NODE_ID[i];	//���E�ߓ_�ԍ�
					if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z��
					else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z��
					else if(NODE[n].r[A_X]<CON->get_XL()+err) NODE[n].boundary_condition=1;		//-X��
					else if(NODE[n].r[A_X]>CON->get_XR()-err) NODE[n].boundary_condition=1;		//+X��
					else if(NODE[n].r[A_Y]<CON->get_YD()+err) NODE[n].boundary_condition=1;		//-Y��
					else if(NODE[n].r[A_Y]>CON->get_YU()-err) NODE[n].boundary_condition=1;		//+Y��
					else
					{
						NODE[n].boundary_condition=0;					//��͋��E�ߓ_�ł͂Ȃ��Aremesh�̈�Ƃ̋��E�ߓ_�Ȃ̂ŁA���E�����̓[��
						NODE[n].remesh=ON;
						NODE[n].BD_node=ON;							//remesh�̈�̋��E���ł���Ƃ������邵
						fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
					}
					if(NODE[n].r[A_Z]>Zmax) Zmax=NODE[n].r[A_Z];
					if(NODE[n].r[A_Y]>Ymax) Ymax=NODE[n].r[A_Y];
					if(NODE[n].r[A_X]>Xmax) Xmax=NODE[n].r[A_X];
				}
				cout<<Xmax<<" "<<Ymax<<" "<<Zmax<<endl;
			}
			else if(CON->get_region_shape()==1)				//��͗̈悪�~���Ȃ�
			{
				for(int i=0;i<BD_num;i++)
				{
					int n=BD_NODE_ID[i];	//���E�ߓ_�ԍ�
					double R=sqrt(NODE[n].r[A_X]*NODE[n].r[A_X]+NODE[n].r[A_Y]*NODE[n].r[A_Y]);
					if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z��
					else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z��
					else if(R>CON->get_RU()*0.99) NODE[n].boundary_condition=1;	
					else
					{
						NODE[n].boundary_condition=0;											//��͋��E�ߓ_�ł͂Ȃ��Aremesh�̈�Ƃ̋��E�ߓ_�Ȃ̂ŁA���E�����̓[��
						NODE[n].remesh=ON;
						NODE[n].BD_node=ON;							//remesh�̈�̋��E���ł���Ƃ������邵
						fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
					}
				}
			}
			delete [] BD_flag;
			fs.close();

			//NODE,ELEM�̂����A�����Ȃ��v�f�A�ߓ_������static_NODE,staticELEM�Ɋi�[����
			static_ELEM.clear();
			static_NODE.clear();
			//static_EDGE.clear();
			memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);
			
			if(CON->get_remesh_sw()==ON)
			{
				node=(int) static_NODE.size()-1;
				nelm=(int) static_ELEM.size()-1;
				delaun_flag=REMESH;				//delaun_flag��REMESH�ɂ��邱�ƂŁA����if���ɓ����ē��I�v�f�𐶐�����
				NODE.clear();
				ELEM.clear();
			}
			else cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;

		}
		

		if(delaun_flag==FULL)									//���ׂĂ��f���[�j����
		{
			point3D NODE0;
			element3D ELEM0;
			for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//�z����m��
			int FINE_sw=CON->get_fine();					//�����ĕ����X�C�b�`
			delaun3D_main(CON,NODE,ELEM, KTJ, KTE,&node,&nelm, FINE_sw);

			cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;
			///���b�V���������m�F
			double *val=new double[KTJ+1];
			for(int i=1;i<=node;i++) val[i]=NODE[i].material;
			//data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,val,t);//�f�ʐ}
			data_avs3(node,nelm,NODE,ELEM,CON,t);//�ގ�

			if(CON->get_remesh_sw()==ON)				//remesh�̈��z�肷��Ȃ�A�ÓI�ߓ_�E�v�f�����L�����Ă���
			{
				//NODE,ELEM�̂����A�����Ȃ��v�f�A�ߓ_������static_NODE,staticELEM�Ɋi�[����
				static_ELEM.clear();
				static_NODE.clear();
				memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);

				/*/�`�F�b�N
				int snode=(int) static_NODE.size()-1;
				int snelm=(int) static_ELEM.size()-1;
				for(int i=1;i<=snode;i++) val[i]=static_NODE[i].remesh;
				data_avs2(CON,snode,snelm,static_NODE,static_ELEM,KTJ,val);//�f�ʐ}*/
			}

			delete [] val;
		}
		else if(delaun_flag==REMESH)									//remesh�̈���f���[�j����
		{
			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);			//�ߓ_�z��m�� +8�̓X�[�p�[�{�b�N�X
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);				//�z����m��
			for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);				//�z����m��

			for(int i=1;i<=node;i++)										//���̎��_��node�ɂ͐ÓI�ߓ_�����i�[����Ă���
			{
				for(int D=0;D<3;D++) NODE[i].r[D]=static_NODE[i].r[D];
				NODE[i].boundary_condition=static_NODE[i].boundary_condition;
				NODE[i].material=static_NODE[i].material;
				NODE[i].particleID=static_NODE[i].particleID;
				NODE[i].remesh=static_NODE[i].remesh;
				NODE[i].BD_node=static_NODE[i].BD_node;
			}

			//static_ELEM�������copy
			for(int i=1;i<=nelm;i++)										//���̎��_��nelm�ɂ͐ÓI�v�f�����i�[����Ă���
			{
				for(int D=0;D<3;D++) ELEM[i].r[D]=static_ELEM[i].r[D];
				for(int j=1;j<=4;j++)
				{
					ELEM[i].node[j]=static_ELEM[i].node[j];
					ELEM[i].elm[j]=static_ELEM[i].elm[j];
					//if(ELEM[i].elm[j]==0) cout<<i<<endl;
				}
				ELEM[i].map=static_ELEM[i].map;
				ELEM[i].material=static_ELEM[i].material;
				ELEM[i].RR=static_ELEM[i].RR;
				ELEM[i].volume=static_ELEM[i].volume;
				//�ӂ́H�H�H
			}

			//�ÓI�v�f��remesh�̈�ɐڂ���Ƃ��낪ELEM[i].elm=0�ƂȂ��Ă���.����ȗv�f��T��
			vector<int> BD_static_ELEM;		//���I�v�f�ɐڂ���ÓI�v�f�ԍ��i�[
			vector<int> BD_static_ELEM_elm;	//�ÓI�v�f�����I�v�f�ɐڂ���ʔԍ��i�[
			for(int i=1;i<=nelm;i++)
			{
				int flag=OFF;
				int J;
				for(int j=1;j<=4;j++) if(NODE[ELEM[i].node[j]].remesh==ON) flag=ON; 
				if(flag==ON)
				{
					
					flag=OFF;
					int num=0;
					for(int j=1;j<=4;j++)
					{
						if(ELEM[i].elm[j]==0)
						{
							flag=ON;
							J=j;
							num++;
						}
					}
					//if(num>1) cout<<i<<" �ÓI�v�f�������ʂœ��I�v�f�Ɛڂ��Ă��܂� mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].boundary_condition<<NODE[ELEM[i].node[2]].boundary_condition<<NODE[ELEM[i].node[3]].boundary_condition<<NODE[ELEM[i].node[4]].boundary_condition<<endl;
					if(num>1) cout<<i<<" �ÓI�v�f�������ʂœ��I�v�f�Ɛڂ��Ă��܂� mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].r[A_Z]<<endl;
				}
				if(flag==ON)
				{
					BD_static_ELEM.push_back(i);
					BD_static_ELEM_elm.push_back(J);//�v�fi�͑�J�ʂœ��I�v�f�ɐڂ��Ă���
				}
			}
			cout<<"KK="<<BD_static_ELEM.size()<<endl;

			double Pn[3];				//�N�_�̍��W

			for(int k=0;k<1;k++)
			{
				int i=BD_static_ELEM[k];
				int j=BD_static_ELEM_elm[k];
				int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
				int ib=ELEM[i].node[4-(j-1)/2*2];
				int ic=ELEM[i].node[3-(j/2%2)*2];

				double iaic[3];//ia��ic���޸�ِ����i�[
				double iaib[3];//ia��ib���޸�ِ����i�[
				double iaicL=0;	//ia��ic�̒���
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					iaicL+=iaic[D]*iaic[D];
				}
				iaicL=sqrt(iaicL);
				///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
				double S[3];//��L���޸�ِ����i�[
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				////�ʐ�S�����Ƃ܂���

				double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�@���޸��

				double Gp[3];								//�\�ʎO�p�`�̏d�S
				for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

				for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//�N�_�̍��W

				//////////////////////
				///�N�_�������I�ɍ��W�w��
				
				Pn[A_X]=0;
				Pn[A_Y]=0;
				Pn[A_Z]=0.22;//0.10325;	//0.225;//	//��ڒꕔ�F0.10125

				//���̂̍ő卂�������߁A����̏�����ɐݒu

				/*//���̐ߓ_�̍ő卂���A�ŏ�����
				double Zmax=0;
				double Zmin=10;
				for(int i=1;i<=node;i++)
				{
					if(NODE[i].material==FLUID)
					{
						if(NODE[i].r[A_Z]>Zmax) Zmax=NODE[i].r[A_Z];
						if(NODE[i].r[A_Z]<Zmin) Zmin=NODE[i].r[A_Z];
					}
				}

				Pn[A_X]=0;
				Pn[A_Y]=0;
				Pn[A_Z]=Zmax-CON->get_distancebp();
				*///

				node++;
				for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
				NODE[node].boundary_condition=0;
				NODE[node].material=AIR;
				NODE[node].particleID=-1;
				NODE[node].remesh=ON;

			}
			cout<<"�N�_����"<<endl;
			int countBD=0;
			for(int i=1;i<=node;i++) if(NODE[i].BD_node==ON) countBD++;
			cout<<"���E�̐ߓ_��="<<countBD<<endl;

			int *imen[4];
			for(int D=0;D<4;D++) imen[D]=new int [BD_static_ELEM.size()+1];
			int *jmen =new int [BD_static_ELEM.size()+1];
			int *kmen =new int [BD_static_ELEM.size()+1];
			double *vol =new double [BD_static_ELEM.size()+1];
			int ip=node;

			/*
			if(BD_static_ELEM.size()>100000) cout<<"imen�Ȃǂ̃������������Ă�������"<<endl;
			int ip=node;			
			int imen[100000][3+1];	//���ʑ̕\�ʎO�p�`�̐ߓ_�ԍ��i�[
			int jmen[100000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̔ԍ��i�[
			int kmen[100000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̗̂אږʔԍ� (����͑扽�ʂŎ����Ɛڂ��Ă��邩)
			double vol[100000];		//���ʑ̂̑̐ς̂U�{
			*/

			for(int k=1;k<=BD_static_ELEM.size();k++)
			{
				int kelm=BD_static_ELEM[k-1];//remesh�̈�ɐڂ���ÓI�v�f
				int j=BD_static_ELEM_elm[k-1];
				int ia=ELEM[kelm].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
				int ib=ELEM[kelm].node[4-(j-1)/2*2];
				int ic=ELEM[kelm].node[3-(j/2%2)*2];

				/*
				imen[k][1]=ic;
				imen[k][2]=ib;
				imen[k][3]=ia;
				*/
				imen[1][k]=ic;
				imen[2][k]=ib;
				imen[3][k]=ia;

				jmen[k]=kelm;
				kmen[k]=j;//jelm��ielm�ɐڂ���ʔԍ�
				//vol[k]=volume3D(NODE,ia,ib,ic,ip);
				vol[k]=volume3D(NODE,ic,ib,ia,ip);	//�ߓ_�̏��Ԃɒ���
				//if(vol[k]<0) cout<<"�̐ϕ��̗v�f���� "<<vol[k]<<endl;
			}

			int ibound=(int) BD_static_ELEM.size();//�\�ʂ̐���\��
			int nelm0=nelm;//�ύX�O�̗v�f�����L��
			
			for(int i=1;i<=ibound;i++)//�v�f��񐶐�
			{   
				nelm++;
				int ielm=nelm0+i;
				int ia=imen[1][i];
				int ib=imen[2][i];
				int ic=imen[3][i];
				ELEM[ielm].node[1]=ia;
				ELEM[ielm].node[2]=ib;
				ELEM[ielm].node[3]=ic;
				ELEM[ielm].node[4]=ip;//�V�_�͂S�Ԗڂƒ�`	
				ELEM[ielm].elm[4]=jmen[i];
				if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
				ELEM[ielm].volume=vol[i];
				
				sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
			
				ELEM[ielm].material=AIR;//�ʏ��poly�֐��Ƃ̑���_�B�����ōގ�����C�ƌ��肷��
			}
			///////////////////

			//�v�f-�v�f�֌W�C��/////////��̏����ő�4�ʂŐڂ���v�f�ԍ��͂킩���Ă���̂ŁA�c������߂�
			//						�����ŁA1�`3�ʂ͑��ʑ̂��\������v�f�Ƃ̋��E�ʂł��邱�Ƃɒ���
			int ix=0;
			
			for(int i=1;i<=ibound;i++)
			{
				int ielm=nelm0+i;
				for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]�͂��łɂ��Ƃ܂�������A����ȊO�����Ƃ߂�
				{
					///ELEM[ielm].node[4]=ip�ł���
					int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3�̂Ƃ��A2,3,1�̏�
					int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3�̂Ƃ��A3,1,2�̏�
					int flag=0;
					for(int k=1;k<=ix;k++)
					{
						if(flag==0)
						{
							int ja=imen[1][k];
							int jb=imen[2][k];
							if(ia==ja && ib==jb)//�ߓ_����v������
							{
								ELEM[ielm].elm[j]=jmen[k];//���炩����ؽĂ��Ă����������i�[
								ELEM[jmen[k]].elm[kmen[k]]=ielm;
								imen[1][k]=imen[1][ix];		//k�Ԗڂ̏��͂����s�v�B�Ȃ̂Ŕz��̈�ԍŌ�̏���k�Ԗڂɂ����Ă��āA����܂ł̏��͔j������
								imen[2][k]=imen[2][ix];
								jmen[k]=jmen[ix];
								kmen[k]=kmen[ix];
								ix--;						//�҂��Ӑ�����
								flag=1;						//ELEM[ielm].elm[j]�͂��Ƃ܂����̂ŁA���̃l�X�g�ɓ���K�v�͂Ȃ��̂�flag=1
							}
						}
					}
					if(flag==0)
					{
						ix++;			//�����ł�ix�́A[�אڊ֌W�𖞂����v�f]���܂��Ă���[��]�̐���\���B
						imen[1][ix]=ib;	//�����̐ߓ_�̕��т��L�������A�ʂ̗v�f�����̕��т𖞂����̂�҂Bib��ia�̕��т��t�ɂ��Ă��邱�Ƃɒ���
						imen[2][ix]=ia;
						jmen[ix]=ielm;
						kmen[ix]=j;
					}
				}	
			}///�v�f-�v�f�֌W�C������

			cout<<"�v�f�쐬����"<<endl;

			fill3D(NODE,ELEM,nelm);

			node0=node;
			cout<<"node0="<<node0<<endl;

			/////////
			for(int D=0;D<4;D++) delete [] imen[D];
			delete [] jmen;
			delete [] kmen;
			delete [] vol;


			/*/���E�ʂ���w���`��
			for(int m=1;m<=1;m++)
			{
				for(int k=1;k<BD_static_ELEM.size();k++)
				{
					int i=BD_static_ELEM[k];
					int j=BD_static_ELEM_elm[k];
					int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];

					double iaic[3];//ia��ic���޸�ِ����i�[
					double iaib[3];//ia��ib���޸�ِ����i�[
					double iaicL=0;	//ia��ic�̒���
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
						iaicL+=iaic[D]*iaic[D];
					}
					iaicL=sqrt(iaicL);
					///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
					double S[3];//��L���޸�ِ����i�[
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					////�ʐ�S�����Ƃ܂���

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�@���޸��

					double Gp[3];								//�\�ʎO�p�`�̏d�S
					for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

					//for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//�N�_�̍��W �I���W�i��
					for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*0.0001*m;		//�N�_�̍��W

					node++;
					for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
					NODE[node].boundary_condition=0;
					NODE[node].material=AIR;
					NODE[node].particleID=-1;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;
				}
			}
			*/


			double u0=PI*4E-7;			//��C�̓�����
			double skin_depth=sqrt(1.0/(PI*CON->get_Hz()*CON->get_ele_conduc()*u0));//�\��[��
			cout<<"���̕\��[��="<<skin_depth<<endl;
			
			ofstream fp("rr.dat");
			/////���̗��q�𓮓I�ߓ_�Ƃ��Ċi�[
			
			if(CON->get_thinout_fluid()==0)
			{
				/////
				for(int i=0;i<fluid_number;i++)
				{
					//if(PART[i].surface==ON  || i%4==0)
					node++;
					for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
					NODE[node].boundary_condition=0;
					NODE[node].material=FLUID;
					NODE[node].particleID=i;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;
					fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				}
				
			}
			
			if(CON->get_thinout_fluid()>0)
			{
				double Xf=0.0;
				double Yf=0.0;
				double Zf=0.0;
				for(int i=0;i<fluid_number;i++)
				{
					Xf+=PART[i].r[A_X];
					Yf+=PART[i].r[A_Y];
					Zf+=PART[i].r[A_Z];
				}
				Xf/=fluid_number;
				Yf/=fluid_number;
				Zf/=fluid_number;
				cout<<"���̏d�S=("<<Xf<<","<<Yf<<","<<Zf<<")"<<endl;
				
				for(int i=0;i<fluid_number;i++)
				{
					double rx=PART[i].r[A_X]-Xf;
					double ry=PART[i].r[A_Y]-Yf;
					double rz=PART[i].r[A_Z]-Zf;
					if(PART[i].surface==ON || sqrt(rx*rx+ry*ry+rz*rz)>CON->get_fluidwidth()*0.001-1.1*skin_depth)
					{
						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=FLUID;
						NODE[node].particleID=i;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
						fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					}
					else if(i%CON->get_thinout_fluid()==0)
					{
						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=FLUID;
						NODE[node].particleID=i;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
						fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					}
				}
			}
			fp.close();
			/////*/



			/*///////
			//�\�ʂ̐��x����̂��߁A���̐ߓ_��ǉ�
			double Rz=0.13125-0.003;
			//�܂��͔��������B���̂��߂ɂ͔����\�ʂ��쐬����K�v������B
			for(int Ri=1;Ri<=5;Ri++)
			{
				double R=0.025-5*0.0001+Ri*0.0002;
				double le=CON->get_distancebp();

				double A=sqrt(3.0)/2;				//�悭�g���W��
				double B=sqrt(2.0/3);						////�悭�g���W��
				int half_WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
				int half_WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
				int half_WZ=(int)(R/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
				double R2=R-0.5*le;				//���������߂̔��a��ݒ�

				///////////�����\��
				int Nt;						//���\�ʂ́A�ƕ����̕�����
				double Lt;					//���\�ʂ́A�ƕ����̕�������
				//calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//���~�̕������͋����E��ǂ���ł��悢
				double temp_N=PI/2*R/(le*A);			//���̕������Ble�Ŋ���؂ꂽ���Ԃ������ǁA�����������Ȃ��Ƃ�������
				int Ns=(int) temp_N;				//�^�̕�����
				double difference=temp_N-Ns;		//���Ɛ^�̍�
				if(difference>0.5) Ns++;
				Lt=PI/2*R/Ns;			//���q�̋���
				Nt=Ns;

				double d_theta=Lt/R;		//�ʂ̒�����Lt�ɂȂ�p�x

				for(int k=0;k<=Nt;k++)//loop��k<Nt�ŏI��点��BNt�ɊY������Ƃ���͂��łɐݒu�ς�
				{
					double THETA=k*d_theta;	//��
					double r=R*sin(THETA);	//���̍����ɂ�����~�̔��a
					double round=2*PI*r;//���̍����ɂ�����~��

					//int Nf=calc_division_N_circle(round,le);//���\�ʂ́A�ƕ����̕�����
					//�Ώ̐����l�������A�~�����L�q���闱�q�͋����łȂ���΂Ȃ�Ȃ��B�����瑼�̕ӕ������Ƃ͈�������������
					//dis:�������鋗��(�~��)
					double temp_num=round/le;		//�~�O���ɐݒu����w���́x���q���B�������O�������܂�le�Ŋ���؂��Ƃ͌���Ȃ�

					int N1=(int)(temp_num/2);
					N1*=2;							
					int N2=N1+2;					//temp_num��N1��N2�̊Ԃɂ���B������N1,N2�͋���

					double dif1=temp_num-N1;		//�eN�Ƃ̍�
					double dif2=N2-temp_num;
					int N=N1;						//������������
					if(dif2<dif1) N=N2;				//���̏���������N�Ƃ��č̗p����B

					int Nf=N;
					
					
					double Lf=round/Nf;						//���\�ʂ́A�ƕ����̕�������
					double d_fai=Lf/r;						//�ʂ̒�����Lf�ɂȂ�p�x
					
					for(int i=0;i<Nf;i++)
					{
						double fai=d_fai*i;
						if(Nt%2==0)
						{
							if(k%2!=0) fai+=0.5*d_fai;//Nt�������Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͊�ԖځB����Ċ�����炷
						}
						else
						{
							if(k%2==0) fai+=0.5*d_fai;//Nt����Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͋����ԖځB����Ċ�����炷
						}
						double xf=r*cos(fai);
						double yf=r*sin(fai);
						double zf=R*cos(THETA);
						
						node++;
						NODE[node].r[A_X]=xf;
						NODE[node].r[A_Y]=yf;
						NODE[node].r[A_Z]=zf+Rz;
						NODE[node].boundary_condition=0;
						//NODE[node].material=FLUID;
						if(Ri<=3) NODE[node].material=FLUID;
						if(Ri>3) NODE[node].material=AIR;
						NODE[node].particleID=0;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;

						if(k!=Nt)
						{
							node++;
							NODE[node].r[A_X]=xf;
							NODE[node].r[A_Y]=yf;
							NODE[node].r[A_Z]=-zf+Rz;
							NODE[node].boundary_condition=0;
							//NODE[node].material=FLUID;
							if(Ri<=3) NODE[node].material=FLUID;
							if(Ri>3) NODE[node].material=AIR;
							NODE[node].particleID=0;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
						}
					}
				}
				if(Nt%2!=0)//Nt����̂Ƃ��́A����ɗ��q���u����Ȃ���΂Ȃ�Ȃ��B���������loop�͂��ꂪ�s�\�B����Ă����Œǉ�
				{
					node++;
					NODE[node].r[A_X]=0;
					NODE[node].r[A_Y]=0;
					NODE[node].r[A_Z]=R+Rz;
					NODE[node].boundary_condition=0;
					//NODE[node].material=FLUID;
					if(Ri<=3) NODE[node].material=FLUID;
					if(Ri>3) NODE[node].material=AIR;
					NODE[node].particleID=0;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;

					
						node++;
						NODE[node].r[A_X]=0;
						NODE[node].r[A_Y]=0;
						NODE[node].r[A_Z]=-R+Rz;
						NODE[node].boundary_condition=0;
						//NODE[node].material=FLUID;
						if(Ri<=3) NODE[node].material=FLUID;
						if(Ri>3) NODE[node].material=AIR;
						NODE[node].particleID=0;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
				
				}
				//////////////////////////////////////
			}

			//*////
			unsigned int timeD=GetTickCount();
			
			/*
			int checkmax=0;
			int checkmin=node;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].BD_node==ON)
				{
					if(i>checkmax) checkmax=i;
					if(i<checkmin) checkmin=i;
				}
			}
			cout<<"���E�ߓ_�ԍ��̍ő�l�A�ŏ��l��"<<checkmax<<","<<checkmin<<endl;

			int checkbig=0;
			int checksmall=0;
			double check=node/2.0;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].BD_node==ON)
				{
					if(i>check) checkbig++;
					if(i<check) checksmall++;
				}
			}
			cout<<"���E�ߓ_�ԍ��F�����̔����Ɣ�r�����召��"<<checkbig<<","<<checksmall<<endl;

			*/
			
			cout<<"node="<<node<<" �f���[�j����"<<endl;

			for(int i=1;i<=nelm;i++)
			{
				ELEM[i].remesh=OFF; //������
				ELEM[i].volume=volume3D(NODE,ELEM[i].node[1],ELEM[i].node[2],ELEM[i].node[3],ELEM[i].node[4]);
				if(ELEM[i].volume<0) cout<<"�̐ϕ��̗v�f�ԍ�="<<i<<" "<<"�̐�="<<ELEM[i].volume<<endl;

				int R1=NODE[ELEM[i].node[1]].remesh;
				int R2=NODE[ELEM[i].node[2]].remesh;
				int R3=NODE[ELEM[i].node[3]].remesh;
				int R4=NODE[ELEM[i].node[4]].remesh;
				
				if(R1==ON && R2==ON && R3==ON && R4==ON) ELEM[i].remesh=ON;//NODE.remesh�́A�ق��̊֐��̓s���ナ���b�V�����E��ON�ɂȂ��Ă���

				sphere3D(NODE,ELEM,ELEM[i].node[1],ELEM[i].node[2],ELEM[i].node[3],ELEM[i].node[4],i);
			}

			for(int i=1;i<KTE;i++) ELEM[i].map=0;//������
			err=1e-14;
			///�����ߓ_�𓱓����Ă���
			int *kv=new int[KTE];//�V�ߓ_���O�ڋ��ɂӂ��ޗv�f�Q
			int *istack=new int[KTE];//�ꎞ�z��

			//////////�ӂ��̏��Ԃɗ��̂�ǉ�
			if(CON->get_defer_f()==OFF)
			{
				for(int i=node0+1;i<=node;i++)
				{   
					//cout<<i<<endl;
					int ip=i;
					double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
					double yp=NODE[ip].r[A_Y];
					double zp=NODE[ip].r[A_Z];
			
					///�V�ߓ_���܂ޗv�f�̒T��
					int loc=locate3D(NODE,ELEM,nelm,xp,yp,zp);
				
					//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
					int iv=0;
					int msk=0;
			
					iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
					kv[iv]=loc;
					ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނ��ǂ����������ς݂Ƃ�������
					msk++;
					istack[msk]=loc;
					if(loc==0)cout<<"loc==0"<<endl;
					
					if(CON->get_CDT_sw()==OFF)//�ʏ�̃f���[�j
					{
						while(msk!=0)
						{   
							int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
								if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
								{
									if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
									{   
					         
										double rad=ELEM[jelm].RR*(1.000000+err);//�O�ڋ����a�̂Q��
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
									
										if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//�O�ڋ����Ɋ܂�
												{
													iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
													kv[iv]=jelm;//���X�g�ɂ����
													ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
													msk++;
													istack[msk]=jelm;
												}
											}
										}
									}
								}
							}
						}
					}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���

					if(CON->get_CDT_sw()==ON)//����t���f���[�j�Bremesh���E���󂳂Ȃ��悤�ɂ�����
					{
						while(msk!=0)
						{   
							int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
								if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
								{
									if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
									{   
										double rad=ELEM[jelm].RR*(1.000000+err);//�O�ڋ����a�̂Q��
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
								
										if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//�O�ڋ����Ɋ܂�
												{
													if(ELEM[jelm].remesh==ON)//isk�Ɛڂ���v�f���\������ߓ_�����ׂ�remesh�ߓ_�A���Ȃ킿remesh�̈���̗v�f�Ȃ�
													{
														iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
														kv[iv]=jelm;//���X�g�ɂ����
														ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
														msk++;
														istack[msk]=jelm;
													}
												}
											}
										}
									}
								}
							}
						}
					}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���
					
					/////////////////
					
					////����ꂽ���ʑ̂��l�ʑ̂ɕ�������
					poly3D2(NODE,ELEM,&iv,kv,ip,&nelm,CON); 
				}
			}

			/////////���ʑ̌`�������܂������Ȃ��������̐ߓ_�̒ǉ�����񂵂ɂ���
			if(CON->get_defer_f()==ON)
			{
				int *flag=new int [node+1];
				for(int i=1;i<=node;i++) flag[i]=ON;
				int *jnb3=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[ ���̐ߓ_���ē������邩�̏���(jnb=0)
				int num_ON=node;//�ǉ��ł��Ă��Ȃ��ߓ_�̐�
				int count=0;
				for(int i=1;i<=node0;i++)//�����b�V���ߓ_�i�N�_�������j�݂̂��V�����ǉ�����ߓ_
				{
					flag[i]=OFF;
					num_ON--;
				}
				
				while(num_ON!=0)
				{
					if(count>0) cout<<"�v�f�`���Ɏ��s�������̐ߓ_�̍Ĕz�u�@num="<<num_ON<<endl;
					//if(count>0) cout<<"���̐ߓ_�̍Ĕz�u�J�n count="<<count+1<<endl;
					count++;
					for(int i=node0+1;i<=node;i++)
					{
						//cout<<i<<endl;
						if(flag[i]==ON)
						{
							int ip=i;
							double xp=NODE[ip].r[A_X];//��������ߓ_�̍��W
							double yp=NODE[ip].r[A_Y];
							double zp=NODE[ip].r[A_Z];
						
							///�V�ߓ_���܂ޗv�f�̒T��
							//unsigned int timeA=GetTickCount();
							int loc=locate3D(NODE,ELEM,nelm,xp,yp,zp);
							//cout<<"locate�����|�|time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
						
							//////////�O�ڋ����ɐV�ߓ_���܂ޗv�f�̒��o
							int iv=0;
							int msk=0;
					
							iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f��
							kv[iv]=loc;
							ELEM[loc].map=1;//map��1�̗v�f�́A�O�ڋ��ɐߓ_i���܂ނ��ǂ����������ς݂Ƃ�������
							msk++;
							istack[msk]=loc;
							if(loc==0)cout<<"loc==0"<<endl;
							
							if(CON->get_CDT_sw()==OFF)//�ʏ�̃f���[�j
							{
								while(msk!=0)
								{   
									int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
										if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
										{
											if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
											{   
							         
												//double rad=ELEM[jelm].RR*(1.000000+err);//�O�ڋ����a�̂Q��
												double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//�O�ڋ����a�̂Q��
												double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
											
												if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
												{
													dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
													if(dst<rad)
													{
														dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
														if(dst<rad)//�O�ڋ����Ɋ܂�
														{
															iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
															kv[iv]=jelm;//���X�g�ɂ����
															ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
															msk++;
															istack[msk]=jelm;
														}
													}
												}
											}
										}
									}
								}
							}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���

							if(CON->get_CDT_sw()==ON)//����t���f���[�j�Bremesh���E���󂳂Ȃ��悤�ɂ�����
							{
								while(msk!=0)
								{   
									int isk=istack[msk];//���ܒ��ڂ��Ă���v�f�̔ԍ�
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//isk�Ɛڂ���v�f
										if(ELEM[jelm].remesh==ON)//isk�Ɛڂ���v�f���\������ߓ_�����ׂ�remesh�ߓ_�A���Ȃ킿remesh�̈���̗v�f�Ȃ�
										{
											if(jelm!=0)//���ꂪ�\�ʂłȂ��Ȃ�
											{
												if(ELEM[jelm].map!=1) //�܂��������ĂȂ��Ȃ�
												{   
													//double rad=ELEM[jelm].RR*(1.000000+err);//�O�ڋ����a�̂Q��
													double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//�O�ڋ����a�̂Q��
													double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///�O�ڋ����S�ƐV�ߓ_�̋���
											
													if(dst<rad)///������dst>rad�Ȃ��΂ɊO�ڋ��Ɋ܂܂Ȃ�
													{
														dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
														if(dst<rad)
														{
															dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
															if(dst<rad)//�O�ڋ����Ɋ܂�
															{
																iv++;//�O�ڋ����ɐV�ߓ_���܂ޗv�f�����{�P
																kv[iv]=jelm;//���X�g�ɂ����
																ELEM[jelm].map=1;//�O�ڋ��ɐV�_�܂ނƂ������邵
																msk++;
																istack[msk]=jelm;
															}
														}
													}
												}
											}
										}
									}
								}
							}//�O�ڋ����ɐV�ߓ_���܂ޗv�f��iv�ƁA���̗v�f�ԍ�kv[iv]�����Ƃ܂���


							
							/////////////////
							
							////����ꂽ���ʑ̂��l�ʑ̂ɕ�������
							//unsigned int timeB=GetTickCount();
							poly3D2(NODE,ELEM,&iv,kv,ip,&nelm,CON);
							//cout<<"poly�����|�|time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

							
						}
					}
				
					set_jnb3D(NODE,ELEM,node,nelm,jnb3);
					for(int i=node0+1;i<=node;i++)
					{
						if(flag[i]==ON) num_ON--;
						flag[i]=OFF;
						
						if(jnb3[i]==0)
						{
							flag[i]=ON;
							num_ON++;
						}
						/*
						if(flag[i]==ON)
						{
							if(jnb3[i]>0)
							{
								num_ON--;
								flag[i]=OFF;
							}
						}
						*/
					}
				}
				delete [] jnb3;
				delete [] flag;
			}
			
			delete [] kv;
			delete [] istack;//*/

			cout<<"�f���[�j���������|�|time="<<(GetTickCount()-timeD)*0.001<<"[sec]"<<endl;

			/////////

			////�ގ�����
			//set_material(CON,NODE,ELEM,node,nelm);
			
			int countmate=0;
			////�ގ�����
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					
					int M1=NODE[ELEM[i].node[1]].material;
					int M2=NODE[ELEM[i].node[2]].material;
					int M3=NODE[ELEM[i].node[3]].material;
					int M4=NODE[ELEM[i].node[4]].material;

					int B1=NODE[ELEM[i].node[1]].BD_node;
					int B2=NODE[ELEM[i].node[2]].BD_node;
					int B3=NODE[ELEM[i].node[3]].BD_node;
					int B4=NODE[ELEM[i].node[4]].BD_node;

					int F1=0;
					int F2=0;
					int F3=0;
					int F4=0;

					if(M1==FLUID || B1==ON) F1=ON;
					if(M2==FLUID || B2==ON) F2=ON;
					if(M3==FLUID || B3==ON) F3=ON;
					if(M4==FLUID || B4==ON) F4=ON;
			
					///4���_���ׂĂ������ގ��Ȃ�v�f������ɂȂ炤�B
					///�ЂƂł��قȂ��Ă������C�ƒ�`
					if(M1==M2 && M2==M3 && M3==M4)
					{
						if(M1==FLUID) ELEM[i].material=FLUID;
					}

					////�����b�V�����E�ߓ_�Ɨ��̐ߓ_�݂̂ō\�������v�f�𗬑̂Ƃ���
					if(CON->get_BD_fluid()==ON)
					{
						if(F1==F2 && F2==F3 && F3==F4)
						{
							if(F1==ON) ELEM[i].material=FLUID;
						}
					}
				}
			}
					//*//////
			for(int i=1;i<=nelm;i++)
			{
				//�ő�Ӓ��������l�𒴂����ꍇ�A���̗v�f�ł͂Ȃ��Ɣ��f����C�v�f�ɕς���
				if( ELEM[i].material==FLUID)
				{
					if(CON->get_material_judge()>0)
					{
						int ia=ELEM[i].node[1];
						int ib=ELEM[i].node[2];
						int ic=ELEM[i].node[3];
						int id=ELEM[i].node[4];
							
						double iaib=0;//�_ia,ib�̋���
						double iaic=0;//�_ia,ic�̋���
						double iaid=0;//�_ia,id�̋���
						double ibic=0;//�_ib,ic�̋���
						double ibid=0;//�_ib,id�̋���
						double icid=0;//�_ic,id�̋���

						for(int D=0;D<3;D++)
						{
							iaib+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);
							iaic+=(NODE[ia].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
							iaid+=(NODE[ia].r[D]-NODE[id].r[D])*(NODE[ia].r[D]-NODE[id].r[D]);
							ibic+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);
							ibid+=(NODE[ib].r[D]-NODE[id].r[D])*(NODE[ib].r[D]-NODE[id].r[D]);
							icid+=(NODE[ic].r[D]-NODE[id].r[D])*(NODE[ic].r[D]-NODE[id].r[D]);
						}
						iaib=sqrt(iaib);
						iaic=sqrt(iaic);
						iaid=sqrt(iaid);
						ibic=sqrt(ibic);
						ibid=sqrt(ibid);
						icid=sqrt(icid);

						double maxL=iaib;//�ő�Ӓ���
						if(iaic>maxL) maxL=iaic;
						if(iaid>maxL) maxL=iaid;
						if(ibic>maxL) maxL=ibic;
						if(ibid>maxL) maxL=ibid;
						if(icid>maxL) maxL=icid;

						if(maxL>CON->get_material_judge()*CON->get_distancebp())
						{
							countmate++;
							ELEM[i].material=AIR;
						}
					}
				}
			}
			if(countmate>0) cout<<"�ő�Ӓ�����"<<CON->get_material_judge()<<"le�𒴂��A��C�ɕϊ����ꂽ���̗v�f�̐�="<<countmate<<endl;


			/////���b�V���̍ו���
			FINE3D(NODE,ELEM,KTJ,KTE,&node,&nelm,CON,1,nelm0);
			//*/
		}
		

		////�v�f�����܂���������Ă��邩�`�F�b�N����
		fill3D(NODE,ELEM,nelm);

		cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;
		/*//���b�V���������m�F
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
		//for(int i=1;i<=dnode;i++) val[i]=dy_NODE[i].material;
		//data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
		//data_avs2(CON,dnode,dnelm,dy_NODE,dy_ELEM,KTJ,val);//�f�ʐ}
		
		if(TIME==0)
		{
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,val,t);//�f�ʐ}
			data_avs3(node,nelm,NODE,ELEM,CON,t);//�ގ�
		}
		delete [] val;
		/////*/
		
	}
	else if(CON->get_mesher()==1)	//TetGen�ɂ�郁�b�V������
	{
		if(CON->get_dimention()==3)
		{
			//delaun3D�֌W�̕ϐ��͑S�ĕs�v

			vector<int> TRANS;		//i�͐ߓ_�ԍ�  �ߓ_i��TRANS[i]�Ԗڂ̗��q�ɑ���
			TRANS.push_back(-1);	//�ŏ��̃f�[�^�͎g��Ȃ��̂ŋl�߂Ă���

			//TetGen�p�z��쐬
			vector<tetgen_node> NODEall;
			vector<tetgen_facet> FACEall;
			vector<tetgen_element> ELEMall;

			///TetGen�ɂ�郁�b�V������////////////////////////////////////////////////////

			tetgen_function TETFUNC;

			if(CON->get_mesh_input()==0)//���ׂĎ����ŗp��
			{
				TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
			}

			//�O������̓Ǎ��͌��ݑΉ����Ă��Ȃ�
			if(CON->get_mesh_input()==1)//magnet�̃��b�V����Ǎ�
			{
				int delaun_flag;			//�f���[�j�������s�����A�s��Ȃ���
				int node0=0;

				if(CON->get_mesh_input()==0)							//MPSTOFEM�ɂ��ߓ_�A�v�f����
				{
					if(CON->get_remesh_sw()==OFF) delaun_flag=FULL;		//remesh�̈��z�肹���A��ɂ��ׂĂ��f���[�j����
					else if(CON->get_remesh_sw()==ON)
					{
						if(t==1) delaun_flag=FULL;	//�S���f�����f���[�j����
						else delaun_flag=REMESH;	//remesh�̈�̂݃f���[�j����
					}
				}
				else if(CON->get_mesh_input()==1)						//Magnet���ǂݍ���
				{
					if(CON->get_remesh_sw()==OFF) delaun_flag=FULL_INPORT;		//���Magnet�̗v�f��ǂݍ��݉�� �Ǘ��җp�H
					else if(CON->get_remesh_sw()==ON)
					{
						if(t==1) delaun_flag=FULL_INPORT;	//�S���f����Magnet�t�@�C�����ǂݍ���
						else delaun_flag=REMESH;	//remesh�̈�̂݃f���[�j����
					}
				}

				if(delaun_flag==FULL) cout<<"FULL �f���[�j�������s"<<endl;
				else if(delaun_flag==REMESH) cout<<"remesh�̈�̂݃f���[�j�������s(tetgen)"<<endl;
				else if(delaun_flag==FULL_INPORT) cout<<"Magnet�����t�@�C�����v�f��񓙓ǂݍ���(tetgen)"<<endl;

				if(delaun_flag==FULL)
				{
					MPS_TO_FEM3Dmain(CON,&node,NODE,PART,  fluid_number,  particle_number);//���q�z�u���ߓ_�z�u�����
					KTJ=node;	
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;
				}
				else if(delaun_flag==REMESH)
				{
					node=(int) static_NODE.size()-1;	//�ÓI�ߓ_��
					nelm=(int) static_ELEM.size()-1;
					KTJ=node+particle_number;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;
				}
				else if(delaun_flag==FULL_INPORT)
				{
					ifstream fin("input_from_MAGNET.dat");
					if(!fin) cout<<"cannot open input_from_MAGNET.dat"<<endl;
					fin.unsetf(ifstream::dec);
					fin.setf(ifstream::skipws);

					fin>>node;			//�ߓ_���ǂݍ���
					fin>>nelm;			//�v�f���ǂݍ���
					KTJ=node+particle_number;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;

					point3D NODE0;
					element3D ELEM0;
					edge3D EDGE0;
					int ID;

					for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);
					for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//�z����m��
					for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);	//�z����m��
					
					//�ߓ_���Ǎ�
					for(int i=1;i<=node;i++)
					{
						fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
						for(int D=0;D<3;D++) fin>>NODE[i].r[D];
						fin>>NODE[i].material;
						if(NODE[i].material==21) NODE[i].material=CRUCIBLE;
						NODE[i].boundary_condition=0;			//�Ƃ肠�����[�����i�[
						NODE[i].particleID=-1;					//�Ή����闱�q�͑��݂��Ȃ�
						NODE[i].remesh=OFF;						//non-remesh�̈�̐ߓ_�ł���Bremesh�̈�Ƃ̋��E�Ɉʒu����ߓ_�Ɋւ��Ă͌�ɏ������{��
						NODE[i].BD_node=OFF;		
					}

					/*//////�ߓ_����static.node�ɏo��
					ofstream fout("static.node");
	
					fout<<"#node"<<endl;
					fout<<node<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=node;i++)
					{
						fout<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].material<<" "<<NODE[i].material<<endl;
					}

					fout.close();
					//////*/

					
					//�v�f-�ߓ_���ǂݍ���
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
						for(int j=1;j<=4;j++) fin>>ELEM[i].node[j];
					}
					//�v�f-�v�f���ǂݍ���
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
						for(int j=1;j<=4;j++) fin>>ELEM[i].elm[j];
					}
					//�v�f�ގ����ǂݍ���
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//�ǂ��ł��������Ȃ̂ŉ��ɂ��g�p���Ȃ��B�������t�@�C���̍��[�ɂ͂��̂悤�ɐߓ_�ԍ������Ă����悤�ɁB���̂ق����`�F�b�N���y
						fin>>ELEM[i].material;
						//////
						if(ELEM[i].material==21) ELEM[i].material=CRUCIBLE;
						/////
						ELEM[i].map=0;			//������
						for(int D=0;D<3;D++)ELEM[i].r[D]=0;
						ELEM[i].RR=0;
						ELEM[i].volume=0;		
					}
					/*//////�v�f����static.ele�ɏo��
					ofstream fout2("static.ele");
	
					fout2<<nelm<<" "<<"4"<<" "<<"1"<<endl;
					for(int i=1;i<=nelm;i++)
					{
						fout2<<i<<" "<<ELEM[i].node[1]<<" "<<ELEM[i].node[2]<<" "<<ELEM[i].node[3]<<" "<<ELEM[i].node[4]<<" "<<ELEM[i].material<<endl;
					}

					fout2.close();

					fin.close();
					/////////////*/
					/*///////////���̐ߓ_��fluid.node�ɏo��
					ofstream fp("rr.dat");
					ofstream fl("fluid.node");
					/////
					fl<<"#node"<<endl;
					fl<<fluid_number<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=fluid_number;i++)
					{
						fl<<i<<" "<<PART[i-1].r[A_X]<<" "<<PART[i-1].r[A_Y]<<" "<<PART[i-1].r[A_Z]<<" "<<WATER<<" "<<WATER<<endl;
					}
					

					fp.close();
					fl.close();
					///////*/

					//////static.node,ele���烁�b�V���č\��
					//TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
					////////

					//���E�����ݒ�
					int *BD_flag=new int[node+1];			//BD_flag=ON�Ȃ狫�E�ߓ_(��͋��E��������Ȃ����Aremesh���E��������Ȃ�)
					for(int i=0;i<=node;i++) BD_flag[i]=OFF;//������
					for(int i=1;i<=nelm;i++)
					{
						for(int j=1;j<=4;j++)
						{
							if(ELEM[i].elm[j]==0)
							{
								int ia=ELEM[i].node[j%4+1];
								int ib=ELEM[i].node[4-(j-1)/2*2];
								int ic=ELEM[i].node[3-(j/2%2)*2];
								BD_flag[ia]=ON;						//���E�ߓ_�Ƃ�����
								BD_flag[ib]=ON;
								BD_flag[ic]=ON;
							}
						}
					}
					vector <int> BD_NODE_ID;						//���E�ߓ_�ԍ��i�[
					for(int i=1;i<=node;i++)
					{
						if(BD_flag[i]==ON) BD_NODE_ID.push_back(i);
					}
					int BD_num=(int) BD_NODE_ID.size();
					ofstream fs("remesh.dat");
					if(CON->get_region_shape()==0)				//��͗̈悪�����̂Ȃ�
					{
						double Xmax=0;
						double Ymax=0;
						double Zmax=0;

				

						for(int i=0;i<BD_num;i++)
						{
							int n=BD_NODE_ID[i];	//���E�ߓ_�ԍ�
							if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z��
							else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z��
							else if(NODE[n].r[A_X]<CON->get_XL()+err) NODE[n].boundary_condition=1;		//-X��
							else if(NODE[n].r[A_X]>CON->get_XR()-err) NODE[n].boundary_condition=1;		//+X��
							else if(NODE[n].r[A_Y]<CON->get_YD()+err) NODE[n].boundary_condition=1;		//-Y��
							else if(NODE[n].r[A_Y]>CON->get_YU()-err) NODE[n].boundary_condition=1;		//+Y��
							else
							{
								NODE[n].boundary_condition=0;					//��͋��E�ߓ_�ł͂Ȃ��Aremesh�̈�Ƃ̋��E�ߓ_�Ȃ̂ŁA���E�����̓[��
								NODE[n].remesh=ON;
								NODE[n].BD_node=ON;							//remesh�̈�̋��E���ł���Ƃ������邵
								fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
							}
							if(NODE[n].r[A_Z]>Zmax) Zmax=NODE[n].r[A_Z];
							if(NODE[n].r[A_Y]>Ymax) Ymax=NODE[n].r[A_Y];
							if(NODE[n].r[A_X]>Xmax) Xmax=NODE[n].r[A_X];
						}
						cout<<Xmax<<" "<<Ymax<<" "<<Zmax<<endl;
					}
					else if(CON->get_region_shape()==1)				//��͗̈悪�~���Ȃ�
					{
						for(int i=0;i<BD_num;i++)
						{
							int n=BD_NODE_ID[i];	//���E�ߓ_�ԍ�
							double R=sqrt(NODE[n].r[A_X]*NODE[n].r[A_X]+NODE[n].r[A_Y]*NODE[n].r[A_Y]);
							if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z��
							else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z��
							else if(R>CON->get_RU()*0.99) NODE[n].boundary_condition=1;	
							else
							{
								NODE[n].boundary_condition=0;											//��͋��E�ߓ_�ł͂Ȃ��Aremesh�̈�Ƃ̋��E�ߓ_�Ȃ̂ŁA���E�����̓[��
								NODE[n].remesh=ON;
								NODE[n].BD_node=ON;							//remesh�̈�̋��E���ł���Ƃ������邵
								fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
							}
						}
					}
					delete [] BD_flag;
					fs.close();

					//NODE,ELEM�̂����A�����Ȃ��v�f�A�ߓ_������static_NODE,staticELEM�Ɋi�[����
					static_ELEM.clear();
					static_NODE.clear();
					memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);
			
					if(CON->get_remesh_sw()==ON)
					{
						node=(int) static_NODE.size()-1;
						nelm=(int) static_ELEM.size()-1;
						delaun_flag=REMESH;				//delaun_flag��REMESH�ɂ��邱�ƂŁA����if���ɓ����ē��I�v�f�𐶐�����
						NODE.clear();
						ELEM.clear();
					}
					else cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;

				}
		
				if(delaun_flag==REMESH)									//remesh�̈���f���[�j����
				{
					point3D NODE0;
					element3D ELEM0;
					for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);			//�ߓ_�z��m�� +8�̓X�[�p�[�{�b�N�X
					for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);				//�z����m��

					for(int i=1;i<=node;i++)										//���̎��_��node�ɂ͐ÓI�ߓ_�����i�[����Ă���
					{
						for(int D=0;D<3;D++) NODE[i].r[D]=static_NODE[i].r[D];
						NODE[i].boundary_condition=static_NODE[i].boundary_condition;
						NODE[i].material=static_NODE[i].material;
						NODE[i].particleID=static_NODE[i].particleID;
						NODE[i].remesh=static_NODE[i].remesh;
						NODE[i].BD_node=static_NODE[i].BD_node;
					}

					//static_ELEM�������copy
					for(int i=1;i<=nelm;i++)										//���̎��_��nelm�ɂ͐ÓI�v�f�����i�[����Ă���
					{
						for(int D=0;D<3;D++) ELEM[i].r[D]=static_ELEM[i].r[D];
						for(int j=1;j<=4;j++)
						{
							ELEM[i].node[j]=static_ELEM[i].node[j];
							ELEM[i].elm[j]=static_ELEM[i].elm[j];
							//if(ELEM[i].elm[j]==0) cout<<i<<endl;
						}
						ELEM[i].map=static_ELEM[i].map;
						ELEM[i].material=static_ELEM[i].material;
						ELEM[i].RR=static_ELEM[i].RR;
						ELEM[i].volume=static_ELEM[i].volume;
						//�ӂ́H�H�H
					}

					//�ÓI�v�f��remesh�̈�ɐڂ���Ƃ��낪ELEM[i].elm=0�ƂȂ��Ă���.����ȗv�f��T��
					vector<int> BD_static_ELEM;		//���I�v�f�ɐڂ���ÓI�v�f�ԍ��i�[
					vector<int> BD_static_ELEM_elm;	//�ÓI�v�f�����I�v�f�ɐڂ���ʔԍ��i�[
					for(int i=1;i<=nelm;i++)
					{
						int flag=OFF;
						int J;
						for(int j=1;j<=4;j++) if(NODE[ELEM[i].node[j]].remesh==ON) flag=ON; 
						if(flag==ON)
						{
					
							flag=OFF;
							int num=0;
							for(int j=1;j<=4;j++)
							{
								if(ELEM[i].elm[j]==0)
								{
									flag=ON;
									J=j;
									num++;
								}
							}
							//if(num>1) cout<<i<<" �ÓI�v�f�������ʂœ��I�v�f�Ɛڂ��Ă��܂� mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].boundary_condition<<NODE[ELEM[i].node[2]].boundary_condition<<NODE[ELEM[i].node[3]].boundary_condition<<NODE[ELEM[i].node[4]].boundary_condition<<endl;
							if(num>1) cout<<i<<" �ÓI�v�f�������ʂœ��I�v�f�Ɛڂ��Ă��܂� mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].r[A_Z]<<endl;
						}
						if(flag==ON)
						{
							BD_static_ELEM.push_back(i);
							BD_static_ELEM_elm.push_back(J);//�v�fi�͑�J�ʂœ��I�v�f�ɐڂ��Ă���
						}
					}
					cout<<"KK="<<BD_static_ELEM.size()<<endl;

					double Pn[3];				//�N�_�̍��W

					for(int k=0;k<1;k++)
					{
						int i=BD_static_ELEM[k];
						int j=BD_static_ELEM_elm[k];
						int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
						int ib=ELEM[i].node[4-(j-1)/2*2];
						int ic=ELEM[i].node[3-(j/2%2)*2];

						double iaic[3];//ia��ic���޸�ِ����i�[
						double iaib[3];//ia��ib���޸�ِ����i�[
						double iaicL=0;	//ia��ic�̒���
						for(int D=0;D<3;D++)
						{
							iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
							iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
							iaicL+=iaic[D]*iaic[D];
						}
						iaicL=sqrt(iaicL);
						///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
						double S[3];//��L���޸�ِ����i�[
						S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
						S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
						S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
						double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
						////�ʐ�S�����Ƃ܂���

						double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�@���޸��

						double Gp[3];								//�\�ʎO�p�`�̏d�S
						for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

						for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//�N�_�̍��W

						//////////////////////
						///�N�_�������I�ɍ��W�w��
				
						Pn[A_X]=0;
						Pn[A_Y]=0;
						Pn[A_Z]=0.225;//0.10325;		//��ڒꕔ�F0.10125

						//���̂̍ő卂�������߁A����̏�����ɐݒu

						/*//���̐ߓ_�̍ő卂���A�ŏ�����
						double Zmax=0;
						double Zmin=10;
						for(int i=1;i<=node;i++)
						{
							if(NODE[i].material==FLUID)
							{
								if(NODE[i].r[A_Z]>Zmax) Zmax=NODE[i].r[A_Z];
								if(NODE[i].r[A_Z]<Zmin) Zmin=NODE[i].r[A_Z];
							}
						}

						Pn[A_X]=0;
						Pn[A_Y]=0;
						Pn[A_Z]=Zmax-CON->get_distancebp();
						*///

						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=AIR;
						NODE[node].particleID=-1;
						NODE[node].remesh=ON;

					}
					cout<<"�N�_����"<<endl;
					int countBD=0;
					for(int i=1;i<=node;i++) if(NODE[i].BD_node==ON) countBD++;
					cout<<"���E�̐ߓ_��="<<countBD<<endl;

					int *imen[4];
					for(int D=0;D<4;D++) imen[D]=new int [BD_static_ELEM.size()+1];
					int *jmen =new int [BD_static_ELEM.size()+1];
					int *kmen =new int [BD_static_ELEM.size()+1];
					double *vol =new double [BD_static_ELEM.size()+1];
					int ip=node;

					/*
					if(BD_static_ELEM.size()>100000) cout<<"imen�Ȃǂ̃������������Ă�������"<<endl;
					int ip=node;			
					int imen[100000][3+1];	//���ʑ̕\�ʎO�p�`�̐ߓ_�ԍ��i�[
					int jmen[100000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̔ԍ��i�[
					int kmen[100000];		//���ʑ̕\�ʎO�p�`�ɗאڂ���l�ʑ̗̂אږʔԍ� (����͑扽�ʂŎ����Ɛڂ��Ă��邩)
					double vol[100000];		//���ʑ̂̑̐ς̂U�{
					*/

					for(int k=1;k<=BD_static_ELEM.size();k++)
					{
						int kelm=BD_static_ELEM[k-1];//remesh�̈�ɐڂ���ÓI�v�f
						int j=BD_static_ELEM_elm[k-1];
						int ia=ELEM[kelm].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
						int ib=ELEM[kelm].node[4-(j-1)/2*2];
						int ic=ELEM[kelm].node[3-(j/2%2)*2];

						/*
						imen[k][1]=ic;
						imen[k][2]=ib;
						imen[k][3]=ia;
						*/
						imen[1][k]=ic;
						imen[2][k]=ib;
						imen[3][k]=ia;

						jmen[k]=kelm;
						kmen[k]=j;//jelm��ielm�ɐڂ���ʔԍ�
						//vol[k]=volume3D(NODE,ia,ib,ic,ip);
						vol[k]=volume3D(NODE,ic,ib,ia,ip);	//�ߓ_�̏��Ԃɒ���
						//if(vol[k]<0) cout<<"�̐ϕ��̗v�f���� "<<vol[k]<<endl;
					}

					int ibound=(int) BD_static_ELEM.size();//�\�ʂ̐���\��
					int nelm0=nelm;//�ύX�O�̗v�f�����L��
			
					for(int i=1;i<=ibound;i++)//�v�f��񐶐�
					{   
						nelm++;
						int ielm=nelm0+i;
						int ia=imen[1][i];
						int ib=imen[2][i];
						int ic=imen[3][i];
						ELEM[ielm].node[1]=ia;
						ELEM[ielm].node[2]=ib;
						ELEM[ielm].node[3]=ic;
						ELEM[ielm].node[4]=ip;//�V�_�͂S�Ԗڂƒ�`	
						ELEM[ielm].elm[4]=jmen[i];
						if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
						ELEM[ielm].volume=vol[i];
				
						sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
			
						ELEM[ielm].material=AIR;//�ʏ��poly�֐��Ƃ̑���_�B�����ōގ�����C�ƌ��肷��
					}
					///////////////////

					//�v�f-�v�f�֌W�C��/////////��̏����ő�4�ʂŐڂ���v�f�ԍ��͂킩���Ă���̂ŁA�c������߂�
					//						�����ŁA1�`3�ʂ͑��ʑ̂��\������v�f�Ƃ̋��E�ʂł��邱�Ƃɒ���
					int ix=0;
			
					for(int i=1;i<=ibound;i++)
					{
						int ielm=nelm0+i;
						for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]�͂��łɂ��Ƃ܂�������A����ȊO�����Ƃ߂�
						{
							///ELEM[ielm].node[4]=ip�ł���
							int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3�̂Ƃ��A2,3,1�̏�
							int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3�̂Ƃ��A3,1,2�̏�
							int flag=0;
							for(int k=1;k<=ix;k++)
							{
								if(flag==0)
								{
									int ja=imen[1][k];
									int jb=imen[2][k];
									if(ia==ja && ib==jb)//�ߓ_����v������
									{
										ELEM[ielm].elm[j]=jmen[k];//���炩����ؽĂ��Ă����������i�[
										ELEM[jmen[k]].elm[kmen[k]]=ielm;
										imen[1][k]=imen[1][ix];		//k�Ԗڂ̏��͂����s�v�B�Ȃ̂Ŕz��̈�ԍŌ�̏���k�Ԗڂɂ����Ă��āA����܂ł̏��͔j������
										imen[2][k]=imen[2][ix];
										jmen[k]=jmen[ix];
										kmen[k]=kmen[ix];
										ix--;						//�҂��Ӑ�����
										flag=1;						//ELEM[ielm].elm[j]�͂��Ƃ܂����̂ŁA���̃l�X�g�ɓ���K�v�͂Ȃ��̂�flag=1
									}
								}
							}
							if(flag==0)
							{
								ix++;			//�����ł�ix�́A[�אڊ֌W�𖞂����v�f]���܂��Ă���[��]�̐���\���B
								imen[1][ix]=ib;	//�����̐ߓ_�̕��т��L�������A�ʂ̗v�f�����̕��т𖞂����̂�҂Bib��ia�̕��т��t�ɂ��Ă��邱�Ƃɒ���
								imen[2][ix]=ia;
								jmen[ix]=ielm;
								kmen[ix]=j;
							}
						}	
					}///�v�f-�v�f�֌W�C������

					cout<<"�v�f�쐬����"<<endl;

					fill3D(NODE,ELEM,nelm);

					node0=node;
					cout<<"node0="<<node0<<endl;

					///////�ߓ_����static.node�ɏo��
					ofstream fout("static.node");
	
					fout<<"#node"<<endl;
					fout<<node<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=node;i++)
					{
						fout<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].material<<" "<<NODE[i].material<<endl;
					}

					fout.close();
					//////*/

					
					///////�v�f����static.ele�ɏo��
					ofstream fout2("static.ele");
	
					fout2<<nelm<<" "<<"4"<<" "<<"1"<<endl;
					for(int i=1;i<=nelm;i++)
					{
						fout2<<i<<" "<<ELEM[i].node[1]<<" "<<ELEM[i].node[2]<<" "<<ELEM[i].node[3]<<" "<<ELEM[i].node[4]<<" "<<ELEM[i].material<<endl;
					}

					fout2.close();

					/////////////*/

					////////////���̐ߓ_��fluid.node�ɏo��
					ofstream fl("static-a.node");
					/////
					fl<<"#node"<<endl;
					fl<<fluid_number<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=fluid_number;i++)
					{
						fl<<i<<" "<<PART[i-1].r[A_X]<<" "<<PART[i-1].r[A_Y]<<" "<<PART[i-1].r[A_Z]<<" "<<AIR<<" "<<AIR<<endl;
					}
					
					fl.close();
					////////

					//////static.node,ele���烁�b�V���č\��
					TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
					////////

					/////////
					for(int D=0;D<4;D++) delete [] imen[D];
					delete [] jmen;
					delete [] kmen;
					delete [] vol;


					/*/���E�ʂ���w���`��
					for(int m=1;m<=1;m++)
					{
						for(int k=1;k<BD_static_ELEM.size();k++)
						{
							int i=BD_static_ELEM[k];
							int j=BD_static_ELEM_elm[k];
							int ia=ELEM[i].node[j%4+1];//���_���̌����Ɉʒu����O�p�`�̒��_ ia��ib��ic�͂�����݂Ĕ����v�܂��
							int ib=ELEM[i].node[4-(j-1)/2*2];
							int ic=ELEM[i].node[3-(j/2%2)*2];

							double iaic[3];//ia��ic���޸�ِ����i�[
							double iaib[3];//ia��ib���޸�ِ����i�[
							double iaicL=0;	//ia��ic�̒���
							for(int D=0;D<3;D++)
							{
								iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
								iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
								iaicL+=iaic[D]*iaic[D];
							}
							iaicL=sqrt(iaicL);
							///0.5(iaic�~iaib)�͎O�p�`ia,ib,ic�̖ʐς̑傫���������A�����͊O�������޸�قƂȂ�
							double S[3];//��L���޸�ِ����i�[
							S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
							S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
							S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
							double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
							////�ʐ�S�����Ƃ܂���

							double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//�@���޸��

							double Gp[3];								//�\�ʎO�p�`�̏d�S
							for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

							//for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//�N�_�̍��W �I���W�i��
							for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*0.0001*m;		//�N�_�̍��W

							node++;
							for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
							NODE[node].boundary_condition=0;
							NODE[node].material=AIR;
							NODE[node].particleID=-1;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
						}
					}
					*/


					double u0=PI*4E-7;			//��C�̓�����
					double skin_depth=sqrt(1.0/(PI*CON->get_Hz()*CON->get_ele_conduc()*u0));//�\��[��
					cout<<"���̕\��[��="<<skin_depth<<endl;
			
					/////���̗��q�𓮓I�ߓ_�Ƃ��Ċi�[
					ofstream fp("rr.dat");
					if(CON->get_thinout_fluid()==0)
					{
						/////
						for(int i=0;i<fluid_number;i++)
						{
							//if(PART[i].surface==ON  || i%4==0)
							node++;
							for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
							NODE[node].boundary_condition=0;
							NODE[node].material=FLUID;
							NODE[node].particleID=i;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
							fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
						}
						////*/
					}
					if(CON->get_thinout_fluid()>0)
					{
						double Xf=0.0;
						double Yf=0.0;
						double Zf=0.0;
						for(int i=0;i<fluid_number;i++)
						{
							Xf+=PART[i].r[A_X];
							Yf+=PART[i].r[A_Y];
							Zf+=PART[i].r[A_Z];
						}
						Xf/=fluid_number;
						Yf/=fluid_number;
						Zf/=fluid_number;
						cout<<"���̏d�S=("<<Xf<<","<<Yf<<","<<Zf<<")"<<endl;
				
						for(int i=0;i<fluid_number;i++)
						{
							double rx=PART[i].r[A_X]-Xf;
							double ry=PART[i].r[A_Y]-Yf;
							double rz=PART[i].r[A_Z]-Zf;
							if(PART[i].surface==ON || sqrt(rx*rx+ry*ry+rz*rz)>CON->get_fluidwidth()*0.001-1.1*skin_depth)
							{
								node++;
								for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
								NODE[node].boundary_condition=0;
								NODE[node].material=FLUID;
								NODE[node].particleID=i;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
								fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
							}
							else if(i%CON->get_thinout_fluid()==0)
							{
								node++;
								for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
								NODE[node].boundary_condition=0;
								NODE[node].material=FLUID;
								NODE[node].particleID=i;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
								fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
							}
						}
					}
					fp.close();

					/*///////
					//�\�ʂ̐��x����̂��߁A���̐ߓ_��ǉ�
					double Rz=0.13125-0.003;
					//�܂��͔��������B���̂��߂ɂ͔����\�ʂ��쐬����K�v������B
					for(int Ri=1;Ri<=5;Ri++)
					{
						double R=0.025-5*0.0001+Ri*0.0002;
						double le=CON->get_distancebp();

						double A=sqrt(3.0)/2;				//�悭�g���W��
						double B=sqrt(2.0/3);						////�悭�g���W��
						int half_WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
						int half_WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
						int half_WZ=(int)(R/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
						double R2=R-0.5*le;				//���������߂̔��a��ݒ�

						///////////�����\��
						int Nt;						//���\�ʂ́A�ƕ����̕�����
						double Lt;					//���\�ʂ́A�ƕ����̕�������
						//calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//���~�̕������͋����E��ǂ���ł��悢
						double temp_N=PI/2*R/(le*A);			//���̕������Ble�Ŋ���؂ꂽ���Ԃ������ǁA�����������Ȃ��Ƃ�������
						int Ns=(int) temp_N;				//�^�̕�����
						double difference=temp_N-Ns;		//���Ɛ^�̍�
						if(difference>0.5) Ns++;
						Lt=PI/2*R/Ns;			//���q�̋���
						Nt=Ns;

						double d_theta=Lt/R;		//�ʂ̒�����Lt�ɂȂ�p�x

						for(int k=0;k<=Nt;k++)//loop��k<Nt�ŏI��点��BNt�ɊY������Ƃ���͂��łɐݒu�ς�
						{
							double THETA=k*d_theta;	//��
							double r=R*sin(THETA);	//���̍����ɂ�����~�̔��a
							double round=2*PI*r;//���̍����ɂ�����~��

							//int Nf=calc_division_N_circle(round,le);//���\�ʂ́A�ƕ����̕�����
							//�Ώ̐����l�������A�~�����L�q���闱�q�͋����łȂ���΂Ȃ�Ȃ��B�����瑼�̕ӕ������Ƃ͈�������������
							//dis:�������鋗��(�~��)
							double temp_num=round/le;		//�~�O���ɐݒu����w���́x���q���B�������O�������܂�le�Ŋ���؂��Ƃ͌���Ȃ�

							int N1=(int)(temp_num/2);
							N1*=2;							
							int N2=N1+2;					//temp_num��N1��N2�̊Ԃɂ���B������N1,N2�͋���

							double dif1=temp_num-N1;		//�eN�Ƃ̍�
							double dif2=N2-temp_num;
							int N=N1;						//������������
							if(dif2<dif1) N=N2;				//���̏���������N�Ƃ��č̗p����B

							int Nf=N;
					
					
							double Lf=round/Nf;						//���\�ʂ́A�ƕ����̕�������
							double d_fai=Lf/r;						//�ʂ̒�����Lf�ɂȂ�p�x
					
							for(int i=0;i<Nf;i++)
							{
								double fai=d_fai*i;
								if(Nt%2==0)
								{
									if(k%2!=0) fai+=0.5*d_fai;//Nt�������Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͊�ԖځB����Ċ�����炷
								}
								else
								{
									if(k%2==0) fai+=0.5*d_fai;//Nt����Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͋����ԖځB����Ċ�����炷
								}
								double xf=r*cos(fai);
								double yf=r*sin(fai);
								double zf=R*cos(THETA);
						
								node++;
								NODE[node].r[A_X]=xf;
								NODE[node].r[A_Y]=yf;
								NODE[node].r[A_Z]=zf+Rz;
								NODE[node].boundary_condition=0;
								//NODE[node].material=FLUID;
								if(Ri<=3) NODE[node].material=FLUID;
								if(Ri>3) NODE[node].material=AIR;
								NODE[node].particleID=0;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;

								if(k!=Nt)
								{
									node++;
									NODE[node].r[A_X]=xf;
									NODE[node].r[A_Y]=yf;
									NODE[node].r[A_Z]=-zf+Rz;
									NODE[node].boundary_condition=0;
									//NODE[node].material=FLUID;
									if(Ri<=3) NODE[node].material=FLUID;
									if(Ri>3) NODE[node].material=AIR;
									NODE[node].particleID=0;
									NODE[node].remesh=ON;
									NODE[node].BD_node=OFF;
								}
							}
						}
						if(Nt%2!=0)//Nt����̂Ƃ��́A����ɗ��q���u����Ȃ���΂Ȃ�Ȃ��B���������loop�͂��ꂪ�s�\�B����Ă����Œǉ�
						{
							node++;
							NODE[node].r[A_X]=0;
							NODE[node].r[A_Y]=0;
							NODE[node].r[A_Z]=R+Rz;
							NODE[node].boundary_condition=0;
							//NODE[node].material=FLUID;
							if(Ri<=3) NODE[node].material=FLUID;
							if(Ri>3) NODE[node].material=AIR;
							NODE[node].particleID=0;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;

					
								node++;
								NODE[node].r[A_X]=0;
								NODE[node].r[A_Y]=0;
								NODE[node].r[A_Z]=-R+Rz;
								NODE[node].boundary_condition=0;
								//NODE[node].material=FLUID;
								if(Ri<=3) NODE[node].material=FLUID;
								if(Ri>3) NODE[node].material=AIR;
								NODE[node].particleID=0;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
				
						}
						//////////////////////////////////////
					}

					//*////
					unsigned int timeD=GetTickCount();
			
					/*
					int checkmax=0;
					int checkmin=node;
					for(int i=1;i<=node;i++)
					{
						if(NODE[i].BD_node==ON)
						{
							if(i>checkmax) checkmax=i;
							if(i<checkmin) checkmin=i;
						}
					}
					cout<<"���E�ߓ_�ԍ��̍ő�l�A�ŏ��l��"<<checkmax<<","<<checkmin<<endl;

					int checkbig=0;
					int checksmall=0;
					double check=node/2.0;
					for(int i=1;i<=node;i++)
					{
						if(NODE[i].BD_node==ON)
						{
							if(i>check) checkbig++;
							if(i<check) checksmall++;
						}
					}
					cout<<"���E�ߓ_�ԍ��F�����̔����Ɣ�r�����召��"<<checkbig<<","<<checksmall<<endl;

					*/
			
					cout<<"node="<<node<<" �f���[�j����"<<endl;
				}
			}
			
			///////////////////////////////////////////////////////////////////////////////

			int N=(int)TRANS.size()-1;			//FEM�ߓ_�Ɋ܂܂�闱�q��(0�Ԗڂ�����TRANS[]�̒���)
			node=(int)NODEall.size();		//�ߓ_��
			nelm=(int)ELEMall.size();		//�v�f��
			KTJ=node;						//�ő�ߓ_��
			KTE=nelm;						//�ő�v�f��
			//int *depth=new int [KTE+1];			//�e�v�f�̐[���i�[

			cout<<"tetgen���node="<<node<<endl;
			cout<<"tetgen���nelm="<<nelm<<endl;

			///�z��m��
			//point3D *NODE2=new point3D [NODEall.size()+1];
			//element3D *ELEM2=new element3D [ELEMall.size()+1];
			
			//NODE.clear();
			//ELEM.clear();
			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			//int ID;
			
			for(int i=0;i<=node+1;i++) NODE.push_back(NODE0);
			for(int i=0;i<=nelm+1;i++) ELEM.push_back(ELEM0);	//�z����m��
			for(int i=0;i<2*KTE;i++) EDGE.push_back(EDGE0);	//�z����m��
			//NODE.resize(NODEall.size()+1);
			//ELEM.resize(ELEMall.size()+1);
			
			for(int i=0;i<node;i++) if(NODEall[i].attribute==WATER) NODEall[i].attribute=FLUID;//tetgen�ł͗��̂�WATER�Ƃ��Ă���̂ŁAFLUID�ɕύX
			for(int i=0;i<nelm;i++) if(ELEMall[i].attribute==WATER) ELEMall[i].attribute=FLUID;//tetgen�ł͗��̂�WATER�Ƃ��Ă���̂ŁAFLUID�ɕύX

			for(int i=0;i<node;i++) if(NODEall[i].attribute==CONDUCT) NODEall[i].attribute=CRUCIBLE;
			for(int i=0;i<nelm;i++) if(ELEMall[i].attribute==CONDUCT) ELEMall[i].attribute=CRUCIBLE;
			
			//TetGen�̐ߓ_�E�v�f�f�[�^���擾	�ߓ_�ԍ���1���炷
			//�ߓ_�f�[�^	
			for(int i=0;i<node;i++)
			{
				NODE[i+1].r[A_X]=NODEall[i].r[A_X];
				NODE[i+1].r[A_Y]=NODEall[i].r[A_Y];
				NODE[i+1].r[A_Z]=NODEall[i].r[A_Z];
				NODE[i+1].material=(int)NODEall[i].attribute;
			}
			//�v�f�f�[�^
			for(int i=0;i<nelm;i++)
			{
				ELEM[i+1].material=ELEMall[i].attribute;								//�ގ�
				for(int n=0;n<4;n++)	ELEM[i+1].node[n+1]=ELEMall[i].node[n]+1;		//�\���ߓ_
				for(int n=0;n<4;n++)													//�v�f-�v�f�֌W 
				{
					if(ELEMall[i].nei_elem[n]==-1)	ELEM[i+1].elm[n+1]=0;	//�ߗחv�f�Ȃ��̏ꍇTetGen��-1��Ԃ��Ă��邽��0�ɏC������
					else							ELEM[i+1].elm[n+1]=ELEMall[i].nei_elem[n]+1;
				}

			}

			//�̐ς���ъO�ڋ��p�����[�^�v�Z
			for(int i=1;i<=nelm;i++)
			{
				//4�̐ߓ_
				int ia=ELEM[i].node[1];
				int ib=ELEM[i].node[2];
				int ic=ELEM[i].node[3];
				int ip=ELEM[i].node[4];

				//�v�f�̐όv�Z
				ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);
				
				///�̐σ`�F�b�N
				if(ELEM[i].volume<0)			{cout<<"i="<<i<<" �̐ς����ł�"<<endl;				cout<<"volume="<<ELEM[i].volume<<endl;}
				else if(ELEM[i].volume==0)		{cout<<"i="<<i<<" �̐ς�0�ł�"<<endl;				cout<<"volume="<<ELEM[i].volume<<endl;}
				//else if(ELEM[i].volume<1e-20)	{cout<<"i="<<i<<" �̐ς����������܂� 1e-20"<<endl;	cout<<"volume="<<ELEM[i].volume<<endl;}
				//*/
				//if(ELEM[i].volume==0)		{cout<<"i="<<i<<" �̐ς�0�ł�"<<endl;	ELEM[i].volume=1e-50;}
				//if(ELEM[i].volume<min_volume && ELEM[i].volume!=0)	min_volume=ELEM[i].volume;

				//�O�ڋ����S���W����є��a�v�Z
				sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
			}

			////���f�����̐ݒ�(���E�����Ȃ�)
			//�Ód����
			if(CON->get_model_number()==14)
			{
				for(int i=0;i<node;i++)
				{
					//i��1����Ă��邱�Ƃɒ���
					if(NODEall[i].boundary==ELECTRODE1)			NODE[i+1].boundary_condition=1;	//�~���d�ɂ���ѓy��
					else if(NODEall[i].boundary==ELECTRODE2)	NODE[i+1].boundary_condition=2;	//���d��
					else										NODE[i+1].boundary_condition=0;	//���̑��̐ߓ_�͖��m��
				}
			}
			if(CON->get_model_number()==25)
			{
				for(int i=1;i<=node;i++)
				{
							
					if(NODE[i].r[A_Z]<CON->get_ZD()+err) NODE[i].boundary_condition=1;		//-Z��
					else if(NODE[i].r[A_Z]>CON->get_ZU()-err) NODE[i].boundary_condition=1;		//+Z��
					else if(NODE[i].r[A_X]<CON->get_XL()+err) NODE[i].boundary_condition=1;		//-X��
					else if(NODE[i].r[A_X]>CON->get_XR()-err) NODE[i].boundary_condition=1;		//+X��
					else if(NODE[i].r[A_Y]<CON->get_YD()+err) NODE[i].boundary_condition=1;		//-Y��
					else if(NODE[i].r[A_Y]>CON->get_YU()-err) NODE[i].boundary_condition=1;		//+Y��
					else
					{
						NODE[i].boundary_condition=0;					//��͋��E�ߓ_�ł͂Ȃ�
					}
				}
			}
	//cout<<"�ŏ��̐�="<<min_volume<<endl;
		}
	}
	/////
	int countf=0;
	int countco=0;
	int countcru=0;
	int counta=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].material==FLUID) countf++;
		if(NODE[i].material==COIL) countco++;
		if(NODE[i].material==CRUCIBLE) countcru++;
		if(NODE[i].material==AIR) counta++;
	}

	if(countf>=0) cout<<"���̐ߓ_="<<countf<<endl;
	if(countco>=0) cout<<"�R�C���ߓ_="<<countco<<endl;
	if(countcru>=0) cout<<"��ڐߓ_="<<countcru<<endl;
	if(counta>=0) cout<<"��C�ߓ_="<<counta<<endl;
	///*/

	/*///��ڍő卂��
	double hmax=0;
	double hmin=10000;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].material==CRUCIBLE)
		{
			if(NODE[i].r[A_Z]>hmax) hmax=NODE[i].r[A_Z];
			if(NODE[i].r[A_Z]<hmin) hmin=NODE[i].r[A_Z];
		}
	}
	cout<<"���Z���W �ő�="<<hmax<<" �ŏ�="<<hmin<<endl;
	///*/

	
	
	
	
	/////�ߓ_-�v�f�֌W
    int *jnb=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//���b�V���X���[�W���O
	//if(CON->get_mesh_smth()>0) for(int i=0;i<CON->get_mesh_smth();i++) laplacian_smoothing2(NODE,ELEM,&node,&nelm,CON,jnb,nei,node0);
	
	/////���b�V���������m�F
	//if(TIME==0)
	//if(t==1 || (t-1)%(CON->get_EM_interval()*CON->get_mesh_output_interval())==0) 
	if(t==1)
	{
		data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//�f�ʐ}�B����t�@�C����������
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
		data_avs2node(CON,node,nelm,NODE,ELEM,val,t);
		data_avs3(node,nelm,NODE,ELEM,CON,t);//�ގ�
		delete [] val;
	}
	//*/

	countf=0;
	countco=0;
	countcru=0;
	counta=0;
	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<endl;
			NODE[i].boundary_condition=1;//���E�������ިظڌ^�ɂ��邱�ƂŁAICCG�ɎQ�������Ȃ�
			if(NODE[i].material==FLUID)
			{
				NODE[i].particleID=-1;//���̐ߓ_����������ꍇ�A�K�v�Ȓl���v�Z����Ȃ����߁A���q�Ƀt�B�[�h�o�b�N�ł��Ȃ�(���Ă͂����Ȃ�)
				countf++;
			}
			if(NODE[i].material==COIL) countco++;
			if(NODE[i].material==CRUCIBLE) countcru++;
			if(NODE[i].material==AIR) counta++;
		}
	}

	if(countf>0) cout<<"�v�f������Ȃ����̐ߓ_���� ����"<<countf<<endl;
	if(countco>0) cout<<"�v�f������Ȃ��R�C���ߓ_���� ����"<<countco<<endl;
	if(countcru>0) cout<<"�v�f������Ȃ���ڐߓ_���� ����"<<countcru<<endl;
	if(counta>0) cout<<"�v�f������Ȃ���C�ߓ_���� ����"<<counta<<endl;

	int node_sta=(int) static_NODE.size()-1;
	cout<<"static_NODE.size()="<<node_sta<<endl;
	//���̂̑̐ϒ���
	double volume_f=0;//���̗v�f�̑̐ς̍��v
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID)
		{
			int ia=ELEM[i].node[1];
			int ib=ELEM[i].node[2];
			int ic=ELEM[i].node[3];
			int ip=ELEM[i].node[4];
		
			double volme=volume3D(NODE,ia,ib,ic,ip);
			//if(volme<=0) cout<<"�̐ς��[������"<<endl;
			volume_f+=volme;
			
			//ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//�̐ς�6�{�ł��邱�Ƃɒ���
			//sphere3D(NODE,ELEM,ia,ib,ic,ip,i);//�O�ڋ��̒��S�i�{���m�C�_)�Ɣ��a�̓����v�Z
		}
	}
	volume_f/=6;
	cout<<"���̑̐�="<<volume_f<<endl;
	
	if(t==1)
	{
		ofstream fv("volume.dat");
		fv.close();
	}
	ofstream vo("volume.dat",ios :: app);
	vo<<t<<" "<<volume_f<<endl;
	vo.close();

	//////�ߓ_�ԍ��̕��ёւ�
	if(CON->get_node_sort()==1) node_sorting(CON,node,nelm,NODE,ELEM,jnb,nei);
	else if(CON->get_node_sort()==2) node_sorting2(CON,node,nelm,NODE,ELEM,jnb,nei);
	//////
	

	/////�ߓ_-�v�f�֌W �ߓ_�ԍ����ύX���ꂽ�̂ł�����x���߂�
    int *jnb3=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb3);
    int **nei3=new int* [node+1];
    for(int i=1;i<=node;i++) nei3[i]=new int [jnb[i]+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    set_nei3D(NODE,ELEM,node,nelm,jnb3,nei3);


	//�ȉ~������
	if(CON->get_EM_calc_type()==1 || CON->get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb3, TIME,PART, fluid_number,nei3);
	if(CON->get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,nedge,NODE,ELEM,EDGE,jnb3, dt, TIME,PART, fluid_number,F,t,nei3,particle_node,NODE_jw,ELEM_jw,node_sta,static_EDGE);
		
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
	delete [] nei;

	delete [] jnb3;
    for(int i=1;i<=node;i++) delete [] nei3[i];
    delete [] nei3;
	return 0;
}

//�Q�d�����l���ɂ��ꂽ�ߓ_�v�f���޸�����ݼ�ٌv�Z�֐�ver.2 ���ݹް�ނ�K�p���č��ӂ����׼�݂Ƃ���
void Avector3D_node_eddy2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t)
{   
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	double j0x,j0y,j0z;					//�d�����x
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///���m��
    int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)//�s���A1x,A1y,A1z,A1��,A2x,A2y�E�E�E�̏��Ɋi�[
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;

			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-3;
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-3;
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///�s��̍ő啝�v�Z  �����Ƃ����̑傫���̂̓ӂ̗v�f���낤�Ɖ��肵�Ă���B�m���ł͂Ȃ����Ƃɒ���
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=4*mat_w;//X,Y,Z,�Ӑ����Ƃ���̂Ł~4 //���ӂ����׼�݂ɂ��Ă��A�ӂɊւ��Ă͂S�����̂܂܂ɒ���

	////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    //*////

	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=4*mat_w;//
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////�z��m��
	double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

	int width_max=0;
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
		width_max+=width_mat[i];
    }

	cout<<"���v��="<<width_max<<endl;

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;

	delete [] width_mat;	
	delete [] width_node;	

	//*/
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
   
    
    /////////�S�̍s����쐬����
	cout<<"�S�̍s��쐬�J�n"<<endl;
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=delta/4*j0x;//x,y,z�͂P�������
				B[I]+=delta/4*j0y;//
				B[I+1]+=delta/4*j0z;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flag=1;
							}
						}
						if(flag==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=J;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y����
				I=npp[N[n]]+1+1;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z����
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						J3=J2+1;		//Aiz�̍�
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J3;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	//////�Q�d�����v�Z
	int J4,flag4;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sigma;//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sigma;//   ��/(36V)
			//double co3=sigma*dt*delta6;

			double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*dt*delta6;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aix�̍�
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag=1;
								}
							
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flag4=1;
								}///
							}
							if(flag==0)//��������i�����݂��Ȃ���������(�����ł͂���Ȃ͂��Ȃ�����)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==J) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=J4;
							}

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J) Sx+=2*Ax;
							else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����

					I=npp[N[n]]+1+1;/////Y����
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J4=J+3;			//�ӂ̍�
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J2) G[I][h]+=co*2;//Aiy�̍�
									else G[I][h]+=co;
									flag2=1;
								}

								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(I==J2) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J2;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=J4;
							}//*/

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J2) Sy+=2*Ay;
							else Sy+=Ay;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sy;//�x�z��������X�������瓾����B�̒l

					/////Z����
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==J3) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    if(I==J3) G[I][H]+=co*2;
								else G[I][H]+=co;
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=J4;
							}//*/
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							if(I==J3) Sz+=2*Az;
							else Sz+=Az;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sz;//�x�z��������X�������瓾����B�̒l

					/////��
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aix�̍�
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flag=1;
								}
								//Aiy�̍�
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flag2=1;
								}
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=J4;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

	//�Ίp�����̒l�`�F�b�N
	cout<<"�Ίp�����̒l�`�F�b�N"<<endl;
	int count_plus=0;
	int count_minus=0;
	int count_zero=0;
	for(int i=1;i<=pn;i++)
	{
		int flag=OFF;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				if(G[i][j]==0)
				{
					count_zero++;
					if(ppn[i-1]==-6) cout<<"�Ίp��������(��i�j i="<<i<<endl;
					else if(ppn[i-1]==-7) cout<<"�Ίp��������(��r�j i="<<i<<endl;
					else cout<<"�Ίp��������(A�j i="<<i<<endl;
				}
				if(G[i][j]>0) count_plus++;
				if(G[i][j]<0) count_minus++;
				flag=ON;

			}
		}
		if(flag==OFF) count_zero++;
	}
	cout<<"���̑Ίp���̐�="<<count_plus<<endl;
	cout<<"���̑Ίp���̐�="<<count_minus<<endl;
	cout<<"��̑Ίp���̐�="<<count_zero<<endl;
	///

	//�Ίp�����̕��z
	diagonal_distribution(CON,node,nelm,NODE,ELEM,pn,NUM,ROW,G);

	 double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/

	cout<<"��됔="<<number<<endl;


	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//*///

	//�s��̎��o��
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	

	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) A[A_X][i]=XX[n];
		else if(i==-1) A[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) A[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) V[ppn[n-3]]=XX[n];//�d�ʃ�
	}	
	
	delete [] XX;

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///�Q�d���Ƌ����d���𑫂����l���o�͂�����
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_A�̍X�V 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_A��̧�ُo��
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

////////////
////////////
////////////
//�Q�d�����l���ɂ��ꂽ�ߓ_�v�f���޸�����ݼ�ٌv�Z�֐� j�֖@ �ߓ_���ƂɎ����A���������ꂼ��i�[
void Avector3D_node_eddy2_jw(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i���֖@�j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///���m�� ���f���Ȃ̂Ŋe�����Ɏ����Ƌ���������
    int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)//�s���Re(A1x),Im(A1x),Re(A1y),Im(A1y),Re(A1z),Re(��1),Im(��1),Re(A2x),�E�E�E�̏��Ɋi�[
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			ppn[num+3]=-3;
			ppn[num+4]=-4;
			ppn[num+5]=-5;
			npp[i]=num;
			num+=6;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;
						ppn[num+1]=-7;
						num+=2;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;
						ppn[num+1]=-7;
						num+=2;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///�s��̍ő啝�v�Z  �����Ƃ����̑傫���̂̓ӂ̗v�f���낤�Ɖ��肵�Ă���B�m���ł͂Ȃ����Ƃɒ���
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8
    *//////

	/*///�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //*/
    

	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=8;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=8;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=6;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=6;
			}
		}
	}	
    //////

	////�z��m��
	double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;

	delete [] width_mat;	
	delete [] width_node;	

    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X�����A����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=delta/4*j0x.real();//Rex
				B[I]+=delta/4*j0x.imag();//Imx
				B[I+1]+=delta/4*j0y.real();//Rey
				B[I+2]+=delta/4*j0y.imag();//Imy
				B[I+3]+=delta/4*j0z.real();//Rez
				B[I+4]+=delta/4*j0z.imag();//Imz
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxr=1;
							}
						}
						if(flagxr==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X�����A����

				I=npp[N[n]]+2;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxi=1;
							}
						}
						if(flagxi==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y����,����
				I=npp[N[n]]+3;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+3;			//Re(Aiy)�̍�
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y����,����
				I=npp[N[n]]+4;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+4;			//Im(Aiy)�̍�
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z����,����
				I=npp[N[n]]+5;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+5;	//Aix�̍�
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z����,����
				I=npp[N[n]]+6;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+6;	//Aix�̍�
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	//////�Q�d�����v�Z
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����,����

					int I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxi)
									{
										//cout<<"�Ίp���ɑ�����Ă���H"<<endl;
										G[I][h]-=co*2;//�N����Ȃ��͂�
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////X����,����

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							//Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����,����

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////�x����,����

					I=npp[N[n]]+4;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							//Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�

					/////Z����,����

					I=npp[N[n]]+5;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////�y����,����

					I=npp[N[n]]+6;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							//Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������(���łɂ���͂��j
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�

					/////��,���� //P77��(4.27)��dt�������� co3�������鍀�̋����ɒ���
					
					I=npp[N[n]]+7;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+8;	//Im�ӂ̍�

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///B�̌v�Z
							//���̍��Ɠ������Aold_A���K�v�Ȃ��̂�B�͍X�V����Ȃ�
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////��,����
					
					I=npp[N[n]]+8;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							//Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							//Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiy�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aiz�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Re�ӂ̍� G=��(��i-j��r)
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	cout<<"G,ROW�̏��ԕ��ёւ��J�n"<<endl;
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N //��Ώ̍s��Ȃ̂Ń`�F�b�N���Ȃ�
	//check_matrix_symmetry(pn,NUM,ROW,G);

	 double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//*///

	//�s��̎��o��
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	

	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AR[A_X][i]=XX[n];
		else if(i==-1) AI[A_X][ppn[n-1]]=XX[n];
		else if(i==-2) AR[A_Y][ppn[n-2]]=XX[n];
		else if(i==-3) AI[A_Y][ppn[n-3]]=XX[n];
		else if(i==-4) AR[A_Z][ppn[n-4]]=XX[n];
		else if(i==-5) AI[A_Z][ppn[n-5]]=XX[n];
		else if(i==-6) VR[ppn[n-6]]=XX[n];//�d�ʃ�
		else if(i==-7) VI[ppn[n-6]]=XX[n];//�d�ʃ�
	}	

	delete [] XX;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	double Am=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*t+phi);
		}
	}

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///�Q�d���Ƌ����d���𑫂����l���o�͂�����
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_A�̍X�V 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_A��̧�ُo��
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

//////////
//j�֖@,�����������̏��Ɋi�[
void Avector3D_node_eddy2_jw2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i���֖@-2�j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///���m�� ���f���Ȃ̂Ŋe�����Ɏ����Ƌ���������
    int pnI=3*(node-NN)+conducter_num;//���m���̂����������̐�
	int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)//�s���Im(A1x),Im(A1y),Im(A1z),Im(��1),Im(A2x),�E�E�E,Im(��n),Re(A1x),Re(A1y)�̏��Ɋi�[
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;//Imx
			ppn[num+1]=-1;//Imy
			ppn[num+2]=-2;//Imz
			ppn[num+pnI]=-3;//Rex
			ppn[num+1+pnI]=-4;//Rey
			ppn[num+2+pnI]=-5;//Rez
			npp[i]=num;
			num+=3;//4�߈ȍ~�͎����Ƃ��ĕ�������Ă���̂Œ���
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;//Im��
						ppn[num+pnI]=-7;
						num+=1;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;//Im��
						ppn[num+pnI]=-7;
						num+=1;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///�s��̍ő啝�v�Z  �����Ƃ����̑傫���̂̓ӂ̗v�f���낤�Ɖ��肵�Ă���B�m���ł͂Ȃ����Ƃɒ���
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8
    //////

	////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //*/
    

	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=4;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=3;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////�z��m��
	double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	//delete [] width_mat;	
	//delete [] width_node;	

    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //*/
    
	cout<<"�S�̍s��쐬�J�n pn="<<pn<<" pnI="<<pnI<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//if(je%10==0) cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X�����A����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=delta/4*j0x.imag();//Imx
				B[I]+=delta/4*j0y.imag();//Imy
				B[I+1]+=delta/4*j0z.imag();//Imz
				B[I-1+pnI]+=delta/4*j0x.real();//Rex
				B[I+pnI]+=delta/4*j0y.real();//Rey
				B[I+1+pnI]+=delta/4*j0z.real();//Rez
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxi=1;
							}
						}
						if(flagxi==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y����,����
				I=npp[N[n]]+2;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+2;			//Im(Aiy)�̍�
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z����,����
				I=npp[N[n]]+3;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+3;	//Aix�̍�
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X�����A����

				I=npp[N[n]]+1+pnI;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxr=1;
							}
						}
						if(flagxr==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y����,����
				I=npp[N[n]]+2+pnI;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+2+pnI;			//Re(Aiy)�̍�
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				
				/////*/

				/////Z����,����
				I=npp[N[n]]+3+pnI;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+3+pnI;	//Aix�̍�
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	cout<<"�x�N�g���|�e���V����������"<<endl;

	//////�Q�d�����v�Z
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����,����

					int I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							//Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////�x����,����

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							//Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;

					/////�y����,����

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							//Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�


					/////��,����//co3�������鍀�̋����ɒ���
					I=npp[N[n]]+4;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiy�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aiz�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Re�ӂ̍� G=��(��i-j��r)
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////X����,����

					 I=npp[N[n]]+1+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxi)
									{
										//cout<<"�Ίp���ɑ�����Ă���H"<<endl;
										G[I][h]-=co*2;//�N����Ȃ��͂�
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					

					/////Y����,����

					I=npp[N[n]]+2+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�

					/////Z����,����

					I=npp[N[n]]+3+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					
					/////��,���� //co3�������鍀�̋����ɒ���
					
					I=npp[N[n]]+4+pnI;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///B�̌v�Z
							//���̍��Ɠ������Aold_A���K�v�Ȃ��̂�B�͍X�V����Ȃ�
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
					
				}
			}
		}
	}//////*/
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N //��Ώ̍s��Ȃ̂Ń`�F�b�N���Ȃ�
	//check_matrix_symmetry(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	/*////
	int countImIm=0;//Im�sIm��
	int countImRe=0;
	int countReRe=0;
	int countReIm=0;
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						if(i<=pnI && ROW[i][j]<=pnI) countImIm++;
						if(i>pnI && ROW[i][j]>pnI) countReRe++;

						if(i<=pnI && ROW[i][j]>pnI) countImRe++;
						if(i>pnI && ROW[i][j]<=pnI) countReIm++;

						//cout<<"�Ώ̐��װ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						//cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						if(ppn[i-1]>0 && ppn[J-1]==-3) cout<<"�Ώ̐��װ (AxI,AxR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-1 && ppn[J-1]==-4) cout<<"�Ώ̐��װ (AyI,AyR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-2 && ppn[J-1]==-5) cout<<"�Ώ̐��װ (AzI,AzR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-6 && ppn[J-1]==-7) cout<<"�Ώ̐��װ(��I,��R) ("<<i<<","<<J<<")="<<G[i][j]<<" ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}
	if(countImIm>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countImIm<<endl;
	if(countImRe>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countImRe<<endl;
	if(countReIm>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countReIm<<endl;
	if(countReRe>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countReRe<<endl;
	////*/

	
	/*////���ی��݂̂̑Ώ̐��𒲂ׂ�
	cout<<"���ی��̑Ώ̐��`�F�b�N"<<endl;
	double **G2=new double *[pnI+1];///�S�̍s��
    for(int i=1;i<=pnI;i++) G2[i]=new double [width_mat[i]+1];
	int **ROW2=new int *[pnI+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pnI;i++) ROW2[i]=new int [width_mat[i]+1];
	int *NUM2=new int [pnI+1]; ///�e�s�́A��[���v�f��

	 for(int i=1;i<=pnI;i++)//������
    {
        NUM2[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G2[i][j]=0;
			ROW2[i][j]=0;
		}
    }

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			if(ROW[i][j]>pnI)
			{
				NUM2[i]=NUM2[i]+1;
				ROW2[i][NUM2[i]]=ROW[i][j]-pnI;
				G2[i][NUM2[i]]=G[i][j];
			}
		}
	}

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM2[i];j++)
		{
			int J=ROW2[i][j];
			int flag=0;
			for(int k=1;k<=NUM2[J];k++)
			{
				if(ROW2[J][k]==i)
				{
					flag=1;
					if(G2[i][j]!=G2[J][k])
					{
						cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G2[i][j]<<"  ("<<J<<","<<ROW2[J][k]<<")="<<G2[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}

	for(int i=0;i<=pnI;i++) delete [] G2[i];
    delete [] G2;
	for(int i=0;i<=pnI;i++) delete [] ROW2[i];
    delete [] ROW2;
	delete [] NUM2;
	////////*/

	double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//*///

	//�s��̎��o��
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
  

	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AI[A_X][i]=XX[n];
		else if(i==-1) AI[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) AI[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) AR[A_X][ppn[n-pnI]]=XX[n];
		else if(i==-4) AR[A_Y][ppn[n-pnI-1]]=XX[n];
		else if(i==-5) AR[A_Z][ppn[n-pnI-2]]=XX[n];
		else if(i==-6) VI[ppn[n-3]]=XX[n];//�d�ʃ�
		else if(i==-7) VR[ppn[n-3-pnI]]=XX[n];//�d�ʃ�
	}	

	delete [] XX;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*t+phi);
		}
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);

	}

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///�Q�d���Ƌ����d���𑫂����l���o�͂�����
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_A�̍X�V 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_A��̧�ُo��
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] width_mat;	
	delete [] width_node;	
}

//////////
//////////
//j�֖@3,���f���s����쐬 //FEM�����邽�сA�ÓI�v�f���܂߂đS�̈�Ńx�N�g���|�e���V�������v�Z����
void Avector3D_node_eddy2_jw3(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw)
{ 
	/////
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i���֖@-3�j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//�f�B���N���^���l

    ///�f�B���N���^���E��������

	ifstream old("old_A.dat");
	if(!old) cout<<"cannot open old_A.dat"<<endl;
	old.unsetf(ifstream::dec);
	old.setf(ifstream::skipws);

	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];

	old.close();	
			
    for(int i=1;i<=node;i++)
    {
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		if(CON->get_static_dirichlet()==ON)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
    }//////////////
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���

	
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///���m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����ƕς��Ȃ�
    //int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    //int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
	vector <int> ppn;
	vector <int> npp;
	npp.reserve(node+1);
    int num=0; 

    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			//ppn[num]=i;
			//ppn[num+1]=-1;
			//ppn[num+2]=-2;
			ppn.push_back(i);
			ppn.push_back(-1);
			ppn.push_back(-2);
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }

	cout<<"num="<<num<<"pn="<<pn<<endl;
	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=4*mat_w;//
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////�z��m��
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aix�̍�
							    flag=1;
							}
						}
						if(flag==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aix�̍�
						    ROW[I][H]=J;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}

				/////Y����
				I=npp[N[n]]+1+1;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiy�̍�
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}
				/////*/

				/////Z����
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						J3=J2+1;		//Aiz�̍�
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiz�̍�
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aiz�̍�
						    ROW[I][H]=J3;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}
			}
		}
	}

	//////�Q�d�����v�Z
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aix�̍�
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//��������i�����݂��Ȃ���������(�����ł͂���Ȃ͂��Ȃ�����)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����

					I=npp[N[n]]+1+1;/////Y����
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J4=J+3;			//�ӂ̍�
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//�x�z��������X�������瓾����B�̒l

					/////Z����
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sz;//�x�z��������X�������瓾����B�̒l

					/////��
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aix�̍�
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiy�̍�
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"��됔="<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	
	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	ofstream ar("old_AR.dat");
	ofstream ai("old_AI.dat");
	ofstream vr("old_VR.dat");
	ofstream vi("old_VI.dat");

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D];//�g�����o��
				//if(CON->get_jw_Faverage()==ON) A[D][i]=Am[D]*sqrt(2.0);
				//else A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
		ar<<AR[A_X][i]<<" "<<AR[A_Y][i]<<" "<<AR[A_Z][i]<<endl;
		ai<<AI[A_X][i]<<" "<<AI[A_Y][i]<<" "<<AI[A_Z][i]<<endl;
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}
	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}

	/*/////old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Je�ɂ͔g�̍���������

	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//�{���̑̐�

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//

	/*//old_A��̧�ُo��
	if(t==1)
	{
		cout<<"�����X�e�b�v�̃x�N�g���|�e���V�����L��"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	///*/
	
	
	//���̃X�e�b�v�ō쐬����NODE,ELEM��NODE_jw,ELEM_jw�ɋL��������B�����ŋ��߂��g���ƒx����ȍ~�̃X�e�b�v�ŗ��p����B�����߂邽��
	NODE_jw.clear();
	ELEM_jw.clear();
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    //delete [] npp;
    //delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

//////////
//j�֖@4,�s�����Ɨ�����̊i�[���Ԃ�ς���(�v�A���܂������Ȃ��j
void Avector3D_node_eddy2_jw4(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i���֖@-4�j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<3*NN<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///���m�� ���f���Ȃ̂Ŋe�����Ɏ����Ƌ���������
    int pnI=3*(node-NN)+conducter_num;//���m���̂����������̐�
	int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)//�s���Im(A1x),Im(A1y),Im(A1z),Im(��1),Im(A2x),�E�E�E,Im(��n),Re(A1x),Re(A1y)�̏��Ɋi�[
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;//Imx
			ppn[num+1]=-1;//Imy
			ppn[num+2]=-2;//Imz
			ppn[num+pnI]=-3;//Rex
			ppn[num+1+pnI]=-4;//Rey
			ppn[num+2+pnI]=-5;//Rez
			npp[i]=num;
			num+=3;//4�߈ȍ~�͎����Ƃ��ĕ�������Ă���̂Œ���
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;//Im��
						ppn[num+pnI]=-7;
						num+=1;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;//Im��
						ppn[num+pnI]=-7;
						num+=1;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///�s��̍ő啝�v�Z  �����Ƃ����̑傫���̂̓ӂ̗v�f���낤�Ɖ��肵�Ă���B�m���ł͂Ȃ����Ƃɒ���
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8
    //////

	////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //*/
    

	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=8*mat_w;//X,Y,Z,�Ӑ����ɂ��ꂼ��Re,Im������̂Ł~8

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=4;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=3;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////�z��m��
	double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	//delete [] width_mat;	
	//delete [] width_node;	

    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    //*/
    
	cout<<"�S�̍s��쐬�J�n pn="<<pn<<" pnI="<<pnI<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//if(je%10==0) cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X�����A����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=delta/4*j0x.real();//Rex
				B[I]+=delta/4*j0y.real();//Rey
				B[I+1]+=delta/4*j0z.real();//Rez
				B[I-1+pnI]+=delta/4*j0x.imag();//Imx
				B[I+pnI]+=delta/4*j0y.imag();//Imy
				B[I+1+pnI]+=delta/4*j0z.imag();//Imz
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxr=1;
							}
						}
						if(flagxr==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y����,����
				I=npp[N[n]]+2;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+2+pnI;			//Re(Aiy)�̍�
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				
				/////*/

				/////Z����,����
				I=npp[N[n]]+3;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+3+pnI;	//Aix�̍�
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X�����A����

				 I=npp[N[n]]+1+pnI;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
							    flagxi=1;
							}
						}
						if(flagxi==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aix�̍�
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y����,����
				I=npp[N[n]]+2+pnI;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+2;			//Im(Aiy)�̍�
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiy�̍�
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z����,����
				I=npp[N[n]]+3+pnI;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+3;	//Aix�̍�
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiz�̍�
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	cout<<"�x�N�g���|�e���V����������"<<endl;

	//////�Q�d�����v�Z
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����,����

					int I=npp[N[n]]+1+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							//Jxi=npp[N[m]]+2;	//Im(Aix)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////�x����,����

					I=npp[N[n]]+2+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							//Jyi=npp[N[m]]+4;	//Im(Aiy)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;

					/////�y����,����

					I=npp[N[n]]+3+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							//Jzi=npp[N[m]]+6;	//Im(Aiz)�̍�
							//Jvr=npp[N[m]]+7;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzr) G[I][h]+=co*2;//�N����Ȃ��͂�
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//�N����Ȃ��͂�
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�


					/////��,����//co3�������鍀�̋����ɒ���
					I=npp[N[n]]+4+pnI;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�
							

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiy�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aiz�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Re�ӂ̍� G=��(��i-j��r)
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////X����,����

					 I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)�̍�
								if(Jxi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jxi)
									{
										//cout<<"�Ίp���ɑ�����Ă���H"<<endl;
										G[I][h]-=co*2;//�N����Ȃ��͂�
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					

					/////Y����,����

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)�̍�
								if(Jyi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jyi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					//�x�z��������X�������瓾����B�̒l//���֖@����B�������Ȃ�

					/////Z����,����

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							//Jvi=npp[N[m]]+8;	//Im�ӂ̍�
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)�̍�
								if(Jzi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									if(I==Jzi) G[I][h]-=co*2;//�N����Ȃ��͂�
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Re�ӂ̍�
								if(Jvr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//�N����Ȃ��͂�
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_A�͂��֖@�̏ꍇ�K�v�Ȃ��̂ŁA������B�͕ύX����Ȃ�

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					
					/////��,���� //co3�������鍀�̋����ɒ���
					
					I=npp[N[n]]+4;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)�̍�
							Jxi=npp[N[m]]+1;	//Im(Aix)�̍�
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)�̍�
							Jyi=npp[N[m]]+2;	//Im(Aiy)�̍�
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)�̍�
							Jzi=npp[N[m]]+3;	//Im(Aiz)�̍�
							Jvr=npp[N[m]]+4+pnI;	//Re�ӂ̍�
							Jvi=npp[N[m]]+4;	//Im�ӂ̍�

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)�̍�
								if(Jxr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)�̍�
								if(Jyr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)�̍�
								if(Jzr==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Im�ӂ̍�
								if(Jvi==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///B�̌v�Z
							//���̍��Ɠ������Aold_A���K�v�Ȃ��̂�B�͍X�V����Ȃ�
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
					
				}
			}
		}
	}//////*/
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�Ώ̐��`�F�b�N"<<endl;
	check_matrix_symmetry(pn,NUM,ROW,G);

	//�Ίp�����̒l�`�F�b�N
	cout<<"�Ίp�����̒l�`�F�b�N"<<endl;
	int count_plus=0;
	int count_minus=0;
	int count_zero=0;
	for(int i=1;i<=pn;i++)
	{
		int flag=OFF;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				if(G[i][j]==0)
				{
					count_zero++;
					if(ppn[i-1]==-6) cout<<"�Ίp��������(��i�j i="<<i<<endl;
					else if(ppn[i-1]==-7) cout<<"�Ίp��������(��r�j i="<<i<<endl;
					else cout<<"�Ίp��������(A�j i="<<i<<endl;
				}
				if(G[i][j]>0) count_plus++;
				if(G[i][j]<0) count_minus++;
				flag=ON;

			}
		}
		if(flag==OFF) count_zero++;
	}
	cout<<"���̑Ίp���̐�="<<count_plus<<endl;
	cout<<"���̑Ίp���̐�="<<count_minus<<endl;
	cout<<"��̑Ίp���̐�="<<count_zero<<endl;
	///

	//�Ώ̐��`�F�b�N
	/*////
	int countImIm=0;//Im�sIm��
	int countImRe=0;
	int countReRe=0;
	int countReIm=0;
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						if(i<=pnI && ROW[i][j]<=pnI) countImIm++;
						if(i>pnI && ROW[i][j]>pnI) countReRe++;

						if(i<=pnI && ROW[i][j]>pnI) countImRe++;
						if(i>pnI && ROW[i][j]<=pnI) countReIm++;

						//cout<<"�Ώ̐��װ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						//cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						if(ppn[i-1]>0 && ppn[J-1]==-3) cout<<"�Ώ̐��װ (AxI,AxR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-1 && ppn[J-1]==-4) cout<<"�Ώ̐��װ (AyI,AyR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-2 && ppn[J-1]==-5) cout<<"�Ώ̐��װ (AzI,AzR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-6 && ppn[J-1]==-7) cout<<"�Ώ̐��װ(��I,��R) ("<<i<<","<<J<<")="<<G[i][j]<<" ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}
	if(countImIm>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countImIm<<endl;
	if(countImRe>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countImRe<<endl;
	if(countReIm>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countReIm<<endl;
	if(countReRe>0) cout<<"�����s,������̍��Ŕ�Ώ� num="<<countReRe<<endl;
	////*/

	
	/*////���ی��݂̂̑Ώ̐��𒲂ׂ�
	cout<<"���ی��̑Ώ̐��`�F�b�N"<<endl;
	double **G2=new double *[pnI+1];///�S�̍s��
    for(int i=1;i<=pnI;i++) G2[i]=new double [width_mat[i]+1];
	int **ROW2=new int *[pnI+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pnI;i++) ROW2[i]=new int [width_mat[i]+1];
	int *NUM2=new int [pnI+1]; ///�e�s�́A��[���v�f��

	 for(int i=1;i<=pnI;i++)//������
    {
        NUM2[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G2[i][j]=0;
			ROW2[i][j]=0;
		}
    }

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			if(ROW[i][j]>pnI)
			{
				NUM2[i]=NUM2[i]+1;
				ROW2[i][NUM2[i]]=ROW[i][j]-pnI;
				G2[i][NUM2[i]]=G[i][j];
			}
		}
	}

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM2[i];j++)
		{
			int J=ROW2[i][j];
			int flag=0;
			for(int k=1;k<=NUM2[J];k++)
			{
				if(ROW2[J][k]==i)
				{
					flag=1;
					if(G2[i][j]!=G2[J][k])
					{
						cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G2[i][j]<<"  ("<<J<<","<<ROW2[J][k]<<")="<<G2[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}

	for(int i=1;i<=pnI;i++) delete [] G2[i];
    delete [] G2;
	for(int i=1;i<=pnI;i++) delete [] ROW2[i];
    delete [] ROW2;
	delete [] NUM2;
	////////*/

	double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//*///

	//�s��̎��o��
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
  

	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AI[A_X][i]=XX[n];
		else if(i==-1) AI[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) AI[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) AR[A_X][ppn[n-pnI]]=XX[n];
		else if(i==-4) AR[A_Y][ppn[n-pnI-1]]=XX[n];
		else if(i==-5) AR[A_Z][ppn[n-pnI-2]]=XX[n];
		else if(i==-6) VI[ppn[n-3]]=XX[n];//�d�ʃ�
		else if(i==-7) VR[ppn[n-3-pnI]]=XX[n];//�d�ʃ�
	}	

	delete [] XX;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*TIME+phi);
		}
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*TIME+phi);

	}

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///�Q�d���Ƌ����d���𑫂����l���o�͂�����
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_A�̍X�V 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_A��̧�ُo��
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] width_mat;	
	delete [] width_node;	
}

//j�֖@5,���f���s����쐬 //�Q�X�e�b�v�ڈȍ~�A�ÓI�v�f�̃x�N�g���|�e���V�����l�͂P�X�e�b�v�ڂ̒l���f�B���N���l�Ƃ��ė��p����
void Avector3D_node_eddy2_jw5(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw)
{ 
	/////
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i���֖@-5 �j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
	int NN_V=0;//�f�B���N���^���E�ߓ_���̂����A�Q�d�����l����ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	
	complex<double> *PHAT_A[3];
	for(int D=0;D<3;D++) PHAT_A[D]=new complex<double> [node+1];//A�̃f�B���N���^���l
	complex<double> *PHAT_V = new complex<double> [node+1];//V�̃f�B���N�����E�l
	

    ///�f�B���N���^���E��������
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//������
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	/////////1�X�e�b�v�ڂ̉�ǂݍ���
	cout<<"�ÓI�v�f�̃f�B���N���l�ǂݍ���"<<endl;
	//AR
	ifstream ar("old_AR.dat");
	if(!ar) cout<<"cannot open old_AR.dat"<<endl;
	ar.unsetf(ifstream::dec);
	ar.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ar>>AR[D][i];
	ar.close();

	//AI
	ifstream ai("old_AI.dat");
	if(!ai) cout<<"cannot open old_AI.dat"<<endl;
	ai.unsetf(ifstream::dec);
	ai.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ai>>AI[D][i];
	ai.close();

	//VR
	ifstream vr("old_VR.dat");
	if(!vr) cout<<"cannot open old_VR.dat"<<endl;
	vr.unsetf(ifstream::dec);
	vr.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vr>>VR[i];
	vr.close();

	//VI
	ifstream vi("old_VI.dat");
	if(!vi) cout<<"cannot open old_VI.dat"<<endl;
	vi.unsetf(ifstream::dec);
	vi.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vi>>VI[i];
	vi.close();
	/////
	cout<<"�f�B���N���l�ǂݍ��݊���"<<endl;		
	int count_r=0;
    for(int i=1;i<=node;i++)
    {
		/*///���̊֐��ɂ�static_dirichlet=OFF�̂Ƃ����ŗ��Ȃ�
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		////*/
		if(CON->get_static_dirichlet()==ON)
		{
			/*///
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			///*/
			////else
			
			{
				if(NODE[i].remesh==OFF)//�����X�e�b�v�ȊO�̓����b�V�����Ȃ��ߓ_�ɂ��čŏ��ɋ��߂��x�N�g���|�e���V�������f�B���N���l�Ƃ��ė^����
				{
					count_r++;
					NODE[i].boundary_condition=3;//�ÓI�v�f�̋��E���������Ƃ�0����3�ɏ���������B�Œ苫�E(1,2)���ς��邪���Ȃ��͂�
					dn[i]=NN;
					PHAT_A[A_X][NN]=complex<double>(AR[A_X][i],AI[A_X][i]);
					PHAT_A[A_Y][NN]=complex<double>(AR[A_Y][i],AI[A_Y][i]);
					PHAT_A[A_Z][NN]=complex<double>(AR[A_Z][i],AI[A_Z][i]);
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;

					
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) NN_V++;
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				
				}
				else
				{
					dn[i]=node+1;
				}
			}
		}
    }//////////////
    cout<<"�ިظڐߓ_����"<<NN<<"�ިظڐ���"<<3*NN+NN_V<<" NN_V="<<NN_V<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	cout<<"non_remesh�ߓ_��="<<count_r<<endl;

	
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			//if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==FLUID &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID) conducter_num++;
			//if(NODE[i].material==CRUCIBLE) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num+100;///���m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����ƕς��Ȃ� //�Ȃ��������Ŏ~�܂�̂ŏ����]�T������Ĕz�������Ă݂�
	cout<<"pn="<<pn<<endl;
    int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
    int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 

	cout<<"pn="<<pn<<" ppn,npp�̊���U��-----";
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-3;
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-3;
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }
	cout<<"ok"<<endl;

	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=4*mat_w;//

	cout<<"calc_matrix�I��"<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].boundary_condition==0)
		{
			cout<<i<<endl;
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID)
				{
					flagi=ON;
					
				}
				if(NODE[i].material==CRUCIBLE)
				{
					cout<<"error! ��ڂ�boundary"<<endl;				
					flagi=ON;
				}
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////�z��m��
	cout<<"G,ROW,NUM�̊m�ۊJ�n"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_dir;//���s��ɓ���f�B���N���l
	complex<double> B_N;//��ɂ�����`��֐��R���̒l
	B_dir=complex<double> (0.0,0.0);
	B_N=complex<double> (0.0,0.0);

    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		//if(ELEM[je].material==COIL)
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aix�̍�
							    flag=1;
							}
						}
						if(flag==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aix�̍�
						    ROW[I][H]=J;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double>(B_temp,0); 
					    B[I-1]-=PHAT_A[A_X][NN]*B_N;
					}
				}

				/////Y����
				I=npp[N[n]]+1+1;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiy�̍�
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B_N=complex<double> (B_temp,0);
						B[I-1]-=PHAT_A[A_Y][NN]*B_N;
					}
				}
				/////*/

				/////Z����
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						J3=J2+1;		//Aiz�̍�
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiz�̍�
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aiz�̍�
						    ROW[I][H]=J3;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double> (B_temp,0);
					    B[I-1]-=PHAT_A[A_Z][NN]*B_N;
					}
				}
			}
		}
	}

	//////�Q�d�����v�Z
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aix�̍�
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//��������i�����݂��Ȃ���������(�����ł͂���Ȃ͂��Ȃ�����)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co;
							if(I==J) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����

					I=npp[N[n]]+1+1;/////Y����
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J4=J+3;			//�ӂ̍�
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ay
							B_temp=co;
							if(I==J2) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//�x�z��������X�������瓾����B�̒l

					/////Z����
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Az
							B_temp=co;
							if(I==J3) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
						}
					}////////*/
					//B[I-1]+=co*Sz;//�x�z��������X�������瓾����B�̒l

					/////��
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aix�̍�
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiy�̍�
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//Ay
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//Az
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//��
							B_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
							B_N=complex<double>(0,-B_temp);
							B[I-1]-=PHAT_V[NN]*B_N;

							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"��됔="<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);

	
	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];

	/*////
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	/////*/
	
	cout<<"���f�����̊���U��----";
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;

	cout<<"ok"<<endl;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D]/sqrt(2.0);//�����l�Ōv�Z
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	a.close();
	p.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}
	

	/*/////old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Je�ɂ͔g�̍���������

	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//�{���̑̐�

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//

	///old_A��̧�ُo��
	if(t==1)
	{
		cout<<"�����X�e�b�v�̃x�N�g���|�e���V�����L��"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	
	
	/////���̃X�e�b�v�ō쐬����NODE,ELEM��NODE_jw,ELEM_jw�ɋL��������B�����ŋ��߂��g���ƒx����ȍ~�̃X�e�b�v�ŗ��p����B�����߂邽��
	
	NODE_jw.clear();
	ELEM_jw.clear();
	/*////
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];
	//////*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT_A[D];
	delete [] PHAT_V;
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

//j�֖@5,���f���s����쐬 //�Q�X�e�b�v�ڈȍ~�A�ÓI�v�f�̃x�N�g���|�e���V�����l�͂P�X�e�b�v�ڂ̒l���f�B���N���l�Ƃ��ė��p����//ver2�ł͖��m���Ɋ֌W����z���vector��
void Avector3D_node_eddy2_jw5_ver2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta)
{ 
	/////
	cout<<"�Q�d�����l���ɂ��ꂽ�ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n�i�O���̈拫�E��ver2 �j"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�ߓ_��
	int NN_V=0;//�f�B���N���^���E�ߓ_���̂����A�Q�d�����l����ߓ_��
    int *dn=new int [node+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
	
	complex<double> *PHAT_A[3];
	for(int D=0;D<3;D++) PHAT_A[D]=new complex<double> [node+1];//A�̃f�B���N���^���l
	complex<double> *PHAT_V = new complex<double> [node+1];//V�̃f�B���N�����E�l
	

    ///�f�B���N���^���E��������
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//������
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	/////////1�X�e�b�v�ڂ̉�ǂݍ���
	cout<<"�ÓI�v�f�̃f�B���N���l�ǂݍ���"<<endl;
	//AR
	ifstream ar("old_AR.dat");
	if(!ar) cout<<"cannot open old_AR.dat"<<endl;
	ar.unsetf(ifstream::dec);
	ar.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ar>>AR[D][i];
	ar.close();

	//AI
	ifstream ai("old_AI.dat");
	if(!ai) cout<<"cannot open old_AI.dat"<<endl;
	ai.unsetf(ifstream::dec);
	ai.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ai>>AI[D][i];
	ai.close();

	//VR
	ifstream vr("old_VR.dat");
	if(!vr) cout<<"cannot open old_VR.dat"<<endl;
	vr.unsetf(ifstream::dec);
	vr.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vr>>VR[i];
	vr.close();

	//VI
	ifstream vi("old_VI.dat");
	if(!vi) cout<<"cannot open old_VI.dat"<<endl;
	vi.unsetf(ifstream::dec);
	vi.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vi>>VI[i];
	vi.close();
	/////
	cout<<"�f�B���N���l�ǂݍ��݊���"<<endl;		
	int count_r=0;
    for(int i=1;i<=node;i++)
    {
		/*///���̊֐��ɂ�static_dirichlet=OFF�̂Ƃ����ŗ��Ȃ�
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		////*/
		if(CON->get_static_dirichlet()==ON)
		{
			/*///
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			///*/
			////else
			
			{
				//if(NODE[i].remesh==OFF)//�����X�e�b�v�ȊO�̓����b�V�����Ȃ��ߓ_�ɂ��čŏ��ɋ��߂��x�N�g���|�e���V�������f�B���N���l�Ƃ��ė^����
				if(i<=node_sta)
				{
					count_r++;
					NODE[i].boundary_condition=3;//�ÓI�v�f�̋��E���������Ƃ�0����3�ɏ���������B�Œ苫�E(1,2)���ς��邪���Ȃ��͂�
					dn[i]=NN;
					PHAT_A[A_X][NN]=complex<double>(AR[A_X][i],AI[A_X][i]);
					PHAT_A[A_Y][NN]=complex<double>(AR[A_Y][i],AI[A_Y][i]);
					PHAT_A[A_Z][NN]=complex<double>(AR[A_Z][i],AI[A_Z][i]);
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;

					
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) NN_V++;
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				
				}
				else
				{
					dn[i]=node+1;
				}
			}
		}
    }//////////////
    cout<<"�ިظڐߓ_����"<<NN<<"�ިظڐ���"<<3*NN+NN_V<<" NN_V="<<NN_V<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
	cout<<"non_remesh�ߓ_��="<<count_r<<endl;

	
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���

	int conducter_num=0;//����(����)�ߓ_��
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			//if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==FLUID &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID) conducter_num++;
			//if(NODE[i].material==CRUCIBLE) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///���m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����ƕς��Ȃ� //�����グ���s���Ă�H//�v�f�`���Ɏ��s�����ߓ_�̐��~3��������Ă�
	cout<<"pn="<<pn<<endl;
    //int *ppn=new int [pn];		///�s���n�Ԗڂ̐ߓ_�͐ߓ_�ԍ�ppn[n]
	//int *npp=new int [node+1];	///�e�ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
	
	vector<int> ppn;

	vector<int> npp;
	npp.reserve(node+1); 
    
    int num=0; 

	cout<<"pn(�Q�l�l)="<<pn<<" ppn,npp�̊���U��-----";
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			//ppn[num]=i;
			//ppn[num+1]=-1;
			//ppn[num+2]=-2;
			ppn.push_back(i);
			ppn.push_back(-1);
			ppn.push_back(-2);
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }
	cout<<"ok"<<"num="<<num<<endl;

	pn=num;//num���^��pn�̂͂�

	///////
	////�s��̕��v�Z �e�s���Ƃɕ������߁A���������팸����
	int mat_w=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
	mat_w=4*mat_w;//

	cout<<"calc_matrix�I��"<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].boundary_condition==0)
		{
			
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID)
				{
					flagi=ON;
					
				}
				if(NODE[i].material==CRUCIBLE)
				{
					cout<<"error! ��ڂ�boundary"<<endl;				
					flagi=ON;
				}
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    /////*/

	////�z��m��
	cout<<"G,ROW,NUM�̊m�ۊJ�n"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_dir;//���s��ɓ���f�B���N���l
	complex<double> B_N;//��ɂ�����`��֐��R���̒l
	B_dir=complex<double> (0.0,0.0);
	B_N=complex<double> (0.0,0.0);

    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		//if(ELEM[je].material==COIL)
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X����

				int I=npp[N[n]]+1;
				//�e�a�͂����ł܂Ƃ߂Čv�Z����
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aix�̍�
							    flag=1;
							}
						}
						if(flag==0)//��������i�����݂��Ȃ���������
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aix�̍�
						    ROW[I][H]=J;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double>(B_temp,0); 
					    B[I-1]-=PHAT_A[A_X][NN]*B_N;
					}
				}

				/////Y����
				I=npp[N[n]]+1+1;/////Y����
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiy�̍�
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B_N=complex<double> (B_temp,0);
						B[I-1]-=PHAT_A[A_Y][NN]*B_N;
					}
				}
				/////*/

				/////Z����
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aix�̍�
						J2=J+1;			//Aiy�̍�
						J3=J2+1;		//Aiz�̍�
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiz�̍�
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aiz�̍�
						    ROW[I][H]=J3;
						}
					}
					else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double> (B_temp,0);
					    B[I-1]-=PHAT_A[A_Z][NN]*B_N;
					}
				}
			}
		}
	}

	//////�Q�d�����v�Z
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//�{���̑̐�
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//�Q�d�����̂����A�P�ï�ߑO���޸�����ݼ�ق��l�����邱�Ƃœ�����A���޸��
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			//double co2=delta6*sig[je];//   ��/(36V)
			//double co3=sig[je]*dt*delta6;

			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=(delta/20.0)*omega*sig[je];//��V/(20dt)�@���x���v�Z���邱�ƂɂȂ�̂ŌW��������	
			double co2=delta6*sig[je];//   ��/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X����

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aix�̍�
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//��������i�����݂��Ȃ���������(�����ł͂���Ȃ͂��Ȃ�����)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co;
							if(I==J) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//�x�z��������X�������瓾����B�̒l

					/////Y����

					I=npp[N[n]]+1+1;/////Y����
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J4=J+3;			//�ӂ̍�
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ay
							B_temp=co;
							if(I==J2) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//�x�z��������X�������瓾����B�̒l

					/////Z����
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Az
							B_temp=co;
							if(I==J3) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//��
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
						}
					}////////*/
					//B[I-1]+=co*Sz;//�x�z��������X�������瓾����B�̒l

					/////��
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aix�̍�
							J2=J+1;			//Aiy�̍�
							J3=J2+1;		//Aiz�̍�
							J4=J+3;			//�ӂ̍�
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aix�̍�
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])�Ƃ����ӂ���c[n]��O�Ɏ����Ă���ƁA�ł��؂�덷�̊֌W�Ŕ�Ώ̂ƔF�������̂ŁA��̎��Ɠ��l�ɍŌ�ɂ����Ă���
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiy�̍�
								if(J2==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aiz�̍�
								if(J3==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//�ӂ̍�
								if(J4==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//��������i�����݂��Ȃ���������
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3��(1/j)���������Ă���B1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///B�̌v�Z //����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//Ay
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//Az
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//��
							B_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
							B_N=complex<double>(0,-B_temp);
							B[I-1]-=PHAT_V[NN]*B_N;

							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"��됔="<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);
	
	//XX�Ɋi�[���ꂽ�����e�ߓ_�ɐU��
	//complex<double> *Ac[3];									//�ߓ_�ɂ������޸�����ݼ�فi���f���j
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];

	/*////
	double *AR[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//�x�N�g���|�e���V��������
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	/////*/
	
	cout<<"���f�����̊���U��----";
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;

	cout<<"ok"<<endl;
	
	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D]/sqrt(2.0);//�����l�Ōv�Z
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	a.close();
	p.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}
	

	/*/////old_A�𗱎q�ɓn��
	cout<<"old_A�𗱎q�ɋL��";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]
	int *count_e=new int[node+1];//�e���_�܂��̉Q�d�����������񂳂ꂽ�v�f�̐����i�[�B

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Je�ɂ͔g�̍���������

	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_loss_n[N[j]]+=Je_loss_e[je];
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			Je_loss_n[i]/=count_e[i];
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	delete [] count_e;
	//

	///old_A��̧�ُo��
	if(t==1)
	{
		cout<<"�����X�e�b�v�̃x�N�g���|�e���V�����L��"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	
	
	/////���̃X�e�b�v�ō쐬����NODE,ELEM��NODE_jw,ELEM_jw�ɋL��������B�����ŋ��߂��g���ƒx����ȍ~�̃X�e�b�v�ŗ��p����B�����߂邽��
	
	NODE_jw.clear();
	ELEM_jw.clear();
	/*////
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];
	//////*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT_A[D];
	delete [] PHAT_V;
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    //delete [] npp;
    //delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

///�޸�����ݼ�ٌv�Z�֐�(�ӗv�f�p)
void Avector3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP)
{
	cout<<"�޸�����ݼ�ٌv�Z�J�n---";

	double u0=4*PI*0.0000001;	//�^��̓�����
    double v0=1/u0;				//���C��R��
	double j0x,j0y,j0z;			//�d�����x[A/m^3]
	unsigned timeA=GetTickCount();

	//���΂̒�������������
	double MA=0;//=CON->get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON->get_magnet_B()*magnet_direction[A_X];
	double My=CON->get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON->get_magnet_B()*magnet_direction[A_Z];
	

	//////////////////////////////////////////////////////*/

	///���͕ǂɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].edge[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [side_num+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num+1���i�[
    double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"�ިظڐ�="<<NN;
	/////////////*/
    
	    
    int pn=side_num-NN;				///���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [side_num+1];	///�e�ӂ��s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]���܂ޗv�f��
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].edge[j];
			nume[edge]=nume[edge]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=side_num;i++)
	{	//�l����:����ӎ���̗v�f�Q���l����B�܂��A���̓��̈����Ԃɂ������i�K�ŁA�ӂ̐���4�B�Ȍ�A�v�f�������邲�Ƃ�3�ӂ�������B����Ď����ƂȂ�
		int width=4+3*(nume[i]-1);//�ꍇ�ɂ���Ă͂����菬�����l�ɂȂ邩������Ȃ��B���Ǒ������������m�ۂ���Ԃ�ɂ͖��Ȃ�
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

	
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
	for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//cout<<"�s��쐬�J�n ";
    for(int je=1;je<=nelm;je++)
    {   
		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		
		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON->get_RP();//�䓧����
		////�v�f��ظ��쐬�J�n
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
			if(SIDE[iside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				int I1=SIDE[iside].node[1];//iside���\������2�_
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jside���\������2�_
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				///B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///�s��̎��ۂ̍ő啝�����߂�
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	

	///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<" �s��쐬 ���F"<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

	///////////////////////�s��v�Z�J�n
	double *XX=new double [pn];//�s��̓����i�[
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

void Avector3D_edge_eddy(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:�v�f���\������ӂU�@EDGE.node:�ӂ�����_�Q�@
	/////
	int flageddy=OFF;
	cout<<"�Q�d�����l���ɂ��ꂽ�ӗv�f�ɂ���޸�����ݼ�ٌv�Z�J�n"<<endl;

	double u0=PI*4e-7;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	double j0x,j0y,j0z;					//�d�����x
	double Sx=0;
	double Sy=0;
	double Sz=0;
	//complex<double> j0x;
	//complex<double> j0y;
	//complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [nedge+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�nedge+1���i�[
	double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
	
	///��͗̈�̋��E�ɌŒ苫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"�ӗv�f�̌Œ苫�E�ݒ芮��"<<endl;

    ///�f�B���N���^���E��������
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }//////////////
    cout<<"�ިظڐ���"<<NN<<endl;
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���
	//�ӗv�f�ł��ӂɂ��Ă͐ߓ_�v�f��p����

	int conducter_num=0;//���̐ߓ_��
	int *conducter_flag=new int [node+1];//���̐ߓ_�Ȃ�ON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///�ӗv�f�̖��m��
	int pn;
	if(flageddy==ON) pn=pn_e+conducter_num;//�ӂ��܂߂��S���m��
	if(flageddy==OFF)pn=pn_e;

	int *ppn=new int [pn];		///�s���n�Ԗڂ͕Ӕԍ�ppn[n]
    int *npp=new int [nedge+node+1];	///�e��,���̐ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
   
	for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//���m��
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
	
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//��
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				ppn[num]=nedge+i;//
				npp[nedge+i]=num;//�ߓ_�v�f��npp�́A�ӗv�f��npp�����ׂĊi�[�������Ƃ̔z���p����
				num++;
			}
			else npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////�s��̕��v�Z �ӗv�f
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]���܂ޗv�f��
	int *width_edge= new int [pn+1];

	for(int i=1;i<=pn;i++) width_edge[i]=0;

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		width_edge[npp[i]+1]=width;	
	}

	
	cout<<"A-A�̍s�񕝌v�Z�I��"<<endl;

	////�s��̕��v�Z �ߓ_�v�f
	int mat_wn=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�

	for(int i=0;i<=pn;i++) width_mat[i]=0;
	mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_node�����܂�
	
	//�e�s�̐ߓ_�v�f�R���̔�됔�����߂� width_mat�Ɋi�[
	
	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	cout<<"nume_max="<<nume_max<<endl;
	
	if(flageddy==ON)
	{
		//A-��
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*5+1];
		for(int i=1;i<=pn;i++)//������
		{
			for(int j=1;j<=nume_max*5;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-�ӂ̍s�񕝌v�Z�I��"<<endl;
		
		//��-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]�̌v�Z

								}
							}		
						}
					}
				}
			}
		}

		cout<<"A-��,��-A�̍s�񕝌v�Z�I��"<<endl;
		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;


		//��-��
		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				width_mat[npp[i+nedge]+1]+=width_node[i];
			}
		}
		cout<<"��-�ӂ̍s�񕝌v�Z�I��"<<endl;
	}
	
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	

	////�z��m��
	//cout<<"�S�̍s��p�z��錾"<<endl;
	double **G=new double *[pn+1];///�S�̍s��
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��


    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0.0;
			ROW[i][j]=0;
		}
    }

		
	delete [] width_node;
	delete [] width_edge;
	delete [] width_mat;

    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0.0;//������
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//int J,J2,J3,flag,flag2,flag3;
	//double G_temp=0.0;
	//double B_temp=0.0;

	//�Î��ꍀ
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�d�S���W
		double Ys=0;
		double Zs=0;

		double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
			
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
    
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////�v�f��ظ��쐬�J�n
		//A-A
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
			if(EDGE[iside].boundary_condition==0)///���m��
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				//int I1=EDGE[iside].node[1];//iside���\������2�_
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jside���\������2�_
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///���m��
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6*v;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6*v;
							ROW[I][H]=J;
						}
					}
					////
					else //jside���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				///B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//�Q�d����
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//�Q�d�����̌v�Z�ΏۂƂȂ�v�f
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;

			double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//��A/��t
		
			double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			
			//cout<<"���̗v�f "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				if(EDGE[iside].boundary_condition==0)///���m��
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					//int I1=EDGE[iside].node[1];//iside���\������2�_
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jside���\������2�_
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///���m��
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									
									//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

									///////////
									G[I][h]+=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G[I][h]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G[I][h]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G[I][h]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G[I][h]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									G[I][h]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									//(Nk)y�E(Nu)y
									//////
									///////
									G[I][h]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G[I][h]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G[I][h]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G[I][h]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G[I][h]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G[I][h]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////
									////////
									//(Nk)z�E(Nu)z
									G[I][h]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G[I][h]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G[I][h]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G[I][h]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G[I][h]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G[I][h]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									//////

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"�N����Ȃ��͂�"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y�E(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z�E(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								//////
							    
								ROW[I][H]=J;
							}
							///B[I-1]���v�Z����
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
							int n=dn[jside];
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}///////////
					}///////
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-��
				double co2=sig[je]*delta*delta6*delta6*delta6;//��V*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						int I1=EDGE[iside].node[1];//iside���\������2�_
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G[I][h]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G[I][h]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G[I][h]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G[I][H]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G[I][H]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G[I][H]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									ROW[I][H]=J;
								}
							}
						}
					}
				}
				//cout<<"A-�� eddy"<<endl;

				//��-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"��-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
							//cout<<"��-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//iside���\������2�_
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G[I][h]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G[I][h]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G[I][h]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G[I][H]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G[I][H]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G[I][H]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
								}
								//B[I-1]�̌v�Z
							}
							else //jside���ިظڌ^���E�ߓ_�Ȃ�
							{
								int n=dn[iside];
								B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}///////////
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"��-A eddy"<<endl;

				//��-��
				double co3=sig[je]*dt*delta*delta6*delta6;
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G[I][h]+=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G[I][H]+=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								int NN=dn[N[j]];
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////
					}
				}
				//cout<<"��-�� eddy"<<endl;
			}
		}
	}

	///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"��[���v�f��"<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N----";
	check_matrix_symmetry(pn,NUM,ROW,G);
	cout<<"����"<<endl;;
	//���s��̒l�`�F�b�N
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/////////////
	///�o���h�������߂�
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//���O�p���̃o���h��
		m2=abs(ROW[i][NUM[i]]-i);//��O�p���̃o���h��
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//�傫������long int�ł�����Ȃ�
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"�ő�o���h��="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"���v�o���h��="<<band_sum<<endl;
	delete [] band;
	///////
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	double *XX=new double [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==4) COCG(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==5) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ӂɐU��

	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			A[i]=XX[n];
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			V[i]=XX[n];
		}
	}	

	delete [] XX;
	
	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/////
	//old_A�̍X�V 
	for(int i=1;i<=node;i++)
	{
		old_A[i]=A[i];
	}///

	//old_A�𗱎q�ɓn��
	//cout<<"old_A�𗱎q�ɋL��";
	//for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	//cout<<" ok"<<endl;
	
	///old_A��̧�ُo��
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[i]<<endl;
	g.close();
	
	//if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//�ӗv�f�ɂ��Q�d�����l���������C�x�N�g���|�e���V�����v�Z�֐�,jw�@
void Avector3D_edge_eddy_jw(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:�v�f���\������ӂU�@EDGE.node:�ӂ�����_�Q�@
	/////
	int flageddy=CON->get_A_phi();
	cout<<"�Q�d�����l���ɂ��ꂽ�ӗv�f�ɂ���޸�����ݼ�ٌv�Z�J�n(jw)"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [nedge+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�nedge+1���i�[
	double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
	
	///��͗̈�̋��E�ɌŒ苫�E������ݒ�

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"�ӗv�f�̌Œ苫�E�ݒ芮��"<<endl;


	///�f�B���N���^���E��������
	set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    
	/*//�f�B���N���^���E��������
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<NN<<endl;
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���
	//�ӗv�f�ł��ӂɂ��Ă͐ߓ_�v�f��p����

	int conducter_num=0;//���̐ߓ_��
	int *conducter_flag=new int [node+1];//���̐ߓ_�Ȃ�ON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///�ӗv�f�̖��m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����̂Q�{�ɂ͂Ȃ�Ȃ�
	int pn=pn_e;//�ӂ��܂߂��S���m��
	if(flageddy==ON) pn=pn_e+conducter_num;//�ӂ��܂߂��S���m��
	//int pn=pn_e;

	//int *ppn=new int [pn];		///�s���n�Ԗڂ͕Ӕԍ�ppn[n] �s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n](���̂̂݊i�[)
    //int *npp=new int [nedge+node+1];	///�e��,���̐ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    vector <int> ppn;
	vector <int> npp;
	npp.push_back(0);//0�Ԗڔz��𖄂߂�

	//npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//���m��
		{
			//ppn[num]=i;
			ppn.push_back(i);
			//npp[i]=num;
			npp.push_back(num);
			num++;
		}
		else npp.push_back(pn+1);   //npp[i]=pn+1;
    }
	cout<<"pn_e="<<pn_e<<" num="<<num<<endl;
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//��
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//�s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n]
				//npp[nedge+i]=num;//�ߓ_�v�f��npp�́A�ӗv�f��npp�����ׂĊi�[�������Ƃ̔z���p����
				npp.push_back(num);
				num++;
			}
			else npp.push_back(pn+1); //npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////�s��̕��v�Z �ӗv�f
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]���܂ޗv�f��
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		if(npp[i]<pn_e) width_edge[npp[i]+1]=width;	
		//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	//cout<<"nume_max="<<nume_max<<endl;

	
	//cout<<"A-A�̍s�񕝌v�Z�I��"<<endl;

	//�ߓ_�v�f�֌W�̍s�񕝕ϐ��̐錾
	int mat_wn=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *N_E_num=new int [node+1]; /////NODE[i]���܂ޗv�f��
	int *width_mat=new int [pn+1]; ///�e�s�̕�



	for(int i=1;i<=node;i++) 
	{
		width_node[i]=0;
		N_E_num[i]=0;
	}

	for(int i=0;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////�s��̕��v�Z
		
		for(int i=1;i<=nedge;i++) nume[i]=0;
		//�e�s�̐ߓ_�v�f�R���̔�됔�����߂� width_mat�Ɋi�[

		//A-��
		//int width_ap=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=6;j++)
				{
					int side=ELEM[je].edge[j];
					nume[side]=nume[side]+1;//
				}
			}
		}
		for(int i=1;i<=nedge;i++)
		{
			int width=0;//�s�̖��m����
			if(nume[i]>0) width=3+nume[i];//������nume[]��0�łȂ��Ƃ������Ƃ́A���̕ӂ͓��̗v�f�Ɋ܂܂�Ă���B�v�f��1�ǉ�����邽�сA���̕ӂɉe�����铱�̐ߓ_�͍ő�1������
			if(npp[i]<pn_e) width_mat[npp[i]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"A-�ӂ̍s�񕝌v�Z�I��"<<endl;

		//��-A
		//int width_pa=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					int N=ELEM[je].node[j];
					N_E_num[N]=N_E_num[N]+1;//
				}
			}
		}
		for(int i=1;i<=node;i++)
		{
			int width=0;//�s�̖��m����
			if(N_E_num[i]>0) width=1+3*N_E_num[i];//������nume[]��0�łȂ��Ƃ������Ƃ́A���̕ӂ͓��̗v�f�Ɋ܂܂�Ă���B�v�f��1�ǉ�����邽�сA���̐ߓ_�ɉe�����铱�̕ӂ͍ő�3������
			if(npp[i+nedge]>=pn_e && npp[i+nedge]<pn) width_mat[npp[i+nedge]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"��-A�̍s�񕝌v�Z�I��"<<endl;
		

		//��-��

		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_node�����܂�	

		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				width_mat[npp[i+nedge]+1]=width_mat[npp[i+nedge]+1]+width_node[i];
			}
		}
		//cout<<"��-�ӂ̍s�񕝌v�Z�I��"<<endl;
		/*////
		//A-��
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*4+1];
		for(int i=1;i<=pn;i++)//������
		{
			for(int j=1;j<=nume_max*4;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-�ӂ̍s�񕝌v�Z�I��"<<endl;
		
		//��-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]�̌v�Z

								}
							}		
						}
					}
				}
			}
		}



		cout<<"��-A�̍s�񕝌v�Z�I��"<<endl;

		

		
		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;
		cout<<"1"<<endl;
		*/////

		
	}
	
	//�s�񕝂̓���
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	
	for(int i=1;i<=pn;i++) 
	{
		if(width_mat[i]<=0)
		{
			cout<<"���m������0�ȉ��̍s���� i="<<i<<" "<<width_mat[i]<<endl;
			width_mat[i]=100;
		}
	}
	delete [] nume;
	
	////�z��m��
	cout<<"�S�̍s��p�z��錾"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
	int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��

    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	
	delete [] width_node;
	delete [] N_E_num;
	delete [] width_edge;
	delete [] width_mat;

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//�f�B���N���l�Ƃ��ĉ��s��ɉ��Z����鍀�̌`��֐���
	B_N=complex<double> (0.0,0.0);

	//�Î��ꍀ
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�d�S���W
		double Ys=0;
		double Zs=0;

		double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////�v�f��ظ��쐬�J�n
		//A-A(�Î���)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
			if(EDGE[iside].boundary_condition==0)///���m��
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				//int I1=EDGE[iside].node[1];//iside���\������2�_
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jside���\������2�_
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///���m��
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jside���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
						B[I-1]-=complex<double> (B_temp,0);
					}//////////*/
				}
				///�����d�����Ɋւ���B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//�Q�d����
	//if(flageddy==ON)
	{
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//�Q�d�����̌v�Z�ΏۂƂȂ�v�f
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;

			double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//��A/��t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"���̗v�f "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				if(EDGE[iside].boundary_condition==0)///���m��
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					//int I1=EDGE[iside].node[1];//iside���\������2�_
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jside���\������2�_
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///���m��
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									
									//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y�E(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z�E(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"�N����Ȃ��͂�"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y�E(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z�E(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]���v�Z����
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//jw�@�ɂ����ẮA�O�̃X�e�b�v�̒l�Ɋ�Â��������݂��Ȃ����߁A�����͕K�v�Ȃ�

							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
							int n=dn[jside];

							//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y�E(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z�E(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"�f�B���N���l�����m���̂͂��̍s�ԍ��ɑ���"<<endl;
								B_N=complex<double>(0,PHAT[n]*B_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�//�N����Ȃ��H
							}
							else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//���̋L�q�͑��v�H
							B[I-1]-=B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-��
				double co2=sig[je]*delta*delta6*delta6*delta6;//��V*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						int I1=EDGE[iside].node[1];//iside���\������2�_
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//���̂łȂ��Ƃ������Ƃ̓�=0�Ȃ̂ŁA������A�̍��̂悤�ȃf�B���N���l�Ɋւ��čl����K�v�͂Ȃ�
						}
					}
				}
				//cout<<"A-�� eddy"<<endl;

				//��-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"��-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
							//cout<<"��-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//iside���\������2�_
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]�̌v�Z
							}
							else //jside���ިظڌ^���E�ߓ_�Ȃ�
							{
								int n=dn[iside];
								B_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHAT�͕��f�f�B���N���ɔ����Ȃ����ׂ�
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"��-A eddy"<<endl;

				//��-��
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//�����1/j���������Ă��邽�߁A���f���Ƃ��Ă͕�����ς��ċ������ɑ�����邱�ƂɂȂ�(1/j=-j)
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								int NN=dn[N[j]];
								B_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}//////*/
					}
				}
				//cout<<"��-�� eddy"<<endl;
			}
		}
	}
	}

	///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"��[���v�f��"<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	//���s��̒l�`�F�b�N
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ӂɐU��
	//cout<<"���f�����̊���U��"<<endl;

	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int i=1;i<=nedge;i++)
	{
		
		AR[i]=0;
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	ofstream ar("old_AR_e.dat");
	ofstream ai("old_AI_e.dat");
	ofstream vr("old_VR_e.dat");
	ofstream vi("old_VI_e.dat");
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	
	delete [] XX;

	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	ofstream a("Am_e.dat");
	ofstream p("phi_e.dat");
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			Am=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am*cos(omega*TIME+phi);
		}
		a<<Am<<endl;
		p<<phi<<endl;
		ar<<AR[i]<<endl;
		ai<<AI[i]<<endl;
	}
	for(int i=1;i<=node;i++)
	{
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}

	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();
	//cout<<"���f���������ϊ�����"<<endl;
	/*///
	for(int i=1;i<=node;i++)
	{
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);
	}
	*/
	
	
	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//�{���̑̐�

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//�ӗv�f�ɂ��Q�d�����l���������C�x�N�g���|�e���V�����v�Z�֐�,jw�@,�ߓ_�v�f�ł���ӂ���2���v�f��
void Avector3D_edge_eddy_jw_with_parabolic_node_element(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:�v�f���\������ӂU�@EDGE.node:�ӂ�����_�Q�@
	/////
	int flageddy=CON->get_A_phi();
	cout<<"�Q�d�����l���ɂ��ꂽ1����2���ߓ_�v�f�ɂ���޸�����ݼ�ٌv�Z�J�n(jw)"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [nedge+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�nedge+1���i�[
	double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
	
	///��͗̈�̋��E�ɌŒ苫�E������ݒ�

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"�ӗv�f�̌Œ苫�E�ݒ芮��"<<endl;


	///�f�B���N���^���E��������
	set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    
	/*//�f�B���N���^���E��������
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }/////////////*/
    cout<<"�ިظڐ���"<<NN<<endl;
	
	//2���v�f�̂��߁A���̕ӂ𐔂���
	int *conduct_edge_flag=new int [nedge+1];//���̕ӂȂ�ON
	for(int i=0;i<=nedge;i++) conduct_edge_flag[i]=OFF;

	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=6;j++)
			{
				int side=ELEM[je].edge[j];
				conduct_edge_flag[side]=ON;//���̗v�f�Ɋ܂܂��ӂ��`�F�b�N����
			}
		}
	}
	int para_node_num=0;
	
	for(int i=1;i<=nedge;i++)
	{
		EDGE[i].para_node_num=0;//������
		if(conduct_edge_flag[i]==ON)
		{
			para_node_num++;
			EDGE[i].para_node_num=node+para_node_num;
		}
	}

	cout<<"���̕Ӑ�="<<para_node_num<<endl;

	int node2=node+para_node_num;//2���ߓ_���܂߂��S�ߓ_��

	
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���
	//�ӗv�f�ł��ӂɂ��Ă͐ߓ_�v�f��p����

	int conducter_num=0;//���̐ߓ_�� ���Ƃ��Ƃ�FEM���b�V���ő��݂��Ă����ߓ_�̂ݐ����Ă���
	int *conducter_flag=new int [node2+1];//���̐ߓ_�Ȃ�ON
	for(int i=1;i<=node2;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
		else if(CON->get_Je_crucible()==-1)
		{
		}
	}
	for(int i=node+1;i<=node2;i++)//2���v�f�Œǉ������ߓ_
	{
		conducter_flag[i]=ON;
	}


	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	cout<<"conducter_num2="<<conducter_num+para_node_num<<endl;
	int pn_e=nedge-NN;///�ӗv�f�̖��m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����̂Q�{�ɂ͂Ȃ�Ȃ�
	int pn=pn_e;//�ӂ��܂߂��S���m��
	if(flageddy==ON) pn=pn_e+conducter_num+para_node_num;//�ӂ��܂߂��S���m��
	//int pn=pn_e;

	//int *ppn=new int [pn];		///�s���n�Ԗڂ͕Ӕԍ�ppn[n] �s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n](���̂̂݊i�[)
    //int *npp=new int [nedge+node+1];	///�e��,���̐ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    vector <int> ppn;
	vector <int> npp;
	npp.push_back(0);//0�Ԗڔz��𖄂߂�

	//npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//���m��
		{
			//ppn[num]=i;
			ppn.push_back(i);
			//npp[i]=num;
			npp.push_back(num);
			num++;
		}
		else npp.push_back(pn+1);   //npp[i]=pn+1;
    }
	cout<<"pn_e="<<pn_e<<" num="<<num<<endl;
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node2;i++)//��
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//�s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n]
				//npp[nedge+i]=num;//�ߓ_�v�f��npp�́A�ӗv�f��npp�����ׂĊi�[�������Ƃ̔z���p����
				npp.push_back(num);
				num++;
			}
			else npp.push_back(pn+1); //npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////�s��̕��v�Z �ӗv�f
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]���܂ޗv�f��
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		if(npp[i]<pn_e) width_edge[npp[i]+1]=width;	
		//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	//cout<<"nume_max="<<nume_max<<endl;

	
	//cout<<"A-A�̍s�񕝌v�Z�I��"<<endl;

	//�ߓ_�v�f�֌W�̍s�񕝕ϐ��̐錾
	int mat_wn=0;
	int *width_node=new int [node2+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *N_E_num=new int [node2+1]; /////NODE[i]���܂ޗv�f��
	int *width_mat=new int [pn+1]; ///�e�s�̕�



	for(int i=1;i<=node2;i++) 
	{
		width_node[i]=0;
		N_E_num[i]=0;
	}

	for(int i=1;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////�s��̕��v�Z
		
	
		//�e�s�̐ߓ_�v�f�R���̔�됔�����߂� width_mat�Ɋi�[

		//A-��
		//int width_ap=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=6;j++)
				{
					int side=ELEM[je].edge[j];
					nume[side]=nume[side]+1;//
				}
			}
		}
		for(int i=1;i<=nedge;i++)
		{
			int width=0;//�s�̖��m����
			//if(nume[i]>0) width=3+nume[i];//������nume[]��0�łȂ��Ƃ������Ƃ́A���̕ӂ͓��̗v�f�Ɋ܂܂�Ă���B�v�f��1�ǉ�����邽�сA���̕ӂɉe�����铱�̐ߓ_�͍ő�1������
			if(nume[i]>0) width=6+4*nume[i];//������nume[]��0�łȂ��Ƃ������Ƃ́A���̕ӂ͓��̗v�f�Ɋ܂܂�Ă���B�v�f��1�ǉ�����邽�сA���̕ӂɉe�����铱�̐ߓ_��2���ߓ_���܂ߍő�4������
			if(npp[i]<pn_e) width_mat[npp[i]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"A-�ӂ̍s�񕝌v�Z�I��"<<endl;

		//��-A
		//int width_pa=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					int N=ELEM[je].node[j];
					N_E_num[N]=N_E_num[N]+1;//
				}
				for(int j=1;j<=6;j++)//2���ߓ_
				{
					int side=ELEM[je].edge[j];
					int N=EDGE[side].para_node_num;
					if(N>0) N_E_num[N]=N_E_num[N]+1;//
				}
			}
		}
		for(int i=1;i<=node2;i++)
		{
			int width=0;//�s�̖��m����
			if(N_E_num[i]>0) width=1+3*N_E_num[i];//������nume[]��0�łȂ��Ƃ������Ƃ́A���̕ӂ͓��̗v�f�Ɋ܂܂�Ă���B�v�f��1�ǉ�����邽�сA���̐ߓ_�ɉe�����铱�̕ӂ͍ő�3������
			if(npp[i+nedge]>=pn_e) width_mat[npp[i+nedge]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"��-A�̍s�񕝌v�Z�I��"<<endl;
		

		//��-��

		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_node�����܂�	

		for(int i=node+1;i<=node2;i++) width_node[i]=2+1;//2���ߓ_�͕ӂ̒��_�Ȃ̂ŁA�������g���܂ߗׂ荇���̂�3�_�̂͂��H

		for(int i=1;i<=node2;i++)
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				width_mat[npp[i+nedge]+1]+=width_node[i];
			}
		}
		//cout<<"��-�ӂ̍s�񕝌v�Z�I��"<<endl;
		/////
	}
		
	
	//�s�񕝂̓���
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	
	////�z��m��
	cout<<"�S�̍s��p�z��錾"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��

    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	delete [] width_node;
	delete [] N_E_num;
	delete [] width_edge;
	delete [] width_mat;
	

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//�f�B���N���l�Ƃ��ĉ��s��ɉ��Z����鍀�̌`��֐���
	B_N=complex<double> (0.0,0.0);

	//�Î��ꍀ
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�d�S���W
		double Ys=0;
		double Zs=0;

		double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////�v�f��ظ��쐬�J�n
		//A-A(�Î���)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
			if(EDGE[iside].boundary_condition==0)///���m��
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				//int I1=EDGE[iside].node[1];//iside���\������2�_
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jside���\������2�_
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///���m��
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jside���ިظڌ^���E�ߓ_�Ȃ�
					{
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
						B[I-1]-=complex<double> (B_temp,0);
					}//////////*/
				}
				///�����d�����Ɋւ���B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//�Q�d����
	//if(flageddy==ON)
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//�Q�d�����̌v�Z�ΏۂƂȂ�v�f
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;

			double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//��A/��t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"���̗v�f "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				if(EDGE[iside].boundary_condition==0)///���m��
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					//int I1=EDGE[iside].node[1];//iside���\������2�_
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jside���\������2�_
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///���m��
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									
									//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y�E(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z�E(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"A�̎��Ԕ������v�Z�ŐV���ɔz��m�ہH"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y�E(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z�E(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]���v�Z����
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//jw�@�ɂ����ẮA�O�̃X�e�b�v�̒l�Ɋ�Â��������݂��Ȃ����߁A�����͕K�v�Ȃ�

							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
							int n=dn[jside];

							//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y�E(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z�E(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"�f�B���N���l�����m���̂͂��̍s�ԍ��ɑ���"<<endl;
								B_N=complex<double>(0,PHAT[n]*B_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�//�N����Ȃ��H
							}
							else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//���̋L�q�͑��v�H
							B[I-1]-=B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-��
				double co2=sig[je]*delta*delta6*delta6*delta6;//��V*(1/6V)^3
				double co4=sig[je]*delta6*delta6*delta6*delta6*delta*4;//(1/6V)^2 * 1/9V^2
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						int I1=EDGE[iside].node[1];//iside���\������2�_
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNi=Ni*(2*Ni-1)
										//gradNi��ꍀ 2Ni^2 
										//(Nk)x * ��N��x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*c[j];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*c[j];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*c[j];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*d[j];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*d[j];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*d[j];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*e[j];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*e[j];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*e[j];

										///////////////										
										//gradNi��� -Ni
										G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j]*(-1);
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j]*(-1);
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j]*(-1);


										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNi=Ni*(2*Ni-1)
									//gradNi��ꍀ 2Ni^2 
									//(Nk)x * ��N��x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*c[j];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*c[j];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*c[j];
										
									///////////

									//(Nk)y * ��N��y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*d[j];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*d[j];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*d[j];

									//(Nk)z * ��N��z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*e[j];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*e[j];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*e[j];

									///////////////										
									//gradNi��� -Ni
									G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j]*(-1);
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j]*(-1);
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j]*(-1);
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//���̂łȂ��Ƃ������Ƃ̓�=0�Ȃ̂ŁA������A�̍��̂悤�ȃf�B���N���l�Ɋւ��čl����K�v�͂Ȃ�
						}

						for(int j=1;j<=6;j++)//�ǉ�2���ߓ_
						{	
							int side=ELEM[je].edge[j];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int J=npp[N+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;	

								int I1=EDGE[iside].node[1];//iside���\������2�_
								int I2=EDGE[iside].node[2];
								int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
								int u2=table[j][2];
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNij=4*NiNj
										//gradNi��ꍀ 2Ni^2 
										//(Nk)x * ��N��x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///Ni��Nj���t�]�����̍�
										////
										//(Nk)x * ��N��x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////										

										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNij=4*NiNj
									//gradNi��ꍀ 2Ni^2 
									//(Nk)x * ��N��x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
									///////////

									//(Nk)y * ��N��y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

									//(Nk)z * ��N��z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

									//
									///Ni��Nj���t�]�����̍�
									////
									//(Nk)x * ��N��x
									G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
									///////////

									//(Nk)y * ��N��y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

									//(Nk)z * ��N��z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//���̂łȂ��Ƃ������Ƃ̓�=0�Ȃ̂ŁA������A�̍��̂悤�ȃf�B���N���l�Ɋւ��čl����K�v�͂Ȃ�
						}
					}
				}
				//cout<<"A-�� eddy"<<endl;

				//��-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"��-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
							//cout<<"��-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//iside���\������2�_
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNi=Ni*(2*Ni-1)
										//gradNi��ꍀ 2Ni^2 
										//(Nk)x * ��N��x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

										///////////////										
										//gradNi��� -Ni
										G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNi=Ni*(2*Ni-1)
									//gradNi��ꍀ 2Ni^2 
									//(Nk)x * ��N��x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
									///////////

									//(Nk)y * ��N��y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

									//(Nk)z * ��N��z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

									///////////////										
									//gradNi��� -Ni
									G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]�̌v�Z
							}
							else //jside���ިظڌ^���E�ߓ_�Ȃ�
							{
								int n=dn[iside];
								//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNi=Ni*(2*Ni-1)
								//gradNi��ꍀ 2Ni^2 
								//(Nk)x * ��N��x
								B_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
								///////////

								//(Nk)y * ��N��y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

								//(Nk)z * ��N��z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

								///////////////										
								//gradNi��� -Ni
								B_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHAT�͕��f�f�B���N���ɔ����Ȃ����ׂ�
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				for(int i=1;i<=6;i++)//�ǉ�2���ߓ_
				{	
					int side=ELEM[je].edge[i];
					int N=EDGE[side].para_node_num;
					if(conducter_flag[N]==ON)
					{
						int I=npp[N+nedge]+1;
						int I1=EDGE[side].node[1];//iside���\������2�_
						int I2=EDGE[side].node[2];
						int u1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[i][2];
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
							//cout<<"��-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//iside���\������2�_
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNij=4*NiNj
										//gradNi��ꍀ 2Ni^2 
										//(Nk)x * ��N��x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///Ni��Nj���t�]�����̍�
										////
										//(Nk)x * ��N��x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////						

										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNij=4*NiNj
										//gradNi��ꍀ 2Ni^2 
										//(Nk)x * ��N��x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///Ni��Nj���t�]�����̍�
										////
										//(Nk)x * ��N��x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ��N��y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ��N��z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////					

									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]�̌v�Z
							}
							else //jside���ިظڌ^���E�ߓ_�Ȃ�
							{
								int n=dn[iside];
								//gradNi���A�񎟂̌`��֐����̗p���邱�Ƃɂ���đ傫���ω�����BNij=4*NiNj
								//gradNi��ꍀ 2Ni^2 
								//(Nk)x * ��N��x
								B_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
								///////////

								//(Nk)y * ��N��y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

								//(Nk)z * ��N��z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

								//
								///Ni��Nj���t�]�����̍�
								////
								//(Nk)x * ��N��x
								B_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
								///////////

								//(Nk)y * ��N��y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

								//(Nk)z * ��N��z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
								///////////////				

								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHAT�͕��f�f�B���N���ɔ����Ȃ����ׂ�
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"��-A eddy"<<endl;

				//��-��
				/////
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//�����1/j���������Ă��邽�߁A���f���Ƃ��Ă͕�����ς��ċ������ɑ�����邱�ƂɂȂ�(1/j=-j)
				double co5=sig[je]*delta*delta6*delta6*delta6*delta6*16/omega;
				double co6=sig[je]*delta*delta6*delta6*delta6*4/omega;
				
				for(int i=1;i<=4;i++)//���Ƃ̐ߓ_
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{

										G_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										G_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										G_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								int NN=dn[N[j]];
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
								B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////

						for(int j=1;j<=6;j++)//�񎟐ߓ_
						{
							int side=ELEM[je].edge[i];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int flag=0;
								int J=npp[N+nedge]+1;
								int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
								int u2=table[j][2];
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{

										G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								//int NN=dn[N[j]];
								cout<<"�ӂ����m���łȂ�?"<<endl;
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
								B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}
					}
				}

				for(int i=1;i<=6;i++)//�ǉ��񎟐ߓ_
				{			
					int side=ELEM[je].edge[i];
					int Ni=EDGE[side].para_node_num;	
					int I=npp[Ni+nedge]+1;
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[Ni]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{

										G_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
										G_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));

										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
									G_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));

									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								//cout<<"�ǉ�2���v�f�ߓ_���f�B���N���H"<<endl;
								int NN=dn[N[j]];
								B_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

								B_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
								B_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////

						for(int j=1;j<=6;j++)//�񎟐ߓ_
						{
							int side=ELEM[je].edge[i];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int flag=0;
								int J=npp[N+nedge]+1;
								int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
								int u2=table[j][2];
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{

										G_temp = co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(b[k2]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(c[k2]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(d[k2]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(e[k2]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(b[k2]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(c[k2]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(d[k2]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(e[k2]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(b[k1]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(c[k1]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(d[k1]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(e[k1]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(b[k1]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(c[k1]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(d[k1]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(e[k1]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										
										
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								int NN=dn[N];
								cout<<"�ӂ��f�B���N���l?"<<endl;
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}
					}
				}
				////*/
				//cout<<"��-�� eddy"<<endl;
			}
		}
	}
	

	///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"��[���v�f��"<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	cout<<"����"<<endl;;
	//���s��̒l�`�F�b�N
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ӂɐU��
	//cout<<"���f�����̊���U��"<<endl;

	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int i=1;i<=nedge;i++)
	{
		
		AR[i]=0;
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	ofstream ar("old_AR_e.dat");
	ofstream ai("old_AI_e.dat");
	ofstream vr("old_VR_e.dat");
	ofstream vi("old_VI_e.dat");
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	
	delete [] XX;

	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	ofstream a("Am_e.dat");
	ofstream p("phi_e.dat");
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			Am=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am*cos(omega*TIME+phi);
		}
		a<<Am<<endl;
		p<<phi<<endl;
		ar<<AR[i]<<endl;
		ai<<AI[i]<<endl;
	}
	for(int i=1;i<=node;i++)
	{
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}

	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();
	//cout<<"���f���������ϊ�����"<<endl;
	/*///
	for(int i=1;i<=node;i++)
	{
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);
	}
	*/
	
	
	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//�{���̑̐�

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	delete [] conduct_edge_flag;
	
	//////
	/////*/
}


//�ӗv�f�ɂ��Q�d�����l���������C�x�N�g���|�e���V�����v�Z�֐�,jw�@ �ÓI�v�f�̃f�B���N���l�ǂݍ���
void Avector3D_edge_eddy_jw2(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega, int node_sta, vector<edge3D> &static_EDGE,double *Am,double *phi)
{
	//ELEM.edge:�v�f���\������ӂU�@EDGE.node:�ӂ�����_�Q�@
	/////
	int flageddy=OFF;//A-�Ӗ@�ɂ��邩���Ȃ���

	cout<<"�Q�d�����l���ɂ��ꂽ�ӗv�f�ɂ���޸�����ݼ�ٌv�Z�J�n(jw_static�ǂݍ���)"<<endl;

	double u0=4*PI*0.0000001;			//��C�̓�����
    double v0=1/u0;						//���C��R��
	//double j0x,j0y,j0z;					//�d�����x
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//�d�C�`����
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	int NN=0;//�f�B���N���^���E�Ӑ�
	int NN_V=0;//�f�B���N���^���E�ߓ_��
    int *dn=new int [nedge+1]; //�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�nedge+1���i�[
	double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
	
	complex<double> *PHAT_A=new complex<double> [nedge+1];//A�̃f�B���N���^���l
	complex<double> *PHAT_V = new complex<double> [node+1];//V�̃f�B���N�����E�l

	///��͗̈�̋��E�ɌŒ苫�E������ݒ�

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
					int p=ELEM[i].node[j];//���E�O�p�`�ɑ����Ȃ��ߓ_
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//�ߓ_p���܂܂Ȃ��ӂ͋��E��
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"�ӗv�f�̌Œ苫�E�ݒ芮��"<<endl;


	///�f�B���N���^���E��������
	//set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    

	//�ÓI�v�f�̃f�B���N���l�ǂݍ���
	double *AR = new double [nedge+1];//�x�N�g���|�e���V��������
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//������
	for(int i=1;i<=nedge;i++)
	{
		//AR[i]=A[i];//���������0�łȂ����E����(��l����Ȃ�)�������Ă���ꍇ���l���A�����œ���Ă���
		AR[i]=0;	
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	if(CON->get_static_dirichlet()==ON)
	{
		/////////1�X�e�b�v�ڂ̉�ǂݍ���
		cout<<"�ÓI�v�f�̃f�B���N���l�ǂݍ���"<<endl;

		//AR
		ifstream ar("old_AR_e.dat");
		if(!ar) cout<<"cannot open old_AR_e.dat"<<endl;
		ar.unsetf(ifstream::dec);
		//ar.setf(ifstream::skipws);
		for(int i=1;i<=nedge;i++) ar>>AR[i];
		ar.close();

		//AI
		ifstream ai("old_AI_e.dat");
		if(!ai) cout<<"cannot open old_AI_e.dat"<<endl;
		ai.unsetf(ifstream::dec);
		//ai.setf(ifstream::skipws);
		for(int i=1;i<=nedge;i++) ai>>AI[i];
		ai.close();

		//VR
		ifstream vr("old_VR_e.dat");
		if(!vr) cout<<"cannot open old_VR.dat"<<endl;
		//vr.unsetf(ifstream::dec);
		//vr.setf(ifstream::skipws);
		for(int i=1;i<=node;i++) vr>>VR[i];
		vr.close();

		//VI
		ifstream vi("old_VI_e.dat");
		if(!vi) cout<<"cannot open old_VI.dat"<<endl;
		//vi.unsetf(ifstream::dec);
		//vi.setf(ifstream::skipws);
		for(int i=1;i<=node;i++) vi>>VI[i];
		vi.close();
		/////
		cout<<"�f�B���N���l�ǂݍ��݊���"<<endl;		
		int count_r=0;

		
		////

		for(int i=1;i<=nedge;i++)
		{
			/*///���̊֐��ɂ�static_dirichlet=OFF�̂Ƃ����ŗ��Ȃ�
			if(CON->get_static_dirichlet()==OFF)
			{
				if(NODE[i].boundary_condition==2)
				{    
					dn[i]=NN;
					PHAT[A_X][NN]=0;
					PHAT[A_Y][NN]=0;
					PHAT[A_Z][NN]=0;
					A[A_X][i]=0;
					A[A_Y][i]=0;
					A[A_Z][i]=0;
					NN++;
				}
				else if(NODE[i].boundary_condition==1)
				{    
					dn[i]=NN;
					PHAT[A_X][NN]=0;
					PHAT[A_Y][NN]=0;
					PHAT[A_Z][NN]=0;
					A[A_X][i]=0;
					A[A_Y][i]=0;
					A[A_Z][i]=0;
					NN++;
				}
				else
				{
					dn[i]=node+1;
				}
			}
			////*/

			
			{
				//if(NODE[i].remesh==OFF)//�����X�e�b�v�ȊO�̓����b�V�����Ȃ��ߓ_�ɂ��čŏ��ɋ��߂��x�N�g���|�e���V�������f�B���N���l�Ƃ��ė^����
				if(EDGE[i].stat==ON)
				{
					count_r++;
					EDGE[i].boundary_condition=3;//�ÓI�v�f�̋��E���������Ƃ�0����3�ɏ���������B�Œ苫�E(1,2)���ς��邪���Ȃ��͂�
					dn[i]=NN;
					PHAT_A[NN]=complex<double>(AR[EDGE[i].static_num],AI[EDGE[i].static_num]);
					//PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;
				}
				else
				{
					dn[i]=nedge+1;
				}
			}
			//*/
		}//////////////
		
			for(int i=1;i<=node;i++)
			{
				if(i<=node_sta)
				{
					NODE[i].boundary_condition=3;
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!���̂��ÓI�ߓ_"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0)
						{
							NN_V++;
						}
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				}
				else
				{
					//dn[i]=node+1;
				}
			}
		
		cout<<"�ިظڕӐ���"<<NN<<"�ިظڐ���"<<NN+NN_V<<" NN_V="<<NN_V<<endl;//���ۂɌv�Z���Ȃ��Ă����ިظڌ^�ߓ_����3*NN�ł���
		//cout<<"non_remesh�Ӑ�="<<count_r<<endl;
	}

    //cout<<"�ިظڐ���"<<NN<<endl;
	    
	///�`�Ӗ@�ɂ͓���(����)�ߓ_���̐������d�ʂɊւ��関�m�������݂���B����𐔂���
	//�ӗv�f�ł��ӂɂ��Ă͐ߓ_�v�f��p����

	int conducter_num=0;//���̐ߓ_��
	int *conducter_flag=new int [node+1];//���̐ߓ_�Ȃ�ON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///�ӗv�f�̖��m�� ���f���̂܂܊i�[����̂Ŏ��ԍ����̂Q�{�ɂ͂Ȃ�Ȃ�
	int pn=pn_e;//�ӂ��܂߂��S���m��
	if(flageddy==ON) pn=pn_e+conducter_num;//�ӂ��܂߂��S���m��
	//int pn=pn_e;

	//int *ppn=new int [pn];		///�s���n�Ԗڂ͕Ӕԍ�ppn[n] �s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n](���̂̂݊i�[)
    //int *npp=new int [nedge+node+1];	///�e��,���̐ߓ_���s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    vector <int> ppn;
	vector <int> npp;
	

	npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//���m��
		{
			//ppn[num]=i;
			ppn.push_back(i);
			npp[i]=num;
			//npp.push_back(num);
			num++;
		}
		else npp[i]=pn+1;// npp.push_back(pn+1);   //
    }
	
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//��
		{
			if(NODE[i].boundary_condition==0 && conducter_flag[i]==ON)//���̐ߓ_
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//�s���n(�������An�͕Ӑ����傫��)�Ԗڂ͐ߓ_�ԍ�ppn[n]
				npp[nedge+i]=num;//�ߓ_�v�f��npp�́A�ӗv�f��npp�����ׂĊi�[�������Ƃ̔z���p����
				//npp.push_back(num);
				num++;
			}
			else npp[nedge+i]=pn+1;// npp.push_back(pn+1); //
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////�s��̕��v�Z �ӗv�f
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]���܂ޗv�f��
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]�����Ƃ܂���
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			int width=7+3*(nume[i]-2);
			if(width>mat_we) mat_we=width;
			if(npp[i]<pn) width_edge[npp[i]+1]=width;
		}
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	cout<<"nume_max="<<nume_max<<endl;

	
	cout<<"A-A�̍s�񕝌v�Z�I��"<<endl;

	//�ߓ_�v�f�֌W�̍s�񕝕ϐ��̐錾
	int mat_wn=0;
	int *width_node=new int [node+1]; ///�e�ߓ_�ׂ̗荇���ߓ_�̐��{�P
	int *width_mat=new int [pn+1]; ///�e�s�̕�

	for(int i=0;i<=node;i++) 
	{
		width_node[i]=0;
	}

	for(int i=0;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////�s��̕��v�Z �Q�d����
		//��-��
		
		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_node�����܂�
	
		//�e�s�̐ߓ_�v�f�R���̔�됔�����߂� width_mat�Ɋi�[
	
	
	
	
		//A-��
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*4+1];
		for(int i=1;i<=pn;i++)//������
		{
			for(int j=1;j<=nume_max*4;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-�ӂ̍s�񕝌v�Z�I��"<<endl;
		
		//��-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]�̌v�Z

								}
							}		
						}
					}
				}
			}
		}



		cout<<"A-��,��-A�̍s�񕝌v�Z�I��"<<endl;

		//��-��
		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//���̐ߓ_
			{
				width_mat[npp[i]+1]+=width_node[i];
			}
		}
		cout<<"��-�ӂ̍s�񕝌v�Z�I��"<<endl;

		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;
	}
	
	//�s�񕝂̓���
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	
	
	////�z��m��
	//cout<<"�S�̍s��p�z��錾"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///�S�̍s��
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��


    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	
	delete [] width_node;
	delete [] width_edge;
	delete [] width_mat;

    complex<double> *B=new complex<double> [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//������
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"�S�̍s��쐬�J�n"<<endl;
    /////////�S�̍s����쐬����
    unsigned int time=GetTickCount();
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//�f�B���N���l�Ƃ��ĉ��s��ɉ��Z����鍀�̌`��֐���
	B_N=complex<double> (0.0,0.0);

	//�Î��ꍀ
    for(int je=1;je<=nelm;je++)
    {
		//if(t>1 && je%1000==0) cout<<je<<endl;
		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�d�S���W
		double Ys=0;
		double Zs=0;

		double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
	
		double delta6=ELEM[je].volume;//�̐ς�6�{
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1�X�e�b�v�ځA���Ȃ킿j0���ő�ɂȂ��Ă���Ƃ��̂ݐ��藧�B���Ƃ�current�̒�`���C�����邱��
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///�䓧����
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////�v�f��ظ��쐬�J�n
		//A-A(�Î���)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
			if(EDGE[iside].boundary_condition==0)///���m��
			{   
				int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
				//int I1=EDGE[iside].node[1];//iside���\������2�_
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jside���\������2�_
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///���m��
					{   
						int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jside���ިظڌ^���E�ӂȂ�
					{
						//if(EDGE[jside].boundary_condition!=3) cout<<"bc��3�ȊO�H"<<endl;
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
						B_N=complex<double>(B_temp,0); 
						B[I-1]-=PHAT_A[n]*B_N;
					}//////////*/
				}
				///�����d�����Ɋւ���B[I-1]���v�Z����
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
					
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}
	cout<<"�Î��ꍀ�̍쐬����"<<endl;

	//�Q�d����
	//if(flageddy==ON)
	{
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//�Q�d�����̌v�Z�ΏۂƂȂ�v�f
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;

			double XXs=0;//���a�@x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//�Ϙa�@x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//��A/��t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//���֖@�̏ꍇ�A1/dt�����ւɒu��������΂悢�B���͊e�����쐬���ɂ��܂����킹��
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"���̗v�f "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				if(EDGE[iside].boundary_condition==0)///���m��
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					//int I1=EDGE[iside].node[1];//iside���\������2�_
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jside���\������2�_
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///���m��
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									
									//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y�E(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z�E(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"�N����Ȃ��͂�"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y�E(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z�E(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]���v�Z����
							//����N[m]���Œ苫�E�ゾ�����ꍇ��Sx�̌v�Z�ɑg�ݍ��܂Ȃ��Ƃ����Ȃ����A�f�B���N���l�͂O�Ɖ��肵�Čv�Z�����O
							//jw�@�ɂ����ẮA�O�̃X�e�b�v�̒l�Ɋ�Â��������݂��Ȃ����߁A�����͕K�v�Ȃ�

							//double Ax=old_A[A_X][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Ay=old_A[A_Y][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//double Az=old_A[A_Z][N[m]];//�ߓ_N[m]�̂P�ï�� �O���޸�����ݼ��X����
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
							int n=dn[jside];

							//����zdV=V*Zs�����A����z*zdV��V*Zs*Zs�ł͂Ȃ��B(x,y�����l) ���ȏ�p84
									//(Nk)x�E(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y�E(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z�E(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"�f�B���N���l�����m���̂͂��̍s�ԍ��ɑ���"<<endl;
								//B_N=complex<double>(0,PHAT[n]*B_temp*2);//co�ɂ�j���������Ă���̂ŋ������֒ǉ�//�N����Ȃ��H
								B_N=complex<double>(0,B_temp*2);
							}
							else B_N=complex<double>(0,B_temp);
							//else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//���̋L�q�͑��v�H
							B[I-1]-=PHAT_A[n]*B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-��
				double co2=sig[je]*delta*delta6*delta6*delta6;//��V*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
					if(EDGE[iside].boundary_condition==0)///���m��
					{   
						int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
						int I1=EDGE[iside].node[1];//iside���\������2�_
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//���̂łȂ��Ƃ������Ƃ̓�=0�Ȃ̂ŁA������A�̍��̂悤�ȃf�B���N���l�Ɋւ��čl����K�v�͂Ȃ�
						}
					}
				}
				//cout<<"A-�� eddy"<<endl;

				//��-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"��-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//�v�fje�̕Ӕԍ�
							int J=npp[iside]+1;///��iside�͍s���I�Ԗ�
							//cout<<"��-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//iside���\������2�_
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///���m��
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]�̌v�Z
							}
							else //jside���ިظڌ^���E�ߓ_�Ȃ�
							{
								int n=dn[iside];
								B_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHAT�͕��f�f�B���N���ɔ����Ȃ����ׂ�
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"��-A eddy"<<endl;

				//��-��
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//�����1/j���������Ă��邽�߁A���f���Ƃ��Ă͕�����ς��ċ������ɑ�����邱�ƂɂȂ�(1/j=-j)
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
									{
										G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//��������i�����݂��Ȃ���������
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]���ިظڌ^���E�ߓ_�Ȃ�
							{
								int NN=dn[N[j]];
								B_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //�ӂ��̋��E�����ł̓ӂ��l�������Ȃ��̂ŁA�ЂƂ܂��폜
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}//////*/
					}
				}
				//cout<<"��-�� eddy"<<endl;
			}
		}
	}
	}

	///////
	cout<<"G,ROW�쐬����"<<endl;

	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"��[���v�f��"<<number<<endl;

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix_complex(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	cout<<"�s��̑Ώ̐��`�F�b�N----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	cout<<"����"<<endl;;
	//���s��̒l�`�F�b�N
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    /////////////////////

	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG�@���s
	complex<double> *XX=new complex<double> [pn];//�s��̓����i�[
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XX�Ɋi�[���ꂽ�����e�ӂɐU��
	cout<<"���f�����̊���U�肨��яo��----";
	
	/*
	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	*/
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
			if(EDGE[i].stat==ON) cout<<"�ÓI�ӗv�f�����m�������H"<<endl;
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	delete [] XX;


	//�i�[�����l����g���A�ʑ��������߂ăx�N�g���|�e���V���������߂�
	///Am,�ӂ�̧�ُo�� �ӂ͓d�ʂł͂Ȃ��ʑ��̒x��
	//ofstream a("Am_e.dat");
	//ofstream p("phi_e.dat");
	//double Am=0.0;
	double Vm=0.0;
	//double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		//if(EDGE[i].boundary_condition==0)
		{
			Am[i]=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi[i]=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am[i]*cos(omega*TIME+phi[i]);
		}
		//a<<Am<<endl;
		//p<<phi<<endl;
	}
	//cout<<"ok"<<endl;
	//a.close();
	//p.close();

	
	/////�Q�d������ 
	double *Je[3];//�Q�d��
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//�Q�d����[W]
	double *Je_loss_n=new double[node+1];//�Q�d����[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//�Q�d������v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//�d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//�{���̑̐�

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//�ߓ_�̉Q�d������Ή����闱�q�֓n��
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//���E���������Ƃɖ߂�(���̽ï�߂ōĂѓd�����v�Z���邩������Ȃ�����)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
	delete [] PHAT;
    delete [] PHAT_A;
	delete [] PHAT_V;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//�d�ʁA���ʂȂǂ̃|�e���V�����v�Z�֐�
void potential_calculation(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei)
{
	double *V=new double [node+1];	//potential
    
    double *Ee[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Ee[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓������܂��͗U�d���i�[�i���݂ł͗U�d���͎g���Ă��Ȃ��j

	double fluid_rp=CON->get_r_perm();		//�U�d��
	if(CON->get_EM_calc_type()==4) fluid_rp=CON->get_RP();	//������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=fluid_rp;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL)
		{
			if(CON->get_EM_calc_type()==1) RP[i]=1000;//�{���͓��̂����灇
			else if(CON->get_EM_calc_type()==4) RP[i]=1;//������
		}
		//else cout<<"error:�ގ��ɑ΂��A�U�d�����邢�͓��������s��"<<endl;
	}

	//�d�ʉ���
	VOLT3D(CON,NODE,ELEM, node, nelm,V,jnb, TIME,PART, fluid_number,nei,RP);


	delete [] V;
    for(int D=0;D<3;D++) delete [] Ee[D];
	delete [] RP;
}

///�d�ʌv�Z�֐�
void VOLT3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei,double *RP)
{
    double V1=0;
    double V2=CON->get_V();			//�d��
	double u0=0.000001256;			//�^��̓�����
	double ep0=8.854e-12;			//�^��̗U�d��
	double le=CON->get_distancebp();//���q�ԋ���
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

	////////�d�ʌ���
	if(CON->get_EM_calc_type()==1)//�d��v�Z�Ȃ�
	{
		double tau=CON->get_dt()*CON->get_V_step();//���萔
		if(CON->get_V_con()==2)
		{
			V2=CON->get_V()*(1-exp(-TIME/tau));//�d�ʂ��w���֐��I�ɑ���������
		}
		else if(CON->get_V_con()==1)
		{
			V2=(CON->get_V()-CON->get_initial_V())*TIME/(CON->get_dt()*CON->get_V_step())+CON->get_initial_V();//�d�ʂ𒼐��I�ɑ���������
		}
		if(V2==0) V2=1;//0���ƃG���[�ɂȂ邩��1�ɂ���
		else if(V2>CON->get_V()) V2=CON->get_V();
		cout<<"�d�ʌv�Z�J�n V="<<V2<<" ";
	}
	else if(CON->get_EM_calc_type()==4)//���ʌv�Z
	{
		double R=1;//�䗦
		V1=0;
		if(CON->get_uniform_B_sw()==OFF) V2=CON->get_magnet_B()/u0*CON->get_magnet_H();
		if(CON->get_uniform_B_sw()==ON)  V2=CON->get_uniform_B()/u0*(CON->get_ZU()-CON->get_ZD());//��l����
		if(CON->get_V_con()==1) 
		{
			V2=(V2)*TIME/(CON->get_dt()*CON->get_V_step());
			R=TIME/(CON->get_dt()*CON->get_V_step());
			if(V2>CON->get_magnet_B()/u0*CON->get_magnet_H()) {V2=CON->get_magnet_B()/u0*CON->get_magnet_H();R=1;}
		}
		if(V2==0) V2=1;//0���ƃG���[�ɂȂ邩��1�ɂ���
		
		cout<<"���ʌv�Z�J�n V2="<<V2<<"("<<R*CON->get_magnet_B()<<"T)"<<endl;
	}
	////////////////////////////

	///���ʌv�Z�̏ꍇ�A��͗̈�̋��E���������R���E�����ɂ���B�i�����ŁAMPSTOFEM�̒i�K�ł��̂悤�ɂ�����FINE3D�����܂������Ȃ��j
	if(CON->get_EM_calc_type()==4)//���ʌv�Z
	{
		if(CON->get_uniform_B_sw()==OFF)//�ʏ�͂�����
		{
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					for(int k=1;k<=4;k++)
					{
						int kelm=ELEM[i].elm[k];
						if(kelm==0)
						{
							int j=k%4+1;//ielm��jelm�̐ڂ���O�p�`���\������ߓ_�ԍ� 
							int m=4-(k-1)/2*2;
							int n=3-(k/2%2)*2;
							NODE[ELEM[i].node[j]].boundary_condition=0;//���E�����̖���
							NODE[ELEM[i].node[m]].boundary_condition=0;
							NODE[ELEM[i].node[n]].boundary_condition=0;
						}
					}
				}
			}
		}
		if(CON->get_uniform_B_sw()==ON)//��l����̏ꍇ
		{
			double ZU=CON->get_ZU();		//��͗̈��[
			double ZD=CON->get_ZD();		//��͗̈扺�[
			double err=1e-14;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].r[A_Z]>ZU-err) NODE[i].boundary_condition=2;//��[
				else if(NODE[i].r[A_Z]<ZD+err) NODE[i].boundary_condition=1;//���[
				else NODE[i].boundary_condition=0;
			}
		}
	}////////*/


    int NN=0;					//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1];	//�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
    double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
    
    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1)//��ŋ��E���������������Ă邩��A���̕���else if�ɂ��邱��
		{    
	        dn[i]=NN;
	        PHAT[NN]=V1;
	        V[i]=V1;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{   
	        dn[i]=NN;
	        PHAT[NN]=V2;
	        V[i]=V2;
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"�ިظڐ���"<<NN<<" ";
	    
    int pn=node-NN;///���m��
    int *ppn=new int [pn];		//�s���n�Ԗڂ̖��m���͐ߓ_�ԍ�ppn[n]�ɑ���
    int *npp=new int [node+1];	//i�Ԗڂ̐ߓ_�͍s���npp[i]�Ԗڂɑ����B�ިظڌ^�̏ꍇ�͍s��ɓ����Ă��Ȃ��̂�pn+1���i�[
    int num=0; 
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)	//���m���Ȃ�
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//������������������Ƃ߂�B�����ǎ኱�̌v�Z�R�X�g�ɂȂ�B
    //////
    
    ////�z��m��
    double **G=new double *[pn+1];						///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1];							///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1];							///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
	
	/*//���s��a�ɐߓ_�̓d�ׂ���
	if(CON->get_charge()==1 && CON->get_FEM_calc_type()==1)
	{
		for(int k=1;k<=WATER_N;k++)//���߂�WATER_N�̐ߓ_��TRANS[k]�̗��q�ɑ���
		{
			if(NODE[k].boundary_condition==0)//���m���Ȃ�B�ɂ��̐ߓ_�̗̈悪���邩����
			{
				int i=TRANS[k];//�ߓ_k�ɑ������闱�q�ԍ�
				for(int n=1;n<=jnb[k];n++)
				{
					int jelm=nei[k][n];//�ߓ_k�̗אڂ���v�f
					int N[5];
					for(int j=1;j<=4;j++) N[j]=ELEM[jelm].node[j];
		
					///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
					ELEM[jelm].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
					B[k]+=PART[i].angle/ep0/4*ELEM[jelm].volume;
					//if(PART[i].angle!=0)cout<<"EE"<<endl;
				}
				
			}
		}
	}///////////*/
	
	double *charge=new double [nelm+1];
	for(int n=1;n<=nelm;n++) charge[n]=0;
	
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
		double delta6=ELEM[je].volume;//�̐ς�6�{
		
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		/////��U�d�����`
		double ep=RP[je];
		/*if(ELEM[je].material==FLUID)
		{
			if(CON->get_EM_calc_type()==1) ep=CON->get_r_perm();
			else if(CON->get_EM_calc_type()==4) ep=CON->get_RP();//���ʌv�Z������䓧����
		}*/
	
		////�v�f��ظ��쐬�J�n
		for(int i=1;i<=4;i++)
		{
			if(NODE[N[i]].boundary_condition==0)///���m�Ȃ�
			{   
				int I=npp[N[i]]+1;///�ߓ_N[i]�͍s���I�Ԗ�
				for(int j=1;j<=4;j++)
				{					
					int J=npp[N[j]]+1;///�ߓ_N[j]�͍s���J�Ԗ�
					if(NODE[N[j]].boundary_condition==0)///���m�Ȃ�
					{   
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
							{
								G[I][h]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								//if(I==85839 && h==1) cout<<h<<" "<<G[I][h]<<" "<<delta<<" "<<(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])<<endl;
								flag=1;
							}
						}
						if(flag==0)
						{   
							NUM[I]=NUM[I]+1;
							int H=NUM[I];
			    
							G[I][H]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
							ROW[I][H]=J;
						}
					}
					else //N[j]���ިظڌ^���E�ߓ_�Ȃ�
					{
						int n=dn[N[j]];
						B[I-1]-=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
					}
				}
				//if(CON->get_charge()==ON) B[I-1]+=charge[je]/4*delta;
			}
		}   
    }
    ///////////////////////*/
    
	///���ۂ̍ő啝���v�Z
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	arrange_matrix(pn,NUM,ROW,G);

	//�Ώ̐��`�F�b�N
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	cout<<"�s��쐬�I�� ��:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    
	///////////////////////�s��v�Z�J�n
	
	double *XX=new double [pn];//�s��̓����i�[
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		V[i]=XX[n];//�d�ʃ�
	}
	delete [] XX;
	////////////////////////////
    
    ///////�d�ʂ�̧�ُo��
	ofstream fp("V.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<2*le && NODE[i].r[A_Y]>-2*le) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<V[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] charge;
}

//�s��̕��v�Z�֐�(�V)
int calc_matrix_width(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
	//�W���s��̂Ȃ��ŁAi�s�ڂɊ܂܂��񂺂�v�f���́A�v�Z�_i����̂т�ӂ̐�(���Ȃ킿�ߗ׌v�Z�_��)+1(���g)�ł���B����̍ő�l�����߂�
	///�z��m��
		
	int *temp_check=new int[node+1];	//�ꎞ�I�������z��
		
	///������
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON->get_FEM_calc_type()==1 || CON->get_FEM_calc_type()==4 || CON->get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].boundary_condition==0)
			{
				int width=0;
				vector<int> NEI2;//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�ANEI2�͐ߓ_-�ߓ_)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width++;
							NEI2.push_back(p);
							temp_check[p]=1;//��������
						}
					}
				}
				for(int k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//������
				if(width>maxwidth) maxwidth=width;
			}
		}
	}

	delete [] temp_check;

	return maxwidth+1;//width�͎�������L�т�ӂ̐��ɓ������B�s��̕��͎������g�������width+1
}

int calc_matrix_width2(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,int *width)
{
	//�W���s��̂Ȃ��ŁAi�s�ڂɊ܂܂��񂺂�v�f���́A�v�Z�_i����̂т�ӂ̐�(���Ȃ킿�ߗ׌v�Z�_��)+1(���g)�ł���B����̍ő�l�����߂�
	///�z��m��
		
	int *temp_check=new int[node+1];	//�ꎞ�I�������z��
		
	///������
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON->get_FEM_calc_type()==1 || CON->get_FEM_calc_type()==4 || CON->get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			width[i]=0;//������
			
			if(NODE[i].boundary_condition==0)
			{
				int width_temp=0;
				vector<int> NEI2;//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�ANEI2�͐ߓ_-�ߓ_)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width_temp++;
							NEI2.push_back(p);
							temp_check[p]=1;//��������
						}
					}
				}
				for(int k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//������
				if(width_temp>maxwidth) maxwidth=width_temp;
				
				width[i]=width_temp+1;
			}
			
			
		}
	}

	delete [] temp_check;

	return maxwidth+1;//width�͎�������L�т�ӂ̐��ɓ������B�s��̕��͎������g�������width+1
}

//�s����ёւ��֐�
void arrange_matrix(int pn,int *NUM,int **ROW,double **G)
{
	///G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}
//���ёւ��A���f��
void arrange_matrix_complex(int pn,int *NUM,int **ROW,complex<double>**G)
{
	///G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
	
    for(int i=1;i<=pn;i++)
    {
        complex<double> tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}

///ICCG�@
void ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	int count;						//�����グ�ϐ�
	double accel=CON->get_CGaccl();	//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *matrix=new double[num2];//�W���s���ۑ�(�Ίp�������܂މ��O�p�s��) ��[���v�f������1���z��Ƃ��ĕۑ�
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				matrix[num2]=val[m];
				ind2[num2]=ind[m];
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];			//������ɂ݂��A�e��̗v�f��
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//�[���v�f�̒l VAL[i][k]��i���k�Ԗڂ̔�[���v�f
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=1;//�덷
    double *r=new double[pn];
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *L=new double[num2];//�s���S����X�L�[������̉��O�p�s��L�i�[
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�����
	Incomplete_Cholesky_Decomposition(CON,L,D1,matrix,ptr2,ind2,pn,num2);//L��D1�ɒl���������܂��
	
	delete [] matrix;

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=L[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/L[0]; //���i3.77�j 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=L[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
			y[i]=(r[i]-sum)/L[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"ICCG�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	count=0;
	double ep=CON->get_FEMCGep();//��������
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alp�����߂�
		rLDLt_r=0;

		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
		}
		double PAP=0;
		for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) X[n]+=alp*P[n];
		
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		//////////////////�덷
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/L[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
					sum+=L[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
				y[i]=(r[i]-sum)/L[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;

	delete [] L;
}

//COCG�@�@���f���p
void COCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
	}
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	for(int n=0;n<pn;n++) BB= BB+norm(B[n]);
	BB=sqrt(BB);
	//cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	
	ofstream c("convergence.dat");

	 cout<<"COCG�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;
			//cout<<count<<"���"<<endl;
			complex<double> rr_be;
			rr_be=rr;
			//cout<<"rr="<<rr<<endl;
			//////////////alp�����߂�
			for(int n=0;n<pn;n++)
			{    
				AP[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
			}
			//cout<<"AP[0]="<<AP[0]<<" AP[1]="<<AP[1]<<endl;//AP[0]=4+i
			complex<double> PAP=(0.0,0.0);
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
			//for(int n=0;n<pn;n++)  PAP+=conj(P[n])*conj(AP[n]);
			alp=rr/PAP;
			//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
			//////////////////////
		
			//////////////// ���X�V�@X(k+1)=X(k)+alp*P
			for(int n=0;n<pn;n++) X[n]+=alp*P[n];
			//cout<<"X[0]="<<X[0]<<" X[1]="<<X[1]<<endl;
			//////////////////////////////
			
			//////////////// r=r-alp*AP
			for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
			/////////////////////////////
			
			///////////////////////beta
			rr=complex<double> (0.0,0.0);
			//#pragma omp parallel for
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
			//cout<<"rrk+1="<<rr<<endl;

			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			beta=rr/rr_be;
			//cout<<"beta="<<beta<<endl;
			///////////////////////

			//////////////////�덷
			//#pragma omp parallel for reduction(+:E)
			E=0.0;
			for(int n=0;n<pn;n++) E+=norm(r[n]);
			E=sqrt(E);
			E/=BB;
			c<<count<<" "<<E<<endl;
			//cout<<"E="<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
			//cout<<"P[0]="<<P[0]<<" P[1]="<<P[1]<<endl;
		}
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] AP;
	delete [] P;
	c.close();
}

//�O������COCG
void ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	//int count=0;
	complex<double> rr;
	rr=complex<double> (0,0);
	double E=1;//�덷
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];

	//�����ߒ���������
	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0,0);
		r[n]=B[n];
		P[n]=r[n];
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	BB=sqrt(BB);
	cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];

	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//�����t�@�N�^
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//���f�V�t�g�p�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				//if(ind[m]==k) val2[num2]+=accel2*Im;//���f�V�t�g
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	complex<double> rLDLt_r;
	complex<double> *y=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///�s���S�ڽ����������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"ICCOCG�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	//cout<<"ICCG�@:���m��="<<pn<<" ---";
	//unsigned int time=GetTickCount();
	int count=0;
	rLDLt_r=complex<double> (0,0);
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//�ŏ���rLDLt_r���������ŋ��߂�
	//while(E>CON->get_FEMCGep())
	double ep=CON->get_FEMCGep();
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		if(count%10000==0) ep*=10;
		//////////////alp�����߂�
		complex<double> PAP;
		PAP=complex<double> (0,0);
		//#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=complex<double>(0,0);
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
		//#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
			r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
			E+=norm(r[n]);						//�X�V��̌덷
		}
		E=sqrt(E);
		E/=BB;
		//cout<<"E="<<E<<endl;
		c<<count<<" "<<E<<endl;
		if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
		if(E<CON->get_FEMCGep()) cout<<"E<�� E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		
		beta=one/rLDLt_r;
		rLDLt_r=complex<double>(0,0);
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j �V
			else
			{
			   sum= complex<double> (0,0);
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
		for(int i=pn-1;i>=0;i--)
		{
		   sum= complex<double> (0,0);
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		//for(int n=0;n<pn;n++) rLDLt_r+=conj(r[n])*conj(LDLt_r[n]);
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	
	delete [] r;
	delete [] AP;
	delete [] P;
}

//openMP�ɂ���ĕ��񉻂���iccocg
void parallel_ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	//int count=0;
	complex<double> rr;
	rr=complex<double> (0,0);
	double E=1;//�덷
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];
	complex<double> *sum3 = new complex<double> [pn];

	//�����ߒ���������
	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0,0);
		r[n]=B[n];
		P[n]=r[n];
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	BB=sqrt(BB);
	cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];

	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//�����t�@�N�^
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//���f�V�t�g�p�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				//if(ind[m]==k) val2[num2]+=accel2*Im;//���f�V�t�g
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	complex<double> rLDLt_r;
	complex<double> *y=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///�s���S�ڽ����������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"parallel_ICCOCG�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	//cout<<"ICCG�@:���m��="<<pn<<" ---";
	//unsigned int time=GetTickCount();
	int count=0;
	rLDLt_r=complex<double> (0,0);
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//�ŏ���rLDLt_r���������ŋ��߂�
	//while(E>CON->get_FEMCGep())
	double ep=CON->get_FEMCGep();

	//#ifdef _OPENMP
	//#pragma omp parallel
	//#endif
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		if(count%10000==0) ep*=10;
		//////////////alp�����߂�
		complex<double> PAP=0.;

		#pragma omp parallel shared(PAP)
		{
			complex< double > priv_PAP=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				AP[n]=complex<double>(0,0);
				for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
				priv_PAP+=P[n]*AP[n];
			}
			#pragma omp critical
			{
			  PAP += priv_PAP;
			}
		}
		#pragma omp parallel// private(sum) 
		{
			//#pragma omp for reduction(+:PAP)
			//for(int n=0;n<pn;n++) PAP+=P[n]*AP[n];
			
			#pragma omp single
			{
				alp=rLDLt_r/PAP;
				E=0;
			}
			
			#pragma omp for reduction(+:E)
			for(int n=0;n<pn;n++)
			{
				X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
				r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
				//E+=norm(r[n]);						//�X�V��̌덷
				E+=abs(r[n])*abs(r[n]);
			}
		//}
			//#pragma omp for reduction(+:E)
			//for(int n=0;n<pn;n++) E+=norm(r[n]);						//�X�V��̌덷
			#pragma omp single
			{
				E=sqrt(E);
				E/=BB;	
				//cout<<"E="<<E<<endl;
				c<<count<<" "<<E<<endl;
				if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
				if(E<CON->get_FEMCGep()) cout<<"E<�� E="<<E<<endl;
				///////////////////////beta
				beta=one/rLDLt_r;
				rLDLt_r=complex<double>(0,0);

				y[0]=r[0]/val2[0];
			}
		}	
			/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
			//#pragma omp parallel for 
			for(int i=1;i<pn;i++)
			{
				sum3[i]=complex<double>(0,0);
			    //complex<double> priv_sum1=0.; 
				//#pragma omp parallel for 
				for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
				{
					sum3[i]+=val2[m]*y[ind2[m]];//���i3.78�j
				}
				y[i]=(r[i]-sum3[i])/val2[ptr2[i+1]-1];
				//int m2=ptr2[i+1]-1;
				//y[i]=(r[i]-sum3[i])/val2[ptr2[i+1]-1];
			}////y[i]�����Ƃ܂����B
			
			//#pragma omp parallel for 

		//}
			/////////LDLt_r[i]�����߂�
		//	#pragma omp for 
			for(int i=pn-1;i>=0;i--)
			{
				sum= complex<double> (0,0);
				//#pragma omp parallel for   
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y[i]-D1[i]*sum;	
			}
			/////////////////*/
		//}
		
			//for(int n=0;n<pn;n++) rLDLt_r+=conj(r[n])*conj(LDLt_r[n]);
		#pragma omp parallel shared(rLDLt_r)
		{
			complex< double > priv_rLDLt_r=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				priv_rLDLt_r+=r[n]*LDLt_r[n];
			}
			#pragma omp critical
			{
			  rLDLt_r += priv_rLDLt_r;
			}
		}

		beta=beta*rLDLt_r;
			/////////////////*/
			
			///////////////////// P=r+beta*P
		#pragma omp parallel for
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}

	

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	
	delete [] r;
	delete [] AP;
	delete [] P;
	delete [] sum3;
}

//////mps115-1-test1���炻�̂܂܈ڐA��������iccg
void parallel_ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	//�����ߒ���������
	ofstream c("convergence.dat");

	double accel=CON->get_CGaccl();//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	double BB=0;
	for(int k=0;k<pn;k++)
	{
		BB+=B[k]*B[k];
	}
	BB=sqrt(BB);

	cout<<"BB="<<BB<<endl;

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=BB;//�덷
    double *r=new double[pn];
	for(int n=0;n<pn;n++) X[n]=0;
	
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
	        if(i==0)
			{
				val2[m]=val2[m];
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				double sum=0;
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				D1[k]=1/val2[m];
            }
	    }
	}    
	///�s���S�ڽ����������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"parallel_ICCG�@�X�^�[�g  -----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep;
	ep=CON->get_FEMCGep();//��������
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count>pn E="<<E<<endl;
		if(count==10000) cout<<"count=10000 E="<<E<<endl;
		if(count%10000==0) ep*=10;
		//////////////alp�����߂�
		rLDLt_r=0;
		double PAP=0;
		//for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//sirial
		#pragma omp parallel for reduction(+:rLDLt_r) reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{	
			//printf("%d\n",omp_get_thread_num());
			rLDLt_r+=r[n]*LDLt_r[n];
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];//sirial
		alp=rLDLt_r/PAP;
		
//		cout<<"alp="<<alp<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++) 
		{
			X[n]+=alp*P[n];	// X(k+1)=X(k)+alp*P
			r[n]-=alp*AP[n];// r=r-alp*AP
			E+=r[n]*r[n];	//�덷
		}
		E=sqrt(E)/BB;
		if(count%10==0) c<<count<<" "<<E<<endl;
		//c<<count<<" "<<E<<endl;
		if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
		//////////////////////////////
		
		//////////////// r=r-alp*AP
		//for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		/*/////////////////�덷
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		///////////////////////*/
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
}

//MRTR�@�@���f���p// "MRTR�@�̕��f�Ώ̐��^�������ւ̊g��(2007)"
void cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(���f�������o�[�ŕ\�L)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	/*/////1�X�e�b�v�ڂ̌v�Z///////////////

	for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
	{    
		Ar[n]=complex<double>(0.0,0.0);
		cAr[n]=complex<double>(0.0,0.0);
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			Ar[n]+=val[j]*r[ind[j]];
			cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
		}
	}

	cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
	//#pragma omp parallel for
	for(int n=0;n<pn;n++)  
	{
		cAr_r+=conj(cAr[n])*r[n];
		cAr_Ar+=conj(cAr[n])*Ar[n];
	}

	zeta=cAr_r/cAr_Ar;
	eta=complex<double>(0.0,0.0);

	//////v�̌v�Z//////////////
	for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

	//////p�̌v�Z//////////////
	for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
	zeta_old=zeta;

	///////X�̌v�Z////////////
	for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
	//////y�̌v�Z////////////
	for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

	//////r�̌v�Z////////////
	for(int n=0;n<pn;n++) r[n]-=y[n];
	///////*/


	
	//for(int n=0;n<pn;n++) rr+=norm(r[n]);
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	//E=sqrt(rr/BB);

	//cout<<"������="<<1<<" E1="<<E<<endl;
	//count++;
	

	 cout<<"cs_MRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
			{    
				Ar[n]=complex<double>(0.0,0.0);
				cAr[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
					cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
				}
			}

			cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
			cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				cAr_r+=conj(cAr[n])*r[n];
				cAr_Ar+=conj(cAr[n])*Ar[n];
				cy_Ar+=y[n]*Ar[n];
				cAr_y+=conj(cAr[n])*y[n];
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////v�̌v�Z//////////////
			v=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//�덷�]��
			rr=0;
			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//�Ίp�X�P�[�����O��MRTR�@�@���f���p// "MRTR�@�̕��f�Ώ̐��^�������ւ̊g��(2007)"
void DS_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(���f�������o�[�ŕ\�L)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���

	//�Ίp�X�P�[�����O
	for(int n=0;n<pn;n++)
	{    
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			val[j]=val[j]*r[ind[j]];
		}
	}


	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	/*/////1�X�e�b�v�ڂ̌v�Z///////////////

	for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
	{    
		Ar[n]=complex<double>(0.0,0.0);
		cAr[n]=complex<double>(0.0,0.0);
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			Ar[n]+=val[j]*r[ind[j]];
			cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
		}
	}

	cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
	//#pragma omp parallel for
	for(int n=0;n<pn;n++)  
	{
		cAr_r+=conj(cAr[n])*r[n];
		cAr_Ar+=conj(cAr[n])*Ar[n];
	}

	zeta=cAr_r/cAr_Ar;
	eta=complex<double>(0.0,0.0);

	//////v�̌v�Z//////////////
	for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

	//////p�̌v�Z//////////////
	for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
	zeta_old=zeta;

	///////X�̌v�Z////////////
	for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
	//////y�̌v�Z////////////
	for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

	//////r�̌v�Z////////////
	for(int n=0;n<pn;n++) r[n]-=y[n];
	///////*/


	
	//for(int n=0;n<pn;n++) rr+=norm(r[n]);
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	//E=sqrt(rr/BB);

	//cout<<"������="<<1<<" E1="<<E<<endl;
	//count++;
	

	 cout<<"cs_MRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
			{    
				Ar[n]=complex<double>(0.0,0.0);
				cAr[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
					cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
				}
			}

			cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
			cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				cAr_r+=conj(cAr[n])*r[n];
				cAr_Ar+=conj(cAr[n])*Ar[n];
				cy_Ar+=y[n]*Ar[n];
				cAr_y+=conj(cAr[n])*y[n];
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////v�̌v�Z//////////////
			v=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//�덷�]��
			rr=0;
			for(int n=0;n<pn;n++) rr+=norm(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//����MRTR�@�@���f���p// "MRTR�@�̕��f�Ώ̐��^�������ւ̊g��(2007)"
void parallel_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(���f�������o�[�ŕ\�L)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	
	/*////
	#pragma omp parallel shared(PAP)
		{
			complex< double > priv_PAP=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				AP[n]=complex<double>(0,0);
				for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
				priv_PAP+=P[n]*AP[n];
			}
			#pragma omp critical
			{
			  PAP += priv_PAP;
			}
		}
		#pragma omp parallel// private(sum) 
		{
			//#pragma omp for reduction(+:PAP)
			//for(int n=0;n<pn;n++) PAP+=P[n]*AP[n];
			
			#pragma omp single
			{
				alp=rLDLt_r/PAP;
				E=0;
			}
			
			#pragma omp for reduction(+:E)
			for(int n=0;n<pn;n++)
			{
				X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
				r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
				E+=norm(r[n]);						//�X�V��̌덷
			}
		//}
			//#pragma omp for reduction(+:E)
			//for(int n=0;n<pn;n++) E+=norm(r[n]);						//�X�V��̌덷
			#pragma omp single
	////////*/
	cout<<"parallel_cs_MRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			#pragma omp parallel 
			{
				#pragma omp for
				for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
				{    
					Ar[n]=complex<double>(0.0,0.0);
					cAr[n]=complex<double>(0.0,0.0);
					for(int j=ptr[n];j<ptr[n+1];j++)
					{
						Ar[n]+=val[j]*r[ind[j]];
						cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
					}
				}
			
				#pragma omp single
				{
					cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
					cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
					cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
					cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
				}
			}
			//#pragma omp parallel for
			#pragma omp parallel shared(cAr_r,cAr_Ar,cy_Ar,cAr_y)
			{
				#pragma omp for
				for(int n=0;n<pn;n++)  
				{
					cAr_r+=conj(cAr[n])*r[n];
					cAr_Ar+=conj(cAr[n])*Ar[n];
					cy_Ar+=y[n]*Ar[n];
					cAr_y+=conj(cAr[n])*y[n];
				}
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////v�̌v�Z//////////////
			v=complex<double>(0.0,0.0);

			#pragma omp parallel shared(v)
			{
				#pragma omp for
				for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];
			}

			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//�덷�]��
			rr=0;
			for(int n=0;n<pn;n++) rr+=norm(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//ic�t��MRTR�@
void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> nu;

	complex<double> cAr_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

	complex<double> cw_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cv_w;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cy_w;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	complex<double> cw_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(���f�������o�[�ŕ\�L)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];
	complex<double> *u = new complex<double> [pn];
	complex<double> *v = new complex<double> [pn];
	complex<double> *w = new complex<double> [pn];
	complex<double> *z = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	nu=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
		u[n]=complex<double> (0.0,0.0);
		v[n]=complex<double> (0.0,0.0);
		w[n]=complex<double> (0.0,0.0);
		z[n]=complex<double> (0.0,0.0);
	}

	//////�O����/////
	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//�����t�@�N�^
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//���f�V�t�g�p�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				//if(ind[m]==k) val2[num2]+=accel2*Im;//���f�V�t�g
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	complex<double> rLDLt_r;
	complex<double> *y2=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///�s���S�ڽ����������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y2[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y2[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y2[i]-D1[i]*sum;	
	}
	/////////////////*///LDLt��M-1�ɑ����H�܂�ALDLt_r��u0�Ɠ�����
	
	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
		/////////////////*/
	//////////////////////

	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;


	 cout<<"cs_ICMRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			/////v�̌v�Z
			for(int n=0;n<pn;n++)
			{    
				v[n]=complex<double>(0.0,0.0);//v�ɑ���
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					v[n]+=val[j]*u[ind[j]];
				}
			}

			/*///
			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
			{    
				Ar[n]=complex<double>(0.0,0.0);//v�ɑ���
				cAr[n]=complex<double>(0.0,0.0);//conj(v)�ɑ���
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*u[ind[j]];
					cAr[n]+=conj(val[j]) * conj(u[ind[j]]);
				}
			}
			///*/

			////w�̌v�Z

			/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
			for(int i=0;i<pn;i++)
			{
				if(i==0) y2[0]=v[0]/val2[0]; //���i3.77�j 
				else
				{
					sum=complex<double> (0,0);
					/////////        
					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//���i3.78�j
					int m=ptr2[i+1]-1;
					y2[i]=(v[i]-sum)/val2[m];
				}
			}////y[i]�����Ƃ܂����B
			for(int i=pn-1;i>=0;i--)
			{
				sum=complex<double> (0,0);
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y2[i]-D1[i]*sum;	
			}
			/////////////////*///LDLt��M-1�ɑ����H�܂�ALDLt_r��u0�Ɠ����� ���l�ɁALDL_v��w�H
	
			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];

			cw_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
			cv_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cy_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			cw_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
			//���ς̌v�Z
			for(int n=0;n<pn;n++)  
			{
				cw_r+=w[n]*r[n];
				cv_w+=v[n]*w[n];
				cy_w+=y[n]*w[n];
				cw_y+=w[n]*y[n];
			}

			zeta=nu*cw_r/(nu*cv_w-cy_w*cw_y);
			eta=-cy_w*cw_r/(nu*cv_w-cy_w*cw_y);

			//////nu�̌v�Z//////////////
			nu=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];

			//////p�̌v�Z//////////////
			for(int n=0;n<pn;n++) P[n]=u[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////X�̌v�Z////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////y�̌v�Z////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*v[n];

			//////r�̌v�Z////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//////z�̌v�Z////////////
			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];

			//////u�̌v�Z////////////
			for(int n=0;n<pn;n++) u[n]-=z[n];

			//�덷�]��
			rr=0;
			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	delete [] y2;

	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	delete [] u;
	delete [] v;
	delete [] w;
	delete [] z;

	c.close();
}

//�s���S�ڽ������
void Incomplete_Cholesky_Decomposition(mpsconfig *CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:�W���s��̂����A�Ίp���܂މ��O�p���i�[����Ă���
	//ptr2[k]:matrix[i]�ɑ΂��鱸���z��
	//L[i],D[i]:�s���S�ڽ�������i�[�z��
	//pn:���m��
	//num2:matrix[i]�̗v�f��

	//cout<<"�s���S�ڽ�������J�n---";
	unsigned int timeC=GetTickCount();

	int UP=1;							//
	int DOWN=2;
	int flag=DOWN;						//flag==UP�Ȃ�A�����W����1���傫�����Ă����BDOWN�Ȃ珬�������Ă���
	int loopsw;							//while���̃X�C�b�`

	double accel=1;//CON->get_CGaccl();		//�����W��
	double maxvalue=0;					//matrix���̍ő�l
	int Line=1;							//L�s��̂Ȃ��ōő�l�����s�ԍ�
	int Column=1;						//L�s��̂Ȃ��ōő�l������ԍ�
	
	
	//1���:�܂��͕��ʂɕs���S�ڽ���������s��()�B������D1[i]<0�Ȃ�����W���𑝉�����(UP)�BD1[i]>0�Ȃ�����W��������������(DOWN)
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
	if(D1[0]<0) cout<<"error in �s���S�ڽ������ D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
		{
			int i=ind2[m];//��ԍ�
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
			L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
		double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//�����W����������
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
		D1[k]=1/L[m];

		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
		
		if(L[m]<0)						//L[m]<0�̂Ƃ�D1[m]<0�ƂȂ�A�������ɒ[�ɒx���Ȃ�.�����h�����߂ɉ����W����傫������
		{
			flag=UP;
		}
	}
	

	if(flag==UP)		//�����W�����P���傫�����Ă����ꍇ
	{
		loopsw=ON;
		while(loopsw==ON)
		{
			accel+=0.05;
			L[0]=matrix[0]*accel;
			D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
			maxvalue=D1[0]*D1[0];
			for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
			{	
				for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
				{
					int i=ind2[m];//��ԍ�
		
					double sum=0;		//accel�̕ύX�ɂ��L��D1�̒l���ς�����̂ŁAsum���v�Z���Ȃ����B
					for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�

					if(L[m]*L[m]>maxvalue)
					{
						maxvalue=L[m]*L[m];
						Line=k; Column=i;
					}
				}
				int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
				double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
				L[m]=matrix[m]*accel-sum;		//�����W����������
				if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
				D1[k]=1/L[m];
				if(L[m]*L[m]>maxvalue)
				{
					maxvalue=L[m]*L[m];
					Line=k; Column=k;
				}
			}
			loopsw=OFF;
			for(int k=0;k<pn;k++) if(D1[k]<0) loopsw=ON;	//1�ł�����D������΂�����x���s����
		}

	}
	else if(flag==DOWN)		//�����W�����P��菬�������Ă����ꍇ
	{
		loopsw=ON;
		int flag2=OFF;
		while(loopsw==ON)
		{
			if(flag2==OFF) accel-=0.05;
			else if(flag2==ON)
			{
				accel+=0.05;
				loopsw=OFF;
			}
			
			L[0]=matrix[0]*accel;
			D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
			maxvalue=D1[0]*D1[0];
			for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
			{	
				for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
				{
					int i=ind2[m];//��ԍ�
		
					double sum=0;		//accel�̕ύX�ɂ��L��D1�̒l���ς�����̂ŁAsum���v�Z���Ȃ����B
					for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
					L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�

					if(L[m]*L[m]>maxvalue)
					{
						maxvalue=L[m]*L[m];
						Line=k; Column=i;
					}
				}
				int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
				double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
				L[m]=matrix[m]*accel-sum;		//�����W����������
				if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
				D1[k]=1/L[m];
				if(L[m]*L[m]>maxvalue)
				{
					maxvalue=L[m]*L[m];
					Line=k; Column=k;
				}
			}
			//loopsw=OFF;
			for(int k=0;k<pn;k++) if(D1[k]<0) flag2=ON;	//1�ł�����D������΁A���������������Ɣ��f���āA1�O�̉����W����p���čēx�v�Z����
			if(Line!=Column) flag2=ON;
		}
	}
	//if(Line!=Column) cout<<"value="<<sqrt(maxvalue)<<" Line="<<Line<<" Column="<<Column<<" D["<<Line<<"]="<<D1[Line]<<endl;

	/*/�Ίp�D�ʐ����`�F�b�N
	for(int k=0;k<pn;k++)
	{	
		double val=0;
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++) val+=fabs(L[m]);
		val*=2;
		if(val>L[ptr2[k+1]-1]) cout<<k<<" "<<val<<" "<<D1[k]<<endl;
	}//*/

	//cout<<"�����W��="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	cout<<"��="<<accel<<" ";
	
}

/*/�s���S�ڽ������(���j
void Incomplete_Cholesky_Decomposition(mpsconfig *CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:�W���s��̂����A�Ίp���܂މ��O�p���i�[����Ă���
	//ptr2[k]:matrix[i]�ɑ΂��鱸���z��
	//L[i],D[i]:�s���S�ڽ�������i�[�z��
	//pn:���m��
	//num2:matrix[i]�̗v�f��

	cout<<"�s���S�ڽ�������J�n---";
	unsigned int timeC=GetTickCount();
	double accel=CON->get_CGaccl();		//�����W��
	double accel2=0;					//�����v�Z���ꂽ�����W��
	double maxvalue=0;					//matrix���̍ő�l
	int Line=1;							//�ő�l�����s�ԍ�
	int Column=1;						//�ő�l������ԍ�
	
	//1���:�܂��͕��ʂɕs���S�ڽ���������s���B������D1[i]<0�Ȃ�����W�����C�����A2��ڂ��s��
	int onemoreflag=OFF;
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
	if(D1[0]<0) cout<<"error in �s���S�ڽ������ D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
		{
			int i=ind2[m];//��ԍ�
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
			L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
		double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//�����W����������
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
		D1[k]=1/L[m];
		
		if(L[m]<0)						//L[m]<0�̂Ƃ�D1[m]<0�ƂȂ�A�������ɒ[�ɒx���Ȃ�.�����h�����߂ɉ����W����傫������
		{
			accel2=sum/akk*1.1;
			if(accel2>accel) accel=accel2;
			onemoreflag=ON;				//���̃X�C�b�`��ON�Ȃ�A�s���S�ڽ��������������x�s���B
		}
		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
	}
	if(Line!=Column) cout<<"value="<<maxvalue<<" Line="<<Line<<" Column="<<Column<<endl;
	if(onemoreflag==ON)
	{
		L[0]=matrix[0]*accel;
		D1[0]=1/L[0];			//L[0]�͌W���s��̒l�ƕς��Ȃ��̂ŁA�[���ł��Ȃ���Ε��ł��Ȃ��Ɨ\�z�����
		for(int k=1;k<pn;k++)	//k=0�͂����ς񂾁Bk>=1����loop
		{	
			for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0����k-1�܂ł̗v�f
			{
				int i=ind2[m];//��ԍ�
		
				double sum=0;		//accel�̕ύX�ɂ��L��D1�̒l���ς�����̂ŁAsum���v�Z���Ȃ����B
				for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
				L[m]=matrix[m]-sum;//i=0�̂Ƃ���sum=0�ƂȂ莩���I��L[m]=matrix[m]�ƂȂ�
			}
			int m=ptr2[k+1]-1;				//����m�͕K���Ίp�����ɊY������
				
			double akk=matrix[m];			//�W���s��̑Ίp����.�l��ۑ�.
			double sum=0;
			for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
			L[m]=matrix[m]*accel-sum;		//�����W����������
			if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkk�����ɏ������Ƃ��A���̂悤�ɂ���B�����ŁA�ǂ݂̂�Lkk<0�Ȃ�������Â炢����ALkk�����̋ɏ��l�ł����̒l�ɂ��Ă��܂��Ă悢�̂ł́H
			D1[k]=1/L[m];
		}
	}	
	cout<<"�����W��="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	///�s���S�ڽ����������//////////
}
///*/

//�s��̑Ώ̐������֐�
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G)
{
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(abs((G[i][j]-G[J][k])/G[i][j])>1e-8)
					{
						//cout<<"�Ώ̐��װ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						cout<<"�Ώ̐��װ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0) cout<<"DDD i="<<i<<endl;
		}
	}
}

void check_matrix_symmetry_complex(int pn,int *NUM,int **ROW,complex<double> **G)
{
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						cout<<"�Ώ̐��װ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0) cout<<"DDD i="<<i<<endl;
		}
	}
	cout<<"�Ώ̐��`�F�b�N����"<<endl;
}

void diagonal_distribution(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int pn,int *NUM,int **ROW,double **G)
{
	ofstream di("dia.dat");
	ofstream xd("exdia.dat");

	double *dia=new double [pn];	//�Ίp����
	double *exdia=new double [pn];	//�Ίp�����ȊO�̘a

	double sum_di=0;
	double sum_xd=0;

	for(int i=1;i<=pn;i++)
	{
		double sum=0;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				dia[i]=G[i][j];
				di<<i<<" "<<dia[i]<<endl;
				sum_di+=dia[i];
			}
			else
			{
				sum+=abs(G[i][j]);
			}
		}
		exdia[i]=sum;
		xd<<i<<" "<<sum<<endl;
		sum_xd+=sum;
	}
	cout<<"�W���s��̃f�[�^�쐬"<<endl;
	di.close();
	xd.close();

	cout<<"r_sum="<<sum_xd<<endl;
	sum_di/=pn; sum_xd/=pn; //���ϒl

	//���U�̌v�Z
	double s2di=0;
	double s2xd=0;
	for(int i=1;i<=pn;i++)
	{
		s2di+=(sum_di-dia[i])*(sum_di-dia[i]);
		s2xd+=(sum_xd-exdia[i])*(sum_xd-exdia[i]);
	}

	s2di/=pn;
	s2xd/=pn;
	
	cout<<"s2D="<<s2di<<" s2r="<<s2xd<<" aveD="<<sum_di<<" aver="<<sum_xd<<endl;

	//�l���傫�ȂƂ��낾���𔲐����� //���ϒl��10�{�ȏゾ�����o
	ofstream di2("dia_p.dat");
	ofstream xd2("exdia_p.dat");

	for(int i=1;i<=pn;i++)
	{
		if(dia[i]>sum_di*5) di2<<i<<" "<<dia[i]<<endl;
		if(exdia[i]>sum_xd*5) xd2<<i<<" "<<exdia[i]<<endl;
	}

	di2.close();
	xd2.close();

	delete [] dia;
	delete [] exdia;

}

//��Ώ̍s���@
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :��[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��
	
	int count=0;
	double pk=1;
	double E=1;//�덷
	double alp,beta,rr,w,ita;

	double *r=new double [pn];	//�c��
	double *P=new double [pn];	//�T���x�N�g��
	double *Pj=new double [pn];
	double *rj=new double [pn];
	double *AP=new double [pn];
	double *AtPj=new double [pn];
	double *e=new double [pn];
	double *Ae=new double [pn];
	double *y=new double [pn];
	double *u=new double [pn];
	double *W=new double [pn];
	double *Z=new double [pn];
	double *e_pre=new double [pn];
	beta=0;
	w=0;
	ita=0;

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0;
		r[n]=B[n];
	}

	for(int n=0;n<pn;n++)
	{
		P[n]=r[n];//������
		rj[n]=r[n];
		Pj[n]=rj[n];
		AP[n]=0;
		y[n]=0;
		u[n]=0;
		W[n]=0;
		Z[n]=0;
		e[n]=0;
		e_pre[n]=0;//1�X�e�b�v�O��e[]
	}
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"rr0="<<rr0<<endl;
	 cout<<"BiCGstab2�@�X�^�[�g  -----���m��="<<pn<<"  ---";
	 unsigned int time=GetTickCount();
	double ep=CON->get_FEMCGep();//��������
	while(E>ep)// EP=CON->get_CGep();//��������(convergence test)
	{
		count++;

		for(int n=0;n<pn;n++)
		{
			//P[n]=r[n]+beta*(P[n]-w*AP[n]);
			P[n]=r[n]+beta*(P[n]-u[n]);
		} 

		////pk(r��rj�̓���)�����߂�
		pk=0;
		for(int n=0;n<pn;n++) pk+=r[n]*rj[n];

		//////////////alp�����߂�
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			//for(int m=0;m<pn;m++) AP[n]+=A[n][m]*P[m];
		}
		double APrj=0;
		for(int n=0;n<pn;n++)  APrj+=rj[n]*AP[n];
		alp=pk/APrj;
		//cout<<"alp="<<alp<<" APrj="<<APrj<<endl;
		//////////////////////

		for(int n=0;n<pn;n++) y[n]=e[n]-r[n]-alp*W[n]+alp*AP[n];

		for(int n=0;n<pn;n++)
		{
			e_pre[n]=e[n];//�ύX�O�̒l���L��
			e[n]=r[n]-alp*AP[n];
		}

		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m];
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]];
		}

		if(count%2!=0)
		{
			double e_dot_ae  = 0;
			double ae_dot_ae = 0;
			for(int n=0;n<pn;n++)
			{
				e_dot_ae+=e[n]*Ae[n];
				ae_dot_ae+=Ae[n]*Ae[n];
			}
			w = e_dot_ae / ae_dot_ae;
			ita=0;
		}
		else
		{
			double e_dot_ae  = 0;
			double y_dot_y=0;
			double y_dot_e=0;
			double Ae_dot_y=0;
			double Ae_dot_Ae=0;
			for(int n=0;n<pn;n++) e_dot_ae+=e[n]*Ae[n];
			for(int n=0;n<pn;n++) y_dot_y+=y[n]*y[n];
			for(int n=0;n<pn;n++) y_dot_e+=y[n]*e[n];
			for(int n=0;n<pn;n++) Ae_dot_y+=Ae[n]*y[n];
			for(int n=0;n<pn;n++) Ae_dot_Ae+=Ae[n]*Ae[n];

			double co=Ae_dot_Ae*y_dot_y-Ae_dot_y*Ae_dot_y;

			w=y_dot_y*e_dot_ae-y_dot_e*Ae_dot_y;
			w/=co;

			ita=Ae_dot_Ae*y_dot_e-Ae_dot_y*e_dot_ae;
			ita/=co;

		}
		

		for(int n=0;n<pn;n++) 
		{
			u[n]=w*AP[n]+ita*(e_pre[n]-r[n]+beta*u[n]);
			Z[n]=w*r[n]+ita*Z[n]-alp*u[n];
			X[n]+=alp*P[n]+Z[n];
			r[n]=e[n]-ita*y[n]-w*Ae[n];
		}

		///beta
		beta=0;
		for(int n=0;n<pn;n++) beta+=r[n]*rj[n];
		beta/=w*APrj;
		
		//cout<<"beta="<<beta<<endl;
		//////////////////////

		//W[n]
		for(int n=0;n<pn;n++) W[n]=Ae[n]+beta*AP[n];
		
		//////////////////�덷
		rr=0;
		for(int n=0;n<pn;n++) rr+=r[n]*r[n];
		E=rr/rr0;
		//E=sqrt(rr);
		//cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;

	delete [] r;
	delete [] Pj;
	delete [] rj;
	delete [] AP;
	delete [] AtPj;
	delete [] P;
	delete [] e;
	delete [] Ae;
	delete [] y;
	delete [] u;
	delete [] W;
	delete [] Z;
	delete [] e_pre;
}

//������v�Z�֐�
void calc_transitional_EM_field(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,int *jnb,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int **nei,int particle_node,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta,vector<edge3D> &static_EDGE)
{
	cout<<"�������͊J�n"<<endl;
	double *V=new double[node+1];
	for(int i=1;i<=node;i++) V[i]=0;				//�������@V�͉Q�d���̓d�ʊi�[�Ɏg��
	double *current[3];								//�e�v�f�̓d�����x[A/m2]�i�[
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];
	double *Be[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//�e�v�f�̓�����
	double *sigma=new double [nelm+1];	//�e�v�f�̓��d��

	//�䓧��������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:�ގ��ɑ΂��A���������s��"<<endl;
	}

	//���d������
	for(int i=1;i<=nelm;i++)
	{
		if(CON->get_temperature_depend()==OFF)
		{
			if(ELEM[i].material==FLUID) sigma[i]=CON->get_ele_conduc();
			else if(ELEM[i].material==AIR) sigma[i]=0;//�g�����Ƃ͂Ȃ��͂������A�ꉞ������
			else if(ELEM[i].material==COIL) sigma[i]=CON->get_ele_conduc2();
			else if(ELEM[i].material==CRUCIBLE) sigma[i]=CON->get_ele_conduc2();
			else cout<<"error:�ގ��ɑ΂��A���d�����s��"<<endl;
		}
	}

	if(CON->get_temperature_depend()==ON)
	{
		cout<<"���x�ˑ��̓��d���v�Z";
		double *sig_f=new double[fluid_number];
		double *sig_n=new double[node+1];
		calc_physical_property(CON,PART,fluid_number,sig_f,fluid_number,4);//��R���̉��x�ˑ�
		for(int i=0;i<fluid_number;i++) sig_f[i]=1/sig_f[i];
		for(int in=1;in<=node;in++) 
		{
			if(NODE[in].particleID>=0)
			{
				sig_n[in]=sig_f[NODE[in].particleID];
				if(NODE[in].material!=FLUID) cout<<"���̗��q���Ή����Ă��Ȃ��_�ŉ��x��v�Z"<<endl;
			}
		}

		//���d����v�f����ߓ_�ɕ��z����
		for(int je=1;je<=nelm;je++)
		{
			double dis=0;
			double dis_sum=0;
			int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
			double X[4+1];
			double Y[4+1];
			double Z[4+1];
			

			if(ELEM[je].material==FLUID)
			{	
				for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
	    
				double Xs=0;//�d�S���W
				double Ys=0;
				double Zs=0;
				for(int j=1;j<=4;j++)
				{
					X[j]=NODE[N[j]].r[A_X];
					Y[j]=NODE[N[j]].r[A_Y];
					Z[j]=NODE[N[j]].r[A_Z];
					Xs+=X[j]*0.25;
					Ys+=Y[j]*0.25;
					Zs+=Z[j]*0.25;
				}
		
				double delta=ELEM[je].volume/6;//�{���̑̐�

				for(int j=1;j<=4;j++)
				{
					dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
					dis_sum+=dis;
				}
				for(int j=1;j<=4;j++)
				{
					dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
					sigma[je]+=sig_n[N[j]]*dis/dis_sum;
				}
			}
			if(ELEM[je].material==AIR) sigma[je]=0;//�g�����Ƃ͂Ȃ��͂������A�ꉞ������
			if(ELEM[je].material==COIL) sigma[je]=CON->get_ele_conduc2();
			if(ELEM[je].material==CRUCIBLE) sigma[je]=CON->get_ele_conduc2();

			//�`�F�b�N
			if(ELEM[je].material==FLUID)
			{
				if(sigma[je]>1e8) cout<<"���d���傫�����H"<<endl;
				if(sigma[je]<1e6) cout<<"���d�����������H"<<endl;
			}

		}

		
		delete [] sig_f;
		delete [] sig_n;
		cout<<" ok"<<endl;
	}

	

	cout<<"�d���C��"<<endl;
	double f=CON->get_Hz();	//�𗬂̎��g��[Hz]
	double omega=2*PI*f;//�p���g��[rad/sec]
	if(CON->get_J_input_way()==1)					//�d�����x�𑼂̃\�t�g����ǂݍ���
	{
		inport_J0_density(CON, node, nelm,NODE,ELEM,current);		
		
		/*///
		////�ϊ�
		double I0=CON->get_I0();			//�d���̐U��
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				current[D][i]*=I0/900;//900��magnet���Ŏw�肵�Ă��鋭���d���B
			}
		}
		///*/

		/*///
		//�𗬕ϊ�
		//double I0=CON->get_I0();			//�d���̐U��
		//double II=I0*cos(omega*TIME);	
		//double II=I0*sin(2*PI*f*(TIME+CON->get_dt()));	
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				//current[D][i]*=II/I0;
				current[D][i]*=cos(omega*TIME);
			}
		}
		////*/

		////
		////�傫���ϊ�
		//double I0=CON->get_I0();
		double I0=CON->get_I0()*sqrt(2.0);	//�d���̐U�� //�����ł͎����l�𑪒肵�Ă���̂ŁA��2��������
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				current[D][i]*=I0/900;//900��magnet���Ŏw�肵�Ă��鋭���d���̑傫���B
			}
		}
		///*/

		check_J0(CON, node, nelm,NODE,ELEM,current);

		//�����d���̔��U���`�F�b�N ��(gradNi�EJ0)dv=0
		double div=0;
		int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		//double b[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		double j0x, j0y, j0z;
		for(int je=1;je<=nelm;je++)
		{   
			if(ELEM[je].material==COIL)
			{
				for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		    
				for(int j=1;j<=4;j++)
				{
					X[j]=NODE[N[j]].r[A_X];
					Y[j]=NODE[N[j]].r[A_Y];
					Z[j]=NODE[N[j]].r[A_Z];
				}
				
				///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
				ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
			
				double delta6=ELEM[je].volume;//�̐ς�6�{
			
				delta6=1/delta6/6;//�v�Z�ɕK�v�Ȃ̂�1/(36V)�Ȃ̂ł����Ŕ��]���Ă���(1/36V^2)*V=1/36V
			
				double delta=ELEM[je].volume/6;//�{���̑̐�
			
				for(int i=1;i<=4;i++)
				{
					int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
					int m=j%4+1;
					int n=m%4+1;
			    
					c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
					d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
					e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
					if(i%2!=0)//i����Ȃ�
					{
						c[i]*=-1;
						d[i]*=-1;
						e[i]*=-1;
					}
				}
				/////////
				

				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
				}
				else
				{
					j0x=0;
					j0y=0;
					j0z=0;
				}
				
				for(int n=1;n<=4;n++)
				{
					div+=c[n]*j0x+d[n]*j0y+e[n]*j0z;
				}
			}
		}
		cout<<"�����d�����U="<<div<<endl;
	}
	else
	{
		//calc_current(CON,NODE,ELEM,EDGE,node,nelm,nedge,jnb,nei,branch_num,current);//�C������K�v����(boundary_condition)
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	cout<<"�d�����x�̐ݒ芮��"<<endl;

	if(CON->get_FEM_elm_type()==0)//�ߓ_�v�f
	{
		double *A[3];									//�ߓ_�ɂ������޸�����ݼ��
		for(int D=0;D<3;D++) A[D]=new double [node+1];
		double *old_A[3];		
		for(int D=0;D<3;D++) old_A[D]=new double [node+1]; //1step�O���޸�����ݼ��

		for(int i=1;i<=node;i++) for(int D=0;D<3;D++)
		{
			A[D][i]=0; old_A[D][i]=0;
		}
				
		/////old_A�ɒl���i�[
		if(t==1 && CON->get_restart()==OFF)				//�ŏ��̽ï�߂̓[���ŏ�����
		{
			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old_A[D][i]=0;	//������
		}
		else//����ȊO��̧�ق���ǂݍ���
		{	
			ifstream old("old_A.dat");
			if(!old) cout<<"cannot open old_A.dat"<<endl;
			old.unsetf(ifstream::dec);
			old.setf(ifstream::skipws);

			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];
					
			old.close();

			//�ǂݍ��񂾒l�ł̓����b�V���̈�����̓��̂̐ߓ_�ԍ�������Ă���̂ŁA���q�ɋL���������l��n��
			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) old_A[D][i]=PART[NODE[i].particleID].old_A[D];
		}

		//�x�N�g���|�e���V��������
		if(CON->get_m_A()==0) Avector3D_node_eddy2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t);
		//if(CON->get_m_A()==1) Avector3D_node_eddy2_jw(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�ߓ_���ƂɎ���������
		//if(CON->get_m_A()==2) Avector3D_node_eddy2_jw2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�����S���������S��
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==OFF) Avector3D_node_eddy2_jw3(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//���f���A
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==ON)
		{
			int dir=1;
			if(dir==0)
			{
				if(t==1) Avector3D_node_eddy2_jw3(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//���f��
				//if(t>1) Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//���f��,�ÓI�v�f�f�B���N����
				if(t>1) Avector3D_node_eddy2_jw5_ver2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw,node_sta);//���f��,�ÓI�v�f�f�B���N����
			}
			else
			{
				 //Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//�ȑO�̉�͂Ōv�Z�����ÓI�v�f�̒l���ŏ�����f�B���N���l�Ƃ��ė��p
				 Avector3D_node_eddy2_jw5_ver2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw,node_sta);//���f��,�ÓI�v�f�f�B���N����
			}

		}
		//if(CON->get_m_A()==4) Avector3D_node_eddy2_jw4(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�s�����A������œ���ւ�
	
		//�������x���� 
		Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be,t,ON);

		int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
		int *count_e=new int [node+1];	//�e���_�܂��ɂ����v�Z�Ώۂ̗v�f�����邩�B�􉽓I�Ȋ֌W�͂��łɋ��߂��Ă��邪�A���̋��E���̊O�Ɠ��̂ǂ���ŕ��ς��邩�ŕς���Ă���̂ł����ł��Ƃ߂�
		for(int i=1;i<=node;i++) count_e[i]=0;
		double B_sum=0;

		//�������x��v�f����ߓ_�ɕ��z����
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					N[j]=ELEM[je].node[j];
					if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
					{
						B_sum=sqrt((Be[A_X][je]*Be[A_X][je]+Be[A_Y][je]*Be[A_Y][je]+Be[A_Z][je]*Be[A_Z][je])/2.0);//�����l�ŏo��
						NODE[N[j]].B+=B_sum;
						count_e[N[j]]=count_e[N[j]]+1;
					}
				}
			}
		}
		for(int i=1;i<=node;i++)
		{
			if(count_e[i]!=0)
			{
				NODE[i].B/=count_e[i];
			}
		}
				
		//////

		//�ߓ_�͖@
		NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

		
		////
		if(CON->get_m_A()==1 && CON->get_jw_Faverage()==ON)//�P����������ɓ����͐ς�]�����邽�߁A�d���͂̕��ϒl�����߂�
		{
			cout<<"�P����������̓d���͂̕��όv�Z"<<endl;
			double TIME2;//�d���͂̎����͎����d���̎����̔����BTIME2��TIME����d���͂̎����̔��������i�߂��l
			//TIME2=TIME+0.000025/4;
			TIME2=TIME+1/(4*CON->get_Hz());
			cout<<"TIME2="<<TIME2<<endl;
			double *F2[DIMENTION];//���������ꂽ���Ԃ̓d���͂��i�[				
			for(int D=0;D<DIMENTION;D++) F2[D]=new double [fluid_number];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F2[D][i]=0;//������
			
			calc_jw_field_node(CON,NODE,ELEM,dt,TIME2,PART,fluid_number,F2,t,node,nelm);

			double Fsum=0;
			double Fsum2=0;
			 for(int i=0;i<fluid_number;i++)
			{
				Fsum+=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]+F[A_Z][i]*F[A_Z][i]);
				Fsum2+=sqrt(F2[A_X][i]*F2[A_X][i]+F2[A_Y][i]*F2[A_Y][i]+F2[A_Z][i]*F2[A_Z][i]);
			 }

			if(t==1)
			{
				ofstream fout("Fsum1.dat");
				fout.close();
			}

			ofstream fout2("Fsum1.dat",ios :: app);
			fout2<<Fsum<<endl;
			fout2.close();

			if(t==1)
			{
				ofstream fouts("Fsum2.dat");
				fouts.close();
			}

			ofstream fouts2("Fsum2.dat",ios :: app);
			fouts2<<Fsum2<<endl;
			fouts2.close();

			
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=(F[D][i]+F2[D][i])/2;//�d���͂̕��ϒl
			
			

			for(int D=0;D<DIMENTION;D++) delete [] F2[D];
		}
		////*/
				

		//�d���̓X���[�W���O�@fem_smn��-1�Ȃ�\�ʁA0�Ȃ�g��Ȃ��A1�Ȃ�S��
		smoothingF3D(CON,PART,fluid_number,F,t);


		delete [] count_e;
		for(int D=0;D<3;D++) delete [] A[D];
		for(int D=0;D<3;D++)delete [] old_A[D];
	}

	if(CON->get_FEM_elm_type()==1)//�ӗv�f
	{
		/////
		//�ӗv�f�쐬
		int *branch_num=new int[node+1];	//�e�ߓ_���אڂ���ߓ_��(�e�ߓ_���牄�т�ӂ̐�)
		int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_��
		int **nei2=new int* [node+1];		//�e�ߓ_�̗אڂ���ߓ_�ԍ��i�[(nei�͐ߓ_-�v�f�ł��邪�Anei2�͐ߓ_-�ߓ_)
		for(int i=1;i<=node;i++) nei2[i]=new int [max];
		
		////�ӗv�f���� (�ߓ_�v�f���g�p����ꍇ�ł��A�d�����x�����߂�Ƃ��ɕӗv�f���ق���)
		int KTJ=node;				//���̂��Ɠ��I�ߓ_(����)���i�[���Ȃ��Ƃ����Ȃ�����AKTJ�𑝉�
		KTJ+=CON->get_add_points();
		int KTE=12*KTJ;
		nedge=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,EDGE,branch_num,nei2,KTE,static_EDGE,t,node_sta);

		//if(ELEM[i].edge

		for(int i=1;i<=node;i++) delete [] nei2[i];
		delete [] nei2;
		delete [] branch_num;
		//////*/

		double *A=new double [nedge+1];	//�ӂɂ������޸�����ݼ��
		double *Am=new double [nedge+1];	//�ӂɂ������޸�����ݼ�ق̏u���l(�g��)	
		double *phi=new double [nedge+1];	//�ӂɂ������޸�����ݼ�ق̈ʑ��x��
		double *old_A=new double [nedge+1]; //1step�O���޸�����ݼ��

		for(int i=1;i<=nedge;i++) 
		{
			A[i]=0;	//�ӂɂ������޸�����ݼ��
			Am[i]=0;//�ӂɂ������޸�����ݼ�ق̏u���l(�g��)	
			phi=0;	//�ӂɂ������޸�����ݼ�ق̈ʑ��x��
			old_A=0; //1step�O���޸�����ݼ��
		}
				
		/*////old_A�ɒl���i�[
		if(t==1 && CON->get_restart()==OFF)				//�ŏ��̽ï�߂̓[���ŏ�����
		{
			for(int i=1;i<=node;i++) old_A[i]=0;	//������
		}
		else//����ȊO��̧�ق���ǂݍ���
		{	
			//�P�X�e�b�v�O�̕ӂ��ĉ��H //�����b�V�����Ă��Q�d���v�Z�����ό`���Ȃ��ꍇ�͍l�����邪�A���̌v�Z���ƁE�E�E
			ifstream old("old_A.dat");
			if(!old) cout<<"cannot open old_A.dat"<<endl;
			old.unsetf(ifstream::dec);
			old.setf(ifstream::skipws);

			for(int i=1;i<=nedge;i++) old>>old_A[i];
					
			old.close();

			//for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) old_A[D][i]=PART[NODE[i].particleID].old_A[D];
		}
		*/

		//Avector3D(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,jnb,branch_num,current,RP);
		//�x�N�g���|�e���V��������
		if(CON->get_m_A()==0) Avector3D_edge_eddy(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�ӗv�f
		//else if(CON->get_m_A()==1) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�ӗv�f�A����

		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==OFF) 
		{
			if(CON->get_parabolic_node_element()==OFF) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�ӗv�f�A����
			else if(CON->get_parabolic_node_element()==ON) Avector3D_edge_eddy_jw_with_parabolic_node_element(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);
		}
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==ON)
		{
			int dir=1;
			if(dir==0)
			{
				if(t==1) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//�ӗv�f�A����
				//if(t>1) Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//���f��,�ÓI�v�f�f�B���N����
				if(t>1) Avector3D_edge_eddy_jw2(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,node_sta,static_EDGE,Am,phi);//�ӗv�f�A���֐ÓI�v�f�f�B���N����
			}
			else
			{
				 //Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//�ȑO�̉�͂Ōv�Z�����ÓI�v�f�̒l���ŏ�����f�B���N���l�Ƃ��ė��p
				 Avector3D_edge_eddy_jw2(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,node_sta,static_EDGE,Am,phi);//�ӗv�f�A���֐ÓI�v�f�f�B���N����
			}

		}
		//�������x���� 
		Bflux3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,Be,t,ON);
		
		double *val=new double[nelm+1];
		for(int i=1;i<=nelm;i++)
		{
			double B=sqrt(Be[A_X][i]*Be[A_X][i]+Be[A_Y][i]*Be[A_Y][i]+Be[A_Z][i]*Be[A_Z][i]);
			val[i]=B;
		}

		/*
		if(t==1)
		{		
			int flux=0;//�Z�O�����g�f��
			data_avs2flux(CON,node,nelm,NODE,ELEM,val,t,flux);//�f�ʐ},����
			flux=1;//�X���b�g�f��
			data_avs2flux(CON,node,nelm,NODE,ELEM,val,t,flux);//�f�ʐ},����
		}
		*/

		delete [] val;



		//�ߓ_�͖@
		NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

		/////
		if(CON->get_m_A()==1 && CON->get_jw_Faverage()==ON)//�P����������ɓ����͐ς�]�����邽�߁A�d���͂̕��ϒl�����߂�
		{
			cout<<"�P����������̓d���͂̕��όv�Z"<<endl;
			double TIME2;//�d���͂̎����͎����d���̎����̔����BTIME2��TIME����d���͂̎����̔��������i�߂��l
			TIME2=TIME+1/(4*CON->get_Hz());
			cout<<"TIME2="<<TIME2<<endl;
			double *F2[DIMENTION];//���������ꂽ���Ԃ̓d���͂��i�[				
			for(int D=0;D<DIMENTION;D++) F2[D]=new double [fluid_number];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F2[D][i]=0;//������
			
			calc_jw_field_edge(CON,NODE,ELEM,EDGE,dt,TIME2,PART,fluid_number,F2,t,node,nelm,nedge,Am,phi);

			double Fsum=0;
			double Fsum2=0;
			 for(int i=0;i<fluid_number;i++)
			{
				Fsum+=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]+F[A_Z][i]*F[A_Z][i]);
				Fsum2+=sqrt(F2[A_X][i]*F2[A_X][i]+F2[A_Y][i]*F2[A_Y][i]+F2[A_Z][i]*F2[A_Z][i]);
			 }

			if(t==1)
			{
				ofstream fout("Fsum1.dat");
				fout.close();
			}

			ofstream fout2("Fsum1.dat",ios :: app);
			fout2<<Fsum<<endl;
			fout2.close();

			if(t==1)
			{
				ofstream fouts("Fsum2.dat");
				fouts.close();
			}

			ofstream fouts2("Fsum2.dat",ios :: app);
			fouts2<<Fsum2<<endl;
			fouts2.close();

			
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=(F[D][i]+F2[D][i])/2;//�d���͂̕��ϒl
			
			

			for(int D=0;D<DIMENTION;D++) delete [] F2[D];
		}
		////*/


		//�d���̓X���[�W���O�@fem_smn��-1�Ȃ�\�ʁA0�Ȃ�g��Ȃ��A1�Ȃ�S��
		smoothingF3D(CON,PART,fluid_number,F,t);

		delete [] A;
		delete [] Am;
		delete [] phi;
		delete [] old_A;

		
		
	}

	for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//���߂��d���͂�PART�Ɋi�[
	output_F_scalar_with_AVS_for_linear(CON,NODE,ELEM,t,PART,node);

	delete [] V;
	delete [] RP;
	delete [] sigma;
	

	for(int D=0;D<3;D++) delete [] current[D];
	for(int D=0;D<3;D++) delete [] Be[D];
}

//jw�@�p�v�Z�֐�(�ߓ_�v�f)
void calc_jw_field_node(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm)
{
	int node11=(int) NODE.size()-1;
	//int nelm=(int) ELEM.size()-1;
	//cout<<"node="<<node<<" NODE.size()="<<node11<<endl;

	ifstream a("Am.dat");
	if(!a) cout<<"cannot open Am.dat"<<endl;
	ifstream p("phi.dat");
	if(!p) cout<<"cannot open phi.dat"<<endl;

	double Am[3];
	double phi[3];

	double *A[3];									//�ߓ_�ɂ������޸�����ݼ��
	for(int D=0;D<3;D++) A[D]=new double [node+1];
	double *Be[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	for(int i=0;i<=node;i++) for(int D=0;D<3;D++) A[D][i]=0.0;

	double omega=2*PI*CON->get_Hz();

	a.unsetf(ifstream::dec);
	a.setf(ifstream::skipws);
	
	p.unsetf(ifstream::dec);
	p.setf(ifstream::skipws);			

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0;
			phi[D]=0;
		}
		for(int D=0;D<3;D++)
		{
			a>>Am[D];
			p>>phi[D];
			A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
		}
	}
	a.close();
	p.close();


	//�������x���� 
	Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be,t,OFF);

	double *RP=new double [nelm+1];	//�e�v�f�̓�����

	//����������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:�ގ��ɑ΂��A���������s��"<<endl;
	}

	int *jnb=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//�ߓ_�͖@
	NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

	//�d���̓X���[�W���O�@fem_smn��-1�Ȃ�\�ʁA0�Ȃ�g��Ȃ��A1�Ȃ�S��
	//smoothingF3D(CON,PART,fluid_number,F,t);

	for(int D=0;D<3;D++) delete [] A[D];
	delete [] RP;
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	for(int D=0;D<3;D++) delete [] Be[D];
}

//jw�@�p�v�Z�֐�(�ӗv�f)
void calc_jw_field_edge(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm,int nedge,double *Am, double *phi)
{
	//int node11=(int) NODE.size()-1;
	//int nelm=(int) ELEM.size()-1;
	//cout<<"node="<<node<<" NODE.size()="<<node11<<endl;

	//ifstream a("Am_e.dat");
	//if(!a) cout<<"cannot open Am_e.dat"<<endl;
	//ifstream p("phi_e.dat");
	//if(!p) cout<<"cannot open phi_e.dat"<<endl;

	//double Am;
	//double phi;

	
	double *A=new double [nedge+1];	//�ӂɂ������޸�����ݼ��
	double *Be[3];					//�v�f���d�E or �v�f�����E
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	for(int i=0;i<=node;i++) A[i]=0.0;

	double omega=2*PI*CON->get_Hz();

	//a.unsetf(ifstream::dec);
	//a.setf(ifstream::skipws);
	
	//p.unsetf(ifstream::dec);
	//p.setf(ifstream::skipws);			

	for(int i=1;i<=nedge;i++)
	{
		//Am=0;
		//phi=0;
		
		//a>>Am;
		//p>>phi;
		A[i]=Am[i]*cos(omega*TIME+phi[i]);
		
	}
	//a.close();
	//p.close();


	//�������x���� 
	Bflux3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,Be,t,OFF);

	double *RP=new double [nelm+1];	//�e�v�f�̓�����

	//����������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:�ގ��ɑ΂��A���������s��"<<endl;
	}

	int *jnb=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//�ߓ_�͖@
	NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

	//�d���̓X���[�W���O�@fem_smn��-1�Ȃ�\�ʁA0�Ȃ�g��Ȃ��A1�Ȃ�S��
	//smoothingF3D(CON,PART,fluid_number,F,t);

	delete [] A;
	delete [] RP;
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	for(int D=0;D<3;D++) delete [] Be[D];
}


///�d�����x�v�Z�֐�(�ӗv�f�p)
void calc_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,int *jnb,int **nei,int *branch_num,double **current)
{
	//II:�d��[A]

	cout<<"�d�����x�v�Z�J�n"<<endl;
	
	double p=1.68e-8;//���̓d�C��R��[��m]

	int side_num2=0;//�R�C�����\������Ӑ�
	for(int i=1;i<=nedge;i++)
	{
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL) side_num2++;
	}///side_num2�����Ƃ܂���

	int *side_id=new int [side_num2+1];//�R�C�����\������Ӕԍ��i�[
	side_num2=0;
	for(int i=1;i<=nedge;i++)
	{
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL)
		{
			side_num2++;
			side_id[side_num2]=i;//�Ӕԍ�i���i�[
		}
	}///�R�C�����\������Ӕԍ����Ǘ�


	//////////////////////////��ٕӂ̓d���Ɋւ��鋫�E������ݒ�
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//���̕\�ʂɗאڂ����C�v�f
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///��j�ʂɗאڂ���O�p�`�ƌ����������Ă��钸�_�́A��j�Ԑߓ_�ł���
						int p=ELEM[i].node[j];//���E�O�p�ɑ����Ȃ��ߓ_
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].edge[k];
							
							int ia=EDGE[iside].node[1];
							int ib=EDGE[iside].node[2];
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///����ŕK��ia<ib�ƂȂ���
							
							if(ia!=p && ib!=p)//��iside�̓R�C�����E��Ƃ�������
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									EDGE[iside].boundary_condition=11;//11���ЂƂł��܂񂾕ӂ͎��R���E����
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									EDGE[iside].boundary_condition=10;//T=0�ƂȂ��
								}
								//T=0�̖ʂ���ю��R���E�����ʂɌŒ苫�E�ӂ�ݒu����
								if(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) EDGE[iside].boundary_condition=21;//�ӂ�21��22�̕���
								else if(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) EDGE[iside].boundary_condition=22;//�ӂ�22��21�̕���
								
							}
						}
					}
				}
			}
		}
	}//�R�C���[�ʂɎ��R���E�A���ʂ�T=0�̌Œ苫�E��~�����B
	////////////////*/

	///���E�����o�́@���܂������Ȃ��Ƃ��ɂ݂�
	//data_avs_J_boundary(node,nelm,NODE,ELEM);

	for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
		//���R���E�����𖢒m������������
        if(EDGE[i].boundary_condition==11) EDGE[i].boundary_condition=0;
	}
	

	double II=CON->get_J0();
	double *T=new double [nedge+1];//�d���x�N�g�����ݼ��
    int NN=0;//�f�B���N���^���E�Ӑ�
    int *dn=new int [nedge+1]; //�e�ӂ��f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�side_num2+1���i�[
    double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
    
	for(int k=1;k<=nedge;k++)T[k]=0;
    ///�f�B���N���^���E��������
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(EDGE[i].boundary_condition==10)
		{
			dn[i]=NN;//i�Ԗڂ̕ӂ�NN�Ԗڂ̃f�B���N�����E��
	        PHAT[NN]=0;
	        T[i]=0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==21)
		{    
			int ia=EDGE[i].node[1];
			int ib=EDGE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==22)
		{   
			int ia=EDGE[i].node[1];
			int ib=EDGE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//�ӂ̒���
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else dn[i]=side_num2+1;
    }
	cout<<"�ިظڐ���"<<NN<<endl;
	/////////////*/

	
    //////int pn=side_num-NN;				///���m��
	int pn=side_num2-NN;				///���m��
    int *ppn=new int [pn];			//�s���n�Ԗڂ͕�ppn[n]
    int *npp=new int [nedge+1];	///�e�ӂ��s��̉��Ԗڂɂ����邩�B�ިظڌ^�̏ꍇ��pn+1���i�[
    int num=0; 
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(EDGE[i].boundary_condition==0)//���m�� 
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////�s��̍ő啝�v�Z
    int mat_w=0;
	for(int k=1;k<=side_num2;k++)
	{
		int i=side_id[k];
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		int width=branch_num[ia]+branch_num[ib];//�s��̕�
		if(width>mat_w) mat_w=width;
	}
	//mat_w*=5;
	////////////
	
    ////�z��m��
    double **G=new double *[pn+1];///�S�̍s��
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///�e�s�́A��[���v�f�̗�ԍ��L��
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///�e�s�́A��[���v�f��
    
    for(int i=1;i<=pn;i++)//������
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//���s��
    for(int i=0;i<pn;i++) B[i]=0;//������
    ////
    
    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
	
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
		
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//�v�f�̏d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;Ys/=4;Zs/=4;
			////////////////////////////
	
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
	        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//i�������̂Ƃ�
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
		
			////�v�f��ظ��쐬�J�n
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				if(EDGE[iside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
				{   
					int I=npp[iside]+1;///��iside�͍s���I�Ԗ�
					int I1=EDGE[iside].node[1];//iside���\������2�_
					int I2=EDGE[iside].node[2];
					int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
					int k2=table[i][2];
				    for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						
						int u1=table[j][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
						int u2=table[j][2];
							
						if(EDGE[jside].boundary_condition==0)///�ިظڂłȂ��B�܂薢�m�Ȃ�
						{   
							int J=npp[jside]+1;///��jside�͍s���J�Ԗ�
							int flag=0;
							//if(J<=I){
							int J1=EDGE[jside].node[1];//jside���\������2�_
							int J2=EDGE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//���łɓ���J���i�[����Ă���Ȃ�
								{
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
								    flag=1;
								}
							}
							if(flag==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
							    ROW[I][H]=J;
							}
							//}
						}
						////
						else //jside���ިظڌ^���E�ߓ_�Ȃ�
						{
						    int n=dn[jside];
						    B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6*PHAT[n];
						}//////////*/
					}
				}
			}
		}
    }
    ///////////////////////*/
	

    int number=0;//�s��̔�[���v�f��
    for(int i=1;i<=pn;i++) number+=NUM[i];
   
    ///���̂܂܂ł�G[i][j]��[j]�̏��������ɕ���ł��Ȃ��̂ŁAG��ROW����ёւ���
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////

	///�Ώ̐��`�F�b�N
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++) if(ROW[J][k]==i) if(G[i][j]!=G[J][k])
			{
				cout<<"matrix isn't symmetric   "<<G[i][j]<<" "<<G[J][k]<<endl;
			}
		}
	}///////////*/

    double *val = new double [number];
    int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
    int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
    
    /////////////////////val,ind ,ptr�ɒl���i�[
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
    cout<<"�s��쐬�I��  "<<endl;
    
    
	//CG3D(val,ind,ptr,pn,ppn,B,T);//CG�@���s
	//ICCG3D(val,ind,ptr,pn,ppn,B,T,number);//ICCG�@���s
	ICCG3D2(CON,val,ind,ptr,pn,B,number,T);
    ///////////
	
	denryu_side(CON,NODE,ELEM,EDGE,node,nelm,nedge,T,current);
    
	delete [] side_id;
	delete [] T;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

////�ӗv�f�d�����x�v�Z�֐�
void denryu_side(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *T,double **current)
{
	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W

	ofstream fp("j.dat");
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//�v�f�̏d�S���W
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;
			Ys/=4;
			Zs/=4;
			////////////////////////////
	
			double delta6=ELEM[je].volume;//�̐ς�6�{(�������̐ς̒l�͂��łɃx�N�g���|�e���V���������߂�ۂɌv�Z���Ă���)
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;
		    
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++) current[D][je]=0;//������
	
			for(int i=1;i<=6;i++)
			{
				int s=ELEM[je].edge[i];
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];
	
				current[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*T[s];
				current[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*T[s];
				current[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*T[s];
				
			}
	
			for(int D=0;D<3;D++) current[D][je]*=delta6*delta6*2;
			
			//if(Zs>0 && Zs<0.0001) fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;
			if(Zs>0.0005 && Zs<0.001)
			fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;		
		}
		
	}	
	fp.close();
}

//�ӗv�f�쐬�֐�
int make_edge_element(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,vector<edge3D> &EDGE,int *branch_num,int **nei2,int KTE,vector<edge3D> &static_EDGE,int t,int node_sta)
{
	unsigned timeA=GetTickCount();		//�v�Z�J�n����

	cout<<"�ӗv�f�����J�n"<<endl;
	cout<<"node="<<node<<" ele="<<nelm<<endl;
	////�ӗv�f����
	int *check=new int [node+1];		//�e�ߓ_������������������ǂ���
	int *flag=new int [node+1];		//�e�ߓ_������������������ǂ���
	int *temp_check=new int[node+1];	//�ꎞ�I�������z��
	int max=1000;						//�ߓ_�ɗאڂ���ő�ߓ_��
	int *ele_side_num=new int[nelm+1];	//�e�v�f�̑扽�ӂ܂ł����Ƃ܂��Ă��邩
	int side_num=0;						//�S�Ӑ��i�[
		
	///������
	for(int i=1;i<=node;i++) 
	{
		check[i]=0;
		temp_check[i]=0;
		branch_num[i]=0;//������
		flag[i]=0;
	}
	for(int i=1;i<=nelm;i++) ele_side_num[i]=0;
	/////////////////////

	//�Ӕԍ��ƕ�-�ߓ_��񐶐�
	for(int i=1;i<=node;i++)
	{
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//�ߓ_i�̗אڂ���ߓ_�����v���X
					temp_check[p]=1;//��������
					nei2[i][branch_num[i]]=p;
					if(check[p]==0)
					{	
						side_num++;
						EDGE[side_num].node[1]=i;//���̱ٺ�ؽ�щ��ɂ����ẮA���i<p�ł���.�Ȃ��Ȃ�i��菬���Ȕԍ��̐ߓ_�͂��ł�check=1������
						EDGE[side_num].node[2]=p;

						///�ߓ_�x�[�X�̋��E������Ӄx�[�X�Ɋg��
						if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//���[���������E�����Ȃ�A���̕ӂ����̋��E�����ɏ]���B���ɗ������m���ł����Ȃ�
						else EDGE[side_num].boundary_condition=0;//����ȊO�͖��m��
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//������
	}

	/*///
	//�Ӕԍ��ƕ�-�ߓ_��񐶐�//���̎��_�ŐÓI�v�f���\������ߓ_�͌��肵�Ă���A�ߓ_�ԍ����Ⴂ�����珇�ɐU���Ă���B�܂��A���̔ԍ��̓X�e�b�v���i�s���Ă��ω����Ȃ��͂��B�ÓI�ߓ_�݂̂ō\�������ӂ��珇�ɔԍ���U��
	for(int i=1;i<=node;i++)
	{		
		if(i<=node_sta) 
		{
			//flag[i]=ON;
			for(int k=1;k<=jnb[i];k++)
			{
				int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
				for(int j=1;j<=4;j++)
				{
					int p=ELEM[jelm].node[j];
					if(p<=node_sta)
					{
						//flag[p]=ON;
				
						if(p!=i && temp_check[p]==0)
						{
							branch_num[i]=branch_num[i]+1;//�ߓ_i�̗אڂ���ߓ_�����v���X
							temp_check[p]=1;//��������
							nei2[i][branch_num[i]]=p;
							if(check[p]==0)
							{	
								//if(flag[i]==ON && flag[p]==ON)
								{
									side_num++;
									EDGE[side_num].node[1]=i;//���̱ٺ�ؽ�щ��ɂ����ẮA���i<p�ł���.�Ȃ��Ȃ�i��菬���Ȕԍ��̐ߓ_�͂��ł�check=1������
									EDGE[side_num].node[2]=p;

									///�ߓ_�x�[�X�̋��E������Ӄx�[�X�Ɋg��
									if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//���[���������E�����Ȃ�A���̕ӂ����̋��E�����ɏ]���B���ɗ������m���ł����Ȃ�
									else EDGE[side_num].boundary_condition=0;//����ȊO�͖��m��
								}
							}
						}
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//������
	}

	for(int i=1;i<=node;i++) check[i]=0;

	for(int i=1;i<=node;i++)
	{		
		if(i<=node_sta) flag[i]=ON;
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p<=node_sta) flag[p]=ON;
				
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//�ߓ_i�̗אڂ���ߓ_�����v���X
					temp_check[p]=1;//��������
					nei2[i][branch_num[i]]=p;
					if(check[p]==0)
					{	
						if(flag[i]==OFF || flag[p]==OFF)
						{
							side_num++;
							EDGE[side_num].node[1]=i;//���̱ٺ�ؽ�щ��ɂ����ẮA���i<p�ł���.�Ȃ��Ȃ�i��菬���Ȕԍ��̐ߓ_�͂��ł�check=1������
							EDGE[side_num].node[2]=p;

							///�ߓ_�x�[�X�̋��E������Ӄx�[�X�Ɋg��
							if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//���[���������E�����Ȃ�A���̕ӂ����̋��E�����ɏ]���B���ɗ������m���ł����Ȃ�
							else EDGE[side_num].boundary_condition=0;//����ȊO�͖��m��
						}
					}
				}
				flag[p]=OFF;
			}
		}
		flag[i]=OFF;
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//������
	}
	///*/
	

	for(int i=1;i<=side_num;i++) EDGE[i].static_num=0;

	
	
	///////////*/

	///�v�f-�ӏ�񐶐�
	for(int i=1;i<=nelm;i++) for(int j=1;j<=6;j++) ELEM[i].edge[j]=0;
	for(int i=1;i<=side_num;i++)
	{
		int ia=EDGE[i].node[1];//��i���\������ߓ_�ԍ�(ia<ib)
		int ib=EDGE[i].node[2];
		for(int k=1;k<=jnb[ia];k++)
		{
			int jelm=nei[ia][k];//�ߓ_ia�̗אڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				if(ELEM[jelm].node[j]==ib)
				{
					ele_side_num[jelm]=ele_side_num[jelm]+1;
					int a=ele_side_num[jelm];
					ELEM[jelm].edge[a]=i;
				}
			}
		}
	}////////////

	//�v�f��remesh��񂩂�A�ӗv�f���ÓI�����I���𔻒f����
	//cout<<"�ÓI�ӗv�f�̐ݒ�J�n"<<endl;
	
	for(int i=1;i<=side_num;i++) EDGE[i].stat=OFF;//������
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].remesh==OFF)//�����b�V�����Ȃ��v�f�̎��ӂ̓����b�V������Ȃ� = �ÓI
		{
			for(int j=1;j<=6;j++)
			{
				EDGE[ELEM[i].edge[j]].stat=ON;
			}
		}
		
	}////////////
	
	
	/*/1step�̏ꍇ�Astatic_edge�ɋL�^���Ă���
	if(t==1)
	{
		int count=0;
		static_EDGE.resize(side_num+1);
		for(int i=1;i<=side_num;i++)
		{
			EDGE[i].static_num=i;
			static_EDGE[i]=EDGE[i];
			if(static_EDGE[i].stat==ON) count++;
		}
		static_EDGE[0].static_num=count;//�ÓI�ӗv�f�̐����i�[
		//cout<<"side_num="<<side_num<<" static_EDGE.size()="<<static_EDGE.size()-1<<endl;
	}

	if(t>1)
	{
		int count=0;
		
		for(int i=1;i<=static_EDGE.size()-1;i++)
		{
			if(static_EDGE[i].stat==ON)
			{
				//if(i%50000==0) cout<<i<<endl;
				if((NODE[EDGE[i].node[1]].r!=NODE[static_EDGE[i].node[1]].r) || (NODE[EDGE[i].node[2]].r!=NODE[static_EDGE[i].node[2]].r) )
				{
					//cout<<i<<endl;
					for(int j=1;j<=side_num;j++)
					{
						if(NODE[EDGE[j].node[1]].r==NODE[static_EDGE[i].node[1]].r && NODE[EDGE[j].node[2]].r==NODE[static_EDGE[i].node[2]].r)//1�X�e�b�v�ڂƓ���̕ӁB�ÓI�v�f�����ɂ������݂��Ȃ��͂�
						{
							EDGE[j].static_num=i;
						}
					}
				}
				else  EDGE[i].static_num=i;
				count++;
				if(count==static_EDGE[0].static_num) break;
			}
		}
	}

	*/

	delete [] check;
	delete [] flag;
	delete [] temp_check;
	delete [] ele_side_num;

	//cout<<"�Ӑ�="<<side_num<<" �ő�Ӑ�="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	cout<<"�Ӑ�="<<side_num<<" �ő�Ӑ�="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	return side_num;
}

///���E�����K�p�֐�(�R�c�ӗv�f�p)
void set_boundary_condition3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A)
{
	int N=0;	//�����グ�ϐ�

	///�f�B���N���^���E��������
	if(CON->get_uniform_B_sw()==OFF)
	{
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{    
		        dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else if(SIDE[i].boundary_condition==2)
			{   
				/*int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				A[i]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        dn[i]=NN;
		        PHAT[NN]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        //A[i]=0;
		        NN++;*/
				dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else dn[i]=side_num+1;
		}
	}
	if(CON->get_uniform_B_sw()==ON)	//��͗̈�S�̂Ɉ�l�����^����ꍇ
	{
		double B=CON->get_uniform_B();//��l����̑傫��[�s�n
		double R[3];					//�ӂ̒����i�[(X,Y,Z����)
		double r[3];					//�ӂ̒��_���W�i�[
		double err=1e-14;
		
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{   
				int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				for(int D=0;D<3;D++)
				{
					R[D]=NODE[ib].r[D]-NODE[ia].r[D];
					r[D]=(NODE[ib].r[D]+NODE[ia].r[D])*0.5;//���_
				}

			//	if(r[A_Z]<CON->get_ZU()-err && r[A_Z]>CON->get_ZD()+err)
				{

					double L=sqrt(R[A_X]*R[A_X]+R[A_Y]*R[A_Y]+R[A_Z]*R[A_Z]);//�ӂ̒���
				
					A[i]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);	//�Ȃ������Ȃ�̂��̓X�g�[�N�X�̒藝�𗘗p����΂킩��B
					dn[i]=N;
					PHAT[N]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);
					N++;
				}
				//else 
				//{
				////	SIDE[i].boundary_condition=0;//���R���E����
				//	dn[i]=side_num+1;
				//}
			}
			else dn[i]=side_num+1;
		}
	}

	*NN=N;
}

//�d�����x�ǂݍ��݊֐�
void inport_J0_density(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	cout<<"current_density.txt��苭���d�����z��ǂݍ���--";
	int id;

	ifstream fp("current_density.dat");
	if(!fp) cout<<"cannot open current_density.dat"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);

	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==COIL)						//current_density.dat�ɂ́A�R�C���v�f�̂ݏo�͂���Ă���d�l�ɂ��Ă�������
		{												//���̏ꍇ�A�R�C���v�f�͐ÓI�v�f�łȂ���΂Ȃ�Ȃ��B���I�Ȃ瑼�̃\�t�g����ǂݍ��߂Ȃ�
			fp>>id;
			for(int D=0;D<3;D++)
			{
				fp>>current[D][i];		//
			}
		}
		else for(int D=0;D<3;D++) current[D][i]=0;
	}
	fp.close();
	cout<<"ok"<<endl;
}

//�d�����x�o�͊֐�
void check_J0(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//�R�C���v�f��
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) coil_num++;

	ofstream fout2("J0.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0 filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0 filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0 filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0 filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0 filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0 filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1;//1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}

void check_J0Je(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//�R�C���v�f��
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL ||ELEM[i].material==FLUID) coil_num++;

	ofstream fout2("J0Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0Je");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1;//1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL ||ELEM[i].material==FLUID)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}


////�ߓ_�v�f�������x�v�Z�֐�
void Bflux3D_node(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B,int t, int flag)
{
	cout<<"�������x�v�Z----------";

	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	double times=CON->get_B_times();
	double le=CON->get_distancebp();
	int plot_type=CON->get_plot_B_type();//1:�޸�ف@2:�X�J���[
	//ofstream fp2("test.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) B[D][je]=0;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		double delta6=ELEM[je].volume;//�̐ς�6�{

		delta6=1/delta6;

		double Xs=0;
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];Ys+=Y[j];Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		//�W���쐬
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int i=1;i<=4;i++)
		{
			B[A_X][je]+=d[i]*A[A_Z][N[i]]-e[i]*A[A_Y][N[i]];
			B[A_Y][je]+=e[i]*A[A_X][N[i]]-c[i]*A[A_Z][N[i]];
			B[A_Z][je]+=c[i]*A[A_Y][N[i]]-d[i]*A[A_X][N[i]];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6;

		//if(Ys>-5*le && Ys<5*le) fp2<<Xs<<" "<<Zs<<" "<<B[A_X][je]*100<<" "<<B[A_Z][je]*100<<endl;

	}
	cout<<"ok"<<endl;
//	fp2.close();

	//�������x�o��
	if(flag==ON)
	{
		int flagB=OFF;
		if(CON->get_m_A()==0) flagB=ON;
		if(CON->get_m_A()==1)
		{
			flagB=ON;
		}
		cout<<"�������x�o�͊J�n----";
		if(flagB==ON)
		{
			//�������x�o��
			ofstream fp("Bflux.dat");
			
			//double Xmin=CON->get_XL()+le; double Xmax=CON->get_XR()-le;//��͗̈�
			//double Zmin=CON->get_ZD()+le; double Zmax=CON->get_ZU()-le;

			double Xmin=CON->get_XL()/2+le; double Xmax=CON->get_XR()/2-le;//��͗̈�
			double Zmin=0.05+le; double Zmax=CON->get_ZU()-le;
			
			double Rmax=CON->get_RU()-le;
			
			double dx; 
			if(plot_type==1) dx=1*le;
			if(plot_type==2) dx=0.01;
			//double dx=0.01;
			int Nx=(int)((Xmax-Xmin)/dx);//�e�����̕�����
			int Nr=(int)(2*Rmax/dx);
			int Nz;
			if(plot_type==1) Nz=(int)((Zmax-Zmin)/dx);
			if(plot_type==2) Nz=15;
			int serch=nelm;//locate�֐��ōŏ��ɒT������v�f�ԍ�
			
			if(CON->get_region_shape()==1)		//�~���̈�̂Ƃ��͕ϐ������������ď���
			{
				Nx=Nr;
				Xmin=-Rmax;
			}


			if(plot_type==1)		//�޸�ٕ\��
			{
				for(int n=0;n<Nx;n++)
				{
					for(int m=0;m<Nz;m++)
					{
						double xp=dx*n+Xmin;//�o�͂���_�̍��W
						double yp=0;
						double zp=dx*m+Zmin;
						int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
						fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
						serch=loc;
						//if(loc==0) cout<<"EE"<<endl;
					}
				}
			}
			else if(plot_type==2)	//�X�J���[�\��
			{
				for(int n=0;n<1;n++)
				{
					for(int m=0;m<Nz;m++)
					{
						double xp=0;//�o�͂���_�̍��W
						//double xp=dx*n+Xmin;//�o�͂���_�̍��W
						double yp=0;
						double zp=dx*m+0.101;
						int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
						double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc])/sqrt(2.0);//�����l�o��
						fp<<xp<<" "<<zp<<" "<<BB<<endl;
						serch=loc;
					}
				}
			}
			fp.close();//*/
		

			/////microAVS�p�̎������x�o��

			int count=0;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)// &&ELEM[i].remesh==OFF)
				{
					count++;
				}
			}

			ofstream fout2("Bflux.fld");
			fout2 << "# AVS field file" << endl;
			fout2 << "ndim=1" << endl;
			//fout2 << "dim1=" << fluid_number <<endl;
			fout2 << "dim1=" << count<<endl;
			fout2 << "nspace=3" << endl;
			fout2 << "veclen=3" << endl;
			fout2 << "data=float" << endl;
			fout2 << "field=irregular" << endl;
			fout2 << "label=e-x e-y e-z" << endl << endl;
			fout2 << "variable 1 file=./Bflux filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout2 << "variable 2 file=./Bflux filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout2 << "variable 3 file=./Bflux filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout2 << "coord    1 file=./Bflux filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout2 << "coord    2 file=./Bflux filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout2 << "coord    3 file=./Bflux filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout2.close();

			ofstream fout("Bflux");
			fout<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)// && ELEM[i].material==OFF)
				{
					fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout.close();

			//�X���b�g��
			int count2=0;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
				{
					count2++;
				}
			}

			ofstream fout3("Bfluxs.fld");
			fout3 << "# AVS field file" << endl;
			fout3 << "ndim=1" << endl;
			//fout2 << "dim1=" << fluid_number <<endl;
			fout3 << "dim1=" << count2<<endl;
			fout3 << "nspace=3" << endl;
			fout3 << "veclen=3" << endl;
			fout3 << "data=float" << endl;
			fout3 << "field=irregular" << endl;
			fout3 << "label=e-x e-y e-z" << endl << endl;
			fout3 << "variable 1 file=./Bfluxs filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout3 << "variable 2 file=./Bfluxs filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout3 << "variable 3 file=./Bfluxs filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout3 << "coord    1 file=./Bfluxs filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout3 << "coord    2 file=./Bfluxs filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout3 << "coord    3 file=./Bfluxs filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout3.close();

			ofstream fout4("Bfluxs");
			fout4<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
				if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
				{
					fout4<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout4.close();

			int flag2=0;
			if(CON->get_EM_interval()>1) flag2=ON;
			//else if(t==1 || t%10==0) flag=ON;


			///////////////////////////////
			if(flag2==ON)
			{
				char filename[25];
				sprintf_s(filename,"Bflux%d.fld", t);
				ofstream fout2(filename);
				fout2 << "# AVS field file" << endl;
				fout2 << "ndim=1" << endl;
				fout2 << "dim1=" << count <<endl;
				//fout2 << "dim1=" << BOnum <<endl;
				fout2 << "nspace=3" << endl;
				fout2 << "veclen=3" << endl;
				fout2 << "data=float" << endl;
				fout2 << "field=irregular" << endl;
				fout2 << "label=e-x e-y e-z" << endl << endl;
				fout2 << "variable 1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
				fout2 << "variable 2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
				fout2 << "variable 3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
				fout2 << "coord    1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
				fout2 << "coord    2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
				fout2 << "coord    3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
				fout2.close();

				char filename2[25];
				sprintf_s(filename2,"Bflux%d", t);
				ofstream fout(filename2);
				fout<<"e-x e-y e-z x y z"<<endl;
				for(int i=1;i<=nelm;i++)
				{
					double r[3]={0,0,0};
					for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
					if(r[A_Y]>-0.1*le && r[A_Y]<0.1*le)
					//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
					{
						fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
					}
				}
				fout.close();
			}
		}
		cout<<"ok"<<endl;
	}

}

//�ӗv�f�p�������x�v�Z
void Bflux3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *A,double **B,int t,int flag)
{
	cout<<"�������x�v�Z�J�n----";
	unsigned timeA=GetTickCount();
	int plot_type=CON->get_plot_B_type();	//̧�ُo�͌`���@1=�޸�� 2=�X�J���[
	double times=CON->get_B_times();		//̧�ُo�͎��̔{��
	double u0=4*PI*1e-7;					//�^��̓�����
	double le=CON->get_distancebp();

	double *Xg=new double [nelm+1];			//�v�f�̏d�S���W
	double *Yg=new double [nelm+1];
	double *Zg=new double [nelm+1];


	//#pragma omp parallel for
	for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<" "<<omp_get_thread_num()<<endl;//�ei�̌v�Z��S�����Ă���گ�ޔԍ��o��
		int N[4+1];								//�v�f�̊e�ߓ_�ԍ��i�[
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		int table[6+1][2+1];					//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W

		//�Ӂ|�ߓ_ð��ٍ쐬
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//�v�f�̏d�S���W
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
		}
		Xg[je]=Xs; Yg[je]=Ys; Zg[je]=Zs;	//�d�S���
		////////////////////////////

		double delta6=ELEM[je].volume;//�̐ς�6�{(�������̐ς̒l�͂��łɃx�N�g���|�e���V���������߂�ۂɌv�Z���Ă���)
	
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
	
		double delta=ELEM[je].volume/6;//�{���̑̐�
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
			if(i%2!=0)//i����Ȃ�
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int D=0;D<3;D++) B[D][je]=0;//������

		for(int i=1;i<=6;i++)
		{
			int s=ELEM[je].edge[i];
			int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
			int k2=table[i][2];

			B[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s];
			B[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
			B[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6*delta6*2;
	}

	if(flag==ON)
	{
		//�������x�o��
		ofstream fp("Bflux.dat");
	
		//double Xmin=CON->get_XL()+le; double Xmax=CON->get_XR()-le;//��͗̈�
		//double Zmin=CON->get_ZD()+le; double Zmax=CON->get_ZU()-le;

		double Xmin=CON->get_XL()/2+le; double Xmax=CON->get_XR()/2-le;//��͗̈�
		double Zmin=0.05+le; double Zmax=CON->get_ZU()-le;
	
		double Rmax=CON->get_RU()-le;
	
		//double dx=1*le;
		//int Nx=(int)((Xmax-Xmin)/dx);//�e�����̕�����
		//int Nr=(int)(2*Rmax/dx);
		//int Nz=(int)((Zmax-Zmin)/dx);
		//int serch=nelm;//locate�֐��ōŏ��ɒT������v�f�ԍ�
		double dx; 
		if(plot_type==1) dx=1*le;
		if(plot_type==2) dx=0.01;
		//double dx=0.01;
		int Nx=(int)((Xmax-Xmin)/dx);//�e�����̕�����
		int Nr=(int)(2*Rmax/dx);
		int Nz;
		if(plot_type==1) Nz=(int)((Zmax-Zmin)/dx);
		if(plot_type==2) Nz=15;
		int serch=nelm;//locate�֐��ōŏ��ɒT������v�f�ԍ�
	
		if(CON->get_region_shape()==1)		//�~���̈�̂Ƃ��͕ϐ������������ď���
		{
			Nx=Nr;
			Xmin=-Rmax;
		}


		if(plot_type==1)		//�޸�ٕ\��
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//�o�͂���_�̍��W
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
					serch=loc;
					//if(loc==0) cout<<"EE"<<endl;
				}
			}
		}
		else if(plot_type==2)	//�X�J���[�\��
		{
			for(int n=0;n<1;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=0;//�o�͂���_�̍��W
					//double xp=dx*n+Xmin;//�o�͂���_�̍��W
					double yp=0;
					double zp=dx*m+0.101;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc])/sqrt(2.0);//�����l�o��
					fp<<xp<<" "<<zp<<" "<<BB<<endl;
					serch=loc;
				}
			}
		}
		/*//////
		else if(plot_type==2)	//�X�J���[�\��
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//�o�͂���_�̍��W
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
					fp<<xp<<" "<<zp<<" "<<BB<<endl;
					serch=loc;
				}
			}
		}
		/////*/
		fp.close();//*/

		/////microAVS�p�̎������x�o��

		int count=0;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
			//if(ELEM[i].material==CRUCIBLE)
			if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			//if(r[A_Y]>-0.01 && r[A_Y]<0.01)
			{
				count++;
			}
		}

		ofstream fout2("Bflux.fld");
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count<<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Bflux filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Bflux filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Bflux filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Bflux filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Bflux filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Bflux filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		ofstream fout("Bflux");
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			//if(ELEM[i].material==CRUCIBLE)
			//if(r[A_Y]>-0.01 && r[A_Y]<0.01)
			{
				fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			}
		}
		fout.close();

		//�X���b�g��
		int count2=0;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
			if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
			{
				count2++;
			}
		}

		ofstream fout3("Bfluxs.fld");
		fout3 << "# AVS field file" << endl;
		fout3 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout3 << "dim1=" << count2<<endl;
		fout3 << "nspace=3" << endl;
		fout3 << "veclen=3" << endl;
		fout3 << "data=float" << endl;
		fout3 << "field=irregular" << endl;
		fout3 << "label=e-x e-y e-z" << endl << endl;
		fout3 << "variable 1 file=./Bfluxs filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout3 << "variable 2 file=./Bfluxs filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout3 << "variable 3 file=./Bfluxs filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout3 << "coord    1 file=./Bfluxs filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout3 << "coord    2 file=./Bfluxs filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout3 << "coord    3 file=./Bfluxs filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout3.close();

		ofstream fout4("Bfluxs");
		fout4<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
			{
				fout4<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			}
		}
		fout4.close();

		int flag2=0;
		if(CON->get_EM_interval()>1) flag2=ON;
		else if(t==1 || t%10==0) flag2=ON;


		///////////////////////////////
		if(flag2==ON)
		{
			char filename[25];
			sprintf_s(filename,"Bflux%d.fld", t);
			ofstream fout2(filename);
			fout2 << "# AVS field file" << endl;
			fout2 << "ndim=1" << endl;
			fout2 << "dim1=" << count <<endl;
			//fout2 << "dim1=" << BOnum <<endl;
			fout2 << "nspace=3" << endl;
			fout2 << "veclen=3" << endl;
			fout2 << "data=float" << endl;
			fout2 << "field=irregular" << endl;
			fout2 << "label=e-x e-y e-z" << endl << endl;
			fout2 << "variable 1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout2 << "variable 2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout2 << "variable 3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout2 << "coord    1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout2 << "coord    2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout2 << "coord    3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout2.close();

			char filename2[25];
			sprintf_s(filename2,"Bflux%d", t);
			ofstream fout(filename2);
			fout<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
				{
					fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout.close();
		}
	}

	delete [] Xg;
	delete [] Yg;
	delete [] Zg;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		
}

///�ߓ_�͖@�v�Z�֐�
void NODE_F3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpsparticle> &PART,double **F,int fluid_number,int t)
{
    cout<<"�ߓ_�͖@�ɂ��d���͌v�Z--------";
    double ep0=8.854e-12;	//�^��̗U�d���B
    double u0=12.5e-7;		//�^��̓�����
    int N[4+1];				//�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double Fz=0;			//Z�����̍���[N]
	double Fsum=0;			//����[N]
	unsigned timeA=GetTickCount();//�v�Z�J�n����

	double *Fn[3];
	for(int D=0;D<3;D++) Fn[D]=new double [node+1];//NN��node�ɕύX�B
    
    //for(int i=1;i<=node;i++) 
	if(CON->get_EM_calc_type()==3) //������
    {
        for(int I=1;I<=node;I++)
        {
			if(NODE[I].material==FLUID)
			{
			 for(int D=0;D<3;D++) Fn[D][I]=0;//������
			
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//�ߓ_i���אڂ���v�f�ԍ�
					
					//if(ELEM[jelm].material==AIR){
					///�}�N�X�E�F���̉��̓e���\��
					double Txx=(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=2*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=2*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=2*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
				
					/////�W��c,d,e�v�Z
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
	    				X[j]=NODE[N[j]].r[A_X];
	    				Y[j]=NODE[N[j]].r[A_Y];
	    				Z[j]=NODE[N[j]].r[A_Z];
					}
					int i=0;///�ߓ_i�͗v�fjelm�̑�J�Ԗڂ̐ߓ_
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
					int m=j%4+1;
					int n=m%4+1;
					//delta6�͑��E�����̂ł���Ȃ�
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
					if( i & 1 )//i����Ȃ�
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		    
					double u=RP[jelm];
						
					Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e)/(2*u);
					Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e)/(2*u);
					Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e)/(2*u);
					//}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/(6.00000000000*u0);
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;
				Fz+=Fn[A_Z][I];
				Fsum+=sqrt(Fn[A_X][I]*Fn[A_X][I]+Fn[A_Y][I]*Fn[A_Y][I]+Fn[A_Z][I]*Fn[A_Z][I]);
				
			}
		}
	}
    //cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*//////
	if(t==1)
	{
		ofstream fout("Fz.dat");
		fout.close();
	}

	ofstream fout2("Fz.dat",ios :: app);
	fout2<<Fz<<endl;
	fout2.close();

	if(t==1)
	{
		ofstream fouts("Fsum.dat");
		fouts.close();
	}

	ofstream fouts2("Fsum.dat",ios :: app);
	fouts2<<Fsum<<endl;
	fouts2.close();
    /////*/

	/*
	int *bound=new int[node+1];
	for(int n=1;n<=node;n++)
	{
		bound[n]=OFF;
		if(NODE[n].material==FLUID)
		{
			for(int k=1;k<=jnb[n];k++)
			{
				int jelm=nei[n][k];//�ߓ_i���אڂ���v�f�ԍ�
				if(ELEM[jelm].material==AIR) bound[n]=ON;
			}
		}
	}
	*/
	
	///F�X�V
	for(int n=1;n<=node;n++)
	{
		if(NODE[n].material==FLUID)
		{
			//if(bound[n]==ON)
			{
				int i=NODE[n].particleID;//i�Ԗڂ̗��̗��q�́An�Ԗڂ̐ߓ_�ɑ���
				if(i>=0) for(int D=0;D<3;D++) F[D][i]=Fn[D][n]; //-1�Ԗڂ̔z��ɃA�N�Z�X���邨���ꂪ����̂ŏ����t��
				
			}
		}
	}
	///�ߓ_�͖@�̗͂̒P�ʂ�[N]�Ȃ̂ŁA�������炳��ɗ��q���ߓ_�̗��q�Ɋւ��Ă��͂����Ƃ߂Ă��K�v�͂Ȃ�

	//delete [] bound;

	/*///�t�@�C���o��
	ofstream fp("Fn.dat");
	double le=CON->get_distancebp();
	double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
		}
	}
	//�}��o��
	fp<<0.015<<" "<<0.155" "<<1.0e-002*times<<" "<<0<<endl;

    fp.close();///////////////////

	////�t�@�C���o��//(�X���b�g��ʂ�f��)
	ofstream fps("Fnslit.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Y]>-le*0.5+sin(PI/24)*NODE[i].r[A_X] && NODE[i].r[A_Y]<le*0.5+sin(PI/24)*NODE[i].r[A_X])
			{
				fps<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<(Fn[A_X][i]*cos(PI/24)+Fn[A_Y][i]*sin(PI/24))*times<<" "<<Fn[A_Z][i]*times<<endl;
			}
		}
	}
    fps.close();///////////////////

	////�t�@�C���o��(Z�f��)
	ofstream fpz("Fnz.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Z]>-le*0.5+0.14125 && NODE[i].r[A_Z]<le*0.5+0.14125) fpz<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Y][i]*times<<endl;	
		}
	}
    fpz.close();///////////////////

	//if(t=1 || t%10==0)
	if(CON->get_EM_interval()>1)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
		{
			if(NODE[i].material==FLUID)
			{
				if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
			}
		}
		fp.close();///////////////////
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=1;i<=node;i++)//���̐ߓ_�̂ݏo��
		{
			if(NODE[i].material==FLUID)
			{
				if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
			}
		}
		fp.close();///////////////////
	}
	///*/

	for(int D=0;D<3;D++) delete [] Fn[D];
	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

}

void calc_eddy_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **old_A,double dt,double *V,double **Je,int t, double *sigma)
{
	//double sigma=CON->get_ele_conduc();//�d�C�`����
	double H=CON->get_height();			//̧�ُo�͗p�����p�����[�^
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-��(dA/dt+grad��) ---(��)

	cout<<"�Q�d�����x�v�Z----------";

	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double dA[3];//�޸�����ݼ�ق̍������i�[�iX,Y,Z)
	double gradV[3];//grad�ӊi�[

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) Je[D][je]=0;

		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			//�W���쐬
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++)
			{
				dA[D]=0;
				gradV[D]=0;
			}
			for(int j=1;j<=4;j++)
			{
				for(int D=0;D<3;D++) dA[D]+=(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(A[D][N[j]]-old_A[D][N[j]]);
			}
			for(int D=0;D<3;D++) dA[D]*=delta6/dt;//���i���j�̑�1���v�Z����

			for(int j=1;j<=4;j++)//grad�ӂ̌v�Z
			{
				gradV[A_X]+=c[j]*V[N[j]];
				gradV[A_Y]+=d[j]*V[N[j]];
				gradV[A_Z]+=e[j]*V[N[j]];
			}	
			for(int D=0;D<3;D++) gradV[D]*=delta6;//���i���j�̑�2���v�Z����

			//for(int D=0;D<3;D++) Je[D][je]=-sigma*(dA[D]+gradV[D]);
			for(int D=0;D<3;D++) Je[D][je]=-sigma[je]*(dA[D]+gradV[D]);

			//if(Zs>H && Zs<H+le) fout<<Xs<<" "<<Ys<<" "<<Je[A_X][je]*times<<" "<<Je[A_Y][je]*times<<endl;
			//sqrt(jx*jx+jy*jy)
		}
	}
	
	//�o��
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			count2++;
		}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count/*nelm*/ <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2/*nelm*/ <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///�ǂݍ��ݗp�t�@�C���쐬
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count/*nelm*/ <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count/*nelm*/ <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}



	cout<<"ok"<<endl;
	//fout.close();

}

void calc_node_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **AR,double **AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t,double TIME, double *sigma,double omega)
{
	//double sigma=CON->get_ele_conduc();//�d�C�`����
	double H=CON->get_height();			//̧�ُo�͗p�����p�����[�^
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-��(dA/dt+grad��) ---(��)

	cout<<"�Q�d�����x�v�Z(jw�@)----------";
	
	complex<double> *Jec[3];
	for(int D=0;D<3;D++) Jec[D]=new complex<double> [nelm+1];
	
    for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Jec[D][i]=complex<double>(0.0,0.0);//������

	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
//	double dA[3];//�޸�����ݼ�ق̍������i�[�iX,Y,Z)
	double temp_AR[3];
	double temp_AI[3];
	double temp_VR[3];
	double temp_VI[3];

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) Je[D][je]=0;

		if(ELEM[je].material==FLUID)// || ELEM[je].material==CRUCIBLE)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double co=omega*sigma[je]*delta6;
			double co2=sigma[je]*delta6;
			//�W���쐬
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++)
			{
				 temp_AR[D]=0;
				 temp_AI[D]=0;
				 temp_VR[D]=0;
				 temp_VI[D]=0;
			}
			for(int j=1;j<=4;j++)
			{

				for(int D=0;D<3;D++)
				{
					temp_AR[D]+=co*((b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(AR[D][N[j]]));
					temp_AI[D]+=co*((b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(AI[D][N[j]]));
				}
				
				temp_VR[A_X]+=co2*c[j]*VR[N[j]];
				temp_VR[A_Y]+=co2*d[j]*VR[N[j]];
				temp_VR[A_Z]+=co2*e[j]*VR[N[j]];
				temp_VI[A_X]+=co2*c[j]*VI[N[j]];
				temp_VI[A_Y]+=co2*d[j]*VI[N[j]];
				temp_VI[A_Z]+=co2*e[j]*VI[N[j]];
				
				/*
				temp_VR[A_X]+=co*c[j]*VR[N[j]];
				temp_VR[A_Y]+=co*d[j]*VR[N[j]];
				temp_VR[A_Z]+=co*e[j]*VR[N[j]];
				temp_VI[A_X]+=co*c[j]*VI[N[j]];
				temp_VI[A_Y]+=co*d[j]*VI[N[j]];
				temp_VI[A_Z]+=co*e[j]*VI[N[j]];
				*/
			}
				
			//for(int D=0;D<3;D++) gradV[D]*=delta6;//���i���j�̑�2���v�Z����
			for(int D=0;D<3;D++)
			{
				double t_R=-temp_AI[D]+temp_VR[D];
				double t_I=temp_AR[D]+temp_VI[D];
				//double t_R=-temp_AI[D]-temp_VI[D];
				//double t_I=temp_AR[D]+temp_VR[D];
				Jec[D][je]=complex<double>(-t_R,-t_I);
			}

		}
			
	}

	double Jem[3];
	double phi[3];

	for(int je=1;je<=nelm;je++)
	{
		Je_loss[je]=0.0;
		for(int D=0;D<3;D++)
		{
			Jem[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		if(ELEM[je].material==FLUID)
		{
			double Jem_sum=0.0;
			for(int D=0;D<3;D++)
			{
				double Re=Jec[D][je].real();
				double Im=Jec[D][je].imag();
				Jem[D]=sqrt(Re*Re+Im*Im);
				Jem_sum+=Jem[D]*Jem[D];
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(Im,Re);
				Je[D][je]=Jem[D]*cos(omega*TIME+phi[D]);
			}
			Jem_sum=sqrt(Jem_sum);
			//�Q�d�����̌v�Z
			Je_loss[je]=Jem_sum*Jem_sum*ELEM[je].volume/(2*sigma[je]*6);
		}
		//a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		//p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	//a.close();
	//p.close();

	cout<<"ok"<<endl;

	//
	//�Q�d����v�f����ߓ_�ɕ��z����(�o�͗p)
	//int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
	//double X[5];
	//double Y[5];
	//double Z[5];
	double Je_sum=0;

	
	int *count_e=new int [node+1];	//�e���_�܂��ɂ����v�Z�Ώۂ̗v�f�����邩�B�􉽓I�Ȋ֌W�͂��łɋ��߂��Ă��邪�A���̋��E���̊O�Ɠ��̂ǂ���ŕ��ς��邩�ŕς���Ă���̂ł����ł��Ƃ߂�

	//�Q�d�����x��v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_sum=sqrt((Je[A_X][je]*Je[A_X][je]+Je[A_Y][je]*Je[A_Y][je]+Je[A_Z][je]*Je[A_Z][je])/2.0);//�����l�ŏo��
					NODE[N[j]].Je+=Je_sum;
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			NODE[i].Je/=count_e[i];
		}
	}
	/*/////////
	//�o��
	cout<<"�Q�d�����x�o��------";
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			count2++;
		}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count<<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2 <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///�ǂݍ��ݗp�t�@�C���쐬
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	cout<<"ok"<<endl;
	/////*/



	
	//fout.close();
	for(int D=0;D<DIMENTION;D++) delete [] Jec[D];
	delete [] count_e;

}

//�ӗv�f���g��������͗p�Q�d���v�Z�֐��@�����@�O�����L���v�f�@�@p79�̎�(4.36)��jw���������̊O�ɏo�Ă��邪�ԈႢ�B�{����Ak�̍������ɂ����邱�Ƃɒ���
void calc_edge_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE, int node,int nelm,double *AR,double *AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t,double TIME, double *sigma,double omega)
{
	//double sigma=CON->get_ele_conduc();//�d�C�`����
	double H=CON->get_height();			//̧�ُo�͗p�����p�����[�^
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-��(dA/dt+grad��) ---(��)

	cout<<"�ӗv�f�Q�d�����x�v�Z(jw�@)----------";
	
	complex<double> *Jec[3];
	for(int D=0;D<3;D++) Jec[D]=new complex<double> [nelm+1];
	
    for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Jec[D][i]=complex<double>(0.0,0.0);//������

	int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//�v�f���\������ӂƂ��̗v�f���ߓ_�ԍ��֌W
//	double dA[3];//�޸�����ݼ�ق̍������i�[�iX,Y,Z)
	double temp_AR[3];
	double temp_AI[3];
	double temp_VR[3];
	double temp_VI[3];

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		
		for(int D=0;D<3;D++)
		{
			Je[D][je]=0;  temp_AR[D]=0; temp_AI[D]=0; temp_VR[D]=0; temp_VI[D]=0;
		}

		
		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		{
			//�Ӂ|�ߓ_ð��ٍ쐬
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//��i�̒[�ɂ���ߓ_�́A�v�fje��j�Ԗڂ̐ߓ_
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//�̐ς�6�{
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double co=omega*sigma[je]*delta6*delta6;//j���������Ă��邪�B�e���̌v�Z�ł��܂����킹��
			double co2=sigma[je]*delta6;
			//�W���쐬
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,n�͏z���鐮��(�O�����L���v�f�@ ���E��͋Z�p�̊�b P14�Q�� ��������1.39��xm-zn��xm-xn�̊ԈႢ) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//i�������̂Ƃ�
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//i�������̂Ƃ�
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//i�������̂Ƃ�
				if(i%2!=0)//i����Ȃ�
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////

			for(int i=1;i<=6;i++)//div(��A��t)
			{	
				int iside=ELEM[je].edge[i];//�v�fje�̕Ӕԍ�
				int k1=table[i][1];//�e�_�͗v�f�̑扽�Ԑߓ_��(�O�����L���v�f�@ ���E��͋Z�p�̊�b p.42(3.11)�Ɠ����L��)
				int k2=table[i][2];

				temp_AR[A_X]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*(-AI[iside]);//�Q�d����x���������Adiv(��A)��
				temp_AI[A_X]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*AR[iside];//

				temp_AR[A_Y]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*(-AI[iside]);
				temp_AI[A_Y]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*AR[iside];

				temp_AR[A_Z]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*(-AI[iside]);
				temp_AI[A_Z]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*AR[iside];
				
			}

			for(int j=1;j<=4;j++)//grad��
			{	
				temp_VR[A_X]+=co2*c[j]*VR[N[ j]];
				temp_VI[A_X]+=co2*c[j]*VI[N[ j]];
				
				temp_VR[A_Y]+=co2*d[j]*VR[N[ j]];
				temp_VI[A_Y]+=co2*d[j]*VI[N[ j]];

				temp_VR[A_Z]+=co2*e[j]*VR[N[ j]];
				temp_VI[A_Z]+=co2*e[j]*VI[N[ j]];
			}
				
			//for(int D=0;D<3;D++) gradV[D]*=delta6;//���i���j�̑�2���v�Z����
			for(int D=0;D<3;D++)
			{
				double t_R=temp_AR[D]+temp_VR[D];
				double t_I=temp_AI[D]+temp_VI[D];
				//double t_R=-temp_AI[D]-temp_VI[D];
				//double t_I=temp_AR[D]+temp_VR[D];
				Jec[D][je]=complex<double>(-t_R,-t_I);
			}

		}
			
	}

	double Jem[3];
	double phi[3];

	for(int je=1;je<=nelm;je++)
	{
		Je_loss[je]=0.0;
		for(int D=0;D<3;D++)
		{
			Jem[D]=0.0;
			phi[D]=0.0;
		}
		
		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		//if(ELEM[je].material==FLUID)
		{
			double Jem_sum=0.0;
			for(int D=0;D<3;D++)
			{
				double Re=Jec[D][je].real();
				double Im=Jec[D][je].imag();
				Jem[D]=sqrt(Re*Re+Im*Im);
				Jem_sum+=Jem[D]*Jem[D];
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(Im,Re);
				Je[D][je]=Jem[D]*cos(omega*TIME+phi[D]);
			}
			Jem_sum=sqrt(Jem_sum);
			//�Q�d�����̌v�Z
			Je_loss[je]=Jem_sum*Jem_sum*ELEM[je].volume/(2*sigma[je]*6);
		}
		//a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		//p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	//a.close();
	//p.close();

	cout<<"ok"<<endl;

	//
	//�Q�d����v�f����ߓ_�ɕ��z����(�o�͗p)
	//int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
	//double X[5];
	//double Y[5];
	//double Z[5];
	double Je_sum=0;

	
	int *count_e=new int [node+1];	//�e���_�܂��ɂ����v�Z�Ώۂ̗v�f�����邩�B�􉽓I�Ȋ֌W�͂��łɋ��߂��Ă��邪�A���̋��E���̊O�Ɠ��̂ǂ���ŕ��ς��邩�ŕς���Ă���̂ł����ł��Ƃ߂�

	//�Q�d�����x��v�f����ߓ_�ɕ��z����
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//���ڂ��Ă���v�f�ŉQ�d�������v�Z���邩���Ȃ���
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_sum=sqrt((Je[A_X][je]*Je[A_X][je]+Je[A_Y][je]*Je[A_Y][je]+Je[A_Z][je]*Je[A_Z][je])/2.0);//�����l�ŏo��
					NODE[N[j]].Je+=Je_sum;
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			NODE[i].Je/=count_e[i];
		}
	}
	//////////
	//�o��
	cout<<"�Q�d�����x�o��------";
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
		//if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		//if(r[A_Y]<=0)
		//if(r[A_Y]>=0 && r[A_X]>=0)
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//for(int j=1;j<=4;j++)
			{
				//int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>-0.01 && r[A_Z]<0.01)
			//if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
			{
				count2++;
			}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count<<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//�n�Z�����̒��S����̒f��
		//if(r[A_Y]>=0 && r[A_X]>=0)
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
			//for(int j=1;j<=4;j++)
			{
				//int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2 <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>-0.01 && r[A_Z]<0.01)
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///�ǂݍ��ݗp�t�@�C���쐬
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//�n�Z�����̒��S����̒f��
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	cout<<"ok"<<endl;
	/////*/

	if(CON->get_model_number()==25)
	{
		//�Q�d�����_��
		//cout<<"�S�d���̗��_���v�Z"<<endl;//�����u�O�����L���v�f�@�vp.89 100*100*10�̓���
		double Ir=0; //�S�Q�d��
		double sum=0;
		for(int i=1;i<=100000;i++)//�������̖����a���v�Z
		{
			sum+=pow(-1.0,i)/(pow((2*i+1),3.0)*cosh(2*i+1)*PI*0.5);	
		}
		Ir=-0.5*CON->get_uniform_B()*(-2*PI*CON->get_Hz()*1*CON->get_ele_conduc2()*(0.1*0.5)*(0.1*0.5)*0.01*(1-32*sum/(pow(PI,3.0))));//*sin(2*PI*CON->get_Hz()*TIME)
		cout<<"�S�Q�d��(���_��)[A]="<<Ir<<endl;

		//��͉�
		Ir=0;
		double Ir2=0; //�S�Q�d��
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				if(-0.01<r[A_X]<0.01 && r[A_Y]>=0)
				{
					Ir+=abs(Je[A_X][i]);
					Ir2+=sqrt(Je[A_X][i]*Je[A_X][i]+Je[A_Y][i]*Je[A_Y][i]+Je[A_Z][i]*Je[A_Z][i]);
				}
			}
		}
		Ir*=0.05*0.01/2;
		Ir2*=0.05*0.01/2;
		cout<<"�S�Q�d��(��͉�,X�����̂�)[A]="<<Ir<<"�@(��͉�,�S����)[A]="<<Ir2<<endl;
	}


	
	//fout.close();
	for(int D=0;D<DIMENTION;D++) delete [] Jec[D];
	delete [] count_e;

}

//���b�V���X���[�W���O
void laplacian_smoothing(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int node0,int nelm0)
{
	cout<<"���v���V�A���X���[�W���O "<<*node<<endl;
	int *flag=new int[*node+1];
	for(int i=0;i<=*node;i++) flag[i]=OFF;
	int *flag2=new int[*nelm+1];
	for(int j=0;j<=*nelm;j++) flag2[j]=OFF;
	int count=0;

	for(int i=node0;i<=*node;i++)
	{
		
		if(NODE[i].material==AIR && NODE[i].BD_node==OFF)//���̐ߓ_�ł������b�V�����E�ߓ_�ł��Ȃ��ߓ_�����X���[�W���O
		{
			for(int D=0;D<3;D++) NODE[i].r[D]=0;
			for(int j=1;j<=*nelm;j++)
			{
				for(int k=1;k<=4;k++) if(ELEM[j].node[k]==i) flag2[j]=ON; //���ڂ��Ă���ߓ_��p���č\�������v�f�Ȃ�ON

				if(flag2[j]==ON) for(int k=1;k<=4;k++) flag[ELEM[j].node[k]]=ON;
				flag2[j]=OFF;
			}
			flag[i]=OFF; //���ڂ��Ă���_��flag��ON�ɂȂ��Ă��邪�A�~�����̂͂܂��̓_�Ȃ̂�flag��OFF
				
			for(int l=1;l<=*node;l++)
			{
				if(flag[l]==ON)
				{
					for(int D=0;D<3;D++) NODE[i].r[D]+=NODE[l].r[D];
					count++;
				}
			}
				
			for(int D=0;D<3;D++) NODE[i].r[D]/=count;
			count=0;
			for(int l=1;l<=*node;l++) flag[l]=OFF;
		}
	}
    delete [] flag;
	delete [] flag2;
}

//���b�V���X���[�W���O
void laplacian_smoothing2(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int *jnb,int **nei,int node0)
{
	//cout<<"���v���V�A���X���[�W���O"<<endl;
	int *flag=new int[*node+1];
	for(int i=0;i<=*node;i++) flag[i]=OFF;
	int *flag2=new int[*nelm+1];
	for(int j=0;j<=*nelm;j++) flag2[j]=OFF;
	int count=0;
	double *r[3];
	for(int D=0;D<3;D++) r[D]=new double [*node+1];

	for(int i=node0;i<=*node;i++)
	{
		for(int D=0;D<3;D++) r[D][i]=NODE[i].r[D];
		if(NODE[i].material==AIR)//���̐ߓ_�ł������b�V�����E�ߓ_�ł��Ȃ��ߓ_�����X���[�W���O�Bnode0�ȍ~�̓����b�V���̈�ɓ��������ߓ_�Q�̂͂�
		{
			for(int D=0;D<3;D++) r[D][i]=0;
			for(int j=1;j<=jnb[i];j++)
			{
				int jelm=nei[i][j];//�ߓ_i�ɗאڂ���v�f�ԍ�
				for(int k=1;k<=4;k++) if(ELEM[jelm].node[k]==i) flag2[jelm]=ON; //���ڂ��Ă���ߓ_��p���č\�������v�f�Ȃ�ON

				if(flag2[jelm]==ON) for(int k=1;k<=4;k++) flag[ELEM[jelm].node[k]]=ON;
				flag2[jelm]=OFF;
			}
			flag[i]=OFF; //���ڂ��Ă���_��flag��ON�ɂȂ��Ă��邪�A�~�����̂͂܂��̓_�Ȃ̂�flag��OFF
				
			for(int l=1;l<=*node;l++)
			{
				if(flag[l]==ON)
				{
					for(int D=0;D<3;D++) r[D][i]+=NODE[l].r[D];
					count++;
					
				}
			}
				
			for(int D=0;D<3;D++) r[D][i]/=count;
			count=0;
			for(int l=1;l<=*node;l++) flag[l]=OFF;
		}
	}
	for(int i=node0;i<=*node;i++) for(int D=0;D<3;D++) NODE[i].r[D]=r[D][i];
    delete [] flag;
	delete [] flag2;
	for(int D=0;D<3;D++) delete [] r[D];
}

//�ߓ_�ԍ��̕��ёւ��֐�
void node_sorting(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei)
{
	unsigned int timeA=GetTickCount();
	cout<<"�ߓ_�ԍ����ёւ�"<<endl;

	////�ߓ_-�ߓ_�֌W
	vector<vector<int>> NEI2;//�ߓ_-�ߓ_�֌W NEI2[i][j]�@�ߓ_�ԍ�i�ɁAj�ԖڂɌ������Ă���ߓ_�ԍ����i�[�����
	NEI2.resize(node+1);
	int *temp=new int[node+1];//�`�F�b�N�p�z��

	for(int i=1;i<=node;i++)//���ڂ��Ă���ߓ_�Ɍq�����Ă���ߓ_�𒲂ׂ�
	{
		//cout<<i<<endl;
		for(int k=1;k<=jnb[i];k++)
		{
			
			int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp[p]==0)
				{
					NEI2[i].push_back(p);
					temp[p]=1;//��������
				}
			}
		}
		//for(int k=0;k<=node;k++) temp[k]=0;//������
		for(int k=0;k<NEI2[i].size();k++) temp[NEI2[i][k]]=0;//������
	}
	cout<<"�ߓ_-�ߓ_�֌W�쐬"<<endl;

	int nei2max=0;
	for(int i=1;i<=node;i++) if(NEI2[i].size()>nei2max) nei2max=(int) NEI2[i].size();
	cout<<"�ő匋����="<<nei2max;

	/*//
	int nei2min=node;
	for(int i=1;i<=node;i++) if(NEI2[i].size()<nei2min) nei2min=NEI2[i].size();
	int L1ID=0;
	int L1flag=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NEI2[i].size()==nei2min)
		{
			if(L1flag==OFF)
			{
				L1ID=i;
				L1flag=ON;
			}
		}
	}
	//*/
	
	////���x���ݒ�
	int *L=new int[node+1];//�ߓ_�̃��x��
	for(int i=0;i<=node;i++) L[i]=0;	

	int count=node;
	for(int i=1;i<=node;i++) if(jnb[i]==0) count--;//jnb=0�̐ߓ_�̓��x����������0�Ȃ̂ŁA���x���Z�b�g����ߓ_�����珜�O����

	L[1]=1; //�ߓ_�ԍ�1�̐ߓ_�����x��1�Ƃ���
	//L[L1ID]=1;
	count--;

	int level=1;
	while(count!=0)
	{
		for(int i=1;i<=node;i++)
		{
			if(L[i]==level)
			{
				for(int k=0;k<NEI2[i].size();k++)
				{
					if(L[NEI2[i][k]]==0)//�ׂ̐ߓ_�̃��x�����܂����܂��ĂȂ�
					{
						L[NEI2[i][k]]=L[i]+1;
						count--;
					}
				}
			}
		}
		level++;
	}

	//���x������̐ߓ_�Ƃ̌����������߂�i�\�[�g�p�j
	int *nei_U=new int[node+1];//��������̃��x���̐ߓ_�Ƃ̌�����
	for(int i=0;i<=node;i++) nei_U[i]=0;
	for(int i=1;i<=node;i++)
	{
		for(int k=0;k<NEI2[i].size();k++)
		{
			if(L[i]<L[NEI2[i][k]]) nei_U[i]+=1;
		}
	}
	int nei_Umax=0;
	for(int i=1;i<=node;i++) if(nei_Umax<nei_U[i]) nei_Umax=nei_U[i];

	int Lmax=0;
	for(int i=1;i<=node;i++) if(Lmax<L[i]) Lmax=L[i];

	//���x�����Ƃ̐ߓ_�������߂�
	int *N=new int[Lmax+1];//���x�����Ƃ̐ߓ_��
	for(int i=0;i<=Lmax;i++) N[i]=0;
	for(int i=0;i<=Lmax;i++) for(int j=1;j<=node;j++) if(L[j]==i) N[i]+=1;

	int Nmax=0;
	for(int i=1;i<=Lmax;i++) if(Nmax<N[i]) Nmax=N[i];
	cout<<"Lmax="<<Lmax<<" Nmax="<<Nmax<<endl;

	vector<vector<int>> Ln;//Ln[i][j] i:���x���@j:���x��i�ɏ�������ߓ_�ԍ������������̂���i�[
	Ln.resize(Lmax+1);
	for(int i=0;i<(int)Ln.size();i++) Ln[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++)
	{
		int k=0;
		for(int j=1;j<=node;j++)
		{
			if(L[j]==i)
			{
				k++;
				Ln[i][k]=j;		
			}
		}
	}

	/*//�`�F�b�N
	for(int i=1;i<=node;i++)
	{
		int countL=0;
		int countR=0;
		if(L[i]==0 && jnb[i]>0) cout<<"���x���Z�b�g����Ă��Ȃ��ߓ_���� i="<<i<<endl;
		for(int k=0;k<NEI2[i].size();k++)
		{
			//if(L[NEI2[i][k]]==L[i]) cout<<"�������x���̐ߓ_���א� i="<<i<<endl;
			if(L[NEI2[i][k]]==L[i]) countL++;
			if(L[NEI2[i][k]]>=L[i]) countR++;
			if(L[NEI2[i][k]]>L[i]+1 ||L[NEI2[i][k]]<L[i]-1) cout<<"2�ȏ㗣�ꂽ���x���̐ߓ_���א� i="<<i<<endl;
		}
		if(countL==NEI2[i].size() && jnb[i]>0) cout<<"�������x���̐ߓ_�̂ݗא� i="<<i<<" L[i]="<<L[i]<<endl;
		if(countR==NEI2[i].size() && jnb[i]>0 && L[i]>1) cout<<"���ʃ��x���̐ߓ_�ɗאڂ��Ă��Ȃ� i="<<i<<" L[i]="<<L[i]<<endl;

	}
	//*/

	//���ёւ��J�n
	cout<<"���ёւ��J�n"<<endl;
	int *flag=new int[node+1];
	int *flag2=new int[node+1];
	vector<vector<int>> Ln2;//�V�����ߓ_�ԍ��Ɋ�Â��č����Ln
	Ln2.resize(Lmax+1);
	for(int i=0;i<(int)Ln2.size();i++) Ln2[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++) for(int j=0;j<=N[i];j++) Ln2[i][j]=0;
	
	Ln2[1][1]=1;
	//Ln2[1][1]=L1ID;

	int sum=1;
	for(int w=0;w<=node;w++) flag2[w]=OFF;
	//int i=2;
	for(int i=2;i<=Lmax;i++)
	{
		sum+=N[i];
		//cout<<"���x��"<<i<<"�ɑ�����ߓ_��"<<N[i]<<endl;
		int endj=1;
		int k=1;
		int swap=0;
		
		while(endj<=N[i])
		{
			for(int iflag=0;iflag<=node;iflag++) flag[iflag]=OFF;
	
			int num=0;
			for(int m=1;m<=N[i-1];m++) if(Ln[i-1][m]==Ln2[i-1][k]) swap=m;

			//���x��i-1��k�Ԗڂ̐ߓ_�Ɍ������Ă���ߓ_��I�яo�� num:�I�񂾐ߓ_�̐�
			for(int j=1;j<=N[i];j++)
			{			
				if(flag2[Ln[i][j]]==OFF)
				{
					for(int knei=0;knei<NEI2[Ln[i][j]].size();knei++)
					{
						//if(NEI2[Ln[i][j]][knei]==Ln[i-1][k])
						if(NEI2[Ln[i][j]][knei]==Ln[i-1][swap])//���ёւ�����̂��̂�k�������Ĕ�r����
						{
							flag[Ln[i][j]]=ON;//����̃��[�v�ŕ��ёւ���ߓ_�͂��̃t���O��ON�ɂȂ�
							num++;
							flag2[Ln[i][j]]=ON;//���̃t���O��ON�ɂȂ����ߓ_�͎��̃��[�v�ŕ��ёւ��̑ΏۂƂ��Ȃ�
						}
					}
				}
			}
			int end_old=endj;
			endj+=num;
			//cout<<"num="<<num<<endl;			
			//cout<<"endj_old="<<end_old<<endl;
			//cout<<"N[i-1]="<<N[i-1]<<" k="<<k<<" endj="<<endj<<" N[i]="<<N[i]<<endl;

			//�I�񂾐ߓ_��Ln[i][j]�`Ln[i][endj+num-1]�Ɋ��蓖�Ă�
			
			int je=end_old;

			/*//
			for(int j=1;j<=node;j++)
			{
				if(flag[j]==ON)
				{
					Ln2[j][je]=n;
				}
			}
			//*/
			
			//���蓖�Ă��ߓ_���A�����������ƂɃ\�[�g����i�����������������̂���ԍ��t��) 
			//�������������ꍇ�A�ߓ_�ԍ����������ق����ɂ���
			//���Ȃ��ق��������H
			///
			while(je<endj)
			{
				int nei_Umin=nei_Umax;
				for(int n=1;n<=node;n++) if(flag[n]==ON) if(nei_U[n]<nei_Umin) nei_Umin=nei_U[n];
				int flag_Umin=OFF;
				for(int n=1;n<=node;n++)
				{
					if(nei_U[n]==nei_Umin)
					{
						if(flag_Umin==OFF)//�ŏ��̌����������ߓ_����������Ƃ��͑������̏����ɂ���
						{
							if(flag[n]==ON)
							{
								Ln2[i][je]=n;
								flag_Umin=ON;
								flag[n]=OFF;
							}
						}
					}
				}
				je++;
			}
			//*/
			k++;
		}
		//cout<<"endj="<<endj<<" N[i]="<<N[i]<<endl;
	}

	/*//
	cout<<"Ln"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln[2][i]<<" "<<nei_U[Ln[2][i]]<<endl;
	}
	cout<<"Ln2"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln2[2][i]<<" "<<nei_U[Ln2[2][i]]<<endl;
	}
	//*/

	cout<<"�V�ߓ_�ԍ��̎Z�o����"<<endl;

	//�V�����ԍ������Ă���
	//���x���̒Ⴂ���́A�\�[�g�������x�����Ƃ̔z��̓��̂ق����珇�ɂ���
	int *newID=new int[node+1];//�V�����ߓ_�ԍ��@���ёւ��O�̐ߓ_�ԍ���i
	for(int i=0;i<=node;i++) newID[i]=0;
	int countID=0;
	for(int i=1;i<=Lmax;i++)
	{
		for(int j=1;j<=N[i];j++)
		{
			countID++;
			newID[Ln2[i][j]]=countID;
		}
	}
	//jnb=0�̐ߓ_�͖����ɂ���
	for(int j=1;j<=N[0];j++)
	{
		countID++;
		newID[Ln[0][j]]=countID;
	}

	//�ߓ_�Ɋ֌W����e���̏C��
	cout<<"�ߓ_�ԍ��ύX�ɔ������̏C���J�n"<<endl;
	modify_node_info(CON, node, nelm,NODE,ELEM,jnb,nei,newID);

	delete [] L;
	delete [] flag;
	delete [] flag2;
	delete [] temp;
	delete [] N;
	delete []newID;

	cout<<"���ёւ������|�|time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//�ߓ_�ԍ��̕��ёւ��֐� ���ǔł������ɁH
void node_sorting2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei)
{
	unsigned int timeA=GetTickCount();
	cout<<"�ߓ_�ԍ����ёւ�"<<endl;

	////�ߓ_-�ߓ_�֌W
	vector<vector<int>> NEI2;//�ߓ_-�ߓ_�֌W NEI2[i][j]�@�ߓ_�ԍ�i�ɁAj�ԖڂɌ������Ă���ߓ_�ԍ����i�[�����
	NEI2.resize(node+1);
	int *temp=new int[node+1];//�`�F�b�N�p�z��

	for(int i=1;i<=node;i++)//���ڂ��Ă���ߓ_�Ɍq�����Ă���ߓ_�𒲂ׂ�
	{
		//cout<<i<<endl;
		for(int k=1;k<=jnb[i];k++)
		{
			
			int jelm=nei[i][k];//�ߓ_i�ɗאڂ���v�f�ԍ�
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp[p]==0)
				{
					NEI2[i].push_back(p);
					temp[p]=1;//��������
				}
			}
		}
		//for(int k=0;k<=node;k++) temp[k]=0;//������
		for(int k=0;k<NEI2[i].size();k++) temp[NEI2[i][k]]=0;//������
	}
	cout<<"�ߓ_-�ߓ_�֌W�쐬"<<endl;

	

	/*//
	int nei2min=node;
	for(int i=1;i<=node;i++) if(NEI2[i].size()<nei2min) nei2min=NEI2[i].size();
	int L1ID=0;
	int L1flag=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NEI2[i].size()==nei2min)
		{
			if(L1flag==OFF)
			{
				L1ID=i;
				L1flag=ON;
			}
		}
	}
	//*/
	
	////���x���ݒ�
	int *L=new int[node+1];//�ߓ_�̃��x��
	for(int i=0;i<=node;i++) L[i]=0;	

	int count=node;
	for(int i=1;i<=node;i++) if(jnb[i]==0) count--;//jnb=0�̐ߓ_�̓��x����������0�Ȃ̂ŁA���x���Z�b�g����ߓ_�����珜�O����

	L[1]=1; //�ߓ_�ԍ�1�̐ߓ_�����x��1�Ƃ���
	//L[L1ID]=1;
	count--;

	int level=1;
	while(count!=0)
	{
		for(int i=1;i<=node;i++)
		{
			if(L[i]==level)
			{
				for(int k=0;k<NEI2[i].size();k++)
				{
					if(L[NEI2[i][k]]==0)//�ׂ̐ߓ_�̃��x�����܂����܂��ĂȂ�
					{
						L[NEI2[i][k]]=L[i]+1;
						count--;
					}
				}
			}
		}
		level++;
	}

	//���x������̐ߓ_�Ƃ̌����������߂�i�\�[�g�p�j
	int *nei_U=new int[node+1];//��������̃��x���̐ߓ_�Ƃ̌�����
	for(int i=0;i<=node;i++) nei_U[i]=0;
	for(int i=1;i<=node;i++)
	{
		for(int k=0;k<NEI2[i].size();k++)
		{
			if(L[i]<L[NEI2[i][k]]) nei_U[i]+=1;
		}
	}
	int nei_Umax=0;
	for(int i=1;i<=node;i++) if(nei_Umax<nei_U[i]) nei_Umax=nei_U[i];

	int Lmax=0;
	for(int i=1;i<=node;i++) if(Lmax<L[i]) Lmax=L[i];

	//���x�����Ƃ̐ߓ_�������߂�
	int *N=new int[Lmax+1];//���x�����Ƃ̐ߓ_��
	for(int i=0;i<=Lmax;i++) N[i]=0;
	for(int i=0;i<=Lmax;i++) for(int j=1;j<=node;j++) if(L[j]==i) N[i]+=1;

	int Nmax=0;
	for(int i=1;i<=Lmax;i++) if(Nmax<N[i]) Nmax=N[i];
	cout<<"Lmax="<<Lmax<<" Nmax="<<Nmax<<endl;

	vector<vector<int>> Ln;//Ln[i][j] i:���x���@j:���x��i�ɏ�������ߓ_�ԍ������������̂���i�[
	Ln.resize(Lmax+1);
	for(int i=0;i<(int)Ln.size();i++) Ln[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++)
	{
		int k=0;
		for(int j=1;j<=node;j++)
		{
			if(L[j]==i)
			{
				k++;
				Ln[i][k]=j;		
			}
		}
	}

	/*//�`�F�b�N
	for(int i=1;i<=node;i++)
	{
		int countL=0;
		int countR=0;
		if(L[i]==0 && jnb[i]>0) cout<<"���x���Z�b�g����Ă��Ȃ��ߓ_���� i="<<i<<endl;
		for(int k=0;k<NEI2[i].size();k++)
		{
			//if(L[NEI2[i][k]]==L[i]) cout<<"�������x���̐ߓ_���א� i="<<i<<endl;
			if(L[NEI2[i][k]]==L[i]) countL++;
			if(L[NEI2[i][k]]>=L[i]) countR++;
			if(L[NEI2[i][k]]>L[i]+1 ||L[NEI2[i][k]]<L[i]-1) cout<<"2�ȏ㗣�ꂽ���x���̐ߓ_���א� i="<<i<<endl;
		}
		if(countL==NEI2[i].size() && jnb[i]>0) cout<<"�������x���̐ߓ_�̂ݗא� i="<<i<<" L[i]="<<L[i]<<endl;
		if(countR==NEI2[i].size() && jnb[i]>0 && L[i]>1) cout<<"���ʃ��x���̐ߓ_�ɗאڂ��Ă��Ȃ� i="<<i<<" L[i]="<<L[i]<<endl;

	}
	//*/

	//���ёւ��J�n

	cout<<"���ёւ��J�n"<<endl;
	int *flag=new int[node+1];
	int *flag2=new int[node+1];
	vector<vector<int>> Ln2;//�V�����ߓ_�ԍ��Ɋ�Â��č����Ln
	Ln2.resize(Lmax+1);
	for(int i=0;i<(int)Ln2.size();i++) Ln2[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++) for(int j=0;j<=N[i];j++) Ln2[i][j]=0;
	
	Ln2[1][1]=1;
	//Ln2[1][1]=L1ID;

	int sum=1;
	for(int w=0;w<=node;w++) flag2[w]=OFF;
	//int i=2;
	for(int i=2;i<=Lmax;i++)
	{
		sum+=N[i];
		cout<<"���x��"<<i<<"�ɑ�����ߓ_��"<<N[i]<<endl;
		int endj=1;
		int k=1;
		int swap=0;
		
		while(endj<=N[i])
		{
			for(int iflag=0;iflag<=node;iflag++) flag[iflag]=OFF;
	
			int num=0;
			for(int m=1;m<=N[i-1];m++) if(Ln[i-1][m]==Ln2[i-1][k]) swap=m;

			//���x��i-1��k�Ԗڂ̐ߓ_�Ɍ������Ă���ߓ_��I�яo�� num:�I�񂾐ߓ_�̐�
			for(int j=1;j<=N[i];j++)
			{			
				if(flag2[Ln[i][j]]==OFF)
				{
					for(int knei=0;knei<NEI2[Ln[i][j]].size();knei++)
					{
						//if(NEI2[Ln[i][j]][knei]==Ln[i-1][k])
						if(NEI2[Ln[i][j]][knei]==Ln[i-1][swap])//���ёւ�����̂��̂�k�������Ĕ�r����
						{
							flag[Ln[i][j]]=ON;//����̃��[�v�ŕ��ёւ���ߓ_�͂��̃t���O��ON�ɂȂ�
							num++;
							flag2[Ln[i][j]]=ON;//���̃t���O��ON�ɂȂ����ߓ_�͎��̃��[�v�ŕ��ёւ��̑ΏۂƂ��Ȃ�
						}
					}
				}
			}
			int end_old=endj;
			endj+=num;
			//cout<<"num="<<num<<endl;			
			//cout<<"endj_old="<<end_old<<endl;
			//cout<<"N[i-1]="<<N[i-1]<<" k="<<k<<" endj="<<endj<<" N[i]="<<N[i]<<endl;

			//�I�񂾐ߓ_��Ln[i][j]�`Ln[i][endj+num-1]�Ɋ��蓖�Ă�

			int je=end_old;
			for(int n=1;n<=node;n++)
			{
				if(flag[n]==ON)
				{
					Ln2[i][je]=n;
					je++;
				}
			}
			k++;
		}
		//cout<<"endj="<<endj<<" N[i]="<<N[i]<<endl;
	}

	/*//
	cout<<"Ln"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln[2][i]<<" "<<nei_U[Ln[2][i]]<<endl;
	}
	cout<<"Ln2"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln2[2][i]<<" "<<nei_U[Ln2[2][i]]<<endl;
	}
	//*/

	cout<<"�V�ߓ_�ԍ��̎Z�o����"<<endl;

	//�V�����ԍ������Ă���
	//���x���̒Ⴂ���́A�\�[�g�������x�����Ƃ̔z��̓��̂ق����珇�ɂ���
	int *newID=new int[node+1];//�V�����ߓ_�ԍ��@���ёւ��O�̐ߓ_�ԍ���i
	for(int i=0;i<=node;i++) newID[i]=0;
	int countID=0;
	for(int i=1;i<=Lmax;i++)
	{
		for(int j=1;j<=N[i];j++)
		{
			countID++;
			newID[Ln2[i][j]]=countID;
		}
	}
	//jnb=0�̐ߓ_�͖����ɂ���
	for(int j=1;j<=N[0];j++)
	{
		countID++;
		newID[Ln[0][j]]=countID;
	}

	//�ߓ_�Ɋ֌W����e���̏C��
	cout<<"�ߓ_�ԍ��ύX�ɔ������̏C���J�n"<<endl;
	modify_node_info(CON, node, nelm,NODE,ELEM,jnb,nei,newID);

	delete [] L;
	delete [] flag;
	delete [] flag2;
	delete [] temp;
	delete [] N;
	delete []newID;

	cout<<"���ёւ������|�|time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


//���C���֐�
void modify_node_info(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei,int *newID)
{
	vector<point3D> N_NODE;
	N_NODE.resize(node+1);

	//NODE�N���X�̓��e���ꎞ�z��ɋL��
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) N_NODE[newID[i]].r[D]=NODE[i].r[D];
		N_NODE[newID[i]].material=NODE[i].material;
		N_NODE[newID[i]].boundary_condition=NODE[i].boundary_condition;
		N_NODE[newID[i]].particleID=NODE[i].particleID;
		N_NODE[newID[i]].remesh=NODE[i].remesh;
		N_NODE[newID[i]].BD_node=NODE[i].BD_node;
	}

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) NODE[i].r[D]=N_NODE[i].r[D];
		NODE[i].material=N_NODE[i].material;
		NODE[i].boundary_condition=N_NODE[i].boundary_condition;
		NODE[i].particleID=N_NODE[i].particleID;
		NODE[i].remesh=N_NODE[i].remesh;
		NODE[i].BD_node=N_NODE[i].BD_node;
	}

	//�v�f�͍\������ߓ_�ԍ��̂ݏC��
	for(int i=1;i<=nelm;i++)
	{
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];

		ELEM[i].node[1]=newID[ia];
		ELEM[i].node[2]=newID[ib];
		ELEM[i].node[3]=newID[ic];
		ELEM[i].node[4]=newID[ip];
	}

	fill3D(NODE,ELEM,nelm);
	
	//jnb,nei�͕��ёւ����I��������Ƌ��߂Ȃ����Ă���

}

//�͂�e�����ʂ��R���^�[�}�ŕ\������֐�
void output_F_scalar_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,int node)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	//t=1;//���܂͂킴�Ɩ��X�e�b�v�㏑��
	//int number=int (PART.size());
	//int node=int(NODE.size()-1);
	//int nelm=int(ELEM.size());

	//sprintf_s(filename,"pressure/pressure%d",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����

	cout<<"���q�̕����ʏo�͊J�n"<<endl;

	////�������x
	//sprintf_s(filename,"PART.B%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/B%d",t);
	ofstream fout(filename);
	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID)
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].B;
					fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout.close();


	//sprintf_s(filename,"PART.B%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/B%d.fld",t);
	ofstream fout2(filename);
	if(!fout2)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	fout2 << "dim1=" << n <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=1" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=Bflux" << endl << endl;
	fout2 << "variable 1 file=B" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    1 file=B" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    2 file=B" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    3 file=B" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2.close();

	////�������x�f��
	n=0;
	//sprintf_s(filename,"PART.B%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Bsec%d",t);
	ofstream fouts(filename);
	if(!fouts)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].B;
					fouts << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fouts.close();


	//sprintf_s(filename,"PART.B%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Bsec%d.fld",t);
	ofstream fout2s(filename);
	if(!fout2s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout2s << "# AVS field file" << endl;
	fout2s << "ndim=1" << endl;
	fout2s << "dim1=" << n <<endl;
	fout2s << "nspace=3" << endl;
	fout2s << "veclen=1" << endl;
	fout2s << "data=float" << endl;
	fout2s << "field=irregular" << endl;
	fout2s << "label=Bflux" << endl << endl;
	fout2s << "variable 1 file=Bsec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2s << "coord    1 file=Bsec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2s << "coord    2 file=Bsec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2s << "coord    3 file=Bsec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2s.close();

	///�Q�d�����x
	n=0;
	//sprintf_s(filename,"PART.J%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Jesec%d",t);
	ofstream fout3(filename);
	if(!fout3)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID)
				
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].Je;
					fout3 << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout3.close();

	//sprintf_s(filename,"PART.J%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Jesec%d.fld",t);
	ofstream fout4(filename);
	if(!fout4)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout4 << "# AVS field file" << endl;
	fout4 << "ndim=1" << endl;
	fout4 << "dim1=" << n <<endl;
	fout4 << "nspace=3" << endl;
	fout4 << "veclen=1" << endl;
	fout4 << "data=float" << endl;
	fout4 << "field=irregular" << endl;
	fout4 << "label=eddy " << endl << endl;
	fout4 << "variable 1 file=Je" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4 << "coord    1 file=Je" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4 << "coord    2 file=Je" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4 << "coord    3 file=Je" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4.close();

	///�Q�d�����x �f��
	n=0;
	//sprintf_s(filename,"PART.J%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Jesec%d",t);
	ofstream fout3s(filename);
	if(!fout3s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				//if(PART[p].type==FLUID)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].Je;
					fout3s << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout3s.close();

	//sprintf_s(filename,"PART.J%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Jesec%d.fld",t);
	ofstream fout4s(filename);
	if(!fout4s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout4s << "# AVS field file" << endl;
	fout4s << "ndim=1" << endl;
	fout4s << "dim1=" << n <<endl;
	fout4s << "nspace=3" << endl;
	fout4s << "veclen=1" << endl;
	fout4s << "data=float" << endl;
	fout4s << "field=irregular" << endl;
	fout4s << "label=eddy " << endl << endl;
	fout4s << "variable 1 file=Jesec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4s << "coord    1 file=Jesec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4s << "coord    2 file=Jesec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4s << "coord    3 file=Jesec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout4s.close();

	///�d����
	n=0;
	//sprintf_s(filename,"PART.F%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/F%d",t);
	ofstream fout5(filename);
	if(!fout5)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				if(PART[p].type==FLUID)
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=sqrt(PART[p].F[A_X]*PART[p].F[A_X]+PART[p].F[A_Y]*PART[p].F[A_Y]+PART[p].F[A_Z]*PART[p].F[A_Z]);
					fout5 << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout5.close();

	//sprintf_s(filename,"PART.F%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/F%d.fld",t);
	ofstream fout6(filename);
	if(!fout6)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout6 << "# AVS field file" << endl;
	fout6 << "ndim=1" << endl;
	fout6 << "dim1=" << n <<endl;
	fout6 << "nspace=3" << endl;
	fout6 << "veclen=1" << endl;
	fout6 << "data=float" << endl;
	fout6 << "field=irregular" << endl;
	fout6 << "label=lorentz" << endl << endl;
	fout6 << "variable 1 file=F" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6 << "coord    1 file=F" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6 << "coord    2 file=F" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6 << "coord    3 file=F" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6.close();

	///�d���͒f��
	n=0;
	//sprintf_s(filename,"PART.F%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Fsec%d",t);
	ofstream fout5s(filename);
	if(!fout5s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=sqrt(PART[p].F[A_X]*PART[p].F[A_X]+PART[p].F[A_Y]*PART[p].F[A_Y]+PART[p].F[A_Z]*PART[p].F[A_Z]);
					fout5s << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout5s.close();

	//sprintf_s(filename,"PART.F%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	sprintf_s(filename,"plot_scalar/Fsec%d.fld",t);
	ofstream fout6s(filename);
	if(!fout6s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout6s << "# AVS field file" << endl;
	fout6s << "ndim=1" << endl;
	fout6s << "dim1=" << n <<endl;
	fout6s << "nspace=3" << endl;
	fout6s << "veclen=1" << endl;
	fout6s << "data=float" << endl;
	fout6s << "field=irregular" << endl;
	fout6s << "label=lorentz" << endl << endl;
	fout6s << "variable 1 file=Fsec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6s << "coord    1 file=Fsec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6s << "coord    2 file=Fsec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6s << "coord    3 file=Fsec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout6s.close();

	///
	/*
	//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
	fp<<"6 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	//fp<<"2 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"6 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"speed,m/s"<<endl;
	fp<<"surface_tension,N/m^3"<<endl;
	fp<<"B,T"<<endl;
	fp<<"Je,A/m^2"<<endl;
	//fp<<"Et,V/m"<<endl;
	fp<<"Fn,N/m^3"<<endl;
	//fp<<"P,N/m^2"<<endl;
	fp<<"value1,??"<<endl;

	

	//�e�ߓ_�̏��l����
	for(int i=0;i<node;i++)
	{
		int p=NODE[i].particleID;
		double speed=0;
		double F=0;
		double vn=0;
		double potential=0;
		double Pst=0;//
		
		//���x�v�Z
		if(p>=0)
		{
			for(int D=0;D<3;D++)
			{
				speed+=PART[p].u[D]*PART[p].u[D];
				F+=PART[p].F[D]*PART[p].F[D];
				potential+=PART[p].potential[D]*PART[p].potential[D];
			}
			Pst=PART[p].dir_Pst;
		}
		speed=sqrt(speed);
		F=sqrt(F)*CON->get_density()/CON->get_particle_mass();
		potential=sqrt(F)*CON->get_density();
		////

		//double P=Pst-NODE[i].Fn;
		if(p>=0) fp<<count<<" "<<speed<<" "<<potential<<" "<<NODE[i].B<<" "<<NODE[i].Je<<" "<<NODE[i].F<<" "<<NODE[i].value1<<endl;
	}
	*/
	
	/*fp<<"3 1 1 1"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	
	//�v�f���o�́@�v�f����microAVS�̉������\�b�h�o�[���́A�u�v�f�f�[�^�̓h��Ԃ��v�������Ό����
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	*/
	//cout<<"OK"<<endl;
	//fp.close();
}

//�͂Ȃǂ̕����ʂ��R���^�[�}�ŕ\������֐�
void output_F_scalar_movie_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,double TIME)
{
	/*
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v
	cout<<"Now F_scalar movie writing-----";
	int nelm=int (ELEM.size());		//�v�f��
	int node=int (NODE.size());		//�ߓ_��
	int STEP=CON->get_step()/(CON->get_EM_interval()*CON->get_mesh_output_interval())+1;				//�o�͂��鑍�X�e�b�v��
	
	if(t==1) 
	{
		ofstream fp("F_scalar_movie.inp");
		fp<<STEP<<endl;
		fp<<"data_geom"<<endl;
		fp.close();

		ofstream fq("step_for_F_scalar_movie.dat");			//���̊֐����Ăяo���ꂽ�񐔂��t�@�C���ɋL�����A�Ȍ�A�֐����Ăяo����邽�тɃJ�E���g�𑝂₵�Ă����B���̒l��step�Ŏg�p����B
		fq<<1<<endl;
		fq.close();
	}

	//step_now�̒l���t�@�C����茈��
	int step_now;											//���̊֐����Ăяo���ꂽ��
	ifstream fin("step_for_F_scalar_movie.dat");
	if(!fin) cout<<"cannot open step_for_F_scalar_movie.dat"<<endl;
	fin>>step_now;
	fin.close();

	//main�t�@�C����������
	ofstream fp("F_scalar_movie.inp",ios :: app);
	fp<<"step"<<step_now<<" TIME="<<TIME<<endl;
	fp<<node<<" "<<nelm<<endl;	//�ߓ_���Ɨv�f���o��
	
	//�ߓ_�ԍ��Ƃ��̍��W�̏o�� 
	for(int i=0;i<node;i++) fp<<i+1<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	
	//�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=0;i<nelm;i++)
	{
		fp<<i+1<<"  0 tri ";
		for(int j=0;j<3;j++)	fp<<ELEM[i].node[j]+1<<" ";
		fp<<endl;
	}

	//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
	fp<<"8 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	//fp<<"2 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	fp<<"speed,m/s"<<endl;
	fp<<"surface_tension,N/m^2"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Et,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	fp<<"P,N/m^2"<<endl;
	fp<<"value1,??"<<endl;

	

	//�e�ߓ_�̏��l����
	for(int i=0;i<node;i++)
	{
		int p=NODE[i].particle;
		double speed=0;
		double vn=0;
		double Pst=0;//
		if(p>=0)
		{
			for(int D=0;D<3;D++)
			{
				speed+=PART[p].u[D]*PART[p].u[D];
			}
			Pst=PART[p].dir_Pst;
		}
		speed=sqrt(speed);
		double P=Pst-NODE[i].Fn;
		fp<<i+1<<" "<<speed<<" "<<Pst<<" "<<NODE[i].potential<<" "<<NODE[i].slop1<<" "<<NODE[i].Et<<" "<<NODE[i].Fn<<" "<<P<<" "<<NODE[i].value1<<endl;
	}

	/*fp<<"3 1 1 1"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	
	//�v�f���o�́@�v�f����microAVS�̉������\�b�h�o�[���́A�u�v�f�f�[�^�̓h��Ԃ��v�������Ό����
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	/
	fp.close();

	step_now++;											//�Ăяo���񐔂�����Ɍ����ā{�{���Ă���
	ofstream fq("step_for_F_scalar_movie.dat");			//���̒l���t�@�C���o��
	fq<<step_now<<endl;
	fq.close();

	cout<<"OK"<<endl;
	*/
	

	
}