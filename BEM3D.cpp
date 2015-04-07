#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"	//#define �i�[
//#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
//#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
//#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
#include<omp.h>
#include<vector>
#include"function.h"

void set_BEM3D_static_model(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *s_node_num,int *s_elem_num);
void set_gauss_parameter3D(int *NGS,double ep1[14][3],double *w1,double *ep_ln,double *w_ln);
void set_GHmatrix3D_for_constant(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,double ep[3],double **H,double **G,double w,int type,int n,int m,int n1,int m1);
double calc_Gnn_for_CONSTANT(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,int n1);


void BEM3D(mpsconfig *CON,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number,double **F,int t)
{
	// t=1�̂Ƃ��ɁA�Œ�v�f���쐬����s_NODE,s_ELEM�Ɋi�[����B�ߓ_���Ɨv�f����s_node_num��s_elem_num�Ɋi�[�B�Ȍ�A�Œ�v�f�Ɋւ��Ă͂������g��������

	//�ÓI�v�Z�_�̔z�u�𓾂�
	if(t==1) set_BEM3D_static_model(CON,s_NODE,s_ELEM,s_node_num,s_elem_num);
	cout<<"�Œ�ߓ_��="<<s_NODE.size()<<"�Œ�v�f��="<<s_ELEM.size()<<endl;
	cout<<*s_elem_num<<endl;

	//���I�v�Z�_�̔z�u�𓾂�
	vector<BEMpoint3D> dy_NODE;
	vector<BEMelement3D> dy_ELEM;
	int d_node_num,d_elem_num;
	set_BEM3D_dynaic_model_for_CONSTANT(CON,PART,dy_NODE,dy_ELEM,&d_node_num,&d_elem_num, particle_number, fluid_number);

	//���҂�����
	vector<BEMpoint3D> NODE;
	vector<BEMelement3D> ELEM;
	vector<REGION> region;
	
	couple3D_NODE_and_ELEM(CON,NODE,ELEM,s_NODE,s_ELEM,dy_NODE,dy_ELEM,region);

	//main�֐����s
	BEM3D_main_for_CONSTANT(CON,NODE,ELEM,PART, fluid_number,F,region);

}

//���v�f�p��BEM�֐�
void BEM3D_main_for_CONSTANT(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<mpsparticle> &PART,int fluid_number,double **F,vector<REGION> &region)
{
	int node_num=int (NODE.size());		//�ߓ_��
	int elemnum=int (ELEM.size());		//�v�f��
	int region_num=(int)(region.size());//�̈�̐�
	int uk=0;								//���m��
	int count;

	cout<<"�̈搔��"<<region_num<<endl;
	
	//���m���v�Z
	for(int i=0;i<elemnum;i++)
	{
		if(ELEM[i].boundary_condition==BOTH) uk+=2;	//BOTH�͑��}�����E��̐ߓ_�ł���A���ݼ�فA�@�������̗��������m�ł��邱�Ƃ�����
		else uk++;									//������͒ʏ�BDiric���낤��Neumn���낤�ƁA�Е��͖��m�Ȃ̂�++;
	}
	cout<<"���m��="<<uk<<endl;

	int *Nid=new int [uk];							//i�Ԗڂ̖��m����i�Ԗڂ̌v�Z�_�ł��邱�Ƃ�����
	int *bd_type=new int [uk];
	int *BOTH_column=new int [elemnum];
	count=0;
	
	for(int i=0;i<elemnum;i++)
	{
		
		if(ELEM[i].boundary_condition==BOTH)
		{
			Nid[i]=i;		//i�Ԗڂ̖��m����i�Ԗڂ̌v�Z�_�ł��邱�Ƃ�����
			bd_type[i]=Neumn;//BOTH�̏ꍇ�A���m�������ݼ�فA�@�������̏��ɒ�`����B�����bd_type�͂��̂悤�ɂ��Ă���
			Nid[elemnum+count]=i;
			bd_type[elemnum+count]=Diric;
			BOTH_column[i]=elemnum+count;//i�Ԗڂ̌v�Z�_��BOTH�̂Ƃ��A�@��������\�����m����BOTH_column[i]�Ԗڂ̖��m���ł���B
			count++;
		}
		else
		{
			//cout<<i<<endl;
			Nid[i]=i;				//i�Ԗڂ̖��m����i�Ԗڂ̌v�Z�_�ł��邱�Ƃ�����
			bd_type[i]=ELEM[i].boundary_condition;
			BOTH_column[i]=-1;		//BOTH�łȂ��Ȃ�_�~�[�Ƃ��ā|�P���i�[
			//cout<<i<<endl;
		}
	}
	
	int gauss_N=4;				//�K�E�X�ϕ��ɂ�����]���_�� //3,4,8
	int NGS[14];				//���Ƃ���Gauss�ϕ��_��4�̂Ƃ���for(int i=NGS[4];i<NGS[4]+4;i++) ep1[i]=~~~~�Ƃ����g������������
	double ep1[14][3];			//�K�E�X�ϕ��ɂ�����Ǐ����W�i�[
	double ep_ln[14];			//���R�ΐ��Ɋւ���K�E�X�ϕ��ɂ�����Ǐ����W�i�[
	double w1[14];				//�K�E�X�ϕ��ɂ�����d�݊i�[
	double w_ln[14];			//���R�ΐ��Ɋւ���K�E�X�ϕ��ɂ�����d�݊i�[

	//�K�E�X�ϕ��̏���
	set_gauss_parameter3D(NGS, ep1,w1,ep_ln,w_ln);

	//matrix�쐬
	double **matrixC=new double*[uk];
	for(int n=0;n<uk;n++) matrixC[n]=new double[uk];
	double **matrixC2=new double*[uk];
	for(int n=0;n<uk;n++) matrixC2[n]=new double[uk];
	double **H=new double*[uk];
	for(int n=0;n<uk;n++) H[n]=new double[uk];
	double **G=new double*[uk];
	for(int n=0;n<uk;n++) G[n]=new double[uk];
	double *Bmatrix=new double[uk];//���s��
	
	for(int n=0;n<uk;n++)
	{
		for(int m=0;m<uk;m++)
		{
			matrixC[n][m]=0;				//������
			matrixC2[n][m]=0;
			H[n][m]=0;
			G[n][m]=0;
		}
		Bmatrix[n]=0;
	}
	cout<<"�s��쐬�J�n"<<endl;
	for(int ID=0;ID<1;ID++)				//�܂��͑�P�̈�̂݌v�Z
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=n1;
			int n2=BOTH_column[n1];
			if(n==2036) cout<<n1<<" "<<ELEM[n1].material<<endl;
				
			double SR=0.5*ELEM[n1].S;//�v�f�ʐ�
					
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
					
				int type=0;//type=0�Ȃ�ʏ�@1�Ȃ�A���̗v�f�͐ߓ_n���܂ޓ��ِϕ�
				if(n1==m1) type=1;
					
				if(type==0)//n=m�̏ꍇ�����l�ϕ�����
				{
					double w;
					double ep[3];		//�]���_�̋Ǐ����W
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						for(int j=0;j<3;j++) ep[j]=ep1[i][j];
						w=w1[i];
							
						set_GHmatrix3D_for_constant(ELEM, NODE,ep,  H, G, w, type,n, m,n1,m1);
					}
					if(m2!=-1)//BOTH�Ȃ�
					{
						G[n][m2]=G[n][m];
						G[n][m]=0;
					}
				}
			}
			H[n][n]=0.5;		//H�̂ݑΊp��������������
			G[n][n]=calc_Gnn_for_CONSTANT(ELEM,NODE, n1);
				

			if(n2!=-1)
			{
				G[n][n2]=G[n][n];
				G[n][n]=0;
			}
		}
	}
	//�@���x�N�g���𔽓]
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<3;D++) ELEM[i].direct[D]*=-1.0;
	
	for(int ID=1;ID<region_num;ID++)				//��2�̈�ȍ~���v�Z
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=BOTH_column[n1];					//�v�Z�_n1�́A�@�������̖��m�ϐ���BOTH_column[n1]�Ԗڂ̖��m��
			
			double SR=0.5*ELEM[n1].S;//�v�f�ʐ�
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
				int type=0;//type=0�Ȃ�ʏ�@1�Ȃ�A���̗v�f�͐ߓ_n���܂ޓ��ِϕ�
				if(n1==m1) type=1;
					
				if(type==0)//n=m�̏ꍇ�����l�ϕ�����
				{
					double w;
					double ep[3];		//�]���_�̋Ǐ����W
						
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						for(int j=0;j<3;j++) ep[j]=ep1[i][j];
						w=w1[i];
							
						set_GHmatrix3D_for_constant(ELEM, NODE,ep,  H, G, w, type,n, m,n1,m1);
					}
					G[n][m2]=-G[n][m]/80;
					//G[n][m2]=-G[n][m]/1;
					G[n][m]=0;
				}
			}
			H[n][n1]=0.5;
			G[n][n]=calc_Gnn_for_CONSTANT(ELEM,NODE, n1);
		}
	}
	//�@���x�N�g���𔽓]
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<3;D++) ELEM[i].direct[D]*=-1.0;
	cout<<"H,G�쐬����"<<endl;

	//���ۂɉ����W���s��쐬

	for(int n=0;n<elemnum;n++)	
	{
		//if(n==2036) cout<<ELEM[n].boundary_condition<<endl;
		if(ELEM[n].boundary_condition==Diric)
		{
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=-G[m][n];
				matrixC2[m][n]=-H[m][n];
			}
		}
		else if(ELEM[n].boundary_condition==Neumn)
		{
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=H[m][n];
				matrixC2[m][n]=G[m][n];
			}
		}
		else if(ELEM[n].boundary_condition==BOTH)
		{
			int n2=BOTH_column[n];
			
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n2]=-G[m][n2];
				//if(H[m][n2]!=0) cout<<"H "<<m<<" "<<n2<<" "<<H[m][n2]<<endl;
			}
				
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=H[m][n];
				if(G[m][n]!=0) cout<<"G "<<m<<" "<<n<<" "<<G[m][n]<<endl;
			}
		}
		//cout<<matrixC[2036][2036]<<endl;
	}
	cout<<"�W���s��쐬����"<<endl;
	cout<<matrixC[2036][2036]<<endl;
	cout<<G[2036][2036]<<endl;
	cout<<H[2036][2036]<<endl;

	//�����艺�ŃG���[�B�ǂ��H�H�H
	//���s��쐬
	for(int n=0;n<uk;n++)
	{
		double val=0;//���s��̒l
		for(int m=0;m<uk;m++)
		{
			int elem=Nid[m];			//���m��n�ɑΉ�����v�f�ԍ�
			if(bd_type[m]==Diric)
			{
				int node=ELEM[elem].node[0];	//�v�f���ިظڌ^�Ȃ�A3���_���ިظڌ^�ŁA�����ިظڒl���������B�Ȃ̂�node[0]���ިظڒl���g��
				val+=matrixC2[n][m]*NODE[node].potential;
			}
			else if(bd_type[m]==Neumn)
			{
				int node=ELEM[elem].node[0];	//�v�f���m�C�}���^�̂Ȃ�A�m�C�}���^�̐ߓ_��T���Ă��̔����l���g�p
				if(NODE[node].boundary_condition==Neumn) node=ELEM[elem].node[1];
				if(NODE[node].boundary_condition==Neumn) node=ELEM[elem].node[2];

				val+=matrixC2[n][m]*NODE[node].slop1;
			}
		}
		Bmatrix[n]=val;
	}
	/////matrix�쐬����

	cout<<"�s��쐬���� �K�E�X�̏����@---";

	gauss(matrixC,Bmatrix,uk);//������Bmatrix�Ɋi�[

	cout<<"ok"<<endl;

	//�����i�[
	for(int n=0;n<uk;n++)
	{
		int elem=Nid[n];		//���m��n�ɑΉ�����v�f�ԍ�
		int node=ELEM[elem].node[0];
		if(bd_type[n]==Neumn) NODE[node].potential=Bmatrix[n];
		else if(bd_type[n]==Diric) NODE[node].slop1=Bmatrix[n]; 
	}

	//�����m�F
	ofstream fp3("V.dat");
	ofstream fp4("slop.dat");
	double le=CON->get_distancebp();
	for(int n=0;n<elemnum;n++)
	{
		int node=ELEM[n].node[0];
		double Ys=ELEM[n].r[A_Y];
		if(Ys>0 && Ys<5*le)
		{
			fp3<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Z]<<" "<<NODE[node].potential<<endl;
			fp4<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Z]<<" "<<NODE[node].slop1<<endl;
		}
	}
	fp3.close();
	fp4.close();

	
/*
	ofstream fp5("F.dat");
	double ep0=8.854e-12;
	double times=1e-3;
	double le=CON->get_distancebp();
	for(int n=0;n<elemnum;n++)
	{
		if(ELEM[n].material==FLUID)
		{
			int node=ELEM[n].node[0];//�Ή�����ߓ_�ԍ�
			int i=NODE[node].particle;//�Ή����闱�q�ԍ�
			double E=NODE[node].slop1;
			double Fs=0.5*ep0*E*E;	//����
			double Fn=Fs*ELEM[n].L;	//��[N]
			for(int D=0;D<2;D++) F[D][i]=-Fn*ELEM[n].direct[D];
			fp5<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<F[A_X][i]*times/le<<" "<<F[A_Y][i]*times/le<<endl;
			//fp5<<NODE[node].r[A_X]<<" "<<NODE[node].r[A_Y]<<" "<<F[A_X][i]*times/le<<" "<<F[A_Y][i]*times/le<<endl;
		//	fp5<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<Fs<<endl;
		}
	}
	fp5.close();

	////�Ѱ��ݸ�
	smoothingF3D(CON,PART,fluid_number,F);

	if(CON->get_dir_for_P()==2 ||CON->get_dir_for_P()==3 )
    {
		ofstream bb("electromagnetic_P.dat");
		for(int i=0;i<fluid_number;i++)
		{
			double fs=0;//�\�ʗ�
			if(PART[i].surface==ON)//�������̂̏ꍇ��fs=0�Ƃ���
			{
				double Fn=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]);
				fs=Fn/le;
				for(int D=0;D<3;D++) F[D][i]=0;//���̓fިظڂƂ��ēd���͂��g�p����̂ł����ł͏�����
			}
			bb<<-fs<<endl;
        }
		bb.close();
	}

*/
	for(int n=0;n<node_num;n++) delete[] matrixC[n];
	delete [] matrixC;
	for(int n=0;n<node_num;n++) delete[] matrixC2[n];
	delete [] matrixC2;
	for(int n=0;n<node_num;n++) delete[] H[n];
	delete [] H;
	for(int n=0;n<node_num;n++) delete[] G[n];
	delete [] G;
	delete [] Bmatrix;

	delete [] Nid;
	delete [] bd_type;
	delete [] BOTH_column;
}

//�O�p�`�̈�ɂ�������Gauss�ϕ��̏����֐�
void set_gauss_parameter3D(int *NGS,double ep1[14][3],double *w1,double *ep_ln,double *w_ln)
{
	//�K�E�X�ϕ��̏��� �l��[���E�v�f�@-��b�Ɖ��p- :J.T.�J�`�J�f�[���X]��p222�����p
	for(int i=0;i<32;i++) NGS[i]=0;
	NGS[3]=0; NGS[4]=3; NGS[7]=7;	//���Ƃ���Gauss�ϕ��_��4�̂Ƃ���for(int i=NGS[4];i<NGS[4]+4;i++) ep1[i][0]=~~~~�Ƃ����g������������

	//gauss_N=3�̂Ƃ�
	ep1[0][0]=0.5;		//�]���_1��\���ʐύ��W
	ep1[0][1]=0.5;
	ep1[0][2]=0;
	ep1[1][0]=0;		//�]���_2��\���ʐύ��W
	ep1[1][1]=0.5;
	ep1[1][2]=0.5;
	ep1[2][0]=0.5;		//�]���_3��\���ʐύ��W
	ep1[2][1]=0;
	ep1[2][2]=0.5;
	w1[0]=1.0/3;		//�]���_1�ł̏d��
	w1[1]=1.0/3;		//�]���_2�ł̏d��
	w1[2]=1.0/3;		//�]���_3�ł̏d��

	//gauss_N=4�̂Ƃ�
	ep1[3][0]=1.0/3;	//�]���_1��\���ʐύ��W
	ep1[3][1]=1.0/3;
	ep1[3][2]=1.0/3;
	ep1[4][0]=0.6;		//�]���_2��\���ʐύ��W
	ep1[4][1]=0.2;
	ep1[4][2]=0.2;
	ep1[5][0]=0.2;		//�]���_3��\���ʐύ��W
	ep1[5][1]=0.6;
	ep1[5][2]=0.2;
	ep1[6][0]=0.2;		//�]���_4��\���ʐύ��W
	ep1[6][1]=0.2;
	ep1[6][2]=0.6;
	w1[3]=-0.5625;		//�]���_1�ł̏d��
	w1[4]=25.0/48.0;	//�]���_2�ł̏d��
	w1[5]=25.0/48.0;	//�]���_3�ł̏d��
	w1[6]=25.0/48.0;	//�]���_4�ł̏d��

	//gauss_N=7�̂Ƃ�
	ep1[7][0]=1.0/3;	//�]���_1��\���ʐύ��W
	ep1[7][1]=1.0/3;
	ep1[7][2]=1.0/3;
	ep1[8][0]=0.79742699;		//�]���_2��\���ʐύ��W
	ep1[8][1]=0.10128651;
	ep1[8][2]=0.10128651;
	ep1[9][0]=0.10128651;		//�]���_3��\���ʐύ��W
	ep1[9][1]=0.79742699;
	ep1[9][2]=0.10128651;
	ep1[10][0]=0.10128651;		//�]���_4��\���ʐύ��W
	ep1[10][1]=0.10128651;
	ep1[10][2]=0.79742699;
	ep1[11][0]=0.05971587;		//�]���_5��\���ʐύ��W
	ep1[11][1]=0.47014206;
	ep1[11][2]=0.47014206;
	ep1[12][0]=0.47014206;		//�]���_6��\���ʐύ��W
	ep1[12][1]=0.05971587;
	ep1[12][2]=0.47014206;
	ep1[13][0]=0.47014206;		//�]���_7��\���ʐύ��W
	ep1[13][1]=0.47014206;
	ep1[13][2]=0.05971587;
	w1[7]=0.225;		//�]���_1�ł̏d��
	w1[8]=0.12593918;	//�]���_2�ł̏d��
	w1[9]=0.12593918;	//�]���_3�ł̏d��
	w1[10]=0.12593918;	//�]���_4�ł̏d��
	w1[11]=0.13239415;	//�]���_5�ł̏d��
	w1[12]=0.13239415;	//�]���_6�ł̏d��
	w1[13]=0.13239415;	//�]���_7�ł̏d��

}

//���v�f�̂��߂�H�}�g���N�X�v�Z�֐�
void set_GHmatrix3D_for_constant(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,double ep[3],double **H,double **G,double w,int type,int n,int m,int n1,int m1)
{
	//n1:�\�[�X�_�v�f�ԍ� n:n1�ɊY������s��v�f�ԍ�
	//m1:����v�f�ԍ� m:m1�ɊY������s��v�f�ԍ�
	//type=0�Ȃ�ʏ� 1�Ȃ�v�Z���Ȃ��B
	
	//Xms,Yms:����v�f�̒��_
	//ep[i]:�]���_�̖ʐύ��W
	//w �]���_�ɂ�����d��

	double S=ELEM[n1].S;		//�v�f�ʐ�

	double Xs[3]={ELEM[n1].r[A_X],ELEM[n1].r[A_Y],ELEM[n1].r[A_Z]};//�\�[�X�_���W(�v�f�̒��_)

	double Gc[3];				//����v�f��Gauss�ϕ��]���_���W
	for(int D=0;D<3;D++) Gc[D]=ep[0]*NODE[ELEM[m1].node[0]].r[D]+ep[1]*NODE[ELEM[m1].node[1]].r[D]+ep[2]*NODE[ELEM[m1].node[2]].r[D];

	double RA=0;				//�\�[�X�_�ƕ]���_�̋���
	for(int D=0;D<3;D++) RA+=(Xs[D]-Gc[D])*(Xs[D]-Gc[D]);
	RA=sqrt(RA);

	double RD[3];				//�\�[�X�_����]���_�֌������P�ʃx�N�g��
	for(int D=0;D<3;D++) RD[D]=(Gc[D]-Xs[D])/RA;

	double RDN=0;				//��L�x�N�g���Ɩ@���x�N�g���Ƃ̓���
	for(int D=0;D<3;D++) RDN+=RD[D]*ELEM[n1].direct[D];

	if(type==0)
	{
		H[n][m]+=-1.0/(4*PI*RA*RA)*RDN*w*S;
		G[n][m]+=1.0/(4*PI*RA)*w*S;
	}
}

//���v�f���g�p���邳����G�̓��Ӑϕ��v�Z�֐�
double calc_Gnn_for_CONSTANT(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,int n1)
{
	//�l��[���E�v�f���-���_�Ɖ��p:�c������]��p93���Q�Ƃ���
	//n1:�\�[�X�_�v�f�ԍ� n:n1�ɊY������s��v�f�ԍ�
	//m1:����v�f�ԍ� m:m1�ɊY������s��v�f�ԍ�
	double val=0;
	int ia=ELEM[n1].node[0];
	int ib=ELEM[n1].node[1];
	int ic=ELEM[n1].node[2];

	double Gp[3];				//�O�p�`�̏d�S
	for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;

	double H[3]={0,0,0};//�d�S�Ɗe���_�Ƃ̋���
	for(int D=0;D<3;D++)
	{
		H[0]+=(Gp[D]-NODE[ia].r[D])*(Gp[D]-NODE[ia].r[D]);//�d�S��ia�Ƃ̋���
		H[1]+=(Gp[D]-NODE[ib].r[D])*(Gp[D]-NODE[ib].r[D]);//�d�S��ib�Ƃ̋���
		H[2]+=(Gp[D]-NODE[ic].r[D])*(Gp[D]-NODE[ic].r[D]);//�d�S��ic�Ƃ̋���
	}
	for(int j=0;j<3;j++) H[j]=sqrt(H[j]);

	double r12,r23,r31;		//�O�p�`�̕ӂ̒���
	r12=0; r23=0; r31=0;
	for(int D=0;D<3;D++)
	{
		r12+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);//ia-ib�̋���
		r23+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);//ib-ic�̋���
		r31+=(NODE[ic].r[D]-NODE[ia].r[D])*(NODE[ic].r[D]-NODE[ia].r[D]);//ic-ia�̋���
	}
	r12=sqrt(r12);	r23=sqrt(r23);	r31=sqrt(r31);


	double theta[3];//��`�͋��ȏ��Q��
	//�]���藝���Ƃ����߂�
	theta[0]=(H[1]*H[1]+H[2]*H[2]-r23*r23)/(2*H[1]*H[2]);		//���̎��_�ł�cos�Ƃ��i�[����Ă��邱�Ƃɒ���
	theta[1]=(H[0]*H[0]+H[2]*H[2]-r31*r31)/(2*H[0]*H[2]);
	theta[2]=(H[1]*H[1]+H[0]*H[0]-r12*r12)/(2*H[1]*H[0]);
	for(int j=0;j<3;j++) theta[j]=acos(theta[j]);

	double alpha[3];//��`�͋��ȏ��Q��
	for(int j=0;j<3;j++) alpha[j]=0;
	for(int D=0;D<3;D++)
	{
		alpha[0]+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ic].r[D]-NODE[ia].r[D]);
		alpha[1]+=(NODE[ia].r[D]-NODE[ib].r[D])*(NODE[ic].r[D]-NODE[ib].r[D]);
		alpha[2]+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
	}
	alpha[0]/=r12*r31;		//���̒i�K�ł�cos(2��)���i�[����Ă���
	alpha[1]/=r12*r23;
	alpha[2]/=r23*r31;
	for(int j=0;j<3;j++)
	{
		alpha[j]=acos(alpha[j]);//���̒i�K�ł�2�����i�[����Ă���
		alpha[j]*=0.5;
	}//alpha�����܂���

	double SS=ELEM[n1].S;//�ʐ�

	double termA,termB,termC;
	termA=0;termB=0;termC=0;
	termA=log(tan((theta[0]+alpha[1])*0.5)/(tan(alpha[1]*0.5)))/r23;
	termB=log(tan((theta[1]+alpha[2])*0.5)/(tan(alpha[2]*0.5)))/r31;
	termC=log(tan((theta[2]+alpha[0])*0.5)/(tan(alpha[0]*0.5)))/r12;

	val=2*SS/3*(termA+termB+termC);

	return val;

}
