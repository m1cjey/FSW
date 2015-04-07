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
#include"MMM_CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢


///���C���[�����g�@�v�Z�֐�
void Magnetic_Moment_Method(mpsconfig *CON,vector<mpsparticle> &PART,double **F,double n0,double lamda,int fluid_number,int particle_number)
{
	
	double u0=4*PI*1e-7;
	MMM_config MCON;

	cout<<"���C���[�����g�@---";

	double le=CON->get_distancebp();
	double R=CON->get_re()*le;	//���U�p�e�����a
	int d=CON->get_dimention();						//����
	unsigned int timeA=GetTickCount();				//�v�Z�J�n����
	int count=0;
	int pn=fluid_number*d;							//���m��:���q���~���xD(����)���� 
	double co=1/(4*PI*u0);							//�v�Z�ɂ悭�����W��
	double RP0=CON->get_RP();						//�䓧����
	double kai=1.0-RP0;								//���C����
	double V=CON->get_particle_volume();
	double S=le;
	if(d==3 && CON->get_model_set_way()==0) S=le*le;
	if(d==3 && CON->get_model_set_way()==1) S=sqrt(3.0)/4*le*le;//�f�ʐ�

	//�s��m�� //���m���̏��Ԃ�Mx0,My0,Mz0,Mx1,My1,Mz1,Mx2,My2,Mz2�E�E�E
	double **matrix=new double*[pn];
	for(int i=0;i<pn;i++) matrix[i]=new double[pn];
	double *B   = new double[pn];					//���s��

	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//�O�����@���x�N�g��

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) direct[D][i]=0;				//������
        if(PART[i].surface==ON)
		{
			direct_f(CON,PART,i,direct);
			for(int D=0;D<d;D++) direct[D][i]*=-1;//�O�������~�������甽�]����
		}
	}

	for(int i=0;i<pn;i++)
	{
		B[i]=0;
		for(int j=0;j<pn;j++) matrix[i][j]=0;		//������
	}
	//////////////////////////////////////////////////////////////////

	//���s��B[]�쐬
	if(MCON.get_Hf_type()==0)
	{
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<d;D++) B[i*d+D]=MCON.get_Hf_H()*sign[D];//2�����Ȃ�A_Y,3�����Ȃ�A_Z�̕��������l������B
		}
	}
	else cout<<"�O������̃^�C�v��������"<<endl;
	/////*/

	
	/////////////////////matrix�ɒl���i�[
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<d;D++)
		{
			int n=i*d+D;//�Ή�����matrix���v�f�ԍ�
			matrix[n][n]+=1/(RP0-1);
		}
		for(int j=0;j<fluid_number;j++)
		{
			if(j!=i)
			{
				double rA[3]={PART[i].r[A_X]-PART[j].r[A_X],PART[i].r[A_Y]-PART[j].r[A_Y],PART[i].r[A_Z]-PART[j].r[A_Z]};
				double disA=sqrt(rA[A_X]*rA[A_X]+rA[A_Y]*rA[A_Y]+rA[A_Z]*rA[A_Z]);//�v�Z�_i�ƌv�Z�_j�̋���
	
				double dis3=disA*disA*disA;		//������3��
				
				//if(PART[j].surface==OFF)
				{
					double co2=1/dis3*co*V;	//�悭�g���W��
					for(int D=0;D<d;D++)
					{
						int n=i*d+D;//�Ή�����matrix���v�f�ԍ�
						for(int k=0;k<PART[j].N;k++)
						{
							int j2=PART[j].NEI[k];	//�v�Z�_j�̎��ӗ��q
							double rB[3]={PART[j2].r[A_X]-PART[j].r[A_X],PART[j2].r[A_Y]-PART[j].r[A_Y],PART[j2].r[A_Z]-PART[j].r[A_Z]};
							double disB=sqrt(rB[A_X]*rB[A_X]+rB[A_Y]*rB[A_Y]+rB[A_Z]*rB[A_Z]);//�v�Z�_j�Ƃ��̋ߗח��q�Ƃ̋���
							double w=kernel(R,disB);
							for(int D2=0;D2<d;D2++)
							{
								int m=j2*d+D2;
								int l=j*d+D2;
								matrix[n][m]+=d/n0*co2*rB[D2]/(disB*disB)*w*rA[D];
								matrix[n][l]+=-d/n0*co2*rB[D2]/(disB*disB)*w*rA[D];
							}
						}	
					}
				}
				/*else if(PART[j].surface==ON)
				{
					double co2=1/dis3*co*S;	//�悭�g���W��
					for(int D=0;D<d;D++)
					{
						int n=i*d+D;//�Ή�����matrix���v�f�ԍ�
						for(int D2=0;D2<d;D2++)
						{
							int l=j*d+D2;
							matrix[n][l]+=-co2*direct[D2][j]*rA[D];
						}	
					}
				}///////*/
			}
		}
	}
	cout<<"�s��쐬"<<endl;

	//�s��l�o��
	//output_matrix(matrix,B, pn);

	cout<<"���m��:"<<pn<<"�ɑ΂��K�E�X�̏����@--";
	unsigned int timeB=GetTickCount();
	gauss(matrix,B,pn);//������B�Ɋi�[
	//jacobi(matrix,B,pn);//������B�Ɋi�[
	cout<<"ok time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

	double *M[DIMENTION];					
	for(int D=0;D<DIMENTION;D++) M[D]=new double [fluid_number];
	double *H[DIMENTION];					
	for(int D=0;D<DIMENTION;D++) H[D]=new double [fluid_number];//���q�ʒu�ł̎���H
	
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			M[D][i]=0;
			H[D][i]=0;
		}
		for(int D=0;D<d;D++)
		{
			M[D][i]=B[i*d+D];
			H[D][i]=M[D][i]/kai;
		}
	}
	ofstream fp("M.dat");
	ofstream fq("H.dat");
	double times=1e-11;
	if(d==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<M[A_X][i]*times<<" "<<M[A_Y][i]*times<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i])<<endl;
		}
	}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le)
			{
				fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]*times<<" "<<M[A_Z][i]*times<<endl;
				fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
			}
		}
	}
	fp.close();
	fq.close();

	//�͂����߂�
	
	double *Fs[DIMENTION];
    for(int D=0;D<DIMENTION;D++) Fs[D]=new double [particle_number];//�P�ʖʐς�����̗�
	double *Fv[DIMENTION];
    for(int D=0;D<DIMENTION;D++) Fv[D]=new double [particle_number];//�P�ʑ̐ς�����̗�
	double *Hgrad[DIMENTION];			//H���z�޸�يi�[
	for(int D=0;D<DIMENTION;D++) Hgrad[D]=new double [fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Fs[D][i]=0;		//������
		if(PART[i].surface==ON)
		{
			double Mn=M[A_X][i]*direct[A_X][i]+M[A_Y][i]*direct[A_Y][i]+M[A_Z][i]*direct[A_Z][i];
			double val=0.5*u0*Mn*Mn;		//���͒l
			for(int D=0;D<DIMENTION;D++) Fs[D][i]=val*direct[D][i];
			for(int D=0;D<DIMENTION;D++) F[D][i]=Fs[D][i]*CON->get_distancebp();
		}
	}
	//�̐ϗ�
	H_gradient1(CON,PART, fluid_number,Hgrad,H);//��H�v�Z
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) Fv[D][i]=u0*M[D][i]*Hgrad[D][i];
	
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=Fv[D][i]*V;

	////�Ѱ��ݸ�
	//smoothingF3D(CON,PART,fluid_number,F);

	ofstream fr("Fs.dat");
	times=5e-2/MCON.get_Hf_H();
	if(d==2){ for(int i=0;i<fluid_number;i++) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Y][i]*times<<endl;}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Z][i]*times<<endl;
	}
	fr.close();
	ofstream ft("Fv.dat");
	times=1e-6;
	if(d==2) for(int i=0;i<fluid_number;i++) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fv[A_X][i]*times<<" "<<Fv[A_Y][i]*times<<endl;
	ft.close();

	if(CON->get_dir_for_P()==2 ||CON->get_dir_for_P()==3 )
    {
		ofstream bb("electromagnetic_P.dat");
		for(int i=0;i<fluid_number;i++)
		{
			double fs=0;//�\�ʗ�
			if(PART[i].surface==ON)
			{
				fs=sqrt(Fs[A_X][i]*Fs[A_X][i]+Fs[A_Y][i]*Fs[A_Y][i]+Fs[A_Z][i]*Fs[A_Z][i]);
				for(int D=0;D<DIMENTION;D++) F[D][i]=0;
			}
			bb<<-fs<<endl;
        }
		bb.close();
	}

	for(int i=0;i<pn;i++) delete [] matrix[i];
	delete [] matrix;
	
    delete [] B;
	for(int D=0;D<DIMENTION;D++) delete [] M[D];
	for(int D=0;D<DIMENTION;D++) delete [] H[D];
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] Fs[D];
	for(int D=0;D<DIMENTION;D++) delete [] Fv[D];
	for(int D=0;D<DIMENTION;D++) delete [] Hgrad[D];

	
}

///H���z�v�Z�֐�ver.1
void H_gradient1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Hgrad,double **H)
{
	double le=CON->get_distancebp();//�������q�ԋ���
	double r=CON->get_re()*le;
	int d=CON->get_dimention();

	double *HH=new double[fluid_number];//�e���q�ʒu�ł̎��ꋭ�xH�i�[
	for(int i=0;i<fluid_number;i++) HH[i]=sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i]);
	

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Hgrad[D][i]=0;//������

		double W=0;//���q�����x�@OUT���������肷��̂�PND[i]�͔���

		for(int k=0;k<PART[i].N;k++)
		{       
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
		
				double w=kernel(r,dis);
				W+=w;
				
				Hgrad[A_X][i]+=(HH[j]-HH[i])*X*w/(dis*dis);
				Hgrad[A_Y][i]+=(HH[j]-HH[i])*Y*w/(dis*dis);
				Hgrad[A_Z][i]+=(HH[j]-HH[i])*Z*w/(dis*dis);
			}

		}
		for(int D=0;D<DIMENTION;D++) if(W!=0) Hgrad[D][i]= Hgrad[D][i]*d/W;
	}///////////////Pgrad[D][i]�v�Z�I��

	

	delete [] HH;
}

//�s��l�o�͊֐�(�Ǘ��p)
void output_matrix(double **matrix,double *Bmatrix,int node_num)
{
	//node_num=30;					//���ゾ���\����������
	int *DDN=new int[node_num];		//�Ίp�D�ʁ@diagonally dominant matrix
	ofstream fp2("matrixMMM.dat");
	for(int n=0;n<node_num;n++)
	{
		double val=0;
		DDN[n]=OFF;
		for(int m=0;m<node_num;m++)
		{
			fp2<<matrix[n][m]<<"\t";
			if(n!=m) val+=sqrt(matrix[n][m]*matrix[n][m]);
		}
		fp2<<endl;
		if(val<sqrt(matrix[n][n]*matrix[n][n])) DDN[n]=ON;//�Ίp�D��
	}
	fp2.close();
	ofstream fp4("BmatrixMMM.dat");
	for(int n=0;n<node_num;n++) fp4<<Bmatrix[n]<<endl;
	fp4.close();

	//for(int n=0;n<node_num;n++) if(DDN[n]==ON) cout<<n<<"�͑Ίp�D��"<<endl;

	delete [] DDN;
}

//jacobi�̔����@ ���͍ŏI�I��B�̂Ȃ���
void jacobi(double **matrix,double *B,int N)
{
	double ep=1e-8;//��������
	double E=10;		//�덷
	int count=0;
	
	double *X=new double[N];//��
	for(int k=0;k<N;k++) X[k]=0;//�����l
	while(E>ep)
	{
		E=0;
		for(int i=0;i<N;i++)
		{
			double L=0;
			double U=0;
			for(int j=0;j<i;j++) L+=matrix[i][j]*X[j]; 
			for(int j=i+1;j<N;j++) U+=matrix[i][j]*X[j];
			double Xnew=(B[i]-L-U)/matrix[i][i];
			E+=fabs(Xnew-X[i]);
			X[i]=Xnew;
		}
		E/=N;
		count++;
		cout<<count<<" E="<<E<<endl;
	}
	
	for(int k=0;k<N;k++) B[k]=X[k];//B�ɓ������i�[
	delete [] X;
}