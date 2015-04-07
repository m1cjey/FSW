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

void P3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,int *jnb,vector<mpsparticle> &PART,int fluid_number,int **nei,double dt,double N0);
void reU3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,vector<mpsparticle> &PART,int fluid_number,double dt,int *jnb);


//�A���
void negative1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un)
{
	int d=CON->get_dimention();
	double *reU[DIMENTION];//���x�C����
	int negativeP=CON->get_negativeP();			//���̈��͂������邩�A���Ȃ���
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	/////���͌v�Z
	pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,CON->get_B_of_P(),negativeP);

	//���͌��z�v�Z
	calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,CON->get_minP());
	///////////////
	
	//���x�C���ƈʒu�C��
	modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
	
	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
}

//�A��͂�2��i�P��ڂ�PND�@�Q��ڂ͑��x���U�j
void negative1_twice(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un,double n0_4)
{
	int d=CON->get_dimention();
	double le=CON->get_distancebp();
	double limit=CON->get_Pgrad_limit()*le;		//�ʒu�C���ʂ̌��E�l
	int negativeP=CON->get_negativeP();			//���̈��͂������邩�A���Ȃ���
	double *reU[DIMENTION];//���x�C����
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	//if(t%2==0)
	if(t%CON->get_P_twice()==0)
	{

		//calc_PND_by_minmum_L_main(CON,PART, particle_number);

		/////���͌v�Z
		//negativeP=OFF;
		pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,0,negativeP);//���q�����x�������Ōv�Z
	
		//���͌��z�v�Z
		int minPsw=OFF;
		calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,minPsw);

		//�ʒu�̂ݏC��
		int modified_num=0;						//�ʒu�C���ʂ��C�����ꂽ���q��
		for(int i=0;i<fluid_number;i++)
		{
			double ABS=0;		//�ʒu�C���ʂ̐�Βl
			for(int D=0;D<d;D++) ABS+=dt*reU[D][i]*dt*reU[D][i];
			ABS=sqrt(ABS);
			if(ABS>limit)
			{
				for(int D=0;D<d;D++) reU[D][i]*=limit/ABS;		//�C���ʂ̑傫�����C��
				modified_num++;
			}
		
			for(int D=0;D<d;D++) PART[i].r[D]+=dt*reU[D][i];
		}

		if(t==1 && CON->get_restart()==OFF)
		{
			ofstream fout("Pgrad_modified_num.dat");
	
			fout<<t<<" "<<modified_num<<endl;		
			fout.close();
		}
		else if(modified_num>0)
		{
			ofstream avs("Pgrad_modified_num.dat",ios :: app);
			avs<<t<<" "<<modified_num<<endl;	
			avs.close();
		}

		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].r[D]+PART[i].r[D]!=2*PART[i].r[D]) cout<<"## i="<<i<<endl;
	
		//���x�C���ƈʒu�C��
		//modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
		
		//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].u[D]-=reU[D][i];//���x�������Ƃɖ߂�
	
		//���q�ʒu���ύX���ꂽ�̂ŁAfreeon�����s���Ȃ��Ƃ����Ȃ�
		int *surface=new int[particle_number];		//�������\�ʂ̒�`�͕ύX���Ăق����Ȃ��B�����Œl��ۑ����Ă���
		for(int i=0;i<particle_number;i++) surface[i]=PART[i].surface;
		calc_neighbor_relation(CON,PART,particle_number,n0_4,fluid_number,out);
		for(int i=0;i<particle_number;i++) PART[i].surface=surface[i];
		delete [] surface;
	}	

	//�ēx�v�Z
	negativeP=CON->get_negativeP();		
	pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,CON->get_B_of_P(),negativeP);

	//���͌��z�v�Z
	calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,CON->get_minP());
	///////////////
	
	//���x�C���ƈʒu�C��
	modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
	
	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
	
}
///���͌v�Z�֐�
void pressure(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,double dt,int t,double lamda,double N0,double *PND2,double n0,int B_of_P,int negativeP)
{
	int P_flag=OFF;	//���͌v�Z���s�����ۂ�
	int pn=0;		//���m��
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			P_flag=ON;//���͌v�Z���s��
			pn++;
		}
		else if(PART[i].type==INWALL && PART[i].surface==OFF)
		{
			P_flag=ON;//���͌v�Z���s��
			pn++;
		}
	}
	
	if(P_flag==ON && fluid_number>0)
	{
		if(CON->get_gridless_P()==OFF) calc_P_main(CON,PART,particle_number,lamda,N0,dt,PND2,n0,fluid_number,out,pn,B_of_P);//MPS�ɂ�鈳�͌v�Z
		if(CON->get_gridless_P()==ON) calc_P_main_with_gridless(CON,PART,particle_number,lamda,N0,dt,PND2,n0,fluid_number,out,pn,B_of_P,t);//gridless�@�ɂ�鈳�͌v�Z
	} 

    ///���͂���������
	//if(CON->get_set_P()!=OFF) set_P(CON,PART,particle_number,fluid_number);
	/////*/
	
	////���͂�splot
	if(P_flag==ON) plot_P(CON ,PART,particle_number,t,fluid_number);
	////////////////////////////
}

void calc_P_main(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P)
{

	///���͌v�Z�ɂ����āA�g�p����d�݊֐��͌��z�ȂǂɎg�p���邻��Ƃ͈قȂ���̂�p����B
	//���R�́A���z�Ɏg�p����d�݊֐��͒P�Ȃ�d�ݕ��ςȂ̂ŉ������悤���Ă��悢���A
	//PPE���̗��U���Ɏg�p����d�݊֐��́A���q�����x�喧�x�ƂȂ�悤�ɐݒ肷��K�v�����邩��B
	//�������҂œ���̏d�݊֐����g�p��������΁Akernel2()�̒��g��kernel�Ɠ���ɏ��������邱�ƁB

	//pn:���͂��������߂̘A�����������m��

	double le=CON->get_distancebp();
	double r2=CON->get_re2()*le;
	unsigned int timeA=GetTickCount();					//�v�Z�J�n����
	double dimention=2;
	if(CON->get_dimention()==3) dimention=3;
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();

	//���͉�͂̂��߂�N0��PART[i].PND2��kernel2()�p�ɂ���������
	//set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	//cout<<"N0="<<N0<<endl;
	///////////////////*/

	cout<<"���͖��m��:"<<pn<<" ";

	int *ppn = new int[pn];					//�s��ɂ������n�Ԗڂ̖��m���͗��q�ԍ�ppn[n]�̗��q�ɑ���
	int *link = new int [particle_number];	//���q�ԍ�i��link[i]�Ԗڂ̖��m��
	double *B   = new double[pn];			//���s��
	
	int count=0;
	///ppn�z���link�z��쐬
	
	double real_lamda=lamda;	//lamda�̒l�L��
	if(CON->get_HL_sw()==ON) lamda=2.0*CON->get_dimention()*real_lamda;//�����̗��U���̏ꍇ�Alamda=2d�ƍĒ�`����΃�/2d�̍��͏�����B��������������Ɨ��ӂ̒l���傫���Ȃ肷����̂ŁAreak_lamda�𗼕ӂɂ����ď��������Ă���
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL && PART[i].surface==OFF)
		{
			ppn[count]=i;
			if(B_of_P==0)
			{
				B[count]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//���ȏ�
			}
			else if(B_of_P==1)//���x���U
			{    
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				//cout<<"D="<<div<<" "<<PART[i].PND2<<endl;
				
				//B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
				B[count]=div*lamda*PART[i].PND2/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==2)//(���ȏ�+���x���U)/2
			{
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				
				B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[count]/=a+b;
			}	
			else if(B_of_P==3)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U
			{    
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==4)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U+PND
			{    
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//���ȏ�
				B[count]/=a+b;
			}
			else if(B_of_P==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				//double div=Dndt(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			//if(B[count]>800000) cout<<"B="<<i<<endl;
			link[i]=count;
			count++;
		}
		else link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
	}
	lamda=real_lamda;//lamda�̒l��߂�
	//////*/

	if(CON->get_dir_for_P()!=OFF && B_of_P!=0)//�\�ʗ��q�̈��͂Ƃ��āA�f�B���N���l
	{
		int flag=CON->get_dir_for_P();
		if(flag==2 || flag==3)
		{
			if(CON->get_EM_method()==OFF) flag=1;//BEM�ɂ��v�Z���s��Ȃ��̂�flag��2��3�Ȃ�ԈႢ�Ȃ̂�1�ɖ߂�
		}
		double *Dirichlet_P=new double [particle_number];
		for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

		if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(flag==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}
		
		for(int i=fluid_number;i<particle_number;i++)//�Ǖ\�ʗ��q��Dirichlet_P�v�Z
		{
			//cout<<i<<endl;
			if(PART[i].surface==ON)//�Ǖ\�ʗ��q�Ȃ�
			{
				double P=0;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					if(PART[j].type==FLUID && PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(r2,dis);
						P+=Dirichlet_P[j]*w;
						W+=w;
					}
				}
				if(W!=0) P/=W;
				Dirichlet_P[i]=P;
			}
		}

		ofstream gg2("Dirichlet_error.dat");
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					int m=link[j];
					if(m<pn)//���qj�����m���Ȃ�
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
				    
						double w=kernel2(r2,dis,dimention);
						if(CON->get_HL_sw()==ON) w*=real_lamda;
						B[m]-=w*Dirichlet_P[i];
						PART[i].P=Dirichlet_P[i];
					}
				}
			}
			else if(PART[i].type==FLUID && Dirichlet_P[i]!=0)
			{
				cout<<"�������q�Ȃ̂��ިظڒl���i�[����Ă��܂��B�l��"<<Dirichlet_P[i]<<endl;
				gg2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		gg2.close();

		//̧�ُo��
		ofstream gg("Dirichlet_P.dat");
		if(dimention==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
		else if(dimention==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
		gg.close();///*/

		//�f�B���N���l���x�N�g���\��
		output_dirichlet_vector_files(CON,PART,fluid_number,flag,Dirichlet_P);

		delete [] Dirichlet_P;
	}///////////

	///���s��o��
	ofstream h("Bmatrix_for_P.dat");
	if(CON->get_dimention()==2) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<B[n]<<endl;
		}
	}
	if(CON->get_dimention()==3) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B[n]<<endl;
		}
	}
	h.close();
	/////////////*/

	
	int number=0;			//�W���s��̔�[���v�f��
	for(int n=0;n<pn;n++)
	{   
	    number++;///�������g�𐔂ɂ����
	    int i=ppn[n];//n�Ԗڂ̖��m���͗��qi
	  
	    for(int k=0;k<PART[i].N2;k++)
	    {
	        int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
	    }
	}
	////number�����܂���
	
    double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	
	if(CON->get_HL_sw()==OFF)//�W���I�ȗ��U��
	{
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//�l��ۑ�
			ind[index]=n;
			index++;
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				if(m<pn)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					    
					double w=kernel2(r2,dis,dimention);
					val[index]=w;
					ind[index]=m;
					index++;
					W+=w;
				}
				else if(PART[j].surface==ON)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					
					//double w=kernel(r2,dis);
					double w=kernel2(r2,dis,dimention);
					W+=w; //������w���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^ 
				}
			}
			val[kk]=-W;
		}
		ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
	}
	else if(CON->get_HL_sw()==ON)//�����̗��U�� �d�݊֐���ς����炱���������ς��Ȃ��Ƃ����Ȃ����Ƃɒ���
	{
		if(CON->get_dimention()==2) cout<<"2D�͔�Ή�"<<endl;
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//�l��ۑ�
			//val[index]=-PART[i].PND;
			ind[index]=n;
			index++;
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				if(m<pn)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					    
					//double w=15*r2*r2*r2/(dis*dis*dis*dis*dis);
					double w=3*r2/(dis*dis*dis);				//R/dis-1
					val[index]=w*real_lamda;
					ind[index]=m;
					index++;
					W+=w;
				}
				else if(PART[j].surface==ON)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					//double w=15*r2*r2*r2/(dis*dis*dis*dis*dis);
					double w=3*r2/(dis*dis*dis);				//R/dis-1
					W+=w; //������w���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^ 
				}
			}
			val[kk]=-W*real_lamda;
		}
		ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
	}
	////////////////////*/	 

	//�����ō����s��͑Ίp���������ׂĕ��Ȃ̂ŁA����𐳂ɂ��邽�߂ɁA�W���s��Ɖ��s���-1��������
	for(int n=0;n<pn;n++)
	{
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			val[j]*=-1;
			if(ind[j]==n && val[j]<0) cout<<"�Ίp�������� n="<<n<<endl;
		}
		B[n]*=-1;
	}
        
	/////////////////////////////////////CG�@
    double *r=new double[pn];
	double *X=new double[pn];
	
	double *AP = new double [pn];
	double *P = new double [pn];

	/////////////////////////�����l//////////////////
    if(CON->get_initialP()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initialP()==ON)//�����l�Ƃ��Č��݂̏���^����B�����ɑ����Ȃ�
	{
		for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	//CG�@�ɂ��s�������
	if(CON->get_solution()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_CGep());//CG�@
	else if(CON->get_solution()==1)
	{
		
		for(int n=0;n<pn;n++)//ICCG�@�łƂ��ꍇ�A�W���s���������ƕ��ѕς����Ȃ��Ă͂Ȃ�Ȃ�
		{      
			int num=ptr[n+1]-ptr[n];
			for(int j=ptr[n]+1;j<ptr[n+1];j++)
			{
				for(int m=ptr[n];m<j;m++)
				{
					if(ind[j]<ind[m])
					{
						double temp=val[m];
						int tempR=ind[m];
						val[m]=val[j];
						ind[m]=ind[j];
						val[j]=temp;
						ind[j]=tempR;
					}
				}
			}
		}
		iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_CGep(),&count);
	}
	
	if(CON->get_negativeP()==OFF)//�������l�����Ȃ��ꍇ
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			PART[i].P=X[n];
			if(PART[i].P<0) PART[i].P=0;
		}
	}
	else if(CON->get_negativeP()==ON)//�������l������ꍇ
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			//PART[i].P=X[n];
			if(2*X[n]==X[n]+X[n])PART[i].P=X[n];
			
		}
	}
	
		
    delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

    //////////////////////////////*/

	//CG_GPU(CON,PART,val,ind,ptr,pn,number,ppn,B,&count);

    delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] density;

	cout<<"������:"<<count<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//���͉�͂̂��߂�N0��PART[i].PND2��kernel2()�p�ɂ���������
void set_N0_and_PND2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *N0,int out)
{
	//0=���ȏ� 1=���x���U 2=0+1 3=���x���U2 4=3+PND 5=(ni-nk)+(nk-n0)

	int SW=CON->get_B_of_P();			//PPE���ɂ�������s�� 

	////////////////N0�̌v�Z
	int size = (int)(CON->get_re2()+1);	//�v�Z�̈�
	double dis2;						//����
	int d=CON->get_dimention();			//��͎���
	double R2=CON->get_re2()*CON->get_distancebp();
	double N=0;
	
	if(d==2)
	{
		if(CON->get_model_set_way()==0)//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis2=sqrt((double)(i*i+j*j));
					if(dis2!=0 && dis2<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis2,d);		
				}
			}
		}
		if(CON->get_model_set_way()==1)//�����z�u�Ƃ��čז����Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
					double dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis,d);			
				}
			}
		}
	}
	if(d==3)
	{
		if(CON->get_model_set_way()==0)//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis2=sqrt((double)(i*i+j*j+k*k));
						if(dis2!=0 && dis2<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis2,d);
					}			
				}
			}
		}
		if(CON->get_model_set_way()==1)//�����z�u�Ƃ��čז����Ƃ����ꍇ
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
						double dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis,d);
					}
				}
			}
		}
	}
	//cout<<"N0="<<N<<endl;
	*N0=N;
	///////////////////*/

	if(SW==0 || SW==2 || SW==4 || SW==5)//���q�����x���g�p����ꍇ��
	{
		///PART[i].PND2�̌v�Z���Ȃ���
		for(int i=0;i<particle_number;i++)//outwall��PART[i].N2��0�Ȃ̂ŁAPND2��0�ɂȂ邪�A���͂Ȃ�
		{
			PART[i].PND2=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				    
				double w=kernel2(R2,dis,d);
				PART[i].PND2+=w;
			}
		}
	}///////////*/
}

//�f�B���N���l�Z�b�g�֐�
void set_Dirichlet_P(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P)
{
	//file���J��
	ifstream fin;
	/*/if(flag==1) fin.open("surface_tension_P.dat", ios::in);
	if(flag==2)
	{fin.open("electromagnetic_P.dat", ios::in);
		
		if(!fin) cout<<"cannot open the file for dir_for_P"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
	}
	///////////////////////////////*/
	

	double *val=new double[fluid_number];

	int *err=new int[fluid_number];//err[i]=ON�Ȃ�A���̕\�ʗ��q�͑Ή�����f�B���N���l���[���Ƃ������ƁB��ɓd���͂�FEM_interval()����1�񂵂����s����Ȃ����ƂɋN������G���[�B
    int modify_sw=OFF;		//modify_sw=ON�Ȃ�G���[�����B�Ώ�����B

	if(flag==1)
	{
		for(int i=0;i<fluid_number;i++)
		{
			Dirichlet_P[i]+=PART[i].dir_Pst;	//�\�ʒ��͍�
			val[i]=PART[i].dir_Pst;
			//if(PART[i].surface==1  && abs(val[i])<0.1) cout<<"error tension�L��ŕ\�ʂȂ̂Ƀf�B���N���l0? cood"<<endl;
		}
	}
	
	if(flag==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			err[i]=OFF;
			val[i]=PART[i].dir_Pem;

			//fin>>val[i];
			Dirichlet_P[i]+=val[i];//̧�ق����ިظڒl�ǂݎ��	
			if(flag==2)
			{
				if(PART[i].type==BOFLUID && val[i]==0)
				{
					//cout<<"val["<<i<<"]=0"<<endl;
					err[i]=ON;
					modify_sw=ON;//�C���X�C�b�`ON
				}
			}
		}
	}

	//if(flag==2) fin.close();

	//�G���[�C��
	if(modify_sw==ON)
	{
		cout<<endl<<"set_Dirichlet_P()�ɂăG���[�����m�B�d�ݕ␳���s�BCON.FEM_interval�̏k����v���B"<<endl;
		double le=CON->get_distancebp();
		double R=CON->get_re()*le;//�C���Ɏg�p����e�����a
		ofstream fp("before_val.dat");
		if(CON->get_dimention()==2) {for(int i=0;i<fluid_number;i++) if(PART[i].type==BOFLUID) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<val[i]<<endl;}
		else if(CON->get_dimention()==3) {for(int i=0;i<fluid_number;i++) if(PART[i].type==BOFLUID) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<val[i]<<endl;}
		fp.close();

		ofstream fq("after_val.dat");
		for(int i=0;i<fluid_number;i++)
		{
			if(err[i]==ON)
			{
				double W=0;//�d�݂̑��a
				double value=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID&& PART[i].surface==ON && err[j]==OFF)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(R,dis);
						value+=val[j]*w;
						W+=w;
					}
				}
				if(W>0) value/=W;
				val[i]=value;
				Dirichlet_P[i]+=val[i];//�G���[�������Ƃ������Ƃ́A�C���O��val[i]=0���Ӗ�����B�Ⴄ�Ȃ炱�̍s�͂�����l�����邱�ƁB
			}
		}
		if(CON->get_dimention()==2) {for(int i=0;i<fluid_number;i++) if(PART[i].type==FLUID&& PART[i].surface) fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<val[i]<<endl;}
		else if(CON->get_dimention()==3) {for(int i=0;i<fluid_number;i++) if(PART[i].type==FLUID&& PART[i].surface) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fq<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<val[i]<<endl;}
		fq.close();
	}


	delete [] val;
	delete [] err;
}

//�f�B���N���l�\�ʗ͏o�͊֐�
void output_dirichlet_vector_files(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P)
{
	//�\�ʗ̓x�N�g�����o�͂���B�P�ʂ�N/m^2 flag=1�Ȃ�\�ʒ��́Aflag=2�Ȃ�d����
	
	double le=CON->get_distancebp();
	int dimention=CON->get_dimention();
	double xp=-100;	//������o��X���W
	double yp=-100;	//������o��Y���W
	double maxP=0;	//��Βl�̍ő�f�B���N���l

	//�@���x�N�g��
	double *direct[DIMENTION];
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];//�������@���x�N�g��
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	}

	//maxP���߂�
	for(int i=0;i<fluid_number;i++) if(fabs(Dirichlet_P[i])>maxP) maxP=fabs(Dirichlet_P[i]);

	double times=4*le/maxP;
	if(CON->get_dir_for_P()==3) times*=0.5;//�����̏ꍇ�͔{������
	

	ofstream gg2("dirichlet_vector.dat");
	if(dimention==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			gg2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]*direct[A_X][i]*times<<" "<<Dirichlet_P[i]*direct[A_Y][i]*times<<endl;
			if(PART[i].r[A_X]>xp) xp=PART[i].r[A_X];
			if(PART[i].r[A_Z]>yp) yp=PART[i].r[A_Y];
		}
		xp+=4*le;
		yp+=4*le;
		gg2<<xp<<"\t"<<yp<<"\t"<<maxP*times<<" "<<0<<" ////�}��̒�����"<<maxP<<endl;
	}
	if(dimention==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le)
			{
				gg2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]*direct[A_X][i]*times<<" "<<Dirichlet_P[i]*direct[A_Z][i]*times<<endl;
				if(PART[i].r[A_X]>xp) xp=PART[i].r[A_X];
				if(PART[i].r[A_Z]>yp) yp=PART[i].r[A_Z];
			}
		}
		xp+=4*le;
		yp+=4*le;
		gg2<<xp<<"\t"<<yp<<"\t"<<maxP*times<<" "<<0<<" ////�}��̒�����"<<maxP<<endl;
	}
	gg2.close();

	for(int D=0;D<DIMENTION;D++) delete [] direct[D]; 
}

///���q�����x�̎��Ԕ����v�Z�֐�
double Dndt(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0)
{
	///�����ŕԂ����Dndt�́A�d�݊֐�pow(R,D)/pow(r,D)��p����MPS�ɂ�葬�x�̔��U���v�Z�����l�ɓ�����
    double R=CON->get_distancebp()*CON->get_re2();	//�e�����a
    double Dndt=0;									//���֐��̒l
	double D=2;										//��͎����B�������v�Z�̓s����Adouble�^
	if(CON->get_dimention()==3) D=3;
    
    for(int k=0;k<PART[i].N2;k++)
    {    
        int j=PART[i].NEI2[k]; 
        double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);
		
		double div=(PART[j].u[A_X]-PART[i].u[A_X])*X+(PART[j].u[A_Y]-PART[i].u[A_Y])*Y+(PART[j].u[A_Z]-PART[i].u[A_Z])*Z;
		Dndt-=D*pow(R,D)/pow(dis,D)*div/(dis*dis);	
    }
    return Dndt;
}

///���̓v���b�g�֐��i�X�J���[�j
void plot_P(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int t,int fluid_number)
{
	ofstream fout("P.dat");
    double le=CON->get_distancebp();

    if(CON->get_dimention()==2)
    {
		for(int i=0;i<particle_number;i++) fout<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].P<<endl;
    }
    else if(CON->get_dimention()==3)
    {
        for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fout<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].P<<endl;
    }
    fout.close();

	ofstream fout2("Pf.dat");
    le=CON->get_distancebp();

    if(CON->get_dimention()==2)
    {
		for(int i=0;i<fluid_number;i++) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].P<<endl;
    }
    else if(CON->get_dimention()==3)
    {
        for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].P<<endl;
    }
    fout.close();

	if(CON->get_P_AVS()>0)//������if(CON->get_P_AVS()>0 &&  t%CON->get_P_AVS())�Ƃ����0�̂Ƃ��o�O��
	{
		if(t==1 || t%CON->get_P_AVS()==0) output_pressuer_avs(CON,PART,t,particle_number,fluid_number);
	}
}

//����AVS�t�@�C���o�͊֐�
void output_pressuer_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	t=1;//���܂͂킴�Ɩ��X�e�b�v�㏑��

	//sprintf_s(filename,"pressure/pressure%d",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"pressure%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	ofstream fout(filename);
	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=PART[i].P;
				fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n++;
			}
		}
	}
	else if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			double x=PART[i].r[A_X]*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
			double y=PART[i].r[A_Y]*1.0E+05;
			double z=PART[i].r[A_Z]*1.0E+05;
			double P=PART[i].P;
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
	}
	fout.close();
	//sprintf_s(filename,"pressure/pressure%d.fld",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"pressure%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
	fout2 << "label=pressure" << endl << endl;
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	fout2 << "variable 1 file=pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    1 file=pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    2 file=pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    3 file=pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2.close();
}

///���q�f�[�^�ǂݎ��֐�
void calc_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double **reU,double n0,double dt,int minPsw)
{
	if(minPsw==ON) cout<<"���͌��z�v�Z�J�n ver."<<CON->get_Pgrad()<<" minP=ON ---------";
	if(minPsw==OFF) cout<<"���͌��z�v�Z�J�n ver."<<CON->get_Pgrad()<<" minP=OFF ---------";
	double le=CON->get_distancebp();	//�������q�ԋ���
	unsigned int timeA=GetTickCount();	//�v�Z�J�n����

	double *direct[DIMENTION];			//�������@���޸�يi�[
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];

	double *Pgrad[DIMENTION];			//���͌��z�޸�يi�[
	for(int D=0;D<DIMENTION;D++) Pgrad[D]=new double [fluid_number];
	
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			direct[D][i]=0.0;
			Pgrad[D][i]=0.0;
		}
	}

	//////�@���޸�ٌv�Z�J�n
	if(CON->get_Pgrad()==2 ||CON->get_Pgrad()==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON ) direct_f(CON,PART,i,direct);
			else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
		}
		
	}/////////////////*/
	
	///////////////////���͌��z�v�Z
	
	if(CON->get_Pgrad()==3)//�\�ʂ̂ݖ@�� && minP=0�̂Ƃ��\�ʗ��q����������͌v�Z
	{
		P_gradient3(CON,PART,fluid_number,direct,dt,reU,Pgrad,minPsw);
	}
	else if(CON->get_Pgrad()==4)//�d�݂��ŏ����@(WLSM)
	{
		P_gradient4(CON,PART,dt,fluid_number,reU,Pgrad,minPsw);
	}
	else if(CON->get_Pgrad()==5)//�d�݂��ŏ����@(WLSM)
	{
		P_gradient5(CON,PART,dt,fluid_number,reU,Pgrad);
	}
	///////////////////////*/

	//���x�C���ʌv�Z
	double limit_reU=CON->get_Pgrad_limit()*CON->get_distancebp();
	double *density=new double[fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*Pgrad[D][i]/density[i];
	delete [] density;

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] Pgrad[D];

	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

///���͌��z�v�Z�֐�ver.3
void P_gradient3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *direct[DIMENTION],double dt,double **reU,double **Pgrad,int minPsw)
{
	double le=CON->get_distancebp();//�������q�ԋ���
	double r=CON->get_re()*le;
	int d=CON->get_dimention();
	
	double *minP=new double[fluid_number];	//���͂̍ŏ����͊i�[

	///minP�����߂�
	set_minP(CON,PART,fluid_number,minP,minPsw);
	///////////////

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Pgrad[D][i]=0;//������

		double W=0;//���q�����x�@OUT���������肷��̂�PND[i]�͔���

		if(PART[i].surface==ON)//�\�ʗ��q�̏ꍇ�A���������̈��͂𒲂ׂČ��z���v�Z����B
		{
			double x1=PART[i].r[A_X]+direct[A_X][i]*le;//���������̍��W
			double y1=PART[i].r[A_Y]+direct[A_Y][i]*le;
			double z1=PART[i].r[A_Z]+direct[A_Z][i]*le;
			double P=0;//���������̈���
			
			for(int k=0;k<PART[i].N3;k++)
			{       
	 			int j=PART[i].NEI3[k];
				if(PART[j].type!=OUTWALL)//BDWALL���Ȃ����ق����悭�Ȃ��H
				{
					double X=PART[j].r[A_X]-x1;
					double Y=PART[j].r[A_Y]-y1;
					double Z=PART[j].r[A_Z]-z1;
					double dis=sqrt(X*X+Y*Y+Z*Z);
					if(dis<r)
					{
						double w=(1-dis/r)*(1-dis/r);//w(0)=���ƂȂ�֐��͎g���Ȃ�
						W+=w;
						P+=PART[j].P*w;			
					}
				}
			}

			if(W!=0) P/=W; //(x1,y1,z1)�ł̈���
			Pgrad[A_X][i]=(P-PART[i].P)*direct[A_X][i]/le;
			Pgrad[A_Y][i]=(P-PART[i].P)*direct[A_Y][i]/le;
			Pgrad[A_Z][i]=(P-PART[i].P)*direct[A_Z][i]/le;
		}
		else//�������̗��q�̏ꍇ�A���͂��爳�͌��z���v�Z
		{
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					
					double w=kernel(r,dis);
					W+=w;
					if(minP[i]==0 && PART[j].type==FLUID && PART[j].surface==ON)///���qj�̈��͂𗱎qi�̈��͂ł���������
					{   
						Pgrad[A_X][i]+=(PART[i].P-minP[i])*X*w/(dis*dis);
						Pgrad[A_Y][i]+=(PART[i].P-minP[i])*Y*w/(dis*dis);
						Pgrad[A_Z][i]+=(PART[i].P-minP[i])*Z*w/(dis*dis);
					}
					else
					{   //�ʏ�ǂ���
						Pgrad[A_X][i]+=(PART[j].P-minP[i])*X*w/(dis*dis);
						Pgrad[A_Y][i]+=(PART[j].P-minP[i])*Y*w/(dis*dis);
						Pgrad[A_Z][i]+=(PART[j].P-minP[i])*Z*w/(dis*dis);
					}
				}
			}
			
			for(int D=0;D<DIMENTION;D++) if(W!=0) Pgrad[D][i]= Pgrad[D][i]*d/W;
		}
	}///////////////Pgrad[D][i]�v�Z�I��

	///////////////////////////////���͌��z��\��
	if(CON->get_iteration_count()==1)				//�A��͂𕡐���s���ꍇ�́A����o�͂��Ă��܂��̂ł�߂�
	{
		plot_Pgradient(CON,PART,fluid_number,Pgrad);
	}
	////////////////////////*/

	delete [] minP;
}

//�ŏ����͌v�Z�֐�
void set_minP(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *minP,int minPsw)
{
	///minP�����߂�
	for(int i=0;i<fluid_number;i++) minP[i]=PART[i].P;
	if(minPsw==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					if(PART[j].P<minP[i]) minP[i]=PART[j].P;
				}
			}
		}
	}////////////////////*/
}

///���͌��z�\���֐�
void plot_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Pgrad)
{
	int d=CON->get_dimention();
	double times=0;					//�\���p�{��
	double le=CON->get_distancebp();
	double density=CON->get_density();

	//�ʏ�ʂ�
	times=CON->get_times()*le*le/density*CON->get_Pgrad_times();
	
	ofstream gra("Pgrad.dat");
	ofstream gra2("Pgrad2.dat");//�X�J���[�\��
	if(d==2)//���₷���悤�ɕ����𔽓]�\��
	{
		for(int i=0;i<fluid_number;i++)
		{
			gra<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<-Pgrad[A_X][i]*times<<" "<<-Pgrad[A_Y][i]*times<<endl;
			gra2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(Pgrad[A_X][i]*Pgrad[A_X][i]+Pgrad[A_Y][i]*Pgrad[A_Y][i])<<endl;
		}
	}
	else if(d==3)//���₷���悤�ɕ����𔽓]�\��
	{
		for(int i=0;i<fluid_number;i++) 
		{
			if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)
			{
				gra<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<-Pgrad[A_X][i]*times<<" "<<-Pgrad[A_Z][i]*times<<endl;
				gra2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(Pgrad[A_X][i]*Pgrad[A_X][i]+Pgrad[A_Y][i]*Pgrad[A_Y][i]+Pgrad[A_Z][i]*Pgrad[A_Z][i])<<endl;
			}
		}
	}
	gra.close();
	gra2.close();
}

///���͌��z�v�Z�֐�ver.4 �d�݂��ŏ����@�Ɋ�Â����z(WLSM)
void P_gradient4(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3],double *P_grad[3],int minPsw)
{
	//WLSM�Ɋւ��钍�ӓ_
	//�EminP=OFF�̂ق������܂�����??
	//�E2D��͂�2�����x�̏ꍇ�Are=2.1���Ɠ��ɕ\�ʗ��q�̌v�Z�ɂ����ċߗח��q�������Ȃ��A�v�Z���ʂ����������Ƃ�������Bre=2.5���炢�𐄏�??�B
	//�E���Ɠ������R�ŁA�Ǘ����q�₻��ɋ߂����q�̌v�Z�͔j�]���邨���ꂠ��B

	double le=CON->get_distancebp();
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;							//�W���s��̌�
	int order=CON->get_Pgrad_order();	//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��

	double *minP=new double[fluid_number];	//���͂̍ŏ����͊i�[

	///minP�����߂�
	set_minP(CON,PART,fluid_number,minP,minPsw);
	///////////////

	//�W���s��̑傫���̌���
	if(d==2)
	{
		if(order==1) N=2;
		else if(order==2) N=5;
	}
	else if(d==3)
	{
		if(order==1) N=3;
		else if(order==2) N=9;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *B=new double [N];	//N�̉��s��
	
	if(d==2 && order==1)//�񎟌�
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//������
			for(int n=0;n<N;n++) B[n]=0;

			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X])/le;// le�Ŋ���̂͑ł��؂�덷�h�~
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/le;
					double dis=sqrt(X*X+Y*Y);
					
					double w=1;
					if(dis>1) w=1/(dis*dis*dis*dis);
					
					matrix[0]+=X*X*w;			//��Xjwj
					matrix[1]+=X*Y*w;		//��XjYjwj
					matrix[3]+=Y*Y*w;			//��Yjwj
				
					B[0]+=(PART[j].P-minP[i])*X*w;//��fjXjwj
					B[1]+=(PART[j].P-minP[i])*Y*w;//��fjYjwj
				}
			}
			matrix[2]=matrix[1];		//��XjYjwj
			for(int n=0;n<N;n++) B[n]/=le;//�ł��؂�덷�h�~

			double determinant=matrix[0]*matrix[3]-matrix[1]*matrix[2];//�s��
			
			P_grad[A_X][i]=(B[0]*matrix[3]-matrix[1]*B[1])/determinant;
			P_grad[A_Y][i]=(B[1]*matrix[0]-matrix[2]*B[0])/determinant;

			int flag=OFF;//OFF�Ȃ���Ȃ��BON�Ȃ�MPS�ł�蒼��
			for(int D=0;D<2;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//������Ȃ�ON
			if(flag==ON)
			{
				//cout<<"�G���[�̂���MPS�ɂ�藣�U�� "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//�����ł͖��Ӗ��ȕϐ��Ȃ̂œK���ɂ���Ƃ�
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			
		}
	}
	else if(d==2 && order==2)//�񎟌�2����
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;

			for(int n=0;n<N*N;n++) matrix[n]=0;//������
			for(int n=0;n<N;n++) B[n]=0;
		
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(PART[j].type!=OUTWALL)
				if(PART[j].type==FLUID)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X]);
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
					double dis=sqrt(X*X+Y*Y);
					double dP=PART[j].P-minP[i];
						
					double w=1;
					
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
						
					matrix[0]+=X*X*w;			//��Xjwj
					matrix[1]+=X*Y*w;		//��XjYjwj
					matrix[2]+=X*X*X*w;			//��Xj^3wj
					matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
					matrix[4]+=X*Y*Y*w;			//��XjYj^2wj
	
					matrix[6]+=Y*Y*w;			//��Yj^2wj
					matrix[9]+=Y*Y*Y*w;			//��Yj^3wj
	
					matrix[12]+=X*X*X*X*w;			//��Xj^4wj
					matrix[13]+=X*X*X*Y*w;			//��Xj^3Yjwj
					matrix[14]+=X*X*Y*Y*w;			//��Xj^2Yj^2wj
		
					matrix[19]+=X*Y*Y*Y*w;			//��XjYj^3wj
	
					matrix[24]+=Y*Y*Y*Y*w;			//��Yj^4wj

				
					B[0]+=dP*X*w;//��dPjXjwj
					B[1]+=dP*Y*w;//��dPjYjwj
					B[2]+=dP*X*X*w;//��dPjXj^2wj
					B[3]+=dP*X*Y*w;//��dPjXjYjwj
					B[4]+=dP*Y*Y*w;//��dPjYj^2wj
					//if(i==551) cout<<endl<<dP*Y*w<<" "<<dP<<" "<<Y<<" "<<w<<endl;
				}
			}
			//if(i==551) cout<<"B[1]="<<B[1]<<endl;
			
			matrix[5]=matrix[1];		//��XjYjwj
			matrix[7]=matrix[3];		//��Xj^2Yjwj
			matrix[8]=matrix[4];		//��XjYj^2wj
			matrix[10]=matrix[2];		//��Xj^3Yjwj
			matrix[11]=matrix[3];		//��Xj^2Yjwj
			matrix[15]=matrix[3];		//��Xj^2Yjwj
			matrix[16]=matrix[4];		//��XjYj^2wj
			matrix[17]=matrix[13];		//��Xj^3Yjwj
			matrix[18]=matrix[14];		//��Xj^2Yj^2wj
			matrix[20]=matrix[4];		//��XjYj^2wj
			matrix[21]=matrix[9];		//��Yj^3wj
			matrix[22]=matrix[14];		//��Xj^2Yj^2wj
			matrix[23]=matrix[19];		//��XjYj^3wj

			/*/�ۂߌ덷�h�~
			for(int n=0;n<N;n++)
			{
				for(int m=0;m<N;m++) matrix[n*N+m]*=1e14;
				B[n]*=1e14;
			}//*/

			double dPdx=0;
			double dPdy=0;
			
			return_X_for5N(matrix,N,B,B,&dPdx,&dPdy);//5���A���������́A���P�C�Q��Ԃ��֐�(�����ł�B���Q�n��)
		
			P_grad[A_X][i]=dPdx;
			P_grad[A_Y][i]=dPdy;
		}
			
	}
	else if(d==3 && order==1)//3����1����
	{
		///�W���s���
		///   ����x2    ����x��y  ����x��z  a = ����x��f  
		///  ����x��y    ����y2   ����y��z  b = ����y��f 
		///  ����x��z   ����y��z   ����z2   c = ����z��f 
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//������
			for(int n=0;n<N;n++) B[n]=0;

			//if(PART[i].N>4)//���ӗ��q�����ɒ[�ɏ��Ȃ��ƍs�񎮂��[���ɂȂ��ĉ����Ȃ�
			if(PART[i].N>=3)
			{
				cacl_WLSM_P_D3_order1(CON,PART,matrix,B,P_grad, i,N,minP);//1���ߎ� �{����3�̋ߗח��q�Ōv�Z�\�����ǁA�ی��������Ă����ł�5�ȏ�̌v�Z�_���K�v�Ƃ��Ă���B
				int flag=OFF;//OFF�Ȃ���Ȃ��BON�Ȃ�MPS�ł�蒼��
				for(int D=0;D<DIMENTION;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//������Ȃ�ON
				if(flag==ON)
				{
					//cout<<"�G���[�̂���MPS�ɂ�藣�U�� "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					double grad[3]={0,0,0};
					double n0=1;//�����ł͖��Ӗ��ȕϐ��Ȃ̂œK���ɂ���Ƃ�
					P_gradient_MPS(CON,PART,i,n0,grad);//
					for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
				}
			}
			else//���ӗ��q�̏d�ݕ��ς�. ���q�������Ȃ��̂ŁAMPS�͂���������
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"�x�� �ߗח��q����4�ȉ�("<<PART[i].N<<")�ł�"<<endl;
			}
		}
	}
	else if(d==3 && order==2)//3����2����
	{
		//P=Pi+a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z     a = ����x��P  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z   b = ����y��P
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2     c = ����z��P
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z     d = ����x2��P
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z  e = ����y2��P
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3     f = ����z2��P
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z  g = ����x��y��P
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2  h = ����y��z��P
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2    g = ����x��z��P

		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//������
			for(int n=0;n<N;n++) B[n]=0;

			if(PART[i].N>8) cacl_WLSM_P_D3_order2(CON,PART,matrix,B,P_grad, i,N,minP);//2���ߎ� �{����6�̋ߗח��q�Ōv�Z�\�����ǁA�ی��������Ă����ł�9�ȏ�̌v�Z�_���K�v�Ƃ��Ă���B
			else if(PART[i].N>4) cacl_WLSM_P_D3_order1(CON,PART,matrix,B,P_grad, i,3,minP);//1���ߎ�
			int flag=OFF;//OFF�Ȃ���Ȃ��BON�Ȃ�MPS�ł�蒼��
			for(int D=0;D<DIMENTION;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//������Ȃ�ON
			if(flag==ON)
			{
				cout<<"�G���[�̂���MPS�ɂ�藣�U�� "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//�����ł͖��Ӗ��ȕϐ��Ȃ̂œK���ɂ���Ƃ�
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			if(PART[i].N<=4)
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"�x�� �ߗח��q����4�ȉ�("<<PART[i].N<<")�ł�"<<endl;
				//for(int D=0;D<3;D++) P_grad[D][i]=0;
			}
		}
	}

	/*/���x�C���ʌv�Z
	double *density=new double[fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*P_grad[D][i]/density[i];
	delete [] density;*/

	///////////////////////////////���͌��z��\��
	if(CON->get_iteration_count()==1)				//�A��͂𕡐���s���ꍇ�́A����o�͂��Ă��܂��̂ł�߂�
	{
		plot_Pgradient(CON,PART,fluid_number,P_grad);
	}
	////////////////////////*/

	delete [] matrix;
	delete [] B;
	delete [] minP;
}

///���͌��z�v�Z�֐�ver.5 ���g��ʂ�Ɖ��肵�Ȃ�WLSM
void P_gradient5(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3],double *P_grad[3])
{
	//P_gradient4()��薢�m�����P�����ꍇ�@�ߎ��ǖʂ����g�̒l��ʂ�Ɖ��肵�Ȃ��B
	//���̓s����AminP�͍l�����Ȃ��B

	double le=CON->get_distancebp();
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N0=0;							//�W���s��̌�
	int order0=CON->get_Pgrad_order();	//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
	double threshold=CON->get_threshold();	//���͌덷��臒l(�A�_�v�e�B�u������ۂɎg�p)

	//�W���s��̑傫���̌���
	if(d==2)
	{
		if(order0==1) N0=3;
		else if(order0==2) N0=6;
		else if(order0==3) N0=10;
	}
	else if(d==3)
	{
		if(order0==1) N0=4;
		else if(order0==2) N0=10;
		else if(order0==3) N0=20;
	}
	////////////////////////////////

	double *matrix=new double [N0*N0];	//N�~N�̌W���s��
	double *B=new double [N0];	//N�̉��s��
	double **MAT= new double*[N0];		//N�~N�̌W���s��(�z���2����)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//���x�N�g���i�[
	ofstream fg("error_in_P.dat");
	
	int N=N0;
	int order=order0;
	if(d==2 )//�񎟌�
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			int jnb=0;		//���ӗ��q���B������OUTWALL�Ȃǂ͏���
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL) jnb++;
			}

			order=order0;
			if(order==3) {if(jnb<10) {order=2; N=6;}}	//���ӗ��q�������Ȃ�������2���ߎ��ɂ���
			if(order==2) {if(jnb<6) {order=1; N=3;}}	//���ӗ��q�������Ȃ�������1���ߎ��ɂ���

			if(jnb>=2)		//���ӗ��q���������菭�Ȃ����́A���͌��z�̓[��
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;//������
				for(int n=0;n<N0;n++) B[n]=0;
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//������

				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type!=OUTWALL)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// le�Ŋ���̂͑ł��؂�덷�h�~
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
						double dP=(PART[j].P)/PART[i].L;
						double dis=sqrt(X*X+Y*Y);
					
						double w=1;
						if(dis>1) w=1/(dis*dis*dis*dis);
						//if(dis>1) w=kernel_in_WLSM( dis, CON->get_re());

						if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//���x�N�g��
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}

						//�s��쐬
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B[n]+=base[n]*w*dP;
						}
					}
				}

				MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
				B[N-1]+=PART[i].P/PART[i].L;	//���qi�̊�^�B
				//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�
				
				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��

				gauss(matrix,B,N);//�K�E�X�̏����@�ŉ���
				
				P_grad[A_X][i]=B[0];	//1���ߎ��ł�2���ߎ��ł��A���m���̏��ԓI�ɂ����Ȃ�悤�ɂ��Ă���
				P_grad[A_Y][i]=B[1];
				
				int flag=OFF;//OFF�Ȃ���Ȃ��BON�Ȃ�MPS�ł�蒼��
				for(int D=0;D<2;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//������Ȃ�ON
				if(flag==ON)
				{
					cout<<"�G���[�̂���MPS�ɂ�藣�U�� "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					double grad[3]={0,0,0};
					double n0=1;//�����ł͖��Ӗ��ȕϐ��Ȃ̂œK���ɂ���Ƃ�
					P_gradient_MPS(CON,PART,i,n0,grad);//
					for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
				}
				/*double Q=0;
				int count_for_Q=0;
				if(flag==OFF)
				{
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=(PART[j].r[A_X]-PART[i].r[A_X]);
							double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
							double dis=sqrt(X*X+Y*Y);
					
							double w=1;
							//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
							if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
							Q+=(a*X+b*Y+PART[i].P-PART[j].P)*(a*X+b*Y+PART[i].P-PART[j].P)*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B
					if(Q>threshold) PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B
				}
				//fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(Q>threshold) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;*/
			}
		}
	}
	else if(d==3 )//3����1����
	{
		///�W���s���
		///   ����x2    ����x��y  ����x��z  a = ����x��f  
		///  ����x��y    ����y2   ����y��z  b = ����y��f 
		///  ����x��z   ����y��z   ����z2   c = ����z��f 

		//�Q���ߎ��Ȃ�
		//P=Pi+a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z     a = ����x��P  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z   b = ����y��P
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2     c = ����z��P
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z     d = ����x2��P
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z  e = ����y2��P
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3     f = ����z2��P
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z  g = ����x��y��P
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2  h = ����y��z��P
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2    g = ����x��z��P

		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			int jnb=0;		//���ӗ��q���B������OUTWALL�Ȃǂ͏���
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL) jnb++;
			}

			order=order0;
			if(order==3) if(jnb<20) {order=2; N=10;}	//���ӗ��q�������Ȃ�������2���ߎ��ɂ���
			if(order==2) if(jnb<10) {order=1; N=4;}	//���ӗ��q�������Ȃ�������1���ߎ��ɂ���

			if(jnb>=3)
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;//������
				for(int n=0;n<N0;n++) B[n]=0;
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//������
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type!=OUTWALL)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double dP=(PART[j].P);//��fj�ɑ���
					
						double w=1;
						if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
						//if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
						
						if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//���x�N�g��
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
						//�s��쐬
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B[n]+=base[n]*w*dP;
						}
					}
				}
				MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
				B[N-1]+=PART[i].P;	//���qi�̊�^�B
				
				//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�*/

				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��
				
				gauss(matrix,B,N);//�K�E�X�̏����@�ŉ���
				double Px=B[0];//X��������
				double Py=B[1];//Y��������
				double Pz=B[2];//Z��������

				P_grad[A_X][i]=Px;
				P_grad[A_Y][i]=Py;
				P_grad[A_Z][i]=Pz;
				
			}
			else//���ӗ��q�̏d�ݕ��ς�. ���q�������Ȃ��̂ŁAMPS�͂���������
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"�x�� �ߗח��q����4�ȉ�("<<PART[i].N<<")�ł�"<<endl;
			}

			int flag=OFF;//OFF�Ȃ���Ȃ��BON�Ȃ�MPS�ł�蒼��
			for(int D=0;D<3;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//������Ȃ�ON
			if(flag==ON)
			{
				cout<<"�G���[�̂���MPS�ɂ�藣�U�� "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//�����ł͖��Ӗ��ȕϐ��Ȃ̂œK���ɂ���Ƃ�
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			double Q=0;
			int count_for_Q=0;
			if(order==1)
			{
				if(flag==OFF && PART[i].N>=3)
				{
					//double a,b,c,d,e;
					double a,b,c;
					a=B[0];//X��������
					b=B[1];//Y��������
					c=B[2];//Z��������
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=(PART[j].r[A_X]-PART[i].r[A_X]);
							double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							double w=1;
							//if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
							if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
							double error=a*X+b*Y+c*Z+PART[i].P-PART[j].P;
							Q+=error*error*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
				
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B
					if(Q>threshold) 
					{
						PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B
						//cout<<"���͌덷���傫�����ߕ���(threshold)"<<endl;
					}
				}
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(fabs(PART[i].r[A_Y])<le) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Q<<endl;
			}
			else if(order==2)
			{
				if(flag==OFF)
				{
					double a,b,c,d,e,f,g,h,ii;
					a=B[0];//X��������
					b=B[1];//Y��������
					c=B[2];//Z��������
					d=B[3];//X����2�K����
					e=B[4];//Y����2�K����
					f=B[5];//Z����2�K����
					g=B[6];//XY��������
					h=B[7];//YZ��������
					ii=B[8];//ZX��������
					
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							double w=1;
							if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
							double error=a*X+b*Y+c*Z+d*X*X+e*Y*Y+f*Z*Z+g*X*Y+h*Y*Z+ii*X*Z+PART[i].P-PART[j].P;
							Q+=error*error*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
				
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B	
					if(Q>threshold)
					{
						PART[i].division_flag=1;		//�덷���傫���Ɣ��f���āA��������B
						//cout<<"���͌덷���傫�����ߕ���(threshold)"<<endl;
					}
				}
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(fabs(PART[i].r[A_Y])<le) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Q<<endl;
			}
		}
	}

	///////////////////////////////���͌��z��\��
	if(CON->get_iteration_count()==1)				//�A��͂𕡐���s���ꍇ�́A����o�͂��Ă��܂��̂ł�߂�
	{
		plot_Pgradient(CON,PART,fluid_number,P_grad);
	}
	////////////////////////*/
	fg.close();
	delete [] matrix;
	delete [] B;
	for(int n=0;n<N0;n++) delete [] MAT[n];
    delete [] MAT;
	delete [] base;
}


//Pgradient4�ɂ�����A3����1���ߎ����s���֐�
void cacl_WLSM_P_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP)
{
	double le=CON->get_distancebp();
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double dP=(PART[j].P-minP[i]);//��fj�ɑ���
					
			double w=1;
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
						
			matrix[0]+=X*X*w;			//��Xj^2wj
			matrix[1]+=X*Y*w;		//��XjYjwj
			matrix[2]+=X*Z*w;		//��XjZjwj
						
			matrix[4]+=Y*Y*w;			//��Yj^2wj
			matrix[5]+=Y*Z*w;		//��YjZjwj
	
			matrix[8]+=Z*Z*w;			//��Zj^2wj
					
			B[0]+=dP*X*w;//��fjXjwj
			B[1]+=dP*Y*w;//��fjYjwj
			B[2]+=dP*Z*w;//��fjZjwj
		}
	}
	matrix[3]=matrix[1];		//��XjYjwj
	matrix[6]=matrix[2];		//��XjZjwj
	matrix[7]=matrix[5];		//��YjZjwj

	//�v�Z���邩���Ȃ����𔻒�
	if(B[0]==0 && B[1]==0 && B[2]==0)
	{
		//���s�񂪂��ׂă[���̂Ƃ��͉����Ȃ��Ȃ�̂ŃK�E�X�̏����@�ɂ�������ł͂����Ȃ��B
		for(int D=0;D<3;D++) P_grad[D][i]=0;//���s�񂪃[���Ƃ������Ƃ͌��z���[���Ƃ������Ƃ����炻�����
	}
	else
	{
		gauss(matrix,B,N);//�K�E�X�̏����@�ŉ���
		double Px=B[0];//X��������
		double Py=B[1];//Y��������
		double Pz=B[2];//Z��������

		P_grad[A_X][i]=Px;
		P_grad[A_Y][i]=Py;
		P_grad[A_Z][i]=Pz;
	}
	
	///*/

	/*double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//�s��
	//if(i==48155)cout<<"D="<<determinant<<endl;		
	P_grad[A_X][i]=(B[0]*matrix[4]*matrix[8]-B[0]*matrix[5]*matrix[7]-matrix[1]*B[1]*matrix[8]+matrix[1]*matrix[5]*B[2]+matrix[2]*B[1]*matrix[7]-matrix[2]*matrix[4]*B[2])/determinant;
	P_grad[A_Y][i]=(matrix[0]*B[1]*matrix[8]-matrix[0]*matrix[5]*B[2]-B[0]*matrix[3]*matrix[8]+B[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B[2]-matrix[2]*B[1]*matrix[6])/determinant;
	P_grad[A_Z][i]=(matrix[0]*matrix[4]*B[2]-matrix[0]*B[1]*matrix[7]-matrix[1]*matrix[3]*B[2]+matrix[1]*B[1]*matrix[6]+B[0]*matrix[3]*matrix[7]-B[0]*matrix[4]*matrix[6])/determinant;	
	///*/
}

//Pgradient4�ɂ�����A3����2���ߎ����s���֐�
void cacl_WLSM_P_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP)
{
	//P=Pi+a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z     a = ����x��P  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z   b = ����y��P
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2     c = ����z��P
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z     d = ����x2��P
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z  e = ����y2��P
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3     f = ����z2��P
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z  g = ����x��y��P
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2  h = ����y��z��P
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2    g = ����x��z��P

	double le=CON->get_distancebp();
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B
	double B_val[3];						//���s���3�v�f�������ۑ����Ă����B���Ƃ�1���ߎ��ł��Ȃ����ꍇ�ɕK�v
	double matrix_val[9];					//�W���s���9�v�f�������ۑ����Ă����B���Ƃ�1���ߎ��ł��Ȃ����ꍇ�ɕK�v

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double dP=(PART[j].P-minP[i]);//��fj�ɑ���
						
			double w=1;
			//if(dis>le) w=le*le/(dis*dis);
			//if(dis>le) w=le/(dis);
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
			weight[k]=w;

			matrix[0]+=X*X*w;		//��Xj^2wj
			matrix[1]+=X*Y*w;		//��XjYjwj
			matrix[2]+=X*Z*w;		//��XjZjwj
			matrix[3]+=X*X*X*w;		//��Xj^3wj
			matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
			matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
			matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
			matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
			matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
	
			matrix[10]+=Y*Y*w;		
			matrix[11]+=Y*Z*w;		
			matrix[12]+=X*X*Y*w;		
			matrix[13]+=Y*Y*Y*w;		
			matrix[14]+=Y*Z*Z*w;
			matrix[15]+=X*Y*Y*w;
			matrix[16]+=Y*Y*Z*w;
						
			matrix[20]+=Z*Z*w;			
			matrix[23]+=Z*Z*Z*w;		
					
			matrix[30]+=X*X*X*X*w;
			matrix[31]+=X*X*Y*Y*w;
			matrix[32]+=X*X*Z*Z*w;	
			matrix[33]+=X*X*X*Y*w;	
			matrix[34]+=X*X*Y*Z*w;	
			matrix[35]+=X*X*X*Z*w;	
					
			matrix[40]+=Y*Y*Y*Y*w;
			matrix[41]+=Y*Y*Z*Z*w;
			matrix[42]+=X*Y*Y*Y*w;
			matrix[43]+=Y*Y*Y*Z*w;
			matrix[44]+=X*Y*Y*Z*w;

			matrix[50]+=Z*Z*Z*Z*w;	//6�s��
			matrix[51]+=X*Y*Z*Z*w;
			matrix[52]+=Y*Z*Z*Z*w;
			matrix[53]+=X*Z*Z*Z*w;

			//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\


			B[0]+=dP*X*w;		//a
			B[1]+=dP*Y*w;		//b
			B[2]+=dP*Z*w;		//c
			B[3]+=dP*X*X*w;		//d
			B[4]+=dP*Y*Y*w;		//e
			B[5]+=dP*Z*Z*w;		//f
			B[6]+=dP*X*Y*w;		//g
			B[7]+=dP*Y*Z*w;		//h
			B[8]+=dP*X*Z*w;		//i
		}
	}
	matrix[9]=matrix[1];		//��XjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];
		
	for(int L=0;L<3;L++) B_val[L]=B[L];//���s��̒l��ۑ�
	for(int L=0;L<3;L++)
	{
		for(int M=0;M<3;M++)
		{
			matrix_val[L*3+M]=matrix[L*9+M];//���s��̒l��ۑ�
		}
	}
		

	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	gauss(matrix,B,N);

	//�덷�𒲍�
	double Q=0;
	double W=0;//�d�݂̑��a
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dP=(PART[j].P-minP[i]);//��fj�ɑ���
			double w=weight[k];
			W+=w;
			double err=B[0]*X+B[1]*Y+B[2]*Z+B[3]*X*X+B[4]*Y*Y+B[5]*Z*Z+B[6]*X*Y+B[7]*Y*Z+B[8]*X*Z-dP;
			Q+=err*err*w;
		}
	}
	if(W>0) Q/=W;

	if(Q>100)//�덷��100�ȏ�Ȃ�(��������������)
	{
		//cout<<i<<endl;
		for(int n=0;n<3;n++) B[n]=B_val[n];	//Px=B[0]�̂悤�ɁAB[n]�ɓ������i�[������悤�ɂ��Ă�̂ŁA�����ł�B���g��
		N=3;								//N=3�ɂ��Ă����Ȃ��Ƃ����Ȃ�
		gauss(matrix_val,B,N);				//matrix_val�ɂ͍���g�p�����s��̂���[3][3]�������P�����Ɋi�[����Ă���B
	}///*/
	
	double Px=B[0];//X��������
	double Py=B[1];//Y��������
	double Pz=B[2];//Z��������
	
	P_grad[A_X][i]=Px;
	P_grad[A_Y][i]=Py;
	P_grad[A_Z][i]=Pz;/////////*/

	delete [] weight;
}

///���͌��z�v�Z�֐�
void P_gradient_MPS(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0,double *P_grad)
{
    double R=CON->get_re()*CON->get_distancebp();
	double W=0;//���q�����x�@OUT���������肷��̂�PND[i]�͔���
	double minP=PART[i].P;//���͂̍ŏ�����
	int d=CON->get_dimention();
	int J=i;
	for(int k=0;k<PART[i].N;k++)
	{       
	    int j=PART[i].NEI[k];
	    if(PART[j].type!=OUTWALL)
	    {
	        if(PART[j].P<minP)
			{
				minP=PART[j].P;
				J=j;
			}
	    }
	}///minP[i]�����܂���
    //minP=PART[i].P;
	////////
	for(int k=0;k<PART[i].N;k++)
	{       
	    int j=PART[i].NEI[k];
	    if(PART[j].type!=OUTWALL)
	    {
	        double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
		
			double w=kernel(R,dis);
			W+=w;
		
			P_grad[A_X]+=(PART[j].P-minP)*X*w/(dis*dis);
			P_grad[A_Y]+=(PART[j].P-minP)*Y*w/(dis*dis);
			P_grad[A_Z]+=(PART[j].P-minP)*Z*w/(dis*dis);
		
		
			//P_grad[A_X]+=(PART[i].P-minP)*X*w/(dis*dis)/2;
			//P_grad[A_Y]+=(PART[i].P-minP)*Y*w/(dis*dis)/2;
			//P_grad[A_Z]+=(PART[i].P-minP)*Z*w/(dis*dis)/2;
		
	    }
	}////*/
	
	for(int D=0;D<DIMENTION;D++) if(W!=0) P_grad[D]= P_grad[D]*d/W;
	//for(int D=0;D<DIMENTION;D++) if(W!=0) P_grad[D]= P_grad[D]*d/n0;
}

//���͌��z�ɂ�鑬�x�E�ʒu�C���֐�
void modify_u_and_x_after_Pcalc(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **reU,double dt,double **Un)
{
	//���x�C���ɂ��ʒu�C��
	int sw=CON->get_temporary_r();		//ON�Ȃ�z��͌�ɉ��̈ʒu���v�Z���Ă���(�ʏ�̉��)
	int d=CON->get_dimention();

	//���x�X�V
	for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].u[D]+=reU[D][i];

	//�\�ʗ��q�̑��x���A��̂悤�ɋ��߂�̂ł͂Ȃ��A�������q����̋ߎ��ŋ��߂�ꍇ
	if(CON->get_interpolate_surface_u()==ON) reset_surface_velocity_by_WLSM(CON,PART, fluid_number);


	//�\�ʂ������瑬�x�̍X�V����߂�
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].u[D]=0;

	//�ʒu�X�V
	//if(CON->get_Position_by_WLSM()==ON) update_position_by_WLSM(CON,PART, fluid_number,dt,Un);//�ŏ����@�ɂ��ʒu�X�V
	//if(CON->get_Position_by_WLSM()==ON) update_position_by_WLSM_2(CON,PART, fluid_number,dt,Un);//�ŏ����@�ɂ��ʒu�X�V
	//else
	{
		if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*reU[D][i]*0.5;//��`��
		else if(sw==OFF) 
		{
			//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*PART[i].u[D];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*(PART[i].u[D]+Un[D][i])*0.5;//��`��
		}///*/
	}

	//�\�ʂ�������ʒu�̍X�V����߂�
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].r[D]-=dt*reU[D][i]*0.5;//��`��
}

//�A���ver,3
void negative3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double lamda,double N0,vector<point3D> &NODE,vector<element3D> &ELEM,int out)
{
	cout<<"FEM�ɂ��A��͎��s"<<endl;

	if(CON->get_dimention()==2) cout<<"2�����͔�Ή�"<<endl;

	point3D NODE0;
	element3D ELEM0;
	NODE.push_back(NODE0);//0�Ԗڂ̗v�f�����������Ă���
	ELEM.push_back(ELEM0);//0�Ԗڂ̗v�f�����������Ă���

	//INWALL�̐��𐔂���
	int inwall_number=0;
	for(int i=fluid_number;i<particle_number;i++) if(PART[i].type==INWALL && PART[i].surface==OFF) inwall_number++;//BDWALL�͌v�Z�ɂ͑g�ݍ��܂Ȃ�

	int node=fluid_number+inwall_number;		//�ߓ_��=���̗��q��
    int KTJ=node;							//�ő�ߓ_��
    int KTE=12*KTJ;								//�ő�v�f���@3���������Ƃ߂�
	int nelm=0;									//���݂̗v�f��
    double err=1.0e-14;							//�덷����̂������l
	double femfactor=1;//1000;							//�ł��؂�덷�΍�W��
	
	///��۰Ƃ��s�����ۂ��̔���
	int delaun_flag;							//��۰ƕ������s�����s��Ȃ���
	if(t==1 || t%100==0) delaun_flag=ON;
	else delaun_flag=OFF;
	///////////////////////
	
	if(delaun_flag==OFF)
	{
		nelm=NODE[0].material;
		node=NODE[0].particleID;
		//nelm=(int)ELEM.size();
		//node=(int)NODE.size()-(8+1);//�X�[�p�[�{�b�N�̂Ԃ��0�Ԗڂ̗v�f�͎󂯎���Ă͂����Ȃ�
	}


    if(delaun_flag==OFF)	//��۰Ƃ͂��Ȃ��Ă��ߓ_�̍��W�͕ύX����
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;//�ߓ_i�͗��qp
			for(int D=0;D<3;D++) NODE[i].r[D]=PART[p].r[D];
		}
	}
	else if(delaun_flag==ON)
	{
		/////////////input
		int num=1;
		for(int i=0;i<fluid_number;i++)
		{
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].material=FLUID;
			NODE[num].particleID=i;									//�Ή����闱�q�ԍ��i�[
			if(PART[i].surface==ON) NODE[num].boundary_condition=1;	//�Œ苫�E����
			else NODE[num].boundary_condition=0;
			num++;
		}
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL && PART[i].surface==OFF)//���̂Ɛڂ��Ă���ǂȂ�
			{
				NODE.push_back(NODE0);
				for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
				NODE[num].material=FLUID;		//�ʓ|�Ȃ̂ł����ł�WATER�ƒ�`
				NODE[num].particleID=i;			//�Ή����闱�q�ԍ��i�[
				NODE[num].boundary_condition=0;	//���R���E����
				num++;
			}
		}//////////////*/
    
		cout<<"�ߓ_��="<<node<<"  �ő�ߓ_��="<<KTJ<<"   �ő�v�f��="<<KTE<<endl;

		///femfacter�ɂ��g��
		for(int i=1;i<=node;i++) for(int D=0;D<3;D++) NODE[i].r[D]*=femfactor;

		for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);//�X�[�p�[�{�b�N�p�̐ߓ_���m��
		for(int i=1;i<=KTE;i++) ELEM.push_back(ELEM0);
	
	    /////////////�ߓ_���W�̐��K��
		double xmin=NODE[1].r[A_X];
		double ymin=NODE[1].r[A_Y];
		double zmin=NODE[1].r[A_Z];
		double xmax=xmin;
		double ymax=ymin;
		double zmax=zmin;

		///���W�̍ő�A�ŏ��l�����߂�
		for(int i=2;i<=node;i++)
		{
			if(NODE[i].r[A_X]<xmin) xmin=NODE[i].r[A_X];
			else if(NODE[i].r[A_X]>xmax) xmax=NODE[i].r[A_X];
	
			if(NODE[i].r[A_Y]<ymin) ymin=NODE[i].r[A_Y];
			else if(NODE[i].r[A_Y]>ymax) ymax=NODE[i].r[A_Y];
		
			if(NODE[i].r[A_Z]<zmin) zmin=NODE[i].r[A_Z];
			else if(NODE[i].r[A_Z]>zmax) zmax=NODE[i].r[A_Z];
		}
		////

		double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X�������̐��@
		double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y�������̐��@
		double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z�������̐��@
		double rmax=rax;							///�ő吡�@
		if(ray>rmax) rmax=ray;
		if(raz>rmax) rmax=raz;						//������else�ɂ�����_��
	
		///���W�ϊ�
		double rrm=1.000000/rmax;///�������������������邱�ƂŁA���l�덷�����点��E�E�H
		for(int i=1;i<=node;i++)
		{   //   A/B�Ƃ����v�Z�������Ƃ��A�`�̒l�ɂ���Ĕ�����1/B�Ƃ����{�����������Ă���̂ł͂Ȃ����ƍl���āA���̂悤�ȏ������ɂ��Ă���
		    NODE[i].r[A_X]=(NODE[i].r[A_X]-xmin)*rrm;
			NODE[i].r[A_Y]=(NODE[i].r[A_Y]-ymin)*rrm;
			NODE[i].r[A_Z]=(NODE[i].r[A_Z]-zmin)*rrm;
		}
		rax*=rrm;
		ray*=rrm;
		raz*=rrm;
		/////

		///��۰ƕ���
		int FINE_sw=OFF;
		delaun3D(CON,NODE,ELEM,KTJ,KTE,rax,ray,raz,&node,&nelm,FINE_sw,rrm);
	
		////���W�����ɖ߂�
		for(int i=1;i<=node;i++)
		{
			NODE[i].r[A_X]=rmax*NODE[i].r[A_X]+xmin;
			NODE[i].r[A_Y]=rmax*NODE[i].r[A_Y]+ymin;
			NODE[i].r[A_Z]=rmax*NODE[i].r[A_Z]+zmin;
		}
    
		for(int i=1;i<=node;i++) for(int D=0;D<3;D++) NODE[i].r[D]/=femfactor;
		/////ү����������

		NODE[0].material=nelm;//�v�f���Ȃǂ�ۑ�
		NODE[0].particleID=node;
	}
	/////////////////////////

	cout<<"�v�f��="<<nelm<<" �ߓ_����"<<node<<endl;

	double *P=new double[node+1];
	for(int i=0;i<=node;i++) P[i]=1;

	 /////�ߓ_-�v�f�֌W
    int *jnb=new int[node+1];///�e�ߓ_�ɗאڂ���v�f���i�[
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//�e�ߓ_�̎��ӗv�f�ԍ��i�[
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);
	
	//���͉�͂̂��߂�N0��PART[i].PND2��kernel2()�p�ɂ���������
	set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	cout<<"N0="<<N0<<endl;
	///////////////////*/

	///////////�L���v�f�@�v�Z�J�n
	P3D(CON,NODE,ELEM,node,nelm,P,jnb,PART,fluid_number,nei,dt,N0);

	reU3D(CON,NODE,ELEM,node,nelm,P,PART,fluid_number,dt,jnb);

    delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	delete [] P;

}

///���͌v�Z�֐��iFEM)
void P3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,int *jnb,vector<mpsparticle> &PART,int fluid_number,int **nei,double dt,double N0)
{
	cout<<"���͌v�Z�J�n--";
    double P0=0;					//�\�ʈ���
	double le=CON->get_distancebp();//���q�ԋ���
	double density=CON->get_density();
	double CO=density/(dt*dt)/N0;		//�v�Z�ɕK�v�ȌW��
	unsigned timeA=GetTickCount();	//�v�Z�J�n����

    int NN=0;					//�f�B���N���^���E�ߓ_��
    int *dn=new int [node+1];	//�e�ߓ_���f�B���N���^���E�̉��Ԗڂ��B�f�B���N���^���E��łȂ��Ȃ�node+1���i�[
    double *PHAT=new double [CON->get_max_DN()];//�f�B���N���^���l
	
	double *Dirichlet_P=new double [fluid_number];//�e���q�̃f�B���N���l

	if(CON->get_dir_for_P()==OFF) 
	{
		for(int i=0;i<fluid_number;i++) Dirichlet_P[i]=P0;
	}
	else //�\�ʗ��q�̈��͂Ƃ��āA�f�B���N���l
	{
		int flag=CON->get_dir_for_P();
		
		
		if(CON->get_dir_for_P()==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(CON->get_dir_for_P()==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(CON->get_dir_for_P()==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}

		//̧�ُo��
		ofstream gg("Dirichlet_P.dat");
		for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON ) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;
		gg.close();///*/

	}///////////

    ///�f�B���N���^���E��������
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			int PN=NODE[i].particleID;//MPS�ł̗��q�ԍ�
	        if(PN<fluid_number) PHAT[NN]=Dirichlet_P[PN];
			else PHAT[NN]=P0;
	        P[i]=PHAT[NN];
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"�ިظڐ�="<<NN<<" ";
	    
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


    /////////�S�̍s����쐬����
    
    int N[4+1]; //�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double U[4+1];//�e�_�̑��x
    double V[4+1];
    double W[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double dN[4+1];	//���q�����x�̍�
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==FLUID)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				int particle_num=NODE[N[j]].particleID;//MPS�ɂ����闱�q�ԍ�
				U[j]=PART[particle_num].u[A_X];
				V[j]=PART[particle_num].u[A_Y];
				W[j]=PART[particle_num].u[A_Z];
				dN[j]=PART[particle_num].PND2-N0;
			}
			///��۰ƕ����̍ۂɋ��߂��̐ς́A���߰�ޯ�����Ɉڍs���Čv�Z����Ă��邽�߂��͂␳�����l�ł͂Ȃ��B�Ȃ̂ł����ŋ��ߒ����B
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//�̐ς�6�{�ł��邱�Ƃɒ���
			
			double delta6=ELEM[je].volume;//�̐ς�6�{
		
			delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
			double delta=ELEM[je].volume/6;//�{���̑̐�
			double co1=delta/20.0;			//�W��

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
		
			double udiv=0;//���x���U
			udiv+=c[1]*U[1]+c[2]*U[2]+c[3]*U[3]+c[4]*U[4];
			udiv+=d[1]*V[1]+d[2]*V[2]+d[3]*V[3]+d[4]*V[4];
			udiv+=e[1]*W[1]+e[2]*W[2]+e[3]*W[3]+e[4]*W[4];
			udiv*=density/dt;
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
									G[I][h]+=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
				    
								G[I][H]+=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								ROW[I][H]=J;
							}
						}
						else //N[j]���ިظڌ^���E�ߓ_�Ȃ�
						{
							int n=dn[N[j]];
							B[I-1]-=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
						}
					}
					B[I-1]+=-udiv/4*delta;
					/*double A=0;				//���q�����x�̍���v�f�ő̐ϐϕ������l
					for(int j=1;j<=4;j++)
					{
						if(j==i) A+=co1*2*dN[j];
						else A+=co1*dN[j];
					}
					B[I-1]+=CO*A;*/
				}
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
	ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		P[i]=XX[n];
	}
	delete [] XX;
	////////////////////////////
    
    ///////���͂�̧�ُo��
	ofstream fp("P.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<0.5*le && NODE[i].r[A_Y]>-0.5*le) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Z]<<"\t"<<P[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

	delete [] Dirichlet_P;
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

}

///���x�C���֐�(FEM)
void reU3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,vector<mpsparticle> &PART,int fluid_number,double dt,int *jnb)
{
	unsigned timeA=GetTickCount();		//�v�Z�J�n����
    
    int N[4+1];							//�v�f�̊e�ߓ_�ԍ��i�[
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double density=CON->get_density();
	double A=dt/density;				//�悭�g���̂ŌW����

	double *newU[3];
	for(int D=0;D<3;D++) newU[D]=new double [node+1];	//�e�ߓ_�̐V�������x�i�[
	double *reU[3];
	for(int D=0;D<3;D++) reU[D]=new double [nelm+1];	//�e�v�f�̑��x�C���ʊi�[

	int *num=new int[node+1];
	for(int i=1;i<=node;i++) num[i]=0;
	///newU�ɒl���i�[
	for(int i=1;i<=node;i++)
	{
		int particle_id=NODE[i].particleID;	//�ߓ_i��particle_id�̗��q
		for(int D=0;D<3;D++)  newU[D][i]=0;//������
	}////////////

    ///�e�ߓ_�̑��x���X�V����
    for(int je=1;je<=nelm;je++)
    {
		for(int D=0;D<3;D++) reU[D][je]=0;
		
        for(int j=1;j<=4;j++)
		{
			N[j]=ELEM[je].node[j];
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
	
		double delta6=ELEM[je].volume;//�̐ς�6�{  �̐ς͈��͋��߂�Ƃ��Ɍv�Z���Ă���
	
		delta6=1/delta6;//�v�Z�ɕK�v�Ȃ̂͋t���Ȃ̂ł����Ŕ��]���Ă���
		
		///�W��c,d,e�v�Z
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
    
		reU[A_X][je]=c[1]*P[N[1]]+c[2]*P[N[2]]+c[3]*P[N[3]]+c[4]*P[N[4]];
        reU[A_Y][je]=d[1]*P[N[1]]+d[2]*P[N[2]]+d[3]*P[N[3]]+d[4]*P[N[4]];
		reU[A_Z][je]=e[1]*P[N[1]]+e[2]*P[N[2]]+e[3]*P[N[3]]+e[4]*P[N[4]];

		for(int j=1;j<=4;j++)
		{
			newU[A_X][N[j]]-=A*reU[A_X][je];
			newU[A_Y][N[j]]-=A*reU[A_Y][je];
			newU[A_Z][N[j]]-=A*reU[A_Z][je];
			
			num[N[j]]=num[N[j]]+1;//jnb�ƈ�v����̂ł́H
		}////*/

    }///���x�X�V�I��
    
	//���x�X�V�𗱎q�ɕϊ�
	for(int i=1;i<=node;i++)
	{
		int id=NODE[i].particleID;
		if(id<fluid_number)
		{	
			for(int D=0;D<3;D++)
			{
				if(num[i]<=1)cout<<"num<=1 ??"<<endl;
				PART[id].u[D]+=newU[D][i]/num[i];
				if(CON->get_temporary_r()==OFF) PART[id].r[D]+=PART[id].u[D]*dt;
				if(CON->get_temporary_r()==ON)  PART[id].r[D]+=newU[D][i]/num[i]*dt;
			}
		}
	}

	for(int D=0;D<3;D++) delete [] newU[D];
	for(int D=0;D<3;D++) delete [] reU[D];
	delete [] num;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//�A���ver,2
void negativeP_iterative(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double lamda,double n0,double N0,int out,double **Un)
{
	cout<<"�����@�ɂ��A��͎��s"<<endl;

	int d=CON->get_dimention();						//����
	double le=CON->get_distancebp();				//�������q�ԋ���
	double r2=CON->get_re2()*CON->get_distancebp();
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();


	int count2=0;									//������
	
	double error=0;
	unsigned timeA=GetTickCount();					//�v�Z�J�n����


	int calctype=0;									//�v�Z�^�C�v 0:���x�̂ݍX�V 1:�ʒu���X�V
	int pn=0;										//���͂��������߂̘A�����������m��
	
	double div_error=0;								//���U�덷
	int numf=0;										//�������̗��q��
	int divflag=2;									//���U�v�Z�׸� 1:MPS 2:WLSM 3:Dndt
	int Btype=CON->get_B_of_P();					//���s��̎��
	
	if(Btype==3) divflag=3;

	for(int i=0;i<particle_number;i++) if(PART[i].surface==OFF) if(PART[i].type==FLUID || PART[i].type==INWALL) pn++;//�������̐�

	double *reU[DIMENTION];//���x�C����
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	double *old_U[DIMENTION];//���̑��x�L��
	for(int D=0;D<DIMENTION;D++) old_U[D] = new double [fluid_number];

	double *old_r[DIMENTION];//���̑��x�L��
	for(int D=0;D<DIMENTION;D++) old_r[D] = new double [fluid_number];

	double *old_reU[DIMENTION];//1�����O�̑��x�C���ʋL��
	for(int D=0;D<DIMENTION;D++) old_reU[D] = new double [fluid_number];

	double *direct[DIMENTION];			//�������@���޸�يi�[
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];

	double *P_grad[DIMENTION];			//���͌��z�޸�يi�[
	for(int D=0;D<DIMENTION;D++) P_grad[D]=new double [fluid_number];

	double *udiv=new double [particle_number];	//�e���q�̑��x���U(�Ǘ��q���܂�)

	//���͉�͂̂��߂�N0��PART[i].PND2��kernel2()�p�ɂ���������
	set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	cout<<"N0="<<N0<<endl;
	///////////////////*/

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<d;D++)
		{
			old_U[D][i]=PART[i].u[D];//���̑��x���L��
			old_reU[D][i]=0;//������
			old_r[D][i]=PART[i].r[D];
		}
	}

	//////�@���޸�ٌv�Z�J�n
	if(CON->get_Pgrad()==2 ||CON->get_Pgrad()==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].type==BOFLUID ) direct_f(CON,PART,i,direct);
			else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
		}
		
	}/////////////////*/

	///�e���q�̑��x���U�v�Z
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL)
		{
			if(divflag==1) udiv[i]=divergence(CON,PART,i,n0);
			else if(divflag==2)
			{
				udiv[i]=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*udiv[i]!=udiv[i]+udiv[i]) udiv[i]=divergence(CON,PART,i,n0);//udiv[i]���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
			}
			else if(divflag==3) udiv[i]=Dndt(CON,PART,i,n0);
		}
	}///////

	//�ŏ��̑��x���U���z�o��
	ofstream f0("div0.dat");
	if(d==2) for(int i=0;i<fluid_number;i++) f0<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<udiv[i]<<endl;
	if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) f0<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<udiv[i]<<endl;
	f0.close();////////////*/

	///���x�̔��U�덷���v�Z
	div_error=0;
	numf=0;//�������̗��q��
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			div_error+=udiv[i]*udiv[i];
			numf++;
		}
	}
	if(numf!=0) 
	{
		div_error=div_error/numf;
		div_error=sqrt(div_error);	//�W���΍�
	}


	int *ppn = new int[pn];					//�s��ɂ������n�Ԗڂ̖��m���͗��q�ԍ�ppn[n]�̗��q�ɑ���
	int *link = new int [particle_number];	//���q�ԍ�i��link[i]�Ԗڂ̖��m��
	double *B   = new double[pn];			//���s��

	int count=0;
	///ppn�z���link�z��쐬
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			if(PART[i].type==FLUID || PART[i].type==INWALL)
			{
				ppn[count]=i;
				link[i]=count;
				count++;
			}
			else link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
		}
		else link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
	}
	//////*/

	int number=0;			//�W���s��̔�[���v�f��
	for(int n=0;n<pn;n++)
	{   
		number++;///�������g�𐔂ɂ����
		int i=ppn[n];//n�Ԗڂ̖��m���͗��qi
	  
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
		}
	}
	////number�����܂���

	double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	double real_lamda=lamda;	//lamda�̒l�L��
	if(CON->get_HL_sw()==ON) lamda=2.0*CON->get_dimention()*real_lamda;//�����̗��U���̏ꍇ�Alamda=2d�ƍĒ�`����΃�/2d�̍��͏�����B��������������Ɨ��ӂ̒l���傫���Ȃ肷����̂ŁAreak_lamda�𗼕ӂɂ����ď��������Ă���
		
	

	
	///////////////////////////////////�����J�n
	int flag=0;
	while(flag==0)
	{
		unsigned int timeC=GetTickCount();	//�v�Z�J�n����
	
		cout<<"���͖��m��:"<<pn<<" �ɑ΂���";

		int remake_SW=OFF;			//�W���s����ēx�쐬���邩�A���Ȃ���
		if(calctype==0) if(count2==0) remake_SW=ON;	//���x�̂ݏC������ꍇ�͍ŏ������A
		if(calctype==1) remake_SW=ON;//�ʒu���C������ꍇ�͖�����Ȃ���

		if(CON->get_HL_sw()==OFF && remake_SW==ON)//�W���I�ȗ��U��
		{
			int index=0;
			for(int n=0;n<pn;n++)
			{
				ptr[n]=index;
				int i=ppn[n];
				int kk=index;//�l��ۑ�
				ind[index]=n;
				index++;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						    
						double w=kernel2(r2,dis,d);
						val[index]=w;
						ind[index]=m;
						index++;
						W+=w;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel2(r2,dis,d);
						W+=w; //������w���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^ 
					}
				}
				val[kk]=-W;
			}
			ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
		}
		////////////////////*/ 

		//if(count2==CON->get_iteration_count()) Btype=4;		//��ԍŌ�̔��������͖��x�^

		///////////////////////////////////////////���s��B�̌v�Z
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//�s���n�Ԗڂ̗��q��i
			if(Btype==0) B[n]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//���ȏ�
			else if(Btype==1)//���x���U
			{    
				double div=udiv[i];
				B[n]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
			}
			else if(Btype==2)//(���ȏ�+���x���U)/2
			{
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				double div=udiv[i];
				
				B[n]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[n]/=a+b;
			}
			else if(Btype==3)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U
			{    
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(Btype==4)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U+PND
			{    
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//���ȏ�
				B[n]/=a+b;
			}
			else if(Btype==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			else cout<<"���s�񖢉��� check B_of_P"<<endl;
		}
		////���s��B�����Ƃ܂���//*/
	
	
		if(CON->get_dir_for_P()!=OFF)//�\�ʗ��q�̈��͂Ƃ��āA�f�B���N���l
		{
			int flag=CON->get_dir_for_P();
			if(flag==2 || flag==3)
			{
				if(CON->get_EM_method()==OFF) flag=1;//EM�ɂ��v�Z���s��Ȃ��̂�flag��2��3�Ȃ�ԈႢ�Ȃ̂�1�ɖ߂�
			}
			double *Dirichlet_P=new double [particle_number];
			for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

			if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
			else if(flag==3)
			{
				set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
				set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
			}

			for(int i=fluid_number;i<particle_number;i++)//BDWALL��Dirichlet_P�v�Z
			{
				if(PART[i].surface==ON)//�Ǖ\�ʗ��q�Ȃ�
				{
					double P=0;
					double W=0;
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						if(PART[j].type==FLUID && PART[j].surface==ON)
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							double w=kernel(r2,dis);
							P+=Dirichlet_P[j]*w;
							W+=w;
						}
					}
					if(W!=0) P/=W;
					Dirichlet_P[i]=P;
				}
			}

			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].surface==ON)
				{
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						int m=link[j];
						if(m<pn)//���qj�����m���Ȃ�
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
				    
							double w=kernel2(r2,dis,d);
							if(CON->get_HL_sw()==ON) w*=real_lamda;
							B[m]-=w*Dirichlet_P[i];
							PART[i].P=Dirichlet_P[i];
						}
					}
				}
			
			}

			//̧�ُo��
			if(count2==0)
			{
				ofstream gg("Dirichlet_P.dat");
				if(d==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
				else if(d==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
				gg.close();///*/
	
				
			}
			delete [] Dirichlet_P;
		}///////////
	
	
		/////////////////////////////////////CG�@
		double *r=new double[pn];
		double *X=new double[pn];
		double *AP = new double [pn];
		double *P = new double [pn];

		/////////////////////////�����l//////////////////
		if(CON->get_initialP()==OFF)
		{
			for(int n=0;n<pn;n++) 
			{
				 X[n]=0;
				 r[n]=B[n];
				 P[n]=r[n];
			}
		}
		else if(CON->get_initialP()==ON)//�����l�Ƃ��Č��݂̏���^����B�����ɑ����Ȃ�
		{
			for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
			for(int n=0;n<pn;n++)
			{
				double AX=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
				r[n]=B[n]-AX;
				P[n]=r[n];
			}
		}
		//////////////////////////////////////////////

		//CG�@�ɂ��s�������
		if(CON->get_solution()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_CGep());
		else if(CON->get_solution()==1)
		{
			for(int n=0;n<pn;n++)//ICCG�@�łƂ��ꍇ�A�W���s���������ƕ��ѕς����Ȃ��Ă͂Ȃ�Ȃ�
			{      
				int num=ptr[n+1]-ptr[n];
				for(int j=ptr[n]+1;j<ptr[n+1];j++)
				{
					for(int m=ptr[n];m<j;m++)
					{
						if(ind[j]<ind[m])
						{
							double temp=val[m];
							int tempR=ind[m];
							val[m]=val[j];
							ind[m]=ind[j];
							val[j]=temp;
							ind[j]=tempR;
						}
					}
				}
			}
			iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_CGep(),&count);
		}
		
		if(CON->get_negativeP()==OFF)//�������l�����Ȃ��ꍇ
		{
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				PART[i].P=X[n];
				if(PART[i].P<0) PART[i].P=0;
			}
		}
		else if(CON->get_negativeP()==ON)//�������l������ꍇ
		{
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				PART[i].P=X[n];
			}
		}
		
		delete [] r;
		delete [] X;
		delete [] AP;
		delete [] P;

		//////////////////////////////*/

		cout<<"������:"<<count<<"  time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;

		if(count2==0) plot_P(CON ,PART,particle_number,t,fluid_number);//�ŏ��݈̂��͏o��

		///////////////////////////////////////////////////���͌��z�v�Z
		cout<<"���͌��z�v�Z�J�n ver."<<CON->get_Pgrad()<<" ---------";
	
		unsigned int timeB=GetTickCount();	//�v�Z�J�n����

		int minPsw=CON->get_minP();

		if(CON->get_Pgrad()==3)//�\�ʂ̂ݖ@�� && minP=0�̂Ƃ��\�ʗ��q����������͌v�Z
		{
			P_gradient3(CON,PART,fluid_number,direct,dt,reU,P_grad,minPsw);
		}
		else if(CON->get_Pgrad()==4)//�d�݂��ŏ����@(WLSM)
		{
			P_gradient4(CON,PART,dt,fluid_number,reU,P_grad,minPsw);
		}
		else if(CON->get_Pgrad()==5)//���g��ʂ�Ȃ��d�݂��ŏ����@(WLSM)
		{
			P_gradient5(CON,PART,dt,fluid_number,reU,P_grad);
		}
	
		cout<<"ok  time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;
		//////////////////////////////////////////////////////////////////////////////*/

		//���x�C���ʂ��v�Z
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*P_grad[D][i]/density[i];
		
		///���x�C���ʂ̏C��
		int mode=0;
		if(mode==1) modify_reU(CON,PART,fluid_number,particle_number,t,dt,n0,udiv,reU);


		//���x�C��
		if(calctype==0)//���x�̂ݏC��
		{
			for(int i=0;i<fluid_number;i++)
			{
				//double divreU=div_of_reU(CON,PART,reU,i,fluid_number);
				double w=1;//0.8;
				for(int D=0;D<d;D++)
				{
					PART[i].u[D]+=w*reU[D][i];//�Ȃ�0.6�`0.9���炢�̔{���������瑁����������
					//PART[i].u[D]+=w*reU[D][i]+(1-w)*old_reU[D][i];
					old_reU[D][i]=reU[D][i];//�X�V
					
				}
			}
		}
		else if(calctype==1)//�ʒu���C���B���̏ꍇ�͏d�݊֐����Čv�Z�B
		{
			int sw=CON->get_temporary_r();		//ON�Ȃ�z��͌�ɉ��̈ʒu���v�Z���Ă���(�ʏ�̉��)
			for(int i=0;i<fluid_number;i++)
			{
				for(int D=0;D<d;D++)
				{
					PART[i].u[D]+=reU[D][i];
					//PART[i].r[D]+=reU[D][i]*dt;
					if(sw==ON) PART[i].r[D]+=reU[D][i]*dt*0.5;
					else 
					{
						PART[i].r[D]=old_r[D][i]+dt*(PART[i].u[D]+old_U[D][i])*0.5;
						//PART[i].r[D]+=dt*(PART[i].u[D]+old_U[D][i])*0.5;//��`�� old_U�͂��̊֐��ɓ˓������ۂ̑��x�B�����̊Ԃ����ƋL�����Ă�
					}
				}
			}

			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type!=OUTWALL)
				{
					double W=0;
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel2(r2,dis,d);
						W+=w;
					}
					PART[i].PND2=W;
				}
			}
		}////////////////////////////

		///���x�̔��U�덷���v�Z
		div_error=0;
		numf=0;//�������̗��q��
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=OUTWALL && PART[i].surface==OFF)
			{
				if(divflag==1) udiv[i]=divergence(CON,PART,i,n0);
				else if(divflag==2)
				{
					udiv[i]=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
					if(2*udiv[i]!=udiv[i]+udiv[i]) udiv[i]=divergence(CON,PART,i,n0);//udiv[i]���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				}
				else if(divflag==3) udiv[i]=Dndt(CON,PART,i,n0);
			}
			if(PART[i].type==FLUID && PART[i].surface==OFF)
			{
				div_error+=udiv[i]*udiv[i];
				numf++;
			}
		}
		if(numf!=0) 
		{
			div_error=div_error/numf;
			div_error=sqrt(div_error);	//�W���΍�
		}///////*/


		if(count2>0 && div_error>=error) flag=1;//2��ڈȍ~�̌v�Z�Ō덷�������Ă��܂�����v�Z���~
		
		cout<<div_error<<endl;

		error=div_error;//�덷�X�V�i�ǂ�ǂ񏬂����Ȃ�j
		count2++;

		if(count2>CON->get_iteration_count()) flag=1;//�E�o
		if(mode==1) flag=1;
		
	}
	lamda=real_lamda;//lamda�̒l��߂�

	//���x�C���ƈʒu�C��
	if(calctype==0) modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);

	//�ŏI�I�ȑ��x���U���z�o��
	ofstream ff("div.dat");
	if(d==2) for(int i=0;i<fluid_number;i++) ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<udiv[i]<<endl;
	if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<udiv[i]<<endl;
	ff.close();

	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_U[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_r[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_reU[D];
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] P_grad[D];

	delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] val;
	delete [] ind;
	delete [] ptr;

	delete [] udiv;

	delete [] density;

	cout<<"ok count="<<count2-1<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//���x�C���ʂ̏C���֐�
void modify_reU(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double n0,double *udiv, double **reU)
{
	cout<<"reU�̏C���J�n--";
	///�����ł��炤n0�͊g�U��n0�Ȃ̂ł́H�v�`�F�b�N
	double le=CON->get_distancebp();
	double R=CON->get_re()*le;
	int    d=CON->get_dimention();
	double dir_val=1;				//�f�B���N���l�@���݂�1���g�p ���Ƃ��Ε\�ʗ��q�����͌l����Ɍ��߂�����ިظڒl�ɂ�������ǂ��Ȃ�H�������H

	int pn=0;	//���m��
	for(int i=0;i<fluid_number;i++) if(PART[i].surface==OFF) pn++;//���݂̂Ƃ���A�\�ʗ��q�͌Œ苫�E�Ƃ��Ă���


	int *ppn = new int[pn];					//�s��ɂ������n�Ԗڂ̖��m���͗��q�ԍ�ppn[n]�̗��q�ɑ���
	int *link = new int [particle_number];		//���q�ԍ�i��link[i]�Ԗڂ̖��m��
	double *B   = new double[pn];			//���s��

	int count=0;
	///ppn�z���link�z��,����щ��s��a�쐬
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			ppn[count]=i;
			link[i]=count;
			B[count]=-udiv[i];
			count++;
		}
		else link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
	}	
	//////*/

	int number=0;			//�W���s��̔�[���v�f��
	for(int n=0;n<pn;n++)
	{   
		number++;			//�������g�𐔂ɂ����
		int i=ppn[n];		//n�Ԗڂ̖��m���͗��qi
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			int m=link[j];
			if(m<pn) number++;
		}
	}
	////number�����܂���

	double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	
	int index=0;
	for(int n=0;n<pn;n++)
	{
	    ptr[n]=index;
	    int i=ppn[n];
		int kk=index;//�l��ۑ�
	    val[index]=-PART[i].PND;
	    ind[index]=n;
	    index++;
		double W=0;
		double d_val=0;//�Ίp����̒l
		//n0=PART[i].PND;
	    for(int k=0;k<PART[i].N;k++)
	    {
	        int j=PART[i].NEI[k];
			
			int m=link[j];
			if(m<pn)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
				val[index]=d/n0*w*(reU[A_X][j]*X+reU[A_Y][j]*Y+reU[A_Z][j]*Z)/(dis*dis);
				ind[index]=m;
				index++;
				W+=w;
	
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}
			else if(PART[j].type==FLUID && PART[j].surface==ON)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
					
				B[n]-=d/n0*w*(reU[A_X][j]*X+reU[A_Y][j]*Y+reU[A_Z][j]*Z)/(dis*dis)*dir_val;
				W+=w; 
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}
			else	//�Ǘ��q
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
					
				W+=w; 
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}//*/
			
		}
		val[kk]=d_val;
		//val[kk]=-W;
	}
	ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
	cout<<"�s��쐬 ";
	////////////////////*/
	

	/////////////////////////////////////GCR�@
	double alp;
	double E=1;//�덷
	unsigned int timeC=GetTickCount();
	
	double *r=new double[pn];
	double *XX=new double[pn];
	for(int n=0;n<pn;n++) XX[n]=0;
	
	double *AP = new double [pn];
	double *Ar = new double [pn];
	double *P = new double [pn];
	double *nextP = new double [pn];

	int max=500;
	double **oldP=new double *[max];
	for(int m=0;m<max;m++) oldP[m]=new double [pn];
	double **oldAP=new double *[max];
	for(int m=0;m<max;m++) oldAP[m]=new double [pn];
	double *APjAPj=new double [max];

	/////////////////////////�����l//////////////////
	
	for(int n=0;n<pn;n++) XX[n]=dir_val;
	/*for(int n=0;n<pn;n++)
	{
		double AX=0;
		for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*XX[ind[j]];
		r[n]=B[n]-AX;
		P[n]=r[n];
	}/////////*/
		
	//////////////////////////////////////////////

	double EP=1e-2;//CON->get_CGep();//��������(convergence test)

	
	cout<<"GCR�@�X�^�[�g------";//������۸��т́w�V���l�v�Z�x���� �`����P46���Q�l�ɂ��Ă���
	for(int n=0;n<pn;n++)
	{
		double AX=0;
		for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*XX[ind[j]];
		r[n]=B[n]-AX;
		P[n]=r[n];
	}/////////*/
	
	count=0;
	
	while(E>EP && count<max)
	{
		for(int n=0;n<pn;n++) oldP[count][n]=P[n];//�T�������޸�ق̕ۑ�

		//////////////alp�����߂� alp=(ri,APi)/(APi,APi)
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
		}
		for(int n=0;n<pn;n++) oldAP[count][n]=AP[n];//�l��ۑ�


		double rAP=0;
		for(int n=0;n<pn;n++) rAP+=r[n]*AP[n];
		double APAP=0;
		for(int n=0;n<pn;n++) APAP+=AP[n]*AP[n];
		alp=rAP/APAP;
			//cout<<"alp="<<alp<<" "<<rAP<<" "<<APAP<<endl;
		//////////////////////

		//(AP,AP)��ۑ�
		APjAPj[count]=APAP;
		////////////
		
		//////////////// ���X�V�@X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) XX[n]+=alp*P[n];
		//////////////////////////////
			
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
			
		//////////////////�덷
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
			
		////////////////////////�Ƃ肠�������̒T�������Ƃ��Ďc�������i�[
		for(int n=0;n<pn;n++) nextP[n]=r[n];
		///////////////////////////////////

		////////////////////////����P �ߋ��̒l���g���Čv�Z
		double betaP=0;//����P
		for(int n=0;n<pn;n++)
		{      
			Ar[n]=0;		
			for(int j=ptr[n];j<ptr[n+1];j++) Ar[n]+=val[j]*r[ind[j]];
		}

		for(int j=0;j<=count;j++)
		{
			double ArAPj=0;
			for(int n=0;n<pn;n++) ArAPj+=Ar[n]*oldAP[j][n];
			double betaj=-ArAPj/APjAPj[j];
			
			for(int n=0;n<pn;n++) nextP[n]+=betaj*oldP[j][n];
		}//////////*/
			
		///////////////////// P�̍X�V
		for(int n=0;n<pn;n++) P[n]=nextP[n];
		
		count++;
	}
	
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		for(int D=0;D<d;D++) reU[D][i]*=XX[n];

	}

	ofstream fx("kreU.dat");
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fx<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<XX[n]<<endl;

	}
	fx.close();
		
	///�z��폜
	delete [] r;
	delete [] XX;
	delete [] AP;
	delete [] Ar;
	delete [] P;
	delete [] nextP;

	//////////////////////////////*/

	cout<<"������:"<<count<<"  time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;


	delete [] ppn;
	delete [] link;
	delete [] B;

	delete [] val;
	delete [] ind;
	delete [] ptr;

	
	for(int m=0;m<max;m++) delete [] oldP[m];
	delete [] oldP;
	for(int m=0;m<max;m++) delete [] oldAP[m];
	delete [] oldAP;
	delete [] APjAPj;
	//////////

	
}

void calc_P_main_with_gridless(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P,int t)
{
	///���͌v�Z�ɂ����āA�g�p����d�݊֐��͌��z�ȂǂɎg�p���邻��Ƃ͈قȂ���̂�p����B
	//���R�́A���z�Ɏg�p����d�݊֐��͒P�Ȃ�d�ݕ��ςȂ̂ŉ������悤���Ă��悢���A
	//PPE���̗��U���Ɏg�p����d�݊֐��́A���q�����x�喧�x�ƂȂ�悤�ɐݒ肷��K�v�����邩��B
	//�������҂œ���̏d�݊֐����g�p��������΁Akernel2()�̒��g��kernel�Ɠ���ɏ��������邱�ƁB

	//pn:���͂��������߂̘A�����������m��
	cout<<"gridless�ɂ�鈳�͌v�Z"<<endl;

	double le=CON->get_distancebp();
	double r2=CON->get_re2()*le;
	unsigned int timeA=GetTickCount();					//�v�Z�J�n����
	double dimention=2;
	if(CON->get_dimention()==3) dimention=3;
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();
	//double lamda=calclamda(CON); //���׼�ݗp��

	//���͉�͂̂��߂�N0��PART[i].PND2��kernel2()�p�ɂ���������
	//set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	//cout<<"N0="<<N0<<endl;
	///////////////////*/

	//���̊֐��ɓ������i�K��pn�����܂��Ă���d�l�ɂȂ��Ă��邪�Agridless�@�ōs���ꍇ�́A���ӗ��q�������Ȃ��������q�̉��Z���Ȃ������B
	//���������pn���ω����邱�ƂɂȂ�B�����ŁA�ēxpn�����߂�B

	int min_neiber=4;				//���ӗ��q���������菭�Ȃ��ꍇ�͈��͂��v�Z���Ȃ�
	if(dimention==3) min_neiber=6;
	pn=0;
	for(int i=0;i<particle_number;i++) if(PART[i].type!=OUTWALL && PART[i].surface==OFF && PART[i].N2>=min_neiber) pn++;

	cout<<"���͖��m��:"<<pn<<" ";

	int *ppn = new int[pn];					//�s��ɂ������n�Ԗڂ̖��m���͗��q�ԍ�ppn[n]�̗��q�ɑ���
	int *link = new int [particle_number];	//���q�ԍ�i��link[i]�Ԗڂ̖��m��
	double *B   = new double[pn];			//���s��
	
	int count=0;

	///ppn�z���link�z��쐬
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL && PART[i].surface==OFF && PART[i].N2>=min_neiber)
		{
			ppn[count]=i;
			if(B_of_P==0)
			{
				//B[count]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//���ȏ�
				B[count]=-density[i]*(PART[i].PND2-N0)/(dt*dt*N0);//���ȏ�
			}
			else if(B_of_P==1)//���x���U
			{    
				//double div=divergence(CON,PART,i,n0);
				double div;
				if(CON->get_divU_method()==1) div=divergence(CON,PART,i,n0);
				else if(CON->get_divU_method()==2) div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());//�ŏ����@�i���g��ʂ�Ȗʁj
				else if(CON->get_divU_method()==3) div=divergence3(CON,PART,i,CON->get_interpolate_surface_u());//�ŏ����@�i���g��ʂ�Ƃ͌���Ȃ��Ȗʁj
				else if(CON->get_divU_method()==4) div=divergence4(CON,PART,i);//������
				if(2*div!=div+div)
				{
					div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
					//cout<<"div_err i="<<i<<" "<<div<<endl;
				}
				
				//B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
				//B[count]=div*lamda*PART[i].PND2/(dt*2*CON->get_dimention())*density[i];
				B[count]=div/(dt)*density[i];
			}
			else if(B_of_P==2)//(���ȏ�+���x���U)/2
			{
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				
				B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[count]/=a+b;
			}	
			else if(B_of_P==3)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U
			{    
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==4)//�d�݊֐��̒��ڔ����ɂ�鑬�x���U+PND
			{    
				double a=CON->get_w_div();double b=1;///�e��@�̏d��
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//���ȏ�
				B[count]/=a+b;
			}
			else if(B_of_P==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				//double div=Dndt(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//div���G���[�Ŕ�����̂Ƃ���MPS�Ōv�Z����
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			link[i]=count;
			count++;
		}
		else
		{
			//if(PART[i].N2<=6) cout<<i<<" "<<PART[i].N2<<endl;
			link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
		}
	}
	//////*/

	if(CON->get_dir_for_P()!=OFF && B_of_P!=0)//�\�ʗ��q�̈��͂Ƃ��āA�f�B���N���l
	{
		int flag=CON->get_dir_for_P();
		if(flag==2 || flag==3)
		{
			if(CON->get_EM_method()==OFF) flag=1;//BEM�ɂ��v�Z���s��Ȃ��̂�flag��2��3�Ȃ�ԈႢ�Ȃ̂�1�ɖ߂�
		}
		double *Dirichlet_P=new double [particle_number];
		for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

		if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(flag==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}
		
		for(int i=fluid_number;i<particle_number;i++)//�Ǖ\�ʗ��q��Dirichlet_P�v�Z
		{
			if(PART[i].surface==ON)//�Ǖ\�ʗ��q�Ȃ�
			{
				double P=0;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					if(PART[j].type==FLUID && PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(r2,dis);
						P+=Dirichlet_P[j]*w;
						W+=w;
					}
				}
				if(W!=0) P/=W;
				Dirichlet_P[i]=P;
			}
		}

		//�\�ʗ��q��PART[i].P�Ƀf�B���N�����E��������
		ofstream gg2("Dirichlet_error.dat");
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON) PART[i].P=Dirichlet_P[i];
			else if(PART[i].type==FLUID && Dirichlet_P[i]!=0)
			{
				cout<<"�������q�Ȃ̂��ިظڒl���i�[����Ă��܂��Bi="<<i<<" �l��"<<Dirichlet_P[i]<<endl;
				gg2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		gg2.close();

		//̧�ُo��
		ofstream gg("Dirichlet_P.dat");
		if(dimention==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
		else if(dimention==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
		gg.close();////

		//�f�B���N���l���x�N�g���\��
		output_dirichlet_vector_files(CON,PART,fluid_number,flag,Dirichlet_P);

		delete [] Dirichlet_P;
	}//////////*/

	///���s��o��
	ofstream h("Bmatrix_for_P.dat");
	if(CON->get_dimention()==2) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<B[n]<<endl;
		}
	}
	if(CON->get_dimention()==3) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			//cout<<B[n]<<endl;
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B[n]<<endl;
		}
	}
	h.close();
	/////////////*/

	
	int number=0;			//�W���s��̔�[���v�f��
	for(int n=0;n<pn;n++)
	{   
	    number++;///�������g�𐔂ɂ����
	    int i=ppn[n];//n�Ԗڂ̖��m���͗��qi
	  
	    for(int k=0;k<PART[i].N2;k++)
	    {
	        int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
	    }
	}
	////number�����܂���
	
    double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	unsigned int timeC=GetTickCount();
	ofstream ft("gridless_to_MPS.dat");				//���U����@��gridless����MPS�ɕύX�ɂȂ������q�̍��W���o��
	if(CON->get_dimention()==2)//3����
	{
		int N=5;
		double *matrix=new double [N*N];	//N�~N�̌W���s��
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//�l��ۑ�
			ind[index]=n;
			index++;
			double valII=0;						//�Ίp�����̒l
			for(int k=0;k<N*N;k++) matrix[k]=0;//������
			double B0=B[n];			//B[n]�̒l��ۑ�
			int index0=index;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				
				int m=link[j];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
				//if(dis>1) w=1/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
				//if(PART[j].surface==ON) w=10;	//2��A�����̂Ɠ�����
				//if(PART[j].surface==ON) w=1;//�\�ʂ͂�����ƒʂ��ė~��������A�d�݂�傫���Ƃ�

				matrix[0]+=X*X*w;			//��Xjwj
				matrix[1]+=X*Y*w;		//��XjYjwj
				matrix[2]+=X*X*X*w;			//��Xj^3wj
				matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
				matrix[4]+=X*Y*Y*w;			//��XjYj^2wj
	
				matrix[6]+=Y*Y*w;			//��Yj^2wj
				matrix[9]+=Y*Y*Y*w;			//��Yj^3wj
	
				matrix[12]+=X*X*X*X*w;			//��Xj^4wj
				matrix[13]+=X*X*X*Y*w;			//��Xj^3Yjwj
				matrix[14]+=X*X*Y*Y*w;			//��Xj^2Yj^2wj
		
				matrix[19]+=X*Y*Y*Y*w;			//��XjYj^3wj
	
				matrix[24]+=Y*Y*Y*Y*w;			//��Yj^4wj
	
				//B[0]+=dP*X*w;//��dPjXjwj
				//B[1]+=dP*Y*w;//��dPjYjwj
				//B[2]+=dP*X*X*w;//��dPjXj^2wj
				//B[3]+=dP*X*Y*w;//��dPjXjYjwj
				//B[4]+=dP*Y*Y*w;//��dPjYj^2wj
				
				/*if(PART[j].surface==ON)
				{
					X*=2; Y*=2; Z*=2; dis*=2;
					
					//if(dis>le) w=le*le/(dis*dis);
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�

					matrix[0]+=X*X*w;			//��Xjwj
					matrix[1]+=X*Y*w;		//��XjYjwj
					matrix[2]+=X*X*X*w;			//��Xj^3wj
					matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
					matrix[4]+=X*Y*Y*w;			//��XjYj^2wj
	
					matrix[6]+=Y*Y*w;			//��Yj^2wj
					matrix[9]+=Y*Y*Y*w;			//��Yj^3wj
	
					matrix[12]+=X*X*X*X*w;			//��Xj^4wj
					matrix[13]+=X*X*X*Y*w;			//��Xj^3Yjwj
					matrix[14]+=X*X*Y*Y*w;			//��Xj^2Yj^2wj
		
					matrix[19]+=X*Y*Y*Y*w;			//��XjYj^3wj
	
					matrix[24]+=Y*Y*Y*Y*w;			//��Yj^4wj
				}*/
			}

			matrix[5]=matrix[1];		//��XjYjwj
			matrix[7]=matrix[3];		//��Xj^2Yjwj
			matrix[8]=matrix[4];		//��XjYj^2wj
			matrix[10]=matrix[2];		//��Xj^3Yjwj
			matrix[11]=matrix[3];		//��Xj^2Yjwj
			matrix[15]=matrix[3];		//��Xj^2Yjwj
			matrix[16]=matrix[4];		//��XjYj^2wj
			matrix[17]=matrix[13];		//��Xj^3Yjwj
			matrix[18]=matrix[14];		//��Xj^2Yj^2wj
			matrix[20]=matrix[4];		//��XjYj^2wj
			matrix[21]=matrix[9];		//��Yj^3wj
			matrix[22]=matrix[14];		//��Xj^2Yj^2wj
			matrix[23]=matrix[19];		//��XjYj^3wj

			//matrix�̋t�s������߂� ���܂����t�s���matrix�̒����㏑�����Ċi�[�����

			/*if(n==110)
			{
				cout<<endl;
				for(int k=0;k<N;k++)
				{
					for(int k2=0;k2<N;k2++)
					{
						cout<<matrix[k*N+k2]<<" ";
					}
					cout<<endl;
				}
			}*/
			
			calc_inverse_matrix(CON,PART, N, matrix);

			int flag=ON;

			for(int k=0;k<PART[i].N2;k++)
			{
				//if(flag==ON)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
						
					double w=1;
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
					//if(dis>1) w=1/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
					//if(PART[j].surface==ON) w=1;//�\�ʂ͂�����ƒʂ��ė~��������A�d�݂�傫���Ƃ�
				
					//�ϐ�c�̍s�~j�Ɋւ�����s��
					double valX=matrix[10]*X+matrix[11]*Y+matrix[12]*X*X+matrix[13]*X*Y+matrix[14]*Y*Y;
				
					//�ϐ�e�̍s�~j�Ɋւ�����s��
					double valY=matrix[20]*X+matrix[21]*Y+matrix[22]*X*X+matrix[23]*X*Y+matrix[24]*Y*Y;

					if(valX+valY<0) flag=OFF;

				//	if(flag==ON)
					{
						if(m<pn)
						{
							val[index]=2*w*(valX+valY);		//2�K�����l��d,e,f��2�{�������̂ł���_�ɒ���
							ind[index]=m;
							index++;
							valII+=2*w*(valX+valY);
						}
						else if(PART[j].surface==ON)
						{
							valII+=2*w*(valX+valY);//�����Œl���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^
							B[n]-=2*w*(valX+valY)*PART[j].P;		//���łɕ\�ʗ��q�ɂ̓f�B���N���^�̈��͂��������Ă���B
						}
					}
				}
			}
			if(valII<0) flag=OFF;
			if(flag==OFF)
			{
				valII=0;
				B[n]=B0;
				index=index0;
				ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					    
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						val[index]=val2;
						ind[index]=m;
						index++;
						valII+=val2;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						//double w=kernel(r2,dis);
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						valII+=val2; //������w���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^
						B[n]-=val2*PART[j].P;		//���łɕ\�ʗ��q�ɂ̓f�B���N���^�̈��͂��������Ă���B
					}
				}

			}
			val[kk]=-valII;
		}
		ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�

		delete [] matrix;
	}
	else if(CON->get_dimention()==3)//3����
	{
		int N=9;
		double *matrix=new double [N*N];	//N�~N�̌W���s��
		int index=0;
		for(int n=0;n<pn;n++)
		{
			//cout<<n<<"/"<<pn<<endl;
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//�l��ۑ�
			ind[index]=n;
			index++;
			double valII=0;						//�Ίp�����̒l
			for(int k=0;k<N*N;k++) matrix[k]=0;//������
			double B0=B[n];			//B[n]�̒l��ۑ�
			int index0=index;		//index�̒l��ۑ�
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				//double dP=(PART[j].P-minP[i]);//��fj�ɑ���

				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				//if(dis>le) w=le/(dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
				//if(PART[j].surface==ON) w=1;//�\�ʂ͂�����ƒʂ��ė~��������A�d�݂�傫���Ƃ�

				matrix[0]+=X*X*w;		//��Xj^2wj
				matrix[1]+=X*Y*w;		//��XjYjwj
				matrix[2]+=X*Z*w;		//��XjZjwj
				matrix[3]+=X*X*X*w;		//��Xj^3wj
				matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
				matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
				matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
				matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
				matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
	
				matrix[10]+=Y*Y*w;		
				matrix[11]+=Y*Z*w;		
				matrix[12]+=X*X*Y*w;		
				matrix[13]+=Y*Y*Y*w;		
				matrix[14]+=Y*Z*Z*w;
				matrix[15]+=X*Y*Y*w;
				matrix[16]+=Y*Y*Z*w;
							
				matrix[20]+=Z*Z*w;			
				matrix[23]+=Z*Z*Z*w;		
					
				matrix[30]+=X*X*X*X*w;
				matrix[31]+=X*X*Y*Y*w;
				matrix[32]+=X*X*Z*Z*w;	
				matrix[33]+=X*X*X*Y*w;	
				matrix[34]+=X*X*Y*Z*w;	
				matrix[35]+=X*X*X*Z*w;	
					
				matrix[40]+=Y*Y*Y*Y*w;
				matrix[41]+=Y*Y*Z*Z*w;
				matrix[42]+=X*Y*Y*Y*w;
				matrix[43]+=Y*Y*Y*Z*w;
				matrix[44]+=X*Y*Y*Z*w;

				matrix[50]+=Z*Z*Z*Z*w;	//6�s��
				matrix[51]+=X*Y*Z*Z*w;
				matrix[52]+=Y*Z*Z*Z*w;
				matrix[53]+=X*Z*Z*Z*w;
	
				//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\


				/*B[0]+=dP*X*w;		//a
				B[1]+=dP*Y*w;		//b
				B[2]+=dP*Z*w;		//c
				B[3]+=dP*X*X*w;		//d	(X��2�K����)
				B[4]+=dP*Y*Y*w;		//e(Y��2�K����)
				B[5]+=dP*Z*Z*w;		//f(Z��2�K����)
				B[6]+=dP*X*Y*w;		//g
				B[7]+=dP*Y*Z*w;		//h
				B[8]+=dP*X*Z*w;		//i*/
			}

			matrix[9]=matrix[1];		//��XjYjwj
			matrix[17]=matrix[7];

			matrix[18]=matrix[2];
			matrix[19]=matrix[11];
			matrix[21]=matrix[8];
			matrix[22]=matrix[16];
			matrix[24]=matrix[7];
			matrix[25]=matrix[14];
			matrix[26]=matrix[5];

			matrix[27]=matrix[3];
			matrix[28]=matrix[12];
			matrix[29]=matrix[21];

			matrix[36]=matrix[4];
			matrix[37]=matrix[13];
			matrix[38]=matrix[22];
			matrix[39]=matrix[31];

			matrix[45]=matrix[5];
			matrix[46]=matrix[14];
			matrix[47]=matrix[23];
			matrix[48]=matrix[32];
			matrix[49]=matrix[41];

			matrix[54]=matrix[6];
			matrix[55]=matrix[15];
			matrix[56]=matrix[24];
			matrix[57]=matrix[33];
			matrix[58]=matrix[42];
			matrix[59]=matrix[51];
			matrix[60]=matrix[31];
			matrix[61]=matrix[44];
			matrix[62]=matrix[34];

			matrix[63]=matrix[7];
			matrix[64]=matrix[16];
			matrix[65]=matrix[25];
			matrix[66]=matrix[34];
			matrix[67]=matrix[43];
			matrix[68]=matrix[52];
			matrix[69]=matrix[61];
			matrix[70]=matrix[41];
			matrix[71]=matrix[51];

			matrix[72]=matrix[8];
			matrix[73]=matrix[17];
			matrix[74]=matrix[26];
			matrix[75]=matrix[35];
			matrix[76]=matrix[44];
			matrix[77]=matrix[53];
			matrix[78]=matrix[62];
			matrix[79]=matrix[71];
			matrix[80]=matrix[32];

			//matrix�̋t�s������߂� ���܂����t�s���matrix�̒����㏑�����Ċi�[�����
			
			calc_inverse_matrix(CON,PART, N, matrix);

			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
						
				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				//if(dis>le) w=le/(dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
				//if(PART[j].surface==ON) w=1;//�\�ʂ͂�����ƒʂ��ė~��������A�d�݂�傫���Ƃ�

				
				//�ϐ�d�̍s�~j�Ɋւ�����s��
				double valX=matrix[27]*X+matrix[28]*Y+matrix[29]*Z+matrix[30]*X*X+matrix[31]*Y*Y+matrix[32]*Z*Z+matrix[33]*X*Y+matrix[34]*Y*Z+matrix[35]*X*Z;
				
				//�ϐ�e�̍s�~j�Ɋւ�����s��
				double valY=matrix[36]*X+matrix[37]*Y+matrix[38]*Z+matrix[39]*X*X+matrix[40]*Y*Y+matrix[41]*Z*Z+matrix[42]*X*Y+matrix[43]*Y*Z+matrix[44]*X*Z;
				
				double valZ=matrix[45]*X+matrix[46]*Y+matrix[47]*Z+matrix[48]*X*X+matrix[49]*Y*Y+matrix[50]*Z*Z+matrix[51]*X*Y+matrix[52]*Y*Z+matrix[53]*X*Z;
				
				//if(valX+valY+valZ<0) w=0;					//���ʂ́A���ɂȂ�ׂ����̂��Ǝv���B���Ȃ�l�����Ȃ�
				if(m<pn)
				{
					val[index]=2*w*(valX+valY+valZ);		//2�K�����l��d,e,f��2�{�������̂ł���_�ɒ���
					ind[index]=m;
					index++;
					valII+=2*w*(valX+valY+valZ);
				}
				else if(PART[j].surface==ON)
				{
					valII+=2*w*(valX+valY+valZ);//�����Œl���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^
					B[n]-=2*w*(valX+valY+valZ)*PART[j].P;		//���łɕ\�ʗ��q�ɂ̓f�B���N���^�̈��͂��������Ă���B
				}
			}
			int flag=ON;
			if(valII<0) flag=OFF;
			else if(valII+valII!=2*valII) flag=OFF;
			if(flag==OFF)
			{
				valII=0;
				B[n]=B0;
				index=index0;
				ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<i<<endl;
				for(int k=0;k<PART[i].N2;k++)
				{
				
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					    
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						val[index]=val2;
						ind[index]=m;
						index++;
						valII+=val2;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						//double w=kernel(r2,dis);
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						valII+=val2; //������w���v�Z�ɑg�ݍ��ނƂ��A���qj���ިظڌ^�A�g�ݍ��܂Ȃ��ƌ��z�[����ɲ�݌^
						B[n]-=val2*PART[j].P;		//���łɕ\�ʗ��q�ɂ̓f�B���N���^�̈��͂��������Ă���B
					}
				}

			}
			val[kk]=-valII;
		}
		ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�

		delete [] matrix;
	}
	ft.close();

	cout<<"�s��쐬-time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
        
	/////////////////////////////////////CG�@
    double *r=new double[pn];
	double *X=new double[pn];
	
	double *AP = new double [pn];
	double *P = new double [pn];
	/////////////////////////�����l//////////////////
    if(CON->get_initialP()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initialP()==ON)//�����l�Ƃ��Č��݂̏���^����B�����ɑ����Ȃ�
	{
		for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	//BiCGStab2_method_with_D_scale_for_sparse(CON,r, pn,X,&count,1e-12,val,ind,ptr);//�Ȃ񂩁A��Ώ̂ȍs��́A������������Ȃ菬�������Ȃ��Ƃ����Ȃ��݂����B
	BiCGStab2_method_with_D_scale_for_sparse(CON,r, pn,X,&count,1e-12,val,ind,ptr);//�Ȃ񂩁A��Ώ̂ȍs��́A������������Ȃ菬�������Ȃ��Ƃ����Ȃ��݂����B
	//MRTR(CON,r, pn,X,&count,1e-10,val,ind,ptr);
	
	if(CON->get_negativeP()==OFF)//�������l�����Ȃ��ꍇ
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			PART[i].P=X[n];
			if(PART[i].P<0) PART[i].P=0;
		}
	}
	else if(CON->get_negativeP()==ON)//�������l������ꍇ
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(2*X[n]==X[n]+X[n])PART[i].P=X[n];
			//else cout<<"!!! i="<<i<<endl;
		}
	}
	
		
    delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

    //////////////////////////////*/

	//CG_GPU(CON,PART,val,ind,ptr,pn,number,ppn,B,&count);

    delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] density;

	cout<<"������:"<<count<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//�Œᗱ�q�ԋ�����藱�q�����x�𐄎Z����B
void calc_PND_by_minmum_L_main(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number)
{
	int dimention=CON->get_dimention();
	double Re=CON->get_re2();
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].N>0)
		{
			double minL=100;//�Œᗱ�q�ԋ���
			double minL2=110;
			double minL3=110;
			double minL4=110;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(dis<minL) minL=dis;
				else if(dis<minL2) minL2=dis;
				else if(dis<minL3) minL3=dis;
				else if(dis<minL4) minL4=dis;
			}
			//minL�����܂���
			//cout<<minL<<" "<<minL2<<endl;
			double ave=(minL+minL2+minL3+minL4)/4;
			//cout<<i<<" "<<PART[i].PND2<<" "<<ave<<" "<<CON->get_distancebp()<<endl;
			double before=PART[i].PND2;
			PART[i].PND2=calc_PND_by_minmum_L(CON, dimention, Re, ave);
			//if(i==586) cout<<i<<" "<<PART[i].PND2<<" "<<before<<endl;
		}
	}
}

//�Œᗱ�q�ԋ�����藱�q�����x�𐄎Z����B
double calc_PND_by_minmum_L(mpsconfig *CON,int dimention,double Re,double L)
{
	//int size = (int)(r+1);//�v�Z�̈�
	int size = (int)(Re+2);//�v�Z�̈�
	int calc_type=CON->get_model_set_way();
	double le=CON->get_distancebp();
	double dis;//����
	double pnd=0;
	if(dimention==2)
	{
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis*L<=Re*le )
					{
						pnd+=kernel(Re*le,dis*L);
					}			
				}
			}
		}
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis*L<=Re*le )
					{
						pnd+=kernel(Re*le,dis*L);
					}			
				}
			}
		}
	}
	else if(dimention==3)
	{
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis*L<=Re*le )
						{
							pnd+=kernel(Re*le,dis*L);
						}
					}			
				}
			}
		}
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis*L<=Re*le )
						{
							pnd+=kernel(Re*le,dis*L);
						}
					}			
				}
			}
		}
	}
	return pnd;  
}

//�t�s������߂�֐�
void calc_inverse_matrix(mpsconfig *CON,vector<mpsparticle> &PART,int N, double *matrix)
{
	//N:���m��
	double buf;
	int i,j,k;


	double **inv_a=new  double*[N];					//�t�s��i�[
	for(int n=0;n<N;n++) inv_a[n]=new double [N]; 
	double **a=new  double*[N];						//matrix�̒l���R�s�[
	for(int n=0;n<N;n++) a[n]=new double [N]; 

	int count=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			a[i][j]=matrix[count];
			count++;
		}
	}
	//�P�ʍs������
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(i==j) inv_a[i][j]=1;
			else		inv_a[i][j]=0;
		}
	}
	//�|���o���@
	for(i=0;i<N;i++)
	{
		buf=1/a[i][i];
		for(j=0;j<N;j++)
		{
			a[i][j]*=buf;
			inv_a[i][j]*=buf;
		}
		for(j=0;j<N;j++)
		{
			if(i!=j)
			{
				buf=a[j][i];
				for(k=0;k<N;k++)
				{
					a[j][k]-=a[i][k]*buf;
					inv_a[j][k]-=inv_a[i][k]*buf;
				}
			}
		}
	}

	count=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			matrix[count]=inv_a[i][j];
			count++;
		}
	}


	for(int n=0;n<N;n++) delete [] inv_a[n];
	delete [] inv_a;
	for(int n=0;n<N;n++) delete [] a[n];
	delete [] a;
}

//MRTR�@//
void MRTR(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//�덷
	double BB=0;
	double rr=0;
		
	double zeta;
	double zeta_old;
	double eta;
	double v;

	double Ar_r;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)

    //double *r= new double[pn];
	double *Ar = new double [pn];  //   _ _
	
	double *P = new double [pn];
	double *y = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0.0;
	//v=complex<double>(1.0,0.0);//1��ڂ̔����ɂ����āA0�łȂ����Ƃɒ���
	v=1.0;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //�����l
	{
		X[n]=0.0;
		//r[n]=B[n];
		P[n]=r[n];
		y[n]=0.0;
	}

	for(int n=0;n<pn;n++) BB+=r[n]*r[n];
	//BB=sqrt(BB);
	//cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
	
	cout<<"MRTR�@�X�^�[�g-----���m��="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//��������(convergence test)
		//while(E>EP)// EP=CON->get_CGep();//��������(convergence test)
		while(E>EP && count<=pn)// EP=CON->get_CGep();//��������(convergence test)
		{
			//if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c�v�Z
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
				}
			}

			Ar_r=0.0;//cAr_r: (cAr,r)(����)���w���B(u,v)=��((cu)*v)
			Ar_Ar=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			y_Ar=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
			Ar_y=0.0;//cAr_Ar: (cAr,Ar)(����)���w���B(u,v)=��((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				Ar_r+=Ar[n]*r[n];
				Ar_Ar+=Ar[n]*Ar[n];
				y_Ar+=y[n]*Ar[n];
				Ar_y+=Ar[n]*y[n];
			}

			zeta=v*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
			eta=-y_Ar*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);

			//////v�̌v�Z//////////////
			v=0.0;
			for(int n=0;n<pn;n++) v+=zeta*Ar[n]*r[n];

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
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"������="<<count<<" E="<<E<<endl;
	
			
		}
	}
	*countN=count;//�����񐔂�n��
	//cout<<"�I�� ������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	//delete [] r;
	delete [] Ar;
	
	delete [] P;
	delete [] y;
	c.close();
}

//�a�s��p�̑Ίp�X�P�[�����O��bicgstab2�@
void BiCGStab2_method_with_D_scale_for_sparse(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr)
{
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double pk=1;
	double E=1;//�덷
	double alp,beta,rr,w,ita;

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
	double *D=new double [pn];			//D^(-1)
	beta=0;
	w=0;
	ita=0;

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
		for(int m=ptr[n];m<ptr[n+1];m++)
		{
			if(ind[m]==n)
			{
				if(val[m]==0) cout<<"RR"<<endl;
				D[n]=1.0/val[m];
			}
		}
		
	}

	ofstream c("convergenceP.dat");
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"BiCGstab2�@:���m��="<<pn<<" ---";
	while(E>EP)// EP=CON->get_CGep();//��������(convergence test)
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
	//	#pragma omp parallel for
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]]*D[ind[m]];	//new
			/*for(int m=0;m<pn;m++)
			{
				AP[n]+=A[n][m]*P[m]*D[m];
			}*/
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

	//	#pragma omp parallel for
		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]]*D[ind[m]];//new
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m]*D[m];
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
			X[n]+=alp*P[n]*D[n]+Z[n]*D[n];
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
		E=sqrt(rr/rr0);
		c<<count<<" "<<E<<endl;
		//cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	*countN=count;//�����񐔂�n��
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

	delete [] D;
	c.close();
}

//�ŏ����@�ɂ��\�ʑ��x��`�֐�
void reset_surface_velocity_by_WLSM(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number)
{
	//�P���ɓ����v�Z�_�Ɋւ��Ĉ��͂̕������������������ł́A�\�ʌv�Z�_�̑��x���U�̓[���ɂȂ�Ȃ��B����āA�\�ʂ̑��x�͉����B�����ŁA�����������q�Ɋւ���ŏ����@�ɂ��⊮����

	double le=CON->get_distancebp();
	
	int d=CON->get_dimention();
	int N0=0;					//�W���s��̌�
	int order=1;//CON->get_divU_order();				//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
    
	//�W���s��̑傫���̌���
	if(d==2)
	{
		if(order==1) N0=3;
		else if(order==2) N0=6;
		else if(order==3) N0=10;
		N0=9+1;	//�ő�l
	}
	else if(d==3)
	{
		if(order==1) N0=4;
		else if(order==2) N0=10;
		else if(order==3) N0=20;
		N0=20;
	}
	////////////////////////////////

	double *matrix=new double [N0*N0];	//N�~N�̌W���s��
	double *matrix2=new double [N0*N0];	//matrix�̃o�b�N�A�b�v
	double *B1=new double [N0];			//N�̉��s��
	double *B2=new double [N0];			//N�̉��s��
	double *B3=new double [N0];			//N�̉��s��
	double **MAT= new double*[N0];		//N�~N�̌W���s��(�z���2����)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//���x�N�g���i�[

	if(d==2)				//�񎟌�
	{
		for(int i=0;i<fluid_number;i++)
		{
			int N=N0;
			int jnb=0;			//���ӓ������q��
			int J=i;
			double dis0=100;
			if(PART[i].surface==ON) 
			{
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID && PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double dis=sqrt(X*X+Y*Y);
						if(dis<dis0)
						{
							dis0=dis; J=j;
						}
					}
				}
				jnb=0;
				//J=i;
				for(int k=0;k<PART[J].N;k++)
				{
					int j=PART[J].NEI[k]; 
					if(PART[j].surface==OFF) jnb++;
					//if(PART[j].type==FLUID && PART[j].surface==OFF) jnb++;
				}
			}


			if(jnb>=3 &&  PART[i].surface==ON) //���ӗ��q�������菭�Ȃ��ꍇ�͑��x��Ԃ��s��Ȃ�
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;	//������
				for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0;}
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//������

				//if(PART[i].N>=15) {order=3; N=10;}
				if(jnb>=15) {order=2; N=6;}		//�\�ʂ�����Q���ߎ��ŏ\���ȋC�͂���
				else if(jnb>=6) {order=2; N=6;}
				else {order=1; N=3;}
				double R=PART[J].L*CON->get_re3();
				for(int k=0;k<PART[J].N3;k++)
				{
					int j=PART[J].NEI3[k];
					//if(PART[j].type==FLUID && PART[j].surface==OFF)
					if(PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[J].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[J].r[A_Y]);
						double U=(PART[j].u[A_X]);
						double V=(PART[j].u[A_Y]);
						double dis=sqrt(X*X+Y*Y);
					
						double w=kernel_in_WLSM(dis,R);	//�X�v���C��
						if(dis>R) w=0;

						if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//���x�N�g��
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y;  base[5]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y;  base[9]=1;}

						//�s��쐬
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B1[n]+=base[n]*w*U;
							B2[n]+=base[n]*w*V;
						}
					}
				}
				if(J!=i)		//���qJ�̊�^�B�������qJ��i�ɓ������Ȃ�l������K�v�͂Ȃ�
				{
					MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
					B1[N-1]+=PART[J].u[A_X];	//���qi�̊�^�B
					B2[N-1]+=PART[J].u[A_Y];	//���qi�̊�^�B
					//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�*/
				}
				
				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��

				for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

				gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

				//���x�̒l����
				double X=(PART[i].r[A_X]-PART[J].r[A_X]);
				double Y=(PART[i].r[A_Y]-PART[J].r[A_Y]);
				if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//���x�N�g��
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y;  base[5]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y;  base[9]=1;}

				for(int D=0;D<d;D++) PART[i].u[D]=0;
				for(int n=0;n<N;n++)
				{
					PART[i].u[A_X]+=B1[n]*base[n];
					PART[i].u[A_Y]+=B2[n]*base[n];
				}
			}
		}
	}
	else if(d==3)				//3����
	{
		for(int i=0;i<fluid_number;i++)
		{
			int N=N0;
			int jnb=0;			//���ӓ������q��
			int J=i;			//�ŋߐړ������̗��q
			double dis0=100;
			if(PART[i].surface==ON) 
			{
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(PART[j].type==FLUID && PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<dis0) {dis0=dis; J=j;}
					}
				}
				jnb=0;
				for(int k=0;k<PART[J].N;k++)
				{
					int j=PART[J].NEI[k]; 
					if(PART[j].surface==OFF)jnb++;
					//if(PART[j].type==FLUID && PART[j].surface==OFF) jnb++;
				}
			}


			if(jnb>=4 &&  PART[i].surface==ON) //���ӗ��q�������菭�Ȃ��ꍇ�͑��x��Ԃ��s��Ȃ�
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;	//������
				for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0; B3[n]=0;}
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//������

				if(jnb>=25) {order=3; N=20;}
				else if(jnb>=10) {order=2; N=10;}
				else {order=1; N=4;}
				double R=PART[J].L*CON->get_re3();

				for(int k=0;k<PART[J].N3;k++)
				{
					int j=PART[J].NEI3[k];
					//if(PART[j].type==FLUID && PART[j].surface==OFF)
					if(PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[J].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[J].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[J].r[A_Z]);
						double U=(PART[j].u[A_X]);
						double V=(PART[j].u[A_Y]);
						double W=(PART[j].u[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=kernel_in_WLSM(dis,R);	//�X�v���C��
						if(dis>R) w=0;								//�X�v���C��

						if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//���x�N�g��
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
						//�s��쐬
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B1[n]+=base[n]*w*U;
							B2[n]+=base[n]*w*V;
							B3[n]+=base[n]*w*W;
						}
					}
				}
				if(J!=i)		//���qJ�̊�^�B�������qJ��i�ɓ������Ȃ�l������K�v�͂Ȃ�
				{
					MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
					B1[N-1]+=PART[J].u[A_X];	//���qi�̊�^�B
					B2[N-1]+=PART[J].u[A_Y];	//���qi�̊�^�B
					B3[N-1]+=PART[J].u[A_Z];	//���qi�̊�^�B
					//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�*/
				}

				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��

				for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

				gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B3,N);//�K�E�X�̏����@�ŉ���

				//���x�̒l����
				double X=(PART[i].r[A_X]-PART[J].r[A_X]);
				double Y=(PART[i].r[A_Y]-PART[J].r[A_Y]);
				double Z=(PART[i].r[A_Z]-PART[J].r[A_Z]);
				if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//���x�N�g��
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
				for(int D=0;D<d;D++) PART[i].u[D]=0;
				for(int n=0;n<N;n++)
				{
					PART[i].u[A_X]+=B1[n]*base[n];
					PART[i].u[A_Y]+=B2[n]*base[n];
					PART[i].u[A_Z]+=B3[n]*base[n];
				}
			}
		}
	}
	
	delete [] matrix;
	delete [] matrix2;
	delete [] B1;
	delete [] B2;
	delete [] B3;
	for(int n=0;n<N0;n++) delete [] MAT[n];
	delete [] MAT;
	delete [] base;
}

