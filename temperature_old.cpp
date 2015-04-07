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

///���x��v�Z�֐�
void calc_Temperature(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	cout<<"���x��v�Z------";

	//�����łʹ����߰����̂ɂ��Čv�Z����B���Ȃ킿�A
	//Dh/Dt=k��T+Q k:�M�`����
	//���Ȃ݂Ɏ��ό`�����
	//DT/Dt=k/(��C)��T+Q/(��C)�ƂȂ�

	unsigned int timeA=GetTickCount();					//�v�Z�J�n����
	int dim=CON->get_dimention();						//��͎���
	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;							//���v���V�A���e�����a

	double V=get_volume(CON);				//���q�̑̐�
	double *density=new double [particle_number];//�e���q�̖��x�i�[
	double *Cp=new double [particle_number];		//��M
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)
		{
			density[i]=CON->get_density();
			Cp[i]=CON->get_Cp();
		}
		else if(PART[i].materialID==2) 
		{
			density[i]=CON->get_density2();
			Cp[i]=CON->get_Cp2();
		}
	}
	for(int i=fluid_number;i<particle_number;i++)
	{
		density[i]=CON->get_wall_density();
		Cp[i]=CON->get_wall_Cp();
	}
	double *mass=new double [fluid_number];	//�e���q�̎��ʊi�[
	for(int i=0;i<fluid_number;i++) mass[i]=density[i]*V;
		
	double *MP=new double [fluid_number];		//�Z�_
	double *latent_H=new double [fluid_number];		//���M
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)
		{	
			MP[i]=CON->get_MP();
			latent_H[i]=CON->get_latent_H();
		}
		else if(PART[i].materialID==2)
		{
			MP[i]=CON->get_MP2();
			latent_H[i]=CON->get_latent_H2();
		}
	}

	double roomT=CON->get_roomT();						//����[K]
	double wall_mass=V*CON->get_wall_density();			//�Ǘ��q�̎���

    double *T_laplacian = new double [particle_number];
    double *T=new double [particle_number];//���x
    double *k=new double [particle_number];//�e���q�̔M�`����
   
    ///�G���^���s�[����e���q�̉��xT[i]�����߂�
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+latent_H[i]*mass[i];			//�Z���I���_�̃G���^���s�[

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//�ő�
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];	//�Z�_
		else if(hs1<PART[i].h)											//�t��
		{       
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if( PART[i].type==SOLID ||PART[i].type==BOSOLID ||PART[i].type>CFD)//���ϑ�
			{
				PART[i].type=FLUID;
				if(PART[i].PND>n0*CON->get_beta()) PART[i].surface=OFF;
				else PART[i].surface=ON;
			}  
		}
    }
	for(int i=fluid_number;i<particle_number;i++)//�Ǘ��q
    {
        T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//�ő�
    }
    ///////////////
	
    ///////�e���q�̔M�`����k[i]�����߂�
	double maxk=0;//�ő�M�`����
	double densityCp=CON->get_density()*CON->get_Cp();
    for(int i=0;i<fluid_number;i++)
    {
		if(PART[i].materialID==1) k[i]=CON->get_k();
		else if(PART[i].materialID==2) k[i]=CON->get_k2();
		if(k[i]>maxk) maxk=k[i];
    }
	for(int i=fluid_number;i<particle_number;i++)
    {
		k[i]=CON->get_wall_k();//�ǂ̔M�`����
		if(k[i]>maxk) maxk=k[i];
    }
	maxk/=densityCp;//�M�g�U��
	if(maxk*dt/(le*le)>0.2) cout<<"�M�`�����Ɋւ���dt���傫������ dt<"<<0.2*le*le/maxk<<"�ɂ��Ă�������"<<endl;

    ////���x�̊g�U�v�Z
    if(CON->get_insulate()==0)//�f�M
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//������
			if(PART[i].type!=INWALL && PART[i].type!=OUTWALL)
			{
				double lam=0;//���m�ȃ�
				double W=0;//���q�����x
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
					if(PART[j].type!=INWALL && PART[j].type!=OUTWALL)
					{ 
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
						double w=kernel(R,dis);//�d�݊֐�
						W+=w;
						lam+=dis*dis*w;
						double dT=T[j]-T[i];//���x��
						if(CON->get_T_laplacian()==0 || CON->get_T_laplacian()==1)
						{
							T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;
						}
						else if(CON->get_T_laplacian()==2)
						{ 
							T_laplacian[i]+=dT*w/(dis*dis)*k[i]*k[j]/(k[i]+k[j])*2;
						}  
					}
				} 
	    
				if(W!=0) lam/=W;
				else if(W==0) lam=lamda;
				int d=CON->get_dimention();
				if(CON->get_T_laplacian()==0)  T_laplacian[i]=T_laplacian[i]*2*d/(n0*lamda);
				else if(CON->get_T_laplacian()==1 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/(W*lam);
				else if(CON->get_T_laplacian()==2 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/W;
			}
		}
    }
    else if(CON->get_insulate()==1)//��f�M
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//������
			//if(PART[i].type!=INWALL && PART[i].type!=OUTWALL)
			{
				double lam=0;//���m�ȃ�
				double W=0;//���q�����x
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
		
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
					double w=kernel(R,dis);//�d�݊֐�
					W+=w;
					lam+=dis*dis*w;
					double dT=T[j]-T[i];//���x��

					if(CON->get_T_laplacian()==0 || CON->get_T_laplacian()==1)
					{
						T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;///(density[i]*Cp[i]);
					}
					else if(CON->get_T_laplacian()==2)
					{ 
						T_laplacian[i]+=dT*w/(dis*dis)*k[i]*k[j]/(k[i]+k[j])*2;
					}
				}
				if(W!=0) lam/=W;
				else if(W==0) lam=lamda;
				int d=CON->get_dimention();
				if(CON->get_T_laplacian()==0)  T_laplacian[i]=T_laplacian[i]*2*d/(n0*lamda);
				else if(CON->get_T_laplacian()==1 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/(W*lam);
				else if(CON->get_T_laplacian()==2 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/W;
				
			}
		}
    }
    //////////////*/
  
    for(int i=0;i<particle_number;i++) PART[i].h+=T_laplacian[i]*dt*V;//�G���^���s�[�X�V

	if(CON->get_air_cool()==ON)//��C�Ƃ̔M�`�B���l������ꍇ
	{
		double h=50;			//��C�̔M�`�B���@�P�ʂ�[W/m^2/K]. �l�͂�����ƒ��ׂ邱�ƁB�k�b�Z�����Ƃ��Ɗ֌W����H
		double S;//���q�̒f�ʐ�
		if(dim==2) S=le;
		else
		{
			if(CON->get_model_set_way()==0) S=le*le;//�����i�q�̏ꍇ
			else if(CON->get_model_set_way()==1) S=sqrt(3.0)/2*le*le;//�ז��i�q�̏ꍇ
		}

		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON) PART[i].h+=h*(roomT-T[i])*S*dt;
		}
	}//////////*/
   
    
    ////���x�v���b�g
	double height=CON->get_TplotZ();
    plot_T(CON ,PART,particle_number,T,height);
	if(CON->get_T_AVS()>0)
	{
		if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs(CON,PART, t, particle_number, fluid_number,T, height);
	}


	if(CON->get_dimention()==3)//XZ���ʂ�T���o��
	{
		ofstream fh("T_XZ.dat");
		for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fh<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<T[i]<<endl;
		fh.close();
	}
    
    delete [] T;
    delete [] k;
    delete [] T_laplacian;

	delete [] density;
	delete [] mass;
	delete [] Cp;
	delete [] MP;
	delete [] latent_H;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//���x�v���b�g�֐�
void plot_T(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,double *T,double height)
{
	ofstream t("T_XY.dat");
	double le=CON->get_distancebp();
	double L=le*0.5;

	if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			double x=PART[i].r[A_X];
			double y=PART[i].r[A_Y];
			double Z=T[i];
			t<<x<<" "<<y<<" "<<Z<<endl;
		}
	}
	else if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			double z=PART[i].r[A_Z];
			if(z>height-L && z<height+L)
			{
				double x=PART[i].r[A_X];
				double y=PART[i].r[A_Y];
				double Z=T[i];
				t<<x<<" "<<y<<" "<<Z<<endl;
			}
		}
	}

	t.close();
}

//���xAVS�t�@�C���o�͊֐�
void output_temperature_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number,double *T,double height)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	t=1;//���܂͂킴�Ɩ��X�e�b�v�㏑��

	//sprintf(filename,"pressure/pressure%d",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf(filename,"T_XZ%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
				double P=T[i];
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
			double P=T[i];
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
	}
	fout.close();
	//sprintf(filename,"pressure/pressure%d.fld",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf(filename,"T_XZ%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
	fout2 << "label=temperature" << endl << endl;
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	fout2 << "variable 1 file=T_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    1 file=T_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    2 file=T_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    3 file=T_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2.close();
}
