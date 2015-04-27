#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"	//#define �i�[
//#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
//#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
//#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
//#include<omp.h>
//#include<vector>
#include"function.h"

///���x��v�Z�֐�
void calc_Temperature(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	int counth=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].heat_generation>0) counth++;
	}
	//cout<<"heat��0�łȂ����q��="<<counth<<endl;
	//for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) PART[i].heat_generation=10;
	cout<<"���x��v�Z------";

	//�����łʹ����߰����̂ɂ��Čv�Z����B���Ȃ킿�A
	//Dh/Dt=k��T+Q k:�M�`����
	//���Ȃ݂Ɏ��ό`�����
	//DT/Dt=k/(��C)��T+Q/(��C)�ƂȂ�

	unsigned int timeA=GetTickCount();					//�v�Z�J�n����
	int dim=CON->get_dimention();						//��͎���
	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;							//���v���V�A���e�����a

	double V=CON->get_particle_volume();			//���q�̑̐�
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


	
	//���M�ʂ��v�Z����
	for(int i=0;i<particle_number;i++)
	{
		PART[i].h+=0.5*dt*(PART[i].heat_generation+PART[i].heat_gene_before1)*V;		//1�X�e�b�v�O�̔��M�ƍ��킹�đ�`���ŋߎ�
		//PART[i].h+=0.5*dt*(PART[i].heat_generation+PART[i].heat_gene_before1)*1;//*V;		//1�X�e�b�v�O�̔��M�ƍ��킹�đ�`���ŋߎ�
	}

    double *T_laplacian = new double [particle_number];
    double *T=new double [particle_number];//���x
	int *insulate=new int [particle_number];//�Ǘ��q���f�M����Ă��邩�ǂ����@�Ǐ��I�ɒf�M��Ԃɂ������Ƃ��Ȃǂɗp����
    double *k=new double [particle_number];//�e���q�̔M�`����

	for(int i=0;i<particle_number;i++) insulate[i]=OFF;
	
   
    ///�G���^���s�[����e���q�̉��xT[i]�����߂�
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+latent_H[i]*mass[i];			//�Z���I���_�̃G���^���s�[

		if(PART[i].h<hs0)
		{
			//cout<<"�ő�"<<endl;
			T[i]=PART[i].h/mass[i]/Cp[i];			//�ő�
		}
		else if(hs0<=PART[i].h && PART[i].h<hs1)
		{
			//cout<<"�Z��?"<<endl;
			T[i]=MP[i];	//�Z�_
			
		}
		else if(hs1<=PART[i].h)											//�t��
		{
			//cout<<"�t��"<<endl;
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if(CON->get_material()==H2O && T[i]>=393) T[i]=393;//���ݕ����͍l���Ȃ��̂ŁA�������x��100�x�z�������100�x�ɖ߂�(�C���M�̎����̕����挈��?�j
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
		if(CON->get_model_number()==26)//���̏ꍇ�A�ǂ̏ꏊ�ɉ����ĕ����n��������
		{
			double air_mass=V*1.205;//��C�̏d��
			double air_Cp=1006;//��C�̔�M
			if(CON->get_dimention()==2)
			{
				if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//���
				{
					T[i]=PART[i].h/air_mass/air_Cp;//�ő�
				}
				else T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//�ő�
			}
			else if(CON->get_dimention()==3)
			{
			}
		}
        else T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//�ő�
    }
    ///////////////

	////fsw�ɂ����āA�ꕔ�ȊO�̕ǖʂ̉��x������ς���
	if(CON->get_model_number()==19)
	{
		double width0=18*1e-3;//25*1e-3;//18*1e-3;//���̂���߂镝
		double Width=width0;		//��18mm+�Ǘ��q�����E��4���q��
		double Height=6*1e-3;	//����6mm+�Ǘ��q������4���q��
		double Depth=30*1e-3;//27*1e-3+6*le*A;	//���s��27mm+�Ǘ��q�����E��4���q��
		double depth0=9e-3;				//�c�[�����S�ƁA��O�̕ǂƂ̋���
	
		//double Z_mod=le*B+(1-sqrt(2.0)/2)*le;//�c�[�������̂ւ߂荞�ނ̂�h�����߂̒����ʁ@���̒l�����c�[���ȊO�̕��̂�������
		double Z_mod=le+(1-sqrt(2.0)/2)*le;
			
		///////////////////box���q�̍ގ�����

		for(int i=fluid_number;i<particle_number;i++) 

		{
			if(PART[i].r[A_X]>(Width*0.5)+0.2*le)//�E�̕�
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
			else if(PART[i].r[A_X]<-(Width*0.5)-0.2*le)//���̕�
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
				
			else if(PART[i].r[A_Y]<-(depth0)-0.2*le)//��O�̕�
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
			else if(PART[i].r[A_Y]>(Depth-depth0)+0.2*le)//���̕�
			{
					T[i]=CON->get_roomT();
					insulate[i]=OFF;
			}
			if(PART[i].r[A_X]<=(Width*0.5)+0.2*le && PART[i].r[A_X]>=-(Width*0.5)-0.2*le)
			{
				if(PART[i].r[A_Y]>=-(depth0)-0.2*le && PART[i].r[A_Y]<=(Depth-depth0)+0.2*le)//���̕�
				{
					if(PART[i].r[A_Z]<-0.2*le-Z_mod)
					{
						insulate[i]=ON;	
					}
				}
			}
		}
	}

	//���E�̉��x���̂��̂������ȂǂŊ��ɂ��Ƃ܂��Ă���ꍇ�A�����Ή����鋫�E�Ɉʒu���闱�q�ɗ^����
	//set_temperature_boundary(CON, PART, fluid_number, particle_number, T,dt,t);
	
	//////////////////////*/
	
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
		if(CON->get_model_number()==26)//���̏ꍇ�A�ǂ̏ꏊ�ɉ����ĕ����n��������
		{
			if(CON->get_dimention()==2)
			{
				if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//���
				{
					k[i]=0.0257;//��C�̔M�`����
				}
				else k[i]=CON->get_wall_k();//�ǂ̔M�`����
			}
			else if(CON->get_dimention()==3)
			{
			}
		}
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
	else if(CON->get_insulate()==2)//�Ǐ��I�ɒf�M
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//������
			if(insulate[i]==OFF)
			{
				double lam=0;//���m�ȃ�
				double W=0;//���q�����x
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
					if(insulate[j]==OFF)
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
							T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;///(density[i]*Cp[i]);
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
    //////////////*/
  
	//�G���^���s�[�X�V
    for(int i=0;i<particle_number;i++) PART[i].h+=T_laplacian[i]*dt*V;

	if(CON->get_air_cool()==ON)//��C�Ƃ̔M�`�B���l������ꍇ
	{
		double h=5;//50;			//��C�̔M�`�B���@�P�ʂ�[W/m^2/K]. �l�͂�����ƒ��ׂ邱�ƁB�k�b�Z�����Ƃ��Ɗ֌W����H
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
			else if(CON->get_model_number()==26 && PART[i].heat_generation>0) PART[i].h+=h*(roomT-T[i])*S*dt;//���ɂ����āA���M�ʂ����݂���(=���̍ŊO��)�����M������Ƃ���
		}
	}//////////*/
    
	double Tmax=0;
	for(int i=0;i<particle_number;i++) if(T[i]>Tmax) Tmax=T[i];
	double Tminw=100000;
	for(int i=0;i<particle_number;i++) if(T[i]<Tminw) if(PART[i].type!=FLUID) Tminw=T[i];
	//cout<<"�ő剷�x="<<Tmax<<"[K]"<<endl;
	double Tmaxf=0;
	for(int i=0;i<fluid_number;i++) if(T[i]>Tmaxf) Tmaxf=T[i];
	cout<<"�ő剷�x="<<Tmax<<"[K]"<<" ���̍ő剷�x"<<Tmaxf<<"[K]"<<endl;

	//���M�ʂ�������
	for(int i=0;i<particle_number;i++)
	{
		PART[i].heat_gene_before1=PART[i].heat_generation;	//1step�O�̏��ւƊi�[���Ȃ���
		PART[i].heat_generation=0;							//�������@���M�ʂ𑝂₷�Ƃ���+=�̋L�q�ɂ���ƁA�����Ȋ֐����Ŕ��M���N�������Ƃ��ɂ��Ώ��ł���//���݂͉Q�d�����݂̂ŁA���X�e�b�v���߂Ă�
	}
	
	//���q�։��x���L��
	for(int i=0;i<particle_number;i++)
	{
		//if(PART[i].type==FLUID) PART[i].T=T[i];
		 PART[i].T=T[i];
		//if(PART[i].T<CON->get_initialT()) cout<<"�������x���Ⴂ?"<<endl; 
	}

	////���x�v���b�g
	double height=CON->get_TplotZ();
    plot_T(CON ,PART,particle_number,T,height);
	if(CON->get_T_AVS()>0)
	{
		if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs(CON,PART, t, particle_number, fluid_number,T, height);
	}

	///microavs�p_�S���o��
	
	if(t==1 || t%10==0)
	{
		int n=0;
		char filename[30];
		sprintf_s(filename,"PART.T%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
		ofstream fout5(filename);
		if(!fout5)
		{
			cout << "cannot open" << filename << endl;
			exit(EXIT_FAILURE);
		}
		if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type!=FLUID)
				{
					double x=PART[i].r[A_X];	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
					double y=PART[i].r[A_Y];//*1.0E+05;
					double z=PART[i].r[A_Z];//*1.0E+05;
					double T=PART[i].T;
					fout5 << T << "\t" << x << "\t" << y << "\t" << z << endl;
					//fout5 << insulate[i] << "\t" << x << "\t" << y << "\t" << z << endl; ���������ƁA�f�M�A��f�M�̏����ݒ�������ł���
					n++;
				}
			}
		}
		fout5.close();
		
		sprintf_s(filename,"PART.T%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
		fout6 << "label=T" << endl << endl;
		fout6 << "variable 1 file=PART.T" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
		fout6 << "coord    1 file=PART.T" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
		fout6 << "coord    2 file=PART.T" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
		fout6 << "coord    3 file=PART.T" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
		fout6.close();
	}
	
	if(CON->get_dimention()==3)//XZ���ʂ�T���o��
	{
		
		ofstream fh("T_XZ.dat");
		for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fh<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<T[i]<<endl;
		fh.close();
	}

    
    delete [] T;
	delete [] insulate;
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
	double cross_section=CON->get_speed_face_p();
	//t=1;//���܂͂킴�Ɩ��X�e�b�v�㏑��

	//sprintf_s(filename,"pressure/pressure%d",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"T_XZ%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
			if(PART[i].type==FLUID)
			{
				if(PART[i].r[A_Y]<cross_section+0.5*le && PART[i].r[A_Y]>cross_section-0.5*le)	
				//if(PART[i].r[A_Y]<0.006+0.5*le && PART[i].r[A_Y]>0.006-0.5*le)	
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
	}
	else if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=FLUID)
			{
			double x=PART[i].r[A_X]*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
			double y=PART[i].r[A_Y]*1.0E+05;
			double z=PART[i].r[A_Z]*1.0E+05;

			//double x=PART[i].r[A_X];//*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
			//double y=PART[i].r[A_Y];//*1.0E+05;
			//double z=PART[i].r[A_Z];//*1.0E+05;
			double P=T[i];
			//double P=PART[i].heat_gene_before1;
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
		}
	}
	fout.close();
	//sprintf_s(filename,"pressure/pressure%d.fld",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"T_XZ%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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


///���x��A�I�v�Z�֐�
void calc_temperature_implicity(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	//���݂�3�c�̂ݑΉ�
	//�����ׂ�����	hi(t+1)-hi(t)=��t*k*(2d/��n0)��(hj-hi)w+q*��t
	//				�̃�t*k*(2d/��n0)��(hj-hi)w-hi(t+1)=-hi(t)-q*��t
	//				�̃�(hj-hi)w-��n0/(2kd��t)hi(t+1)=-��n0/(2kd��t)hi(t)-��n0/(2kd)q
	//				��k��(hj-hi)w-��n0/(2d��t)hi(t+1)=-��n0/(2d��t)hi(t)-��n0/(2d)q

	cout<<"���x��v�Z(�A�I)--";

	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;						//���׼�ݗp�e�����a
	int d=CON->get_dimention();						//����
	unsigned int timeA=GetTickCount();				//�v�Z�J�n����
	int count=0;
	
	double co=lamda*n0/(2*d*dt);					//�v�Z�ɂ悭�����W��
	double co2=lamda*n0/(2*d);

	double V=CON->get_particle_volume();
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
		else cout<<"���q"<<i<<"�̍ގ�ID���s��?"<<endl;
	}
	for(int i=fluid_number;i<particle_number;i++)
	{
		density[i]=CON->get_wall_density();
		Cp[i]=CON->get_wall_Cp();
	}
	double *mass=new double [particle_number];	//�e���q�̎��ʊi�[
	for(int i=0;i<particle_number;i++) mass[i]=density[i]*V;
		
	
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

   
    double *T=new double [particle_number];//���x
    double *K=new double [particle_number];//�e���q�̔M�`����

   
    ///�G���^���s�[����e���q�̉��xT[i]�����߂�
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+latent_H[i]*mass[i];			//�Z���I���_�̃G���^���s�[

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//�ő�
		else if(hs0<=PART[i].h && PART[i].h<=hs1)
		{
			T[i]=MP[i];	//�Z�_
			//cout<<"�n�Z���̗��q������"<<endl;
		}
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
        T[i]=PART[i].h/mass[i]/Cp[i];//�ő�
    }
    ///////////////
	
    ///////�e���q�̔M�`����k[i]�����߂�
    for(int i=0;i<fluid_number;i++)
    {
		double densityCp=density[i]*Cp[i];
		if(PART[i].materialID==1) K[i]=CON->get_k()/densityCp;
		else if(PART[i].materialID==2) K[i]=CON->get_k2()/densityCp;
    }
	for(int i=fluid_number;i<particle_number;i++)
    {
		double densityCp=density[i]*Cp[i];
		K[i]=CON->get_wall_k()/densityCp;//�ǂ̔M�`����
    }

	//�����_step(k)�ł̉��x���v���V�A�����v�Z�@�����ߎ�����ۂɎg�p����
	double *T_lap_before1=new double [particle_number];
	for(int i=0;i<particle_number;i++)
	{
		double lam=0;					//�e���q�̐��m�ȃ�
		double W=0;						//�e���q�̐��m�ȗ��q�����x
		T_lap_before1[i]=0.0;			//������
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k]; 
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
					
			double w=kernel(R,dis);//�d�݊֐�
			W+=w;
			lam+=dis*dis*w;
			double dT=T[j]-T[i];
			//T_lap_before1[i]+=dT*w;
			T_lap_before1[i]+=dT*w*K[i]*K[j]/(K[i]+K[j])*2;
		}
		T_lap_before1[i]=T_lap_before1[i]*2*d/(n0*lamda);
	}
   ///////*/

	

	int pn=0;	//���m��
	for(int i=0;i<particle_number;i++) if(PART[i].surface==OFF) pn++; 

	int *ppn = new int[pn];					//�s��ɂ������n�Ԗڂ̖��m���͗��q�ԍ�ppn[n]�̗��q�ɑ���
	int *link = new int [particle_number];	//���q�ԍ�i��link[i]�Ԗڂ̖��m��

	double *B   = new double[pn];					//���s��

	count=0;
	//���s��B[]�쐬 -��n0/(2d��t)hi(t)
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			ppn[count]=i;
			B[count]=-co*T[i];
			//B[count]=-co*PART[i].h/V;
			
			//B[count]-=PART[i].heat_generation*co2/(density[i]*Cp[i]);	//���M��[W/m3]
			B[count]-=0.5*(PART[i].heat_generation+PART[i].heat_gene_before1)*co2/(density[i]*Cp[i]);	//���M��[W/m3] 1�X�e�b�v�O�̔��M�ʂ�p���đ�`���ŋߎ�

			B[count]-=0.5*T_lap_before1[i]*co2;

			link[i]=count;
			count++;
		}
		else link[i]=pn+1;//�s��Ɋ܂܂�Ȃ����q�ɂ���а�Ƃ���(pn+1)���i�[
	}

	int number=0;			//�W���s��̔�[���v�f��
	for(int n=0;n<pn;n++)	//pn�͒f�M�����Ȃ�fluid_number,��f�M�Ȃ�particle_number�Ɉ�v
	{
		int i=ppn[n];		//n�Ԗڂ̖��m���̗��q�ԍ�
		int num=1;//���qi���ӂ̗��̗��q�� �����l��1�Ȃ͎̂������g���Ă��Ă��邩��
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			if(link[j]<pn) num++;
		}
		number+=num;
	}///number�����Ƃ܂���
	
    double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	int index=0;
	for(int n=0;n<pn;n++)
	{   
		
		ptr[n]=index;
		int KK=index;		//matrix�̑Ίp�������i�[�����ꏊ���L��
		double AA=0;
	    ind[index]=n;
	    index++;

		int i=ppn[n];

	    for(int k=0;k<PART[i].N2;k++)
	    {
			
	        int j=PART[i].NEI2[k]; 
			
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);

			if(link[j]<pn)//���qj�����m���Ȃ�
			{
				val[index]=w*K[i]*K[j]/(K[i]+K[j])*2*0.5;
				AA+=val[index];
				ind[index]=link[j];
				index++;
			}
			else
			{//�m�C�}���^�Ȃ�A�����͋L�q�Ȃ��H
				//B[n]-=w*PART[j].h/V*K[i];//���qj�����̂łȂ��Ȃ�
				//AA+=w*K[i];
			}
	    }
		val[KK]=-AA-co;
	}
	ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
	////////////////////*/

	//�����ō����s��͑Ίp���������ׂĕ��Ȃ̂ŁA����𐳂ɂ��邽�߂ɁA�W���s��Ɖ��s���-1��������
	for(int n=0;n<pn;n++)
	{
		for(int j=ptr[n];j<ptr[n+1];j++) val[j]*=-1;
		B[n]*=-1;
	}

	double *r=new double[pn];
	double *X=new double[pn];		//�s��̓����i�[
	double *AP = new double [pn];
	double *P = new double [pn];

	/////////////////////////�����l//////////////////
	if(CON->get_initial_u()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initial_u()==ON) //�����l��^����ꍇ�@����ȂɌ��ʂȂ�����
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//���Ԗڂ̖��m���ɊY�����闱�q�ԍ�
			X[n]=T[i];
		}
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	if(CON->get_T_CG()==0)//CG�@
	{
		//cout<<CON->get_T_CGep()<<endl;
		CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_T_CGep()); //CG�@�ɂ��s�������

		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//���Ԗڂ̖��m���ɊY�����闱�q�ԍ�
			PART[i].h+=(X[n]-T[i])*V;//�G���^���s�[�X�V
	//		laplacian[D][i]=(X[n]-PART[i].u[D])/(dt*vis0);//����������renewal�֐��Ȃ��ő��x�X�V����`�ɂȂ�
		}
		cout<<"������:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	}
	//////////*/

	//cout<<mass[191480]<<" "<<Cp[191480]<<endl;

	if(CON->get_T_CG()==1)//ICCG�@
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

		iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_T_CGep(),&count);
		
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//���Ԗڂ̖��m���ɊY�����闱�q�ԍ�
			//cout<<i<<endl;
			//PART[i].h+=(X[n]-T[i])*V;//�G���^���s�[�X�V
			//PART[i].h=X[n]*V;//�G���^���s�[�X�V
			PART[i].h=X[n]*mass[i]*Cp[i];//�G���^���s�[�X�V
		}
		cout<<"������:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	}
	///////////////////////////*/

	//���M�ʂ�������
	for(int i=0;i<particle_number;i++)
	{
		PART[i].heat_gene_before1=PART[i].heat_generation;	//1step�O�̏��ւƊi�[���Ȃ���
		PART[i].heat_generation=0;							//�������@���M�ʂ𑝂₷�Ƃ���+=�̋L�q�ɂ���ƁA�����Ȋ֐����Ŕ��M���N�������Ƃ��ɂ��Ώ��ł���
	}


	///�G���^���s�[����e���q�̉��xT[i]�����߂�(���x���v���b�g���邽��)
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+latent_H[i]*mass[i];			//�Z���I���_�̃G���^���s�[

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//�ő�
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];	//�Z�_
		else if(hs1<PART[i].h)											//�t��
		{       
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if(CON->get_material()==H2O && T[i]>=393) T[i]=393;//���ݕ����͍l���Ȃ��̂ŁA�������x��100�x�z�������100�x�ɖ߂�(�C���M�̎����̕����挈��?�j
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
    //////////////*/

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

	double Tmax=0;
	for(int i=0;i<particle_number;i++) if(T[i]>Tmax) Tmax=T[i];
	cout<<"�ő剷�x="<<Tmax<<"[K]"<<endl;

	delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;

	delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

	delete [] T;
    delete [] K;

	delete [] density;
	delete [] mass;
	delete [] Cp;
	delete [] MP;
	delete [] latent_H;

	delete [] ppn;
	delete [] link;

	delete [] T_lap_before1;
	
}

///���x���E�����ݒ�֐� //���E�̉��x���̂��̂������ȂǂŊ��ɂ��Ƃ܂��Ă���ꍇ�A�����Ή����鋫�E�Ɉʒu���闱�q�ɗ^����
void set_temperature_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double *T, double dt,int t)
{
	//cout<<"���x���E�����̓ǂݍ���"<<endl;
	////�����œǂ߂�悤�ɂ��邱��
	int b_num=16;//����_�̐�
	int t_num=469;//���ԏ����̐�
	double time=dt*t;//���݂̎���
	////
	double **bound_T=new double* [t_num];//���x�̋��E����
	for(int i=0;i<t_num;i++) bound_T[i]=new double [b_num];
	
	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]=0.0; 
	//�f�[�^�̓ǂݍ���
    

	//�t�@�C������f�[�^��ǂݍ��ށB
	int x=0;


	ifstream fin("Tbound.prn");//�G�N�Z���̌��f�[�^����A�X�y�[�X��؂�ŏo�͂����t�@�C��
	if(!fin) cout<<"cannot open Tbound"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);
	for(int i=0;i<t_num;i++)
	{ 
		for(int j=0;j<b_num;j++)
		{ 
			fin>>bound_T[i][j];
		}
	}
	fin.close();	

	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]+=273; //���f�[�^�́��Ȃ̂ŁAK�ɏC��
	//cout<<"x="<<x<<" boundT "<<bound_T[t_num-1][b_num-1]<<endl;

	if(CON->get_model_number()==26)//IH�����f��
	{

		
		if(CON->get_dimention()==2)//2����
		{

			for(int i=0; i<particle_number; i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL) //�ǋ��E���̓f�[�^�𗘗p���� ���W�n�ŕ�Ԃ�����A���Ԍn�ŕ�Ԃ���
				{
					double X= abs(PART[i].r[A_X]);
					double Y= PART[i].r[A_Y]+0.04; //model�z�u���͗��̉����Ə㕔��y���W��Βl����v���Ă��邽�߁A�킩��₷���悤�ɉ����̍��W��0�ɂȂ�悤�ɕ␳����
					double temp_t1;
					double temp_t2;
					int t1=(int) floor(time);
					int t2=(int) ceil(time);

					int x1=(int) floor(abs(X*100));//����_��10mm����
					int x2=(int) ceil(abs(X*100));
					int y1=(int) floor(Y*100);
					int y2=(int) ceil(Y*100);

					if(time<t_num-1)
					{
						if(Y<=0.0)//����
						{
							if(abs(X)<=0.08) //���ǂ̂����A����l�����݂���͈�
							{
								if(t1!=t2)
								{
								
									if(x1!=x2)
									{
										temp_t1= bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
										temp_t2= bound_T[t2][x1]+(bound_T[t2][x2]-bound_T[t2][x1])*(X-0.01*x1)/(0.01*(x2-x1));
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
									else
									{
										temp_t1= bound_T[t1][x1];
										temp_t2= bound_T[t2][x1];
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
								}
								else
								{
									if(x1!=x2)
									{
										T[i] = bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
									}
									else
									{
										T[i]= bound_T[t1][x1];
										
									}
								}
							}
							else //���ǁA���E������outwall�����B�Q�Ƃ���l�������̂łƂ肠����X=0.08�̂Ƃ��̒l�𗘗p
							{
								if(t1!=t2)
								{
									temp_t1= bound_T[t1][8];
									temp_t2= bound_T[t2][8];
									T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
								}
								else
								{
									T[i]= bound_T[t1][8];
								}
							}
						
						}
						else//���E�̕�
						{
							if(Y<0.07) //���ǂ̂����A����l�����݂���͈�
							{
								//y1+=8; //0����8�Ԗڂ̑���l�͒��(8���p)�Ȃ̂ŁA8����15�܂ł�y���W�ɉ����Q�Ƃ���悤�ɏC������
								//y2+=8;
								if(t1!=t2)
								{
									if(y1!=y2)
									{
										temp_t1= bound_T[t1][y1+8]+(bound_T[t1][y2+8]-bound_T[t1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
										temp_t2= bound_T[t2][y1+8]+(bound_T[t2][y2+8]-bound_T[t2][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1));
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
									else
									{
										temp_t1= bound_T[t1][y1+8];
										temp_t2= bound_T[t2][y1+8];
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
								}
								else
								{
									if(y1!=y2)
									{
										T[i]= bound_T[t1][y1+8]+(bound_T[t1][y2+8]-bound_T[t1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
									}
									else
									{
										T[i] = bound_T[t1][y1+8];
									}
								}

							}
							else //���E�ǁA������̑���l�����݂��Ȃ������B�Ƃ肠����Y=0.07�̂Ƃ��̒l�𗘗p
							{
								if(t1!=t2)
								{
									temp_t1= bound_T[t1][15];
									temp_t2= bound_T[t2][15];
									T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
								}
								else
								{
									T[i]= bound_T[t1][15];
								}
							}
						
						}
					}
					else//���ԏ���������l���z���Ă���ꍇ�A�ŏI�����̂��̂𗘗p����
					{
						if(Y<=0.0)//����
						{
							if(abs(X)<=0.08) //���ǂ̂����A����l�����݂���͈�
							{
								if(x1!=x2)
								{
									T[i]= bound_T[t_num-1][x1]+(bound_T[t_num-1][x2]-bound_T[t_num-1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									T[i]= bound_T[t_num-1][x1];
								}
							}
							else //���ǁA���E������outwall�����B�Q�Ƃ���l�������̂łƂ肠����X=0.08�̂Ƃ��̒l�𗘗p
							{
								T[i]= bound_T[t_num-1][8];
							}
						}
						else//���E�̕�
						{
							if(Y<0.07) //���ǂ̂����A����l�����݂���͈�
							{
								//y1+=8; //0����8�Ԗڂ̑���l�͒��(8���p)�Ȃ̂ŁA8����15�܂ł�y���W�ɉ����Q�Ƃ���悤�ɏC������
								//y2+=8;
								if(y1!=y2)
								{
									T[i]= bound_T[t_num-1][y1+8]+(bound_T[t_num-1][y2+8]-bound_T[t_num-1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									T[i]= bound_T[t_num-1][y1+8];
								}
							}
							else //���E�ǁA������̑���l�����݂��Ȃ������B�Ƃ肠����Y=0.07�̂Ƃ��̒l�𗘗p
							{
								T[i]= bound_T[t_num-1][15];
							}
						}
					}
				}
			}
		
		}

		if(CON->get_dimention()==3)//3����
		{
		}



	}

	for(int i=0;i<t_num;i++) delete [] bound_T[i];
    delete [] bound_T;
}

///���M���E�����ݒ�֐� //���E�̉��x�Ȃǂ������ȂǂŊ��ɂ��Ƃ܂��Ă���ꍇ�A�����Ή����鋫�E�Ɉʒu���闱�q�ɗ^����
void set_Q_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double dt,int t)
{
	//cout<<"���M�����̓ǂݍ���"<<endl;
	////�����œǂ߂�悤�ɂ��邱��
	int b_num=101;//����_�̐�
	int t_num=1;//���ԏ����̐�
	double time=dt*t;//���݂̎���
	////
	double **bound_T=new double* [t_num];//���x�̋��E����
	for(int i=0;i<t_num;i++) bound_T[i]=new double [b_num];
	
	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]=0.0; 
	//�f�[�^�̓ǂݍ���
    

	//�t�@�C������f�[�^��ǂݍ��ށB
	int x=0;


	ifstream fin("q0_f.txt");//�G�N�Z���̌��f�[�^����A�X�y�[�X��؂�ŏo�͂����t�@�C��
	if(!fin) cout<<"cannot open q0_f"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);
	for(int i=0;i<t_num;i++)
	{ 
		for(int j=0;j<b_num;j++)
		{ 
			//fin>>x;
			fin>>bound_T[i][j];
		}
	}
	fin.close();	

	//cout<<"x="<<x<<" boundT "<<bound_T[t_num-1][b_num-1]<<endl;

	if(CON->get_model_number()==26)//IH�����f��
	{
		if(CON->get_dimention()==2)//2����
		{

			for(int i=0; i<particle_number; i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL) //�ǋ��E���̓f�[�^�𗘗p���� ���W�n�ŕ�Ԃ�����A���Ԍn�ŕ�Ԃ���
				{
					double X= abs(PART[i].r[A_X]);
					double Y= PART[i].r[A_Y]; //model�z�u���͗��̉����Ə㕔��y���W��Βl����v���Ă��邽�߁A�킩��₷���悤�ɉ����̍��W��0�ɂȂ�悤�ɕ␳����
					double temp_t1;
					double temp_t2;
					//int t1=(int) floor(time);
					//int t2=(int) ceil(time);
					int t1=0;
					int t2=0;

					int x1=(int) floor(abs(X/0.0008));//����_��0.8mm����
					int x2=(int) ceil (abs(X/0.0008));
					int y1=(int) floor(Y/0.0008);
					int y2=(int) ceil(Y/0.0008);

					//if(t1==t2) cout<<"?"<<endl;

					//if(time<t_num-1)
					{
						//if(Y<=-0.05-0.0001)//���� INWALL��OUTWALL�S��
						if(Y<=-0.058-0.0001)//���� OUTWALL�̍ŊO�������ɂȂ�
						{
							if(abs(X)<=0.08) //���ǂ̂����A����l�����݂���͈�
							{
								if(x1!=x2)
								{
									PART[i].heat_generation += bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.0008*x1)/(0.0008*(x2-x1)); //��ʓI�Ȓ����̕������By-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									PART[i].heat_generation += bound_T[t1][x1];
								}
							}
							else //���ǁA���E������outwall�����B�Q�Ƃ���l�������̂łƂ肠����X=0.08�̂Ƃ��̒l�𗘗p
							{
								PART[i].heat_generation += bound_T[t1][100];
							}
						
						}
					}
					
				}
				//���ۂ̔��M�����̌��݂��l�����A���M�ʖ��x��ύX����
				PART[i].heat_generation/=CON->get_particle_volume()/(0.00055*CON->get_distancebp());
			}
		
		}

		if(CON->get_dimention()==3)//3����
		{
		}



	}

	for(int i=0;i<t_num;i++) delete [] bound_T[i];
    delete [] bound_T;
}