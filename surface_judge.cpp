#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"	//#define �i�[
//include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
//#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
//#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
//#include<omp.h>
//#include<vector>
#include"function.h"

//�ߗח��q�T���֐�
void calc_neighbor_relation(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int fluid_number,int out)
{
	//�ߗח��q�T���̂��߂̏���
	int *INDEX=new int[CON->get_number_of_mesh()];	//�e�i�q�Ɋ܂܂�闱�q�����i�[
	reload_INDEX(CON,PART,particle_number,INDEX);//�i�q���̗��q���X�V

	int **MESH = new int *[CON->get_number_of_mesh()];
	for(int i=0;i<CON->get_number_of_mesh();i++) MESH[i]=new int [INDEX[i]];
		
	reload_INDEX2(CON,PART,particle_number,MESH);
	////////////////////////*/

		
	unsigned int timeA=GetTickCount();
	if(CON->get_freeon()==1) freeon(CON,PART,particle_number,n0_4,INDEX,MESH,fluid_number,out);//�\�ʔ���
	else if(CON->get_freeon()==2) freeon2(CON,PART,particle_number,n0_4,INDEX,MESH,fluid_number,out);//�\�ʔ���
	//else if(CON->get_freeon()==4) freeon4(CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,n0,fluid_number,e0,out,t);
	else cout<<"�\�ʔ��薢����"<<endl;

	if(CON->get_surface_judge2()==ON && CON->get_adaptive_sw()==OFF) 
	{
		//surface_judge2(CON,PART,fluid_number,particle_number);
		surface_judge2_old(CON,PART,fluid_number,particle_number);
		//surface_judge2_new(CON,PART,fluid_number,particle_number);
	}
	
	delete [] INDEX;
	for(int i=0;i<CON->get_number_of_mesh();i++) delete [] MESH[i];
	delete [] MESH;
	cout<<"���q�ˑ��֌W�v�Z�I��--time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


///�\�ʔ���֐�
void freeon(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out)
{
    ///����:freeon�֐��ł͊e���q�̗��q�����x�����߂�ƂƂ��ɁA����ꂽ���x����\�ʔ�����s���B�܂��A�Œᗱ�q�ԋ��������łɂ��Ƃ߂Ă���
	///freeon2���x�����A���̂Ԃ���񉻂��e�ՁB
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	
	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{//��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
		out=particle_number;
	}
	
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	//omp_set_num_threads(8);//�گ�ސ��w��
	#pragma omp parallel for
	for(int i=0;i<out;i++)//OUTWALL�ȊO�̗��q�B//OUTWALL�̗��q�����x�Ȃǂ͂���Ȃ�
	{    
		//printf("%d %d\n",i,omp_get_thread_num());//�ei�̌v�Z��S�����Ă���گ�ޔԍ��o��
	       
		PART[i].PND=0;//������
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;
		////���q�����x����
		double pnd=0;//���q�����x
		double pnd2=0;//���׼�ݗp���q�����x
		double pnd4=0;//�\�ʔ���p
		int N=0;
		int N2=0;
		int N3=0;
		if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
		{       
			for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
			{       
				for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];
				     
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							if(dis<=CON->get_re()*le && j!=i)//���z�E���U
							{       
								double r=CON->get_re()*le;
								double w=kernel(r,dis);
								pnd+=w;
								PART[i].NEI[N]=j;
								N++;
							}
							if(dis<=CON->get_re2()*le && j!=i)
							{       
								double r=CON->get_re2()*le;
								double w=kernel(r,dis);
								pnd2+=w;
								PART[i].NEI2[N2]=j;
								N2++;
							}
							if(dis<=CON->get_re3()*le && j!=i)//�\�ʒ���re3
							{       
								PART[i].NEI3[N3]=j;
								N3++;
							}
							if(dis<=CON->get_re4()*le && j!=i)
							{       
							    double r=CON->get_re4()*le;
								pnd4+=kernel(r,dis);
							}
						}
					}
				}
			}
		}
		PART[i].PND=pnd;
		PART[i].PND2=pnd2;
		PART[i].N=N;
		PART[i].N2=N2;
		PART[i].N3=N3;
		if(PART[i].N3>800) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==FLUID)
		{
			if(pnd4<n0_4*CON->get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�\�ʗ��q�Ƃ���
				//PART[i].P=0;
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==INWALL)
		{
			if(pnd4<n0_4*CON->get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�Ǖ\�ʗ��q�Ƃ���
				//PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	    
	}

	/*///�Œᗱ�q�ԋ��������Ƃ߂�
	double min0=CON->get_distancebp();//�Œᗱ�q�ԋ���
	int type1,type2,surface1,surface2;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////�Œᗱ�q�ԋ��������Ƃ܂���
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
	*mindis=min0;
	///*/
}

///�\�ʔ���֐�ver.2
void freeon2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out)
{
	///freeon2�֐��̐���:flag1[i]�̓����ɂ�荂�����B�������A1CPU�Ȃ瑁�����A���CPU�ɂ����񉻂͌�����
	//cout<<"�\�ʔ���(freeon2)"<<endl;
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	double d=2;
	if(CON->get_dimention()==3) d=3;

	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{
		out=particle_number;	//��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	}

	/////�ǋ��E�|���S�����̃e�X�g
	int wall_poly_number=0;//�ǋ��E�̃|���S����
	wall_poly_number=3;//�ǂݍ��݂Ȃǂł��Ƃ߂邱��
	double wall_a[3];//2�����̒���������ax+by+c=0�@3DFEM���Q�l��node��edge�̂悤�ȏ����������邱��
	double wall_b[3];
	double wall_c[3];
	if(CON->get_wall_poly()==1)
	{
		//�ǋ��E�̐ݒ� //���b�V���Ȃǂ���ǂݍ��ނ��Ƃ��l����ƁA���̈ʒu�͂܂�������
		if(CON->get_model_number()==1)
		{
			if(CON->get_dimention()==2)
			{
				wall_a[0]=1; wall_b[0]=0; wall_c[0]=0.011;//��������
				wall_a[1]=1; wall_b[1]=0; wall_c[1]=-0.011;//�E������
				wall_a[2]=0; wall_b[2]=1; wall_c[2]=0.021;//��������
			}
		}
	}
	////�Ǐd�݊֐��̌v�Z
	double wallZ_pnd=0;//1E-05x-1.636
	double wallZ_pnd2=0;
	double wallZ_pnd4=0;//4E-06x-1.979
	///////////////
	double *PND4=new double [out];//�\�ʔ���p���q�����x
	int *flag1=new int [out];		//�����t���O�B0:�������@1:�����ς�
	///������
	#pragma omp parallel for
	for(int i=0;i<out;i++)
	{
		PART[i].PND=0;//������
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;

		PND4[i]=0;
		flag1[i]=0;
	}
	
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	for(int i=0;i<out;i++)//OUTWALL�ȊO��CFD���q�B//OUTWALL�̗��q�����x�Ȃǂ͂���Ȃ�
	{    
		if(CON->get_wall_poly()==0)
		{
			////���q�����x����
			if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
			{       
				for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
				{       
					for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
					{
						for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
						{
							for(int L=0;L<INDEX[I+J+K];L++)
							{       
								int j=MESH[I+J+K][L];
								if(j<out)
								{
									if(flag1[j]==0 && j!=i)//�܂��������ĂȂ��Ȃ�
									{
										double X=PART[j].r[A_X]-PART[i].r[A_X];
										double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
										double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
										double dis=sqrt(X*X+Y*Y+Z*Z);
							
										if(dis<=CON->get_re()*le)//���z�E���U
										{       
											double r=CON->get_re()*le;
											double w=kernel(r,dis);
											PART[i].PND+=w;
											PART[j].PND+=w;
											PART[i].NEI[PART[i].N]=j;
											PART[j].NEI[PART[j].N]=i;
											PART[i].N++;
											PART[j].N++;
										}
										if(dis<=CON->get_re2()*le)
										{       
											double r=CON->get_re2()*le;
											double w=kernel(r,dis);
							
											PART[i].PND2+=w;
											PART[j].PND2+=w;
										
											PART[i].NEI2[PART[i].N2]=j;
											PART[j].NEI2[PART[j].N2]=i;
											PART[i].N2++;
											PART[j].N2++;
										}
										if(dis<=CON->get_re3()*le)//�\�ʒ���re3
										{       
											PART[i].NEI3[PART[i].N3]=j;
											PART[j].NEI3[PART[j].N3]=i;
											PART[i].N3++;
											PART[j].N3++;
										}
										if(dis<=CON->get_re4()*le)
										{       
											double r=CON->get_re4()*le;
											double w=kernel(r,dis);
											PND4[i]+=w;
											PND4[j]+=w;
										}
									}
								}
								else
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON->get_re()*le)//���z�E���U
									{       
										double r=CON->get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[i].N++;
									}
									if(dis<=CON->get_re2()*le)
									{       
										double r=CON->get_re2()*le;
										double w=kernel(r,dis);
										PART[i].PND2+=w;
										PART[i].NEI2[PART[i].N2]=j;
										PART[i].N2++;
									}
									if(dis<=CON->get_re3()*le)//�\�ʒ���re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[i].N3++;
									}
									if(dis<=CON->get_re4()*le)
									{       
										double r=CON->get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
									}
								}
							}
						}
					}
				}
			}
		}

		if(CON->get_wall_poly()==1)//�|���S���ǁB�����ł͗��̗��q���m�̑��ݍ�p�����v�Z�����
		{
			////���q�����x����
			if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
			{       
				for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
				{       
					for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
					{
						for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
						{
							for(int L=0;L<INDEX[I+J+K];L++)
							{       
								int j=MESH[I+J+K][L];
								if(j<out)
								{
									if(flag1[j]==0 && j!=i)//�܂��������ĂȂ��Ȃ�
									{
										double X=PART[j].r[A_X]-PART[i].r[A_X];
										double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
										double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
										double dis=sqrt(X*X+Y*Y+Z*Z);
							
										if(dis<=CON->get_re()*le)//���z�E���U
										{       
											double r=CON->get_re()*le;
											double w=kernel(r,dis);
											PART[i].PND+=w;
											PART[j].PND+=w;
											PART[i].NEI[PART[i].N]=j;
											PART[j].NEI[PART[j].N]=i;
											PART[i].N++;
											PART[j].N++;
										}
										if(dis<=CON->get_re2()*le)
										{       
											double r=CON->get_re2()*le;
											double w=kernel(r,dis);
							
											PART[i].PND2+=w;
											PART[j].PND2+=w;
										
											PART[i].NEI2[PART[i].N2]=j;
											PART[j].NEI2[PART[j].N2]=i;
											PART[i].N2++;
											PART[j].N2++;
										}
										if(dis<=CON->get_re3()*le)//�\�ʒ���re3
										{       
											PART[i].NEI3[PART[i].N3]=j;
											PART[j].NEI3[PART[j].N3]=i;
											PART[i].N3++;
											PART[j].N3++;
										}
										if(dis<=CON->get_re4()*le)
										{       
											double r=CON->get_re4()*le;
											double w=kernel(r,dis);
											PND4[i]+=w;
											PND4[j]+=w;
										}
									}
								}
								else
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON->get_re()*le)//���z�E���U
									{       
										double r=CON->get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[i].N++;
									}
									if(dis<=CON->get_re2()*le)
									{       
										double r=CON->get_re2()*le;
										double w=kernel(r,dis);
										PART[i].PND2+=w;
										PART[i].NEI2[PART[i].N2]=j;
										PART[i].N2++;
									}
									if(dis<=CON->get_re3()*le)//�\�ʒ���re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[i].N3++;
									}
									if(dis<=CON->get_re4()*le)
									{       
										double r=CON->get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
									}
								}
							}
						}
					}
				}
			}
		}
		flag1[i]=1;//�����I��

		if(CON->get_wall_poly()==1)//�|���S���ǁB�����ł͕ǂ���̍�p���v�Z�����
		{
			double riw=100*le;//�ǋ��E�Ɨ��q�̋���
			double dis_w=0;
			for(int p=0; p<3;p++)//�e�ǋ��E����̗��q�Ƃ̋��������߂�Bp�̒l�̓|���S�������Q�Ƃł���悤�ɂ��邱��
			{
				dis_w=abs(wall_a[p]*PART[i].r[A_X]+wall_b[p]*PART[i].r[A_Y]+wall_c[p])/(sqrt(wall_a[p]*wall_a[p]+wall_b[p]*wall_b[p]));
				riw=dis_w;
					/*/
				if(dis_w<=riw)
				{
					riw=dis_w;

				}
				/*/
				if(riw<=CON->get_re()*le) PART[i].PND+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
				if(riw<=CON->get_re2()*le) PART[i].PND2+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
				if(riw<=CON->get_re4()*le)//pow(riw,-1.979)*pow(10.0,-6)*4;//4E-06x-1.979//= -9.245ln(x) - 55.708
				{
					if(-9.245*log(riw)-55.708>=0) PND4[i]+=-9.245*log(riw)-55.708;
				}
			}
			/*////
			if(riw<=CON->get_re()*le) PART[i].PND+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
			if(riw<=CON->get_re2()*le) PART[i].PND2+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
			if(riw<=CON->get_re4()*le)//pow(riw,-1.979)*pow(10.0,-6)*4;//4E-06x-1.979//= -9.245ln(x) - 55.708
			{
				if(-9.245*log(riw)-55.708>=0) PND4[i]+=-9.245*log(riw)-55.708;
			}
			///*/

		}

		if(PART[i].N3>450) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==FLUID)
		{
			if(PND4[i]<n0_4*CON->get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�\�ʗ��q�Ƃ���
				PART[i].P=0;
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON->get_beta())//���ȉ��Ȃ�
			{
				PART[i].surface=ON;//�Ǖ\�ʗ��q�Ƃ���
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	}

	/*///�Œᗱ�q�ԋ��������Ƃ߂�
	double min0=CON->get_distancebp();//�Œᗱ�q�ԋ���
	int type1,type2,surface1,surface2;
	type1=0;
	type2=0;
	surface1=0;
	surface2=0;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////�Œᗱ�q�ԋ��������Ƃ܂���
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	//////////
	*mindis=min0;
	

	if(t==1)
	{
		ofstream fouts("mindis.dat");
		fouts.close();
	}

	ofstream fouts2("mindis.dat",ios :: app);
	fouts2<<min0<<endl;
	fouts2.close();
	/////*/

	/*/���q�����x�o��
	plot_PND(CON,PART,fluid_number,PND4,t);
	if(CON->get_PND_interval()>0)
	{
		if(t==1 || t%CON->get_PND_interval()==0) plot_PND_each(CON,PART,fluid_number,PND4,t);
	}
	///*/


	delete [] PND4;
	delete [] flag1;
}

//�\�ʔ���֐�ver.2
void surface_judge2_old(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//�]���̕\�ʔ���ɉ����āA�����ł���ɂӂ邢�ɂ�����B
	//���qi�̖@���x�N�g���ƁA���qj�����޸�قƂ̊p�x������l�𒴂��Ă�����\�ʂł͂Ȃ��Ɣ��f
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	for(int i=0;i<particle_number;i++)
    {
       if(PART[i].type==FLUID ||PART[i].type==INWALL) PART[i].surface==ON;//���q�����x�ɂ��\�ʔ����p���Ȃ�
	}
	
    //////�@���޸�ٌv�Z
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
					if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
				}
			}
			PART[i].surface=flag;
		}
	}

	//��
	//int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM�p�̃��b�V�����쐬����ۂɁA���ׂĂ̕Ǘ��q���o�͂���ꍇ���AOUTWALL�̕\�ʔ��������K�v�����邽�߁Aout=particle_number�Ƃ���B

	//if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
			{
				int flag=ON;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j>=fluid_number)							//�Ǘ��q��surface_jusge2�͕Ǘ��q�����ōs��
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
						double ny=Y/dis;
						double nz=Z/dis;
						double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
						if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
					}
				}
				PART[i].surface=flag;
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
}

//�\�ʔ���֐�ver.2
void surface_judge2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//�]���̕\�ʔ���ɉ����āA�����ł���ɂӂ邢�ɂ�����B
	//���qi�̖@���x�N�g���ƁA���qj�����޸�قƂ̊p�x������l�𒴂��Ă�����\�ʂł͂Ȃ��Ɣ��f
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	
    //////�@���޸�ٌv�Z
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		//if(PART[i].surface==ON)  direct_f2(CON,PART,i,direct);//���݂̂̂����݂���������ŋ��߂�
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}
	double le=CON->get_distancebp();
	//�@���x�N�g�����o��
	ofstream nn("normal.dat");
	if(CON->get_dimention()==2) for(int i=0;i<particle_number;i++) nn<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<direct[A_X][i]*le<<" "<<direct[A_Y][i]*le<<endl;
	nn.close();
	
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(j<fluid_number)
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(PART[j].type==FLUID || PART[j].type==INWALL)
				//if(PART[j].type==FLUID || (PART[j].type==INWALL && dis<1.5*CON->get_distancebp()))
				//if(PART[j].type==FLUID)// �e�X�g�@�Ǖt�߂̔��肪�}�V�ɂȂ邩�H
				{
					
					double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
					if(inp<-0.6) 
					{
						flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
						//cout<<"judge2�Ŕ���ϊ� i="<<i<<endl;
					}
				}
			}
			PART[i].surface=flag;
		}
	}
	
	//��
	int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM�p�̃��b�V�����쐬����ۂɁA���ׂĂ̕Ǘ��q���o�͂���ꍇ���AOUTWALL�̕\�ʔ��������K�v�����邽�߁Aout=particle_number�Ƃ���B

	////
	if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
			{
				if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
				{
					int flag=ON;
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(j>=fluid_number)							//�Ǘ��q��surface_jusge2�͕Ǘ��q�����ōs�� �Ȃ��H
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
							double ny=Y/dis;
							double nz=Z/dis;
							double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
							if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
						}
					}
					PART[i].surface=flag;
				}
			}
		}
	}
	////*/

	/*//////
	//�Ǘ����ĕ\�ʗ��q�Ȃ��̂́A�������q�Ɣ��肷��
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=OFF;					//ON�Ȃ�\�ʁ@OFF�Ȃ����
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==FLUID &&PART[j].surface==ON && dis<CON->get_re()*CON->get_distancebp()) flag=ON;	//���͂ɕ\�ʗ��q������΂���ł悵
			}
			PART[i].surface=flag;
			//if(flag==OFF) cout<<"judge2�ŌǗ��\�ʔ���ϊ� i="<<i<<endl; //���������Ǝ~�܂�H

		}
	}//////////////*/

	/*/�Ǘ����ĕ\�ʗ��q��INWALL�́A�������q�Ɣ��肷��
	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			int flag=OFF;					//ON�Ȃ�\�ʁ@OFF�Ȃ����
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==INWALL &&PART[j].surface==ON) flag=ON;	//���͂ɕ\�ʗ��q������΂���ł悵
			}
			PART[i].surface=flag;
			
		}
	}//////////////*/
	
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];

}

////
//�\�ʔ���֐�
void surface_judge2_new(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{

	//���qi�̖@���x�N�g���ƁA���qj�����޸�قƂ̊p�x������l�𒴂��Ă�����\�ʂł͂Ȃ��Ɣ��f
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	for(int i=0;i<particle_number;i++)
    {
       if(PART[i].type==FLUID ||PART[i].type==INWALL) PART[i].surface==ON;//���q�����x�ɂ��\�ʔ����p���Ȃ�
	}

    //////�@���޸�ٌv�Z
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}
	
	int d=CON->get_dimention();

	double *L=new double [fluid_number];		//���̊֐��Ŏg�p�������q�ԋ���

	for(int i=0;i<fluid_number;i++)  L[i]=PART[i].L;

	/*////
	for(int i=0;i<fluid_number;i++)
	{
		int jnb=PART[i].N;
		if(jnb>4)
		{
			for(int n=0;n<4;n++)
			{
				int temp_k=n;
				int A=PART[i].NEI[n];
				int temp_j=A;
				double X=PART[i].r[A_X]-PART[A].r[A_X];
				double Y=PART[i].r[A_Y]-PART[A].r[A_Y];
				double Z=PART[i].r[A_Z]-PART[A].r[A_Z];
				double mindis=X*X+Y*Y+Z*Z;
					
				for(int k=n+1;k<jnb;k++) 
				{
					int j=PART[i].NEI[k];
					X=PART[j].r[A_X]-PART[i].r[A_X];
					Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=X*X+Y*Y+Z*Z;
					if(dis<mindis)
					{
						mindis=dis; 
						temp_k=k;
						temp_j=j;
					}
				}
				PART[i].NEI[temp_k]=A;
				PART[i].NEI[n]=temp_j;
			}
			double ave_dis=0;
			for(int n=0;n<4;n++)
			{
				int j=PART[i].NEI[n];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				ave_dis+=dis*0.25;
			}
			L[i]=ave_dis;
			//L[i]=PART[i].L;
			//cout<<i<<" "<<L[i]/PART[i].L<<endl;
		}
	}////////////*/

	
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=ON;					//ON�Ȃ�\�ʁ@OFF�Ȃ����
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(PART[j].type==FLUID || (PART[j].type==INWALL && dis<1.2*PART[i].L))		//�Ǘ��q�̏ꍇ�͐ڋ߂̒�`�����߂�
				{
					//if(dis<CON->get_re()*L[i])
					{
					double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP�ʑ��΋����޸�قƂ̓���
					if(inp<-0.6) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
					}
				}
			}
			PART[i].surface=flag; 
		}
	}


	//�Ǘ����ĕ\�ʗ��q�Ȃ��̂́A�������q�Ɣ��肷��
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=OFF;					//ON�Ȃ�\�ʁ@OFF�Ȃ����
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==FLUID &&PART[j].surface==ON && dis<CON->get_re()*L[i]) flag=ON;	//���͂ɕ\�ʗ��q������΂���ł悵
			}
			PART[i].surface=flag;
		}
	}//////////////*/
	
	

	//��
	//int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM�p�̃��b�V�����쐬����ۂɁA���ׂĂ̕Ǘ��q���o�͂���ꍇ���AOUTWALL�̕\�ʔ��������K�v�����邽�߁Aout=particle_number�Ƃ���B

	//if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
			{
				int flag=ON;
				double direct2[3];
				for(int D=0;D<d;D++) direct2[D]=direct[D][i]*(-1);	//�O�����@�� 
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j<fluid_number && PART[j].surface==OFF)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<1.5*CON->get_distancebp())
						{
							X/=dis; Y/=dis; Z/=dis;
							double COS=X*direct2[A_X]+Y*direct2[A_Y]+Z*direct2[A_Z];
							if(COS>sqrt(3.0)*0.5) flag=OFF;
							//flag=OFF;	//�������̗��q�̋߂��ɂ���ǂ͓���
						}
					}
				}
				PART[i].surface=flag;
			}
			//else if(PART[i].type==OUTWALL) PART[i].surface=ON;	//OUTWALL�͂Ƃ肠�����\�ʂɂ��Ă���
		}
	}

	//�Ǘ����ĕ\�ʗ��q��INWALL�́A�������q�Ɣ��肷��
	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			int flag=OFF;					//ON�Ȃ�\�ʁ@OFF�Ȃ����
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==INWALL &&PART[j].surface==ON) flag=ON;	//���͂ɕ\�ʗ��q������΂���ł悵
			}
			PART[i].surface=flag;
		}
	}//////////////*/

	//�������q�Ɣ��肳�ꂽ��A�@���x�N�g�����[���ɂ��A�����͂��[���ɂ��Ă���
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			//for(int D=0;D<d;D++) PART[i].direct[D]=0;		//adaptive2�̂��Ƃ��l����ƁA�����Ŗ@���������Ȃ��ق����悢
			PART[i].P=0;
		}
	}

	delete [] L;
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
}
///*/

//�\�ʔ���֐�ver.3
void surface_judge3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//�]���̕\�ʔ���ɉ����āA�����ł���ɂӂ邢�ɂ�����B
	//���qi�̖@���x�N�g���ƁA���qj�����޸�قƂ̊p�x������l�𒴂��Ă�����\�ʂł͂Ȃ��Ɣ��f
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	
    //////�@���޸�ٌv�Z
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
					if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
				}
			}
			PART[i].surface=flag;
		}
	}

	//��
	int flag=OFF;
	if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM�p�̃��b�V�����쐬����ۂɁA���ׂĂ̕Ǘ��q���o�͂���ꍇ���AOUTWALL�̕\�ʔ��������K�v�����邽�߁Aout=particle_number�Ƃ���B

	if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].surface==ON)//���̕\�ʗ��q���{���ɕ\�ʗ��q���肤�邩���肷��
			{
				int flag=ON;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j>=fluid_number)							//�Ǘ��q��surface_jusge2�͕Ǘ��q�����ōs��
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double nx=X/dis;//i���q����j���q�ւƌ������P���޸��
						double ny=Y/dis;
						double nz=Z/dis;
						double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//�@���޸�قƒP���޸�قƂ̓���
						if(inp<-0.5) flag=OFF;//���ς�-0.5�̂��̂��ЂƂł�����Γ������̗��q�B���ς�-0.5�Ƃ������Ƃ́A���̊p�x��120�x�����e�͈͂Ƃ�������
					}
				}
				PART[i].surface=flag;
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];

}

///freeon�֐� ver.3 ���q�����x�̂ݍČv�Z�B���q�\���q�֌W�͕ω����Ȃ��Ɖ��肵�Ă���
void freeon3(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int out)
{
	double d=2;
	if(CON->get_dimention()==3) d=3;
	double R=CON->get_re()*CON->get_distancebp();
	double R2=CON->get_re2()*CON->get_distancebp();
	for(int i=0;i<out;i++)
	{
		double W=0;
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);
			//double w=kernel2(R,dis,d);
			W+=w;
		}
		PART[i].PND=W;
	}
	if(CON->get_re()!=CON->get_re2())
	{
		for(int i=0;i<out;i++)
		{
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R2,dis);
				//double w=kernel2(R2,dis,d);
				W+=w;
			}
			PART[i].PND2=W;
		}
	}
	else for(int i=0;i<out;i++) PART[i].PND2=PART[i].PND;
}

//���q�����x�o��
void plot_PND(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t)
{
	double le=CON->get_distancebp();

	ofstream fp("PND.dat");
	if(CON->get_dimention()==2) for(int i=0;i<fluid_number;i++)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PND4[i]<<endl;
	if(CON->get_dimention()==3) for(int i=0;i<fluid_number;i++)	if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PND4[i]<<endl;
	fp.close();///////////*/
}

//���q�����x�v���b�g�֐�(���X�e�b�v�o��)
void plot_PND_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t)
{
	double le=CON->get_distancebp();

	char filename[30];
	sprintf_s(filename,"PND%d.dat", t);

	ofstream fp(filename);
	if(CON->get_dimention()==2) for(int i=0;i<fluid_number;i++)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PND4[i]<<endl;
	if(CON->get_dimention()==3) for(int i=0;i<fluid_number;i++)	if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PND4[i]<<endl;
	fp.close();///////////*/
}