#ifndef MPSPARTICLE
#define MPSPARTICLE

#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B

class mpsparticle
{
public:
	int id;
	int firstID;
	double r[DIMENTION];
	double r0[DIMENTION];
	double u[DIMENTION];
	double n[DIMENTION];	//�P�ʖ@���x�N�g��
	double old_A[DIMENTION];
	double F[DIMENTION]; //���q�ɓ����d����
	double potential[DIMENTION];//�\�ʒ��͍�//m/s2
	


	double P;//����
	double h;//�G���^���s�[
	double PND;//re��p�������q�����x
	double PND2;//re2��p�������q�����x 
	double val;//���̓s�x�K���Ȓl�̊i�[�ɗ��p
	double T; //���x[K]
	double L; //����� �ω𑜓x���q�@�֌W�ŗp����B���𑜓x�̏ꍇ��le�Ɠ��l
	double vis;//�S��
	double ep;//�����Ђ��ݗ�
	double sigma;//������������
	
	int type;	//FLUID INWALL OUTWALL
	int materialID;	//�ގ��ԍ� 1�Ȃ�density,2�Ȃ�density2�Ȃǂ��g�p
	int surface;//0:���� 1:�\��
	int index;//�i�[����Ă���i�q�̔ԍ�
	int N;    //re���ɑ��݂�����ӗ��q��
	int N2;   //re2���ɑ��݂�����ӗ��q��
	int N3;   //re3���ɑ��݂�����ӗ��q��
	int NEI[300];//re���ɑ��݂�����ӗ��q�ԍ�
	int NEI2[300];
	int NEI3[450]; 

	int toBEM;	//FEM�ɓ]�����邩�ǂ����BON or OFF
	int color;	//�F�@1�Ȃ��2�Ȃ�@FSW�ł̍P�v�I�ȐF�t���p

	double heat_gene_before1;
	double heat_generation;

	double dir_Pst;
	double dir_Pem;
	int trace;

	int division_flag;
};


#endif