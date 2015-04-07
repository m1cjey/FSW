#ifndef FEM3DCLASS
#define FEM3DCLASS

#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B

class point3D//�ߓ_�N���X
{
public:
	double r[DIMENTION];
	int material; //�ގ�
	int boundary_condition; //���E���� 0=���ʁ@1,2=�Œ苫�E
	int particleID;			//�Ή����闱�q�ԍ� ���݂��Ȃ��Ƃ���-1���i�[
	int remesh;				//�����b�V���̈�ɑ����邩�A���Ȃ���
	int BD_node;			//���E�ߓ_���ǂ���
	
	double F;				//�͂̑傫��
	double B;	
	double Je;
	double value1;
};

class element3D//�v�f�N���X
{
public:
    int node[5];//�v�f���\�z����node�̔ԍ�  �����v�܂�聨�Ă��؂�̏��Ɋi�[
	int elm[5];//�v�f�Ɨאڂ���v�f�ԍ� //�\�ʂ�0�Ɗi�[
	int edge[6+1];//�v�f���\������Ӕԍ�
	double volume;//�̐ς̂U�{
	double r[DIMENTION];//�{���m�C�_���W(�O�ڋ����S)
	double RR;//�O�ڋ����a�̓��
	int map;//�}�b�s���O
	int material;//�ގ�
	int remesh;
};

class edge3D//�ӃN���X
{
public:
	int node[2+1];//�ӂ��\������ߓ_�ԍ��i�[
	int boundary_condition; //���E���� 0=���ʁ@1,2=�Œ苫�E
	int stat;//�ÓI�v�f�Ɉʒu����ӂ��ǂ��� 
	int static_num;//�ÓI�v�f���ɂ�����A�Ή�����1�X�e�b�v�ڂ̕Ӕԍ����i�[
	int para_node_num;//2���ߓ_�v�f���l����ꍇ�A�ӂ̒��_�ɐV���Ȑߓ_���ǉ������B�ǉ����ꂽ�ߓ_�͂��Ƃ��Ƃ������ߓ_�̌�̔ԍ��ɒǉ������B���̔ԍ����i�[
};

#endif
