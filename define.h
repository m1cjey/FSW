#define DIMENTION 3
#define A_X 0
#define A_Y 1
#define A_Z 2
#define ON 1
#define OFF 0
#define PI 3.14159265358979323846
#define POSITIVE 0
#define NEGATIVE 1

#define FLUID   0
#define WALL 1
#define FRFLUID 0
#define BOFLUID 1
#define INWALL  2
#define OUTWALL 3
#define BDWALL 4
#define SOLID 5 //�ő�
#define BOSOLID 6//�ő̕\��
#define AIR 7   //��C
#define AIRWALL 8 //�Œ肳�ꂽ��C
#define NOCALC 9//�v�Z���Ȃ����q
#define CFD 10
#define ELASTIC 11//�e����
#define BOELASTIC 12//�e���̕\��
#define INELASTIC 13//�e���̊E��
#define ELASWALL 14//�e���̕�
#define BOELASWALL 15//�e���̕Ǖ\��
#define INELASWALL 16//�e���̕ǊE��
#define MOVE 2			//�ړ����q

#define ERR 1.0e-14//ERR 1.0e-14
#define WATER 17
#define ELECTRODE 18
#define MAGNET 19
#define COIL 20
#define CRUCIBLE 21
#define CONDUCT  30

#define ELECTRODE1 22	//�d��1(nanoe�~���d��)
#define ELECTRODE2 23	//�d��2(nanoe���d��)

//GPU CG�@�ɂ�����W���s��̊i�[���@
#define CSR_scl 0
#define CSR_vec 1
#define ELL 2

#define Diric 0
#define Neumn 1
#define BOTH  2
#define CONSTANT 0
#define LINER 1

#define UNDEFINED -1	//����`

#define Fe 0
#define Al 1
#define H2O 2
#define FSWA1 3