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
#define SOLID 5 //固体
#define BOSOLID 6//固体表面
#define AIR 7   //空気
#define AIRWALL 8 //固定された空気
#define NOCALC 9//計算しない粒子
#define CFD 10
#define ELASTIC 11//弾性体
#define BOELASTIC 12//弾性体表面
#define INELASTIC 13//弾性体界面
#define ELASWALL 14//弾性体壁
#define BOELASWALL 15//弾性体壁表面
#define INELASWALL 16//弾性体壁界面
#define MOVE 2			//移動粒子

#define ERR 1.0e-14//ERR 1.0e-14
#define WATER 17
#define ELECTRODE 18
#define MAGNET 19
#define COIL 20
#define CRUCIBLE 21
#define CONDUCT  30

#define ELECTRODE1 22	//電極1(nanoe円柱電極)
#define ELECTRODE2 23	//電極2(nanoe平板電極)

//GPU CG法における係数行列の格納方法
#define CSR_scl 0
#define CSR_vec 1
#define ELL 2

#define Diric 0
#define Neumn 1
#define BOTH  2
#define CONSTANT 0
#define LINER 1

#define UNDEFINED -1	//未定義

#define Fe 0
#define Al 1
#define H2O 2
#define FSWA1 3