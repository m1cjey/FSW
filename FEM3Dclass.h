#ifndef FEM3DCLASS
#define FEM3DCLASS

#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。

class point3D//節点クラス
{
public:
	double r[DIMENTION];
	int material; //材質
	int boundary_condition; //境界条件 0=普通　1,2=固定境界
	int particleID;			//対応する粒子番号 存在しないときは-1を格納
	int remesh;				//リメッシュ領域に属するか、しないか
	int BD_node;			//境界節点かどうか
	
	double F;				//力の大きさ
	double B;	
	double Je;
	double value1;
};

class element3D//要素クラス
{
public:
    int node[5];//要素を構築するnodeの番号  反時計まわり→てっぺんの順に格納
	int elm[5];//要素と隣接する要素番号 //表面は0と格納
	int edge[6+1];//要素を構成する辺番号
	double volume;//体積の６倍
	double r[DIMENTION];//ボロノイ点座標(外接球中心)
	double RR;//外接球半径の二乗
	int map;//マッピング
	int material;//材質
	int remesh;
};

class edge3D//辺クラス
{
public:
	int node[2+1];//辺を構成する節点番号格納
	int boundary_condition; //境界条件 0=普通　1,2=固定境界
	int stat;//静的要素に位置する辺かどうか 
	int static_num;//静的要素部における、対応する1ステップ目の辺番号を格納
	int para_node_num;//2次節点要素を考える場合、辺の中点に新たな節点が追加される。追加された節点はもともとあった節点の後の番号に追加される。その番号を格納
};

#endif
