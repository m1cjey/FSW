#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
#include"FEM3Dclass.h"
#include<omp.h>
#include<vector>
#include"function.h"

void box3D(vector<point3D> &NODE,vector<element3D> &ELEM,double alpha,int *nelm,double rax,double ray,double raz,int KTJ);
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip);
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);
int iface3D(vector<element3D> &ELEM,int ielm,int jelm);
void qsorti3D(vector<point3D> &NODE,vector<element3D> &ELEM,int n,int *list);
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void remove3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,int iv,int *kv);
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm);
void set_material(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm);
int poly3D_for_FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig *CON,double rrm,int startID);
int make_air_layer(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,mpsconfig *CON,int *node_num,int KTE,double rrm,int KTJ);



////デローニ分割main関数
void delaun3D_main(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node_num,int *nelm,int FINE_sw)
{
	/////////////節点座標の正規化
    double xmin=NODE[1].r[A_X];
    double ymin=NODE[1].r[A_Y];
    double zmin=NODE[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    ///座標の最大、最小値を求める
    for(int i=2;i<=*node_num;i++)
    {
        if(NODE[i].r[A_X]<xmin) xmin=NODE[i].r[A_X];
		else if(NODE[i].r[A_X]>xmax) xmax=NODE[i].r[A_X];
	
		if(NODE[i].r[A_Y]<ymin) ymin=NODE[i].r[A_Y];
		else if(NODE[i].r[A_Y]>ymax) ymax=NODE[i].r[A_Y];
	
		if(NODE[i].r[A_Z]<zmin) zmin=NODE[i].r[A_Z];
		else if(NODE[i].r[A_Z]>zmax) zmax=NODE[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X軸方向の寸法
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y軸方向の寸法
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z軸方向の寸法
    double rmax=rax;		///最大寸法
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //ここはelseにしたらダメ

    ///座標変換
    double rrm=1.000000/rmax;///こういう書き方をすることで、数値誤差を減らせる・・？
    for(int i=1;i<=*node_num;i++)
    {   //   A/Bという計算をしたとき、Ａの値によって微妙に1/Bという倍率がちがってくるのではないかと考えて、下のような書き方にしている
        NODE[i].r[A_X]=(NODE[i].r[A_X]-xmin)*rrm;
		NODE[i].r[A_Y]=(NODE[i].r[A_Y]-ymin)*rrm;
		NODE[i].r[A_Z]=(NODE[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////
	

	delaun3D(CON,NODE,ELEM, KTJ, KTE, rax, ray, raz,node_num,nelm, FINE_sw,rrm);
	/////ﾒｯｼｭ生成完了

	////座標を元に戻す
    for(int i=1;i<=*node_num;i++)
    {
        NODE[i].r[A_X]=rmax*NODE[i].r[A_X]+xmin;
		NODE[i].r[A_Y]=rmax*NODE[i].r[A_Y]+ymin;
		NODE[i].r[A_Z]=rmax*NODE[i].r[A_Z]+zmin;
    }
}

///デローニ分割
void delaun3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,double rax,double ray,double raz,int *node_num,int *nelm,int FINE_sw,double rrm)
{   
	cout<<"ﾃﾞﾛｰﾆ分割開始----";

	unsigned int timeA=GetTickCount();	//計算開始時刻
	double alpha=CON->get_boxalpha();	//ｽｰﾊﾟｰﾎﾞｯｸｽの大きさを決める係数α
    int node=*node_num;
    
    for(int i=1;i<KTE;i++) ELEM[i].map=0;//初期化
    
    ////ｽｰﾊﾟｰﾎﾞｯｸｽへの移動量を計算(ｽｰﾊﾟｰﾎﾞｯｸｽの大きさはrax*alpha,ray*alpha,raz*alpha)
    double xbar=0.500000*(alpha-1.000000)*rax;
    double ybar=0.500000*(alpha-1.000000)*ray;
    double zbar=0.500000*(alpha-1.000000)*raz;
    
    ///ｽｰﾊﾟｰﾎﾞｯｸｽへ移動 たんなる平行移動であって、倍率は変化していない
    for(int i=1;i<=node;i++)
    {
        NODE[i].r[A_X]+=xbar;
		NODE[i].r[A_Y]+=ybar;
		NODE[i].r[A_Z]+=zbar;
    }
    ///////////
    
    ///ｽｰﾊﾟｰﾎﾞｯｸｽを6つの四面体に分割する
    box3D(NODE,ELEM,alpha,nelm,rax,ray,raz,KTJ);//要素が6つ生成され、各要素の体積・外接球中心座標がもとまった
    
    ///順次節点を導入していく
    int *kv=new int[KTE];//新節点を外接球にふくむ要素群
    int *istack=new int[KTE];//一時配列
	
	//int flag[node]=ON;
	//int num_ON=node;
	//while(num_ON!=0){
		for(int i=1;i<=node;i++)
		{   
			//if(flag[i]==ON{
			int ip=i;
			double xp=NODE[ip].r[A_X];//導入する節点の座標
			double yp=NODE[ip].r[A_Y];
			double zp=NODE[ip].r[A_Z];
		
			///新節点を含む要素の探索
			int loc=locate3D(NODE,ELEM,*nelm,xp,yp,zp);
		
			//////////外接球内に新節点を含む要素の抽出
			int iv=0;
			int msk=0;
	
			iv++;//外接球内に新節点を含む要素数
			kv[iv]=loc;
			ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むかどうかを検査済みということ
			msk++;
			istack[msk]=loc;
		
			while(msk!=0)
			{   
				int isk=istack[msk];//いま注目している要素の番号
				msk--;
				for(int j=1;j<=4;j++)
				{
					int jelm=ELEM[isk].elm[j];//iskと接する要素
					if(jelm!=0)//それが表面でないなら
					{
						if(ELEM[jelm].map!=1) //まだ検査してないなら
						{   
		         
							double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
							double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
						
							if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
							{
								dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
								if(dst<rad)
								{
									dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
									if(dst<rad)//外接球内に含む
									{
										iv++;//外接球内に新節点を含む要素数を＋１
										kv[iv]=jelm;//リストにいれる
										ELEM[jelm].map=1;//外接球に新点含むというしるし
										msk++;
										istack[msk]=jelm;
									}
								}
							}
						}
					}
				}
			}//外接球内に新節点を含む要素数ivと、その要素番号kv[iv]がもとまった
		/////////////////
		
			int n0=*nelm;
			////得られた多面体を四面体に分割する
			poly3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
		
		//if(succes) flag[i]=OFF

    }
    //count ON


    ///ｽｰﾊﾟｰﾎﾞｯｸｽを頂点にもつ四面体を除去する
    ///ここでivとkvの役目を変える。ivは除去される四面体数であり、kvはそのリスト。
    int iv=0;
    for(int i=1;i<=*nelm;i++)  kv[i]=0;//初期化
    
    for(int i=1;i<=*nelm;i++)
    {   //どれかひとつでもｽｰﾊﾟｰﾎﾞｯｸｽの頂点が含まれていたら
        if(ELEM[i].node[1]>KTJ || ELEM[i].node[2]>KTJ||ELEM[i].node[3]>KTJ||ELEM[i].node[4]>KTJ)
		{
			iv++;
			kv[iv]=i;
		}
    }
	int before_nelm=*nelm;//除去前要素数記憶

	
	cout<<"ok 除去前要素数="<<before_nelm<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

    remove3D(NODE,ELEM,nelm,iv,kv);
    ///////////////////
    
    ////要素がうまく生成されているかチェックする
    fill3D(NODE,ELEM,*nelm);

	/////要素材質決定
    set_material(CON,NODE,ELEM,node,*nelm);

	//set_depth(CON,NODE,ELEM,node,*nelm,depth,KTE);//要素深さ決定

    ///座標を元に戻す  ただし完全に戻すのはmainで
    for(int i=1;i<=node;i++)
    {
        NODE[i].r[A_X]-=xbar;
		NODE[i].r[A_Y]-=ybar;
		NODE[i].r[A_Z]-=zbar;
    }

    delete [] kv;
    delete [] istack;

	////////////////////////////////////////////ﾒｯｼｭの自動細分割
	if(FINE_sw==ON) FINE3D(NODE,ELEM,KTJ,KTE,node_num,nelm,CON,rrm,1);

}

///ｽｰﾊﾟｰﾎﾞｯｸｽ作成関数
void box3D(vector<point3D> &NODE,vector<element3D> &ELEM,double alpha,int *nelm,double rax,double ray,double raz,int KTJ)
{
	//ｽｰﾊﾟｰﾎﾞｯｸｽの大きさはrax*alpha,ray*alpha,raz*alpha

    *nelm=6;
    double xone=alpha*rax;//最大座標
    double yone=alpha*ray;
    double zone=alpha*raz;
    
    
    /////四面体情報作成
    //節点座標
    NODE[KTJ+1].r[A_X]=0.000000;
    NODE[KTJ+1].r[A_Y]=0.000000;
    NODE[KTJ+1].r[A_Z]=0.000000;
    
    NODE[KTJ+2].r[A_X]=xone;
    NODE[KTJ+2].r[A_Y]=0.000000;
    NODE[KTJ+2].r[A_Z]=0.000000;
    
    NODE[KTJ+3].r[A_X]=xone;
    NODE[KTJ+3].r[A_Y]=yone;
    NODE[KTJ+3].r[A_Z]=0.000000;
    
    NODE[KTJ+4].r[A_X]=0.000000;
    NODE[KTJ+4].r[A_Y]=yone;
    NODE[KTJ+4].r[A_Z]=0.000000;
    
    
    NODE[KTJ+5].r[A_X]=0.000000;
    NODE[KTJ+5].r[A_Y]=0.000000;
    NODE[KTJ+5].r[A_Z]=zone;
    
    NODE[KTJ+6].r[A_X]=xone;
    NODE[KTJ+6].r[A_Y]=0.000000;
    NODE[KTJ+6].r[A_Z]=zone;
    
    NODE[KTJ+7].r[A_X]=xone;
    NODE[KTJ+7].r[A_Y]=yone;
    NODE[KTJ+7].r[A_Z]=zone;
    
    NODE[KTJ+8].r[A_X]=0.000000;
    NODE[KTJ+8].r[A_Y]=yone;
    NODE[KTJ+8].r[A_Z]=zone;
    ///////
    
    ////要素―節点座標
    ELEM[1].node[1]=KTJ+2;
    ELEM[1].node[2]=KTJ+7;
    ELEM[1].node[3]=KTJ+5;
    ELEM[1].node[4]=KTJ+6;
    
    ELEM[2].node[1]=KTJ+1;
    ELEM[2].node[2]=KTJ+2;
    ELEM[2].node[3]=KTJ+3;
    ELEM[2].node[4]=KTJ+5;
    
    ELEM[3].node[1]=KTJ+2;
    ELEM[3].node[2]=KTJ+3;
    ELEM[3].node[3]=KTJ+5;
    ELEM[3].node[4]=KTJ+7;
    
    ELEM[4].node[1]=KTJ+5;
    ELEM[4].node[2]=KTJ+4;
    ELEM[4].node[3]=KTJ+8;
    ELEM[4].node[4]=KTJ+7;
    
    ELEM[5].node[1]=KTJ+1;
    ELEM[5].node[2]=KTJ+3;
    ELEM[5].node[3]=KTJ+4;
    ELEM[5].node[4]=KTJ+5;
    
    ELEM[6].node[1]=KTJ+4;
    ELEM[6].node[2]=KTJ+3;
    ELEM[6].node[3]=KTJ+7;
    ELEM[6].node[4]=KTJ+5;
    /////////
    
    ////要素―要素情報
    ELEM[1].elm[1]=0;
    ELEM[1].elm[2]=0;
    ELEM[1].elm[3]=0;
    ELEM[1].elm[4]=3;
    
    ELEM[2].elm[1]=3;
    ELEM[2].elm[2]=5;
    ELEM[2].elm[3]=0;
    ELEM[2].elm[4]=0;
    
    ELEM[3].elm[1]=6;
    ELEM[3].elm[2]=1;
    ELEM[3].elm[3]=0;
    ELEM[3].elm[4]=2;
    
    ELEM[4].elm[1]=0;
    ELEM[4].elm[2]=0;
    ELEM[4].elm[3]=6;
    ELEM[4].elm[4]=0;
    
    ELEM[5].elm[1]=6;
    ELEM[5].elm[2]=0;
    ELEM[5].elm[3]=2;
    ELEM[5].elm[4]=0;
    
    ELEM[6].elm[1]=3;
    ELEM[6].elm[2]=4;
    ELEM[6].elm[3]=5;
    ELEM[6].elm[4]=0;
    /////
    
    for(int i=1;i<=6;i++)
    {
        int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//体積の6倍であることに注意
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
    }
    
}

///体積計算関数
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip)
{
    ////4頂点ia,ib,ic,ipから体積をもとめる。ただし体積の6倍の値を返す(どうせ計算に必要なのはＶ＊６だから)
    
    double xa=NODE[ia].r[A_X];
    double ya=NODE[ia].r[A_Y];
    double za=NODE[ia].r[A_Z];
    
    double xb=NODE[ib].r[A_X];
    double yb=NODE[ib].r[A_Y];
    double zb=NODE[ib].r[A_Z];
    
    double xc=NODE[ic].r[A_X];
    double yc=NODE[ic].r[A_Y];
    double zc=NODE[ic].r[A_Z];
    
    double xp=NODE[ip].r[A_X];
    double yp=NODE[ip].r[A_Y];
    double zp=NODE[ip].r[A_Z];
    
    double va=xb*yc*zp+xa*ya*zp+xb*ya*za+xa*yc*za-(xb*yc*za+xa*ya*za+xb*ya*zp+xa*yc*zp);
    double vb=yb*zc*xp+ya*za*xp+yb*za*xa+ya*zc*xa-(yb*zc*xa+ya*za*xa+yb*za*xp+ya*zc*xp);
    double vc=zb*xc*yp+za*xa*yp+zb*xa*ya+za*xc*ya-(zb*xc*ya+za*xa*ya+zb*xa*yp+za*xc*yp);
    
    double wa=xb*zc*ya+xa*za*ya+xb*za*yp+xa*zc*yp-(xb*zc*yp+xa*za*yp+xb*za*ya+xa*zc*ya);
    double wb=yb*xc*za+ya*xa*za+yb*xa*zp+ya*xc*zp-(yb*xc*zp+ya*xa*zp+yb*xa*za+ya*xc*za);
    double wc=zb*yc*xa+za*ya*xa+zb*ya*xp+za*yc*xp-(zb*yc*xp+za*ya*xp+zb*ya*xa+za*yc*xa);
    
    double volume=va+vb+vc+wa+wb+wc;
    
    return volume;
}

///外接球ﾊﾟﾗﾒｰﾀ計算関数
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i)
{
    double xa=NODE[ia].r[A_X];
    double ya=NODE[ia].r[A_Y];
    double za=NODE[ia].r[A_Z];
    
    double xb=NODE[ib].r[A_X];
    double yb=NODE[ib].r[A_Y];
    double zb=NODE[ib].r[A_Z];
    
    double xc=NODE[ic].r[A_X];
    double yc=NODE[ic].r[A_Y];
    double zc=NODE[ic].r[A_Z];
    
    double xp=NODE[ip].r[A_X];
    double yp=NODE[ip].r[A_Y];
    double zp=NODE[ip].r[A_Z];
    
    double p11=yc*zp+ya*za+yp*za+ya*zc-(yc*za+ya*zp+yp*zc+ya*za);
    double p12=xp*zc+xa*za+xc*za+xa*zp-(xp*za+xa*zc+xc*zp+xa*za);
    double p13=xc*yp+xa*ya+xp*ya+xa*yc-(xc*ya+xa*yp+xp*yc+xa*ya);
    double p21=yp*zb+ya*za+yb*za+ya*zp-(yp*za+ya*zb+yb*zp+ya*za);
    double p22=xb*zp+xa*za+xp*za+xa*zb-(xb*za+xa*zp+xp*zb+xa*za);
    double p23=xp*yb+xa*ya+xb*ya+xa*yp-(xp*ya+xa*yb+xb*yp+xa*ya);
    double p31=yb*zc+ya*za+yc*za+ya*zb-(yb*za+ya*zc+yc*zb+ya*za);
    double p32=xc*zb+xa*za+xb*za+xa*zc-(xc*za+xa*zb+xb*zc+xa*za);
    double p33=xb*yc+xa*ya+xc*ya+xa*yb-(xb*ya+xa*yc+xc*yb+xa*ya);
    
    double xyza=xa*xa+ya*ya+za*za;
    double aa=0.5000000*(xb*xb+yb*yb+zb*zb-xyza);
    double bb=0.5000000*(xc*xc+yc*yc+zc*zc-xyza);
    double cc=0.5000000*(xp*xp+yp*yp+zp*zp-xyza);
    
    double xx=p11*aa+p21*bb+p31*cc;
    double yy=p12*aa+p22*bb+p32*cc;
    double zz=p13*aa+p23*bb+p33*cc;
    
    double determ=ELEM[i].volume;//体積の６倍
    double xv=xx/determ; 
    double yv=yy/determ;
    double zv=zz/determ;
    
    ELEM[i].r[A_X]=xv; //外接球中心座標
    ELEM[i].r[A_Y]=yv;
    ELEM[i].r[A_Z]=zv;
    
    ELEM[i].RR=xa*xa+xv*xv+ya*ya+yv*yv+za*za+zv*zv-2.000000*(xa*xx+ya*yy+za*zz)/determ;
}

///節点内包要素探索関数
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp)
{
    int itet=nelm;//一番最後に生成された要素番号
    int flag=1;
    while(flag==1)
    {
        flag=0;
        for(int n=1;n<=4;n++)//第１〜第４面を調べる
        {
            if(flag==0)
			{
                int i=ELEM[itet].node[n%4+1];
				int j=ELEM[itet].node[4-(n-1)/2*2];
				int k=ELEM[itet].node[3-n/2%2*2];
	
				double xi=NODE[i].r[A_X];//節点iの座標
				double yi=NODE[i].r[A_Y];
				double zi=NODE[i].r[A_Z];
	
				double xj=NODE[j].r[A_X];//節点jの座標
				double yj=NODE[j].r[A_Y];
				double zj=NODE[j].r[A_Z];
	
				double xk=NODE[k].r[A_X];//節点kの座標
				double yk=NODE[k].r[A_Y];
				double zk=NODE[k].r[A_Z];
	
				double a=yi*zj+yj*zk+yk*zi-(yi*zk+yj*zi+yk*zj);
				double b=zi*xj+zj*xk+zk*xi-(zi*xk+zj*xi+zk*xj);
				double c=xi*yj+xj*yk+xk*yi-(xi*yk+xj*yi+xk*yj);
				double d=-a*xi-b*yi-c*zi;
				if(a*xp+b*yp+c*zp+d<-ERR)
				{   //外側
					itet=ELEM[itet].elm[n];//itetの第n面に隣接している要素を次の検査対象とする
					flag=1;
					if(itet==0) return 0;	//領域の外にでた場合は0を返す。通常のデローニではありえない状況だが、たとえばremesh領域における	層生成時には、法線ベクトル上の節点位置が静的要素領域にはいるとこのようなエラーとなる
				}
			}
        }
    }
    return itet;
}

//多面体分割関数
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON)
{   
    ////この関数の前の段階で、新節点を外接球に含む四面体数ivとその要素番号kv[iv]がもとまっている。
	////また、新節点を外接球に含む四面体はELEM[i].map=1となっている
    
	int memory=CON->get_poly_memory();		//動的に確保するメモリ数
	
    int ix=0;				//表面の数 一般にiv=ixとはならない.ひとつの要素が複数の表面を担うので。(iv<=ixである)
	int *imen[3+1];
	for(int i=1;i<=3;i++) imen[i]=new int [memory]; //多面体表面三角形の節点番号格納
    int *jmen=new int [memory];		//多面体表面三角形に隣接する四面体番号格納
    int *kmen=new int [memory];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
    double *vol=new double [memory];		//多面体の体積の６倍
    
    ///////多面体表面数ixとそれを構成する頂点をもとめる
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*ivは新点を外接球内に含む要素の数
        {
			if(flag==1)
			{
        		int ielm=kv[i];//新点を外接球内に含む要素番号
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];//ielm要素に隣接する要素番号
						int ia=ELEM[ielm].node[j%4+1];//ielmとjelmの接する三角形を構成する節点番号  ia,ib,ic,ipを頂点として新しい要素が作られる
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
	                
						if(ix>49998) cout<<"pinch"<<endl;
						if(jelm==0)//表面ならielmは表面要素であるとわかる
						{
							ix++;
							imen[1][ix]=ia;	//多面体表面三角形の頂点番号
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=0;		//多面体表面三角形に隣接する要素番号
							kmen[ix]=0;		//多面体表面三角形に隣接する要素の隣接面番号  相手は0だから、ここでは0を代入
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
						}
						else if(ELEM[jelm].map==0)//ielmに接するjelmは新点を外接球に含まない。つまりielmとjelmの境界は多面体表面ということ
						{
							ix++;
						
							imen[1][ix]=ia;
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelmがielmに接する面番号
							
							if(ix>=memory)cout<<" ix>memory"<<endl;
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
							if(vol[ix]<=ERR)//教科書図3.10のように、不適切な四面体を想定してしまい、結果体積が負になる。この場合、教科書に書いてあるとおり、最初からやりなおす(goto 10)
							{   
								*iv=(*iv)+1;//リストを増やす *iv++はだめ？
								
								kv[*iv]=jelm;//本当は外接球内に新節点を含まないが、多面体分割を可能にするためにリストにいれる
								
								ELEM[jelm].map=1;
								flag=0;
							}
						}
					}
				}
			}
        }
    }
    //////多面体表面数ixとそれを構成する頂点がもとまった
   
    /////////////表面の頂点と新節点をつなげる.すると多面体からix個の要素が生成される

    int ibound=ix;//ixの代わり。つまり表面の数を表す
    
    for(int i=*(iv)+1;i<=ibound;i++) ///iv個の要素が破壊されて新たにibound個の要素が生成されるから、増える要素数はibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++という書き方ではだめ？
        kv[i]=*nelm;		//新節点を外接球に含む要素リストも増える。
		ELEM[*nelm].map=1;	//生成される要素は必ず新節点を外接球に含む(としている).(この行意味なくない？)
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//マッピングの初期化
    
    for(int i=1;i<=ibound;i++)//要素情報生成
    {   
        int ielm=kv[i];
		double determ=vol[i];
		
		int ia=imen[1][i];
		int ib=imen[2][i];
		int ic=imen[3][i];
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;		//新点は４番目と定義	
		ELEM[ielm].elm[4]=jmen[i];	
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		///これでいい？
		ELEM[ielm].volume=vol[i];
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//外接球の中心（ボロノイ点)と半径の二乗を計算
    }
    ///////////////////
    
    //要素―要素関係修正/////上の処理で第4面で接する要素番号はわかっているので、残りを求める
	//						ここで、1〜3面は多面体を構成する要素との境界面であることに注意
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]はすでにもとまったから、それ以外をもとめる
		{
			///ELEM[ielm].node[4]=ipである
			ELEM[ielm].elm[j]=-1;				//初期化
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3のとき、2,3,1の順
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3のとき、3,1,2の順
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					int ja=imen[1][k];
					int jb=imen[2][k];
					if(ia==ja && ib==jb)//節点が一致したら
					{
					    ELEM[ielm].elm[j]=jmen[k];//あらかじめﾘｽﾄしてあった情報を格納
					    ELEM[jmen[k]].elm[kmen[k]]=ielm;
					    imen[1][k]=imen[1][ix];//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
					    imen[2][k]=imen[2][ix];
					    jmen[k]=jmen[ix];
					    kmen[k]=kmen[ix];
					    ix--;	//待ち辺数減少
						flag=1;	//ELEM[ielm].elm[j]はもとまったので、下のネストに入る必要はないのでflag=1
					}
				}
			}
			if(flag==0)
			{
			    ix++;			//ここでのixは、[隣接関係を満たす要素]をまっている[辺]の数を表す。
				imen[1][ix]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
				imen[2][ix]=ia;
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///要素-要素関係修正完了

    /////////新たに作られた四面体の個数が、古いものより少なくなった場合(多面体を構成したとき、どの面も多面体の境界ではない内部要素が存在するとき)
    if(*iv>ibound)
    {   
	
        int ir=*(iv)-ibound;	//ir個の要素が削除されなくてはならない.ただし普通に削除したのでは、要素配列に[穴]が生じることになる

		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];//mapとしてとりあえずその要素番号を格納
		}
		///上で代入したmapの値(要素番号)を小さい順にならびかえる。その都合上、上ではkv[1],kv[2]・・・と値をつめている
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//iが増えるにつれ、ir-i+1はirから1ずつ小さくなっていく
			ELEM[ielm].map=0;	//上のほうでmapの初期化を行っているが、iv>iboundの場合、ir個の要素は初期化されていないので、ここで初期化
			
			if(ielm!=*nelm)		//ielm番目の要素に、nelm番目の要素情報を上書き
			{   
			    ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				for(int j=1;j<=4;j++)
				{   
				    ELEM[ielm].node[j]=ELEM[*nelm].node[j];
				    int jelm=ELEM[*nelm].elm[j];
				    ELEM[ielm].elm[j]=jelm;
				    if(jelm!=0)
				    {   
				        int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
				    }
				}
			}
			*nelm=*nelm-1;
		}
    }///////////////////*/

	delete [] jmen;
	delete [] kmen;
	delete [] vol;

	for(int i=1;i<=3;i++) delete [] imen[i];
	
}		

////隣接面番号計算関数
int iface3D(vector<element3D> &ELEM,int ielm,int jelm)
{
    int iface=-1;//returnする変数
    for(int n=1;n<=4;n++)
    {
        if(ELEM[jelm].elm[n]==ielm)
		{
			iface=n;
			// break;
		}
    }
    if(iface<0) cout<<"error in inface"<<endl;
    return iface;
}

///並べ替え関数
void qsorti3D(vector<point3D> &NODE,vector<element3D> &ELEM,int n,int *list)
{
    ///list=kv[i]
    int maxstk=100;//パラメータ
    
    int ll=1;
    int lr=n;//並び替える個数
    int istk=0;
    
    int *ilst=new int[maxstk];
    int *irst=new int[maxstk];
    
    int flag2=0;
    while(flag2==0)  //10 continue
    {
        flag2=1;
        while(ll<lr)
        {
            int nl=ll;
			int nr=lr;
			int lm=(ll+lr)/2;
			int iguess=ELEM[list[lm]].map;
	
			int flag=0;
			while(flag==0)
			{   
				flag=1;
				while(ELEM[list[nl]].map<iguess)
				{
					nl++;
				}
	
				while(iguess<ELEM[list[nr]].map)
				{
					nr--;
				}
	
				if(nl<nr-1)
				{
					int ltemp=list[nl];
					list[nl]=list[nr];
					list[nr]=ltemp;
					nl++;
					nr--;
				flag=0;
				}
			}
	
			if(nl<=nr)
            {
                if(nl<nr)
				{
					int ltemp=list[nl];
					list[nl]=list[nr];
					list[nr]=ltemp;
				}
				nl++;
				nr--;
            }
    
            istk++;
			if(istk>maxstk) cout<<"ERROR(qsorti) istk>maxstk"<<endl;
	
			if(nr<lm)
			{
				ilst[istk]=nl;
				irst[istk]=lr;
				lr=nr;
			}
			else
			{
				ilst[istk]=ll;
				irst[istk]=nr;
				ll=nl;
			}
	
        }
        if(istk!=0)
        {
            ll=ilst[istk];
			lr=irst[istk];
			istk--;
			flag2=0;
        }
    }
    delete [] ilst;
    delete [] irst;
}

///不要要素除去関数
void remove3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,int iv,int *kv)
{
	//iv:除去要素数　kv[i]:除去要素番号

    int m=0;///除去されない要素数 計算の都合上2種類必要
    int n=0;///除去されない要素数
    for(int i=1;i<=*nelm;i++) ELEM[i].map=1;	//初期化 ここではmapは消去後の要素番号を表す
    
    ///////各要素のmap(新しい要素番号)をもとめる

    for(int i=1;i<=iv;i++) ELEM[kv[i]].map=0;//消去されるのだから0番目となる
    
    for(int i=1;i<=*nelm;i++)
    {
        if(ELEM[i].map!=0)//除去されないなら
		{
			m++;
			ELEM[i].map=m;
		}
    }/////////
    
    ////各情報のひきつぎ　
    for(int i=1;i<=*nelm;i++)
    {
        if(ELEM[i].map!=0)//除去されないなら
		{
			n++;
			ELEM[n].r[A_X]=ELEM[i].r[A_X];
			ELEM[n].r[A_Y]=ELEM[i].r[A_Y];
			ELEM[n].r[A_Z]=ELEM[i].r[A_Z];
			ELEM[n].RR=ELEM[i].RR;
			for(int ia=1;ia<=4;ia++)
			{
				//要素-節点情報のコピー
				ELEM[n].node[ia]=ELEM[i].node[ia];

				//要素-要素情報のコピー
				if(ELEM[i].elm[ia]==0) ELEM[n].elm[ia]=0;	//表面と接しているならその情報を単純にコピー
				else //そうでないなら、要素番号の変化を考慮にいれたコピーを行わないといけない
				{
					ELEM[n].elm[ia]=ELEM[ELEM[i].elm[ia]].map;//ここで、ｽｰﾊﾟｰﾎﾞｯｸｽ頂点を含む要素(不要要素)はmap=0だから、不要要素に隣接する要素は隣接要素としてゼロが格納される
				}
			}    
		}
    }///////////
    
	//不要要素の情報初期化
    for(int i=n+1;i<=*nelm;i++)
    {
        ELEM[i].r[A_X]=0.00000;
		ELEM[i].r[A_Y]=0.00000;
		ELEM[i].r[A_Z]=0.00000;
		ELEM[i].RR=0.00000;
		for(int ia=1;ia<=4;ia++)
		{
			ELEM[i].node[ia]=0;
			ELEM[i].elm[ia]=0;
		}
    }
    
    *nelm=*nelm-iv; 
}

///要素生成確認関数
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm)
{
    for(int i=1;i<=nelm;i++)
    {
        int ielm=i;
		for(int j=1;j<=4;j++)
		{
			int ia=ELEM[ielm].node[j%4+1];
			int ib=ELEM[ielm].node[4-(j-1)/2*2];
			int ic=ELEM[ielm].node[3-j/2%2*2];
			int jelm=ELEM[ielm].elm[j];
			
			if(jelm!=0 && ielm<jelm)
			{
				int k=iface3D(ELEM,ielm,jelm);
				int ja=ELEM[jelm].node[k%4+1];
				int jb=ELEM[jelm].node[4-(k-1)/2*2];
				int jc=ELEM[jelm].node[3-k/2%2*2];
				int flag=0;
				if(ia==jc && ib==jb && ic==ja) flag=1;
				if(ib==jc && ic==jb && ia==ja) flag=1;
				if(ic==jc && ia==jb && ib==ja) flag=1;
				if(flag==0) cout<<"ERROR IN FILL"<<endl;
			}
		}
    }
}

///要素材質決定関数
void set_material(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm)
{
    for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==AIR)//リメッシュ領域内部の要素は、ひとまずAIRと判定されている
		{
			int M1=NODE[ELEM[i].node[1]].material;
			int M2=NODE[ELEM[i].node[2]].material;
			int M3=NODE[ELEM[i].node[3]].material;
			int M4=NODE[ELEM[i].node[4]].material;
		
			///4頂点すべてが同じ材質なら要素もそれにならう。
			///ひとつでも異なっていたら空気と定義
			if(M1==M2 && M2==M3 && M3==M4) ELEM[i].material=M1;
			else if(M1!=AIR && M2!=AIR && M3!=AIR && M4!=AIR) 
			{
				if(CON->get_model_number()==15) ELEM[i].material=AIR;
				else ELEM[i].material=FLUID;
			}
			else ELEM[i].material=AIR;
		}
    }
}

//ﾒｯｼｭﾃﾞｰﾀ出力関数
void data_avs(int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double *val,mpsconfig *CON)
{
	ofstream fq("mesh1.inp");

    cout<<"Now writing-----";	
	
	//mesh1
	fq<<"# Micro AVS"<<endl;
	fq<<"1"<<endl;
	fq<<"data"<<endl;
	fq<<"step1"<<endl;
	//fq<<node_number+8<<" "<<nelm<<endl;//こっちはｽｰﾊﾟｰﾎﾞｯｸｽも表示するとき
	fq<<node_number<<" "<<nelm<<endl;

	/*for(i=ktj+1;i<=ktj+8;i++) //こっちはｽｰﾊﾟｰﾎﾞｯｸｽも表示するとき
	{
		fq<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
		//fprintf(fp, "\t%d\t%f\t%f\t%f\n", i, NODE[i].r[A_X], NODE[i].r[A_Y], NODE[i].r[A_Z]);
		
	} */

	//節点番号とその座標の出力
	for(int i=1;i<=node_number;i++) fq<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=1;i<=nelm;i++)
	{
		fq<<i<<" t0 tet ";
		for(int j=1;j<=4;j++) fq<<ELEM[i].node[j]<<" ";
		fq<<endl;
	}
	
	fq<<"1 0"<<endl;
	fq<<"1 1"<<endl;
	fq<<"element, e"<<endl;

    /*for(int i=ktj+1;i<=ktj+8;i++) //こっちはｽｰﾊﾟｰﾎﾞｯｸｽも表示するとき
	{
		fq<<i<<" 0"<<endl;
	}*/

	//各節点の値代入
	for(int i=1;i<=node_number;i++) fq<<i<<" "<<val[i]<<endl;		//出力したいスカラー値をvalに入力しておくこと
	
	cout<<"OK"<<endl;
	fq.close();
}

//ﾒｯｼｭﾃﾞｰﾀ出力関数ver.2 ある断面のﾒｯｼｭをみたいときに使う
void data_avs2(mpsconfig *CON,int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,int t)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	cout<<"Now ver.2a writing-----";
	int step=CON->get_step()/(CON->get_EM_interval()*CON->get_mesh_output_interval())+1;				//出力する総ステップ数
	int nelm2=0;			//出力する要素数
	int *node_flag=new int[node+1];			//各節点を出力するか、しないか
	int node_on=0;							//出力する節点数
	
	for(int i=1;i<=nelm;i++) ELEM[i].map=0;//ﾏｯﾋﾟﾝｸﾞ初期化(mapは要素作成の際に利用していたが、ここでもそのﾒﾓﾘを使う)
	for(int i=0;i<=node;i++) node_flag[i]=OFF;		//初期化

	for(int i=1;i<=nelm;i++)
	{
		double Y=0;//要素の重心のy座標
		for(int j=1;j<=4;j++) Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
		if(Y<0)
		//if(Y<0 && ELEM[i].remesh==ON)
		{
			nelm2++;
			ELEM[i].map=1;//出力するという印
			for(int j=1;j<=4;j++) node_flag[ELEM[i].node[j]]=ON;//関係する節点のフラグをON
		}
	}
	///出力する節点と要素、およびその数がもとまった

	for(int i=1;i<=node;i++) if(node_flag[i]=ON) node_on++;

	if(t==1)
	{
		ofstream fout("mesh2.inp");

		fout<<step<<endl;
		fout<<"data_geom"<<endl;
		
		fout.close();
	}

	ofstream fp("mesh2.inp",ios :: app);
	if(CON->get_EM_interval()==1) fp<<"step"<<t/(CON->get_EM_interval()*CON->get_mesh_output_interval())<<endl;
	else if(CON->get_EM_interval()>1) fp<<"step"<<t/(CON->get_EM_interval()*CON->get_mesh_output_interval())+1<<endl;
	
	fp<<node_on<<" "<<nelm2<<endl;	//節点は無関係のものも出力していいけど、それだとファイルが重くなるから、必要な節点だけ出力 
	
	
	
	//節点番号とその座標の出力
	for(int i=1;i<=node;i++)//節点は無関係のものも出力していいけど、それだとファイルが重くなるから、必要な節点だけ出力 
	{
		//fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
		if(node_flag[i]==ON) fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==1)
		{
			fp<<i<<"  0 tet ";
			for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
			fp<<endl;
		}
	}

	fp<<"0 1"<<endl;//節点の情報量がゼロで、要素の情報量がnum_infoということ。
	fp<<"1 1"<<endl;
	fp<<"material, material"<<endl;

	//各節点の値代入
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	//for(int i=1;i<=node;i++) fp<<i<<" "<<val[i]<<endl;

	int mate=0;
	
		for(int i=1;i<=nelm;i++)
		{
			if(ELEM[i].material==10000) mate=AIR;
			else mate=ELEM[i].material;
			if(ELEM[i].map==1) fp<<i<<"  "<<mate<<endl;
		}
	//for(int i=1;i<=nelm;i++) if(ELEM[i].map==1) fp<<i<<"  "<<ELEM[i].material<<endl;
	
	cout<<"OK"<<endl;
	fp.close();
	delete [] node_flag;
}


void data_avs2flux(mpsconfig *CON,int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,double *val,int t,int flux)
{
	int nelm2=0;//出力する要素数
	double Ymin=CON->get_YD();

	for(int i=1;i<=nelm;i++) ELEM[i].map=0;//ﾏｯﾋﾟﾝｸﾞ初期化(mapは要素作成の際に利用していたが、ここでもそのﾒﾓﾘを使う)
	
	for(int i=1;i<=nelm;i++)
	{
		double X=0;//要素の重心のx座標
		for(int j=1;j<=4;j++) X+=NODE[ELEM[i].node[j]].r[A_X]*0.25;
		double Y=0;//要素の重心のy座標
		for(int j=1;j<=4;j++) Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
		//if(Y<0 && Y>Ymin*0.1)
		//if(Y<0 && Y>-1.5*CON->get_distancebp())
		//if(Y<0)
		if(flux==0)
		{
			if(X<sin(PI/24)*Y)
			{
				nelm2++;
				ELEM[i].map=1;//出力するという印
			}
		}
		else if(flux==1)
		{
			if(Y<sin(PI/24)*X)
			{
				nelm2++;
				ELEM[i].map=1;//出力するという印
			}
		}
	}
	///出力する節点と要素、およびその数がもとまった


	cout<<"Now ver.2 writing-----";
	char filename[30];
	sprintf_s(filename,"mesh2/mesh2flux%d_%d.inp",flux,t);
	ofstream fp(filename);
	//ofstream fp("mesh2.inp");
	
	fp<<"# Micro AVS"<<endl;
	fp<<"1"<<endl;
	fp<<"data"<<endl;
	fp<<"step1"<<endl;
	fp<<node_number<<" "<<nelm2<<endl;
	
	
	
	//節点番号とその座標の出力
	for(int i=1;i<=node_number;i++)//節点は無関係のものも出力すればいい 
	{
		fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==1)
		{
			fp<<i<<"  0 tet ";
			for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
			fp<<endl;
		}
	}
	
	//fp<<"1 0"<<endl;
	//fp<<"1 1"<<endl;
	//fp<<"element, e"<<endl;

	fp<<"0 1"<<endl;
	fp<<"1 1"<<endl;
	fp<<"material, material"<<endl;

	//各節点の値代入
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<val[i]<<endl;

	for(int i=1;i<=nelm;i++) if(ELEM[i].map==1) fp<<i<<"  "<<val[i]<<endl;
	
	cout<<"OK"<<endl;
	fp.close();
}

//ﾒｯｼｭﾃﾞｰﾀ出力関数ver.2 ある断面のﾒｯｼｭをみたいときに使う
void data_avs2node(mpsconfig *CON,int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,double *val,int t)
{
	int nelm2=0;//出力する要素数
	double Ymin=CON->get_YD();

	for(int i=1;i<=nelm;i++) ELEM[i].map=0;//ﾏｯﾋﾟﾝｸﾞ初期化(mapは要素作成の際に利用していたが、ここでもそのﾒﾓﾘを使う)
	
	for(int i=1;i<=nelm;i++)
	{
		double Y=0;//要素の重心のy座標
		for(int j=1;j<=4;j++) Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
		//if(Y<0 && Y>Ymin*0.1)
		//if(Y<0 && Y>-1.5*CON->get_distancebp())
		if(Y<0)
		{
			nelm2++;
			ELEM[i].map=1;//出力するという印
		}
	}
	///出力する節点と要素、およびその数がもとまった


	cout<<"Now ver.2 writing-----";
	char filename[30];
	sprintf_s(filename,"mesh2n/mesh2node_%d.inp",t);
	ofstream fp(filename);
	//ofstream fp("mesh2.inp");
	
	fp<<"# Micro AVS"<<endl;
	fp<<"1"<<endl;
	fp<<"data"<<endl;
	fp<<"step1"<<endl;
	fp<<node_number<<" "<<nelm2<<endl;
	
	
	
	//節点番号とその座標の出力
	for(int i=1;i<=node_number;i++)//節点は無関係のものも出力すればいい 
	{
		fp<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}
	
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==1)
		{
			fp<<i<<"  0 tet ";
			for(int j=1;j<=4;j++)	fp<<ELEM[i].node[j]<<" ";
			fp<<endl;
		}
	}
	
	fp<<"1 0"<<endl;
	fp<<"1 1"<<endl;
	fp<<"element, e"<<endl;

	//各節点の値代入
	//for(int i=1;i<=node_number;i++) fp<<i<<" "<<NODE[i].material<<endl;
	for(int i=1;i<=node_number;i++) fp<<i<<" "<<val[i]<<endl;
	
	cout<<"OK"<<endl;
	fp.close();
}

//ﾒｯｼｭﾃﾞｰﾀ出力関数(材質ver）
void data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig *CON,int t)
{
	//ofstream fout("mesh3.mgf");
	char filename[30];
	sprintf_s(filename,"mesh3/mesh3_%d.mgf",t);
	ofstream fout(filename);
	cout<<"Now ver.3 writing-----";
	int step=1;
	fout <<"# Micro AVS Geom:2.00" << endl;
	//fout << "step" << step<< endl;
	fout << "disjoint polygon" << endl;
	fout <<"element" << endl;
	fout <<"facet" << endl;
	fout <<"color" << endl;

	int count=0;//物質表面数
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		//if(ELEM[i].material==FLUID)
		//if(ELEM[i].material==COIL)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR) count++;
			}
		}
	}///表面数がもとまった
	
	fout <<count<< endl;//要素数出力
	double red, green, blue;
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material!=AIR)
		//if(ELEM[i].material==FLUID)
		//if(ELEM[i].material==COIL)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)
				{
					if(ELEM[i].material==FLUID)
					{
						red   = 0.0;
						green = 0.0;
						blue  = 1.0;
					}
					else if(ELEM[i].material==ELECTRODE || ELEM[i].material==MAGNET)
					{
						red   = 1.0;
						green = 0.0;
						blue  = 0.0;
					}
					else
					{
						red   = 0.0;
						green = 1.0;
						blue  = 0.0;

						red   = 1.0;
						green = 1.0;
						blue  = 1.0;
					}
					
					fout <<"3"<< endl;
					int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];
					for(int D=0;D<3;D++) fout << NODE[ia].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ib].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;

					for(int D=0;D<3;D++) fout << NODE[ic].r[D]<<" ";
					fout << red << " " << green << " " << blue;
					fout << endl;
				}
			}
		}
	}////*/

	fout.close();
	cout<<"OK"<<endl;
}

///自動メッシュ生成関数
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig *CON,double rrm,int startID)
{
	
	cout<<"要素の再分割実行 節点数:"<<*node<<" 要素数:"<<*nelm<<endl;

	unsigned int timeA=GetTickCount();	//計算開始時刻

	/////要素の体積とﾎﾞﾛﾉｲ点をもとめておく(スーパーボックスのそとで関数が呼び出されても対応できるように)
	for(int i=1;i<=*nelm;i++)
    {   
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//体積の6倍であることに注意
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);//外接球の中心（ボロノイ点)と半径の二乗を計算
    }
	//////////////*/

	//空気層生成
	int layer_node_num=0;
	if(CON->get_layer_num()>0) layer_node_num=make_air_layer(NODE,ELEM,nelm,CON,node, KTE,rrm,KTJ);	//layer_node_numはこの関数により増加した節点数

	int Mnum=10;								//歪だと判定され、再分割された要素数
	int newnode=0;								//新しく導入した節点数
	int limit_num=CON->get_add_points()-layer_node_num;//導入できる最大節点数
	//if(limit_num<0) limit_num=0;
	int *kv=new int[KTE];						//新節点を外接球にふくむ要素群
    int *istack=new int[KTE];					//一時配列

	int MESH;
	double limit=CON->get_co_fine();	//辺ﾍﾞｰｽの場合は最少辺と最大辺の比率のしきい値
	int remesh;

	//重心ベースで使う値を求める
	double Gp[3];
	for(int D=0;D<3;D++) Gp[D]=0;
	double volf=0;
	int countf=0;
	for(int i=1;i<=*nelm;i++)
	{
		if(ELEM[i].material==FLUID)
		{
			countf++;
			volf+=ELEM[i].volume;
		}
	}
	volf/=countf;
	cout<<"目標体積="<<volf<<endl;

	countf=0;
	for(int i=1;i<=*node;i++)
	{
		if(NODE[i].material==FLUID)
		{
			countf++;
			for(int D=0;D<3;D++) Gp[D]+=NODE[i].r[D];
		}
	}
	for(int D=0;D<3;D++) Gp[D]/=countf;
	cout<<"流体重心座標="<<Gp[A_X]<<" "<<Gp[A_Y]<<" "<<Gp[A_Z]<<endl;

	
	////////////
	if(CON->get_fine()==1 && limit_num>0)//辺ﾍﾞｰｽ
	{
		while(Mnum>0)
		{
			Mnum=0;			
		
			///初期化
			for(int i=1;i<KTE;i++) ELEM[i].map=0;
	
			///順次節点を導入していく
	    
			MESH=*nelm;
		
			for(int je=startID;je<=MESH;je++)
			{
				//cout<<je<<endl;
				if(newnode<limit_num)
				{
					//if(ELEM[je].material==AIR && depth[je]>=1)
					//if(ELEM[je].material==AIR)
					if(ELEM[je].material==AIR && ELEM[je].volume>volf/10)//執拗にメッシュを細分化しているところがあるので、目標体積程度で打ち切ってみる
					//if(ELEM[je].material==AIR && ELEM[je].volume>1e-9)//執拗にメッシュを細分化しているところがあるので、目標体積程度で打ち切ってみる
					{
						int flag3=0;//1なら新節点を導入する箇所が決まったということ
						double r[3];//新節点の座標格納
						for(int j=1;j<=4;j++)
						{
							if(flag3==0 && ELEM[je].elm[j]!=0)
							{
								int ia=ELEM[je].node[j%4+1];//ielmとjelmの接する三角形を構成する節点番号  ia,ib,ic,ipを頂点として新しい要素が作られる
								int ib=ELEM[je].node[4-(j-1)/2*2];
								int ic=ELEM[je].node[3-(j/2%2)*2];
		
								double iaib=0;//辺ia,ibの距離
								double iaic=0;//辺ia,icの距離
								double ibic=0;//辺ib,icの距離
								for(int D=0;D<3;D++)
								{
									iaib+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);
									iaic+=(NODE[ia].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
									ibic+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);
								}
								iaib=sqrt(iaib);
								iaic=sqrt(iaic);
								ibic=sqrt(ibic);
								double minL=iaib;//最小辺長さ
								double maxL=iaib;//最大辺長さ
								int N[2]={ia,ib};//最大辺を構成する節点番号
								if(iaic<minL) minL=iaic;
								else if(iaic>maxL) 
								{
									maxL=iaic;
									N[0]=ia;
									N[1]=ic;
								}
								if(ibic<minL) minL=ibic;
								else if(ibic>maxL)
								{
									maxL=ibic;
									N[0]=ib;
									N[1]=ic;
								}
								if(maxL>limit*minL) //最大長さが最小のlimit倍以上あれば
								{
									int n1=N[0];//最大辺を構成する節点番号
									int n2=N[1];
									if(NODE[n1].boundary_condition==0 && NODE[n2].boundary_condition==0)//未知数なら
									{
										flag3=1;
										for(int D=0;D<3;D++)  r[D]=0.5*(NODE[n1].r[D]+NODE[n2].r[D]);	//新点は最大辺の中点に設置
										
										if(NODE[n1].remesh==NODE[n2].remesh) remesh=NODE[n1].remesh;	//remeshが同じならそれを引き継ぐ。異なるならそれは非remesh領域内の節点なのでremesh=OFF
										else remesh=OFF;
									}
								}
							}
						}
						if(flag3==1)
						{
							Mnum++;
							newnode++;		//新点数増加
							*node=*node+1;
							
							int ip=*node;
							for(int D=0;D<3;D++) NODE[ip].r[D]=r[D];
							NODE[ip].boundary_condition=0;
							NODE[ip].material=AIR;
							NODE[ip].remesh=remesh;
							NODE[ip].BD_node=OFF;
							NODE[ip].particleID=-1;		//対応する粒子は存在しない
							double xp=NODE[ip].r[A_X];//導入する節点の座標
							double yp=NODE[ip].r[A_Y];
							double zp=NODE[ip].r[A_Z];
			
							///新節点を含む要素の探索
							
							int loc=je;//明らかにjeは新点を外接球に含む
							
							//////////外接球内に新節点を含む要素の抽出
							int iv=0;
							int msk=0;
			
							iv++;//外接球内に新節点を含む要素数
							kv[iv]=loc;
							ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むということ
							msk++;
							istack[msk]=loc;
							
							if(CON->get_CDT_sw()==OFF)//通常のデローニ
							{
								while(msk!=0)
								{	   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素
									
										if(jelm!=0)//それが表面でないなら
										{
											if(ELEM[jelm].map==0) //まだ検査してないなら
											{   
												//double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
												double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
												double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
												if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
												{
													dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
													if(dst<rad)
													{
														dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
														if(dst<rad)//外接球内に含む
														{
															if(ELEM[jelm].material==FLUID) cout<<"流体要素がFINEの対象となりました"<<endl;
															if(ELEM[jelm].remesh==OFF) cout<<"静的要素がFINEの対象となりました"<<endl;
															iv++;//外接球内に新節点を含む要素数を＋１
															kv[iv]=jelm;//リストにいれる
															msk++;
															istack[msk]=jelm;
															ELEM[jelm].map=1;//jelmは外接球内に新節点を含む
														}
													}
												}
											}
										}
									}
								}
							}//新点を外接球内に含む要素数ivとその要素番号kv[iv]がもとまった

							if(CON->get_CDT_sw()==ON)//制約付きデローニ
							{
								while(msk!=0)
								{	   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素

										if(ELEM[jelm].remesh==ON)
										{
											if(jelm!=0)//それが表面でないなら
											{
												if(ELEM[jelm].map==0) //まだ検査してないなら
												{   
													//double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
													double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
													double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
													if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
													{
														dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
														if(dst<rad)
														{
															dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
															if(dst<rad)//外接球内に含む
															{
																if(ELEM[jelm].remesh==ON && ELEM[jelm].material==AIR) //後で多面体を分割するときに足される要素が流体だと弾けない
																{
																	iv++;//外接球内に新節点を含む要素数を＋１
																	kv[iv]=jelm;//リストにいれる
																	msk++;
																	istack[msk]=jelm;
																	ELEM[jelm].map=1;//jelmは外接球内に新節点を含む
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}//新点を外接球内に含む要素数ivとその要素番号kv[iv]がもとまった





							/////////////////
							int NN=*nelm;//分割前の要素数
		
							////得られた多面体を四面体に分割し、新しい要素を生成する。また、関数内で新要素の材質を決定する
							int FLAG=poly3D_for_FINE3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
							if(FLAG==OFF)		//polyが失敗した場合
							{
								*node=*node-1;	//新節点導入をあきらめる
								Mnum--;
								newnode--;
							}
						}
					}
				}
				else Mnum=0;//newnodeがlimit_numに等しくなったらMnum=0にしてﾙｰﾌﾟから脱出
			}
		}
	}
	////////////

	if(CON->get_fine()==2 && limit_num>0)//重心ﾍﾞｰｽ
	{
		while(Mnum>0)
		{
			Mnum=0;			
		
			///初期化
			for(int i=1;i<KTE;i++) ELEM[i].map=0;
	
			///順次節点を導入していく
	    
			MESH=*nelm;
		
			for(int je=startID;je<=MESH;je++)
			{
				//cout<<je<<endl;
				if(newnode<limit_num)
				{
					//if(ELEM[je].material==AIR && depth[je]>=1)
					if(ELEM[je].material==AIR)
					{
						int flag3=0;//1なら新節点を導入する箇所が決まったということ
						double r[3];//新節点の座標格納
						for(int j=1;j<=4;j++)
						{
							if(flag3==0 && ELEM[je].elm[j]!=0)
							{
								int ia=ELEM[je].node[j%4+1];
								int ib=ELEM[je].node[4-(j-1)/2*2];
								int ic=ELEM[je].node[3-(j/2%2)*2];
								int id=ELEM[je].node[j];
								for(int D=0;D<3;D++) r[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D]+NODE[id].r[D])/4;	//新点は重心に
								double dx=r[A_X]-Gp[A_X];
								double dy=r[A_Y]-Gp[A_Y];
								double dz=r[A_Z]-Gp[A_Z];
								double dis=sqrt(dx*dx+dy*dy+dz*dz);
								double w =dis/0.03;
								//if(w<1) w=1;

								//if(ELEM[je].volume>volf*4*w*w)//2*w*w
								if(ELEM[je].volume>5e-9)//2*w*w
								{
									flag3=1;
									int remesh=ON;
								}
							}
						}
						if(flag3==1)
						{
							Mnum++;
							newnode++;		//新点数増加
							*node=*node+1;
							
							int ip=*node;
							for(int D=0;D<3;D++) NODE[ip].r[D]=r[D];
							NODE[ip].boundary_condition=0;
							NODE[ip].material=AIR;
							NODE[ip].remesh=remesh;
							NODE[ip].BD_node=OFF;
							NODE[ip].particleID=-1;		//対応する粒子は存在しない
							double xp=NODE[ip].r[A_X];//導入する節点の座標
							double yp=NODE[ip].r[A_Y];
							double zp=NODE[ip].r[A_Z];
			
							///新節点を含む要素の探索
							
							int loc=je;//明らかにjeは新点を外接球に含む
							
							//////////外接球内に新節点を含む要素の抽出
							int iv=0;
							int msk=0;
			
							iv++;//外接球内に新節点を含む要素数
							kv[iv]=loc;
							ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むということ
							msk++;
							istack[msk]=loc;
							
							if(CON->get_CDT_sw()==OFF)//通常のデローニ
							{
								while(msk!=0)
								{	   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素
									
										if(jelm!=0)//それが表面でないなら
										{
											if(ELEM[jelm].map==0) //まだ検査してないなら
											{   
												//double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
												double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
												double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
												if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
												{
													dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
													if(dst<rad)
													{
														dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
														if(dst<rad)//外接球内に含む
														{
															if(ELEM[jelm].material==FLUID) cout<<"流体要素がFINEの対象となりました"<<endl;
															if(ELEM[jelm].remesh==OFF) cout<<"静的要素がFINEの対象となりました"<<endl;
															iv++;//外接球内に新節点を含む要素数を＋１
															kv[iv]=jelm;//リストにいれる
															msk++;
															istack[msk]=jelm;
															ELEM[jelm].map=1;//jelmは外接球内に新節点を含む
														}
													}
												}
											}
										}
									}
								}
							}//新点を外接球内に含む要素数ivとその要素番号kv[iv]がもとまった

							if(CON->get_CDT_sw()==ON)//制約付きデローニ
							{
								while(msk!=0)
								{	   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素

										if(ELEM[jelm].remesh==ON)
										{
											if(jelm!=0)//それが表面でないなら
											{
												if(ELEM[jelm].map==0) //まだ検査してないなら
												{   
													//double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
													double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
													double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
													if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
													{
														dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
														if(dst<rad)
														{
															dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
															if(dst<rad)//外接球内に含む
															{
																if(ELEM[jelm].material==AIR && ELEM[jelm].remesh==ON) //後で多面体を分割するときに足される要素が流体だと弾けない
																{
																	iv++;//外接球内に新節点を含む要素数を＋１
																	kv[iv]=jelm;//リストにいれる
																	msk++;
																	istack[msk]=jelm;
																	ELEM[jelm].map=1;//jelmは外接球内に新節点を含む
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}//新点を外接球内に含む要素数ivとその要素番号kv[iv]がもとまった

							/////////////////
							int NN=*nelm;//分割前の要素数
		
							////得られた多面体を四面体に分割し、新しい要素を生成する。また、関数内で新要素の材質を決定する
							int FLAG=poly3D_for_FINE3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
							if(FLAG==OFF)		//polyが失敗した場合
							{
								*node=*node-1;	//新節点導入をあきらめる
								Mnum--;
								newnode--;
							}
						}
					}
				}
				else Mnum=0;//newnodeがlimit_numに等しくなったらMnum=0にしてﾙｰﾌﾟから脱出
			}
		}
	}	

	cout<<"節点+="<<newnode<<endl;
	cout<<"FINE完了 time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	

	delete [] kv;
	delete [] istack;
}

//FINE用多面体分割関数
int poly3D_for_FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON)
{   
    ////この関数の前の段階で、新節点を外接球に含む四面体数ivとその要素番号kv[iv]がもとまっている。
    ////また、新節点を外接球に含む四面体はELEM[i].map=1となっている

	int ix=0;				//表面の数. 一般にiv=ixとはならない.ひとつの要素が複数の表面を担うので。(iv<=ixである)
    int imen[20000][3+1];	//多面体表面三角形の節点番号格納
    int jmen[20000];		//多面体表面三角形に隣接する四面体番号格納
    int kmen[20000];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
    double vol[20000];		//多面体の体積の６倍
	
	/*
	int memory=CON->get_poly_memory();		//動的に確保するメモリ数
	int ix=0;				//表面の数 一般にiv=ixとはならない.ひとつの要素が複数の表面を担うので。(iv<=ixである)
	int *imen[3+1];
	for(int i=1;i<=3;i++) imen[i]=new int [memory]; //多面体表面三角形の節点番号格納
    int *jmen=new int [memory];		//多面体表面三角形に隣接する四面体番号格納
    int *kmen=new int [memory];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
    double *vol=new double [memory];		//多面体の体積の６倍
    */

    ///////多面体表面数ixとそれを構成する頂点をもとめる
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*ivは新点を外接球内に含む要素の数
        {
			if(flag==1)
			{
        		int ielm=kv[i];//新点を外接球内に含む要素番号
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];		//ielm要素に隣接する要素番号
						int ia=ELEM[ielm].node[j%4+1];	//ielmとjelmの接する三角形を構成する節点番号  ia,ib,ic,ipを頂点として新しい要素が作られる
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
							
						if(jelm==0)//表面ならielmは表面要素であるとわかる
						{
							ix++;
							
							/*
							imen[1][ix]=ia;	//多面体表面三角形の頂点番号
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							*/
						
							imen[ix][1]=ia;	//多面体表面三角形の頂点番号
							imen[ix][2]=ib;
							imen[ix][3]=ic;
				
							jmen[ix]=0;		//多面体表面三角形に隣接する要素番号
							kmen[ix]=0;		//多面体表面三角形に隣接する要素の隣接面番号  相手は0だから、ここでは0を代入
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);

							//if(vol[ix]<ERR)//この行を消したほうがうまくいく？気のせい？
							if(vol[ix]<=1.0e-18)
							{
								
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
								return OFF;//新節点導入をあきらめる
								
							}
						}
						else if(ELEM[jelm].map==0)//ielmに接するjelmは新点を外接球に含まない。つまりielmとjelmの境界は多面体表面ということ
						{
							ix++;
							/*
							imen[1][ix]=ia;
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							*/
							
							imen[ix][1]=ia;
							imen[ix][2]=ib;
							imen[ix][3]=ic;
							
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelmがielmに接する面番号

							if(ix>=40000)cout<<"ix>40000"<<endl;
							//if(ix>=memory)cout<<"ix>memory"<<endl;

							if(kmen[ix]==-1) 
							{
								cout<<"error ielm/jelm="<<ielm<<"/"<<jelm<<"  ELEM[jelm].elm="<<ELEM[jelm].elm[1]<<","<<ELEM[jelm].elm[2]<<","<<ELEM[jelm].elm[3]<<","<<ELEM[jelm].elm[4]<<endl;
								//なんでifaceエラーになるかよくわからない
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
								return OFF;//新節点導入をあきらめる
							}
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
							//if(ELEM[jelm].material!=AIR || depth[jelm]==1)//物体は破壊してほしくない+深さ1もいや。ただしこれするとたまにifaceがエラーになる。なんで？
							
							//深さ1の空気もきってよいときは下を消す(nanoeはＯＮ ferrofluidはＯＦＦ）
							//if(CON->get_model_number()!=15)
							{
								//if(ELEM[jelm].material==COIL)//物体は破壊してほしくないから
								//if(ELEM[jelm].material==COIL || ELEM[jelm].material==FLUID)
								//{	
								//	for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
								//	return OFF;//新節点導入をあきらめる
								//}
							}//////*/

							/*if(ELEM[jelm].material==MAGNET)//物体は破壊してほしくないから
							{	//なんでifaceエラーになるかよくわからない
								for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
								return OFF;//新節点導入をあきらめる
							}//////*/
							//if(vol[ix]<=ERR)//教科書図3.10のように、不適切な四面体を想定してしまい、結果体積が負になる。この場合、教科書に書いてあるとおり、最初からやりなおす(goto 10)
							if(vol[ix]<=1.0e-18)
							{
								//if(vol[ix]<0) cout<<"polyFINE vol[ix]<0 ielm="<<ielm<<endl;
								if(CON->get_CDT_sw()==OFF)
								{
									if(ELEM[jelm].material==FLUID) cout<<"流体要素がFINEの対象となりました(vol[ix]<ERR)"<<endl;
									*iv=*iv+1;//リストを増やす *iv++はだめ？		
									kv[*iv]=jelm;//本当は外接球内に新節点を含まないが、多面体分割を可能にするためにリストにいれる
									ELEM[jelm].map=1;
									flag=0;
								}
								else if(CON->get_CDT_sw()==ON)
								{
									if(ELEM[jelm].material==FLUID || ELEM[jelm].remesh==OFF)
									{
										for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
										return OFF;//新節点導入をあきらめる
									}
									else
									{
										*iv=*iv+1;//リストを増やす *iv++はだめ？		
										kv[*iv]=jelm;//本当は外接球内に新節点を含まないが、多面体分割を可能にするためにリストにいれる
										//if(ELEM[jelm].material==FLUID) cout<<"流体要素がFINEで巻き込まれました"<<endl;
										//if(ELEM[jelm].remesh==OFF) cout<<"静的要素がFINEで巻き込まれました"<<endl;
										ELEM[jelm].map=1;
										flag=0;
									}
								}
							}
						}
					}
				}
			}
        }
    }
    //////多面体表面数ixとそれを構成する頂点がもとまった
   
    /////////////表面の頂点と新節点をつなげる.すると多面体からix個の要素が生成される

    int ibound=ix;//ixの代わり。つまり表面の数を表す
    
	///体積が0になるならやめる
	for(int i=1;i<=ibound;i++)
	{
		if(vol[i]<=0)
		{
			if(vol[i]==0) cout<<"vol=0"<<" ip="<<ip<<" i="<<i<<endl;
			else if(vol[i]<0) cout<<"vol<=0"<<endl;
			//なんで0になるかよくわからない
			for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;//ここで初期化しておかないと次の節点のときバグる
			return OFF;//新節点導入をあきらめる
		}
	}////

	int nelm0=*nelm;//変更前の要素数を記憶

    for(int i=*(iv)+1;i<=ibound;i++) ///iv個の要素が破壊されて新たにibound個の要素が生成されるから、増える要素数はibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++という書き方ではだめ？
        kv[i]=*nelm;		//新節点を含む要素リストも増える。
		ELEM[*nelm].map=1;	//生成される要素は必ず新節点を外接球に含む(としている).(この行意味なくない？)
		
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//マッピングの初期化

	/*///新しいibound個の要素の深さをもとめる
	int aveDEP=0;//平均深さ
	for(int i=1;i<=ibound;i++)
	{
		aveDEP+=depth[kv[i]];
	}
	aveDEP/=ibound;		*/							////平均深さがもとまった

    for(int i=1;i<=ibound;i++)//要素情報生成
    {   
        int ielm=kv[i];
	//	double determ=vol[i];

		/*
		int ia=imen[1][i];
		int ib=imen[2][i];
		int ic=imen[3][i];
		*/
		
		int ia=imen[i][1];
		int ib=imen[i][2];
		int ic=imen[i][3];
		
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;//新点は４番目と定義	
		ELEM[ielm].elm[4]=jmen[i];
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		ELEM[ielm].volume=vol[i];
		ELEM[ielm].remesh=ON;
		
	
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//外接球の中心（ボロノイ点)と半径の二乗を計算
		
		ELEM[ielm].material=AIR;//通常のpoly関数との相違点。ここで材質を空気と決定する

		//depth[ielm]=aveDEP;//仮の深さとして、平均深さを代入
    }
    ///////////////////
    
    //要素-要素関係修正/////////上の処理で第4面で接する要素番号はわかっているので、残りを求める
	//						ここで、1〜3面は多面体を構成する要素との境界面であることに注意
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]はすでにもとまったから、それ以外をもとめる
		{
			///ELEM[ielm].node[4]=ipである
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3のとき、2,3,1の順
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3のとき、3,1,2の順
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					/*
					int ja=imen[1][k];
					int jb=imen[2][k];
					*/
					
					int ja=imen[k][1];
					int jb=imen[k][2];
					
					if(ia==ja && ib==jb)//節点が一致したら
					{
						ELEM[ielm].elm[j]=jmen[k];//あらかじめﾘｽﾄしてあった情報を格納
						ELEM[jmen[k]].elm[kmen[k]]=ielm;

						/*
						imen[1][k]=imen[1][ix];		//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
						imen[2][k]=imen[2][ix];
						*/
						
						imen[k][1]=imen[ix][1];		//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
						imen[k][2]=imen[ix][2];
						
						jmen[k]=jmen[ix];
						kmen[k]=kmen[ix];
						ix--;						//待ち辺数減少
						flag=1;						//ELEM[ielm].elm[j]はもとまったので、下のネストに入る必要はないのでflag=1
					}
				}
			}
			if(flag==0)
			{
				ix++;			//ここでのixは、[隣接関係を満たす要素]をまっている[辺]の数を表す。

				/*
				imen[1][ix]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
				imen[2][ix]=ia;
				*/
				
				imen[ix][1]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
				imen[ix][2]=ia;
				
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///要素-要素関係修正完了

	/*/////////////////////////////////////////////////新要素の深さ決定(深さ変化が連続になるように、収束するまで計算を繰り返す)
	int count2=10;	//深さが修正された要素数
	while(count2>0)
	{
		count2=0;	//初期化
		for(int i=1;i<=ibound;i++)
		{
		    int ielm=kv[i];
			int flag=0;		//これが1なら物体節点をもっているということ。その場合は深さを1に設定する
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[ielm].node[j];
				if(NODE[ia].material!=AIR)
				{
					flag=1;			//境界要素であるとわかる
					if(depth[ielm]!=1) count2++;
					depth[ielm]=1;	//境界要素は深さ1に設定
				}
			}
			if(flag==0)						//境界要素でないのなら
			{
				int smallDEP=depth[ielm];	//隣接する要素のなかでの最少深さをもとめる
				for(int j=1;j<=4;j++)		//ELEM[ielm].elm[4]は多面体の外の要素である
				{
					int jelm=ELEM[ielm].elm[j];
					if(jelm!=0) if(depth[jelm]<smallDEP) smallDEP=depth[jelm];
				}
				///要素の深さを、周囲の最少深さ+1と定義。物体に隣接する場合は下の式より深さは1となる
				if(depth[ielm]!=smallDEP+1) count2++;
				depth[ielm]=smallDEP+1;//多面体内要素のdepthをいっぺんにもとめて大丈夫か？(物体に隣接する場合は式より深さは1となるから)
			}
		}
	}

	//多面体のひとつ外側の要素の深さを修正
	for(int i=1;i<=ibound;i++)
	{
		int ielm=kv[i];
		int jelm=ELEM[ielm].elm[4];//ELEM[ielm].elm[4]は多面体の外の要素である
		///jelmがielmより深い場合、深さをielmより1だけ大きくする。もともとそうだったなら下の式を計算しても変化しない.
		//また、jelmのほうが浅かったり(その場合は-1のはず)、深さがielmと同じである場合はそのままでよい
		if(jelm!=0) if(depth[jelm]>depth[ielm]) depth[jelm]=depth[ielm]+1;
	}
	////////////////////深さ修正完了*/

	//ﾗﾌﾟﾗｼｱﾝ法実行
	//LAPLAS1(NODE,ELEM,ip,ip,kv,ibound);


    /////////新たに作られた四面体の個数が、古いものより少なくなった場合(多面体を構成したとき、どの面も多面体の境界ではない内部要素が存在するとき)
    if(*iv>ibound)
    {   
        int ir=*(iv)-ibound;		//ir個の要素が削除されなくてはならない.ただし普通に削除したのでは、要素配列に[穴]が生じることになる
		
		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];	//mapとしてとりあえずその要素番号を格納
		}
		///上で代入したmapの値を小さい順にならびかえる。その都合上、上ではkv[1],kv[2]・・・と値をつめている
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//iが増えるにつれ、ir-i+1はirから1ずつ小さくなっていく
			ELEM[ielm].map=0;	//上のほうでmapの初期化を行っているが、iv>iboundの場合、ir個の要素は初期化されていないので、ここで初期化
	    
			if(ielm!=*nelm)		//ielm番目の要素に、nelm番目の要素情報を上書き
			{   
				ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				ELEM[ielm].material=ELEM[*nelm].material;	//普通のpolyとの相違点。材質も引き継ぐ
				//depth[ielm]=depth[*nelm];					//普通のpolyとの相違点。深さも引き継ぐ
				ELEM[ielm].remesh=ELEM[*nelm].remesh;
				for(int j=1;j<=4;j++)
				{   
					ELEM[ielm].node[j]=ELEM[*nelm].node[j];
					int jelm=ELEM[*nelm].elm[j];
					ELEM[ielm].elm[j]=jelm;
					if(jelm!=0)
					{   
						int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
					}
				}
			}
			*nelm=*nelm-1;
		}
    } ///////////////////*/

	//for(int D=0;D<4;D++) delete [] imen[D];
	//delete [] jmen;
	//delete [] kmen;
	//delete [] vol;
	return ON;//polyが成功したというしるしを返す
}		


//空気層生成関数
int make_air_layer(vector<point3D> &NODE,vector<element3D> &ELEM,int *nelm,mpsconfig *CON,int *node_num,int KTE,double rrm,int KTJ)
{
	//空気層を作成する関数 return として節点増加数を返す
	cout<<"空気層生成--";
	int node=*node_num;		//現在の節点数
	int nelm0=*nelm;		//現在の要素数
	int limit_num=CON->get_add_points();		//追加できる最大節点数
	ELEM[0].material=AIR;	//ELEM[0]の材質は定まっていないから、アクセスしたらエラーになるかも。なのでここで適当にいれておく。

	int *nei_num=new int[node+1];			//近隣表面要素数
	
	double *direct[DIMENTION];				//法線ベクトル作成
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [node+1];//外向き法線ベクトル

	for(int i=0;i<=node;i++)
	{
		nei_num[i]=0;
		for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	}
	
	//流体節点の最大高さ、最小高さ
	double Zmax=0;
	double Zmin=10;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Z]>Zmax) Zmax=NODE[i].r[A_Z];
			if(NODE[i].r[A_Z]<Zmin) Zmin=NODE[i].r[A_Z];
		}
	}

	for(int je=1;je<=nelm0;je++)
	{
		if(ELEM[je].material==FLUID)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[je].elm[j];
				if(jelm!=0 && ELEM[jelm].material==AIR)//空気と接した面なら
				{
					int ia=ELEM[je].node[j%4+1];//jeとjelmの接する三角形を構成する節点番号  ia,ib,ic,ipを頂点として新しい要素が作られる
					int ib=ELEM[je].node[4-(j-1)/2*2];
					int ic=ELEM[je].node[3-(j/2%2)*2];

					double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
					double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					}
					///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
					double S[3];//上記のﾍﾞｸﾄﾙ成分格納
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);

					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					
					////面積Sがもとまった

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//外向き単位法線ﾍﾞｸﾄﾙ
					
					for(int D=0;D<3;D++)
					{
						direct[D][ia]+=n[D];
						direct[D][ib]+=n[D];
						direct[D][ic]+=n[D];
					}
					nei_num[ia]++;
					nei_num[ib]++;
					nei_num[ic]++;
				}
			}
		}
	}

	//法線ベクトルを正規化
	for(int i=1;i<=node;i++)
	{
		if(nei_num[i]>0)//表面節点なら
		{
			double val=0;				//法線ベクトルの大きさ
			for(int D=0;D<DIMENTION;D++) val+=direct[D][i]*direct[D][i];
			val=sqrt(val);
			for(int D=0;D<DIMENTION;D++) direct[D][i]/=val;		//正規化
		}
	}
	//単位法線ベクトルが求まった
	

	//単位法線ベクトル上の新節点を用いて再分割実行
	
	int newnode=0;							//新しく導入した節点数
	int *kv=new int[KTE];					//新節点を外接球にふくむ要素群
    int *istack=new int[KTE];				//一時配列

	int MESH;								//locate3Dで最初に探索する要素番号
	int add_num=CON->get_layer_num();							//空気層を何層生成するか
	double dL=CON->get_distancebp()*CON->get_layer_thin();	//空気層の幅
		
	///初期化
	for(int i=1;i<KTE;i++) ELEM[i].map=0;
	
	///順次節点を導入していく
	    
	MESH=*nelm;

	for(int i=1;i<=node;i++)
	{
		if(nei_num[i]>0) //表面節点なら
		{
			for(int k=1;k<=add_num;k++)
			{
				if(newnode<=limit_num)
				{
					double r[3];//新節点の座標格納
					for(int D=0;D<3;D++) r[D]=NODE[i].r[D]+direct[D][i]*dL*k;
					
					/*
					if(i<=node) for(int D=0;D<3;D++) r[D]=NODE[i].r[D]+direct[D][i]*dL*k;
					else
					{
						r[A_X]=0;
						r[A_Y]=0;
						if(i==node+1) r[A_Z]=Zmin-dL*k;
						else if(i==node+2) r[A_Z]=Zmax+dL*k;
					}
					/*/

					double xp=r[A_X];
					double yp=r[A_Y];
					double zp=r[A_Z];
					int flag3=OFF;//ONなら新節点を導入するかどうかが決まったということ
					int loc=locate3D(NODE,ELEM,MESH,r[A_X],r[A_Y],r[A_Z]);//新点を含む要素を探索 要素番号[MESH]から探索する
					if(loc>0)	//remesh領域のみデローニ分割時には、新点が静的要素領域に侵入する場合が考えられる。そのときはloc=0となってしまう。
					{
						if(ELEM[loc].material==AIR) flag3=ON;	//新点を含む要素が空気要素ならFINE実行
						MESH=loc;								//loc=0のときにMESH=locすると、次からlocate3Dできないので、この文をこのカッコの外に出さないよう注意
					}
					if(flag3==ON )
					{
						newnode++;		//新点数増加
						*node_num=*node_num+1;
						int ip=*node_num;
						for(int D=0;D<3;D++) NODE[ip].r[D]=r[D];
						NODE[ip].boundary_condition=0;
						NODE[ip].material=AIR;
						NODE[ip].remesh=ON;			//十中八九、remesh領域に存在することになる。
						NODE[ip].particleID=-1;		//対応する粒子は存在しな
								
						//////////外接球内に新節点を含む要素の抽出
						int iv=0;
						int msk=0;
				
						iv++;//外接球内に新節点を含む要素数
						kv[iv]=loc;
						ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むということ
						msk++;
						istack[msk]=loc;
						
						while(msk!=0)
						{   
							int isk=istack[msk];//いま注目している要素の番号
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//iskと接する要素
								if(jelm!=0)//それが表面でないなら
								{
									if(ELEM[jelm].map==0) //まだ検査してないなら
									{   
										double rad=ELEM[jelm].RR*(1.000000+ERR);//外接球半径の２乗
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
										if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//外接球内に含む
												{
													if(ELEM[jelm].material!=FLUID)
													{
														iv++;//外接球内に新節点を含む要素数を＋１
														kv[iv]=jelm;//リストにいれる
														msk++;
														istack[msk]=jelm;
														ELEM[jelm].map=1;//jelmは外接球内に新節点を含む
													}
												}
											}
										}
									}
								}
							}
						}//新点を外接球内に含む要素数ivとその要素番号kv[iv]がもとまった

						/////////////////
						int NN=*nelm;//分割前の要素数
						
						////得られた多面体を四面体に分割し、新しい要素を生成する。また、関数内で新要素の材質を決定する
						int FLAG=poly3D_for_FINE3D(NODE,ELEM,&iv,kv,ip,nelm,CON); 
						
						if(FLAG==OFF)		//polyが失敗した場合
						{
							*node_num=*node_num-1;	//新節点導入をあきらめる
							newnode--;
						}
					}
				}
			}
		}
		
	}
	
	delete [] nei_num;
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	delete [] kv;			
    delete [] istack;
	cout<<"完了  節点数+="<<newnode<<endl;
	return newnode;
}


//静的要素保存関数
void memorize_static_NODE_and_ELEM(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<point3D> &static_NODE,vector<element3D> &static_ELEM,int node,int nelm)
{
	//NODE,ELEMのうち、動かない要素、節点だけをstatic_NODE,staticELEMに格納する
	cout<<"静的節点・要素情報作成--";
	int STATIC=1;
	int DYNAMIC=2;
	int snode=0;		//静的節点数
	int snelm=0;		//静的要素数

	int *node_map=new int[node+1];					//各節点が静的かどうかを判断する
	int *newID=new int[node+1];						//各節点がstaticとして何番目の節点に変更になったか.DYNAMICならダミーとして-1を格納
	int *new_elemID=new int[nelm+1];				//各要素がstaticとして何番目の要素に変更になったか.DYNAMICならダミーとして-1を格納

	for(int i=1;i<=nelm;i++)
	{
		new_elemID[i]=-1;
		ELEM[i].map=DYNAMIC;	//初期化 
	}
	for(int i=0;i<=node;i++)
	{
		node_map[i]=DYNAMIC;	//初期化
		newID[i]=-1;
	}

	//各要素が静的か動的かを判断する。結果は.mapに格納
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int ip=ELEM[i].node[j];
			if(NODE[ip].remesh==OFF) ELEM[i].map=STATIC;//non-remesh節点をひとつでも含んでいれば、それは静的要素
		}
	}////////////

	//各節点の静的・動的を判断
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==STATIC)
		{
			for(int j=1;j<=4;j++) node_map[ELEM[i].node[j]]=STATIC;//静的要素を構成する節点も静的と判断
		}
	}

	//static_NODEに出力
	point3D NODE0;
	static_NODE.push_back(NODE0);		//節点番号は1からなので、ここで一回、[0]の配列を確保だけしておく
	for(int i=1;i<=node;i++)
	{
		if(node_map[i]==STATIC)			//静的な節点情報のみstatic_NODEに格納
		{
			snode++;
			static_NODE.push_back(NODE0);
			for(int D=0;D<3;D++) static_NODE[snode].r[D]=NODE[i].r[D];
			static_NODE[snode].boundary_condition=NODE[i].boundary_condition;
			static_NODE[snode].material=NODE[i].material;
			static_NODE[snode].particleID=NODE[i].particleID;
			static_NODE[snode].remesh=NODE[i].remesh;
			static_NODE[snode].BD_node=NODE[i].BD_node;

			newID[i]=snode;							//節点iは節点snodeになった
		}
	}/////////

	//////static_ELEM作成
	element3D ELEM0;
	static_ELEM.push_back(ELEM0);
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].map==STATIC)
		{
			snelm++;
			static_ELEM.push_back(ELEM0);
			for(int j=1;j<=4;j++)
			{
				static_ELEM[snelm].node[j]=newID[ELEM[i].node[j]];//節点番号が新しくなっていることに注意
				static_ELEM[snelm].elm[j]=ELEM[i].elm[j];			//要素番号は古いまま。あとで修正する
			}
			static_ELEM[snelm].volume=ELEM[i].volume;
			static_ELEM[snelm].material=ELEM[i].material;
			static_ELEM[snelm].RR=ELEM[i].RR;
			static_ELEM[snelm].map=ELEM[i].map;
			for(int D=0;D<3;D++) static_ELEM[snelm].r[D]=ELEM[i].r[D];
			//辺情報はどうする？ 辺もnewID作らないといけない
			new_elemID[i]=snelm;								//i番目の要素はsnelm番目になった
		}
	}

	//要素ー要素情報作成
	for(int i=1;i<=snelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=static_ELEM[i].elm[j];//格納されている要素番号。古いので新しいのにかえる
			if(jelm!=0)						//ゼロのときは書き換える必要なし
			{
				if(new_elemID[jelm]!=-1) static_ELEM[i].elm[j]=new_elemID[jelm];
				else
				{
					static_ELEM[i].elm[j]=0;//動的要素と接している場合、とりあえずゼロを格納
					//cout<<i<<endl;
				}
			}
		}
	}

	delete [] node_map;
	delete [] newID;
	delete [] new_elemID;
	cout<<"ok"<<endl;
}

////節点隣接要素数計算関数
void set_jnb3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb)
{
    for(int i=1;i<=node;i++) jnb[i]=0;//初期化
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++) jnb[ELEM[je].node[j]]=jnb[ELEM[je].node[j]]+1;
    }
}

////節点隣接要素番号格納関数
void set_nei3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
    int *num=new int [node+1];///数えあげ変数
    for(int i=1;i<=node;i++)
    {
        num[i]=0;//初期化
        for(int j=1;j<=jnb[i];j++) nei[i][j]=0;//初期化
    }
    
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++)
		{
			num[ELEM[je].node[j]]=num[ELEM[je].node[j]]+1;///set_jnb3Dでやってるし、二度手間か・・・
			nei[ELEM[je].node[j]][num[ELEM[je].node[j]]]=je;
        }
    }
    /*//check
    for(int i=1;i<=node;i++)
    {
        for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];
			int flag=0;
			for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==i) flag=1;
			if(flag==0) cout<<"EE"<<endl;
		}
    }/////////*/
    
    delete [] num;
}  




int poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON)
{   
    ////この関数の前の段階で、新節点を外接球に含む四面体数ivとその要素番号kv[iv]がもとまっている。
	////また、新節点を外接球に含む四面体はELEM[i].map=1となっている
    
	int memory=CON->get_poly_memory();		//動的に確保するメモリ数
	
    int ix=0;				//表面の数 一般にiv=ixとはならない.ひとつの要素が複数の表面を担うので。(iv<=ixである)
	int *imen[3+1];
	for(int i=1;i<=3;i++) imen[i]=new int [memory]; //多面体表面三角形の節点番号格納
    int *jmen=new int [memory];		//多面体表面三角形に隣接する四面体番号格納
    int *kmen=new int [memory];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
    double *vol=new double [memory];		//多面体の体積の６倍
    
    ///////多面体表面数ixとそれを構成する頂点をもとめる
    int flag=0;
    while(flag==0)
    {   
        ix=0;
        flag=1;
        for(int i=1;i<=*iv;i++)//*ivは新点を外接球内に含む要素の数
        {
			if(flag==1)
			{
        		int ielm=kv[i];//新点を外接球内に含む要素番号
				for(int j=1;j<=4;j++)
				{
					if(flag==1)
					{
						int jelm=ELEM[ielm].elm[j];//ielm要素に隣接する要素番号
						int ia=ELEM[ielm].node[j%4+1];//ielmとjelmの接する三角形を構成する節点番号  ia,ib,ic,ipを頂点として新しい要素が作られる
						int ib=ELEM[ielm].node[4-(j-1)/2*2];
						int ic=ELEM[ielm].node[3-(j/2%2)*2];
	                
						if(ix>49998) cout<<"pinch"<<endl;
						if(jelm==0)//表面ならielmは表面要素であるとわかる
						{
							ix++;
							imen[1][ix]=ia;	//多面体表面三角形の頂点番号
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=0;		//多面体表面三角形に隣接する要素番号
							kmen[ix]=0;		//多面体表面三角形に隣接する要素の隣接面番号  相手は0だから、ここでは0を代入
		
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);
						}
						else if(ELEM[jelm].map==0)//ielmに接するjelmは新点を外接球に含まない。つまりielmとjelmの境界は多面体表面ということ
						{
							ix++;
						
							imen[1][ix]=ia;
							imen[2][ix]=ib;
							imen[3][ix]=ic;
							jmen[ix]=jelm;
							kmen[ix]=iface3D(ELEM,ielm,jelm);//jelmがielmに接する面番号
							
							if(ix>=memory)cout<<" ix>memory"<<endl;
							vol[ix]=volume3D(NODE,ia,ib,ic,ip);

							//if(vol[ix]<=ERR)//教科書図3.10のように、不適切な四面体を想定してしまい、結果体積が負になる。この場合、教科書に書いてあるとおり、最初からやりなおす(goto 10)
							if(vol[ix]<=1.0e-18)
							//if(vol[ix]<=0)
							{ 
								//if(vol[ix]<0) cout<<"poly3D vol[ix]<0 ielm="<<ielm<<endl;
								//if(vol[ix]==0) cout<<"vol[ix]=0"<<endl;
								if(CON->get_defer_f()==OFF)
								{
									*iv=(*iv)+1;//リストを増やす *iv++はだめ？
									if(ELEM[jelm].remesh==OFF) cout<<"静的要素が多面体形成時に巻き込まれました(vol<ERR)"<<endl;
									kv[*iv]=jelm;//本当は外接球内に新節点を含まないが、多面体分割を可能にするためにリストにいれる
									ELEM[jelm].map=1;
									flag=0;
								}
								if(CON->get_defer_f()==ON)
								{
									/*
									if(ELEM[jelm].remesh==OFF)
									{
										for(int ii=1;ii<=*iv;ii++) ELEM[kv[ii]].map=0;
										return 0;
									}
									*/
									//else
									{
										*iv=(*iv)+1;//リストを増やす *iv++はだめ？
										if(ELEM[jelm].remesh==OFF) cout<<"静的要素が多面体形成時に巻き込まれました(vol<ERR)"<<endl;
										kv[*iv]=jelm;//本当は外接球内に新節点を含まないが、多面体分割を可能にするためにリストにいれる
										ELEM[jelm].map=1;
										flag=0;
									}
								}
							}
						}
					}
				}
			}
        }
    }
    //////多面体表面数ixとそれを構成する頂点がもとまった
   
    /////////////表面の頂点と新節点をつなげる.すると多面体からix個の要素が生成される

    int ibound=ix;//ixの代わり。つまり表面の数を表す
    
    for(int i=*(iv)+1;i<=ibound;i++) ///iv個の要素が破壊されて新たにibound個の要素が生成されるから、増える要素数はibound-iv
    {   
        *nelm=*nelm+1;		//*nelm++という書き方ではだめ？
        kv[i]=*nelm;		//新節点を外接球に含む要素リストも増える。
		ELEM[*nelm].map=1;	//生成される要素は必ず新節点を外接球に含む(としている).(この行意味なくない？)
    }
    
    for(int i=1;i<=ibound;i++) ELEM[kv[i]].map=0;//マッピングの初期化
    
    for(int i=1;i<=ibound;i++)//要素情報生成
    {   
        int ielm=kv[i];
		double determ=vol[i];
		
		int ia=imen[1][i];
		int ib=imen[2][i];
		int ic=imen[3][i];
		ELEM[ielm].node[1]=ia;
		ELEM[ielm].node[2]=ib;
		ELEM[ielm].node[3]=ic;
		ELEM[ielm].node[4]=ip;		//新点は４番目と定義	
		ELEM[ielm].elm[4]=jmen[i];	
		if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
		///これでいい？
		ELEM[ielm].volume=vol[i];
		ELEM[ielm].material=AIR;
		ELEM[ielm].remesh=ON;
		sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//外接球の中心（ボロノイ点)と半径の二乗を計算
    }
    ///////////////////
    
    //要素―要素関係修正/////上の処理で第4面で接する要素番号はわかっているので、残りを求める
	//						ここで、1〜3面は多面体を構成する要素との境界面であることに注意
    ix=0;
    for(int i=1;i<=ibound;i++)
    {
        int ielm=kv[i];
		for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]はすでにもとまったから、それ以外をもとめる
		{
			///ELEM[ielm].node[4]=ipである
			ELEM[ielm].elm[j]=-1;				//初期化
			int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3のとき、2,3,1の順
			int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3のとき、3,1,2の順
			int flag=0;
			for(int k=1;k<=ix;k++)
			{
				if(flag==0)
				{
					int ja=imen[1][k];
					int jb=imen[2][k];
					if(ia==ja && ib==jb)//節点が一致したら
					{
					    ELEM[ielm].elm[j]=jmen[k];//あらかじめﾘｽﾄしてあった情報を格納
					    ELEM[jmen[k]].elm[kmen[k]]=ielm;
					    imen[1][k]=imen[1][ix];//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
					    imen[2][k]=imen[2][ix];
					    jmen[k]=jmen[ix];
					    kmen[k]=kmen[ix];
					    ix--;	//待ち辺数減少
						flag=1;	//ELEM[ielm].elm[j]はもとまったので、下のネストに入る必要はないのでflag=1
					}
				}
			}
			if(flag==0)
			{
			    ix++;			//ここでのixは、[隣接関係を満たす要素]をまっている[辺]の数を表す。
				imen[1][ix]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
				imen[2][ix]=ia;
				jmen[ix]=ielm;
				kmen[ix]=j;
			}
		}
    }///要素-要素関係修正完了

    /////////新たに作られた四面体の個数が、古いものより少なくなった場合(多面体を構成したとき、どの面も多面体の境界ではない内部要素が存在するとき)
    if(*iv>ibound)
    {   
	
        int ir=*(iv)-ibound;	//ir個の要素が削除されなくてはならない.ただし普通に削除したのでは、要素配列に[穴]が生じることになる

		for(int i=1;i<=ir;i++)
		{
			kv[i]=kv[ibound+i];
			ELEM[kv[i]].map=kv[i];//mapとしてとりあえずその要素番号を格納
		}
		///上で代入したmapの値(要素番号)を小さい順にならびかえる。その都合上、上ではkv[1],kv[2]・・・と値をつめている
		qsorti3D(NODE,ELEM,ir,kv);
	
		for(int i=1;i<=ir;i++)
		{   
			int ielm=kv[ir-i+1];//iが増えるにつれ、ir-i+1はirから1ずつ小さくなっていく
			ELEM[ielm].map=0;	//上のほうでmapの初期化を行っているが、iv>iboundの場合、ir個の要素は初期化されていないので、ここで初期化
			
			if(ielm!=*nelm)		//ielm番目の要素に、nelm番目の要素情報を上書き
			{   
			    ELEM[ielm].r[A_X]=ELEM[*nelm].r[A_X];
				ELEM[ielm].r[A_Y]=ELEM[*nelm].r[A_Y];
				ELEM[ielm].r[A_Z]=ELEM[*nelm].r[A_Z];
				ELEM[ielm].RR=ELEM[*nelm].RR;
				ELEM[ielm].remesh=ELEM[*nelm].remesh;
				for(int j=1;j<=4;j++)
				{   
				    ELEM[ielm].node[j]=ELEM[*nelm].node[j];
				    int jelm=ELEM[*nelm].elm[j];
				    ELEM[ielm].elm[j]=jelm;
				    if(jelm!=0)
				    {   
				        int N=iface3D(ELEM,*nelm,jelm);
						ELEM[jelm].elm[N]=ielm;
						//材質は？？？
				    }
				}
			}
			*nelm=*nelm-1;
		}
    }///////////////////*/

	delete [] jmen;
	delete [] kmen;
	delete [] vol;

	for(int i=1;i<=3;i++) delete [] imen[i];

	return ON;
	
}		