// MD_fo_mol.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"


/*改変履歴
2008/1/24	鈴木プログラムのRKGを使って、とりあえず完成。
			変数名を整理して、構造化を進める。パラメータファイルの書式を変更。
			分子回転やTOF、PIPICO作成機能は削除。（おそらく他のプログラムに分けた方が良い。）
			相互作用にCoulomb、LJ、IDを含む。
2008/1/25	snapshot出力機能を付加。
			変数の型を整理。（無駄なdouble型をintに変更）
			エネルギー保存と運動量保存の確認ルーチンは再チェックすること！
2012/10/    多価イオン生成を許す。
			入力、出力ファイルフォーマット変更。全ての結果をinformation.txtに入れることにする。
2015/1/28   charge build-up 時間を導入。これに合わせて、、入力ファイルにパラメータ追加	
2015/7/10   charge transfer 時間を導入。これに合わせて、、入力ファイルにパラメータ追加	
*/ //改変履歴終わり


//CoulombExplosionMomentum.cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <direct.h>
using namespace std;

//#define CONST 138.6 
//時間の単位は ps 長さの単位は nm エネルギーの単位は amu nm^2/ps^2
//CONST=6.02*10^23*e^2/(4pi*e0)        鈴木定数

#define MAX_NUM 200 
//最大粒子数

#define MAX_LABEL 20 
//最大label数

#define PI   3.14159282      //Pi
#define E_0  8.854187817e-12 //epsilon_0 [F/m]
#define AMU  1.6605402e-27   //atomic mass unit [kg]
#define K_B  1.380658e-23    //Boltzmann constant[J/K]
#define UNIT_E 1.602177e-19   //elemental charge[C]
#define EPSILON 166.4    //epsilon/k_B for Kr [K]
#define SIGMA 0.365      //LJ parameter sigma [nm] for Kr 
//#define EPSILON 120.0      //LJ parameter epsilon/k_B for Ar
//#define SIGMA 0.341      //LJ parameter sigma [nm] for Ar 
//#define EPSILON 36.7     //LJ parameter epsilon/k_B for Ne  
//#define SIGMA 0.279      //LJ parameter sigma [nm] for Ne 
#define  ALPHA_CGS  2.46e-30       //polarizabiliy [m^3] for Kr, CGS unit,  4*pi*e_0*alpha_CGS[cm^3] =10^6 * alpha_SI[C m^2 / V] 
//#define  ALPHA_CGS  1.62e-30       //polarizabiliy [m^3] for Ar, CGS unit
//#define  ALPHA_CGS  0.390e-30      //polarizabiliy [m^3] for Ne, CGS unit


//電荷の個情報
class denka
{
	private:
	int PN;		                        // number of simulated atom (ion + neutral)
	int CN;			                    // number of charge
	long double position[MAX_NUM][3];   // coordinates of i-th atom    [nm]
	long double velocity[MAX_NUM][3];   // velocity of i-th atom      [nm/ps]
	long double mass[MAX_NUM];          // mass of i-th atom          [amu]
	long double momentum[MAX_NUM][3];   // momentum of i-th atom      [amu nm/ps]
  	int	charge_dest[MAX_NUM];		    // destination charge of i-th atom (distributed roundomly at the beggining)
  	long double	charge[MAX_NUM];        // charge of i-th atom
  	long double	charge_last[MAX_NUM];   // charge of i-th atom at last step
  	long double	charge_uracil;			// charge of uracil ring at t
  	long double	charge_next[MAX_NUM];   // charge of i-th atom at last step
  	int	Q_limit[MAX_NUM];		        // 最大charge of i-th atom

	long double t_charge_bu, t_CT, dt_temp; // 電荷buld-up time [ps], charge transfer time [ps]
  	long double	CT_flag[MAX_NUM];       // charge transfer flags of i-th atom
  	long double	CT_exp[MAX_NUM];        // exponent of charge transfer flags for i-th atom
  	int			CT_center,CT_para;      // label number of charge transfer center atom

	int random_charge[MAX_NUM];         // 電荷をランダムに分配する原子の番号を入れておく
	int num_random_atom;				// 電荷をランダムに分配する原子の数
	int num_fixed_atom;				    // 初期データで電荷が確定されている原子の数
	string label[MAX_NUM];              // label for identifying particles
	int num_label;				        // number of independent label 
	string valid_label[MAX_NUM];        // list of independent label 
	long double t,t0,tn;                // time 
	int n_step;			                //時間の刻み数
	long double interval_t;				//時間の刻み巾 = (tn-t0)/n_step
	long double q[MAX_NUM][3][2];       //丸めの誤差
	long double px,py,pz;				//系全体の運動量の成分。 
	long double refmomentum;			//精度評価で使用。最大の全運動量の和の成分。 
	long double minmomentum;			//精度評価で使用。最小の運動量成分。　
	long double refenergy;				// 精度評価で使用
	int snap_num;						// snapshot is took every snap_num step

	ofstream foutinformation;			//情報ファイル名　fimename[i]にlabel番号iのファイル名　
    char InfoFile[32];					//情報ファイルストリーム　
    char filename[MAX_LABEL*2][32];		//出力ファイル名　fimename[i]にlabel番号iのファイル名　
	ofstream outfile[MAX_LABEL*2];		//出力ファイルストリーム outfile[i]にlabel番号iのファイル
    char filename_2body[MAX_LABEL][MAX_LABEL][32];		//出力ファイル名　fimename[i]にlabel番号iのファイル名　
	ofstream outfile_2body[MAX_LABEL][MAX_LABEL];		//出力ファイルストリーム outfile[i]にlabel番号iのファイル
	char SnapFile[32];			        //スナップショットファイル名　
    ofstream fsnapshot;					//スナップショット用ファイルストリーム　

	char pot_flag;			//相互作用　default: C, C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole
	char pot_label[8];		//相互作用識別ラベル　C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole

public:
	int trial;              // 試行回数
	void keisan(int n_th);
	long double kyori(long double,long double,long double ,long double,long double,long double);
	int dice(int faces);
	void rkg();

	long double force(int i,int j); //粒子間の力を計算。粒子 i の dv/dt の右辺を返す   j は座標ｘｙｚ
	long double potential_energy(); //全ポテンシャルエネルギーを計算 [amu nm^2/ps^2]
	long double kinetic_energy();   //全運動エネルギーを計算 [amu nm^2/ps^2]
	long double PE_plus_KE();       //全エネルギー（PE+KE）を計算 [amu nm^2/ps^2]
	void cal_momentum();    //ランダムに電荷分配する
	void dist_charge();    //ランダムに電荷分配する

    void mk_out_dir(char *makingdirectory);
    void open_output_files(char *makingdirectory);
    void close_output_files(char *makingdirectory);
	int get_data(char *input_file);
	void snap_output(long double snap_time);
	void message();

};




void denka::message(){

cout << "usage: program.exe  parameterfile " << endl; 
cout << " " << endl; 
cout << "パラメータファイルの書式 " << endl; 
cout << "１行目　　開始時間[ps]  終了時間[ps]  MDステップ数[整数]  スナップショット回数" << endl; 
cout << "　　      数字の間はスペースで区切る。" << endl; 
cout << "　　      スナップショット回数=0ならばスナップショットを出力しない。" << endl; 
cout << " " << endl; 
cout << "２行目　　全粒子数　全電価数  相互作用指定 電荷buld-up時間[ps] 電荷移動時間[ps] 電荷移動CTのフラグ(=1で個別)" << endl; 
cout << "　　　　　全粒子数 > 全電価数。" << endl; 
cout << "　　　　　電荷移動CTのフラグ: defaultはウラシルリング全体が同じレート。=1で個別)。" << endl; 
cout << "　　　　　相互作用: default=C: Coulomb, L: LJ, M: C+LJ, N: C+LJ+Induced dipole " << endl; 
cout << " " << endl; 
cout << "３行目　　シミュレーションの試行回数" << endl; 
cout << " " << endl; 
cout << "４行目以降に初期原子データを書く " << endl; 
cout << "書式 " << endl; 
cout << "ラベル  質量[amu]  x[nm] y[nm] z[nm] vx[nm/ps] vy[nm/ps] vz[nm/ps] 初期電荷 最大電荷 電荷移動CTのフラグ" << endl; 
cout << "ラベルは５文字以内。これに従って結果を分類する。 " << endl; 
cout << "初期電荷を置いてシミュレーション。初期電荷は０以上の整数" << endl; 
cout << "最大電荷まで粒子にランダムに電荷を分配する。 " << endl; 
cout << "CTのフラグ：それぞれの原子の電荷移動時間をflag×電荷移動時間にするフラグ。 " << endl; 
exit(0);
}

void denka::mk_out_dir(char *makingdirectory){

char NewDir[32];
int i;

//実行ファイル名と同じ（拡張子除く）ディレクトリ名を作成　・・.exe 実行ファイル　パラメタファイル
    i=0;
	do {
		NewDir[i]=makingdirectory[i];
        i=i+1;
        if(i>31) {cout<<"Filename is too long. Be fewer 32 letters."<<endl;exit (1);}
	} while(makingdirectory[i]!='.');
    NewDir[i]=NULL;
	cout<<"simulation result -> "<<NewDir<<endl;
    _mkdir(NewDir); //ディレクトリ作成完了


}
//出力ディレクトリ作成

void denka::open_output_files(char *makingdirectory){
char DirName[32];
int i,j;

    i=0;
	do {
		DirName[i]=makingdirectory[i];
        i=i+1;
        if(i>31) {cout<<"Filename is too long. Be fewer 32 letters."<<endl;exit (1);}
	} while(makingdirectory[i]!='.');
    DirName[i]=NULL;


         //シミュレーション条件の出力ファイルの作成
		   strcpy_s(InfoFile,DirName);
		   strcat_s(InfoFile,"\\Infomation.txt");
			foutinformation.open(InfoFile);
				cout << InfoFile << endl;
			if(!foutinformation) {cout<<" cannot make Infomation FILE. \n";}
		  foutinformation<<"START TIME T0 "<<t0<<" ps"<<endl;
	      foutinformation<<"STOP TIME Tn "<<tn<<" ps"<<endl;
	      foutinformation<<"kizami n "<< n_step <<endl;
		  foutinformation<<"trial number "<<trial<<endl;
	      foutinformation<<"Number of simulated particles is "<<PN<<endl;
		  foutinformation<<"Number of charge is "<<CN<<endl;
		  foutinformation<<"Interaction : "<< pot_label <<endl;
		  foutinformation<<"Charge buld-up time [ps] : "<< t_charge_bu <<endl;
		  foutinformation<<"New directory name is "<<DirName<<endl;
		  foutinformation<<"N"<<DirName<<endl;
		  foutinformation<<endl;
		  foutinformation<<endl;
		  foutinformation<<endl;


			//ファイルを開く、"label"P.dat (ion) labelN.dat (neutral)  などなど

			if(MAX_LABEL<num_label){cout << "too much label! " <<endl; }

			for(i=0;i<num_label;i++){
				strcpy_s(filename[2*i]  ,DirName);			// label i-th ion
				strcpy_s(filename[2*i+1],DirName);			// label i-th neutral 
				strcat_s(filename[2*i]  ,"\\");
				strcat_s(filename[2*i+1],"\\");
				strcat_s(filename[2*i]  ,valid_label[i].c_str());
				strcat_s(filename[2*i+1],valid_label[i].c_str());
				strcat_s(filename[2*i]  ,"_P.dat");
				strcat_s(filename[2*i+1],"_N.dat");
              //ファイル作成
				outfile[2*i].open(  filename[2*i]);
//				cout << filename[2*i] << endl;
                if(!outfile[2*i]) {cout << "cannot make \n ";}
				outfile[2*i+1].open(filename[2*i+1]);
//				cout << filename[2*i+1] << endl;
				if(!outfile[2*i+1]) {cout << "cannot make \n ";}
				
			}


			for(i=0;i<num_label;i++){
			for(j=i;j<num_label;j++){
				strcpy_s(filename_2body[i][j]  ,DirName);			// label i-th ion
				strcat_s(filename_2body[i][j]  ,"\\");
				strcat_s(filename_2body[i][j]  ,"2body_");
				strcat_s(filename_2body[i][j]  ,valid_label[i].c_str());
				strcat_s(filename_2body[i][j]  ,valid_label[j].c_str());
				strcat_s(filename_2body[i][j]  ,".dat");
              //ファイル作成
				outfile_2body[i][j].open(  filename_2body[i][j]);
//				cout << filename[2*i] << endl;
                if(!outfile_2body[i][j]) {cout << "cannot make \n ";}
			}
			}

			if(snap_num>0){
				strcpy_s(SnapFile ,DirName);		
				strcat_s(SnapFile,"\\Snapshot.dat");
				fsnapshot.open(SnapFile);
				cout << SnapFile << endl;
				if(!fsnapshot) {cout<<" cannot make Infomation FILE. \n";}
			}
}

void denka::close_output_files(char *makingdirectory){
int i;

			for(i=0;i<num_label;i++){outfile[2*i].close();outfile[2*i+1].close();}
			foutinformation.close();
			if(snap_num>0){	fsnapshot.close();}

}

void denka::keisan(int n_th) {
	int i,j,k,temp_i;
	int suffix_1,suffix_2;
	long double firstP[3]={0,0,0};   //初期の全運動量計算
	long double firstenergy=0;
	long double finalenergy=0;
	double Vx1,Vy1,Vz1;
	double Vx2,Vy2,Vz2;
	double InnerProduct;
	long double ev[MAX_NUM]; //各粒子の運動エネルギー
	long double exp_chg; 
	int	file_num;

          cout<<"you entered... "<<endl;
	      cout<<"START TIME T0 "<<t0<<" ps"<<endl;
	      cout<<"STOP TIME Tn "<<tn<<" ps"<<endl;
	      cout<<"kizami n "<< n_step <<endl;
		  cout<<"trial number "<< trial <<endl;
	      cout<<"Number of simulated particles is "<<PN<<endl;
	      cout<<"Number of total charge is "<<CN<<endl;

			cal_momentum();
			cout<<"initial total momentum of system : [amu nm/ps]"<<endl<<
				"px="<<px<<" py="<<py<<" pz="<<pz<<endl;


//		  //初期の全エネルギー計算
			cout<<"first potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換
			cout<<"first kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換
			firstenergy=PE_plus_KE();  ////単位は [amu nm^2/ps^2]
			cout<<"first total energy  "<< (AMU*1e6/UNIT_E)*firstenergy<<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換

			 
            //Runge Kutta Gill で計算　n は時間の刻みの数。その数だけ for で計算

			if(t_charge_bu!=0){	
				dt_temp = interval_t / t_charge_bu; 
				charge_uracil = 0;
				for(k=0;k<PN;k++){	
					charge[k] = 0; 
					charge_last[k] = 0;  
					charge_next[k] = 0;
				} 
			}

			for(i=0;i<=n_step;i++) {
				t=t0+i*interval_t;

			if(t_charge_bu!=0){		
				exp_chg=(1- exp(-1 * t / t_charge_bu));
				if(t_CT==0){
					for(k=0;k<PN;k++){	charge[k] = charge_dest[k] * exp_chg;	} //解析式で、charge build-up のみ。 動作OK
//					for(k=0;k<PN;k++){	
//						charge_next[k] = charge[k] + (dt_temp - pow(dt_temp, 2)/2 + pow(dt_temp, 3) / 6 - pow(dt_temp, 4) / 24) * (charge_dest[k] - charge[k]) ; // Taylor 4th, OK
//					} //差分式で、charge build-up のみ。 動作OK
				}else{
					for(k=0;k<PN;k++){
						if((CT_flag[k] * t_CT) == 0){
							charge_next[k] = charge_dest[k];
						}else{
							if(CT_para == 1){
								if(k != CT_center){charge_next[k] = charge[k] + charge[CT_center] * interval_t * charge[CT_center] / (CT_flag[k] * t_CT) ;}
								if(charge_next[k] > charge_dest[k]){ charge_next[k] = charge_dest[k];} 
							}else{
								charge_uracil = charge_uracil + charge[CT_center] * interval_t / t_CT ;
								if(charge_uracil > (CN - charge_dest[CT_center])){ charge_uracil = CN - charge_dest[CT_center];} 
								for(k=0;k<PN;k++){
									if(k != CT_center){charge_next[k] = charge_uracil * charge_dest[k] / (CN - charge_dest[CT_center]); }
								}
							}
						};
					} //電荷移動に応じてそれぞれの原子で計算 
				charge_next[CT_center]= exp_chg * CN;
				for(k=0;k<PN;k++){
						if(k != CT_center){charge_next[CT_center] = charge_next[CT_center] - charge_next[k]; }
					}				
				}

			}
			
						rkg();
						if( (snap_num!=0) && (i%snap_num==0) ){snap_output(t);}

			if(t_charge_bu!=0){		for(k=0;k<PN;k++){	charge[k] = charge_next[k];} }

			//エネルギー計算
   				for(k=0;k<PN;k++)  {
						ev[k]=0.5*(AMU*1e6/UNIT_E)*mass[k]*(velocity[k][0]*velocity[k][0] + velocity[k][1]*velocity[k][1]+velocity[k][2]*velocity[k][2]);
			  	};    // エネルギー単位 [amu nm^2/ps^2] --> [eV]　へ変換
			}

			if( (snap_num!=0) && (n_step%snap_num!=0) ){snap_output(t);} //最後のstepでにsanpshotを出力する。


			cal_momentum();
			cout<<"total momentum / [amu nm/ps]"<<endl<<
				"px="<<px<<" py="<<py<<" pz="<<pz<<endl;

//		  //初期の全エネルギー計算
			cout<<"final potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換
			cout<<"final kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換
			finalenergy=PE_plus_KE();  ////単位は [amu nm^2/ps^2]
			cout<<"final total energy  "<< (AMU*1e6/UNIT_E)*finalenergy<<"  eV"<<endl;  //最後に [amu nm^2/ps^2] --> [eV]　へ変換

			if(refenergy<(fabs((finalenergy-firstenergy)/finalenergy))) {
				refenergy = (fabs((finalenergy-firstenergy)/finalenergy));}

			if(refmomentum<fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum))  //大きい方をrefmomentumに採用
			{refmomentum=fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum);}
       
           
		   cout<<" seido (first total momentum-final total momentum)/(minmomentum) "<<refmomentum<<endl;
		   cout<<" seido max(first-final/final)energy = "<<refenergy<<endl;
        

           foutinformation<<endl;
		   foutinformation<<"momentum seido "<<refmomentum<<endl;
		   foutinformation<<"energy seido   "<<refenergy<<endl;



		   //labelに従って結果をファイルに出力する
			for(i=0;i<PN;i++)  {
               if(charge_dest[i]==0) {
				   for(j=0;j<num_label;j++){if(label[i]==valid_label[j]){file_num=2*j+1;};}
//						cout << charge_dest[i] << " " << label[i] << " " << filename[file_num] << endl;
			   }
               if(charge_dest[i]>0) {
				   for(j=0;j<num_label;j++){if(label[i]==valid_label[j]){file_num=2*j;};}
//						cout << charge_dest[i] << " " << label[i] << " " << filename[file_num] << endl;
			   }
			outfile[file_num] <<n_th<< " " << charge_dest[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;
			foutinformation  <<n_th<< " "  <<label[i]<< " " << charge_dest[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;


				if(charge_dest[i]>0) {
					for(j=i+1;j<PN;j++){
						if(charge_dest[j]>0) {
							for(k=0;k<num_label;k++){
								if(label[i]==valid_label[k]){suffix_1=k;};
								if(label[j]==valid_label[k]){suffix_2=k;};
							}
							if(suffix_1>suffix_2){temp_i= suffix_1; suffix_1=suffix_2; suffix_2=temp_i; }
							Vx1= momentum[i][0]; Vy1= momentum[i][1]; Vz1= momentum[i][2]; 
							Vx2= momentum[j][0]; Vy2= momentum[j][1]; Vz2= momentum[j][2]; 
							InnerProduct = (Vx1*Vx2 + Vy1*Vy2 + Vz1*Vz2)/sqrt((Vx1*Vx1 + Vy1*Vy1 + Vz1*Vz1)*(Vx2*Vx2 + Vy2*Vy2 + Vz2*Vz2));
							outfile_2body[suffix_1][suffix_2]  <<n_th<< " "  << label[i]  << " "  << charge_dest[i]  << " "  << label[j]  << " "  << charge_dest[j]  << " "  << InnerProduct  << " "  << ev[i]  << " " << ev[j]  << " " << endl;
						}
					}
				}

			}

			

}

//運動量計算
void denka::cal_momentum(){
int i,j;

//運動量を求める
for(i=0;i<PN;i++) {for(j=0;j<=2;j++)  { momentum[i][j]=mass[i]*velocity[i][j];};}

px=0;	py=0;	pz=0;

//全運動量の和
	for(i=0;i<PN;i++) {
		px=px+momentum[i][0];        py=py+momentum[i][1];		pz=pz+momentum[i][2];
		//運動量精度評価　運動量ベクトルの絶対値で評価する 先ず最も小さいイオンの運動量ベクトルを探す
		if(i==0) {minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);
		}else  if(i) {
			if(minmomentum>kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0))   
			  {minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);}
		}
		}
}			//全運動量の和終了


void denka::snap_output(long double snap_time){
	int i;

	fsnapshot << " simlation time : " << snap_time << " [ps]" << endl;
//	fsnapshot << " t_CBU, T_CT : " << t_charge_bu << " [ps]" << t_CT << " [ps]" << endl;
	for(i=0;i<PN;i++){
		fsnapshot	<< label[i] << "     "  << mass[i]  
					<< "     "  << position[i][0] << "     "  << position[i][1] << "     "  << position[i][2]
					<< "     "  << velocity[i][0] << "     "  << velocity[i][1] << "     "  << velocity[i][2]
					<< "     "  << charge_dest[i] << "     "  << charge[i] << endl;
	};
	fsnapshot	<< endl;  

}

long double denka::kyori(long double a1,long double b1,long double a2,long double b2,long double a3,long double b3)  {
return sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3));
}



//２点間の距離を計算

void denka::rkg()
      
      {
       static long double p1=0.2928932188134524, p2=1.707106781186547;   //  p1 = 1 - sqrt(1/2),  p2 = 1 + sqrt(1/2) 
	   long double ck0[MAX_NUM][3][2];
	   long double q1[MAX_NUM][3][2];
	   long double r1[MAX_NUM][3][2];
	   long double y1[MAX_NUM][3][2];
       int i,j,k;
   
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {ck0[i][j][k]=0; q1[i][j][k]=0;  r1[i][j][k]=0;y1[i][j][k]=0;};};}

      
	    //i 粒子を表す j 座標を表す
	   // ４ステップの計算開始
	   //  step 1
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t/2;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t/2;}
	   ;};}  ;}
        
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+ck0[i][j][k]-q[i][j][k];}
	   else    {y1[i][j][k]=position[i][j]+ck0[i][j][k]-q[i][j][k];}
					  };};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

	  
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];
	   };};}

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}
       


       //  step 2
//       xx=t+interval_t*0.5;

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t;}
	   };};}


        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);}
	   else    {y1[i][j][k]=position[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);}
					  };};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

      
	    for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p1;
	   };};}
      
	   
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}
    
	  

       //  step 3
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t;}
	   };};}       

        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);}
	   else    {y1[i][j][k]=position[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);}
					  };};}


       
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

       
        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p2;
	   };};}


	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}

	  
       
       //  step 4
//     xx=t+interval_t/2;
      
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t/2;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t/2;}
	   };};}       
    

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+(ck0[i][j][k]-q[i][j][k])/3;}
	   else    {y1[i][j][k]=position[i][j]+(ck0[i][j][k]-q[i][j][k])/3;}
	   };};}      


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];};};};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];};};}
	   

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];};};};}

          
    
 	  }



//Runge Kutta Gill calculation
int denka::dice(int faces) { 	return (1 + rand()%faces);}  
// n-face dice  返り値1～n

void denka::dist_charge() {
int i,num,max_charge,fixed_charge;         

     srand( (unsigned)time(NULL) );

	max_charge=0; fixed_charge=0;
	for(i=0;i<PN;i++){max_charge=max_charge+Q_limit[i];
			fixed_charge=fixed_charge+charge_dest[i];
		}

	if(CN > max_charge){ cout << " charge number mismatch" << endl; exit(0);}

	for(i=0;i<(CN - fixed_charge);i++){
			num = dice(PN)-1;
	        while( charge_dest[num]>=Q_limit[num]){num = dice(PN)-1; }
				charge_dest[num]=charge_dest[num]+1;
		}						//  電荷(CN)をランダム分配。

}


	  

//力の計算   粒子 i の dv/dt の右辺を返す   j は座標ｘｙｚ  ここのエネルギー単位は [amu nm^2/ps^2]　力は [amu nm/ps^2]
long double denka::force(int i ,int j)
{
  
   int seisu=0;
   long double nagasa=0;
   long double uhen=0;
	     
		   switch(pot_flag){
				case 'C'	: 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						if(seisu-i) {uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);}
				   };  //一価イオン間のクーロン相互作用
					break;
				case 'L'	:
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+(4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))*(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
					   };
					};  //中性原子間のLJ相互作用　イオン間クーロンは無視
					break;
				case 'M' : 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+((UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]/nagasa + 4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))
						   		 *(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
					   };
					};  //クーロン + LJ相互作用
					break;
				case 'N'	:
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
							if(seisu-i) {
								uhen=uhen+ ( (UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]/nagasa 
										+ 4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12))
										- ( ALPHA_CGS * UNIT_E*UNIT_E/(8*PI*E_0) *1e30 )*charge[i]/pow(nagasa,4) )
										* (position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
							};
				   };  //クーロン + LJ + イオン－誘起双極子相互作用
					break;
				default	: 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);
						   };
				   };  //デフォルトはイオン間のクーロン相互作用のみ。
					break;
			}
    return uhen;
}


//全エネルギー計算
long double denka::PE_plus_KE(){

	return potential_energy() + kinetic_energy();

} //全エネルギー（PE+KE）を計算 [amu nm^2/ps^2]
//全ポテンシャルエネルギー計算
long double denka::potential_energy(){

	int i,j;
	long double nagasa=0;
	long double t_energy=0;

switch(pot_flag){
	case 'C'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa;
			};	
		};  //一価イオン間のクーロン相互作用。エネルギー単位は [amu nm^2/ps^2]
		break;
	case 'L'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+4* (EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
			};	
		};  //中性原子間のLJ相互作用　イオン間クーロンは無視。エネルギー単位は [amu nm^2/ps^2]
		break;
	case 'M'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
								 +4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
			};	
		};  ////クーロン + LJ相互作用。エネルギー単位は [amu nm^2/ps^2]
		break;
	case 'N'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
								 +4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) 
								 - (ALPHA_CGS * UNIT_E * UNIT_E * 1e30/(8*PI*E_0*AMU) )*charge[i]/pow(nagasa,4) ;
			};	
		};  ////クーロン + LJ相互作用+ イオン－誘起双極子相互作用。エネルギー単位は [amu nm^2/ps^2]
		break;
	default	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa;
			};	
		};  //イオン間のクーロン相互作用。エネルギー単位は [amu nm^2/ps^2]
		break;
}
return t_energy;
} //全粒子のpotentialエネルギーを計算
					
//全運動エネルギー計算
long double denka::kinetic_energy(){
int i;
long double t_energy=0;

for(i=0;i<PN;i++) {	t_energy=t_energy+0.5*mass[i]*(velocity[i][0]*velocity[i][0]+velocity[i][1]*velocity[i][1] + velocity[i][2]*velocity[i][2]);}
return t_energy;
} //全粒子の運動エネルギーを計算


int denka::get_data(char *input_file){

  char buf[512];
  int  i,j,k,LINE_NUM,flag;

 ifstream fin(input_file);

 if(!fin) {cout << "cannot open \n ";
          return 1;}


 LINE_NUM=0;

    while( fin.getline( buf, sizeof(buf) )){
		LINE_NUM=LINE_NUM+1;
//		cout << buf << endl;
		istringstream is( buf );
		if(LINE_NUM==1){ is >> t0 >> tn >> n_step>> snap_num ;}
		if(LINE_NUM==2){ is >> PN >> CN >> pot_flag >> t_charge_bu >> t_CT >> CT_para;}
		if(LINE_NUM==3){ is >> trial ;}
		if(LINE_NUM >3){ is >> label[LINE_NUM-4] >> mass[LINE_NUM-4] 
		                    >> position[LINE_NUM-4][0] >> position[LINE_NUM-4][1] >> position[LINE_NUM-4][2] 
							>> velocity[LINE_NUM-4][0] >> velocity[LINE_NUM-4][1] >> velocity[LINE_NUM-4][2] 
							>> charge_dest[LINE_NUM-4] >> Q_limit[LINE_NUM-4] >> CT_flag[LINE_NUM-4] ;
		}
	}
fin.close();

interval_t=(tn-t0)/n_step;	
CT_center=0;
for(k=0;k<PN;k++){if(CT_flag[k]==0){ CT_center=k;}} // 電荷移動の中心原子を捜す。入力ファイルで、CT_flag=0となる最後の原子が中心。

switch(pot_flag){
	case 'C'	: strcpy_s(pot_label, "Coulomb");	break;
	case 'L'	: strcpy_s(pot_label, "LJ");	break;
	case 'M'	: strcpy_s(pot_label, "C+LJ");	break;
	case 'N'	: strcpy_s(pot_label, "C+LJ+ID");	break;
	default		: strcpy_s(pot_label, "Coulomb");	break;
}

for(j=0;j<3;j++) {for(i=0;i<PN;i++) {momentum[i][j]=0;} ;}//初期化
for(i=0;i<PN;i++) {for(j=0;j<3;j++)  {for(k=0;k<2;k++) {q[i][j][k]=0;};};}   //丸めの誤差初期化

valid_label[0] = label[0];
num_label=1;
for(i=1;i<PN;i++){
	flag=0;
	for(j=0;j<num_label;j++){ if(label[i] == valid_label[j]){ flag=1;} ;}
	if(flag==0){valid_label[num_label] = label[i]; num_label=num_label+1 ;}
}			//  labelの数と、重複しないラベルをチェック

cout << LINE_NUM-3 << " particle's data. " << PN << " particles are valid." << endl;


num_random_atom=0;
num_fixed_atom=0;
for(i=0;i<PN;i++){
	if(charge_dest[i]==-1){random_charge[num_random_atom]=i;num_random_atom=num_random_atom+1;};
	if(charge_dest[i]!=-1){num_fixed_atom=num_fixed_atom+1;};
	}

//cout << num_random_atom << " 粒子に電荷を分配する。 " << endl;
//cout << num_fixed_atom << " 粒子は電荷が確定。 " << endl;
//cout << CN << " 電荷数。 " << endl;

//for(i=0;i<num_label;i++){ 	cout << label[i] << " " << valid_label[i] << endl;}
//cout << num_label << " label数。 " << endl;


return 0;
}


//argv[1] がデータファイル 初期時刻　最終時刻　時間刻み数　粒子数　質量　電荷　初期位置xyz　初速度xyz の順　
//argy[2] がパラメータファイル　イオン数　確定原子番号　0 ランダム原子番号 0

//int _tmain(int argc, _TCHAR* argv[])
//{
//	return 0;
//}
//
//
int main(int argc,char *argv[])
{
int i;


cout<<"********************************PROGRAM START********************************"<<endl;
cout<<"length:nm,time:ps,energy:eV"<<endl;    

  denka sim_data;				//class 生成
  if(argv[1]==NULL){sim_data.message();}

  sim_data.mk_out_dir(argv[1]); //出力ディレクトリ作成

	sim_data.get_data(argv[1]);	//初期データ取得 ファイル作成のみに使用する
	sim_data.open_output_files(argv[1]);//出力ファイル作成

for(i=0;i<sim_data.trial;i++) {
	sim_data.get_data(argv[1]);	//初期データ取得,電荷分配し直し

	cout<<"******"<<i+1<<"th simulation"<< "****************************************************"<<endl;
	sim_data.dist_charge();       //電荷をランダムに分配。
	sim_data.keisan(i+1);       //ｎstepのシミュレーションをRunge-Kutta-Gillで解く


 }
	sim_data.close_output_files(argv[1]);


}

