//
//  para.h
//  GA2017_08_07
//
//  Created by 西浦直人 on 2017/08/07.
//  Copyright © 2017年 西浦直人. All rights reserved.
//

#ifndef para_h
#define para_h

static const int MAXCELLNUMBER = 100;
static const int MAXPARASITENUMBER =0;
static const int MAXFILENUMBER = 100;//読み取るfileの数
static const double PI=3.14150265358979323846;
const int ProteinVAL =64;
const int Parasite_ProteinVAL =0;
const int check_num =1;//振動を判定する=checkpointの数
const double average_time=500.0;//平均を出す
const int MaxPathNum = 30;
const int MaxTargetPath = 8;//一つの遺伝子がtargetを制御できる上限
const int Species_numb = 1;//hostの種数
const int target_number= 8;//targetの数
const double mutation_rate = 0.0;//変異確率を与える
const int input_number = 8;//inputの数
const double input_strength=-5.0;//input_geneが受け取るmorphgenの大きさ
const double medium_V=1.0;
const double initial_volume=1.0;
const double death_volume=0.2;
const double S_resource=0.5;
const int death_time = 3000;
const double death_rate = 0.1;
const double S_D=0.05;
const double trspD=0.0;
const double Sigma=0.03;
const int EM_h=10;//EMでの分割数
const double growth_const=0.0002; //0.0005:T3-T4 0.00004:T5 0.00003:T6
const int histdiv =20;//ヒストグラムの分割数
const int elite_num=5;//次世代にそのまま残す個体数
const int mutation_num=4;//変異させる個体数
const double Host_mutaion_rate = 0.03;
const double Host_Sigma_mutaion_rate = 0.0;
const int selection_num = 50;//下位の集団は取り除く--> MAXCELLNUMBER-selection_numが消える
const int depth =2;
const int host_generation_times = 10;//hostのカウントを行う時間
const int parasite_generation_times_over = 1;//parasiteの世代時間がhostを越える時使う
const int parasite_generation_times = 100;//parasiteの世代時間 population dynamicsを再度走らせる
const double parasite_population_dynamics_cal_times=0.1;//population dynamicsを何秒間再度走らせるか
const double parasite_dt = 0.01;//population dynamicsを走らせる時の刻み時間
const int Parasite_genome_size = 1;
const int Parasite_species_num = 2;//2**Parasite_genome_size種類
const int Parasite_target_genome = 7;//attackする遺伝子の指定 Parasite_target_genome+Parasite_genome_size=target_number
const double maximum_virulence =3.0;
const int on =0;//使わないときは0にしておく
const int Glp_Num=15;//遺伝子集団の分割数それぞれのVipを求める
const double LAMBDA =5.0;
#endif /* para_h */




