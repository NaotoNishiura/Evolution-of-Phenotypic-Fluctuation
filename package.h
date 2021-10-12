//
//  package.h
//  Host-para_growth_noread2017_11_19
//
//  Created by 西浦直人 on 2018/06/25.
//  Copyright © 2018年 西浦直人. All rights reserved.
//


//
//  package.h
//  input_competition
//
//  Created by 西浦直人 on 2017/04/24.
//  Copyright © 2017年 西浦直人. All rights reserved.
//


#ifndef package_h
#define package_h

#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <new>
#include "para.h"

using namespace std;

class RK
{
public:
    RK(){};
    RK(double t0,int div,double tn0,int _protein_number,int _terget_protein_number,double mutTime,double diviTime,int NumberofCell,int Numberofonecelltype,double _beta,double _eps,double volume,double threshold);
    ~RK(void);
    RK(const RK& r);
    RK& operator=(const RK& r);
    void ShowDate();
    void Show();
    void archive();
    void archive_evo(int gene,int ith);
    void archive_evo_mini(int gene,int ith);
    void archive_Meanfieldval(int gene);
    void archive_GRN(int ith);
    void archive_GRN_every_gene(int ith,int time);
    void archive_GRN_every_gene_one(int ith,int time);
    void archive_YinJ(int ith);
    void archive_evo_mini_stationary(double time,int ith);
    void archive_Meanfieldval_stationary(int gene);
    void archive_GRN_stationary(int ith);
    double func_protein_g(double t,double *x,double *g_rand,int i);
    double func_protein(double t,double *x,int i);
    double func_protein_continuas(double t,double *x,int i);
    double func_protein_MF(double t,double *x,int i);
    double func_protein_signal(double t,double *x,int i);
    double func_protein_delution(double t,double *x,int i);
    void MutationJ(int ith,int lth);
    double MutationG(int ith,double sigma);
    double MutationD(int ith,double sigma);
    void calculateNS_protein();
    void calculateNS_protein_MF();
    void calculateNS_protein_g();
    void calculateNS_protein_g_EM(vector<vector< double> > g);
    void calculateNS_protein_delution_g_EM(vector<vector< double> > g);
    void calculateNS_protein_signal();
    void calculateNS_protein_signal_g_EM(vector<vector< double> > g);
    void failclean();
    void failclean_evo();
    void converterJ();//行列を隣接リストに変換
    
    int rand();
    int rand_onoff();
    int rand_select(int a,int b,int c);
    void calculateNS_volume();
    double func_volume(double t,double x);
    double gauss_rand(double mu,double sigma);
    double gabs_gauss_rand(double mu,double sigma);
    void calculate_Mu();
    void calculate_growth();
    void calculate_growth_dist();
    void calculate_Fitness();
    void Calculate_Fitness_Target();
    void Calculate_Fitness_Target_parasite();
    void Calculate_Stack_fitness_Target_parasite_average();
    void Calculate_Stack_fitness_Target_noparasite_average();
    void calculate_Input_list();
    
    //以下、GeterとSetter
    void SetResult(double x,int lth);
    void SetResult_average(double rl,int lth);
    void SetResult_check(double rc,int lth,int check_t);
    void SetJ(double j,int ith,int lth);
    void SetEL(int j,int ith,int lth);
    void SetWEL(int j,int ith,int lth);
    void SetSEL(int ith,int size);//ith proteinの手の数
    void SetREL(int j,int ith,int lth);
    void SetRSEL(int lth,int size);//lthから出るpathの数
    
    void SetProteinD(double diff,int lth);
    void SetThreshold_g(double g,int lth);
    void SetTime(double nowtime);
    void SetTime1000();
    void SetGene(int i);
    void SetSorting_number(int num);
    void SetMeanfieldval(double meanval,int lth);
    void SetTotalcellnumber(int new_totalnumber);
    void SetVolume(double volume);
    void SetGrowth(double growth);
    void SetK_delta(int onoff,int ith);
    void SetInput_list(int onoff,int ith,int lth);
    void SetTargetPattern(int onoff,int lth);
    void SetInfection_virulence(double vir,int lth);
    void SetQ(double probability);
    void SetCelltype(int type);
    void SetAverage_growth(double averagegrowth);
    void SetInitialSpecies(int species_numb);
    void SetDivision_number(int numb);
    void SetMutation_number(int numb);
    void SetMu(double _mu);
    void SetHammingD(double ha);
    void SetRelativeHammingD(double Rha);
    void SetStationaryResult(double x,int lth);
    void SetStationaryTime(double sta_time);
    void SetStationaryJ(int j,int ith,int lth);
    void SetStationaryGrowth(double growth_sta);
    void SetFitness(double F);
    void SetFitness_p(double F_p);
    void SetCumulativ_p(double cp);
    void SetPara_fitness(double p_fitness);
    void SetFitness_Target_parasite_average(double tpa_fitness);
    void SetStack_fitness_Target_parasite_average(double stpa_fitness);//平均値計算用の前世代の入れ物
    void SetFitness_Target_noparasite_average(double tpa_fitness);
    void SetStack_fitness_Target_noparasite_average(double stpa_fitness);//平均値計算用の前世代の入れ物
    void SetSigma_Str(double str);
    
    double GetResult(int lth) const;//result[lth] protein_x
    double GetResult_average(int lth) const;//平均の発現量
    double GetResult_check(int lth,int check_t) const;
    
    double GetStationaryResult(int lth) const;//定常状態のresult
    double GetJ(int ith,int lth) const;//J[ith][lth] from protein l to i
    int GetStationaryJ(int ith,int lth) const;//定常状態の遺伝子
    int GetK_delt(int lth);//K_delt[lth][lth] input gene
    int GetInput_list(int ith,int lth);
    int GetInput_listSize(int ith);
    int GetTargetPattern(int lth);
    double GetInfection_virulence(int lth);
    inline int GetEL(int ith,int lth) const;
    inline int GetWEL(int ith,int lth) const;
    inline int GetELSize(int ith) const;
    inline int GetSEL(int ith) const;
    inline int GetREL(int ith,int lth) const;//REL(i,l) iから出るpathを記録　ただしREL[i][l]と表記
    inline int GetRELSize(int ith) const;
    inline int GetRSEL(int ith) const;
    double GetDelt() const;
    double GetProD(int lth) const;
    double GetThreshold_g(int lth) const;
    double GetBeta() const;
    int GetNum() const;
    int GetCelltype() const;
    int Getprotein_x_number() const;
    int GetTotalcellnumber() const;
    double GetMeanfieldval(int lth) const;
    double GetTime() const;
    double GetMutationTime() const;
    double GetDiviTime() const;
    int GetTime1000() const;
    int GetGene() const;
    int GetSorting_num() const;
    int GetMutationTime1000() const;
    int GetDiviTimeTime1000() const;
    double GetVolume() const;
    double GetInitialVolume() const;
    double GetGrowthRate() const;
    double GetQ() const;
    double GetAverage_Growth() const;
    int GetInitialSpecies() const;
    int GetMutation_number() const;
    int GetDivision_number() const;
    double GetHammingD() const;
    double GetRelativeHammingD() const;
    double GetMu() const;
    double GetStationaryTime() const;
    double GetStationaryGrowth() const;
    char GetTree_history() const;
    double GetFitness() const;
    double GetFitness_p() const;//fitnessから得た確率
    double GetCumulative_p() const;//fitnessから得たGA用確率
    double GetFinishTime() const;
    double GetPara_fitness() const;
    double GetFitness_Target_parasite_average() const;
    double GetStack_fitness_Target_parasite_average() const;
    double GetFitness_Target_noparasite_average() const;
    double GetStack_fitness_Target_noparasite_average() const;
    double GetSigma_Str() const;
    string tree_history;//種ごとの子孫を把握するため
    string ori_Sstr_J;//初期のネットワークを記憶
    string Sstr_J;//進化後のネットワークを記憶
    
    //以下Medium　input物質
    void SetY(double Y_med,int ith);//input物質が存在するmedium
    double GetY(int ith) const;
    void SetYin_J(int yin_j,int ith,int lth);//input物質のYin_Jのsetter
    int GetYin_J(int ith,int lth) const;//input物質のYin_Jのgetter
    void SetMedium_D(double Y_input_D,int ith);
    double GetMedium_D(int ith) const;
    void SetYin(double Y_in,int ith);//cell中のinput物質
    double GetYin(int ith) const;
    void calculateNS_Yin();
    double func_Yin(double t,double *y,int i);
    void SetX(double X_med,int ith);//拡散物質が存在するmedium
    double GetX(int ith) const;
    
    //parasite用
    int GetParasite_attacked_number() const;
    void SetParasite_attacked_number(int num);
    double GetParasite_population() const;
    void SetParasite_population(double pp);
    double GetParasite_population_rate(int ith);//全parasite情報を与える
    void SetParasite_population_rate(double pp,int ith);//全parasite情報を与える
    
    //VipVg測定用
    void calculateNS_protein_g_EM_VipVg(vector<vector< double> > g);//VipVg計算用　通常の計算とは分けておく
    double func_protein_VipVg(double t,double *x,int i);
    void SetResult_VipVg(double x,int lth);
    inline double GetResult_VipVg(int lth) const;//Result_VipVgを使う
    
    void converterJ_VipVg();//行列を隣接リストに変換
    
    void SetJ_VipVg(int j,int ith,int lth);
    int GetJ_VipVg(int ith,int lth) const;//J[ith][lth] from protein l to i
    
    inline int GetEL_VipVg(int ith,int lth) const;
    inline int GetELSize_VipVg(int ith) const;
    inline int GetWEL_VipVg(int ith,int lth) const;
    inline int GetSEL_VipVg(int ith) const;
    inline int GetREL_VipVg(int ith,int lth) const;//REL(i,l) iから出るpathを記録　ただしREL[i][l]と表記
    inline int GetRELSize_VipVg(int ith) const;
    inline int GetRSEL_VipVg(int ith) const;
    
    void SetEL_VipVg(int j,int ith,int lth);
    void SetWEL_VipVg(int j,int ith,int lth);
    void SetSEL_VipVg(int ith,int size);//ith proteinの手の数
    void SetREL_VipVg(int j,int ith,int lth);
    void SetRSEL_VipVg(int lth,int size);//lthから出るpathの数
    
private:
    
    //以下protein
    vector<double> result;
    vector<double> result_average;//400~500秒の平均
    vector<vector<double > > result_check;//400~500秒の平均
    vector<vector<double> > J;
    vector<vector<int> > EL;//J(i,l) iへ繋がるpathを記録
    vector<vector<int> > WEL;
    vector<int> SEL;
    vector<vector<int> > REL;//J(i,l) jから出るpathを記録　ただしREL[j][i]と表記
    vector<int> RSEL;
    vector<int> k_delta;//クロネッカー
    vector<vector<int> > input_list;//inputの繋ぎ先を記録
    
    vector<double> Meanfieldval;
    vector<double> proteinD;
    vector<double> threshold_g;
    vector<double> StationaryResult;//定常状態の発現量を記録するよう
    vector<vector<int> > StationaryJ;//定常状態の遺伝子を記録するよう
    vector<int> TergetNumber;
    vector<int> TargetPattern;//fitnessをパターンにする時
    vector<double> Parasite_population_rate;//signalを感知するとした時のparasiteの情報
    vector<double> infection_virulence;//パラサイトに対して距離をとったものの総和
    
    int LOOP;
    int protein_x_number;
    int terget_protein_number;
    int number;
    int totalcellnumber;
    static int objnum;
    int celltype;
    int Initial_species;
    int division_number;
    int mutaion_number;
    int Dist;//類似度
    double hammingD;
    double RelativehammingD;//0thCELLから見たHamming距離
    int gene;//世代数
    int sorting_number;//hamming距離でsortした時の番号
    int parasite_attacked_number;//parasiteに対応するhostの種類
    double delt,time,startTime,finishTime,mutationTime,divisionTime,beta,v,eps,initial_volume,mu;
    int time1000,matationTime1000,divisionTime1000;//一秒に対応させる
    double growth_rate,average_growth,q;//q選択確率
    double StationaryTime,StationaryGrowth;//定常状態の時間を記録するよう
    double fitness,fitness_p,cumulative_p,fitness_Target_parasite_average,stack_fitness_Target_parasite_average;
    double fitness_Target_noparasite_average,stack_fitness_Target_noparasite_average;
    double para_fitness;
    double parasite_population,sigma_str;//attackしてくるparasiteの割合
    
    //以下Medium input物質
    vector<vector<int> > Yin_J;//input物質の制御を決める
    vector<double> Yin;//cell中のinput物質の量
    vector<double> Y;//medium中のinput物質の量
    vector<double> X;//medium中の拡散タンパク物質の量
    vector<double> Medium_D;//cell-medium間のinput物質の拡散係数
    
    //VipVg計算用
    vector<double> result_VipVg;//Vip計算用に入れ物を分けておく
    vector<vector<int> > J_VipVg;
    
    vector<vector<int> > EL_VipVg;//J(i,l) iへ繋がるpathを記録
    vector<vector<int> > WEL_VipVg;
    vector<int> SEL_VipVg;
    vector<vector<int> > REL_VipVg;//J(i,l) jから出るpathを記録　ただしREL[j][i]と表記
    vector<int> RSEL_VipVg;
};

int RK::objnum =0;


struct data_t {
    int num=0;
    double hammD=0;
    double Relative_hammD=0;//0thCELLから見たHamming距離
    string str;
    
    // 最後のconstを忘れると"instantiated from here"というエラーが出てコンパイルできないので注意
    bool operator<( const data_t& right ) const {
        //return hammD == right.hammD ? Relative_hammD < right.Relative_hammD : hammD < right.hammD;
        //InitialGRNからhamming距離でsort
        return Relative_hammD == right.Relative_hammD ? hammD < right.hammD : Relative_hammD < right.Relative_hammD;
        //0thCELLから見たHamming距離でsort
        //return hammD == right.hammD ? str < right.str : hammD < right.hammD;//first geneからhamming距離でソート
        //return str == right.str ? hammD < right.hammD : str < right.str;
    }
};

struct data_t_sta {
    int num;
    int hammD=0;
    int Relative_hammD=0;//0thCELLから見たHamming距離
    
    bool operator<( const data_t_sta& right ) const {
        //return hammD == right.hammD ? Relative_hammD < right.Relative_hammD : hammD < right.hammD;//InitialGRNからhamming距離でsort
        return Relative_hammD == right.Relative_hammD ? hammD < right.hammD : Relative_hammD < right.Relative_hammD;//0thCELLから見たHamming距離でsort
    }
};//定常状態の値の記録用

struct data_t_fitness{
    int num=0;
    double fit=0;
    int hammD=0;//0thCELLから見たHamming距離
    
    // 最後のconstを忘れると"instantiated from here"というエラーが出てコンパイルできないので注意
    bool operator<( const data_t_fitness& right ) const {
        return fit == right.fit ? hammD < right.hammD : fit < right.fit;//fitnessでsorts
    }
};

class CELL
{
private:
    double _t0,_tn0,_mutTime,_diviTime,_beta,_eps,_volume,_threshold,statinary_time,genome_variance,fitness_variance,mean_fitness;
    //statinary_time定常状態カウント用の区間
    int _div,_protein_number,_cellnumber,terget,_terget_protein_number,_Numberofonecelltype,generation,flag_of_VipVg;
    
    vector<vector<double> > sim;//phenotypeの類似度
    vector<vector<double> > genotype_sim;//genotypeの類似度
    vector<vector<double> > Jij_variance;//jijの分散
    vector<vector<double> > Jij_initial;//初期のJijを保存する
    vector<data_t> data_array;//sort用に使う入れ物
    vector<data_t_sta> data_array_sta;//定常状態の値の記録用
    vector<data_t_fitness> data_array_fitness;//fitnessごとにsort
    vector<double> variance ;
    vector<double> variance_sta;
    vector<double> variance_sta_Vg;
    vector<double> variance_average;
    vector<double> variance_threshold;
    vector<double> fitness_variance_by_type;//系列ごとの分散
    vector<double> mean_fitness_by_type;//系列ごとの分散
    vector<double> mean_fitness_nontarget_by_type;//攻撃を受けないtargetfitnessの平均
    vector<double> efficincy;//fitness　資源の利用効率
    vector<vector<int> > histogram;//1世代前のヒストグラム
    
    //VipVg測定用
    vector<vector<double> > variance_gloup_Vip;//Glp_Numだけ用意した遺伝子型の分散
    vector<vector<double> > variance_gloup_Vg;//Glp_Numだけ用意した
    vector<vector<double> > variance_gloup_Vip_average;//Glp_Numだけ用意した遺伝子型の分散
    vector<vector<double> > variance_gloup_Vip_average_stack;//Glp_Numだけ用意した遺伝子型の分散
    
    //Parasite用
    vector<double> hosto_num;//マッチするhostoの数
    vector<double> parasite_population;//parasiteの数
    vector<double> parasite_virulence;//parasiteの毒性の強さ
    vector<vector<double> > parasite_genotype;//parasiteの遺伝子型
    vector<vector<int> > parasite_list;//hamming distanceが1になるもの
    
    std::mt19937 create_rand_engine(){
        std::random_device rnd;
        std::vector<std::uint_least32_t> v(10);// 初期化用ベクタ
        std::generate(v.begin(), v.end(), std::ref(rnd));// ベクタの初期化
        std::seed_seq seed(v.begin(), v.end());
        return std::mt19937(seed);// 乱数エンジン
    }
    
    std::vector<int> make_rand_array_unique(const size_t size, int rand_min, int rand_max){
        if(rand_min > rand_max) std::swap(rand_min, rand_max);
        const size_t max_min_diff = static_cast<size_t>(rand_max - rand_min + 1);
        if(max_min_diff < size) throw std::runtime_error("引数が異常です");
        
        std::vector<int> tmp;
        auto engine = create_rand_engine();
        std::uniform_int_distribution<int> distribution(rand_min, rand_max);
        
        const size_t make_size = static_cast<size_t>(size*1.2);
        
        while(tmp.size() < size){
            while(tmp.size() < make_size) tmp.push_back(distribution(engine));
            std::sort(tmp.begin(), tmp.end());
            auto unique_end = std::unique(tmp.begin(), tmp.end());
            
            if(size < std::distance(tmp.begin(), unique_end)){
                unique_end = std::next(tmp.begin(), size);
            }
            tmp.erase(unique_end, tmp.end());
        }
        
        std::shuffle(tmp.begin(), tmp.end(), engine);
        return std::move(tmp);
    }
    
public:
    CELL(double __t0,int __div,double __tn0,int __protein_number,int __terget_protein_number,int __Numberofonecelltype,double __mutTime,double __diviTime,int __cellnumber,double __beta,double __eps,double __medium_D,double __volume,double __threshold)
    :_t0(__t0),_tn0(__tn0),_mutTime(__mutTime),_diviTime(__diviTime),_div(__div),_protein_number(__protein_number),_terget_protein_number(__terget_protein_number),_cellnumber(__cellnumber),_Numberofonecelltype(__Numberofonecelltype),_beta(__beta),_eps(__eps),_volume(__volume),_threshold(__threshold)
    {
        sim.shrink_to_fit();
        sim.resize(MAXCELLNUMBER);
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            sim[i].resize(MAXCELLNUMBER);
        }
        
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            for(int l=0;l<MAXCELLNUMBER;l++)
            {
                sim[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        Jij_variance.shrink_to_fit();
        Jij_variance.resize(ProteinVAL);
        for(int i=0;i<ProteinVAL;i++)
        {
            Jij_variance[i].resize(ProteinVAL);
        }
        
        Jij_initial.shrink_to_fit();
        Jij_initial.resize(ProteinVAL);
        for(int i=0;i<ProteinVAL;i++)
        {
            Jij_initial[i].resize(ProteinVAL);
        }
        
        for(int i=0;i<ProteinVAL;i++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                Jij_variance[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        genotype_sim.shrink_to_fit();
        genotype_sim.resize(MAXCELLNUMBER);
        
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            genotype_sim[i].resize(MAXCELLNUMBER);
        }
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            for(int l=0;l<MAXCELLNUMBER;l++)
            {
                genotype_sim[i][l] =0.0;//genoetype_simを０で初期化
            }
        }
        
        data_array.shrink_to_fit();
        data_array.resize(MAXCELLNUMBER);
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            data_array[i].num = i;
            data_array[i].hammD = 0;
            data_array[i].str = {};
        }
        
        data_array_sta.shrink_to_fit();
        data_array_sta.resize(MAXCELLNUMBER);
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            data_array_sta[i].num = i;
            data_array_sta[i].hammD = 0;
            data_array_sta[i].Relative_hammD = 0;
        }
        
        data_array_fitness.shrink_to_fit();
        data_array_fitness.resize(MAXCELLNUMBER);
        for(int i=0;i<MAXCELLNUMBER;i++)
        {
            data_array_fitness[i].num = i;
            data_array_fitness[i].hammD = 0;
            data_array_fitness[i].fit =0;
        }
        
        variance.shrink_to_fit();
        variance.resize(_protein_number);
        for(int i=0;i<_protein_number;i++)
        {
            variance[i] =0.0;
        }
        
        variance_threshold.shrink_to_fit();
        variance_threshold.resize(_protein_number);
        for(int i=0;i<_protein_number;i++)
        {
            variance_threshold[i] =0.0;
        }
        
        variance_sta.shrink_to_fit();
        variance_sta.resize(_protein_number);
        for(int i=0;i<_protein_number;i++)
        {
            variance_sta[i] =0.0;
        }
        
        
        variance_sta_Vg.shrink_to_fit();
        variance_sta_Vg.resize(_protein_number);
        for(int i=0;i<_protein_number;i++)
        {
            variance_sta_Vg[i] =0.0;
        }
        
        variance_gloup_Vip.shrink_to_fit();
        variance_gloup_Vip.resize(Glp_Num);
        for(int i=0;i<Glp_Num;i++)
        {
            variance_gloup_Vip[i].resize(ProteinVAL);
        }
        
        for(int i=0;i<Glp_Num;i++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                variance_gloup_Vip[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        variance_gloup_Vip_average.shrink_to_fit();
        variance_gloup_Vip_average.resize(Glp_Num);
        for(int i=0;i<Glp_Num;i++)
        {
            variance_gloup_Vip_average[i].resize(ProteinVAL);
        }
        
        for(int i=0;i<Glp_Num;i++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                variance_gloup_Vip_average[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        variance_gloup_Vip_average_stack.shrink_to_fit();
        variance_gloup_Vip_average_stack.resize(Glp_Num);
        for(int i=0;i<Glp_Num;i++)
        {
            variance_gloup_Vip_average_stack[i].resize(ProteinVAL);
        }
        
        for(int i=0;i<Glp_Num;i++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                variance_gloup_Vip_average_stack[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        variance_gloup_Vg.shrink_to_fit();
        variance_gloup_Vg.resize(Glp_Num);
        for(int i=0;i<Glp_Num;i++)
        {
            variance_gloup_Vg[i].resize(ProteinVAL);
        }
        
        for(int i=0;i<Glp_Num;i++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                variance_gloup_Vg[i][l] =0.0;//phenotype_simを０で初期化
            }
        }
        
        fitness_variance_by_type.shrink_to_fit();
        fitness_variance_by_type.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            fitness_variance_by_type[i] =0.0;
        }
        
        mean_fitness_by_type.shrink_to_fit();
        mean_fitness_by_type.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            mean_fitness_by_type[i] =0.0;
        }
        
        mean_fitness_nontarget_by_type.shrink_to_fit();
        mean_fitness_nontarget_by_type.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            mean_fitness_nontarget_by_type[i] =0.0;
        }
        
        efficincy.shrink_to_fit();
        efficincy.resize((target_number));
        for(int i=0;i<target_number;i++)
        {
            efficincy[i] =0.0;
        }
        
        histogram.shrink_to_fit();
        histogram.resize(histdiv);
        
        for(int i=0;i<histdiv;i++)
        {
            histogram[i].resize(target_number+Parasite_genome_size);
        }
        for(int i=0;i<histdiv;i++)
        {
            for(int j=0;j<target_number+Parasite_genome_size;j++)
            {
                histogram[i][j] = 0;
            }
        }
        for(int i=0;i<_cellnumber;i++)
        {
            new(rk + i) RK(_t0,_div,_tn0,_protein_number,_terget_protein_number,_mutTime,_diviTime,_cellnumber,_Numberofonecelltype,_beta,_eps,_volume,_threshold);
        }
        
        //parasiteの変数
        hosto_num.shrink_to_fit();
        hosto_num.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            hosto_num[i] = 0;
        }
        
        parasite_population.shrink_to_fit();
        parasite_population.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            parasite_population[i] = (double)(1.0/Parasite_species_num);//全種均等に与える
            //parasite_population[i] = (double)(MAXPARASITENUMBER/Parasite_species_num);//全種均等に与える
        }

        parasite_virulence.shrink_to_fit();
        parasite_virulence.resize(Parasite_species_num);
        for(int i=0;i<Parasite_species_num;i++)
        {
            parasite_virulence[i] = maximum_virulence;//全種初期値は同じ
        }
        
        parasite_genotype.shrink_to_fit();
        parasite_genotype.resize(Parasite_species_num);
        
        for(int i=0;i<Parasite_species_num;i++)
        {
            parasite_genotype[i].resize(Parasite_genome_size);
        }
        
       parasite_genotype[0][0]=0.0;
        parasite_genotype[0][1]=0.0;
        //parasite_genotype[0][2]=0.0;
        parasite_genotype[1][0]=1.0;
        parasite_genotype[1][1]=1.0;
        //parasite_genotype[1][2]=1.0;
        //parasite_genotype[2][0]=0.0;
        //parasite_genotype[2][1]=1.0;
        //parasite_genotype[2][2]=0.0;
        //parasite_genotype[3][0]=0.0;
        //parasite_genotype[3][1]=1.0;
        //parasite_genotype[3][2]=1.0;
        //parasite_genotype[4][0]=1.0;
        //parasite_genotype[4][1]=0.0;
        //parasite_genotype[4][2]=0.0;
        //parasite_genotype[5][0]=1.0;
        //parasite_genotype[5][1]=0.0;
        //parasite_genotype[5][2]=1.0;
        //parasite_genotype[6][0]=1.0;
        //parasite_genotype[6][1]=1.0;
        //parasite_genotype[6][2]=0.0;
        //parasite_genotype[7][0]=1.0;
        //parasite_genotype[7][1]=1.0;
        //parasite_genotype[7][2]=1.0;
        
        parasite_list.shrink_to_fit();
        parasite_list.resize(Parasite_species_num);
        
        for(int i=0;i<Parasite_species_num;i++)
        {
            parasite_list[i].resize(Parasite_genome_size);
        }
        
        for(int i=0;i<Parasite_species_num;i++)
        {
            for(int l=0;l<Parasite_genome_size;l++)
            {
                parasite_list[i][l] =0;
            }
        }
        
        genome_variance = 0.0;
        fitness_variance = 0.0;
        mean_fitness = 0.0;
    }
    
    ~CELL(){cout<<"Class CELL is destroyed"<<endl;};
    RK  rk[MAXCELLNUMBER];//細胞の分クラス作成
    RK  rk_VipVg[MAXFILENUMBER];//VipVg用のクラスの作成
    
    void newcellmaker();//新しい細胞を一つ作成
    void removecell(int ith);//細胞を一つ削除
    void death_check(int ith);
    void death_event();
    void cell_death(int ith);//ith cellが死亡する
    void death_volume_check();
    void division_time(int ith,double simga);//time1000を使用
    void division_volume(int ith,double sigma);
    void division_volume_mutation(int ith,double sigma);
    void division_volume_mutation_network(int ith,double sigma,int time);
    void division_volume_mutation_network_continuas(int ith,double sigma,int time);
    void division_volume_mutation_network_sigma(int ith,double sigma,int time);
    void division_volume_mutation_network_couple(int ith,double sigma,int time);
    void division_fitness_mutation(double sigma,int time);
    
    void archive_population(int time);
    void archive_target_input(int gene);
    void archive_average(int time);
    void archive_sort_cell(int time);
    void archive_growth_volume(int time);
    void check_tree_history(int mother,int daughter);
    void calculate_meanfield();
    void calculate_NFDS_growth();
    void calculate_Result_average();//result_averageを与える
    void calculate_Fitness_Target_parasite_average();//fitnessの平均を与える
    void calculate_Fitness_Target_noparasite_average();//non attacked genesのfitnessの平均を与える
    void calculate_sim(int ith,int jth);
    void calculate_simJij(int ith,int jth);
    void calculate_EuclidDistanceJij(int ith,int jth);
    double calculate_EuclidDistanceJij_initial(int ith);
    void calculate_sim_sta(int ith,int jth);
    double GetSim(int ith,int lth) const;//sim[ith][lth]
    void SetSim(double e,int ith,int lth);//類似度'e'をset
    void archive_sim(int time);//類似度記録
    void archive_sim_sta(int time);//定常状態の類似度を記録
    void genosim_label_set();//hamming距離でのsortingのlabelをrkにつける
    double GetJij_variance(int ith,int lth) const;
    void SetJij_variance(double v,int ith,int lth);
    double GetJij_initial(int ith,int lth) const;
    void SetJij_initila(double v,int ith,int lth);
    void archive_Jij_variance(int time);//類似度記録
    void caluculate_genotype_sim(int ith,int lth);
    double GetGenotype_sim(int ith ,int lth) const;
    void SetGenotype_sim(double e,int ith,int lth);
    void archive_genotype_sim(int time);//類似度記録
    void SetData_array();
    void SetVariance(double var,int i);
    void SetVariance_thereshold(double var,int i);
    double GetVariance(int i);
    double GetVariance_thereshold(int i);
    void SetFitness_variance_by_type(double var,int i);
    double GetFitness_variance_by_type(int i);
    void SetMean_fitness_by_type(double mean,int i);
    double GetMean_fitness_by_type(int i);
    double GetMean_fitness_nontarget_by_type(int i);
    void SetMean_fitness_nontarget_by_type(double mean,int i);
    void archive_var_genome(int time);
    void SetGenome_Variance(double var);
    double GetGenome_Variance() const;
    void SetFitness_variance(double var);
    double GetFitness_variance() const;
    void SetMean_fitness(double mean);
    double GetMean_fitness() const;
    void Calculate_Variance();
    void Calculate_Variance_sta();
    void Calculate_Variance_threshold();
    void Calculate_genome_Variance();
    void Calculate_Jij_Variance();
    void archive_Jij_Variance(int time);//行列の分散
    void Calculate_fitness_Variance();
    void Calculate_fitness_Variance_by_type();//host系列ごとのfitness_varianceを計算
    void archive_var(int time);
    void archive_var_sta(int time);
    void archive_var_threshold(int time);
    void archive_Sigma_Str_all(int time);
    void archive_Sigma_Str(int time);
    void archive_Parasite_Virulence(int time);
    void archive_threshold(int time);
    void SetData_array_sta();
    void SetData_array_fitness();
    void Calculate_Cumulative_p();//fitnessから親選びのための累積確率を作る
    void Calculate_Cumulative_p_power();//fitnessから親選びのための累積確率を作る 10**fitness
    int  Roulette();//子供を作る親を選ぶ関数
    void elite(int time);//エリート戦略
    void elite_noinput(int time);//エリート戦略inputは変異しない
    void elite_one(int time);//エリート戦略
    void elite_noinput_one(int time);//エリート戦略inputは変異しない
    void Selection_Pressure_one(int time);//選択圧の変更  上位の個体だけを残す
    void Selection_Pressure_noinput_one(int time);//選択圧の変更
    void Selection_Pressure_norule(int time);//nodeの選択の制限を入れない
    void Selection_WF_norule(int time);//WFで選択的に子供作成
    void Selection_WF_self_regulation(int time);//WFでtargetの自己制御を許す
    void Initial_Jij_maker();//初期のネットワークがfitness=0のとき作り直す
    
    void reset();//初期条件に戻す関数
    void Calculate_Efficiency();//利用効率の計算
    void histo(int time);//ヒストグラム作成用
    int GetHistogram(int k,int j);
    void SetHistogram(int his,int k,int j);
    
    //VipVg測定用
    void Calculate_VipVg(int time);//ここでVipVgの計算をまとめて行う
    void Calculate_Variance_random_Vip();//複数の遺伝子型から分散を求める時
    void Calculate_Variance_random_Vg();//複数の遺伝子型から分散を求める時
    void Calculate_Variance_random_Vip_average(int time);//複数の遺伝子型から分散を求める時
    void Calculate_Variance_random_Vip_average_stack();//複数の遺伝子型から分散を求める時
    
    void SetVariance_Sta(double var,int i);
    double GetVariance_Sta(int i);
    void SetVarianceVg(double var,int i);
    double GetVarianceAverage(int i);
    void SetVarianceAverage(double var,int i);
    double GetVarianceVg(int i);
    void SetVariance_Gloup_Vip(int ith,int Protein_jth,double var);
    void SetVariance_Gloup_Vg(int ith,int Protein_jth,double var);
    void SetVariance_gloup_Vip_average(int ith,int Protein_jth,double var);
    void SetVariance_gloup_Vip_average_stack(int ith,int Protein_jth,double var);
    double GetVariance_Gloup_Vip(int ith,int Protein_jth);
    double GetVariance_Gloup_Vg(int ith,int Protein_jth);
    double GetVariance_gloup_Vip_average(int ith,int Protein_jth);
    double GetVariance_gloup_Vip_average_stack(int ith,int Protein_jth);
    void robustness_function(int ith);//ith細胞のrobustnessを調べる
    void robustness_random_cell();//ith細胞のrobustnessを調べる
    void Vg_cell_make();//遺伝子がほんの少し違う集団を用意
    void Vg_cell_make_random();
    void Load_Jij();//50世代ごとの遺伝子型をcopy
    void Load_Jji_random();//複数の遺伝子型を選んでcopyを作成
    void ReLoad_Jji();//Load_Jji_random()で選んだ遺伝子型に振り直す
    void Load_Jij_typical();//集団の典型的な遺伝子型
    void SetGeneration(int gene);
    int  GetGeneration();
    void SetFlag_of_VipVg(int flag);
    int Getflag_of_VipVg();
    void archive_var_Vip(int time);
    void archive_var_Vip_average(int time);
    
    //Parasite用
    double GetParasite_population(int ith) const;//ithはparasiteの種番号
    void SetParasite_population(int ith,double pp);
    double GetParasite_Genotype(int ith,int jth) const;//ithはparasiteの種番号
    void SetParasite_Genotype(int ith,int jth,double pg);//ithはparasiteの種番号
    double GetHosto_num(int ith);//ithはparasiteの種番号
    void SetHosto_num(int ith,int hn);
    double GetParasite_virulence(int ith) const;//ithはparasiteの種番号
    void SetParasite_virulence(int ith,double vir);
    void calculate_parasite_dynamics();//ここで少しpopulationの計算を行う
    void calculate_parasite_dynamics_evo();//ここで少しpopulationの計算を行う
    void calculate_parasite_dynamics_dist_fitness();//ここで少しpopulationの計算を行う 距離を導入
    void calculate_parasite_dist_fitness();//全細胞の感染をrateを計算
    void population_dynamics();
    void population_dynamics_dist_fitness();
    double func_population_dynamics(double *y,int ith);
    double func_population_dynamics_dist_fitness(double *y,int ith);
    void caluculate_hosto_num();
    void archive_parasite_population(int time);
    void SetParasite_list(int list,int ith,int jth);
    int GetParasite_list(int ith,int jth);
    void parasite_list_maker();
    
};

class Medium
{
public:
    vector<double> Y;//栄養成分
    vector<double> X;//細胞由来の成分
    double Medium_D,delt,time,startTime,finishTime,V,S_out;
    int cell_number,LOOP;
    Medium(){};
    Medium(double t0,int div,double tn0,double medium_D,int _cell_number,double Volume,double S);
    ~Medium(){cout<<"Class Medium is destroyed"<<endl;};
    void SetMedium_D(double medium_D);
    double GetMedium_D() const;
    void SetY(double Y_med,int ith);
    double GetY(int ith) const;
    double GetX(int ith) const;
    void SetX(double X_med,int ith);
    void calculateNS_med(RK *rk);
    double func_medium(double t, double *y, int i,RK *rk);//栄養成分
    double func_medium_X(double t, double *x, int i,RK *rk);//細胞由来の成分
    double func_medium_V(double t,double volume);
    void archive_medium(int time);
    void reset();//初期条件に戻す関数
};




#endif /* package_h */




