#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <new>
#include "package.h"
#include "para.h"

//const double Route = sqrtl(ProteinVAL*0.4);//25:25:50の時
//const double Route = 1.0;//25:25:50の時
//const double Route = sqrtl(ProteinVAL);//25:25:50の時
const double Route = sqrtl(ProteinVAL*0.3);//15:15:70の時
//const double Route = sqrtl(ProteinVAL*0.5);//30:30:40の時
//const double Route = sqrtl(ProteinVAL*0.1);//10:10:80の時

using namespace std;

int sstoi(std::string str){
    int ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    return ret;
}

double sstod(std::string str){
    double ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    return ret;
}

int iHammingDist(std::string sStr1,std::string sStr2){
    
    /*ハミング距離カウント用変数*/
    int iDist = 0;
    
    for(int i=0;i<sStr1.size();i++)
    {
        if(sStr1[i] != sStr2[i]){
            ++iDist;
        }
    }
    
    /*ハミング距離を戻り値にする*/
    return iDist;
    
}

void CELL::newcellmaker()
{
    int ith = 0;
    ith = rk[0].GetTotalcellnumber();//総細胞数の次に一つ追加
    new(rk + ith) RK(_t0,_div,_tn0,_protein_number,_terget_protein_number,_mutTime,_diviTime,_cellnumber,_Numberofonecelltype,_beta,_eps,_volume,_threshold);
    ith = ith + 1;//Totallcellnumberに一つ足した数
    for(int i=0;i<ith;i++)
    {
        rk[i].SetTotalcellnumber(ith);//分裂分の細胞数の更新
    }
}

void CELL::removecell(int ith)
{
    int lastnum = 0,p_val=0,cellnum=0;
    double nowtime =0.0;
    
    lastnum = rk[0].GetTotalcellnumber()-1;//総細胞数の次に一つ追加
    cellnum = rk[0].GetTotalcellnumber();
    p_val =rk[0].Getprotein_x_number();
    nowtime = rk[lastnum].GetTime();
    
    rk[ith].SetStationaryTime(rk[lastnum].GetStationaryTime());//分裂時間を更新
    rk[ith].SetStationaryGrowth(rk[lastnum].GetStationaryGrowth());//分裂時間を更新
    for(int j=0;j<p_val;j++)
    {
        rk[ith].SetStationaryResult(rk[lastnum].GetStationaryResult(j),j);//分裂時の発現量を更新
        for(int l=0;l<p_val;l++)
        {
            rk[ith].SetStationaryJ(rk[lastnum].GetStationaryJ(j, l), j, l);//分裂時のネットワークを更新
        }
    }
    for(int l=0;l<p_val;l++)
    {
        rk[ith].SetResult(rk[lastnum].GetResult(l),l);
        rk[ith].SetProteinD(rk[lastnum].GetProD(l), l);
        rk[ith].SetThreshold_g(rk[lastnum].GetThreshold_g(l), l);
        rk[ith].SetMeanfieldval(rk[lastnum].GetMeanfieldval(l),l);
        
        for(int k=0;k<p_val;k++)
        {
            rk[ith].SetJ(rk[lastnum].GetJ(l, k), l, k);
        }
    }
    rk[ith].converterJ();
    rk[ith].SetVolume(rk[lastnum].GetVolume());
    rk[ith].SetCelltype(rk[lastnum].GetCelltype());//ithのcelltype引き継ぎ
    rk[ith].SetDivision_number(rk[lastnum].GetDivision_number());//親のカウントを増やす
    rk[ith].tree_history.clear();
    rk[ith].tree_history = rk[lastnum].tree_history;
    rk[ith].ori_Sstr_J.clear();
    rk[ith].ori_Sstr_J.reserve(p_val*p_val);
    for(int i=0;i<p_val*p_val;i++)
    {
        rk[ith].ori_Sstr_J.push_back(rk[lastnum].ori_Sstr_J[i]);
    }
    for(int k=0;k<input_number;k++)
    {
        rk[ith].SetYin(rk[lastnum].GetYin(k), k);//分裂時に子供のinput物質は0に
        rk[ith].SetY(rk[lastnum].GetY(k), k);//medium中のinput物質は引き継ぎ
        rk[ith].SetX(rk[lastnum].GetX(k), k);//medium中の拡散物質は引き継ぎ
        rk[ith].SetMedium_D(rk[lastnum].GetMedium_D(k), k);
        for(int l=0;l<p_val;l++)
        {
            rk[ith].SetYin_J(rk[lastnum].GetYin_J(k, l), k, l);
        }
    }
    
    for(int i=0;i<cellnum;i++)
    {
        rk[i].SetTotalcellnumber(lastnum);//削除後の細胞数の更新
    }
    rk[lastnum].~RK();
}

void CELL::death_check(int ith)
{
    if(rk[ith].GetTime()-rk[ith].GetStationaryTime()>death_time)
    {
        removecell(ith);
        cout<<"delete cell"<<ith<<endl;
    }
}

void CELL::death_event()
{
    double probability=0.0,rara=0.0; // 確率（1%）
    int cellnumber=0;
    std::random_device rnd;     // 非決定的な乱数生成器
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<> rand01(0, 1);
    cellnumber= rk[0].GetTotalcellnumber();
    probability = death_rate;
    
    for(int i=0;i<cellnumber;i++)
    {
        rara = rand01(mt);
        if ( rara < probability ) {
            removecell(i);
            cout<<"delete_event cell"<<i<<endl;
        }
        cellnumber= rk[0].GetTotalcellnumber();
    }
}

void CELL::cell_death(int ith)//ith cellが死亡する
{
    removecell(ith);
}

void CELL::death_volume_check()
{
    double vol=0.0;
    int cellnumber=0;
    cellnumber= rk[0].GetTotalcellnumber();
    
    for(int i=0;i<cellnumber;i++)
    {
        vol = rk[i].GetVolume();
        if(vol<death_volume)
        {
            removecell(i);
            cout<<"delete["<<i<<"] #"<<rk[i].GetParasite_attacked_number()<<endl;
        }
        cellnumber= rk[0].GetTotalcellnumber();
    }
}

void CELL::division_time(int ith,double sigma)
{
    int time=0;
    int nowtime=0,divtime=0,newcellnumber=0;
    nowtime = rk[ith].GetTime1000();
    divtime = rk[ith].GetDiviTimeTime1000();
    time = nowtime%divtime;
    if(time==0)
    {
        cout<<"division_time="<<rk[ith].GetTime()<<endl;
        cout<<"nowtime1000="<<nowtime<<endl;
        cout<<"divtime1000="<<divtime<<endl;
        newcellmaker();
        newcellnumber=rk[0].GetTotalcellnumber()-1;
        for(int l=0;l<rk[0].Getprotein_x_number();l++)
        {
            //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
            rk[newcellnumber].SetResult(rk[ith].GetResult(l),l);
            //ガウスノイズを与えてcopy
            rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
            rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
            rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
            
            for(int k=0;k<rk[0].Getprotein_x_number();k++)
            {
                rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
            }
        }
        for(int i=0;i<rk[0].GetTotalcellnumber();i++)
        {
            rk[i].Show();
        }
    }
}

void CELL::division_volume(int ith,double sigma)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    uniform_real_distribution<> rand01(0, 1);
    
    int nowtime=0,divtime=0,newcellnumber=0,eliminate_cellnumber=0;
    double volume=0.0,initialvolume=0.0;
    nowtime = rk[ith].GetTime1000();
    divtime = rk[ith].GetDiviTimeTime1000();
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    newcellnumber = rk[0].GetTotalcellnumber();
    
    if(volume>=initialvolume)
    {
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<rk[ith].GetTime()<<endl;
            cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                //ガウスノイズを与えてcopy
                //rk[newcellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                
                for(int k=0;k<rk[0].Getprotein_x_number();k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            rk[ith].SetVolume(rk[ith].GetInitialVolume());//体積を初期化（半分）
            for(int i=0;i<rk[0].GetTotalcellnumber();i++)
            {
                rk[i].Show();
            }
        }
        else
        {
            cout<<"time["<<ith<<"]="<<rk[ith].GetTime()<<endl;
            cout<<"volume["<<ith<<"]="<<volume<<endl;
            eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
            cout<<eliminate_cellnumber<<"--->"<<ith<<endl;
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                rk[eliminate_cellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                for(int k=0;k<rk[0].Getprotein_x_number();k++)
                {
                    rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                }
                rk[eliminate_cellnumber].SetVolume(rk[ith].GetInitialVolume());
                rk[ith].SetVolume(rk[ith].GetInitialVolume());//体積を初期化（半分）
            }
            
        }
    }
}

void CELL::division_volume_mutation(int ith,double sigma)
{
    int nowtime=0,newcellnumber=0,eliminate_cellnumber=0,mother_celltype,choice=0,p_val=0;
    double volume=0.0,initialvolume=0.0,halfvolume=0.0;
    
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    
    if(volume>=initialvolume)
    {
        p_val =rk[ith].Getprotein_x_number();
        nowtime = rk[ith].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        random_device rnd;     // 非決定的な乱数生成器を生成
        mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
        uniform_int_distribution<> randProtein(0, ProteinVAL-1);
        uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
        uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
        uniform_real_distribution<> rand01(0, 1);
        uniform_int_distribution<> randchoice(1,mutation_rate);//mutationの選択用
        choice = randchoice(mt);
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        
        rk[ith].archive_evo_mini_stationary(nowtime,ith);//分裂直前を記録
        rk[ith].SetStationaryTime(nowtime);//分裂時間を更新
        rk[ith].SetStationaryGrowth(rk[ith].GetGrowthRate());//分裂時間を更新
        
        for(int j=0;j<ProteinVAL;j++)
        {
            rk[ith].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量を更新
            for(int l=0;l<ProteinVAL;l++)
            {
                rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
            }
        }
        
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                rk[newcellnumber].SetResult(rk[ith].gauss_rand(0.1,sigma),l);
                //rk[newcellnumber].SetResult(rk[ith].GetResult(l),l);
                //ガウスノイズを与えてcopy
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                
                for(int k=0;k<p_val;k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            rk[newcellnumber].converterJ();
            halfvolume =rk[ith].GetVolume()/2.0 -rk[ith].gauss_rand(0.05,sigma);
            rk[newcellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
            rk[newcellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[newcellnumber].tree_history.clear();
            rk[newcellnumber].tree_history = rk[ith].tree_history;
            rk[newcellnumber].ori_Sstr_J.clear();
            rk[newcellnumber].ori_Sstr_J.reserve(p_val*p_val);
            for(int i=0;i<p_val*p_val;i++)
            {
                rk[newcellnumber].ori_Sstr_J.push_back(rk[ith].ori_Sstr_J[i]);
            }
            for(int k=0;k<input_number;k++)
            {
                rk[newcellnumber].SetYin(0, k);//分裂時に子供のinput物質は0に
                rk[newcellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                rk[newcellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    rk[newcellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
        }
        else
        {
            //cout<<"time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
            while(eliminate_cellnumber==ith)
            {
                eliminate_cellnumber =randcell(mt);//親子が一致しないように振り直し続ける
            }
            if(eliminate_cellnumber==ith)
            {
                cout<<eliminate_cellnumber<<endl;
                cout<<ith<<endl;
            }
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            //cout<<"消す"<<rk[eliminate_cellnumber].GetCelltype()<<"----->"<<"増やす"<<rk[ith].GetCelltype()<<endl;
            
            for(int l=0;l<p_val;l++)
            {
                //rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                rk[eliminate_cellnumber].SetResult(rand01(mt),l);//子供の初期値はランダムに
                //rk[ith].SetResult(rand01(mt),l);//親の初期値もランダムに
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(0.05,sigma),l);//初期値を0付近に
                //rk[ith].SetResult(rk[ith].gauss_rand(0.05,sigma),l);//親の初期値を0付近に
                //rk[eliminate_cellnumber].SetResult(0.0,l);//初期値を0に
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                //rk[eliminate_cellnumber].SetProteinD(rk[ith].gabs_gauss_rand(rk[ith].GetProD(l),sigma),l);
                //rk[eliminate_cellnumber].SetThreshold_g(rk[ith].gauss_rand(rk[ith].GetThreshold_g(l),sigma), l);
                rk[eliminate_cellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                for(int k=0;k<p_val;k++)
                {
                    rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                }
                //cout<<"convert前="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
                rk[eliminate_cellnumber].converterJ();
                //rk[eliminate_cellnumber].SetSEL(l,rk[eliminate_cellnumber].GetELSize(l));//隣接リストの長さを取得
                //cout<<"convert後="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            }
            halfvolume =rk[ith].GetVolume()/2.0-rk[ith].gauss_rand(0.05,sigma);
            rk[eliminate_cellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
            rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
            //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[eliminate_cellnumber].SetDivision_number(0);//子供のカウントを0に
            
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //for(int i=0;i<rk[ith].tree_history.size();i++)
            //{rk[eliminate_cellnumber].tree_history.push_back(rk[ith].tree_history[i]);}//親のhistory引き継ぎ
            rk[eliminate_cellnumber].tree_history.clear();
            rk[eliminate_cellnumber].tree_history = rk[ith].tree_history;
            
            rk[eliminate_cellnumber].ori_Sstr_J.clear();
            rk[eliminate_cellnumber].ori_Sstr_J.reserve(p_val*p_val);
            for(int i=0;i<p_val*p_val;i++)
            {
                rk[eliminate_cellnumber].ori_Sstr_J.push_back(rk[ith].ori_Sstr_J[i]);
            }
            for(int k=0;k<input_number;k++)
            {
                rk[eliminate_cellnumber].SetYin(0, k);//分裂時に子供のinput物質は0に
                //rk[ith].SetYin(0, k);//分裂時に親のinput物質は0に
                rk[eliminate_cellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                rk[eliminate_cellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    rk[eliminate_cellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
        }
    }
}

void CELL::division_volume_mutation_network(int ith,double sigma,int time)
{
    int nowtime=0,newcellnumber=0,eliminate_cellnumber=0,child1,child2,child3,child4,j12,j32,j34,mother_celltype,
    p_val=0,count1=0,count2=0;
    double volume=0.0,initialvolume=0.0,halfvolume=0.0,rara1=0.0,rara2=0.0,probability=0.0;
    int count_limit=20;
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    probability = Host_mutaion_rate;
    
    if(volume>=initialvolume)
    {
        p_val =rk[ith].Getprotein_x_number();
        nowtime = rk[ith].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        random_device rnd;     // 非決定的な乱数生成器を生成
        mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
        uniform_int_distribution<> randProtein(0, ProteinVAL-1);
        uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
        uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
        uniform_real_distribution<> rand01(0, 1);
        uniform_real_distribution<> rand0th(0,_threshold);
        
        rara1 = rand01(mt);
        rara2 = rand01(mt);
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        rk[ith].archive_evo_mini_stationary(time,ith);//分裂直前を記録
        rk[ith].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
        rk[ith].SetStationaryGrowth(rk[ith].GetGrowthRate());//分裂時間を更新
        
        for(int j=0;j<ProteinVAL;j++)
        {
            rk[ith].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量を更新
            for(int l=0;l<ProteinVAL;l++)
            {
                rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
            }
        }
        
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            //cout<<"rara="<<rara<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                rk[newcellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                for(int k=0;k<p_val;k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[newcellnumber].SetStationaryJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            rk[newcellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[newcellnumber].converterJ();
            halfvolume =rk[ith].GetVolume()/2.0 ;
            rk[newcellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//親の体積を初期化（半分）
            
            rk[newcellnumber].SetCelltype(mother_celltype);//親ithのcelltype引き継ぎ
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[newcellnumber].tree_history.clear();
            rk[newcellnumber].tree_history = rk[ith].tree_history;
            rk[newcellnumber].SetTime(rk[ith].GetTime());//０秒を入れる
            rk[newcellnumber].SetStationaryTime(nowtime);//分裂時間を更新
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[newcellnumber].SetTime(0.0);//０秒を入れる
            
            for(int k=0;k<input_number;k++)
            {
                //rk[newcellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma), k);//分裂時に子供のinput物質は引き継ぎ
                //rk[newcellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[newcellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                //cout<<"MTd"<<rara1<<endl;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                while(rk[newcellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                child3 = randProtein(mt);
                while(rk[newcellnumber].GetJ(child3,child2) != 0 )//0のpathを探す
                {
                    child3 = randProtein(mt);
                }
                
                for(int l=0;l<input_number;l++)
                {
                    if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        child3 = randProtein_non_T(mt);//繋ぎ変える先からtarget-geneを除外する
                        child1 = randProtein_non_T(mt);
                    }
                }
                
                for(int i=0;i<target_number;i++)
                {
                    if(child3 == i||child1 == i)//繋ぎ変える先がtargetの時
                    {
                        count1=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[newcellnumber].GetJ(j,child2) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                        }
                        if(count1>=4)
                        {
                            child3 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                            child1 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                        }
                    }
                }
                
                j12 = rk[newcellnumber].GetJ(child1,child2);
                j32 = rk[newcellnumber].GetJ(child3,child2);
                
                if(j12 != j32)
                {
                    
                    //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<j12<<endl;
                    //cout<<"oldJ["<<child3<<"]["<<child2<<"]="<<j32<<endl;
                    //cout<<"time["<<ith<<"]="<<nowtime<<endl;
                    
                    rk[newcellnumber].SetJ(j12,child3,child2);//path切り替え 2--->1を2--->3に
                    rk[newcellnumber].SetJ(j32,child1,child2);//path切り替え 2--->3を2--->1に
                    rk[newcellnumber].converterJ();//隣接リストにmutationを反映
                    //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child1, child2)<<endl;
                    //cout<<"newJ["<<child3<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child3, child2)<<endl;
                    rk[newcellnumber].SetMutation_number(rk[newcellnumber].GetMutation_number()+1);
                    //変異体のカウントを増やす
                    //check_tree_history(ith,newcellnumber);
                    //rk[eliminate_cellnumber].Show();
                    //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
                }
            }
        }
        else
        {
            //cout<<"time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
            while(eliminate_cellnumber==ith)
            {
                eliminate_cellnumber =randcell(mt);//親子が一致しないように振り直し続ける
            }
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            //cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増やす"<<ith<<endl;
            
            for(int j=0;j<ProteinVAL;j++)
            {
                rk[eliminate_cellnumber].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量0に更新
            }
        
            for(int l=0;l<p_val;l++)
            {
                //rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[eliminate_cellnumber].SetResult(rand01(mt),l);//子供の初期値はランダムに
                rk[eliminate_cellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(0.1,sigma),l);//初期値を0付近に
                //rk[ith].SetResult(rk[ith].gauss_rand(0.05,sigma),l);//初期値を0付近に
                //rk[eliminate_cellnumber].SetResult(0.0,l);//初期値を0に
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                //rk[eliminate_cellnumber].SetProteinD(rk[ith].gabs_gauss_rand(rk[ith].GetProD(l),sigma),l);
                //rk[eliminate_cellnumber].SetThreshold_g(rk[ith].gauss_rand(rk[ith].GetThreshold_g(l),sigma), l);
                //rk[eliminate_cellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                for(int k=0;k<p_val;k++)
                {
                    rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[eliminate_cellnumber].SetStationaryJ(rk[ith].GetJ(l, k),l,k);
                }
            }
            rk[eliminate_cellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[eliminate_cellnumber].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
            rk[eliminate_cellnumber].SetStationaryGrowth( rk[ith].GetStationaryGrowth());//定常状態の成長速度を更新
            
            
            //cout<<"convert前="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            rk[eliminate_cellnumber].converterJ();
            //cout<<"convert後="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            
            //halfvolume =rk[ith].GetVolume()/2.0 -rk[ith].gauss_rand(rk[ith].GetVolume()/2.0,sigma);
            halfvolume =rk[ith].GetVolume()/2.0;
            rk[eliminate_cellnumber].SetVolume(halfvolume);
            //rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
            rk[ith].SetVolume(halfvolume);//体積を初期化（半分）
            rk[ith].SetTime(0.0);//０秒を入れる
            rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
            //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[eliminate_cellnumber].SetDivision_number(rk[ith].GetDivision_number() + 1);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[eliminate_cellnumber].SetTime(0.0);//０秒を入れる
            
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //for(int i=0;i<rk[ith].tree_history.size();i++)
            //{rk[eliminate_cellnumber].tree_history.push_back(rk[ith].tree_history[i]);}//親のhistory引き継ぎ
            //rk[eliminate_cellnumber].tree_history.clear();
            //rk[eliminate_cellnumber].tree_history = rk[ith].tree_history;
            
            //rk[eliminate_cellnumber].ori_Sstr_J.clear();
            //rk[eliminate_cellnumber].ori_Sstr_J.reserve(p_val*p_val);
            for(int i=0;i<p_val*p_val;i++)
            {
                //rk[eliminate_cellnumber].ori_Sstr_J.push_back(rk[ith].ori_Sstr_J[i]);
            }
            for(int k=0;k<input_number;k++)
            {
                //rk[eliminate_cellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
                //rk[eliminate_cellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[eliminate_cellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[eliminate_cellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[eliminate_cellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            
            //cout<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            //cout<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<ith<<"].tree_history="<<rk[ith].tree_history<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //cout<<"rk["<<ith<<"].ori_Sstr="<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].ori_Sstr_J="<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                
                cout<<"time["<<ith<<"]="<<nowtime<<endl;
                cout<<"volume["<<ith<<"]="<<volume<<endl;
                int count1=0,count2=0;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
                child3 = randProtein(mt);
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
                
                while(rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);//inputは変化しない
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                    count1++;
                    if(count1>count_limit){
                        //while文の上限を設定
                        break;
                    }
                }
                
                while(rk[eliminate_cellnumber].GetJ(child3,child4) != 0 ||
                      rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                      rk[eliminate_cellnumber].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
                {
                    child3 = randProtein(mt);//inputは変化しない
                    child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                    if(count1>count_limit){
                        break;
                    }
                }
                
                //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
                //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
                
                for(int l=0;l<input_number;l++)
                {
                    if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        //cout<<"Input_child4="<<child4<<endl;
                        //cout<<"old_child3="<<child3<<endl;
                        child3 = randProtein_non_T(mt);
                        while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                              rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                            if(count1>count_limit){
                                break;
                            }
                        }
                        //cout<<"Input_child4="<<child4<<endl;
                        //cout<<"new_child3="<<child3<<endl;
                        break;
                    }
                }
                
                if(on==1)//使わない
                {
                    for(int l=0;l<Parasite_genome_size;l++)//３種
                    {
                        if(child2 == l+Parasite_target_genome||child1 == l)//attcked-geneとtargetを必ず繋いでるpathだった時
                        {
                            while(rk[eliminate_cellnumber].GetSEL(child1)>=MaxPathNum ||
                                  rk[eliminate_cellnumber].GetJ(child1,child2) != 0 )//別の0のpathを探す
                            {
                                child1 = randProtein_non_T(mt);//targetでないところから探す
                                if(count1>count_limit){
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
                
                for(int t=0;t<target_number;t++)
                {
                    if(child3 == t)//繋ぎ変える先がtargetの時
                    {
                        count1=0;count2=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[eliminate_cellnumber].GetJ(j,child4) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                            else if(rk[eliminate_cellnumber].GetJ(j,child4) == -1)
                            {
                                count2 = count2 +1;
                            }
                        }
                        //if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                        if(count1>=MaxTargetPath)
                        {
                            child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                            while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                                  rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                            {
                                child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                                if(count1>count_limit){
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
                
                j12 = rk[eliminate_cellnumber].GetJ(child1,child2);//-1 or 1
                j34 = rk[eliminate_cellnumber].GetJ(child3,child4);//0
                
                if(j12 != j34 && count1<count_limit)
                {
                    cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                    cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                    
                    rk[eliminate_cellnumber].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                    rk[eliminate_cellnumber].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                    
                    cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                    cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                    rk[eliminate_cellnumber].converterJ();//隣接リストにmutationを反映
                    rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                    //変異回数のカウントを増やす
                    //check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
                }
            }
        }
    }
}


void CELL::division_volume_mutation_network_continuas(int ith,double sigma,int time)
{
    int nowtime=0,newcellnumber=0,eliminate_cellnumber=0,child1,child2,child3,child4,j12,j32,j34,mother_celltype,
    p_val=0,count1=0,count2=0;
    double volume=0.0,initialvolume=0.0,halfvolume=0.0,rara1=0.0,rara2=0.0,probability=0.0;
    
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    probability = Host_mutaion_rate;
    
    if(volume>=initialvolume)
    {
        p_val =rk[ith].Getprotein_x_number();
        nowtime = rk[ith].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        random_device rnd;     // 非決定的な乱数生成器を生成
        mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
        uniform_int_distribution<> randProtein(0, ProteinVAL-1);
        uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
        uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
        uniform_real_distribution<> rand01(0, 1);
        uniform_real_distribution<> rand0th(0,_threshold);
        
        rara1 = rand01(mt);
        rara2 = rand01(mt);
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        rk[ith].archive_evo_mini_stationary(time,ith);//分裂直前を記録
        rk[ith].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
        rk[ith].SetStationaryGrowth(rk[ith].GetGrowthRate());//分裂時間を更新
        
        for(int j=0;j<ProteinVAL;j++)
        {
            rk[ith].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量を更新
            for(int l=0;l<ProteinVAL;l++)
            {
                rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
            }
        }
        
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            //cout<<"rara="<<rara<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                rk[newcellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                for(int k=0;k<p_val;k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[newcellnumber].SetStationaryJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            rk[newcellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[newcellnumber].converterJ();
            halfvolume =rk[ith].GetVolume()/2.0 ;
            rk[newcellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//親の体積を初期化（半分）
            
            rk[newcellnumber].SetCelltype(mother_celltype);//親ithのcelltype引き継ぎ
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[newcellnumber].tree_history.clear();
            rk[newcellnumber].tree_history = rk[ith].tree_history;
            rk[newcellnumber].SetTime(rk[ith].GetTime());//０秒を入れる
            rk[newcellnumber].SetStationaryTime(nowtime);//分裂時間を更新
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[newcellnumber].SetTime(0.0);//０秒を入れる
            
            for(int k=0;k<input_number;k++)
            {
                //rk[newcellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma), k);//分裂時に子供のinput物質は引き継ぎ
                //rk[newcellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[newcellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                //cout<<"MTd"<<rara1<<endl;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                while(rk[newcellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                child3 = randProtein(mt);
                while(rk[newcellnumber].GetJ(child3,child2) != 0 )//0のpathを探す
                {
                    child3 = randProtein(mt);
                }
                
                for(int l=0;l<input_number;l++)
                {
                    if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        child3 = randProtein_non_T(mt);//繋ぎ変える先からtarget-geneを除外する
                        child1 = randProtein_non_T(mt);
                    }
                }
                
                for(int i=0;i<target_number;i++)
                {
                    if(child3 == i||child1 == i)//繋ぎ変える先がtargetの時
                    {
                        count1=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[newcellnumber].GetJ(j,child2) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                        }
                        if(count1>=4)
                        {
                            child3 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                            child1 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                        }
                    }
                }
                
                j12 = rk[newcellnumber].GetJ(child1,child2);
                j32 = rk[newcellnumber].GetJ(child3,child2);
                
                if(j12 != j32)
                {
                    
                    //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<j12<<endl;
                    //cout<<"oldJ["<<child3<<"]["<<child2<<"]="<<j32<<endl;
                    //cout<<"time["<<ith<<"]="<<nowtime<<endl;
                    
                    rk[newcellnumber].SetJ(j12,child3,child2);//path切り替え 2--->1を2--->3に
                    rk[newcellnumber].SetJ(j32,child1,child2);//path切り替え 2--->3を2--->1に
                    rk[newcellnumber].converterJ();//隣接リストにmutationを反映
                    //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child1, child2)<<endl;
                    //cout<<"newJ["<<child3<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child3, child2)<<endl;
                    rk[newcellnumber].SetMutation_number(rk[newcellnumber].GetMutation_number()+1);
                    //変異体のカウントを増やす
                    //check_tree_history(ith,newcellnumber);
                    //rk[eliminate_cellnumber].Show();
                    //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
                }
            }
        }
        else
        {
            //cout<<"time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
            while(eliminate_cellnumber==ith)
            {
                eliminate_cellnumber =randcell(mt);//親子が一致しないように振り直し続ける
            }
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            //cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増やす"<<ith<<endl;
            
            for(int l=0;l<p_val;l++)
            {
                //rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[eliminate_cellnumber].SetResult(rand01(mt),l);//子供の初期値はランダムに
                rk[eliminate_cellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                
                for(int k=0;k<p_val;k++)
                {
                    rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[eliminate_cellnumber].SetStationaryJ(rk[ith].GetJ(l, k),l,k);
                }
            }
            rk[eliminate_cellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[eliminate_cellnumber].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
            rk[eliminate_cellnumber].SetStationaryGrowth( rk[ith].GetStationaryGrowth());//定常状態の成長速度を更新
            for(int j=0;j<ProteinVAL;j++)
            {
                rk[eliminate_cellnumber].SetStationaryResult(rk[ith].GetStationaryResult(j),j);//分裂時の発現量0に更新
            }
            
            //cout<<"convert前="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            rk[eliminate_cellnumber].converterJ();
            //cout<<"convert後="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            
            //halfvolume =rk[ith].GetVolume()/2.0 -rk[ith].gauss_rand(rk[ith].GetVolume()/2.0,sigma);
            halfvolume =rk[ith].GetVolume()/2.0;
            rk[eliminate_cellnumber].SetVolume(halfvolume);
            //rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
            rk[ith].SetVolume(halfvolume);//体積を初期化（半分）
            rk[ith].SetTime(0.0);//０秒を入れる
            rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
            //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[eliminate_cellnumber].SetDivision_number(rk[ith].GetDivision_number() + 1);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[eliminate_cellnumber].SetTime(0.0);//０秒を入れる
            
            //cout<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            //cout<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<ith<<"].tree_history="<<rk[ith].tree_history<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //cout<<"rk["<<ith<<"].ori_Sstr="<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].ori_Sstr_J="<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                
                cout<<"time["<<ith<<"]="<<nowtime<<endl;
                cout<<"volume["<<ith<<"]="<<volume<<endl;
                
                for(int l=0;l<p_val;l++)
                {
                    for(int k=0;k<p_val;k++)
                    {
                        rk[eliminate_cellnumber].SetJ(rk[eliminate_cellnumber].gauss_rand(rk[eliminate_cellnumber].GetJ(l, k),sigma), l, k);//変異が加わる
                    }
                }
                
                for(int l=0;l<input_number;l++)
                {
                    for(int t=0;t<target_number;t++)
                    {
                        rk[eliminate_cellnumber].SetJ(0,t,l);
                    }
                }
                
                for(int t=0;t<target_number;t++)
                {
                    for(int l=0;l<p_val;l++)
                    {
                        rk[eliminate_cellnumber].SetJ(0,l,t);
                    }
                }
                rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            }
        }
    }
}


void CELL::division_volume_mutation_network_sigma(int ith,double sigma,int time)
{
    int nowtime=0,newcellnumber=0,eliminate_cellnumber=0,child1,child2,child3,child4,j12,j32,j34,mother_celltype,
    p_val=0,count1=0,count2=0;
    double volume=0.0,initialvolume=0.0,halfvolume=0.0,rara1=0.0,rara2=0.0,probability=0.0;
    
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    probability = Host_mutaion_rate;
    
    if(volume>=initialvolume)
    {
        p_val =rk[ith].Getprotein_x_number();
        nowtime = rk[ith].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        random_device rnd;     // 非決定的な乱数生成器を生成
        mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
        uniform_int_distribution<> randProtein(0, ProteinVAL-1);
        uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
        uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
        uniform_real_distribution<> rand01(0, 1);
        uniform_real_distribution<> rand0th(0,_threshold);
        
        rara1 = rand01(mt);
        rara2 = rand01(mt);
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        rk[ith].archive_evo_mini_stationary(time,ith);//分裂直前を記録
        rk[ith].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
        rk[ith].SetStationaryGrowth(rk[ith].GetGrowthRate());//分裂時間を更新
        
        for(int j=0;j<ProteinVAL;j++)
        {
            rk[ith].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量を更新
            for(int l=0;l<ProteinVAL;l++)
            {
                rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
            }
        }
        
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            //cout<<"rara="<<rara<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                rk[newcellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                
                for(int k=0;k<p_val;k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[newcellnumber].SetStationaryJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            
            rk[newcellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[newcellnumber].converterJ();
            halfvolume =rk[ith].GetVolume()/2.0 ;
            rk[newcellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//親の体積を初期化（半分）
            
            rk[newcellnumber].SetCelltype(mother_celltype);//親ithのcelltype引き継ぎ
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[newcellnumber].tree_history.clear();
            rk[newcellnumber].tree_history = rk[ith].tree_history;
            rk[newcellnumber].SetTime(rk[ith].GetTime());//０秒を入れる
            rk[newcellnumber].SetStationaryTime(nowtime);//分裂時間を更新
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[newcellnumber].SetTime(0.0);//０秒を入れる
            
            for(int k=0;k<input_number;k++)
            {
                //rk[newcellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma), k);//分裂時に子供のinput物質は引き継ぎ
                //rk[newcellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[newcellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                //cout<<"MTd"<<rara1<<endl;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                while(rk[newcellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                child3 = randProtein(mt);
                while(rk[newcellnumber].GetJ(child3,child2) != 0 )//0のpathを探す
                {
                    child3 = randProtein(mt);
                }
                
                for(int l=0;l<input_number;l++)
                {
                    if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        child3 = randProtein_non_T(mt);//繋ぎ変える先からtarget-geneを除外する
                        child1 = randProtein_non_T(mt);
                    }
                }
                
                for(int i=0;i<target_number;i++)
                {
                    if(child3 == i||child1 == i)//繋ぎ変える先がtargetの時
                    {
                        count1=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[newcellnumber].GetJ(j,child2) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                        }
                        if(count1>=4)
                        {
                            child3 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                            child1 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                        }
                    }
                }
                
                j12 = rk[newcellnumber].GetJ(child1,child2);
                j32 = rk[newcellnumber].GetJ(child3,child2);
                
                if(j12 != j32)
                {
                    
                    //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<j12<<endl;
                    //cout<<"oldJ["<<child3<<"]["<<child2<<"]="<<j32<<endl;
                    //cout<<"time["<<ith<<"]="<<nowtime<<endl;
                    
                    rk[newcellnumber].SetJ(j12,child3,child2);//path切り替え 2--->1を2--->3に
                    rk[newcellnumber].SetJ(j32,child1,child2);//path切り替え 2--->3を2--->1に
                    rk[newcellnumber].converterJ();//隣接リストにmutationを反映
                    //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child1, child2)<<endl;
                    //cout<<"newJ["<<child3<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child3, child2)<<endl;
                    rk[newcellnumber].SetMutation_number(rk[newcellnumber].GetMutation_number()+1);
                    //変異体のカウントを増やす
                    //check_tree_history(ith,newcellnumber);
                    //rk[eliminate_cellnumber].Show();
                    //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
                }
                rk[newcellnumber].SetSigma_Str(rk[ith].gabs_gauss_rand(rk[ith].GetSigma_Str(),sigma));//変異が加わる
            }
        }
            else
            {
                //cout<<"time["<<ith<<"]="<<nowtime<<endl;
                //cout<<"volume["<<ith<<"]="<<volume<<endl;
                eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
                while(eliminate_cellnumber==ith)
                {
                    eliminate_cellnumber =randcell(mt);//親子が一致しないように振り直し続ける
                }
                mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
                //cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増やす"<<ith<<endl;
                
                for(int l=0;l<p_val;l++)
                {
                    //rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                    rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                    rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                    //rk[eliminate_cellnumber].SetResult(rand01(mt),l);//子供の初期値はランダムに
                    rk[eliminate_cellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                    rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                    //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(0.1,sigma),l);//初期値を0付近に
                    //rk[ith].SetResult(rk[ith].gauss_rand(0.05,sigma),l);//初期値を0付近に
                    //rk[eliminate_cellnumber].SetResult(0.0,l);//初期値を0に
                    //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                    //rk[eliminate_cellnumber].SetProteinD(rk[ith].gabs_gauss_rand(rk[ith].GetProD(l),sigma),l);
                    //rk[eliminate_cellnumber].SetThreshold_g(rk[ith].gauss_rand(rk[ith].GetThreshold_g(l),sigma), l);
                    //rk[eliminate_cellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                    for(int k=0;k<p_val;k++)
                    {
                        rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                        rk[eliminate_cellnumber].SetStationaryJ(rk[ith].GetJ(l, k),l,k);
                    }
                }
                rk[eliminate_cellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
                rk[eliminate_cellnumber].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
                rk[eliminate_cellnumber].SetStationaryGrowth( rk[ith].GetStationaryGrowth());//定常状態の成長速度を更新
                for(int j=0;j<ProteinVAL;j++)
                {
                    rk[eliminate_cellnumber].SetStationaryResult(rk[ith].GetStationaryResult(j),j);//分裂時の発現量0に更新
                }
                
                //cout<<"convert前="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
                rk[eliminate_cellnumber].converterJ();
                //cout<<"convert後="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
                
                //halfvolume =rk[ith].GetVolume()/2.0 -rk[ith].gauss_rand(rk[ith].GetVolume()/2.0,sigma);
                halfvolume =rk[ith].GetVolume()/2.0;
                rk[eliminate_cellnumber].SetVolume(halfvolume);
                //rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
                rk[ith].SetVolume(halfvolume);//体積を初期化（半分）
                rk[ith].SetTime(0.0);//０秒を入れる
                rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
                //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
                rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
                rk[eliminate_cellnumber].SetDivision_number(rk[ith].GetDivision_number() + 1);//子供のカウントを0に
                rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
                rk[eliminate_cellnumber].SetSigma_Str(rk[ith].GetSigma_Str());//sigmaを引き継ぐ
                rk[eliminate_cellnumber].SetTime(0.0);//０秒を入れる
                
                //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
                //for(int i=0;i<rk[ith].tree_history.size();i++)
                //{rk[eliminate_cellnumber].tree_history.push_back(rk[ith].tree_history[i]);}//親のhistory引き継ぎ
                //rk[eliminate_cellnumber].tree_history.clear();
                //rk[eliminate_cellnumber].tree_history = rk[ith].tree_history;
                
                //rk[eliminate_cellnumber].ori_Sstr_J.clear();
                //rk[eliminate_cellnumber].ori_Sstr_J.reserve(p_val*p_val);
                for(int i=0;i<p_val*p_val;i++)
                {
                    //rk[eliminate_cellnumber].ori_Sstr_J.push_back(rk[ith].ori_Sstr_J[i]);
                }
                for(int k=0;k<input_number;k++)
                {
                    //rk[eliminate_cellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
                    //rk[eliminate_cellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                    //rk[eliminate_cellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                    //rk[eliminate_cellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                    for(int l=0;l<p_val;l++)
                    {
                        //rk[eliminate_cellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                    }
                }
                if (rara2<Host_Sigma_mutaion_rate)
                {
                    rk[eliminate_cellnumber].SetSigma_Str(rk[ith].gabs_gauss_rand(rk[ith].GetSigma_Str(),sigma));//変異が加わる
                }
                if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
                {
                    
                    cout<<"time["<<ith<<"]="<<nowtime<<endl;
                    cout<<"volume["<<ith<<"]="<<volume<<endl;
                    child1 = randProtein(mt);
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
                    child3 = randProtein(mt);
                    child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
                    
                    while(rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                    {
                        child1 = randProtein(mt);//inputは変化しない
                        child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                    }
                    
                    while(rk[eliminate_cellnumber].GetJ(child3,child4) != 0 ||
                          rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                          rk[eliminate_cellnumber].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
                    {
                        child3 = randProtein(mt);//inputは変化しない
                        child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                    }
                    
                    //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
                    //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
                    
                    for(int l=0;l<input_number;l++)
                    {
                        if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                        {
                            //cout<<"Input_child4="<<child4<<endl;
                            //cout<<"old_child3="<<child3<<endl;
                            child3 = randProtein_non_T(mt);
                            while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                                  rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                            {
                                child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                            }
                            //cout<<"Input_child4="<<child4<<endl;
                            //cout<<"new_child3="<<child3<<endl;
                            break;
                        }
                    }
                    
                    if(on==1)//使わない
                    {
                        for(int l=0;l<Parasite_genome_size;l++)//３種
                        {
                            if(child2 == l+Parasite_target_genome||child1 == l)//attcked-geneとtargetを必ず繋いでるpathだった時
                            {
                                while(rk[eliminate_cellnumber].GetSEL(child1)>=MaxPathNum ||
                                      rk[eliminate_cellnumber].GetJ(child1,child2) != 0 )//別の0のpathを探す
                                {
                                    child1 = randProtein_non_T(mt);//targetでないところから探す
                                }
                                break;
                            }
                        }
                    }
                    
                    for(int t=0;t<target_number;t++)
                    {
                        if(child3 == t)//繋ぎ変える先がtargetの時
                        {
                            count1=0;count2=0;
                            for(int j=0;j<target_number;j++)
                            {
                                if(rk[eliminate_cellnumber].GetJ(j,child4) == 1)
                                    //すでに2からtargetがactivateされている時
                                {
                                    count1 = count1 + 1;
                                }
                                else if(rk[eliminate_cellnumber].GetJ(j,child4) == -1)
                                {
                                    count2 = count2 +1;
                                }
                            }
                            //if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                            if(count1>=MaxTargetPath)
                            {
                                child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                                while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                                      rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                                {
                                    child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                                }
                            }
                            break;
                        }
                    }
                    
                    j12 = rk[eliminate_cellnumber].GetJ(child1,child2);//-1 or 1
                    j34 = rk[eliminate_cellnumber].GetJ(child3,child4);//0
                    
                    if(j12 != j34)
                    {
                        cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                        cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                        
                        rk[eliminate_cellnumber].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                        rk[eliminate_cellnumber].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                        
                        cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                        cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                        rk[eliminate_cellnumber].converterJ();//隣接リストにmutationを反映
                        rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                        //変異回数のカウントを増やす
                        //check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
                    }
                }
                
            }
        }
}

void CELL::division_volume_mutation_network_couple(int ith,double sigma,int time)
{
    int nowtime=0,newcellnumber=0,eliminate_cellnumber=0,child1,child2,child3,child4,j12,j32,j34,mother_celltype,
    p_val=0,count1=0,count2=0;
    double volume=0.0,initialvolume=0.0,halfvolume=0.0,rara1=0.0,rara2=0.0,probability=0.0;
    
    volume = rk[ith].GetVolume();
    initialvolume = rk[ith].GetInitialVolume();
    initialvolume = initialvolume*2;//体積２倍で分裂
    probability = Host_mutaion_rate;
    
    if(volume>=initialvolume)
    {
        p_val =rk[ith].Getprotein_x_number();
        nowtime = rk[ith].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        random_device rnd;     // 非決定的な乱数生成器を生成
        mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
        uniform_int_distribution<> randProtein(0, ProteinVAL-1);
        uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
        uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
        uniform_real_distribution<> rand01(0, 1);
        uniform_real_distribution<> rand0th(0,_threshold);
        
        rara1 = rand01(mt);
        rara2 = rand01(mt);
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        rk[ith].archive_evo_mini_stationary(time,ith);//分裂直前を記録
        rk[ith].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
        rk[ith].SetStationaryGrowth(rk[ith].GetGrowthRate());//分裂時間を更新
        
        for(int j=0;j<ProteinVAL;j++)
        {
            rk[ith].SetStationaryResult(rk[ith].GetResult(j),j);//分裂時の発現量を更新
            for(int l=0;l<ProteinVAL;l++)
            {
                rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
            }
        }
        
        if(newcellnumber<MAXCELLNUMBER)//最大cell数を定めておく
        {
            cout<<"division_time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            cout<<newcellnumber<<endl;
            //cout<<"rara="<<rara<<endl;
            newcellmaker();
            newcellnumber=rk[0].GetTotalcellnumber()-1;
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            for(int l=0;l<rk[0].Getprotein_x_number();l++)
            {
                //rk[newcellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                rk[newcellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[newcellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[newcellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[newcellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                
                
                for(int k=0;k<p_val;k++)
                {
                    rk[newcellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[newcellnumber].SetStationaryJ(rk[ith].GetJ(l, k), l, k);
                }
            }
            rk[newcellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[newcellnumber].converterJ();
            halfvolume =rk[ith].GetVolume()/2.0 ;
            rk[newcellnumber].SetVolume(halfvolume);
            rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//親の体積を初期化（半分）
            
            rk[newcellnumber].SetCelltype(mother_celltype);//親ithのcelltype引き継ぎ
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[newcellnumber].tree_history.clear();
            rk[newcellnumber].tree_history = rk[ith].tree_history;
            rk[newcellnumber].SetTime(rk[ith].GetTime());//０秒を入れる
            rk[newcellnumber].SetStationaryTime(nowtime);//分裂時間を更新
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[newcellnumber].SetDivision_number(0);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[newcellnumber].SetTime(0.0);//０秒を入れる
            
            for(int k=0;k<input_number;k++)
            {
                //rk[newcellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma), k);//分裂時に子供のinput物質は引き継ぎ
                //rk[newcellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[newcellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[newcellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                //cout<<"MTd"<<rara1<<endl;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                while(rk[newcellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                child3 = randProtein(mt);
                while(rk[newcellnumber].GetJ(child3,child2) == 0 )//0のpathを探す
                {
                    child3 = randProtein(mt);
                }
                
                for(int l=0;l<input_number;l++)
                {
                    if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        child3 = randProtein_non_T(mt);//繋ぎ変える先からtarget-geneを除外する
                        child1 = randProtein_non_T(mt);
                    }
                }
                
                for(int i=0;i<target_number;i++)
                {
                    if(child3 == i||child1 == i)//繋ぎ変える先がtargetの時
                    {
                        count1=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[newcellnumber].GetJ(j,child2) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                        }
                        if(count1>=4)
                        {
                            child3 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                            child1 = randProtein_non_T(mt);//target-geneへのpathの本数を規制する
                        }
                    }
                }
                
                j12 = rk[newcellnumber].GetJ(child1,child2);
                j32 = rk[newcellnumber].GetJ(child3,child2);
                
                if(j12 != j32)
                {
                    
                    //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<j12<<endl;
                    //cout<<"oldJ["<<child3<<"]["<<child2<<"]="<<j32<<endl;
                    //cout<<"time["<<ith<<"]="<<nowtime<<endl;
                    
                    rk[newcellnumber].SetJ(j12,child3,child2);//path切り替え 2--->1を2--->3に
                    rk[newcellnumber].SetJ(j32,child1,child2);//path切り替え 2--->3を2--->1に
                    rk[newcellnumber].converterJ();//隣接リストにmutationを反映
                    //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child1, child2)<<endl;
                    //cout<<"newJ["<<child3<<"]["<<child2<<"]="<<rk[eliminate_cellnumber].GetJ(child3, child2)<<endl;
                    rk[newcellnumber].SetMutation_number(rk[newcellnumber].GetMutation_number()+1);
                    //変異体のカウントを増やす
                    //check_tree_history(ith,newcellnumber);
                    //rk[eliminate_cellnumber].Show();
                    //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
                }
                rk[newcellnumber].SetSigma_Str(rk[ith].gabs_gauss_rand(rk[ith].GetSigma_Str(),sigma));//変異が加わる
            }
        }
        else
        {
            //cout<<"time["<<ith<<"]="<<nowtime<<endl;
            //cout<<"volume["<<ith<<"]="<<volume<<endl;
            eliminate_cellnumber = randcell(mt);//削除したcellをdaughtercellにする　ithがmother
            while(eliminate_cellnumber==ith)
            {
                eliminate_cellnumber =randcell(mt);//親子が一致しないように振り直し続ける
            }
            mother_celltype = rk[ith].GetCelltype();//親のcelltypeを引き継ぎ
            //cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増やす"<<ith<<endl;
            
            for(int l=0;l<p_val;l++)
            {
                //rk[eliminate_cellnumber].SetResult(rk[ith].GetResult(l),l);
                rk[eliminate_cellnumber].SetProteinD(rk[ith].GetProD(l), l);
                rk[eliminate_cellnumber].SetThreshold_g(rk[ith].GetThreshold_g(l), l);
                //rk[eliminate_cellnumber].SetResult(rand01(mt),l);//子供の初期値はランダムに
                rk[eliminate_cellnumber].SetResult(rand0th(mt),l);//閾値以下で乱数
                rk[ith].SetResult(rand0th(mt),l);//親の初期値もランダムに
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(0.1,sigma),l);//初期値を0付近に
                //rk[ith].SetResult(rk[ith].gauss_rand(0.05,sigma),l);//初期値を0付近に
                //rk[eliminate_cellnumber].SetResult(0.0,l);//初期値を0に
                //rk[eliminate_cellnumber].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
                //rk[eliminate_cellnumber].SetProteinD(rk[ith].gabs_gauss_rand(rk[ith].GetProD(l),sigma),l);
                //rk[eliminate_cellnumber].SetThreshold_g(rk[ith].gauss_rand(rk[ith].GetThreshold_g(l),sigma), l);
                //rk[eliminate_cellnumber].SetMeanfieldval(rk[ith].GetMeanfieldval(l),l);
                for(int k=0;k<p_val;k++)
                {
                    rk[eliminate_cellnumber].SetJ(rk[ith].GetJ(l, k), l, k);
                    rk[eliminate_cellnumber].SetStationaryJ(rk[ith].GetJ(l, k),l,k);
                }
            }
            rk[eliminate_cellnumber].SetSigma_Str(rk[ith].GetSigma_Str());
            rk[eliminate_cellnumber].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
            rk[eliminate_cellnumber].SetStationaryGrowth( rk[ith].GetStationaryGrowth());//定常状態の成長速度を更新
            for(int j=0;j<ProteinVAL;j++)
            {
                rk[eliminate_cellnumber].SetStationaryResult(rk[ith].GetStationaryResult(j),j);//分裂時の発現量0に更新
            }
            
            //cout<<"convert前="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            rk[eliminate_cellnumber].converterJ();
            //cout<<"convert後="<<rk[eliminate_cellnumber].GetSEL(l)<<endl;
            
            //halfvolume =rk[ith].GetVolume()/2.0 -rk[ith].gauss_rand(rk[ith].GetVolume()/2.0,sigma);
            halfvolume =rk[ith].GetVolume()/2.0;
            rk[eliminate_cellnumber].SetVolume(halfvolume);
            //rk[ith].SetVolume(rk[ith].GetVolume()-halfvolume);//体積を初期化（半分）
            rk[ith].SetVolume(halfvolume);//体積を初期化（半分）
            rk[ith].SetTime(0.0);//０秒を入れる
            rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
            //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
            rk[ith].SetDivision_number(rk[ith].GetDivision_number() + 1);//親のカウントを増やす
            rk[eliminate_cellnumber].SetDivision_number(rk[ith].GetDivision_number() + 1);//子供のカウントを0に
            rk[ith].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number());//初期値からの変異数を記録
            rk[eliminate_cellnumber].SetTime(0.0);//０秒を入れる
            
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //for(int i=0;i<rk[ith].tree_history.size();i++)
            //{rk[eliminate_cellnumber].tree_history.push_back(rk[ith].tree_history[i]);}//親のhistory引き継ぎ
            //rk[eliminate_cellnumber].tree_history.clear();
            //rk[eliminate_cellnumber].tree_history = rk[ith].tree_history;
            
            //rk[eliminate_cellnumber].ori_Sstr_J.clear();
            //rk[eliminate_cellnumber].ori_Sstr_J.reserve(p_val*p_val);
            for(int i=0;i<p_val*p_val;i++)
            {
                //rk[eliminate_cellnumber].ori_Sstr_J.push_back(rk[ith].ori_Sstr_J[i]);
            }
            for(int k=0;k<input_number;k++)
            {
                //rk[eliminate_cellnumber].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
                //rk[eliminate_cellnumber].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
                //rk[eliminate_cellnumber].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
                //rk[eliminate_cellnumber].SetMedium_D(rk[ith].GetMedium_D(k), k);
                for(int l=0;l<p_val;l++)
                {
                    //rk[eliminate_cellnumber].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
                }
            }
            
            //cout<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            //cout<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<ith<<"].tree_history="<<rk[ith].tree_history<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].tree_history="<<rk[eliminate_cellnumber].tree_history<<endl;
            //cout<<"rk["<<ith<<"].ori_Sstr="<<rk[ith].ori_Sstr_J.size()<<endl;
            //cout<<"rk["<<eliminate_cellnumber<<"].ori_Sstr_J="<<rk[eliminate_cellnumber].ori_Sstr_J.size()<<endl;
            
            if(rara1<probability)//mutationが起こる確率を設定 pathに制限あり
            {
                
                cout<<"time["<<ith<<"]="<<nowtime<<endl;
                cout<<"volume["<<ith<<"]="<<volume<<endl;
                child1 = randProtein(mt);
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
                child3 = randProtein(mt);
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
                
                while(rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
                {
                    child1 = randProtein(mt);//inputは変化しない
                    child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                
                while(rk[eliminate_cellnumber].GetJ(child3,child4) == 0 ||
                      rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                      rk[eliminate_cellnumber].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
                {
                    child3 = randProtein(mt);//inputは変化しない
                    child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
                }
                
                //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
                //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
                
                for(int l=0;l<input_number;l++)
                {
                    if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                    {
                        //cout<<"Input_child4="<<child4<<endl;
                        //cout<<"old_child3="<<child3<<endl;
                        child3 = randProtein_non_T(mt);
                        while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                              rk[eliminate_cellnumber].GetJ(child3,child4) == 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                        }
                        //cout<<"Input_child4="<<child4<<endl;
                        //cout<<"new_child3="<<child3<<endl;
                        break;
                    }
                }
                
                if(on==1)//使わない
                {
                    for(int l=0;l<Parasite_genome_size;l++)//３種
                    {
                        if(child2 == l+Parasite_target_genome||child1 == l)//attcked-geneとtargetを必ず繋いでるpathだった時
                        {
                            while(rk[eliminate_cellnumber].GetSEL(child1)>=MaxPathNum ||
                                  rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//別の0のpathを探す
                            {
                                child1 = randProtein_non_T(mt);//targetでないところから探す
                            }
                            break;
                        }
                    }
                }
                
                for(int t=0;t<target_number;t++)
                {
                    if(child3 == t)//繋ぎ変える先がtargetの時
                    {
                        count1=0;count2=0;
                        for(int j=0;j<target_number;j++)
                        {
                            if(rk[eliminate_cellnumber].GetJ(j,child4) == 1)
                                //すでに2からtargetがactivateされている時
                            {
                                count1 = count1 + 1;
                            }
                            else if(rk[eliminate_cellnumber].GetJ(j,child4) == -1)
                            {
                                count2 = count2 +1;
                            }
                        }
                        //if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                        if(count1>=MaxTargetPath)
                        {
                            child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                            while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                                  rk[eliminate_cellnumber].GetJ(child3,child4) == 0 )//別の0のpathを探す
                            {
                                child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                            }
                        }
                        break;
                    }
                }
                
                j12 = rk[eliminate_cellnumber].GetJ(child1,child2);//-1 or 1
                j34 = rk[eliminate_cellnumber].GetJ(child3,child4);//0
                
                if(j12 != j34)
                {
                    cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                    cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                    
                    rk[eliminate_cellnumber].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                    rk[eliminate_cellnumber].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                    
                    cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                    cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                    rk[eliminate_cellnumber].converterJ();//隣接リストにmutationを反映
                    rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                    //変異回数のカウントを増やす
                    //check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
                }
            }
            if (rara2<Host_Sigma_mutaion_rate)
            {
                rk[eliminate_cellnumber].SetSigma_Str(rk[ith].gabs_gauss_rand(rk[ith].GetSigma_Str(),sigma));//変異が加わる
            }
        }
    }
}

void CELL::division_fitness_mutation(double sigma,int time)
{
    {
        int nowtime=0,newcellnumber=0,p_val=0,ith;
        
        p_val =rk[0].Getprotein_x_number();
        nowtime = rk[0].GetTime();
        newcellnumber = rk[0].GetTotalcellnumber();
        
        //rk[ith].archive_GRN_stationary(ith);//分裂直前を記録 こっちが先
        //Calculate_Efficiency();//fitnessの計算
        for(ith=0;ith<rk[0].GetTotalcellnumber();ith++)
        {
            //rk[ith].calculate_Fitness();
            //rk[ith].Calculate_Fitness_Target();
            //rk[ith].Calculate_Fitness_Target_parasite();
            rk[ith].archive_evo_mini_stationary(time,ith);//分裂直前を記録
            rk[ith].SetStationaryTime(rk[ith].GetTime());//分裂時間を更新
            for(int j=0;j<ProteinVAL;j++)
            {
                rk[ith].SetStationaryResult(rk[ith].GetResult_average(j),j);//分裂時の発現量を更新
                for(int l=0;l<ProteinVAL;l++)
                {
                    rk[ith].SetStationaryJ(rk[ith].GetJ(j,l), j, l);//分裂時のネットワークを更新
                }
            }
        }
        //histo(time);
        archive_sim_sta(time);//類似度をここで測る
        archive_var_sta(time);//分散をここで測る
        //Calculate_Cumulative_p();//ルーレット用の累積確率を求める
        Calculate_Cumulative_p_power();//ルーレット用の累積確率を求める power版
        Selection_WF_norule(time);//randomに取り除く
        //Selection_Pressure_norule(time);
    }
}

void CELL::check_tree_history(int mother,int daughter)
{
    int type,ch_type,total_number,mother_size,i;
    type = rk[mother].GetCelltype();
    total_number = rk[mother].GetTotalcellnumber();
    mother_size = (int)rk[mother].tree_history.size();
    string tree;
    
    //例えば1-001-001-000-000 ---> tree=1-001-001とする
    // 1-000-000-000-000 --->tree=1 1-001-000-000-000 --->tree=1-001
    
    if(type<10)//cell_typeが一桁の時　*-100-000 size=1+3i
    {
        tree.push_back(rk[mother].tree_history[0]);//先頭だけ
        for(i=1;i<(mother_size+2)/3;i++)
        {
            if(rk[mother].tree_history[3*i-2] != '0')//*-100-000-
            {
                tree.push_back(rk[mother].tree_history[3*i-2]);
                tree.push_back(rk[mother].tree_history[3*i-1]);
                tree.push_back(rk[mother].tree_history[3*i]);
            }
            else if(rk[mother].tree_history[3*i-1] != '0')//-010-000-
            {
                tree.push_back(rk[mother].tree_history[3*i-2]);
                tree.push_back(rk[mother].tree_history[3*i-1]);
                tree.push_back(rk[mother].tree_history[3*i]);
            }
            else//-001-000-
            {
                if(rk[mother].tree_history[3*i] == '0') break;//-000-000-
                else//-001-000-
                {
                    tree.push_back(rk[mother].tree_history[3*i-2]);
                    tree.push_back(rk[mother].tree_history[3*i-1]);
                    tree.push_back(rk[mother].tree_history[3*i]);
                }
            }
        }
    }
    else//**-100-000 size=2+3i
    {
        tree.push_back(rk[mother].tree_history[0]);//先頭だけ
        tree.push_back(rk[mother].tree_history[1]);//先頭だけ
        for(i=1;i<(mother_size+1)/3;i++)
        {
            if(rk[mother].tree_history[3*i-1] != '0')//**-100-000-
            {
                tree.push_back(rk[mother].tree_history[3*i-1]);
                tree.push_back(rk[mother].tree_history[3*i]);
                tree.push_back(rk[mother].tree_history[3*i+1]);
            }
            else if(rk[mother].tree_history[3*i] != '0')//**-010-000-
            {
                tree.push_back(rk[mother].tree_history[3*i-1]);
                tree.push_back(rk[mother].tree_history[3*i]);
                tree.push_back(rk[mother].tree_history[3*i+1]);
            }
            else//-001-000-
            {
                if(rk[mother].tree_history[3*i+1] == '0') break;//**-000-000-
                else//**-001-000-
                {
                    tree.push_back(rk[mother].tree_history[3*i-1]);
                    tree.push_back(rk[mother].tree_history[3*i]);
                    tree.push_back(rk[mother].tree_history[3*i+1]);
                }
            }
        }
        
    }
    
    //cout<<"type="<<type<<endl;
    //cout<<"tree="<<tree<<endl;
    //cout<<"tree_history["<<mother<<"]="<<rk[mother].tree_history<<endl;
    
    rk[daughter].tree_history.clear();
    for(i=0;i<3*4;i++)
    {
        rk[daughter].tree_history.push_back('0');//-000-000-000-000作成　下で頭の文字を与える
    }
    rk[daughter].tree_history = to_string(type) + rk[daughter].tree_history;//親のcelltypeを子供に渡す
    
    mother_size = (int)tree.size();//不要な0要素を消去した親のsize　例えば1-002-000- --->1-002
    
    if(mother_size==1)//'*-000-000-000-000'---->'*-001-000-000-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(1,3)<<endl;
        
        if(rk[(int)maxIndex].tree_history[3]=='9')//*-009-000-000の時
        {
            if(rk[(int)maxIndex].tree_history[2]=='9')//*-099-000-000の時
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1] +1;//新しい系列の番号を与える
                rk[daughter].tree_history[2] = '0';
                rk[daughter].tree_history[3] = '0';
            }
            else//*-189-000-000の時
                
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1] ;
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2] +1;//新しい系列の番号を与える
                rk[daughter].tree_history[3] = '0';
            }
            
        }
        else//*-008-000-000の時
        {
            rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1] ;
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2] ;
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3] +1;
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    else if(mother_size==2)//'**-000-000-000-000'---->'**-001-000-000-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(2,3)<<endl;
        
        if(rk[(int)maxIndex].tree_history[4]=='9')//**-009-000-000の時
        {
            if(rk[(int)maxIndex].tree_history[3]=='9')//**-099-000-000の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2] +1;//新しい系列の番号を与える
                rk[daughter].tree_history[3] = '0';
                rk[daughter].tree_history[4] = '0';
            }
            else//**-189-000-000の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2] ;
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3] +1;//新しい系列の番号を与える
                rk[daughter].tree_history[4] = '0';
            }
            
        }
        else//**-008-000-000の時
        {
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2] ;
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3] ;
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4] +1;
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else if(mother_size==4)//'*-***-000-000-000'---->'*-***-001-000-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,4) != rk[mother].tree_history.substr(0,4))
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(2,3)<<endl;
        
        if(rk[(int)maxIndex].tree_history[6]=='9')//*-***-009-000-000の時
        {
            if(rk[(int)maxIndex].tree_history[5]=='9')//*-***-099-000-000の時
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[5] = '0';
                rk[daughter].tree_history[6] = '0';
            }
            else//*-***-019-000-000の時
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[6] = '0';
            }
        }
        else//*-08-08-00の時
        {
            rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6]+1;
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else if(mother_size==5)//'**-***-000-000-000'---->'**-***-001-000-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,5) != rk[mother].tree_history.substr(0,5))
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(2,3)<<endl;
        
        if(rk[(int)maxIndex].tree_history[7]=='9')//*-***-009-000-000の時
        {
            if(rk[(int)maxIndex].tree_history[6]=='9')//*-***-099-000-000の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[6] = '0';
                rk[daughter].tree_history[7] = '0';
            }
            else//*-***-019-000-000の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[7] = '0';
            }
        }
        else//*-08-08-00の時
        {
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
            rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7]+1;
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else if(mother_size==7)//'*-***-***-000-000'---->'*-***-***-001-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,7) != rk[mother].tree_history.substr(0,7))
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(1,6)<<endl;
        
        if(rk[(int)maxIndex].tree_history[9]=='9')//*-***-***-009-000の時
        {
            if(rk[(int)maxIndex].tree_history[8]=='9')//*-***-***-099-000の時
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[8] = '0';
                rk[daughter].tree_history[9] = '0';
            }
            else
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[9] = '0';
            }
        }
        else//*-**-**-08-00の時
        {
            rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
            rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
            rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
            rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9]+1;//新しい系列の番号を与える
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else if(mother_size==8)//'**-***-***-000-000'---->'**-***-***-001-000'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,8) != rk[mother].tree_history.substr(0,8))
            {
                sortv[i] = 0;
            }
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(2,6)<<endl;
        
        if(rk[(int)maxIndex].tree_history[10]=='9')//*-***-***-009-000の時
        {
            if(rk[(int)maxIndex].tree_history[9]=='9')//*-***-***-099-000の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[9] = '0';
                rk[daughter].tree_history[10] = '0';
            }
            else
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
                rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[10] = '0';
            }
        }
        else//*-**-**-08-00の時
        {
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
            rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
            rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
            rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
            rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10]+1;//新しい系列の番号を与える
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else if(mother_size==10) //'*-***-***-***-000'---->'*-***-***-***-008'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,10) != rk[mother].tree_history.substr(0,10))
            {
                sortv[i] = 0;
            }
            
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(1,9)<<endl;
        
        if(rk[(int)maxIndex].tree_history[12]=='9')//*-***-***-***-009の時
        {
            if(rk[(int)maxIndex].tree_history[11]=='9')//*-***-***-***-099の時
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
                rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
                rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[11] = '0';
                rk[daughter].tree_history[12] = '0';
            }
            else
            {
                rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
                rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
                rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10];
                rk[daughter].tree_history[11] = rk[(int)maxIndex].tree_history[11]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[12] = '0';
            }
        }
        else//*-**-**-08の時
        {
            rk[daughter].tree_history[1] = rk[(int)maxIndex].tree_history[1];
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
            rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
            rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
            rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
            rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10];
            rk[daughter].tree_history[11] = rk[(int)maxIndex].tree_history[11];
            rk[daughter].tree_history[12] = rk[(int)maxIndex].tree_history[12]+1;//新しい系列の番号を与える
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
    else //'*-***-***-***-000'---->'*-***-***-***-008'
    {
        vector<int> sortv;
        sortv.clear();
        sortv.resize(total_number);
        for(i=0;i<total_number;i++)
        {
            //cout<<"rk["<<i<<"].tree_history="<<rk[i].tree_history<<endl;
            //sortv[i]=sstoi(rk[i].tree_history);
            sortv[i]=sstod(rk[i].tree_history)*0.000001;
            //cout<<"sortv["<<i<<"]="<< sortv[i]<<endl;
        }
        for(int i=0;i<total_number;i++)
        {
            ch_type = rk[i].GetCelltype();
            //cout<<"ch_type="<<ch_type<<endl;
            if(ch_type != type)
            {
                sortv[i] = 0;
            }
            else if(rk[i].tree_history.substr(0,11) != rk[mother].tree_history.substr(0,11))
            {
                sortv[i] = 0;
            }
            
            //cout<<"sortv["<<i<<"]="<<sortv[i]<<endl;
        }
        
        vector<int>::iterator maxIt = max_element(sortv.begin(), sortv.end());
        
        //int max = *max_element(sortv.begin(), sortv.end());
        
        //cout<<"max="<<max<<endl;
        
        size_t maxIndex = distance(sortv.begin(), maxIt);
        
        //cout<<"maxIndex="<<maxIndex<<endl;
        //cout<<"max tree_history="<<rk[(int)maxIndex].tree_history<<endl;
        //cout<<"系統番号="<<rk[(int)maxIndex].tree_history.substr(2,9)<<endl;
        
        if(rk[(int)maxIndex].tree_history[13]=='9')//*-***-***-***-009の時
        {
            if(rk[(int)maxIndex].tree_history[12]=='9')//*-***-***-***-099の時
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
                rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
                rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10];
                rk[daughter].tree_history[11] = rk[(int)maxIndex].tree_history[11]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[12] = '0';
                rk[daughter].tree_history[13] = '0';
            }
            else
            {
                rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
                rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
                rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
                rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
                rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
                rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
                rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
                rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
                rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10];
                rk[daughter].tree_history[11] = rk[(int)maxIndex].tree_history[11];
                rk[daughter].tree_history[12] = rk[(int)maxIndex].tree_history[12]+1;//新しい系列の番号を与える
                rk[daughter].tree_history[13] = '0';
            }
        }
        else//*-**-**-08の時
        {
            rk[daughter].tree_history[2] = rk[(int)maxIndex].tree_history[2];
            rk[daughter].tree_history[3] = rk[(int)maxIndex].tree_history[3];
            rk[daughter].tree_history[4] = rk[(int)maxIndex].tree_history[4];
            rk[daughter].tree_history[5] = rk[(int)maxIndex].tree_history[5];
            rk[daughter].tree_history[6] = rk[(int)maxIndex].tree_history[6];
            rk[daughter].tree_history[7] = rk[(int)maxIndex].tree_history[7];
            rk[daughter].tree_history[8] = rk[(int)maxIndex].tree_history[8];
            rk[daughter].tree_history[9] = rk[(int)maxIndex].tree_history[9];
            rk[daughter].tree_history[10] = rk[(int)maxIndex].tree_history[10];
            rk[daughter].tree_history[11] = rk[(int)maxIndex].tree_history[11];
            rk[daughter].tree_history[12] = rk[(int)maxIndex].tree_history[12];
            rk[daughter].tree_history[13] = rk[(int)maxIndex].tree_history[13]+1;//新しい系列の番号を与える
        }
        //cout<<"rk["<<daughter<<"].tree_history="<<rk[daughter].tree_history<<endl;
    }
    
}

void CELL::archive_population(int time)
{
    FILE *pf = nullptr;
    FILE *pg = nullptr;
    
    int l,_number=0,totalcellnum=0,celltypenumer=Species_numb;
    double total=0.0;
    int *n;double *ro,*sumofg,*ag;
    n = new int[celltypenumer];
    ro = new double[celltypenumer];
    sumofg = new double[celltypenumer];
    ag = new double[celltypenumer];
    
    totalcellnum = rk[0].GetTotalcellnumber();
    for(l=0;l<celltypenumer;l++)
    {
        n[l]= 0;
        ro[l] = 0.0;
        sumofg[l] = 0.0;
        ag[l] = 0.0;
    }
    for(l=0;l<totalcellnum;l++)
    {
        _number = rk[l].GetCelltype();
        for(int i=0;i<celltypenumer;i++)
        {
            if(_number%celltypenumer == 0)
            {
                ++n[celltypenumer-1];
                sumofg[celltypenumer-1] = sumofg[celltypenumer-1] + rk[l].GetGrowthRate();
                break;
            }
            else if(_number%celltypenumer == i)
            {
                ++n[i-1];
                sumofg[i-1] = sumofg[i-1] + rk[l].GetGrowthRate();
                break;
            }
        }
        _number = 0;
    }
    
    total = rk[0].GetTotalcellnumber();
    for(l=0;l<celltypenumer;l++)
    {
        //cout<<n[l]<<endl;
        ro[l] = (double)(n[l]/total);
        if(n[l]!=0)
        {
            ag[l] = (double)(sumofg[l]/n[l]);
        }
        else{
            ag[l] = 0.0;
        }
        //cout<<ro[l]<<endl;
    }
    
    pf = fopen("population.txt","a+");
    pg = fopen("averagegrowth.txt","a+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",time);
        for(l=0;l<celltypenumer;l++)
        {
            fprintf(pf," %f ",ro[l]);
        }
        fprintf(pf,"\n");
        fprintf(pg," %d ",time);
        for(l=0;l<celltypenumer;l++)
        {
            fprintf(pg," %f ",ag[l]);
        }
        fprintf(pg,"\n");
    }
    
    fclose(pf);
    fclose(pg);
    delete []n;
    delete []ro;
    delete []sumofg;
    delete []ag;
}

void CELL::archive_target_input(int gene)
{
    FILE *pf = nullptr;
    FILE *pr = nullptr;
    //FILE *pg = nullptr;
    FILE *pp = nullptr;
    
    int l,_number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("target.txt","a+");
    pr = fopen("simi.txt","w+");
    //pg = fopen("Labels.txt","w+");
    pp = fopen("cell_number.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",gene);
            fprintf(pf," %d " ,ith+1);
            fprintf(pf," %d " ,rk[ith].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[ith].tree_history.c_str());//変異数を確認
            for(l=0;l<target_number;l++)
            {
                fprintf(pf," %.8f ",rk[ith].GetResult(l));//targetの発現量を記録
            }
            for(l=ProteinVAL-2*input_number;l<ProteinVAL-input_number;l++)
            {
                fprintf(pf," %.8f",rk[ith].GetResult(l));//input gene を記録
            }
            for(l=ProteinVAL-input_number;l<ProteinVAL;l++)
            {
                fprintf(pf," %.8f",rk[ith].GetResult(l));//cell中のtransporterを記録
            }
            fprintf(pf," %.8f ",rk[ith].GetMu());//発現コストを記録
            
            for(l=0;l<input_number;l++)
            {
                fprintf(pf," %.8f",rk[ith].GetYin(l));//cell中のinput物質を記録
            }
            for(l=0;l<input_number;l++)
            {
                fprintf(pf," %.8f",rk[ith].GetY(l));//medium中のinput物質を記録
            }
            fprintf(pf,"\n");
            
            //fprintf(pr," %d " ,ith+1);//simirarity計算用
            for(l=0;l<ProteinVAL;l++)
            {
                fprintf(pr," %.8f ",rk[ith].GetResult(l));//targetの発現量を記録
                fprintf(pr,",");
            }
            fprintf(pr,"\n");
            //fprintf(pg,"\"%d\"" ,ith+1);//simirarity用ラベル
            //fprintf(pg,"," );
            //fprintf(pg,"\n");
        }
        fprintf(pp," %d ",gene);
        fprintf(pp," %d ",rk[0].GetTotalcellnumber());
        fprintf(pp,"\n");
    }
    fclose(pf);
    fclose(pr);
    //fclose(pg);
    fclose(pp);
}

void CELL::archive_sort_cell(int time)//sortしたものを出力
{
    int totalcellnum=0,ith=0;
    
    char filename[50]={};
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    for(ith=0;ith<rk[0].GetTotalcellnumber();ith++)
    {
        FILE *ph = nullptr;
        sprintf(filename,"Sort_cell%04d.txt",ith+1);
        ph = fopen(filename,"a+");
        fprintf(ph," %d ",time);
        fprintf(ph," #%d " ,rk[ith].GetSorting_num());
        fprintf(ph," %d " ,rk[rk[ith].GetSorting_num()].GetParasite_attacked_number());
        for(int l=0;l<ProteinVAL;l++)
        {
            fprintf(ph," %f ",rk[rk[ith].GetSorting_num()].GetResult(l));
        }
        fprintf(ph,"\n");
        fclose(ph);
    }
    
}

void CELL::archive_growth_volume(int time)
{
    int totalcellnum=0,ith=0;
    
    char filename[50]={};
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    FILE *ph = nullptr;
    sprintf(filename,"volume.txt");
    ph = fopen(filename,"a+");
    
    for(ith=0;ith<totalcellnum;ith++)
    {
        fprintf(ph," %d ",time);
        fprintf(ph," #%d #%d " ,ith,rk[ith].GetParasite_attacked_number());
        fprintf(ph," %f  %f " ,rk[ith].GetGrowthRate(),rk[ith].GetVolume());
        fprintf(ph,"\n");
    }
    fclose(ph);
}

RK::RK(double t0,int div,double tn0,int _protein_number,int _terget_protein_number,double mutTime,double diviTime,int NumberofCell,int Numberofonecelltype,double _beta,double _eps,double volume,double threshold):LOOP(div),totalcellnumber(NumberofCell),protein_x_number(_protein_number),time(t0),startTime(t0),finishTime(tn0),mutationTime(mutTime),divisionTime(diviTime),beta(_beta),eps(_eps),v(volume),initial_volume(volume)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);        // [0, 1] 範囲の一様乱数　rand01(mt)と使う
    uniform_real_distribution<> rand0th(0, threshold);
    
    //vectorを以下で定義
    result.clear();//proteinの発現量
    result.resize(_protein_number);
    fill(result.begin(),result.end(),0.0);
    for(int i =0;i<_protein_number;i++)
    {
        //result.at(i) = rand01(mt);//初期値設定
        result.at(i) = rand0th(mt);//閾値以下の初期値
        //result.at(i) = gauss_rand(0.01,0.001);//0付近からスタート
        //result.at(i) = 0.0;
    }
    result_average.clear();//400~500秒の平均
    result_average.resize(ProteinVAL);
    fill(result_average.begin(),result_average.end(),0.0);
    result_check.clear();
    result_check.resize(ProteinVAL);
    for(int i=0;i<target_number;i++)
    {
        for(int l=0;l<check_num;l++)
        {
            result_check[i].push_back(0.0);
        }
    }
    
    result_VipVg.clear();//400~500秒の平均
    result_VipVg.resize(ProteinVAL);
    fill(result_VipVg.begin(),result_VipVg.end(),0.0);
    for(int i =0;i<_protein_number;i++)
    {
        result_VipVg.at(i) = rand0th(mt);//閾値以下の初期値
    }
    
    StationaryResult.clear();
    StationaryResult.resize(_protein_number);
    fill(StationaryResult.begin(),StationaryResult.end(),0.0);
    
    J.clear();
    J.resize(_protein_number);
    for(int i=0;i<_protein_number;i++)
    {
        for(int l=0;l<_protein_number;l++)
        {
            J[i].push_back(0);
        }
    }
    
    J_VipVg.clear();
    J_VipVg.resize(_protein_number);
    for(int i=0;i<_protein_number;i++)
    {
        for(int l=0;l<_protein_number;l++)
        {
            J_VipVg[i].push_back(0);
        }
    }
    
    StationaryJ.clear();
    StationaryJ.resize(_protein_number);
    for(int i=0;i<_protein_number;i++)
    {
        for(int l=0;l<_protein_number;l++)
        {
            StationaryJ[i].push_back(0);
        }
    }
    
    EL.clear();
    EL.resize(_protein_number);
    
    WEL.clear();
    WEL.resize(_protein_number);
    
    SEL.clear();
    for(int l=0;l<_protein_number;l++)
    {
        SEL.push_back(0);//手の本数を数える
    }
    
    REL.clear();
    REL.resize(_protein_number);
    
    RSEL.clear();
    for(int l=0;l<_protein_number;l++)
    {
        RSEL.push_back(0);//手の本数を数える
    }
    
    EL_VipVg.clear();
    EL_VipVg.resize(_protein_number);
    
    WEL_VipVg.clear();
    WEL_VipVg.resize(_protein_number);
    
    SEL_VipVg.clear();
    for(int l=0;l<_protein_number;l++)
    {
        SEL_VipVg.push_back(0);//手の本数を数える
    }
    
    REL_VipVg.clear();
    REL_VipVg.resize(_protein_number);
    
    RSEL_VipVg.clear();
    for(int l=0;l<_protein_number;l++)
    {
        RSEL_VipVg.push_back(0);//手の本数を数える
    }
    
    input_list.clear();
    input_list.resize(depth);
    for(int i=0;i<depth;i++)
    {
        for(int l=0;l<_protein_number;l++)
        {
            input_list[i].push_back(0);
        }
    }
    
    Meanfieldval.clear();
    Meanfieldval.resize(_protein_number);
    fill(Meanfieldval.begin(),Meanfieldval.end(),0.0);
    proteinD.clear();
    proteinD.resize(_protein_number);
    fill(proteinD.begin(),proteinD.end(),0.0);
    threshold_g.clear();
    threshold_g.resize(_protein_number);
    fill(threshold_g.begin(),threshold_g.end(),threshold);
    k_delta.clear();
    k_delta.resize(_protein_number);
    fill(k_delta.begin(),k_delta.end(),0);
    
    TargetPattern.clear();
    TargetPattern.resize(target_number);
    fill(TargetPattern.begin(),TargetPattern.end(),0.0);
    
    infection_virulence.clear();
    infection_virulence.resize(Parasite_species_num);
    fill(infection_virulence.begin(),infection_virulence.end(),0.0);
    //それ以外
    delt = (double)((tn0-t0)/div);
    time1000 = (int)(time*1000);//時間1秒表示用
    matationTime1000 = (int)(mutationTime*1000);
    divisionTime1000 = (int)(divisionTime*1000);
    objnum++;
    number = objnum;
    average_growth =0.0;
    growth_rate = 0.0;
    division_number=0;
    mutaion_number=0;
    StationaryTime =0.0;
    q = 0.0;
    fitness=0.0;
    fitness_p = 0.0;
    cumulative_p=0.0;
    gene =0;
    sorting_number =0;
    para_fitness = 0.0;
    sigma_str = Sigma;
    mu=0.0;//発現コスト用
    
    for(int i=2;i<=Species_numb+1;i++)//=つけた
    {
        if(number<=1*Numberofonecelltype)
        {
            celltype = 1;
        }
        else if(number<=i*Numberofonecelltype)
        {
            celltype = i;
            break;
        }
    }
    tree_history = to_string(celltype) +tree_history;
    for(int i=0;i<3*4;i++)
    {
        tree_history.push_back('0');
    }
    
    ori_Sstr_J.reserve(_protein_number*_protein_number);//Jを分ここで確保しておく
    Sstr_J.reserve(_protein_number*_protein_number);//Jをf分ここで確保しておく
    hammingD =0;//初期はゼロ
    RelativehammingD =0;
    Yin_J.shrink_to_fit();
    Yin_J.resize(input_number);
    for(int i=0;i<input_number;i++)
    {
        Yin_J[i].resize(ProteinVAL);
    }
    
    for(int i=0;i<input_number;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            Yin_J[i][l] = 0;
        }
    }
    
    Y.shrink_to_fit();//Medium中のinput物質
    Y.resize(input_number);
    fill(Y.begin(),Y.end(),S_resource);
    
    X.shrink_to_fit();//Medium中の拡散物質
    X.resize(ProteinVAL);
    fill(X.begin(),X.end(),0.0);
    
    Yin.shrink_to_fit();//Medium中のinput物質
    Yin.resize(input_number);
    
    //fill(Yin.begin(),Yin.end(),S_resource);
    for(int l=0;l<input_number;l++)
    {
        Yin.at(l)= S_resource*rand01(mt);
        Yin.at(l)= 0.0;
    }
    
    Medium_D.shrink_to_fit();
    Medium_D.resize(input_number);
    fill(Medium_D.begin(),Medium_D.end(),0.0);
    
    //parasite用
    parasite_attacked_number=0;//hostの種類を決める　2**n種
    parasite_population =0.0;
    Parasite_population_rate.shrink_to_fit();//全パラサイトの情報
    Parasite_population_rate.resize(Parasite_species_num);
    fill(Parasite_population_rate.begin(),Parasite_population_rate.end(),0.0);
}

RK::~RK()
{
    cout<<"Class RK is destroyed"<<endl;
}

RK::RK(const RK& r)
{}

void RK::ShowDate()
{
    int i,l;
    cout<<"---------------------"<<endl;
    for(i=0;i<result.size();i++)
    {
        cout<<"protein["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<result.size();i++)
    {
        //cout<<"MeanVal["<<i<<"]="<<Meanfieldval.at(i)<<endl;
    }
    for(i=0;i<J.size();i++)
    {
        for(l=0;l<J.size();l++)
        {
            //cout<<"J["<<i<<"]["<<l<<"]="<<J[i][l]<<endl;
        }
    }
    for(i=0;i<result.size();i++)
    {
        cout<<"proteinD["<<i<<"]="<<proteinD[i]<<endl;
    }
    for(i=0;i<result.size();i++)
    {
        cout<<"g["<<i<<"]="<<threshold_g[i]<<endl;
    }
    for(i=0;i<J.size();i++)
    {
        //cout<<"k_delta["<<i<<"]="<<k_delta[i]<<endl;
    }
    cout<<"celltype="<<celltype<<endl;
    cout<<"tree_history="<<tree_history<<endl;
    cout<<"growth_rate="<<growth_rate<<endl;
    cout<<"volume="<<v<<endl;
    cout<<"beta="<<beta<<endl;
    //cout<<"esi_sub="<<eps<<endl;
    cout<<"delt="<<delt<<endl;
    cout<<"time="<<time<<endl;
    cout<<"num="<<number<<endl;
    cout<<"totalcellnumber="<<totalcellnumber<<endl;
    cout<<"---------------------"<<endl;
}
void RK::Show()
{
    int i;
    cout<<"---------------------"<<endl;
    cout<<"time="<<time<<endl;
    cout<<"celltype="<<celltype<<endl;
    cout<<"tree_history="<<tree_history<<endl;
    cout<<"growth_rate="<<growth_rate<<endl;
    cout<<"volume="<<v<<endl;
    cout<<"num="<<number<<endl;
    for(i=0;i<result.size();i++)
    {
        cout<<"protein["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<result.size();i++)
    {
        cout<<"MeanVal["<<i<<"]="<<Meanfieldval.at(i)<<endl;
    }
    cout<<"---------------------"<<endl;
}

void RK::archive()
{
    FILE *pf = nullptr;
    FILE *pr = nullptr;
    FILE *pg = nullptr;
    
    int l,_number=0,p_val=0;
    double time=0,delttime=0;
    char filename[50]={};
    char filenamer[50] ={};
    char filenameg[50] ={};
    double *x,*g;
    x = new double[protein_x_number];
    g = new double[protein_x_number];
    
    _number  = GetNum();
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    sprintf(filename,"CELL%03d.txt",_number);
    sprintf(filenamer,"GROWTH%03d.txt",_number);
    sprintf(filenameg,"Threshold_g%03d.txt",_number);
    
    pf = fopen(filename,"a+");
    pr = fopen(filenamer,"a+");
    pg = fopen(filenameg,"a+");
    
    if(pf==NULL||pr==NULL||pg==NULL)
    {
        cout << "Could not open file"<<endl;
        cout <<pf<<endl;
        cout <<pr<<endl;
        cout <<pg<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %f ",time);
        for(l=0;l<p_val;l++)
        {
            x[l] = GetResult(l);
            fprintf(pf," %f ",x[l]);
        }
        
        fprintf(pf,"\n");
        
        fprintf(pr," %f " ,time);
        fprintf(pr," %f ",growth_rate);
        fprintf(pr,"\n");
        
        fprintf(pg," %f ",time);
        for(l=0;l<p_val;l++)
        {
            g[l] = GetThreshold_g(l);
            fprintf(pg," %f ",g[l]);
        }
        fprintf(pg,"\n");
    }
    //time = time+delttime;
    //SetTime(time);
    fclose(pf);
    fclose(pr);
    fclose(pg);
    delete []x;
    delete []g;
}


void RK::archive_evo(int gene,int ith)
{
    
    FILE *pg = nullptr;
    FILE *ppg = nullptr;
    FILE *pd = nullptr;
    
    int l,_number=0,p_val=0;
    double time=0,delttime=0;
    
    char filenameg[50] ={};
    char filenamepg[50] ={};
    char filenamepd[50] ={};
    
    _number  = GetCelltype() ;//同じ細胞を同一ファイルに記録
    //_number = ith+1;
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    sprintf(filenameg,"threshold_g%04d.txt",_number);
    sprintf(filenamepg,"p-g%04d.txt",_number);
    sprintf(filenamepd,"proteinD%04d.txt",_number);
    pg = fopen(filenameg,"a+");
    ppg = fopen(filenamepg,"a+");
    pd = fopen(filenamepd,"a+");
    
    if(pg==NULL||ppg==NULL||pd==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pg," %d " ,gene);
        fprintf(pg," %d " ,ith);
        for(l=0;l<p_val;l++)
        {
            fprintf(pg," %f ",GetThreshold_g(l));
        }
        fprintf(pg,"\n");
        
        fprintf(ppg," %d " ,gene);
        fprintf(ppg," %d " ,ith);
        for(l=0;l<p_val;l++)
        {
            fprintf(ppg," %f ",GetThreshold_g(l));
        }
        for(l=0;l<p_val;l++)
        {
            fprintf(ppg," %f ",GetResult(l));
        }
        fprintf(ppg,"\n");
        
        fprintf(pd," %d " ,gene);
        fprintf(pd," %d " ,ith);
        for(l=0;l<p_val;l++)
        {
            fprintf(pd," %f ",GetProD(l));
        }
        fprintf(pd,"\n");
    }
    
    fclose(pg);
    fclose(ppg);
    fclose(pd);
    
}

void RK::archive_evo_mini(int gene,int ith)
{
    FILE *pf = nullptr;
    //FILE *pr = nullptr;
    
    
    int l,_number=0,p_val=0,k,memo=0;
    double time=0,delttime=0;
    char filename[50]={};
    //char filenamer[50] ={};
    
    int iDist1 = 0;
    
    _number  = GetCelltype() ;//同じ細胞を同一ファイルに記録
    
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    sprintf(filename,"CELL%d.txt",_number);
    //sprintf(filenamer,"cell%04d.txt",ith+1);
    
    pf = fopen(filename,"a+");
    //pr = fopen(filenamer,"a+");
    
    Sstr_J.clear();
    Sstr_J.reserve(p_val*p_val);
    for(l=0;l<p_val;l++)
    {
        for(k=0;k<p_val;k++)
        {
            memo = GetJ(l,k);
            if(memo == -1)
            {
                Sstr_J = Sstr_J + 'a' ;//初期のネットワーク構造 GRN -1を文字aで記録
            }
            else if(memo == 0)
            {
                Sstr_J = Sstr_J + 'b' ;//初期のネットワーク構造 GRN 0を文字bで記録
            }
            else
            {
                Sstr_J = Sstr_J + 'c' ;//初期のネットワーク構造 GRN 1を文字cで記録
            }
        }
    }
    
    iDist1 = iHammingDist(Sstr_J,ori_Sstr_J);
    SetHammingD(iDist1);
    
    if(pf==NULL)
        //if(pf==NULL||pr==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",gene);
        fprintf(pf," %d " ,ith+1);
        //fprintf(pf," %s ",tree_history.c_str());
        for(l=0;l<p_val;l++)
        {
            fprintf(pf," %f ",GetResult(l));
        }
        for(l=0;l<input_number;l++)
        {
            fprintf(pf," %.8f",GetYin(l));//cell中のinput物質を記録
        }
        for(l=0;l<input_number;l++)
        {
            fprintf(pf," %.8f",GetY(l));//medium中のinput物質を記録
        }
        fprintf(pf,"\n");
        
        //fprintf(pr," %d ",gene);
        //fprintf(pr," %d " ,ith+1);
        //fprintf(pr," %s ",tree_history.c_str());
        for(l=0;l<p_val;l++)
        {
            //fprintf(pr," %f ",GetResult(l));
        }
        
        //fprintf(pr," %f ",GetMu());//発現コストmu
        for(l=0;l<input_number;l++)
        {
            //fprintf(pr," %.8f",GetYin(l));//cell中のinput物質を記録
        }
        for(l=0;l<input_number;l++)
        {
            //fprintf(pr," %.8f",GetY(l));//medium中のinput物質を記録
        }
        //fprintf(pr,"\n");
    }
    
    fclose(pf);
    //fclose(pr);
}

void RK::archive_Meanfieldval(int gene)
{
    FILE *pf = nullptr;
    
    int l,_number=0,p_val=0;
    double time=0,delttime=0;
    
    _number  = GetCelltype() ;//同じ細胞を同一ファイルに記録
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    pf = fopen("Meanfield.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",gene);
        for(l=0;l<p_val;l++)
        {
            fprintf(pf," %f ",GetMeanfieldval(l));
        }
        fprintf(pf,"\n");
        
    }
    
    fclose(pf);
}

void RK::archive_GRN(int ith)
{
    FILE *pf = nullptr;
    FILE *pr = nullptr;
    
    int l,k,_number=0,p_val=0;
    double time=0,delttime=0;
    char filename[50]={},filenamer[50]={};
    
    _number  = ith+1;
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    sprintf(filename,"GRN%03d.txt",_number);
    sprintf(filenamer,"GRN_Edge%03d.txt",_number);
    pf = fopen(filename,"w+");
    pr = fopen(filenamer,"w+");
    
    if(pf==NULL||pr==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        //fprintf(pf," %#f ",time);
        //fprintf(pf,"\n");
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                fprintf(pf," %f ",GetJ(k,l));
                //fprintf(pf,",");
            }
            fprintf(pf,"\n");
        }
        fprintf(pr,"#num #RSEL #SEL");
        fprintf(pr,"\n");
        for(k=0;k<p_val;k++)
        {
            fprintf(pr," %d ",k);
            fprintf(pr," %d ",GetRSEL(k));
            fprintf(pr," %d ",GetSEL(k));
            fprintf(pr,"\n");
        }
    }
    fclose(pf);
    fclose(pr);
}

void RK::archive_GRN_every_gene(int ith,int time)
{
    FILE *pf = nullptr;

    
    int l,k,_number=0,p_val=0;
    double delttime=0;
    char filename[50]={},filenamer[50]={},filenameg[50]={};
    
    _number  = ith+1;
    p_val = Getprotein_x_number();
    
    delttime = GetDelt();
    
    sprintf(filename,"GRN_gene%d_%03d.txt",time,_number);
    //sprintf(filenamer,"GRN_Edge_gene%d_%03d.txt",time,_number);
    //sprintf(filenameg,"Reverse_GRN_gene%d_%03d.txt",time,_number);
    pf = fopen(filename,"a+");
    //pr = fopen(filenamer,"a+");
    //pg = fopen(filenameg,"a+");
    
    //if(pf==NULL||pr==NULL||pg==NULL)
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        //fprintf(pf," %d ",time/500);
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                fprintf(pf," %f ",GetJ(k,l));//inputしてくるgeneを記録
            }
            fprintf(pf,"\n");
        }
        //fprintf(pr," %d ",time/500);
        //fprintf(pr,"#num #RSEL #SEL");
        //fprintf(pr,"\n");
        for(k=0;k<p_val;k++)
        {
            //fprintf(pr," %d ",k);
            //fprintf(pr," %d ",GetRSEL(k));
            //fprintf(pr," %d ",GetSEL(k));
            //fprintf(pr,"\n");
        }
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                //fprintf(pg," %d ",GetJ(l,k));//outputする先のgeneを記録
            }
            //fprintf(pg,"\n");
        }
    }
    fclose(pf);
    //fclose(pr);
    //fclose(pg);
}
void RK::archive_GRN_every_gene_one(int ith,int time)
{
    //FILE *pf = nullptr;
    FILE *pr = nullptr;
    FILE *pg = nullptr;
    
    int l,k,_number=0,p_val=0;
    double delttime=0;
    char filename[50]={},filenamer[50]={},filenameg[50]={};
    
    _number  = ith+1;
    p_val = Getprotein_x_number();
    
    delttime = GetDelt();
    
    //sprintf(filename,"GRN_gene%d_%03d.txt",time,_number);
    sprintf(filenamer,"GRN_Edge_gene%d_%03d.txt",time,_number);
    sprintf(filenameg,"Reverse_GRN_gene%d_%03d.txt",time,_number);
    //pf = fopen(filename,"a+");
    pr = fopen(filenamer,"a+");
    pg = fopen(filenameg,"a+");
    
    //if(pf==NULL||pr==NULL||pg==NULL)
    if(pr==NULL||pg==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        //fprintf(pf," %d ",time/500);
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                //fprintf(pf," %d ",GetJ(k,l));//inputしてくるgeneを記録
                //fprintf(pf,",");
            }
            //fprintf(pf,"\n");
        }
        //fprintf(pr," %d ",time/500);
        fprintf(pr,"#num #RSEL #SEL");
        fprintf(pr,"\n");
        for(k=0;k<p_val;k++)
        {
            fprintf(pr," %d ",k);
            fprintf(pr," %d ",GetRSEL(k));
            fprintf(pr," %d ",GetSEL(k));
            fprintf(pr,"\n");
        }
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                fprintf(pg," %f ",GetJ(l,k));//outputする先のgeneを記録
            }
            fprintf(pg,"\n");
        }
    }
    //fclose(pf);
    fclose(pr);
    fclose(pg);
}

void RK::archive_YinJ(int ith)
{
    FILE *pf = nullptr;
    
    int k,_number=0,p_val=0,G=0;
    double time=0;
    char filename[50]={};
    
    _number  = ith+1;
    p_val = ProteinVAL;
    
    time = GetTime();
    
    sprintf(filename,"YinJ%03d.txt",_number);
    pf = fopen(filename,"w+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %f ",time);
        fprintf(pf,"\n");
        for(k=0;k<p_val;k++)
        {
            G=GetYin_J(ith,k);
            fprintf(pf," %d ",G);
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
}

void RK::archive_evo_mini_stationary(double time,int ith)
{
    //FILE *pf = nullptr;
    FILE *pr = nullptr;
    FILE *pg = nullptr;
    
    int l,k,_number=0,p_val=0,memo=0;
    int iDist1 = 0;
    //char filename[50]={};
    char filenamer[50] ={};
    
    //_number  = GetCelltype() ;//同じ細胞を同一ファイルに記録
    _number  = GetParasite_attacked_number() ;//同じ細胞を同一ファイルに記録
    p_val = Getprotein_x_number();
    
    Sstr_J.clear();
    Sstr_J.reserve(p_val*p_val);
    for(l=0;l<p_val;l++)
    {
        for(k=0;k<p_val;k++)
        {
            memo = GetJ(l,k);
            if(memo == -1)
            {
                Sstr_J = Sstr_J + 'a' ;//初期のネットワーク構造 GRN -1を文字aで記録
            }
            else if(memo == 0)
            {
                Sstr_J = Sstr_J + 'b' ;//初期のネットワーク構造 GRN 0を文字bで記録
            }
            else
            {
                Sstr_J = Sstr_J + 'c' ;//初期のネットワーク構造 GRN 1を文字cで記録
            }
        }
    }
    
    iDist1 = iHammingDist(Sstr_J,ori_Sstr_J);
    SetHammingD(iDist1);
    
    //sprintf(filename,"cell_stationary%04d.txt",_number);
    sprintf(filenamer,"fitness%04d.txt",_number+1);
    
    //pf = fopen(filename,"a+");
    pr = fopen(filenamer,"a+");
    pg = fopen("Fitness_stationary.txt","a+");
    
    //if(pf==NULL||pg==NULL)
    if(pg==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        //fprintf(pf," %f ",time);
        //fprintf(pf," %d " ,ith+1);
        //fprintf(pf," %s ",tree_history.c_str());
        //fprintf(pf," %d " ,iDist1);
        for(l=0;l<p_val;l++)
        {
            //fprintf(pf," %.8f ",GetResult_average(l));
        }
        //fprintf(pf," %.8f ",GetFitness_Target_parasite_average());
        //fprintf(pf," %f ",GetMu());//発現コストmu
        for(l=0;l<input_number;l++)
        {
            //fprintf(pf," %.8f",GetYin(l));//cell中のinput物質を記録
        }
        for(l=0;l<input_number;l++)
        {
            //fprintf(pf," %.8f",GetY(l));//medium中のinput物質を記録
        }
        
        //fprintf(pf,"\n");
        
        fprintf(pr," %f " ,time);
        fprintf(pr," %d " ,ith+1);
        fprintf(pr," %f ",GetStationaryTime());
        //fprintf(pr," %s ",tree_history.c_str());
        //fprintf(pr," %d " ,iDist1);
        //fprintf(pr," %.8f ",GetFitness_Target_parasite_average());
        //fprintf(pr," %.8f ",GetFitness_Target_noparasite_average());
        fprintf(pr," %.8f ",GetGrowthRate());
        //fprintf(pr," %.8f ",GetVolume());
        //fprintf(pr," %f ",GetMu());//発現コストmu
        fprintf(pr," %d ",division_number);
        fprintf(pr," %d ",mutaion_number);
        fprintf(pr,"\n");
        
        fprintf(pg," %f " ,time);
        fprintf(pg," %d " ,ith+1);
        fprintf(pg," %d " ,_number);
        fprintf(pg," %f " ,GetStationaryTime());
        //fprintf(pg," %s ",tree_history.c_str());
        //fprintf(pg," %d " ,iDist1);
        //fprintf(pg," %.8f ",GetFitness_Target_parasite_average());
        //fprintf(pg," %.8f ",GetFitness_Target_noparasite_average());
        fprintf(pg," %.8f ",GetGrowthRate());
        //fprintf(pg," %.8f ",GetVolume());
        //fprintf(pg," %f ",GetMu());//発現コストmu
        fprintf(pg," %d ",division_number);
        fprintf(pg," %d ",mutaion_number);
        fprintf(pg,"\n");
        
    }
    
    //fclose(pf);
    fclose(pr);
    fclose(pg);
}

void RK::archive_Meanfieldval_stationary(int gene)
{
    FILE *pf = nullptr;
    
    int l,_number=0,p_val=0;
    double time=0,delttime=0;
    
    _number  = GetCelltype() ;//同じ細胞を同一ファイルに記録
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    pf = fopen("Meanfield_stationary.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",gene);
        for(l=0;l<p_val;l++)
        {
            fprintf(pf," %f ",GetMeanfieldval(l));
        }
        fprintf(pf,"\n");
        
    }
    
    fclose(pf);
}

void RK::archive_GRN_stationary(int ith)
{
    FILE *pf = nullptr;
    
    int l,k,_number=0,p_val=0,memo=0;
    double time=0,delttime=0;
    char filename[50]={};
    
    _number  = ith+1;
    p_val = Getprotein_x_number();
    time = GetTime();
    delttime = GetDelt();
    
    sprintf(filename,"GRN_stationary%03d.txt",_number);
    pf = fopen(filename,"w+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %f ",time);
        fprintf(pf,"\n");
        for(l=0;l<p_val;l++)
        {
            for(k=0;k<p_val;k++)
            {
                memo=GetJ(k,l);
                if(memo==-1)
                {
                    fprintf(pf,"\"%d\"",memo+3);
                }
                else
                {
                    fprintf(pf,"\"%d\"",memo);
                }
                fprintf(pf,",");
            }
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
}


inline void RK::converterJ()
{
    int j=0;
    //EL.shrink_to_fit();
    EL.clear();
    EL.resize(ProteinVAL);
    //WEL.shrink_to_fit();
    WEL.clear();
    WEL.resize(ProteinVAL);
    REL.clear();
    REL.resize(ProteinVAL);
    
    for(int i=0;i<ProteinVAL;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            j = GetJ(i,l);//lからiへのpath
            if(j != 0)
            {
                EL[i].push_back(l);//iへのedgeを記録
                WEL[i].push_back(j);//重み(-1,0,1)を記録
            }
        }
        for(int l=0,n=GetELSize(i);l<n;l++)
        {
            //cout<<"EL["<<i<<"]["<<l<<"]="<<EL[i][l]<<endl;
        }
    }
    
    for(int l=0;l<ProteinVAL;l++)
    {
        for(int i=0;i<ProteinVAL;i++)
        {
            j = GetJ(i,l);//lからiへのpath
            if(j != 0)
            {
                REL[l].push_back(i);//iへのedgeを記録
            }
        }
        for(int i=0,n=GetRELSize(l);i<n;i++)
        {
            //cout<<"REL["<<l<<"]["<<i<<"]="<<REL[l][i]<<endl;
        }
    }
    
    SEL.clear();
    SEL.resize(ProteinVAL);
    for(int l=0;l<ProteinVAL;l++)
    {
        SEL[l]=GetELSize(l);
        //cout<<"SEL["<<l<<"]="<<SEL[l]<<endl;
    }
    
    RSEL.clear();
    RSEL.resize(ProteinVAL);
    for(int l=0;l<ProteinVAL;l++)
    {
        RSEL[l]=GetRELSize(l);
        //cout<<"RSEL["<<l<<"]="<<RSEL[l]<<endl;
    }
}

inline void RK::converterJ_VipVg()
{
    int j=0;
    //EL.shrink_to_fit();
    EL_VipVg.clear();
    EL_VipVg.resize(ProteinVAL);
    //WEL.shrink_to_fit();
    WEL_VipVg.clear();
    WEL_VipVg.resize(ProteinVAL);
    REL_VipVg.clear();
    REL_VipVg.resize(ProteinVAL);
    
    for(int i=0;i<ProteinVAL;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            j = GetJ_VipVg(i,l);//lからiへのpath
            if(j != 0)
            {
                EL_VipVg[i].push_back(l);//iへのedgeを記録
                WEL_VipVg[i].push_back(j);//重み(-1,0,1)を記録
            }
        }
        for(int l=0,n=GetELSize_VipVg(i);l<n;l++)
        {
            //cout<<"EL["<<i<<"]["<<l<<"]="<<EL[i][l]<<endl;
        }
    }
    
    for(int l=0;l<ProteinVAL;l++)
    {
        for(int i=0;i<ProteinVAL;i++)
        {
            j = GetJ_VipVg(i,l);//lからiへのpath
            if(j != 0)
            {
                REL_VipVg[l].push_back(i);//iへのedgeを記録
            }
        }
        for(int i=0,n=GetRELSize_VipVg(l);i<n;i++)
        {
            //cout<<"REL["<<l<<"]["<<i<<"]="<<REL[l][i]<<endl;
        }
    }
    
    SEL_VipVg.clear();
    SEL_VipVg.resize(ProteinVAL);
    for(int l=0;l<ProteinVAL;l++)
    {
        SEL_VipVg[l]=GetELSize_VipVg(l);
        //cout<<"SEL["<<l<<"]="<<SEL[l]<<endl;
    }
    
    RSEL_VipVg.clear();
    RSEL_VipVg.resize(ProteinVAL);
    for(int l=0;l<ProteinVAL;l++)
    {
        RSEL_VipVg[l]=GetRELSize_VipVg(l);
        //cout<<"RSEL["<<l<<"]="<<RSEL[l]<<endl;
    }
}

void RK::calculate_Input_list()
{
    int i,j,l,k=0,path_size=0;
    input_list.clear();
    input_list.resize(depth);
    
    i=0;//depth=1 input直下の遺伝子
    path_size = GetRELSize(ProteinVAL-input_number)+GetRELSize(ProteinVAL-input_number+1);//２本
    for(j=0;j<path_size;j++)
    {
        input_list[i].push_back(GetREL(ProteinVAL-input_number+j,0));
        //cout<<"Input_list["<<i<<"]["<<j<<"]="<<GetInput_list(i, j)<<endl;
    }
    
    for(i=1;i<depth;i++)
    {
        path_size=0;
        l=0;
        path_size =  GetRELSize(GetInput_list(i-1,l));
        for(j=0;j<path_size;j++)
        {
            if(GetREL(GetInput_list(i-1,l),j)!=GetInput_list(i-1,l))
            {
                input_list[i].push_back(GetREL(GetInput_list(i-1,l),j));
                //cout<<"Input_list["<<i<<"]["<<j<<"]="<<GetInput_list(i, j)<<endl;
                //cout<<"J["<<GetREL(GetInput_list(i-1,l),j)<<"]["<<GetInput_list(i-1,l)<<"]="<<GetJ(GetREL(GetInput_list(i-1,l),j),GetInput_list(i-1,l))<<endl;
            }
        }
        for(l=1;l<GetInput_listSize(i-1);l++)
        {
            k=0;
            
            for(j=path_size;j< path_size + GetRELSize(GetInput_list(i-1,l));j++)
            {
                if(GetREL(GetInput_list(i-1,l),k)!=GetInput_list(i-1,l))
                {
                    input_list[i].push_back(GetREL(GetInput_list(i-1,l),k));
                    //cout<<"Input_list["<<i<<"]["<<j<<"]="<<GetInput_list(i, j)<<endl;
                    //cout<<"J["<<GetREL(GetInput_list(i-1,l),k)<<"]["<<GetInput_list(i-1,l)<<"]="<<GetJ(GetREL(GetInput_list(i-1,l),k),GetInput_list(i-1,l))<<endl;
                }
                k++;
            }
            path_size = path_size + GetRELSize(GetInput_list(i-1,l));//for文範囲の更新
        }
    }
}

void RK::SetResult(double x,int lth)
{
    if(x<-1||x>100)
    {
        cout<<"result is inappropriate"<<endl;
        cout<<"result["<<lth<<"]="<<x<<endl;
        exit(1);
    }
    else if(x<0.0)
    {
        result.at(lth) = -x;//0以下は境界で反射させる
    }
    else
    {
        result.at(lth) = x;
    }
}

void RK::SetResult_average(double rl,int lth)
{
    if(rl<-0.1||rl>1000000)
    {
        cout<<"rl is inappropriate"<<endl;
        cout<<"result_average["<<lth<<"]="<<rl<<endl;
        exit(1);
    }
    else
    {
        result_average.at(lth) = rl;
    }
}

void RK::SetResult_check(double rc,int lth,int check_t)
{
    if(rc<-0.1||rc>100)
    {
        cout<<"rc is inappropriate"<<endl;
        cout<<"result_average["<<lth<<"]="<<rc<<endl;
        exit(1);
    }
    else
    {
        result_check[lth][check_t] = rc;
    }
}

void RK::SetResult_VipVg(double x,int lth)
{
    if(x<-1||x>100)
    {
        cout<<"result is inappropriate"<<endl;
        cout<<"result["<<lth<<"]="<<x<<endl;
        exit(1);
    }
    else if(x<0.0)
    {
        result_VipVg.at(lth) = -x;//0以下は境界で反射させる
    }
    else
    {
        result_VipVg.at(lth) = x;
    }
}

void RK::SetJ(double j,int ith,int lth)//j is onoff from lth to ith
{
    J[ith][lth] = j;
}

void RK::SetJ_VipVg(int j,int ith,int lth)//j is onoff from lth to ith
{
    if(j==-1||j==0||j==1)
    {
        J_VipVg[ith][lth] = j;
    }
    else{
        cout<<"j is inappropriate"<<endl;
        exit(1);
    }
}

void RK::SetStationaryJ(int j,int ith,int lth)
{
    if(j==-1||j==0||j==1)
    {
        StationaryJ[ith][lth] = j;
    }
    else{
        cout<<"j is inappropriate"<<endl;
        exit(-1);
    }
}

void RK::SetEL(int j,int ith,int lth)
{
    EL[ith][lth] = j;
}

void RK::SetWEL(int j,int ith,int lth)
{
    if(j==-1||j==0||j==1)
    {
        WEL[ith][lth] = j;
    }
    else{
        cout<<"j is inappropriate"<<endl;
        cout<<j<<endl;
        exit(1);
    }
}

void RK::SetSEL(int ith,int size)
{
    if(size<=ProteinVAL)
    {
        SEL[ith] = size;
    }
    else
    {
        cout<<"size is inappropriate"<<endl;
        cout<<size<<endl;
        exit(1);
    }
}

void RK::SetEL_VipVg(int j,int ith,int lth)
{
    EL[ith][lth] = j;
}

void RK::SetWEL_VipVg(int j,int ith,int lth)
{
    if(j==-1||j==0||j==1)
    {
        WEL[ith][lth] = j;
    }
    else{
        cout<<"j is inappropriate"<<endl;
        cout<<j<<endl;
        exit(1);
    }
}

void RK::SetSEL_VipVg(int ith,int size)
{
    if(size<=ProteinVAL)
    {
        SEL[ith] = size;
    }
    else
    {
        cout<<"size is inappropriate"<<endl;
        cout<<size<<endl;
        exit(1);
    }
}
void RK::SetProteinD(double diff,int lth)
{
    if(diff<-1||diff>10)
    {
        cout<<"proteinD is inappropriate"<<endl;
        cout<<diff<<endl;
        exit(1);
    }
    else
    {
        proteinD[lth] = diff;
    }
}

void RK::SetThreshold_g(double g,int lth)
{
    if(g<-100||g>100)
    {
        cout<<"threshold is inappropriate"<<endl;
        cout<<g<<endl;
        exit(1);
    }
    else
    {
        threshold_g[lth] = g;
    }
}

void RK::SetTime(double nowtime)
{
    if(nowtime<startTime||nowtime/100>finishTime+delt)
    {
        cout<<"nowtime is inappropriate"<<endl;
        cout<<"nowtime="<<nowtime<<endl;
        ShowDate();
        exit(1);
    }
    else
    {
        time = nowtime;
    }
}

void RK::SetTime1000()
{
    double nowtime=0;
    nowtime = GetTime();
    time1000 = (int)(nowtime/GetDelt());
}

void RK::SetGene(int i)
{
    gene = i;
}

void RK::SetSorting_number(int num)
{
    sorting_number=num;
}

void RK::SetMeanfieldval(double sum_x,int lth)
{
    Meanfieldval[lth] = sum_x;
}

void RK::SetTotalcellnumber(int new_totalnumber)
{
    totalcellnumber = new_totalnumber;
}

void RK::SetVolume(double volume)
{
    if(volume<0.0)
    {
        v = 0.0;
    }
    else{
        v = volume;
    }
}

void RK::SetGrowth(double growth)
{
    if(growth>1000.0||growth<-maximum_virulence)
    {
        cout<<"growth is inappropriate"<<endl;
        cout<<"growth="<<growth<<endl;
        exit(1);
    }
    else
    {
        growth_rate = growth;
    }
}

void RK::SetK_delta(int onoff,int ith)
{
    if(onoff<0.0||onoff>1.0)
    {
        cout<<"input is inappropriate"<<endl;
        cout<<"input="<<onoff<<endl;
        ShowDate();
        exit(1);
    }
    else
    {
        k_delta[ith] = onoff;
    }
}

void RK::SetInput_list(int onoff,int ith,int lth)//inputから繋がる先　onoffには遺伝子番号がはいる
{
    if(onoff<0||onoff>ProteinVAL)
    {
        cout<<"input 1is inappropriate"<<endl;
        cout<<"input="<<onoff<<endl;
        ShowDate();
        exit(1);
    }
    else
    {
        input_list[ith][lth] = onoff;
    }
}

void RK::SetTargetPattern(int onoff,int lth)
{
    if(onoff<0||onoff>1)
    {
        cout<<"pattern is inappropriate"<<endl;
        cout<<"onoff="<<onoff<<endl;
        ShowDate();
        exit(1);
    }
    else
    {
        TargetPattern[lth] = onoff;
    }
}

void RK::SetInfection_virulence(double vir,int lth)
{
    if(vir<0||vir>maximum_virulence*MAXCELLNUMBER*100)//パラサイトの個体数がホストの100倍になったら計算を終了
    {
        cout<<"vir is inappropriate"<<endl;
        cout<<"vir="<<vir<<endl;
        exit(1);
    }
    else
    {
        infection_virulence[lth] = vir;
    }
}

void RK::SetQ(double probability)
{
    if(probability<0.0||probability>1.1)
    {
        cout<<"probability is inappropriate"<<endl;
        cout<<"probability="<<probability<<endl;
        exit(1);
    }
    q = probability;
}

void RK::SetCelltype(int type)
{
    celltype = type;
}

void RK::SetAverage_growth(double averagegrowth)
{
    if(averagegrowth<0.0||averagegrowth>100)
    {
        cout<<"averagegrowth is inappropriate"<<endl;
        cout<<"averagegrowth="<<averagegrowth<<endl;
        //ShowDate();
        exit(1);
    }
    else
    {
        average_growth = averagegrowth;
    }
}

void RK::SetDivision_number(int numb)
{
    division_number = numb;
}

void RK::SetMutation_number(int numb)
{
    mutaion_number = numb;
}

void RK::SetMu(double _mu)
{
    if(_mu<-0.1||_mu>1000)
    {
        cout<<"_mu="<<_mu<<endl;
        exit(-1);
    }
    else
    {
        mu = _mu;
    }
}

void RK::SetHammingD(double ha)
{
    if(ha<0)
    {
        cout<<"hamming="<<ha;
        exit(-1);
    }
    else{
        hammingD = ha;
    }
}

void RK::SetRelativeHammingD(double Rha)
{
    if(Rha<0)
    {
        cout<<"Rhamming="<<Rha;
        exit(-1);
    }
    else{
        RelativehammingD = Rha;
    }
}

void RK::SetStationaryResult(double x,int lth)
{
    for(int j=0;j<ProteinVAL;j++)
    {
        StationaryResult[lth] = x;
    }
}

void RK::SetStationaryTime(double sta_time)
{
    StationaryTime = sta_time;
}

void RK::SetStationaryGrowth(double growth_sta)
{
    StationaryGrowth = growth_sta;
}
void RK::SetFitness(double F)
{
    if(F<0.0)
    {
        cout<<"fitness is "<<F<<endl;
        exit(-1);
    }
    else
    {
        fitness = F;
    }
}
void RK::SetFitness_p(double F_p)
{
    if(F_p<0.0)
    {
        cout<<"fitness_p is "<<F_p<<endl;
        exit(-1);
    }
    else
    {
        fitness_p = F_p;
    }
}
void RK::SetCumulativ_p(double cp)
{
    if(cp<0.0)
    {
        cout<<"cp is "<<cp<<endl;
    }
    else{
        cumulative_p =cp;
    }
}
void RK::SetPara_fitness(double p_fitness)
{
    para_fitness = p_fitness;
}

void RK::SetFitness_Target_parasite_average(double tpa_fitness)
{
    fitness_Target_parasite_average =tpa_fitness;
}

void RK::SetStack_fitness_Target_parasite_average(double stpa_fitness)
{
    stack_fitness_Target_parasite_average =stpa_fitness;
}

void RK::SetFitness_Target_noparasite_average(double tpa_fitness)
{
    fitness_Target_noparasite_average =tpa_fitness;
}
void RK::SetStack_fitness_Target_noparasite_average(double stpa_fitness)//平均値計算用の前世代の入れ物
{
    stack_fitness_Target_noparasite_average =stpa_fitness;
}

void RK::SetSigma_Str(double str)
{
    if(str>0)
    {
        sigma_str =str;
    }
    else{
        sigma_str =-str;
    }
}
inline double RK::GetResult(int lth) const
{
    return result[lth];
}

inline double RK::GetResult_average(int lth) const
{
    return result_average[lth];
}
inline double RK::GetResult_check(int lth,int check_t) const
{
    return result_check[lth][check_t];
}

inline double RK::GetResult_VipVg(int lth) const//result[lth] protein_x
{
    return result_VipVg[lth];
}
inline double RK::GetStationaryResult(int lth) const
{
    return StationaryResult[lth];
}

inline double RK::GetJ(int ith,int lth) const//from l to i
{
    return J[ith][lth];
}

inline int RK::GetJ_VipVg(int ith,int lth) const//from l to i
{
    return J_VipVg[ith][lth];
}

inline int RK::GetStationaryJ(int ith,int lth) const
{
    return StationaryJ[ith][lth];
}

inline int RK::GetEL(int ith,int lth) const//from l to i
{
    return EL[ith][lth];
}

inline int RK::GetWEL(int ith,int lth) const//from l to i
{
    return WEL[ith][lth];
}

inline int RK::GetELSize(int ith) const
{
    return (int)EL[ith].size();
}

inline int RK::GetELSize_VipVg(int ith) const
{
    return (int)EL_VipVg[ith].size();
}

inline int RK::GetSEL(int ith) const
{
    return SEL[ith];
}

inline int RK::GetREL(int ith,int lth) const
{
    return REL[ith][lth];//ithから出るpath
}

inline int RK::GetRELSize(int ith) const
{
    return (int)REL[ith].size();
}

inline int RK::GetRSEL(int ith) const
{
    return RSEL[ith];
}


inline int RK::GetEL_VipVg(int ith,int lth) const//from l to i
{
    return EL_VipVg[ith][lth];
}

inline int RK::GetSEL_VipVg(int ith) const
{
    return SEL_VipVg[ith];
}

inline int RK::GetWEL_VipVg(int ith,int lth) const//from l to i
{
    return WEL_VipVg[ith][lth];
}

inline int RK::GetREL_VipVg(int ith,int lth) const
{
    return REL_VipVg[ith][lth];//ithから出るpath
}

inline int RK::GetRELSize_VipVg(int ith) const
{
    return (int)REL_VipVg[ith].size();
}

inline int RK::GetRSEL_VipVg(int ith) const
{
    return RSEL_VipVg[ith];
}

int RK::GetK_delt(int lth)
{
    return k_delta[lth];
}

int RK::GetInput_list(int ith,int lth)
{
    return input_list[ith][lth];
}

int RK::GetInput_listSize(int ith)
{
    return (int)input_list[ith].size();//depth =ithのgeneから出るpathの本数
}

int RK::GetTargetPattern(int lth)
{
    return TargetPattern[lth];
}

double RK::GetInfection_virulence(int lth)
{
    return infection_virulence[lth];
}

double RK::GetDelt() const
{
    return delt;
}

double RK::GetProD(int lth) const
{
    return proteinD[lth];
}

double RK::GetThreshold_g(int lth) const
{
    return threshold_g[lth];
}

int RK::GetNum() const
{
    return number;
}
int RK::GetCelltype() const
{
    return celltype;
}
double RK::GetTime() const
{
    return time;
}

int RK::GetGene() const
{
    return gene;
}

int RK::GetSorting_num() const
{
    return sorting_number;
}

double RK::GetMutationTime() const
{
    return mutationTime;
}
double RK::GetBeta() const
{
    return beta;
}

double RK::GetDiviTime() const
{
    return divisionTime;
}

int RK::GetTime1000() const
{
    return time1000;
}
int RK::GetMutationTime1000() const
{
    return matationTime1000;
}
int RK::GetDiviTimeTime1000() const
{
    return divisionTime1000;
}

int RK::Getprotein_x_number() const
{
    return protein_x_number;
}

double RK::GetMeanfieldval(int lth) const
{
    return Meanfieldval[lth];
}

int RK::GetTotalcellnumber() const
{
    return totalcellnumber;
}

double RK::GetVolume() const
{
    return v;
}

double RK::GetInitialVolume() const
{
    return initial_volume;
}
double RK::GetGrowthRate() const
{
    return growth_rate;
}

double RK::GetQ() const
{
    return q;
}

double RK::GetAverage_Growth() const
{
    return average_growth;
}

int RK::GetMutation_number() const
{
    return mutaion_number;
}

double RK::GetHammingD() const
{
    return hammingD;
}

double RK::GetRelativeHammingD() const
{
    return RelativehammingD;
}

int RK::GetDivision_number() const
{
    return division_number;
}

double RK::GetStationaryTime() const
{
    return StationaryTime;
}

double RK::GetStationaryGrowth() const
{
    return StationaryGrowth;
}

double RK::GetMu() const
{
    return mu;
}

double RK::GetFitness() const
{
    return fitness;
}
double RK::GetFitness_p() const//fitnessから得た確率
{
    return fitness_p;
}
double RK::GetCumulative_p() const
{
    return cumulative_p;
}
double RK::GetFinishTime() const
{
    return finishTime;
}

double RK::GetPara_fitness() const
{
    return para_fitness;
}

double RK::GetFitness_Target_parasite_average() const
{
    return fitness_Target_parasite_average;
}

double RK::GetStack_fitness_Target_parasite_average() const
{
    return stack_fitness_Target_parasite_average;
}

double RK::GetFitness_Target_noparasite_average() const
{
    return fitness_Target_noparasite_average;
}
double RK::GetStack_fitness_Target_noparasite_average() const
{
    return stack_fitness_Target_noparasite_average;
}

double RK::GetSigma_Str() const
{
    return sigma_str;
}

void RK::MutationJ(int ith,int lth)
{
    int time=0;
    int nowtime=0,muttime=0;
    int rnd=0;
    nowtime = GetTime1000();
    muttime = GetMutationTime1000();
    
    time = nowtime%muttime;
    if(time==0)
    {
        cout<<"time="<<GetTime()<<endl;
        cout<<"nowtime1000="<<nowtime<<endl;
        cout<<"muttime1000="<<muttime<<endl;
        ShowDate();
        for(int i=0;i<result.size();i++)
        {
            for(int l=0;l<result.size();l++)
            {
                rnd = rand();
                //cout<<"rnd="<<rnd<<endl;
                //SetJ(rnd, i, l);
            }
        }
        ShowDate();
    }
}

double RK::MutationG(int ith,double sigma)//ith成分
{
    double parent_g =0.0,child_g=0.0;
    parent_g = GetThreshold_g(ith);
    //cout<<"parent_g="<<parent_g<<endl;
    child_g  = gauss_rand(parent_g, sigma);
    //cout<<"child_g="<<child_g<<endl;
    
    return child_g;
}

double RK::MutationD(int ith,double sigma)
{
    double parent_D = 0.0,child_D = 0.0;
    parent_D = GetProD(ith);
    cout<<"parent_D="<<parent_D<<endl;
    child_D  = gauss_rand(parent_D, sigma);
    cout<<"child_D="<<child_D<<endl;
    
    return child_D;
}

inline void RK::calculateNS_protein()
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0,h=0;
    int i;
    
    
    double xa[protein_x_number];
    double k1a[protein_x_number];
    double k2a[protein_x_number];
    double k3a[protein_x_number];
    double k4a[protein_x_number];
    double tmpa[protein_x_number];
    
    //vector<double> xa;
    //vector<double> k1a;
    //vector<double> k2a;
    //vector<double> k3a;
    //vector<double> k4a;
    //vector<double> tmpa;
    //xa.reserve(protein_x_number);
    //k1a.reserve(protein_x_number);
    //k2a.reserve(protein_x_number);
    //k3a.reserve(protein_x_number);
    //k4a.reserve(protein_x_number);
    //tmpa.reserve(protein_x_number);
    //x = new double[protein_x_number];
    //k1 = new double[protein_x_number];
    //k2 = new double[protein_x_number];
    //k3 = new double[protein_x_number];
    //k4 = new double[protein_x_number];
    //tmp = new double[protein_x_number];
    
    t=GetTime();
    h=GetDelt();
    
    for(i=0;i<ProteinVAL;i++)
    {
        //x[i] = GetResult(i);
        xa[i] = GetResult(i);
        
        //cout<<"x["<<i<<"]="<<x[i]<<endl;
        //cout<<"result["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k1[i] = func_protein(t,x,i);
        //tmp[i] = x[i] + h*k1[i]/2;
        k1a[i] = func_protein(t,xa,i);
        tmpa[i] =xa[i] + h*k1a[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<ProteinVAL;i++)
    {
        //k2[i]= func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k2[i]/2;
        k2a[i] = func_protein(t+h/2,tmpa,i);
        tmpa[i] =xa[i] + h*k2a[i]/2;
        
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k3[i] = func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k3[i];
        k3a[i] = func_protein(t+h/2,tmpa,i);
        tmpa[i] = xa[i] + h*k3a[i];
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k4[i] = func_protein(t+h,tmp,i);
        //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        k4a[i] = func_protein(t+h,tmpa,i);
        xa[i] = xa[i] + (k1a[i]+2*k2a[i]+2*k3a[i]+k4a[i])*h/6;
        SetResult(xa[i], i);
    }
    
    //delete [] x;
    //delete [] k1;
    //delete [] k2;
    //delete [] k3;
    //delete [] k4;
    //delete [] tmp;
}

inline double RK::func_protein(double t,double *x,int i)//ith proteinの関数
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g=0.0,d=0.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=0,n=GetSEL(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL(i,l)*x[GetEL(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    //z = z + S*k_delta[i] ;//input シグナル無限に入って来ているnonzeroなら
    
    //for(int l=0;l<input_number;l++)
    //{//z = z + GetYin(l)*GetYin_J(l,i);//from ith input成分 to lth protein
    //cout<<"z["<<l<<"]="<<z<<endl;
    //}
    //d = 1+exp_fast(-beta_*(z-g));//分母
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    //r=sigmoid -x[i]+ diff*(MF-x[i]) +eps;
    //r=gamma*(sigmoid -x[i])+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    //r=sigmoid -x[i]+ eps ;
    r=gamma*(sigmoid -x[i] + eps);//相互作用なし
    //r=sigmoid -x[i];//相互作用なし
    //r=sigmoid*v -x[i]*v+ diff*(MF-x[i]) ;//希釈あり 全体に
    //r=sigmoid -x[i]*v + diff*(MF-x[i]) ;//希釈あり
    //r=sigmoid -x[i]+ diff*(GetX(i)-x[i]) ;//希釈なし
    //r=sigmoid -x[i]*_mu+ diff*(GetX(i)-x[i]) +eps;
    
    return r;
}

inline double RK::func_protein_continuas(double t,double *x,int i)//ith proteinの関数
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g=0.0,d=0.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=target_number,n=ProteinVAL;l<n;l++)
    {
        z = z + GetJ(i,l)*x[l];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    r=gamma*(sigmoid -x[i] + eps);//相互作用なし
    
    return r;
}

inline void RK::calculateNS_protein_MF()
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0,h=0;
    int i;
    
    
    double xa[protein_x_number];
    double k1a[protein_x_number];
    double k2a[protein_x_number];
    double k3a[protein_x_number];
    double k4a[protein_x_number];
    double tmpa[protein_x_number];
    
    //vector<double> xa;
    //vector<double> k1a;
    //vector<double> k2a;
    //vector<double> k3a;
    //vector<double> k4a;
    //vector<double> tmpa;
    //xa.reserve(protein_x_number);
    //k1a.reserve(protein_x_number);
    //k2a.reserve(protein_x_number);
    //k3a.reserve(protein_x_number);
    //k4a.reserve(protein_x_number);
    //tmpa.reserve(protein_x_number);
    //x = new double[protein_x_number];
    //k1 = new double[protein_x_number];
    //k2 = new double[protein_x_number];
    //k3 = new double[protein_x_number];
    //k4 = new double[protein_x_number];
    //tmp = new double[protein_x_number];
    
    t=GetTime();
    h=GetDelt();
    
    for(i=0;i<ProteinVAL;i++)
    {
        //x[i] = GetResult(i);
        xa[i] = GetResult(i);
        //cout<<"x["<<i<<"]="<<x[i]<<endl;
        //cout<<"result["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k1[i] = func_protein(t,x,i);
        //tmp[i] = x[i] + h*k1[i]/2;
        k1a[i] = func_protein_MF(t,xa,i);
        tmpa[i] =xa[i] + h*k1a[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<ProteinVAL;i++)
    {
        //k2[i]= func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k2[i]/2;
        k2a[i] = func_protein_MF(t+h/2,tmpa,i);
        tmpa[i] =xa[i] + h*k2a[i]/2;
        
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k3[i] = func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k3[i];
        k3a[i] = func_protein_MF(t+h/2,tmpa,i);
        tmpa[i] = xa[i] + h*k3a[i];
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k4[i] = func_protein(t+h,tmp,i);
        //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        k4a[i] = func_protein_MF(t+h,tmpa,i);
        xa[i] = xa[i] + (k1a[i]+2*k2a[i]+2*k3a[i]+k4a[i])*h/6;
        SetResult(xa[i], i);
    }
    
    //delete [] x;
    //delete [] k1;
    //delete [] k2;
    //delete [] k3;
    //delete [] k4;
    //delete [] tmp;
}

inline double RK::func_protein_MF(double t,double *x,int i)
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g=0.0,d=0.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=0,n=GetSEL(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL(i,l)*x[GetEL(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    //z = z + S*k_delta[i] ;//input シグナル無限に入って来ているnonzeroなら
    
    //for(int l=0;l<input_number;l++)
    //{//z = z + GetYin(l)*GetYin_J(l,i);//from ith input成分 to lth protein
    //cout<<"z["<<l<<"]="<<z<<endl;
    //}
    //d = 1+exp_fast(-beta_*(z-g));//分母
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    r=sigmoid -x[i]+ diff*(MF-x[i]);
    //r=gamma*(sigmoid -x[i])+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    //r=sigmoid -x[i]+ eps ;
    //r=gamma*(sigmoid -x[i] + eps);//相互作用なし
    //r=sigmoid -x[i];//相互作用なし
    //r=sigmoid*v -x[i]*v+ diff*(MF-x[i]) ;//希釈あり 全体に
    //r=sigmoid -x[i]*v + diff*(MF-x[i]) ;//希釈あり
    //r=sigmoid -x[i]+ diff*(GetX(i)-x[i]) ;//希釈なし
    //r=sigmoid -x[i]*_mu+ diff*(GetX(i)-x[i]) +eps;
    
    return r;
}

inline void RK::calculateNS_protein_g()
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0,h=0;
    int i;
    
    double a,z;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    
    double xa[protein_x_number];
    double k1a[protein_x_number];
    double k2a[protein_x_number];
    double k3a[protein_x_number];
    double k4a[protein_x_number];
    double tmpa[protein_x_number];
    double g_rand[protein_x_number];
    
    //vector<double> xa;
    //vector<double> k1a;
    //vector<double> k2a;
    //vector<double> k3a;
    //vector<double> k4a;
    //vector<double> tmpa;
    //xa.reserve(protein_x_number);
    //k1a.reserve(protein_x_number);
    //k2a.reserve(protein_x_number);
    //k3a.reserve(protein_x_number);
    //k4a.reserve(protein_x_number);
    //tmpa.reserve(protein_x_number);
    //x = new double[protein_x_number];
    //k1 = new double[protein_x_number];
    //k2 = new double[protein_x_number];
    //k3 = new double[protein_x_number];
    //k4 = new double[protein_x_number];
    //tmp = new double[protein_x_number];
    
    t=GetTime();
    h=GetDelt();
    
    for(i=0;i<ProteinVAL;i++)
    {
        //x[i] = GetResult(i);
        xa[i] = GetResult(i);
        g_rand[i] = Sigma*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));
        //g_rand[i] = 0.0;
        
        //cout<<"x["<<i<<"]="<<x[i]<<endl;
        //cout<<"result["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k1[i] = func_protein(t,x,i);
        //tmp[i] = x[i] + h*k1[i]/2;
        k1a[i] = func_protein_g(t,xa,g_rand,i);
        tmpa[i] =xa[i] + h*k1a[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<ProteinVAL;i++)
    {
        //k2[i]= func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k2[i]/2;
        k2a[i] = func_protein_g(t+h/2,tmpa,g_rand,i);
        tmpa[i] =xa[i] + h*k2a[i]/2;
        
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k3[i] = func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k3[i];
        k3a[i] = func_protein_g(t+h/2,tmpa,g_rand,i);
        tmpa[i] = xa[i] + h*k3a[i];
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k4[i] = func_protein(t+h,tmp,i);
        //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        k4a[i] = func_protein_g(t+h,tmpa,g_rand,i);
        xa[i] = xa[i] + (k1a[i]+2*k2a[i]+2*k3a[i]+k4a[i])*h/6;
        SetResult(xa[i], i);
    }
    
    //delete [] x;
    //delete [] k1;
    //delete [] k2;
    //delete [] k3;
    //delete [] k4;
    //delete [] tmp;
}

inline double RK::func_protein_g(double t,double *x,double *g_rand,int i)//ith proteinの関数
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g= 0.0,d=0.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=0,n=GetSEL(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL(i,l)*x[GetEL(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    //z = z + S*k_delta[i] ;//input シグナル無限に入って来ているnonzeroなら
    
    //for(int l=0;l<input_number;l++)
    //{//z = z + GetYin(l)*GetYin_J(l,i);//from ith input成分 to lth protein
    //cout<<"z["<<l<<"]="<<z<<endl;
    //}
    //d = 1+exp_fast(-beta_*(z-g));//分母
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    //r=sigmoid -x[i]+ diff*(MF-x[i]) +eps;
    //r=sigmoid -x[i]+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    r=gamma*(sigmoid +eps -x[i]) +g_rand[i];
    //r=sigmoid -x[i];//相互作用なし
    //r=sigmoid*v -x[i]*v+ diff*(MF-x[i]) ;//希釈あり 全体に
    //r=sigmoid -x[i]*v + diff*(MF-x[i]) ;//希釈あり
    //r=sigmoid -x[i]+ diff*(GetX(i)-x[i]) ;//希釈なし
    //r=sigmoid -x[i]*_mu+ diff*(GetX(i)-x[i]) +eps;
    
    return r;
}

const double SH = sqrt(0.05/(double)(EM_h));//EM_h=5:SH=0.1 EM_h=10:SH=0.07 EM_h=20:SH=0.05 EM_h=50:SH=0.033
const double DELT_H = 0.05/(double)(EM_h);//EM_h=5:DELT_H=0.01 EM_h=10:DELT_H=0.005 EM_h=20:DELT_H=0.0025 EM_h=50:DELT_H=0.001

inline void RK::calculateNS_protein_g_EM(vector<vector< double> > g)
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0.0,h=0.0,sh=0.0;
    int i,j;
    
    double xa[protein_x_number];
    double k4a[protein_x_number];
    
    t=GetTime();
    h = DELT_H;
    sh = SH;
    
    for(j=0;j<EM_h;j++)
    {
        for(i=0;i<ProteinVAL;i++)
        {
            xa[i] = GetResult(i);
            //g_rand[i] = 0.0;
            
            //cout<<"x["<<i<<"]="<<x[i]<<endl;
            //cout<<"result["<<i<<"]="<<result[i]<<endl;
        }
        
        for(i=0;i<ProteinVAL;i++)
        {
            //k4[i] = func_protein(t+h,tmp,i);
            //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
            k4a[i] = func_protein(t+h,xa,i);
            //k4a[i] = func_protein_continuas(t+h,xa,i);
            //k4a[i] = func_protein_delution(t+h,xa,i);
            xa[i] = xa[i] + k4a[i]*h + g[j][i]*sh;
            SetResult(xa[i], i);
        }
    }
}

inline void RK::calculateNS_protein_delution_g_EM(vector<vector< double> > g)
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0.0,h=0.0,sh=0.0;
    int i,j;
    
    double xa[protein_x_number];
    double k4a[protein_x_number];
    
    t=GetTime();
    h = DELT_H;
    sh = SH;
    
    for(j=0;j<EM_h;j++)
    {
        for(i=0;i<ProteinVAL;i++)
        {
            xa[i] = GetResult(i);
        }
        for(i=0;i<ProteinVAL;i++)
        {
            k4a[i] = func_protein_delution(t+h,xa,i);
            xa[i] = xa[i] + k4a[i]*h + g[j][i]*sh;
            SetResult(xa[i], i);
        }
    }
}

void RK::calculateNS_protein_g_EM_VipVg(vector<vector< double> > g)//VipVg計算用　通常の計算とは分けておく
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0.0,h=0.0,sh=0.0;
    int i,j;
    
    double xa[protein_x_number];
    double k4a[protein_x_number];
    
    t=GetTime();
    h = DELT_H;
    sh = SH;
    
    for(j=0;j<EM_h;j++)
    {
        for(i=0;i<ProteinVAL;i++)
        {
            xa[i] = GetResult_VipVg(i);
            
            //cout<<"x["<<i<<"]="<<x[i]<<endl;
            //cout<<"result["<<i<<"]="<<result[i]<<endl;
        }
        
        for(i=0;i<ProteinVAL;i++)
        {
            //k4[i] = func_protein(t+h,tmp,i);
            //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
            k4a[i] = func_protein_VipVg(t+h,xa,i);
            //k4a[i] = func_protein_delution(t+h,xa,i);
            xa[i] = xa[i] + k4a[i]*h + g[j][i]*sh;
            SetResult_VipVg(xa[i], i);
        }
    }
}

inline double RK::func_protein_VipVg(double t,double *x,int i)//ith proteinの関数
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g=0.0,d=0.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=0,n=GetSEL_VipVg(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL_VipVg(i,l)*x[GetEL_VipVg(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    
    z = z/Route;
    //z = z + S*k_delta[i] ;//input シグナル無限に入って来ているnonzeroなら
    
    //for(int l=0;l<input_number;l++)
    //{//z = z + GetYin(l)*GetYin_J(l,i);//from ith input成分 to lth protein
    //cout<<"z["<<l<<"]="<<z<<endl;
    //}
    //d = 1+exp_fast(-beta_*(z-g));//分母
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    //r=sigmoid -x[i]+ diff*(MF-x[i]) +eps;
    //r=gamma*(sigmoid -x[i])+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    //r=sigmoid -x[i]+ eps ;
    r=gamma*(sigmoid -x[i] + eps);//相互作用なし
    //r=sigmoid -x[i];//相互作用なし
    //r=sigmoid*v -x[i]*v+ diff*(MF-x[i]) ;//希釈あり 全体に
    //r=sigmoid -x[i]*v + diff*(MF-x[i]) ;//希釈あり
    //r=sigmoid -x[i]+ diff*(GetX(i)-x[i]) ;//希釈なし
    //r=sigmoid -x[i]*_mu+ diff*(GetX(i)-x[i]) +eps;
    
    return r;
}

void RK::calculateNS_protein_signal()
{
    double t=0,h=0;
    int i;
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    
    double xa[protein_x_number];
    double k1a[protein_x_number];
    double k2a[protein_x_number];
    double k3a[protein_x_number];
    double k4a[protein_x_number];
    double tmpa[protein_x_number];
    double g_rand[protein_x_number];
    
    t=GetTime();
    h=GetDelt();
    
    for(i=0;i<ProteinVAL;i++)
    {
        //x[i] = GetResult(i);
        xa[i] = GetResult(i);
        g_rand[i] = Sigma*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));
        //g_rand[i] = 0.0;
        
        //cout<<"x["<<i<<"]="<<x[i]<<endl;
        //cout<<"result["<<i<<"]="<<result[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k1[i] = func_protein(t,x,i);
        //tmp[i] = x[i] + h*k1[i]/2;
        k1a[i] = func_protein_signal(t,xa,i);
        tmpa[i] =xa[i] + h*k1a[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<ProteinVAL;i++)
    {
        //k2[i]= func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k2[i]/2;
        k2a[i] = func_protein_signal(t+h/2,tmpa,i);
        tmpa[i] =xa[i] + h*k2a[i]/2;
        
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k3[i] = func_protein(t+h/2,tmp,i);
        //tmp[i] = x[i] + h*k3[i];
        k3a[i] = func_protein_signal(t+h/2,tmpa,i);
        tmpa[i] = xa[i] + h*k3a[i];
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<ProteinVAL;i++)
    {
        //k4[i] = func_protein(t+h,tmp,i);
        //x[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        k4a[i] = func_protein_signal(t+h,tmpa,i);
        xa[i] = xa[i] + (k1a[i]+2*k2a[i]+2*k3a[i]+k4a[i])*h/6;
        SetResult(xa[i], i);
    }
}

double RK::func_protein_signal(double t,double *x,int i)
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,g= 0.0,d=0.0,S=10.0;//S is input factor
    double MF =0.0,diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    //double v = GetGrowthRate();
    //double _mu = 0.5+0.5*GetMu()/10.0;
    
    for(int l=0,n=GetSEL(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL(i,l)*x[GetEL(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    //z = z + S*k_delta[i] ;//input シグナル無限に入って来ているnonzeroなら
    
    if(i==ProteinVAL-input_number-1)
    {
        z = z - S*GetParasite_population();
    }
    
    //for(int l=0;l<input_number;l++)
    //{//z = z + GetYin(l)*GetYin_J(l,i);//from ith input成分 to lth protein
    //cout<<"z["<<l<<"]="<<z<<endl;
    //}
    //d = 1+exp_fast(-beta_*(z-g));//分母
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    //r=sigmoid -x[i]+ diff*(MF-x[i]) +eps;
    //r=sigmoid -x[i]+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    //r=sigmoid -x[i]+ eps ;
    r=sigmoid -x[i];//相互作用なし
    //r=sigmoid*v -x[i]*v+ diff*(MF-x[i]) ;//希釈あり 全体に
    //r=sigmoid -x[i]*v + diff*(MF-x[i]) ;//希釈あり
    //r=sigmoid -x[i]+ diff*(GetX(i)-x[i]) ;//希釈なし
    //r=sigmoid -x[i]*_mu+ diff*(GetX(i)-x[i]) +eps;
    
    return r;
}

double RK::func_protein_delution(double t,double *x,int i)
{
    double z=0.0,r=0.0,sigmoid=0.0,beta_=0.0,gamma=0.1,g= 0.0,d=0.0;//S is input factor
    double diff= 0.0;
    diff = GetProD(i);
    beta_ = GetBeta();
    g = GetThreshold_g(i);
    
    //MF = GetMeanfieldval(i);
    
    double v = GetGrowthRate();
    double mu = ((v/growth_const) + maximum_virulence)/(double)(target_number);
    
    for(int l=0,n=GetSEL(i);l<n;l++)
    {
        //cout<<GetEL(i,l)<<endl;
        z = z + GetWEL(i,l)*x[GetEL(i,l)];//from l to i EL[i][l]*x[EL[i][l]]
    }
    z = z/Route;
    
    d = 1.0+exp(-beta_*(z-g));//分母
    sigmoid = 1.0/d;//f(xi)=exp(-beta*(z-thita))
    
    //r=sigmoid -x[i]+ diff*(MF-x[i]) +eps;
    //r=sigmoid -x[i]+ diff*(MF-x[i]);
    //r=sigmoid -x[i]+ diff*(MF-x[i]) + eps +g_rand[i];
    
    //r=sigmoid -x[i]+ eps ;
    //r=sigmoid -x[i];//相互作用なし
    r=gamma*mu*(sigmoid -x[i]);//相互作用なし　希釈あり
    //r=sigmoid*v -x[i]*v ;//希釈あり 全体に
    
    return r;
}

void RK::calculateNS_protein_signal_g_EM(vector<vector< double> > g)
{
    //double *x, *k1,*k2,*k3,*k4,*tmp;
    double t=0.0,h=0.0,sh=0.0;
    int i,j;
    
    double xa[protein_x_number];
    double k4a[protein_x_number];
    
    t=GetTime();
    h = DELT_H;
    sh = SH;
    
    for(j=0;j<EM_h;j++)
    {
        for(i=0;i<ProteinVAL;i++)
        {
            xa[i] = GetResult(i);
        }
        for(i=0;i<ProteinVAL;i++)
        {
            //k4a[i] = func_protein_signal(t+h,xa,i);signalが入る時
            k4a[i] = func_protein_delution(t+h,xa,i);//成長で希釈される時
            xa[i] = xa[i] + k4a[i]*h + g[j][i]*sh;
            SetResult(xa[i], i);
        }
    }
}

int  RK::rand()//-1,0,1の乱数発生
{
    random_device rnd;
    int rnt,rnf;
    
    rnt = rnd();
    rnt =abs(rnt);
    rnf = (int)(rnt%3-1);
    return rnf;
}

int RK::rand_select(int a,int b,int c)//a:１の割合、b:0の割合 c:-1の割合
{
    int rnt,rnf;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> randProtein(1, 100);
    
    rnt = randProtein(mt);
    if(rnt>=1 && rnt<=a)
    {
        rnf = 1;
    }
    else if(rnt>=a && rnt<=a+b)
    {
        rnf = 0;
    }
    else
    {
        rnf = -1;
    }
    
    return rnf;
}

int  RK::rand_onoff()//0,1の乱数発生
{
    random_device rnd;
    int rnt,rnf;
    
    rnt = rnd();
    rnt =abs(rnt);
    rnf = (int)(rnt%2);
    return rnf;
}

//標準分布に従う乱数を返す関数
double RK::gauss_rand(double mu,double sigma)
{
    double a,z;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    z = sqrt(-2.0*log(rand01(mt)) )*sin(2.0*PI*rand01(mt) );
    a = mu + sigma*z;
    return a;
}

double RK::gabs_gauss_rand(double mu,double sigma)
{
    double a,z;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    z = sqrt(-2.0*log(rand01(mt)) )*sin(2.0*PI*rand01(mt) );
    a = mu + sigma*z;//絶対値をとり足す
    if(a<0.0)
    {
        a = mu;
    }
    return a;
}
void RK::failclean()
{
    char filename[50]={};
    char filenamer[50]={};
    char filenameg[50] ={};
    
    int _number=0;
    _number  = GetNum();
    
    sprintf(filename,"CELL%03d.txt",_number);
    sprintf(filenamer,"GROWTH%03d.txt",_number);
    sprintf(filenameg,"Threshold_g%03d.txt",_number);
    remove(filename);
    remove(filenamer);
    remove(filenameg);
}

void RK::failclean_evo()
{
    char filename[50]={};
    char filenamer[50]={};
    char filenameg[50] ={};
    char filenamepg[50] ={};
    char filenamepd[50] ={};
    char filenamesta[50] = {};
    char filenamersta[50] ={};
    
    int _number=0;
    _number  = GetNum();
    
    sprintf(filename,"cell%04d.txt",_number);
    sprintf(filenamer,"growth%04d.txt",_number);
    sprintf(filenameg,"threshold_g%04d.txt",_number);
    sprintf(filenamepg,"p-g%04d.txt",_number);
    sprintf(filenamepd,"proteinD%04d.txt",_number);
    sprintf(filenamesta,"cell_stationary%04d.txt",_number);
    sprintf(filenamersta,"growth_stationary%04d.txt",_number);
    
    remove(filename);
    remove(filenamer);
    remove(filenameg);
    remove(filenamepg);
    remove(filenamepd);
}

void RK::calculateNS_volume()
{
    double x=0.0,k1=0.0,k2=0.0,k3=0.0,k4=0.0,tmp=0.0;
    double t=0.0,h=0.0;
    
    t=GetTime();
    h=GetDelt();
    x = GetVolume();
    
    k1 = func_volume(t, x);
    tmp = x + h*k1/2;
    k2= func_volume(t+h/2,tmp);
    tmp = x + h*k2/2;
    k3 = func_volume(t+h/2,tmp);
    tmp = x + h*k3;
    k4 = func_volume(t+h,tmp);
    x = x + (k1+2*k2+2*k3+k4)*h/6;
    SetVolume(x);
}

double RK::func_volume(double t,double x)
{
    double r=0.0,growth=0.0;
    double c = 1.0 ;//volumeの大きさを調整--->1に固定してgrowthを調整
    growth = GetGrowthRate();
    
    return r = c*growth*x;
}

void RK::calculate_Mu()
{
    double sigma=0.0;
    for(int l=0;l<ProteinVAL;l++)
    {
        sigma = sigma + GetResult(l);//発現量の和をコストとする
    }
    
    SetMu(sigma);
}


//mediumの関数
void RK::SetY(double Y_med,int ith)//input物質が存在するmedium
{
    if(Y_med<-0.1)
    {
        
        cout<<"Y_mid="<<Y_med<<endl;
        exit(-1);
    }
    else{
        Y[ith] = Y_med;
    }
}
double RK::GetY(int ith) const
{
    return Y[ith];
}

void RK::SetX(double X_med,int ith)//拡散物質が存在するmedium
{
    if(X_med<-0.1)
    {
        
        cout<<"Y_mid="<<X_med<<endl;
        exit(-1);
    }
    else{
        X[ith] = X_med;
    }
}

double RK::GetX(int ith) const
{
    return X[ith];
}

void RK::SetYin(double Y_in,int ith)//cell中のinput物質
{
    if(Y_in<-0.1)
    {
        cout<<"Y_in="<<Y_in<<endl;
        exit(-1);
    }
    else{
        Yin[ith] = Y_in;
        
    }
}
double RK::GetYin(int ith) const
{
    return Yin[ith];
}

void RK::SetMedium_D(double Y_input_D,int ith)
{
    Medium_D[ith] = Y_input_D;
}

double RK::GetMedium_D(int ith) const
{
    return Medium_D[ith];
}

int RK::GetYin_J(int ith,int lth) const//input物質のYin_Jのgetter
{
    return Yin_J[ith][lth];
}

void RK::SetYin_J(int yin_j,int ith,int lth)//input物質のYin_Jのsetter
{
    if(yin_j>1)
    {
        cout<<"Yib_J=";
        cout<<Yin_J[ith][lth]<<endl;
        exit(-1);
    }
    else
    {
        Yin_J[ith][lth] = yin_j;
    }
}

void RK::calculateNS_Yin()
{
    double *k1,*k2,*k3,*k4,*tmp;
    double t=0.0,h=0.0;
    double *y ;
    int i;
    
    y = new double[input_number];
    k1 = new double[input_number];
    k2 = new double[input_number];
    k3 = new double[input_number];
    k4 = new double[input_number];
    tmp = new double[input_number];
    
    t=GetTime();
    h=GetDelt();
    
    for(i=0;i<input_number;i++)
    {
        y[i] = GetYin(i);
        //cout<<"y["<<i<<"]="<<y[i]<<endl;
    }
    for(i=0;i<input_number;i++)
    {
        k1[i] = func_Yin(t,y,i);
        tmp[i] = y[i] + h*k1[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<input_number;i++)
    {
        k2[i]= func_Yin(t+h/2,tmp,i);
        
        tmp[i] = y[i] + h*k2[i]/2;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<input_number;i++)
    {
        k3[i] = func_Yin(t+h/2,tmp,i);
        tmp[i] = y[i] + h*k3[i];
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<input_number;i++)
    {
        k4[i] = func_Yin(t+h,tmp,i);
        y[i] = y[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        SetYin(y[i], i);
        //cout<<"y["<<i<<"]="<<GetYin(i)<<endl;
    }
    
    delete [] y;
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] tmp;
    y = 0;
    k1 = 0;
    k2 = 0;
    k3 = 0;
    k4 = 0;
    tmp = 0;
    
}

double RK::func_Yin(double t,double *y,int i)//inputのith成分
{
    double r=0.0,D_input=0.0,x=0.0;
    double growth=0.0;
    //double _mu =0.0;
    D_input = trspD;
    x = GetResult(i);//transporter とりあえずラストのタンパクから割り当てておく
    growth = GetGrowthRate();
    
    //r = c*D_input*x*GetY(i) - _mu*y[i];
    r = D_input*x*GetY(i)- growth*y[i];//transporterのみ
    //r = D_input*x*GetY(i) + S_D*(GetY(i)-y[i])-growth*y[i];//拡散あり
    return r;
}

void RK::calculate_growth()
{
    double r=0.0,f=0.0;
    double s=maximum_virulence;
    
    for(int i=0;i<target_number;i++)
    {
        r = GetResult(i);
        f = f + 1 -abs(GetTargetPattern(i)-r);
    }
    f = f - s*GetParasite_population_rate(GetParasite_attacked_number());//自分とマッチングするものを使う
    
    f=growth_const*f;//あとで定数で与える
    SetGrowth(f);
    //SetStack_fitness_Target_parasite_average(f);
}

void RK::calculate_growth_dist()
{
    double r=0.0,f=0.0,sum=0.0;
    double s=maximum_virulence;
    
    for(int i=0;i<target_number;i++)
    {
        r = GetResult(i);
        f = f + 1 -abs(GetTargetPattern(i)-r);//targetのパターンについての一致具合
    }
    for(int i=0;i<Parasite_species_num;i++)
    {
        sum =sum+ GetInfection_virulence(i);
    }
    
    f = f - sum;//自分とマッチングするものを使う
    
    f=growth_const*f;//あとで定数で与える
    SetGrowth(f);
    //SetStack_fitness_Target_parasite_average(f);
}

void RK::calculate_Fitness()
{
    double r=0.0,f=0.0;
    for(int i=0;i<target_number;i++)
    {
        r = GetResult(i);
        f = f + 2*(r-0.5)*(r-0.5);
    }
    
    SetFitness(f);
    
    //cout<<"fitness="<<GetFitness()<<endl;
    //cout<<"protein="<<GetResult(0)<<endl;
}

void RK::Calculate_Fitness_Target()//fitnessがtargetの時
{
    double r=0.0,f=0.0;
    
    for(int i=0;i<target_number;i++)
    {
        r = GetResult(i);
        f = f +1 -abs(GetTargetPattern(i)-r);
    }
    SetFitness(f);
}

void RK::Calculate_Fitness_Target_parasite()//fitnessがtargetの時
{
    double r=0.0,f=0.0,s=0.0;
    s = maximum_virulence;//parasiteの毒性の強さ
    for(int i=0;i<target_number;i++)
    {
        r = GetResult(i);
        f = f +1 -abs(GetTargetPattern(i)-r);
    }
    
    f = f - s*GetParasite_population();
    if(f<0)
    {
        f =0.0;
        SetFitness(f);
    }
    else{
        SetFitness(f);
    }
}


void RK::Calculate_Stack_fitness_Target_parasite_average()//fitnessがtargetの時parasiteあり averageをとる前に記録する
{
    double r=0.0,f=0.0,s=0.0;
    s = maximum_virulence;
    if(GetTime() >GetDiviTime()-average_time)
    {
        for(int i=0;i<target_number;i++)
        {
            r = GetResult(i);
            f = f +1 -abs(GetTargetPattern(i)-r);
        }
        //f  = 10.0;//targetをつける制約を入れない
        
        f = f - s*GetParasite_population_rate(GetParasite_attacked_number());//自分とマッチングするものを使う
        
        if(f<-s||f>10)
        {
            f =0.0;
            SetStack_fitness_Target_parasite_average(f);//値を保存しておく
        }
        else{
            SetStack_fitness_Target_parasite_average(f);
        }
    }
    
}

void RK::Calculate_Stack_fitness_Target_noparasite_average()//fitnessがtargetの時parasiteあり averageをとる前に記録する
{
    double r=0.0,f=0.0,s=0.0;
    s = maximum_virulence;
    if(GetTime() >GetDiviTime()-average_time)
    {
        for(int i=0;i<Parasite_target_genome;i++)
        {
            r = GetResult(i);
            f = f +1 -abs(GetTargetPattern(i)-r);
        }
        
        if(f<0||f>10)
        {
            f =0.0;
            SetStack_fitness_Target_noparasite_average(f);//値を保存しておく
        }
        else{
            SetStack_fitness_Target_noparasite_average(f);
        }
    }
}

//parasite用
int RK::GetParasite_attacked_number() const
{
    return parasite_attacked_number;
}
void RK::SetParasite_attacked_number(int num)
{
    parasite_attacked_number =num;
}

double RK::GetParasite_population() const
{
    return parasite_population;
}

void RK::SetParasite_population(double pp)
{
    parasite_population = pp;
}


double RK::GetParasite_population_rate(int ith)//全parasite情報を与える
{
    return Parasite_population_rate[ith];
}
void RK::SetParasite_population_rate(double pp,int ith)//全parasite情報を与える
{
    Parasite_population_rate[ith] = pp;
}

//その他関数//
void archive_J_3jigen(int *x,int species)
{
    FILE *pf = nullptr;
    
    char filename[50]={};
    
    species = species+1;
    sprintf(filename,"J_3jigen%d.txt",species);
    pf = fopen(filename,"a+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            fprintf(pf," %d ",x[l]);
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void archive_ProteinD(double x,int species)
{
    FILE *pf = nullptr;
    
    char filename[50]={};
    
    species = species+1;
    sprintf(filename,"ProteinD_%d.txt",species);
    pf = fopen(filename,"a+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %f ",x);
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void archive_threshold(double *x,int species)
{
    FILE *pf = nullptr;
    
    char filename[50]={};
    
    species = species+1;
    sprintf(filename,"p_threshold%d.txt",species);
    pf = fopen(filename,"a+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            fprintf(pf," %f ",x[l]);
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
}

void archive_J_3jigen_graph(int x,int species)
{
    FILE *pf = nullptr;
    
    char filename[50]={};
    
    species = species+1;
    sprintf(filename,"J_3jigen_graph%d.txt",species);
    pf = fopen(filename,"a+");
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        int memo=x;
        if(memo==-1)
        {
            fprintf(pf,"\"%d\"",memo+3);
        }
        else
        {
            fprintf(pf,"\"%d\"",memo);
        }
        fprintf(pf,",");
    }
    fclose(pf);
}

void remove_J_3jigen()
{
    char filename[50]={};
    char filenamer[50]={};
    for(int i=1;i<=Species_numb;i++)
    {
        sprintf(filename,"J_3jigen%d.txt",i);
        sprintf(filenamer,"J_3jigen_stationary%d.txt",i);
        remove(filename);
        remove(filenamer);
    }
}

void CELL::calculate_meanfield()
{
    int N=0,i,l;
    double n=0.0;
    N=rk[0].GetTotalcellnumber();
    n = (double)(rk[0].GetTotalcellnumber());
    
    vector<double> sum;
    sum.resize(ProteinVAL);
    for(l=0;l<ProteinVAL;l++)
    {
        sum[l] =0.0;
    }
    for(l=0;l<ProteinVAL;l++)
    {
        if(rk[0].GetProD(l)>0.0)
        {
            for(i=0;i<N;i++)
            {
                sum[l] = sum[l] + rk[i].GetResult(l);
            }
        }
        else{
            sum[l] = 0.0;
        }
    }
    for(l=0;l<ProteinVAL;l++)
    {
        sum[l] =sum[l]/n;
        for(i=0;i<N;i++)
        {
            rk[i].SetMeanfieldval(sum[l], l);
        }
    }
}
void CELL::calculate_NFDS_growth()
{
    double NFDS[MAXCELLNUMBER]={},SumValues[target_number]={};
    
    int cellnumber;
    
    cellnumber=rk[0].GetTotalcellnumber();
    
    //ここからgrowthの計算
    for(int i=0;i<cellnumber;i++)
    {
        NFDS[i] =0;
    }
    for(int l=0;l<target_number;l++)
    {
        SumValues[l] = 0.0;
    }
    for(int l=0;l<target_number;l++)
    {
        SumValues[l] = rk[0].GetResult(l)*rk[0].GetResult(l);//0th cellの発現量
    }
    for(int l=0;l<target_number;l++)
    {
        for(int i=1;i<cellnumber;i++)
        {
            SumValues[l] = SumValues[l]+rk[i].GetResult(l)*rk[i].GetResult(l);//ith cellの発現量の和
        }
    }
    for(int i=0;i<cellnumber;i++)
    {
        for(int l=0;l<target_number;l++)
        {
            NFDS[i] = NFDS[i] +(100*growth_const*rk[i].GetResult(l)*rk[i].GetResult(l))/SumValues[l];//ex/sum of ex
        }
        rk[i].SetGrowth(NFDS[i]);
    }
}


//ここまで

void CELL::calculate_Result_average()
{
    if(rk[0].GetTime() >rk[0].GetDiviTime()-average_time)
    {
        if(rk[0].GetTime() >=rk[0].GetDiviTime())
        {
            for(int l=0;l<rk[0].GetTotalcellnumber();l++)
            {
                for(int j=0;j<ProteinVAL;j++)
                {
                    rk[l].SetResult_average(rk[l].GetResult_average(j)/(average_time/rk[0].GetDelt()), j);
                    //足し合わせたResult_averageをここで規格化
                }
            }
            //cout <<"normalize_time="<<rk[0].GetTime()<<endl;
            //for(int j=0;j<ProteinVAL;j++){cout<<"rusult_average["<<j<<"]="<<rk[0].GetResult_average(j)<<endl;}
        }
        else{
            for(int l=0;l<rk[0].GetTotalcellnumber();l++)
            {
                for(int j=0;j<ProteinVAL;j++)
                {
                    rk[l].SetResult_average((rk[l].GetResult_average(j)+rk[l].GetResult(j)), j);
                }
            }
            if(on==1)
            {
                for(int i=0;i<check_num;i++)
                {
                    if(rk[0].GetTime() >=rk[0].GetDiviTime()-average_time +i*average_time/check_num&&
                       rk[0].GetTime()<rk[0].GetDiviTime()-average_time +rk[0].GetDelt() + i*average_time/check_num)
                    {
                        for(int l=0;l<rk[0].GetTotalcellnumber();l++)
                        {
                            for(int j=0;j<ProteinVAL;j++)
                            {
                                rk[l].SetResult_check(rk[l].GetResult(j), j, i);
                            }
                        }
                        //cout<<"average["<<rk[0].GetTime()<<"]="<<rk[0].GetResult_average(0)<<endl;
                        cout<<"Result_check[0]["<<i<<"]="<<rk[0].GetResult_check(0, i)<<endl;
                        //cout<<"time="<<rk[0].GetTime()<<endl;
                        
                        break;
                    }
                }
            }
        }
    }
}

void CELL::calculate_Fitness_Target_parasite_average()
{
    double fit=0.0;
    int totalcellnum=0;
    totalcellnum =rk[0].GetTotalcellnumber();
    if(rk[0].GetTime() >rk[0].GetDiviTime()-average_time)
    {
        if(rk[0].GetTime() >rk[0].GetDiviTime())
        {
            for(int l=0;l<totalcellnum;l++)
            {
                fit =rk[l].GetFitness_Target_parasite_average();
                if(fit>=0)
                {
                    rk[l].SetFitness_Target_parasite_average(rk[l].GetFitness_Target_parasite_average()/(average_time/rk[0].GetDelt()));
                }
                else{
                    rk[l].SetFitness_Target_parasite_average(0.0);
                }
                //足し合わせたFitness_averageをここで規格化
            }
        }
        else{
            for(int l=0;l<totalcellnum;l++)
            {
                rk[l].SetFitness_Target_parasite_average(rk[l].GetFitness_Target_parasite_average()+rk[l].GetStack_fitness_Target_parasite_average());
            }
        }
    }
}

void CELL::calculate_Fitness_Target_noparasite_average()
{
    if(rk[0].GetTime() >rk[0].GetDiviTime()-average_time)
    {
        if(rk[0].GetTime() >rk[0].GetDiviTime())
        {
            for(int l=0;l<rk[0].GetTotalcellnumber();l++)
            {
                rk[l].SetFitness_Target_noparasite_average(rk[l].GetFitness_Target_noparasite_average()/(average_time/rk[0].GetDelt()));
                //足し合わせたFitness_averageをここで規格化
            }
        }
        else{
            for(int l=0;l<rk[0].GetTotalcellnumber();l++)
            {
                rk[l].SetFitness_Target_noparasite_average(rk[l].GetFitness_Target_noparasite_average()+rk[l].GetStack_fitness_Target_noparasite_average());
            }
        }
    }
}

double CELL::GetSim(int ith,int lth) const
{
    return sim[ith][lth];
}

void CELL::SetSim(double e,int ith,int lth)
{
    sim[ith][lth] = e;
}

void CELL::calculate_sim(int ith,int jth)//細胞X(ith ,Y(jth
{
    int l,cellnum;;//i,j 細胞用 l gene
    double inner=0.0,Xd=0.0,Yd=0.0,E=0.0;
    cellnum = rk[0].GetTotalcellnumber();
    for(l=0;l<ProteinVAL;l++)
    {
        inner = inner + rk[ith].GetResult(l)*rk[jth].GetResult(l);//内積
    }
    for(l=0;l<ProteinVAL;l++)
    {
        Xd = Xd + rk[ith].GetResult(l)*rk[ith].GetResult(l);//ノルム
    }
    Xd = sqrt(Xd);
    for(l=0;l<ProteinVAL;l++)
    {
        Yd = Yd + rk[jth].GetResult(l)*rk[jth].GetResult(l);//ノルム
    }
    Yd = sqrt(Yd);
    E = inner/(Xd*Yd);//cosine距離
    SetSim(E,ith,jth);
}

void CELL::calculate_simJij(int ith,int jth)//細胞X(ith ,Y(jth
{
    int i,l,cellnum;;//i,j 細胞用 l gene
    double inner=0.0,Xd=0.0,Yd=0.0,E=0.0;
    cellnum = rk[0].GetTotalcellnumber();
    for(i=0;i<ProteinVAL;i++)
    {
        inner=0.0;Xd=0.0;Yd=0.0;
        for(l=0;l<ProteinVAL;l++)
        {
            inner = inner + rk[ith].GetJ(i, l)*rk[jth].GetJ(i, l);//内積
        }
        for(l=0;l<ProteinVAL;l++)
        {
            Xd = Xd + rk[ith].GetJ(i, l)*rk[ith].GetJ(i, l);//ノルム
        }
        Xd = sqrt(Xd);
        for(l=0;l<ProteinVAL;l++)
        {
            Yd = Yd + rk[jth].GetJ(i, l)*rk[jth].GetJ(i, l);//ノルム
        }
        Yd = sqrt(Yd);
        E = E+inner/(Xd*Yd);//cosine距離
    }
    E= (double)(E/ProteinVAL);
    SetGenotype_sim(E,ith,jth);
}

void CELL::calculate_EuclidDistanceJij(int ith,int jth)//細胞X(ith ,Y(jth
{
    int i,l;//i,j 細胞用 l gene
    double inner=0.0;
    for(i=0;i<ProteinVAL;i++)
    {
        inner=0.0;
        for(l=0;l<ProteinVAL;l++)
        {
            inner = inner + (rk[ith].GetJ(i, l)-rk[jth].GetJ(i, l))*(rk[ith].GetJ(i, l)-rk[jth].GetJ(i, l));//内積
        }
    }
    inner = sqrt(inner);
    SetGenotype_sim(inner,ith,jth);
}

double CELL::calculate_EuclidDistanceJij_initial(int ith)//細胞X(ith ,Y(jth
{
    int i,l;//i,j 細胞用 l gene
    double inner=0.0;
    for(i=0;i<ProteinVAL;i++)
    {
        inner=0.0;
        for(l=0;l<ProteinVAL;l++)
        {
            inner = inner + (GetJij_initial(i, l)-rk[ith].GetJ(i, l))*(GetJij_initial(i, l)-rk[ith].GetJ(i, l));
        }
    }
    inner = sqrt(inner);
    return inner;
}

void CELL::calculate_sim_sta(int ith,int jth)//細胞X(ith ,Y(jth
{
    int l,cellnum;;//i,j 細胞用 l gene
    double inner=0.0,Xd=0.0,Yd=0.0,E=0.0;
    cellnum = rk[0].GetTotalcellnumber();
    for(l=0;l<ProteinVAL;l++)
    {
        inner = inner + rk[ith].GetResult_average(l)*rk[jth].GetResult_average(l);//内積
    }
    for(l=0;l<ProteinVAL;l++)
    {
        Xd = Xd + rk[ith].GetResult_average(l)*rk[ith].GetResult_average(l);//ノルム
    }
    Xd = sqrt(Xd);
    for(l=0;l<ProteinVAL;l++)
    {
        Yd = Yd + rk[jth].GetResult_average(l)*rk[jth].GetResult_average(l);//ノルム
    }
    Yd = sqrt(Yd);
    E = inner/(Xd*Yd);//cosine距離
    SetSim(E,ith,jth);
}

double CELL::GetGenotype_sim(int ith ,int lth) const
{
    return genotype_sim[ith][lth];
}
void CELL::SetGenotype_sim(double e,int ith,int lth)
{
    if(e<0)
    {
        cout << "genotype_sim = "<<e<<endl;
        exit(1);
    }
    else
    {
        genotype_sim[ith][lth] = e;
    }
}
double CELL::GetJij_variance(int ith,int lth) const
{
    return Jij_variance[ith][lth];
}
void CELL::SetJij_variance(double v,int ith,int lth)
{
    if(v<-10)
    {
        cout << "Jij_variance = "<<v<<endl;
        exit(1);
    }
    else
    {
        Jij_variance[ith][lth] = v;
    }
}

double CELL::GetJij_initial(int ith,int lth) const
{
    return Jij_initial[ith][lth];
}
void CELL::SetJij_initila(double v,int ith,int lth)
{
    if(v<-10)
    {
        cout << "Jij_variance = "<<v<<endl;
        exit(1);
    }
    else
    {
        Jij_initial[ith][lth] = v;
    }
}


void CELL::archive_Jij_Variance(int time)
{
    FILE *pf = nullptr,*pg = nullptr,*ph = nullptr;
    
    char filename[50]={},filenameg[50]={},filenameh[50]={};
    double sum=0.0,sumall=0.0;
    sprintf(filename,"Jij_var%d.txt",time);
    sprintf(filenameg,"Jij_var.txt");
    sprintf(filenameh,"Jij_var_all.txt");
    
    pf = fopen(filename,"w+");
    pg = fopen(filenameg,"a+");
    ph = fopen(filenameh,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        Calculate_Jij_Variance();//ここで分散を計算
        //fprintf(pf,"%d",time);
        //fprintf(pf,"\n");
        fprintf(ph,"%d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            sum = 0.0;
            for(int l=0;l<ProteinVAL;l++)
            {
                sum = sum+GetJij_variance(i,l);
                sumall = sumall+GetJij_variance(i,l);
            }
            sum = (double)(sum/ProteinVAL);
            fprintf(pf,"#%d %.8f",i+1,sum);
            fprintf(pf,"\n");
            fprintf(ph,"%.8f ",sum);
        }
        fprintf(ph,"\n");
        sumall= (double)(sumall/(ProteinVAL*ProteinVAL));
        fprintf(pf,"#all %.8f",sumall);
        fprintf(pf,"\n");
        fprintf(pg,"%d %.8f",time,sumall);
        fprintf(pg,"\n");
    }
    fclose(pf);
    fclose(pg);
    fclose(ph);
}

void CELL::archive_sim(int time)
{
    FILE *pf = nullptr;
    FILE *pg = nullptr;
    FILE *ph = nullptr;
    
    int totalcellnum=0,ith,jth,iDist1=0,iDist2=0,p_val=ProteinVAL,memo;
    double E=0.0;
    double G=0;
    char filename[50]={},filenameg[50]={},filenameh[50]={};
    bool flags=false;
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    sprintf(filename,"sim%d.txt",time);
    pf = fopen(filename,"w+");
    sprintf(filenameg,"genotype_sim%d.txt",time);
    pg = fopen(filenameg,"w+");
    sprintf(filenameh,"cell_sim%d.txt",time);//sortした記録
    ph = fopen(filenameh,"w+");
    
    if (flags){
        for(ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
        {
            data_array[ith].str.clear();
            data_array[ith].hammD =0;
            data_array[ith].Relative_hammD =0;
            data_array[ith].num =0;
        }
        for(int i=0;i<totalcellnum;i++)
        {
            for(int j=0;j<totalcellnum;j++)
            {
                if(i==j)
                {
                    SetGenotype_sim(0,i,j);
                }
                else
                {
                    //Jijが連続変化する時は以下
                    //calculate_EuclidDistanceJij(i,j);
                }
            }
            data_array[i].hammD =GetGenotype_sim(0,i);
            data_array[i].Relative_hammD =calculate_EuclidDistanceJij_initial(i);
        }
        for(int j=0;j<totalcellnum;j++)
        {
            if(j==0)
            {
                SetSim(1.0,0,j);
            }
            else
            {
                calculate_sim(0,j);//i細胞とj細胞の類似度計算 定常状態
            }
        }
        
        for(int i=0;i<totalcellnum;i++)
        {
            rk[i].SetHammingD(data_array[i].hammD);//元のgenetypeと比べた時の絶対hamming距離iDist1を更新
            rk[i].SetRelativeHammingD(data_array[i].Relative_hammD);//相対Hamming距離iDist(0,i)を更新
        }
        
        SetData_array();//ここで絶対hamming距離と相対hamming距離とtreeを更新
        
        for(ith=0;ith<totalcellnum;ith++)//Xith
        {
            data_array[ith].num = ith;//sort前検索用の番号をセット
        }
        
        sort(data_array.begin(), data_array.begin()+totalcellnum);
        
    }
    else{
        //Jijが離散の時
        for(ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
        {
            data_array[ith].str.clear();
            data_array[ith].hammD =0;
            data_array[ith].Relative_hammD =0;
            data_array[ith].num =0;
        }
        
        for(int i=0;i<totalcellnum;i++)
        {
            for(int j=0;j<totalcellnum;j++)
            {
                if(i==j)
                {
                    SetSim(1.0, i,j);
                }
                else
                {
                    calculate_sim(i,j);//i細胞とj細胞の類似度計算
                }
            }
        }
        
        for(int i=0;i<totalcellnum;i++)
        {
            rk[i].Sstr_J.clear();
            rk[i].Sstr_J.reserve(p_val*p_val);
            for(int l=0;l<p_val;l++)
            {
                for(int k=0;k<p_val;k++)
                {
                    memo = rk[i].GetJ(l,k);
                    if(memo == -1)
                    {
                        rk[i].Sstr_J = rk[i].Sstr_J + 'a' ;//今のネットワーク構造 GRN -1を文字aで記録
                    }
                    else if(memo == 0)
                    {
                        rk[i].Sstr_J = rk[i].Sstr_J + 'b' ;//今のネットワーク構造 GRN 0を文字bで記録
                    }
                    else
                    {
                        rk[i].Sstr_J = rk[i].Sstr_J + 'c' ;//今のネットワーク構造 GRN 1を文字cで記録
                    }
                }
            }
            iDist1 = iHammingDist(rk[i].Sstr_J,rk[i].ori_Sstr_J);
            rk[i].SetHammingD(iDist1);//元のgenetypeと比べた時の絶対hamming距離iDist1を更新
        }
        for(int i=0;i<totalcellnum;i++)
        {
            for(int j=0;j<totalcellnum;j++)
            {
                if(i==j)
                {
                    SetGenotype_sim(0,i,j);
                }
                else
                {
                    iDist2 = iHammingDist(rk[i].Sstr_J,rk[j].Sstr_J);
                    SetGenotype_sim(iDist2,i,j);//類似度を計算 相対hamming距離
                }
            }
            rk[i].SetRelativeHammingD(GetGenotype_sim(0, i));//相対Hamming距離iDist(0,i)を更新
        }
        SetData_array();//ここで絶対hamming距離と相対hamming距離とtreeを更新
        for(ith=0;ith<totalcellnum;ith++)//Xith
        {
            data_array[ith].num = ith;//sort前検索用の番号をセット
        }
        
        sort(data_array.begin(), data_array.begin()+totalcellnum);
    }
 
    
    if(pg==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(ith=0;ith<totalcellnum;ith++)//Xith
        {
            //cout<<"data_array["<<ith<<"]="<<data_array[ith].num<<endl;
            //cout<<"rk["<<data_array[ith].num<<"]_type="<<rk[data_array[ith].num].GetCelltype()<<endl;
            //cout<<data_array[ith].hammD<<endl;
            //cout<<data_array[ith].str<<endl;
            for(jth=0;jth<totalcellnum;jth++)//Yjth
            {
                //E = GetSim(ith, jth);
                //fprintf(pf,"%d %d %f",ith+1,jth+1,E);
                E = GetSim(data_array[ith].num,data_array[jth].num);
                G = GetGenotype_sim(data_array[ith].num,data_array[jth].num);
                //G = GetGenotype_sim(ith,jth);
                fprintf(pf,"%d %d %f #%d #%d #%d",ith+1,jth+1,E,data_array[ith].num,data_array[jth].num,rk[data_array[ith].num].GetParasite_attacked_number());
                fprintf(pf,"\n");
                fprintf(pg,"%d %d %f %d %d %d",ith+1,jth+1,G,data_array[ith].num,data_array[jth].num,rk[data_array[jth].num].GetParasite_attacked_number());
                fprintf(pg,"\n");
            }
            fprintf(pf,"\n");
            fprintf(pg,"\n");
            fprintf(ph," %d ",time);
            fprintf(ph," %d " ,ith+1);
            fprintf(ph," %d " ,data_array[ith].num);
            fprintf(ph," %d " ,rk[data_array[ith].num].GetParasite_attacked_number());
            for(int l=0;l<ProteinVAL;l++)
            {
                fprintf(ph," %f ",rk[data_array[ith].num].GetResult(l));
            }
            fprintf(ph," %f ",rk[data_array[ith].num].GetHammingD());
            fprintf(ph," %f ",rk[data_array[ith].num].GetRelativeHammingD());
            fprintf(ph," %f ",GetSim(0,data_array[ith].num));//0番目のcellとの類似度を記録する
            fprintf(ph," %f ",rk[data_array[ith].num].GetGrowthRate());
            fprintf(ph," %f ",rk[data_array[ith].num].GetVolume());
            fprintf(ph,"\n");
        }
    }
    fclose(pf);
    fclose(pg);
    fclose(ph);
}

void CELL::archive_sim_sta(int time)
{
    //FILE *pf = nullptr;
    FILE *pg = nullptr;
    FILE *ph = nullptr;
    
    int totalcellnum=0,ith,jth,iDist1=0,iDist2=0,p_val=ProteinVAL,memo;
    //double E=0.0;
    int G=0;
    //char filename[50]={};
    char filenameg[50]={},filenameh[50]={};
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    //sprintf(filename,"sim_sta%d.txt",time);
    //pf = fopen(filename,"w+");
    sprintf(filenameg,"genotype_sta_sim%d.txt",time);
    pg = fopen(filenameg,"w+");
    sprintf(filenameh,"cell_sta_sim%d.txt",time);//sortした記録
    ph = fopen(filenameh,"w+");
    
    for(ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_sta[ith].hammD =0;
        data_array_sta[ith].Relative_hammD =0;
        data_array_sta[ith].num =0;
    }
    
    for(int i=0;i<totalcellnum;i++)
    {
        for(int j=0;j<totalcellnum;j++)
        {
            if(i==j)
            {
                //SetSim(1.0, i,j);
            }
            else
            {
                //calculate_sim_sta(i,j);//i細胞とj細胞の類似度計算 定常状態
            }
        }
    }
    
    for(int i=0;i<totalcellnum;i++)
    {
        rk[i].Sstr_J.clear();
        rk[i].Sstr_J.reserve(p_val*p_val);
        for(int l=0;l<p_val;l++)
        {
            for(int k=0;k<p_val;k++)
            {
                memo = rk[i].GetStationaryJ(l,k);
                if(memo == -1)
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'a' ;//今のネットワーク構造 GRN -1を文字aで記録
                }
                else if(memo == 0)
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'b' ;//今のネットワーク構造 GRN 0を文字bで記録
                }
                else
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'c' ;//今のネットワーク構造 GRN 1を文字cで記録
                }
            }
        }
        iDist1 = iHammingDist(rk[i].Sstr_J,rk[i].ori_Sstr_J);
        rk[i].SetHammingD(iDist1);//元のgenetypeと比べた時の絶対hamming距離iDist1を更新
    }
    for(int i=0;i<totalcellnum;i++)
    {
        for(int j=0;j<totalcellnum;j++)
        {
            if(i==j)
            {
                SetGenotype_sim(0,i,j);
            }
            else
            {
                iDist2 = iHammingDist(rk[i].Sstr_J,rk[j].Sstr_J);
                SetGenotype_sim(iDist2,i,j);//類似度を計算 相対hamming距離
            }
        }
        rk[i].SetRelativeHammingD(GetGenotype_sim(0, i));//相対Hamming距離iDist(0,i)を更新
    }
    
    SetData_array_sta();//ここで絶対hamming距離と相対hamming距離とtreeを更新
    
    for(ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_sta[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_sta.begin(), data_array_sta.begin()+totalcellnum);
    if(ph==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(ith=0;ith<totalcellnum;ith++)//Xith
        {
            //cout<<"data_array["<<ith<<"]="<<data_array[ith].num<<endl;
            //cout<<"rk["<<data_array[ith].num<<"]_type="<<rk[data_array[ith].num].GetCelltype()<<endl;
            //cout<<data_array[ith].hammD<<endl;
            //cout<<data_array[ith].str<<endl;
            for(jth=0;jth<totalcellnum;jth++)//Yjth
            {
                //E = GetSim(ith, jth);
                //fprintf(pf,"%d %d %f",ith+1,jth+1,E);
                //E = GetSim(data_array_sta[ith].num,data_array_sta[jth].num);//sort済みの類似度
                G = GetGenotype_sim(data_array_sta[ith].num,data_array_sta[jth].num);
                //G = GetGenotype_sim(ith,jth);
                //fprintf(pf,"%d %d %f #%d #%d",ith+1,jth+1,E,data_array_sta[ith].num,data_array_sta[jth].num);
                //fprintf(pf,"\n");
                fprintf(pg,"%d %d %d #%d #%d %d",ith+1,jth+1,G,data_array_sta[ith].num,data_array_sta[jth].num,rk[data_array_sta[jth].num].GetParasite_attacked_number());
                fprintf(pg,"\n");
            }
            //fprintf(pf,"\n");
            fprintf(pg,"\n");
            fprintf(ph," %d ",time);
            fprintf(ph," %d " ,ith+1);
            fprintf(ph," #%d " ,data_array_sta[ith].num);
            fprintf(ph," #%f " ,rk[ith].GetStationaryTime());
            //fprintf(ph," %s ",rk[data_array_sta[ith].num].tree_history.c_str());
            fprintf(ph," #%d ",rk[data_array_sta[ith].num].GetParasite_attacked_number());
            for(int l=0;l<ProteinVAL;l++)
            {
                fprintf(ph," %f ",rk[data_array_sta[ith].num].GetStationaryResult(l));
                //fprintf(ph," %f ",rk[data_array_sta[ith].num].GetResult_average(l));
            }
            //fprintf(ph," %f ",rk[data_array_sta[ith].num].GetFitness_Target_parasite_average());
            fprintf(ph," %f ",rk[data_array_sta[ith].num].GetFitness());
            fprintf(ph," %f ",rk[data_array_sta[ith].num].GetGrowthRate());
            fprintf(ph," %f ",rk[data_array_sta[ith].num].GetVolume());
            fprintf(ph,"\n");
        }
    }
    //fclose(pf);
    fclose(pg);
    fclose(ph);
}

void CELL::genosim_label_set()
{
    int totalcellnum=0,ith,iDist1=0,iDist2=0,p_val=ProteinVAL,memo;
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    for(ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_sta[ith].hammD =0;
        data_array_sta[ith].Relative_hammD =0;
        data_array_sta[ith].num =0;
    }
    
    for(int i=0;i<totalcellnum;i++)
    {
        rk[i].Sstr_J.clear();
        rk[i].Sstr_J.reserve(p_val*p_val);
        for(int l=0;l<p_val;l++)
        {
            for(int k=0;k<p_val;k++)
            {
                memo = rk[i].GetJ(l,k);
                if(memo == -1)
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'a' ;//今のネットワーク構造 GRN -1を文字aで記録
                }
                else if(memo == 0)
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'b' ;//今のネットワーク構造 GRN 0を文字bで記録
                }
                else
                {
                    rk[i].Sstr_J = rk[i].Sstr_J + 'c' ;//今のネットワーク構造 GRN 1を文字cで記録
                }
            }
        }
        iDist1 = iHammingDist(rk[i].Sstr_J,rk[i].ori_Sstr_J);
        rk[i].SetHammingD(iDist1);//元のgenetypeと比べた時の絶対hamming距離iDist1を更新
    }
    for(int i=0;i<totalcellnum;i++)
    {
        for(int j=0;j<totalcellnum;j++)
        {
            if(i==j)
            {
                SetGenotype_sim(0,i,j);
            }
            else
            {
                iDist2 = iHammingDist(rk[i].Sstr_J,rk[j].Sstr_J);
                SetGenotype_sim(iDist2,i,j);//類似度を計算 相対hamming距離
            }
        }
        rk[i].SetRelativeHammingD(GetGenotype_sim(0, i));//相対Hamming距離iDist(0,i)を更新
    }
    
    SetData_array_sta();//ここで絶対hamming距離と相対hamming距離とtreeを更新
    
    for(ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_sta[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_sta.begin(), data_array_sta.begin()+totalcellnum);
    
    for(int i=0;i<totalcellnum;i++)
    {
        rk[i].SetSorting_number(data_array_sta[i].num);//sortした番号を記録
        cout<<"data_array_sta["<<i<<"].num="<<data_array_sta[i].num<<endl;
    }
}

void CELL::SetData_array()
{
    int totalcellnum =0;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(int i=0;i<totalcellnum;i++)
    {
        data_array[i].hammD = rk[i].GetHammingD();
        data_array[i].Relative_hammD = rk[i].GetRelativeHammingD();
        data_array[i].str = rk[i].tree_history;
    }
}

void CELL::SetData_array_sta()
{
    int totalcellnum =0;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(int i=0;i<totalcellnum;i++)
    {
        data_array_sta[i].hammD = rk[i].GetHammingD();
        data_array_sta[i].Relative_hammD = rk[i].GetRelativeHammingD();
    }
}

void CELL::SetData_array_fitness()
{
    int totalcellnum =0;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(int i=0;i<totalcellnum;i++)
    {
        data_array_fitness[i].hammD = rk[i].GetHammingD();
        //data_array_fitness[i].fit = rk[i].GetFitness();//平均取らない
        data_array_fitness[i].fit = rk[i].GetFitness_Target_parasite_average();//時間平均ver
    }
}

void CELL::SetVariance(double var,int i)
{
    if(var<0)
    {
        cout <<"var="<<var;
        exit(-1);
    }
    else
    {
        variance[i] = var;
    }
}

void CELL::SetVariance_thereshold(double var,int i)
{
    if(var<0)
    {
        cout <<"var="<<var;
        exit(-1);
    }
    else
    {
        variance_threshold[i] = var;
    }
}

double CELL::GetVariance(int i)
{
    return variance[i];
}

double CELL::GetVariance_thereshold(int i)
{
    return variance_threshold[i];
}
void CELL::SetFitness_variance_by_type(double var,int i)
{
    if(var<0)
    {
        cout <<"Fitness_variance_by_type="<<var;
        exit(-1);
    }
    else
    {
        fitness_variance_by_type[i] = var;
    }
}

double CELL::GetFitness_variance_by_type(int i)
{
    return fitness_variance_by_type[i];
}

void CELL::SetMean_fitness_by_type(double mean,int i)
{
    if(mean<0)
    {
        mean_fitness_by_type[i] = 0.0;
    }
    else
    {
        mean_fitness_by_type[i] = mean;
    }
}

double CELL::GetMean_fitness_by_type(int i)
{
    return mean_fitness_by_type[i];
}

double CELL::GetMean_fitness_nontarget_by_type(int i)
{
    return mean_fitness_nontarget_by_type[i];
}

void CELL::SetMean_fitness_nontarget_by_type(double mean,int i)
{
    if(mean<0)
    {
        mean_fitness_nontarget_by_type[i] = 0.0;
    }
    else
    {
        mean_fitness_nontarget_by_type[i] = mean;
    }
}

void CELL::SetGenome_Variance(double var)
{
    if(var<0)
    {
        cout <<"Genome_Variance="<<var;
        exit(-1);
    }
    else
    {
        genome_variance = var;
    }
}
double CELL::GetGenome_Variance() const
{
    return genome_variance;
}
void CELL::SetFitness_variance(double var)
{
    if(var<0)
    {
        cout <<"Fitness_variance="<<var;
        exit(-1);
    }
    else
    {
        fitness_variance = var;
    }
}
double CELL::GetFitness_variance() const
{
    return fitness_variance;
}
void CELL::SetMean_fitness(double mean)
{
    if(mean<0)
    {
        mean_fitness=0;
    }
    else
    {
        mean_fitness = mean;
    }
}
double CELL::GetMean_fitness() const
{
    return mean_fitness;
}
void CELL::Calculate_Variance()
{
    int i,j,totalcellnum;
    double sum,var;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(i=0;i<ProteinVAL;i++)
    {
        sum =0.0;
        var=0.0;
        
        for(j=0;j<totalcellnum;j++)
        {
            sum = sum +rk[j].GetResult(i);
        }
        sum = (double)(sum/totalcellnum);//平均
        //cout<<"sum="<<sum<<endl;
        for(j=0;j<totalcellnum;j++)
        {
            var =  var + (sum -rk[j].GetResult(i))*(sum -rk[j].GetResult(i));
        }
        var = (double)(var/(totalcellnum-1));
        //cout<<"var="<<var<<endl;
        SetVariance(var, i);
    }
}

void CELL::Calculate_Variance_sta()
{
    int i,j,totalcellnum;
    double sum,var;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(i=0;i<ProteinVAL;i++)
    {
        sum =0.0;
        var=0.0;
        for(j=0;j<totalcellnum;j++)
        {
            sum = sum +rk[j].GetResult_average(i);
        }
        sum = (double)(sum/totalcellnum);//平均
        //cout<<"sum="<<sum<<endl;
        for(j=0;j<totalcellnum;j++)
        {
            var =  var + (sum -rk[j].GetResult_average(i))*(sum -rk[j].GetResult_average(i));
        }
        var = (double)(var/(totalcellnum-1));
        //cout<<"var="<<var<<endl;
        SetVariance(var, i);
    }
}

void CELL::Calculate_Variance_threshold()
{
    int i,j,totalcellnum;
    double sum,var;
    totalcellnum = rk[0].GetTotalcellnumber();
    for(i=0;i<ProteinVAL;i++)
    {
        sum =0.0;
        var=0.0;
        
        for(j=0;j<totalcellnum;j++)
        {
            sum = sum +rk[j].GetThreshold_g(i);
        }
        sum = (double)(sum/totalcellnum);//平均
        //cout<<"sum="<<sum<<endl;
        for(j=0;j<totalcellnum;j++)
        {
            var =  var + (sum -rk[j].GetThreshold_g(i))*(sum -rk[j].GetThreshold_g(i));
        }
        var = (double)(var/(totalcellnum-1));
        //cout<<"var="<<var<<endl;
        SetVariance_thereshold(var, i);
    }
}

void CELL::Calculate_fitness_Variance()
{
    int j,totalcellnum;
    double sum,var,fit;
    
    sum =0.0;
    var=0.0;
    totalcellnum=0;
    
    for(j=0;j<rk[0].GetTotalcellnumber();j++)
    {
        fit =rk[j].GetFitness_Target_parasite_average();
        if(fit>0&&fit<=10)
        {
            sum = sum +fit;
            totalcellnum =totalcellnum+1;
        }
    }
    sum = (double)(sum/totalcellnum);//平均
    //cout<<"sum="<<sum<<endl;
    for(j=0;j<rk[0].GetTotalcellnumber();j++)
    {
        fit =rk[j].GetFitness_Target_parasite_average();
        if(fit>0&&fit<=10)
        {
            var =  var + (sum -fit)*(sum -fit);
        }
    }
    var = (double)(var/(totalcellnum-1));
    SetFitness_variance(var);
    SetMean_fitness(sum);
    
}

void CELL::Calculate_fitness_Variance_by_type()
{
    int i,j,totalcellnum,paranum=0;
    double sum[Parasite_species_num]={},var[Parasite_species_num]={},sum_nontarget[Parasite_species_num]={},fit=0.0,fit_non_target=0.0;
    int totalcellnum_by_type[Parasite_species_num]={};
    
    totalcellnum = rk[0].GetTotalcellnumber();
    
    for(i=0;i<Parasite_species_num;i++)
    {
        totalcellnum_by_type[i] = 0;
        sum[i] = 0.0;
        var[i] = 0.0;
        sum_nontarget[i] = 0.0;
    }
    
    for(j=0;j<totalcellnum;j++)
    {
        paranum = rk[j].GetParasite_attacked_number();
        fit = rk[j].GetFitness_Target_parasite_average();
        if(fit<=10&&fit>0)
        {
            sum[paranum] = sum[paranum] +fit;
            totalcellnum_by_type[paranum]=totalcellnum_by_type[paranum]+1;
        }
        fit_non_target = rk[j].GetFitness_Target_noparasite_average();
        if(fit_non_target<=10&&fit_non_target>0)
        {
            sum_nontarget[paranum] = sum_nontarget[paranum] +fit_non_target;
        }
        //cout<<"sum["<<paranum<<"]="<<sum[paranum]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //cout<<"sum["<<i<<"]="<<sum[i]<<endl;
        //cout<<"sum_nontarget["<<i<<"]="<<sum_nontarget[i]<<endl;
        if(totalcellnum_by_type[i]>0)
        {
            sum[i] = (double)(sum[i]/totalcellnum_by_type[i]);//平均化
            sum_nontarget[i] = (double)(sum_nontarget[i]/totalcellnum_by_type[i]);//平均化
        }
        else{
            sum[i] =0.0;
            sum_nontarget[i] =0.0;
        }
        //cout<<"sum["<<i<<"]="<<sum[i]<<endl;
        //cout<<"sum_nontarget["<<i<<"]="<<sum_nontarget[i]<<endl;
    }
    
    for(j=0;j<totalcellnum;j++)
    {
        paranum = rk[j].GetParasite_attacked_number();
        fit = rk[j].GetFitness_Target_parasite_average();
        if(fit<=10&&fit>0)
        {
            var[paranum] =  var[paranum] + (sum[paranum] -fit)*(sum[paranum] -fit);
        }
        //cout<<"var["<<paranum<<"]="<<var[paranum]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        if(totalcellnum_by_type[i]>1)
        {
            var[i] = (double)(var[i]/(totalcellnum_by_type[i]-1));//平均化
        }
        else
        {
            var[i] =0.0;
        }
        SetFitness_variance_by_type(var[i],i);
        SetMean_fitness_by_type(sum[i],i);
        SetMean_fitness_nontarget_by_type(sum_nontarget[i],i);
    }
}

void CELL::Calculate_genome_Variance()
{
    int j,totalcellnum;
    double sum,var;
    totalcellnum = rk[0].GetTotalcellnumber();
    
    sum =0.0;
    var=0.0;
    for(j=0;j<totalcellnum;j++)
    {
        sum = sum +rk[j].GetHammingD();
    }
    sum = (double)(sum/totalcellnum);//平均
    //cout<<"sum="<<sum<<endl;
    for(j=0;j<totalcellnum;j++)
    {
        var =  var + (sum -rk[j].GetHammingD())*(sum -rk[j].GetHammingD());
    }
    var = (double)(var/(totalcellnum-1)*sum);
    //cout<<"var="<<var<<endl;
    SetGenome_Variance(var);
    
}

void CELL::Calculate_Jij_Variance()
{
    int i,j,l,totalcellnum;
    double sum,var;
    totalcellnum = rk[0].GetTotalcellnumber();
    
    for(i=0;i<ProteinVAL;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            sum =0.0;
            var=0.0;
            for(j=0;j<totalcellnum;j++)
            {
                sum = sum +rk[j].GetJ(i,l);
            }
            sum = (double)(sum/totalcellnum);//平均
            //cout<<"sum="<<sum<<endl;
            for(j=0;j<totalcellnum;j++)
            {
                var =  var + (sum -rk[j].GetJ(i,l))*(sum -rk[j].GetJ(i,l));
            }
            var = (double)(var/(totalcellnum-1));
            //cout<<"var="<<var<<endl;
            SetJij_variance(var,i,l);
        }
    }
}

void CELL::archive_var(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    
    sprintf(filename,"var0001.txt");
    pf = fopen(filename,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        Calculate_Variance();//ここで分散を計算
        fprintf(pf," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pf,"%.8f ",GetVariance(i));
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_var_sta(int time)
{
    FILE *pf = nullptr;
    //FILE *pg = nullptr;
    FILE *pr = nullptr;
    char filename[50]={};
    //char filenameg[50]={};
    char filenamer[50]={};
    FILE *pb = nullptr;
    char filenameb[50]={};
    
    double sum=0.0;
    sprintf(filename,"var_sta.txt");
    pf = fopen(filename,"a+");
    //sprintf(filenameg,"genome_var0001.txt");
    //pg = fopen(filenameg,"a+");
    
    sprintf(filenamer,"fitness_var.txt");//全部の分散
    pr = fopen(filenamer,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        Calculate_Variance_sta();//ここで分散を計算　平均値を使う
        //Calculate_genome_Variance();
        Calculate_fitness_Variance();
        Calculate_fitness_Variance_by_type();//系列ごとの計算
        fprintf(pf," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pf,"%.8f ",GetVariance(i));
            sum = sum + GetVariance(i);
        }
        sum =(double)(sum/ProteinVAL);
        fprintf(pf,"%.8f ",sum);//分散の平均値
        fprintf(pf,"\n");
        //fprintf(pg," %d ",time);
        //fprintf(pg,"%.8f ",GetGenome_Variance());
        //fprintf(pg,"\n");
        fprintf(pr," %d ",time);
        fprintf(pr,"%f ",GetMean_fitness());//平均fitness
        fprintf(pr,"%.8f ",GetFitness_variance());
        fprintf(pr,"\n");
        for(int i=0;i<Parasite_species_num;i++)
        {
            sprintf(filenameb,"fitness_var%04d.txt",i+1);
            pb = fopen(filenameb,"a+");
            fprintf(pb," %d ",time);
            fprintf(pb,"%f ",GetMean_fitness_by_type(i));//系列ごとの平均fitness
            fprintf(pb,"%.8f ",GetFitness_variance_by_type(i));
            fprintf(pb,"%f ",GetMean_fitness_nontarget_by_type(i));//系列ごとの平均fitness
            fprintf(pb,"\n");
            fclose(pb);
        }
    }
    fclose(pf);
    //fclose(pg);
    fclose(pr);
}

void CELL::archive_var_threshold(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    
    sprintf(filename,"var_threshold.txt");
    pf = fopen(filename,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        Calculate_Variance_threshold();//ここで分散を計算
        fprintf(pf," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pf,"%.8f ",GetVariance_thereshold(i));
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_var_genome(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    
    sprintf(filename,"var_genome.txt");
    pf = fopen(filename,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        Calculate_genome_Variance();//ここで分散を計算
        fprintf(pf," %d ",time);
        fprintf(pf,"%.8f ",GetGenome_Variance());
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_Sigma_Str_all(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    int MaxCellNum=0;
    
    sprintf(filename,"Sigma_Str.txt");
    pf = fopen(filename,"a+");
    MaxCellNum = rk[0].GetTotalcellnumber();
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<MaxCellNum;i++)
        {
            fprintf(pf,"%d %d %.8f \n",time,i,rk[i].GetSigma_Str());
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_Sigma_Str(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    int MaxCellNum=0;
    
    sprintf(filename,"Sigma_Str%d.txt",time);
    
    pf = fopen(filename,"a+");
    MaxCellNum = rk[0].GetTotalcellnumber();
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<MaxCellNum;i++)
        {
            fprintf(pf,"%d %d %.8f \n",time,i,rk[i].GetSigma_Str());
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_Parasite_Virulence(int time)
{
    
    FILE *pf = nullptr;
    char filename[50]={};
    int MaxCellNum=0;
    
    sprintf(filename,"Parasite_Virulence.txt");
    pf = fopen(filename,"a+");
    MaxCellNum = rk[0].GetTotalcellnumber();
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<Parasite_species_num;i++)
        {
            fprintf(pf,"%d %d %f \n",time,i,GetParasite_virulence(i));
        }
        fprintf(pf,"\n");
    }
    fclose(pf);
}

void CELL::archive_threshold(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    
    sprintf(filename,"Threshold_Str%d.txt",time);
    pf = fopen(filename,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int j=0;j<MAXCELLNUMBER;j++)
        {
            fprintf(pf," %d ",j);
            for(int i=0;i<ProteinVAL;i++)
            {
                fprintf(pf,"%.8f ",rk[j].GetThreshold_g(i));
            }
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
}

void CELL::Calculate_Cumulative_p()
{
    double sumoffiness=0.0,sumofp=0.0;
    int i;
    vector<double> probability;
    probability.reserve(MAXCELLNUMBER);
    
    for(i=0;i<rk[0].GetTotalcellnumber();i++)
    {
        //cout<<rk[i].GetFitness_Target_parasite_average()<<endl;
        sumoffiness += rk[i].GetFitness_Target_parasite_average();
        //cout<<"sumoffiness="<<sumoffiness<<endl;
    }
    if(sumoffiness>0.0)
    {
        for(i=0;i<rk[0].GetTotalcellnumber();i++)
        {
            probability.push_back(rk[i].GetFitness_Target_parasite_average()/sumoffiness);
            sumofp = sumofp + probability[i];
            rk[i].SetCumulativ_p(sumofp);
            rk[i].SetFitness_p(probability[i]);
            //cout<<"sumofp="<<sumofp<<endl;
        }
    }
    else{//全fitnessが0の時
        rk[i].SetCumulativ_p((double)(1.0/rk[0].GetTotalcellnumber()));
    }
}

void CELL::Calculate_Cumulative_p_power()//fitnessから親選びのための累積確率を作る 10**fitness
{
    double sumoffiness=0.0,sumofp=0.0,fit=0.0;
    int i;
    vector<double> probability;
    vector<double> power_fitness;
    probability.reserve(MAXCELLNUMBER);
    power_fitness.reserve(MAXCELLNUMBER);
    
    for(i=0;i<rk[0].GetTotalcellnumber();i++)
    {
        //cout<<rk[i].GetFitness_Target_parasite_average()<<endl;
        fit=rk[i].GetFitness_Target_parasite_average();
        if(fit<=10||fit>=0)
        {
            power_fitness.push_back(pow(10.0,fit));
            sumoffiness += power_fitness[i];
            //cout<<"sumoffiness="<<sumoffiness<<endl;
        }
    }
    
    if(sumoffiness>0.0)
    {
        for(i=0;i<rk[0].GetTotalcellnumber();i++)
        {
            probability.push_back(power_fitness[i]/sumoffiness);
            sumofp = sumofp + probability[i];
            rk[i].SetCumulativ_p(sumofp);
            rk[i].SetFitness_p(probability[i]);
            //cout<<"sumofp="<<sumofp<<endl;
        }
    }
    else{
        //全fitnessが0の時
        for(i=0;i<rk[0].GetTotalcellnumber();i++)
        {
            rk[i].SetCumulativ_p((double)(1.0/rk[0].GetTotalcellnumber()));
        }
    }
}

int  CELL::Roulette()
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    
    double rd=0.0;
    int i=0,parent_num=0;
    
    rd = rand01(mt);
    
    while(rd>rk[i].GetCumulative_p())
    {
        parent_num = i;
        i++;
    }
    
    return parent_num;
}

void CELL::elite(int time)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_int_distribution<> randnonelite(0,MAXCELLNUMBER-1-elite_num);//上位以外
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[mutation_num] ={};
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    parentnum[0] =randnonelite(mt);
    
    for(i=1;i<mutation_num;i++)
    {
        while(parentnum[i+1]==parentnum[i])
        {
            parentnum[i] = randnonelite(mt);
        }
    }
    for(i=0;i<mutation_num;i++)
    {
        eliminate_cellnumber = randnonelite(mt);
        while(eliminate_cellnumber==parentnum[i])
        {
            eliminate_cellnumber =randnonelite(mt);//親子が一致しないように振り直し続ける
        }
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        child1 = randProtein(mt);
        child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
        child3 = randProtein(mt);
        child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child1 = randProtein(mt);//inputは変化しない
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)>=10)//0のpathを探す。ただし総本数が10本を超えないように
        {
            child3 = randProtein(mt);//inputは変化しない
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
                      rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                }
            }
        }
        
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=4||count2>=4)
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
                          rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
        }
        
        j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
        j34 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4);//0
        
        if(j12 != j34)
        {
            
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
            rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            //変異体のカウントを増やす
            check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
        }
    }
}

void CELL::elite_noinput(int time)//inputは変異しない
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_int_distribution<> randnonelite(0,MAXCELLNUMBER-1-elite_num);//上位以外
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_TI(target_number+1,ProteinVAL-1-input_number);
    uniform_int_distribution<> randProtein_non_I(0,ProteinVAL-1-input_number);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[mutation_num] ={};
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    parentnum[0] =randnonelite(mt);
    
    for(i=1;i<mutation_num;i++)
    {
        while(parentnum[i+1]==parentnum[i])
        {
            parentnum[i] = randnonelite(mt);
        }
    }
    for(i=0;i<mutation_num;i++)
    {
        eliminate_cellnumber = randnonelite(mt);
        while(eliminate_cellnumber==parentnum[i])
        {
            eliminate_cellnumber =randnonelite(mt);//親子が一致しないように振り直し続ける
        }
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        child1 = randProtein_non_I(mt);
        child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
        child3 = randProtein_non_I(mt);
        child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child1 = randProtein_non_I(mt);
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)>=10)//0のpathを探す。ただし総本数が10本を超えないように
        {
            child3 = randProtein_non_I(mt);//inputは変化しない
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
                      rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_TI(mt);//inputは変化しない,input-targetの直結は禁止
                }
            }
        }
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=4||count2>=4)//targetに繋ぎすぎていれば
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
                          rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_TI(mt);//inputは変化しない
                    }
                }
            }
        }
        
        j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
        j34 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4);//0
        
        if(j12 != j34)
        {
            
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
            rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            //変異体のカウントを増やす
            check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
        }
    }
}

void CELL::elite_one(int time)//エリート戦略
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_int_distribution<> randnonelite(0,MAXCELLNUMBER-1-elite_num);//上位以外
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,count1=0,count2=0,j12=0,j32=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[mutation_num] ={};
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    if(elite_num == 0)
    {
        parentnum[0] =Roulette();
        for(i=1;i<mutation_num;i++)
        {
            while(parentnum[i+1]==parentnum[i])
            {
                parentnum[i] = Roulette();
            }
        }
    }
    else{
        parentnum[0] =randnonelite(mt);
        for(i=1;i<mutation_num;i++)
        {
            while(parentnum[i+1]==parentnum[i])
            {
                parentnum[i] = randnonelite(mt);
            }
        }
    }
    
    
    for(i=1;i<mutation_num;i++)
    {
        while(parentnum[i+1]==parentnum[i])
        {
            parentnum[i] = randnonelite(mt);
        }
    }
    
    for(i=0;i<mutation_num;i++)
    {
        eliminate_cellnumber = randnonelite(mt);
        while(eliminate_cellnumber==parentnum[i])
        {
            eliminate_cellnumber =randnonelite(mt);//親子が一致しないように振り直し続ける
        }
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        child1 = randProtein(mt);//TargetのFBや自己複製を禁止 pathがある方
        child2 = randProtein_non_T(mt);
        child3 = randProtein(mt);//TargetのFBや自己複製を禁止 pathがない方
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)>=10)//0のpathを探す。ただし総本数が10本を超えないように
        {
            child3 = randProtein(mt);//TargetのFBや自己複製を禁止
        }
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child2 = randProtein(mt);
            child1 = randProtein(mt);//TargetのFBや自己複製を禁止
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
                      rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                }
            }
        }
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=4||count2>=4)
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
        }
        
        j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
        j32 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2);//0
        
        if(j12 != j32)
        {
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j32,child1,child2);//path切り替え 2--->1を2--->3に
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child2);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
            rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            //変異体のカウントを増やす
            check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
        }
    }
}

void CELL::elite_noinput_one(int time)//エリート戦略inputは変異しない
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_int_distribution<> randnonelite(0,MAXCELLNUMBER-1-elite_num);//上位以外
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_module(target_number+1, ProteinVAL-1-2);
    uniform_int_distribution<> randProtein_non_T(target_number+1,ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_TI(target_number+1,ProteinVAL-1-input_number);
    uniform_int_distribution<> randProtein_non_I(0,ProteinVAL-1-input_number);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,count1=0,count2=0,j12=0,j32=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[mutation_num] ={};
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    parentnum[0] =randnonelite(mt);
    
    for(i=1;i<mutation_num;i++)
    {
        while(parentnum[i+1]==parentnum[i])
        {
            parentnum[i] = randnonelite(mt);
        }
    }
    for(i=0;i<mutation_num;i++)
    {
        eliminate_cellnumber = randnonelite(mt);
        while(eliminate_cellnumber==parentnum[i])
        {
            eliminate_cellnumber =randnonelite(mt);//親子が一致しないように振り直し続ける
        }
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        child2 = randProtein_non_module(mt);//TargetのFBや自己複製を禁止 pathがある方
        child3 = randProtein_non_I(mt);
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)>=10||rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)<2)//0のpathを探す。ただし総本数が10本を超えないように
        {
            child2 = randProtein_non_module(mt);//targetは選ばない turing-moduleも
        }
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0)//0のpathを探す。
        {
            child3 = randProtein_non_I(mt);//TargetのFBや自己複製を禁止
        }
        
        child1 = randProtein_non_I(mt);
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child1 = randProtein_non_I(mt);//inputは変化しない
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_TI(mt);//inputは変化しない,input-targetの直結は禁止
                }
            }
        }
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=4||count2>=4)//targetに繋ぎすぎていれば
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_TI(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
        }
        
        j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
        j32 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2);//0
        
        if(j12 != j32)
        {
            
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j32,child1,child2);//path切り替え 2--->1を2--->3に
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child2);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
            rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            //変異体のカウントを増やす
            check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
        }
    }
}
void CELL::Selection_Pressure_one(int time)//選択圧の変更  上位の個体だけを残す
{
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,count1=0,count2=0,j12=0,j32=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[MAXCELLNUMBER-selection_num] ={};//下位個体の入れ物
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_int_distribution<> randnonelite(0,selection_num-1);//上位以外
    uniform_int_distribution<> randelite(selection_num,MAXCELLNUMBER-1);//上位個体
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_module(target_number, ProteinVAL-1-2);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_TI(target_number,ProteinVAL-1-input_number);
    uniform_int_distribution<> randProtein_non_I(0,ProteinVAL-1-input_number);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    for(i=0;i<MAXCELLNUMBER-1-selection_num;i++)
    {
        parentnum[i] = randelite(mt);
    }
    
    for(i=0;i<MAXCELLNUMBER-1-selection_num;i++)
    {
        eliminate_cellnumber = i;
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        rk[data_array_fitness[eliminate_cellnumber].num].calculate_Input_list();
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        child1 = randProtein(mt);//TargetのFBや自己複製を禁止 pathがある方
        child2 = randProtein_non_T(mt);
        child3 = randProtein(mt);//TargetのFBや自己複製を禁止 pathがない方
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=10 ||
              rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)>=10)//0のpathを探す。ただし総本数が10本を超えないように
        {
            child3 = randProtein(mt);//TargetのFBや自己複製を禁止
        }
        
        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child2 = randProtein(mt);
            child1 = randProtein(mt);//TargetのFBや自己複製を禁止
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child2 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                }
            }
        }
        
        for(int d=0;d<depth;d++)
        {
            for(int l=0;l<rk[data_array_fitness[eliminate_cellnumber].num].GetInput_listSize(d);l++)
            {
                if(child2 == rk[data_array_fitness[eliminate_cellnumber].num].GetInput_list(d,l))////inputから距離i+1の遺伝子が選ばれた時
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
        }
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=4||count2>=4)
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
        }
        
        j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
        j32 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2);//0
        
        if(j12 != j32)
        {
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j32,child1,child2);//path切り替え 2--->1を2--->3に
            rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child2);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
            rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
            rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
            //変異体のカウントを増やす
            check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
        }
    }
    
}

void CELL::Selection_Pressure_noinput_one(int time)//選択圧の変更
{
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,count1=0,count2=0,j12=0,j32=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[MAXCELLNUMBER-selection_num] ={};//下位個体の入れ物
    double p=0.0;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_real_distribution<> Probability(0,1);
    uniform_int_distribution<> randnonelite(0,selection_num-1);//上位以外
    uniform_int_distribution<> randelite(selection_num,MAXCELLNUMBER-1);//上位個体
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();//ここで時間平均を使うかを与える
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            //fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetParasite_attacked_number());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    for(i=0;i<MAXCELLNUMBER-1-selection_num;i++)
    {
        parentnum[i] = randelite(mt);//親選び
    }
    
    for(i=selection_num-1;i<MAXCELLNUMBER;i++)
    {
        //parentnum[i] = i;//親は上位のみ
    }
    
    for(i=0;i<MAXCELLNUMBER-1-selection_num;i++)
    {
        eliminate_cellnumber = i;//子供は下位のcellを削除して増やす
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増やす"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        //rk[data_array_fitness[eliminate_cellnumber].num].calculate_Input_list();
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        p = Probability(mt);//Hostのmutation_rate
        
        if(p<Host_mutaion_rate)
        {
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
            child3 = randProtein(mt);
            
            while(rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)>=MaxPathNum+1||rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)<2)//0のpathを探す。ただし総本数が9本を超えないように
            {
                child2 = randProtein_non_T(mt);//targetは選ばない turing-moduleも
            }
            
            while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0)//0のpathを探す。
            {
                child3 = randProtein(mt);//TargetのFBや自己複製を禁止
            }
            
            child1 = randProtein(mt);
            while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
            {
                child1 = randProtein(mt);//inputは変化しない
            }
            
            //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
            //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child2)<<endl;
            
            for(int l=0;l<input_number;l++)
            {
                if(child3 == ProteinVAL-1-l)//input-geneが選ばれた時
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=MaxPathNum ||
                          rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child3) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
            
            for(int t=0;t<target_number;t++)
            {
                if(child3 == t)//繋ぎ変える先がtargetの時
                {
                    count1=0;count2=0;
                    for(int j=0;j<target_number;j++)
                    {
                        if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == 1)
                            //すでに2からtargetがactivateされている時
                        {
                            count1 = count1 + 1;
                        }
                        else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child2) == -1)
                        {
                            count2 = count2 +1;
                        }
                    }
                    if(count1>=3||count2>=3)//targetに繋ぎすぎていれば　8/29:4本から3本へ変更
                    {
                        while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2) != 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                        }
                    }
                }
            }
            
            j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
            j32 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child2);//0
            
            if(j12 != j32)
            {
                
                //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j32,child1,child2);//path切り替え 2--->1を2--->3に
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child2);//path切り替え 2--->3を2--->1に
                
                //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
                rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                //変異回数のカウントを増やす
                check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
            }
        }
    }
}

void CELL::Selection_Pressure_norule(int time)
{
    int i=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[MAXCELLNUMBER-selection_num] ={};//下位個体の入れ物
    double p=0.0;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_real_distribution<> Probability(0,1);
    uniform_int_distribution<> randnonelite(0,selection_num-1);//上位以外
    uniform_int_distribution<> randelite(selection_num,MAXCELLNUMBER-1);//上位個体
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0;
        data_array_fitness[ith].fit =0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();//ここで時間平均を使うかを与える
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            //fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetParasite_attacked_number());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            //fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness());
            fprintf(pf,"\n");
        }
    }
    fclose(pf);
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        parentnum[i] = randelite(mt);//親選び
    }
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        eliminate_cellnumber = i;//子供は下位のcellを削除して増やす
        
        mother_celltype = rk[data_array_fitness[parentnum[i]].num].GetCelltype();//親のcelltypeを引き継ぎ
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増す"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(rk[data_array_fitness[parentnum[i]].num].GetJ(l, k), l, k);
            }
        }
        
        rk[data_array_fitness[eliminate_cellnumber].num].SetStationaryTime(rk[data_array_fitness[parentnum[i]].num].GetTime());//分裂時間を更新
        rk[data_array_fitness[eliminate_cellnumber].num].converterJ();
        //rk[data_array_fitness[eliminate_cellnumber].num].calculate_Input_list();
        rk[data_array_fitness[eliminate_cellnumber].num].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[data_array_fitness[parentnum[i]].num].SetDivision_number(rk[data_array_fitness[parentnum[i]].num].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history.clear();
        rk[data_array_fitness[eliminate_cellnumber].num].tree_history = rk[data_array_fitness[parentnum[i]].num].tree_history;
        
        for(int k=0;k<input_number;k++)
        {
            //rk[data_array_fitness[eliminate_cellnumber].num].SetYin(rk[ith].gauss_rand(rk[ith].GetYin(k),sigma),k);
            //rk[data_array_fitness[eliminate_cellnumber].num].SetY(rk[ith].GetY(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetX(rk[ith].GetX(k), k);//medium中のinput物質は引き継ぎ
            //rk[data_array_fitness[eliminate_cellnumber].num].SetMedium_D(rk[ith].GetMedium_D(k), k);
            for(int l=0;l<ProteinVAL;l++)
            {
                //rk[data_array_fitness[eliminate_cellnumber].num].SetYin_J(rk[ith].GetYin_J(k, l), k, l);
            }
        }
        
        p = Probability(mt);//Hostのmutation_rate
        
        if(p<Host_mutaion_rate)
        {
            
            cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増す"<<data_array_fitness[parentnum[i]].num<<"fitness="<<rk[data_array_fitness[parentnum[i]].num].GetFitness_Target_parasite_average()<<endl;
            child1 = randProtein(mt);
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
            child3 = randProtein(mt);
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
            
            while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
            {
                child1 = randProtein(mt);//inputは変化しない
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            while(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 ||
                  rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=MaxPathNum ||
                  rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
            {
                child3 = randProtein(mt);//inputは変化しない
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
            //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
            
            for(int l=0;l<input_number;l++)
            {
                if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                {
                    while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=MaxPathNum ||
                          rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                }
            }
            
            for(int t=0;t<target_number;t++)
            {
                if(child3 == t)//繋ぎ変える先がtargetの時
                {
                    count1=0;count2=0;
                    for(int j=0;j<target_number;j++)
                    {
                        if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == 1)
                            //すでに2からtargetがactivateされている時
                        {
                            count1 = count1 + 1;
                        }
                        else if(rk[data_array_fitness[eliminate_cellnumber].num].GetJ(j,child4) == -1)
                        {
                            count2 = count2 +1;
                        }
                    }
                    if(count1>=4||count2>=4)
                    {
                        while(rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)>=MaxPathNum ||
                              rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4) != 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//繋ぎすぎの時は中間層に
                        }
                    }
                }
            }
            
            j12 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1,child2);//-1 or 1
            j34 = rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3,child4);//0
            
            if(j12 != j34)
            {
                //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                rk[data_array_fitness[eliminate_cellnumber].num].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                
                //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                rk[data_array_fitness[eliminate_cellnumber].num].converterJ();//隣接リストにmutationを反映
                rk[data_array_fitness[eliminate_cellnumber].num].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                //変異回数のカウントを増やす
                check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
            }
        }
    }
}

void CELL::Selection_WF_norule(int time)
{
    int i=0,j=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[MAXCELLNUMBER-selection_num] ={};//下位個体の入れ物
    double rd[MAXCELLNUMBER-selection_num] ={};
    
    double p=0.0;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_real_distribution<> Probability(0,1);
    uniform_int_distribution<> randelite(selection_num,MAXCELLNUMBER-1);//上位個体
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);//target以外
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0.0;
        data_array_fitness[ith].fit =0.0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();//ここで時間平均を使うかを与える
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            //fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetParasite_attacked_number());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness_Target_noparasite_average());
            fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness_p());//子孫を残す確率
            
            fprintf(pf,"\n");
            //cout<<"p["<<ith<<"]="<<rk[ith].GetCumulative_p()<<endl;
        }
    }
    fclose(pf);
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        rd[i] = rand01(mt);
        j=0;
        //cout<<"rd["<<i<<"]="<<rd[i]<<endl;
        
        while(rd[i]>rk[j].GetCumulative_p())
        {
            j++;
        }
        parentnum[i] = j;
        //cout<<"Cumulative_p["<<parentnum[i]-1<<"]="<<rk[parentnum[i]-1].GetCumulative_p()<<endl;
        //cout<<"Cumulative_p["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetCumulative_p()<<endl;
        //cout<<"fitness_p["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetFitness_p()<<endl;
        //cout<<"fitness["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetFitness_Target_parasite_average()<<endl;
    }
    
    const auto arr = make_rand_array_unique(MAXCELLNUMBER-selection_num, MAXCELLNUMBER-1, 0);
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        //eliminate_cellnumber = randcell(mt);//randomに取り除く
        eliminate_cellnumber = arr[i];
        //cout<<"eliminate_cellnumber="<<eliminate_cellnumber<<endl;
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増す"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[eliminate_cellnumber].SetJ(rk[parentnum[i]].GetJ(l, k), l, k);
            }
        }
        
        rk[eliminate_cellnumber].SetStationaryTime(rk[parentnum[i]].GetTime());
        //分裂時間を更新
        rk[eliminate_cellnumber].converterJ();
        //rk[data_array_fitness[eliminate_cellnumber].num].calculate_Input_list();
        mother_celltype = rk[parentnum[i]].GetCelltype();//親のcelltypeを引き継ぎ
        rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[eliminate_cellnumber].SetDivision_number(rk[parentnum[i]].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[eliminate_cellnumber].tree_history.clear();
        rk[eliminate_cellnumber].tree_history = rk[parentnum[i]].tree_history;
        
        p = Probability(mt);//Hostのmutation_rate
        
        if(p<Host_mutaion_rate)
        {
            cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増す"<<parentnum[i]<<"fitness="<<rk[parentnum[i]].GetFitness_Target_parasite_average()<<endl;
            child1 = randProtein(mt);
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
            child3 = randProtein(mt);
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
            
            while(rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
            {
                child1 = randProtein(mt);//inputは変化しない
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            while(rk[eliminate_cellnumber].GetJ(child3,child4) != 0 ||
                  rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                  rk[eliminate_cellnumber].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
            {
                child3 = randProtein(mt);//inputは変化しない
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
            //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
            
            for(int l=0;l<input_number;l++)
            {
                if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                {
                    //cout<<"Input_child4="<<child4<<endl;
                    //cout<<"old_child3="<<child3<<endl;
                    child3 = randProtein_non_T(mt);
                    while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                          rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                    //cout<<"Input_child4="<<child4<<endl;
                    //cout<<"new_child3="<<child3<<endl;
                    break;
                }
            }
            
            if(on==1)//使わない
            {
                for(int l=0;l<Parasite_genome_size;l++)//３種
                {
                    if(child2 == l+Parasite_target_genome||child1 == l)//attcked-geneとtargetを必ず繋いでるpathだった時
                    {
                        while(rk[eliminate_cellnumber].GetSEL(child1)>=MaxPathNum ||
                              rk[eliminate_cellnumber].GetJ(child1,child2) != 0 )//別の0のpathを探す
                        {
                            child1 = randProtein_non_T(mt);//targetでないところから探す
                        }
                        break;
                    }
                }
            }
            
            for(int t=0;t<target_number;t++)
            {
                if(child3 == t)//繋ぎ変える先がtargetの時
                {
                    count1=0;count2=0;
                    for(int j=0;j<target_number;j++)
                    {
                        if(rk[eliminate_cellnumber].GetJ(j,child4) == 1)
                            //すでに2からtargetがactivateされている時
                        {
                            count1 = count1 + 1;
                        }
                        else if(rk[eliminate_cellnumber].GetJ(j,child4) == -1)
                        {
                            count2 = count2 +1;
                        }
                    }
                    //if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                    if(count1>=MaxTargetPath)
                    {
                        child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                        while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                              rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                        }
                    }
                    break;
                }
            }
            
            
            j12 = rk[eliminate_cellnumber].GetJ(child1,child2);//-1 or 1
            j34 = rk[eliminate_cellnumber].GetJ(child3,child4);//0
            
            if(j12 != j34)
            {
                //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                
                rk[eliminate_cellnumber].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                rk[eliminate_cellnumber].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                
                //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                rk[eliminate_cellnumber].converterJ();//隣接リストにmutationを反映
                rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                //変異回数のカウントを増やす
                //check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
            }
        }
    }
}

void CELL::Selection_WF_self_regulation(int time)//WFでtargetの自己制御を許す
{
    int i=0,j=0,totalcellnum=0,eliminate_cellnumber=0,mother_celltype=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    int parentnum[MAXCELLNUMBER-selection_num] ={};//下位個体の入れ物
    double rd[MAXCELLNUMBER-selection_num] ={};
    
    double p=0.0;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0,1);
    uniform_real_distribution<> Probability(0,1);
    uniform_int_distribution<> randelite(selection_num,MAXCELLNUMBER-1);//上位個体
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);//target以外
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    
    for(int ith=0;ith<totalcellnum;ith++)//sort用のdataを初期化
    {
        data_array_fitness[ith].hammD =0.0;
        data_array_fitness[ith].fit =0.0;
        data_array_fitness[ith].num =0;//sort用の入れ物
    }
    SetData_array_fitness();//ここで時間平均を使うかを与える
    
    for(int ith=0;ith<totalcellnum;ith++)//Xith
    {
        data_array_fitness[ith].num = ith;//sort前検索用の番号をセット
    }
    
    sort(data_array_fitness.begin(), data_array_fitness.begin()+totalcellnum);
    
    FILE *pf = nullptr;
    
    int _number=0;
    
    _number  = rk[0].GetTotalcellnumber() ;//同じ細胞を同一ファイルに記録
    
    pf = fopen("fitness_sort.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        for(int ith=0 ;ith<_number;ith++)
        {
            fprintf(pf," %d ",time);
            fprintf(pf," %d " ,ith);
            fprintf(pf," #%d " ,data_array_fitness[ith].num);
            //fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetCelltype());//celltypeを記録
            fprintf(pf," %d " ,rk[data_array_fitness[ith].num].GetParasite_attacked_number());//celltypeを記録
            //fprintf(pf," %s ",rk[data_array_fitness[ith].num].tree_history.c_str());//変異数を確認
            fprintf(pf," %f",data_array_fitness[ith].fit);
            fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness_Target_noparasite_average());
            fprintf(pf," %f",rk[data_array_fitness[ith].num].GetFitness_p());//子孫を残す確率
            
            fprintf(pf,"\n");
            //cout<<"p["<<ith<<"]="<<rk[ith].GetCumulative_p()<<endl;
        }
    }
    fclose(pf);
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        rd[i] = rand01(mt);
        j=0;
        //cout<<"rd["<<i<<"]="<<rd[i]<<endl;
        
        while(rd[i]>rk[j].GetCumulative_p())
        {
            j++;
        }
        parentnum[i] = j;
        //cout<<"Cumulative_p["<<parentnum[i]-1<<"]="<<rk[parentnum[i]-1].GetCumulative_p()<<endl;
        //cout<<"Cumulative_p["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetCumulative_p()<<endl;
        //cout<<"fitness_p["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetFitness_p()<<endl;
        //cout<<"fitness["<<parentnum[i]<<"]="<<rk[parentnum[i]].GetFitness_Target_parasite_average()<<endl;
    }
    
    const auto arr = make_rand_array_unique(MAXCELLNUMBER-selection_num, MAXCELLNUMBER-1, 0);
    
    for(i=0;i<MAXCELLNUMBER-selection_num;i++)
    {
        //eliminate_cellnumber = randcell(mt);//randomに取り除く
        eliminate_cellnumber = arr[i];
        //cout<<"eliminate_cellnumber="<<eliminate_cellnumber<<endl;
        //cout<<"消す"<<data_array_fitness[eliminate_cellnumber].num<<"----->"<<"増す"<<data_array_fitness[parentnum[i]].num<<endl;
        
        for(int l=0;l<ProteinVAL;l++)
        {
            for(int k=0;k<ProteinVAL;k++)
            {
                rk[eliminate_cellnumber].SetJ(rk[parentnum[i]].GetJ(l, k), l, k);
            }
        }
        
        rk[eliminate_cellnumber].SetStationaryTime(rk[parentnum[i]].GetTime());
        //分裂時間を更新
        rk[eliminate_cellnumber].converterJ();
        //rk[data_array_fitness[eliminate_cellnumber].num].calculate_Input_list();
        mother_celltype = rk[parentnum[i]].GetCelltype();//親のcelltypeを引き継ぎ
        rk[eliminate_cellnumber].SetCelltype(mother_celltype);//ithのcelltype引き継ぎ
        //cout<<"type(daughter)="<<rk[eliminate_cellnumber].GetCelltype()<<endl;
        rk[eliminate_cellnumber].SetDivision_number(rk[parentnum[i]].GetDivision_number() + 1);//親のカウントを増やす
        
        rk[eliminate_cellnumber].tree_history.clear();
        rk[eliminate_cellnumber].tree_history = rk[parentnum[i]].tree_history;
        
        p = Probability(mt);//Hostのmutation_rate
        
        if(p<Host_mutaion_rate)
        {
            cout<<"消す"<<eliminate_cellnumber<<"----->"<<"増す"<<parentnum[i]<<"fitness="<<rk[parentnum[i]].GetFitness_Target_parasite_average()<<endl;
            child1 = randProtein(mt);
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
            child3 = randProtein(mt);
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
            
            while(rk[eliminate_cellnumber].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
            {
                child1 = randProtein(mt);//inputは変化しない
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            while(rk[eliminate_cellnumber].GetJ(child3,child4) != 0 ||
                  rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                  rk[eliminate_cellnumber].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
            {
                child3 = randProtein(mt);//inputは変化しない
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
            //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
            
            for(int l=0;l<input_number;l++)
            {
                if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                {
                    //cout<<"Input_child4="<<child4<<endl;
                    //cout<<"old_child3="<<child3<<endl;
                    child3 = randProtein_non_T(mt);
                    while(rk[eliminate_cellnumber].GetSEL(child3)>=MaxPathNum ||
                          rk[eliminate_cellnumber].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                    //cout<<"Input_child4="<<child4<<endl;
                    //cout<<"new_child3="<<child3<<endl;
                    break;
                }
            }
            
            
            j12 = rk[eliminate_cellnumber].GetJ(child1,child2);//-1 or 1
            j34 = rk[eliminate_cellnumber].GetJ(child3,child4);//0
            
            if(j12 != j34)
            {
                //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                
                rk[eliminate_cellnumber].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                rk[eliminate_cellnumber].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                
                //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child1, child2)<<endl;
                //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[data_array_fitness[eliminate_cellnumber].num].GetJ(child3, child4)<<endl;
                rk[eliminate_cellnumber].converterJ();//隣接リストにmutationを反映
                rk[eliminate_cellnumber].SetMutation_number(rk[eliminate_cellnumber].GetMutation_number()+1);
                //変異回数のカウントを増やす
                //check_tree_history(data_array_fitness[parentnum[i]].num,data_array_fitness[eliminate_cellnumber].num);
            }
        }
    }
}

void CELL::reset()
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    
    //uniform_real_distribution<> rand01(0, 1);
    uniform_real_distribution<> rand0th(0,_threshold);
    
    for(int i=0;i<rk[0].GetTotalcellnumber();i++)
    {
        for(int l=0;l<rk[0].Getprotein_x_number();l++)
        {
            //rk[i].SetResult(rk[ith].gauss_rand(rk[ith].GetResult(l),sigma),l);
            rk[i].SetResult(rand0th(mt),l);//閾値以下で乱数
            //rk[i].SetResult(rand01(mt),l);
            rk[i].SetResult_average(0.0,l);
            
        }
        for(int k=0;k<input_number;k++)
        {
            rk[i].SetYin(0,k);//分裂時に子供のinput物質は0に
            rk[i].SetY(S_resource,k);
        }
        rk[i].SetFitness_Target_parasite_average(0.0);
        
        rk[i].SetTime(0.0);
    }
}

void CELL::Calculate_Efficiency()//利用効率の計算とfitnessの計算
{
    int i,j;
    double fit=0.0,S=0.0;
    vector<double> square;
    square.reserve(target_number);
    S = (double)(MAXCELLNUMBER/target_number);
    
    for(j=0;j<target_number;j++)
    {
        square[j] =0.0;
    }
    for(j=0;j<target_number;j++)
    {
        for(i=0;i<rk[0].GetTotalcellnumber();i++)
        {
            square[j] = square[j] + rk[i].GetResult(j)*rk[i].GetResult(j);
        }
    }
    
    for(i=0;i<rk[0].GetTotalcellnumber();i++)
    {
        for(j=0;j<target_number;j++)
        {
            if(square[j]==0)//分母が0の時を除く
            {
                fit = fit + 0.0;
            }
            else{
                fit = fit + S*rk[i].GetResult(j)*rk[i].GetResult(j)/square[j];
            }
        }
        rk[i].SetFitness(fit);
        fit =0.0;
    }
}

int CELL::GetHistogram(int k,int j)
{
    return histogram[k][j];
}

void CELL::SetHistogram(int his,int k,int j)
{
    if(his>MAXCELLNUMBER)
    {
        cout<<"his ="<<his<<endl;
        exit(-1);
    }
    else
    {
        histogram[k][j]=his;
    }
}

void CELL::histo(int time)
{
    int i, j, k,gene=0;
    int  width=histdiv;
    double y_max=1.0,y_min=0.0,dw=0.0;//0から１の範囲
    double u_limit, l_limit;
    
    dw = (double)((y_max-y_min)/width);
    gene = time/500;
    
    vector<vector<int> > histo;
    
    histo.shrink_to_fit();
    histo.resize(histdiv);
    
    for(int i=0;i<histdiv;i++)
    {
        histo[i].resize(target_number+Parasite_genome_size);
    }
    for(k=0;k<width;k++)
    {
        for(j=0;j<target_number+Parasite_genome_size;j++)
        {
            histo[k][j] = 0;
        }
    }
    
    for(j=0;j<target_number;j++)
    {
        u_limit = (y_min + dw);
        l_limit =  y_min;
        for(k=0;k<width;k++)
        {
            for(i=0;i<rk[0].GetTotalcellnumber();i++)
            {
                if(rk[i].GetResult(j) < u_limit && rk[i].GetResult(j) >= l_limit)
                {
                    histo[k][j] = histo[k][j] +1;
                }
            }
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
    }
    
    for(j=target_number;j<target_number + Parasite_genome_size;j++)//inputのhistogram
    {
        u_limit = (y_min + dw);
        l_limit =  y_min;
        for(k=0;k<width;k++)
        {
            for(i=0;i<rk[0].GetTotalcellnumber();i++)
            {
                if(rk[i].GetResult(j) < u_limit && rk[i].GetResult(j) >= l_limit)
                {
                    histo[k][j] = histo[k][j] +1;
                }
            }
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
    }
    
    for(j=0;j<target_number;j++)
    {
        FILE *pf = nullptr;
        char filename[50]={};
        sprintf(filename,"histo%d.txt",j);
        pf = fopen(filename,"a+");
        u_limit = (y_min + dw);
        l_limit =  y_min;
        
        fprintf(pf,"%d %f %d ",gene-1,l_limit,0);//1世代前の
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,l_limit,histogram[0][j]);//1世代前の
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,u_limit,histogram[0][j]);
        fprintf(pf,"\n");
        l_limit = u_limit;
        u_limit = u_limit +dw;
        
        for(k=1;k<width;k++)
        {
            fprintf(pf,"%d %f %d ",gene-1,l_limit,histogram[k][j]);//1世代前の
            fprintf(pf,"\n");
            fprintf(pf,"%d %f %d ",gene-1,u_limit,histogram[k][j]);
            fprintf(pf,"\n");
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
        fprintf(pf,"\n");
        u_limit = (y_min + dw);
        l_limit =  y_min;
        
        fprintf(pf,"%d %f %d ",gene-1,l_limit,0);
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,l_limit,histo[0][j]);
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,u_limit,histo[0][j]);
        fprintf(pf,"\n");
        l_limit = u_limit;
        u_limit = u_limit +dw;
        for(k=1;k<width;k++)
        {
            fprintf(pf,"%d %f %d ",gene-1,l_limit,histo[k][j]);
            fprintf(pf,"\n");
            fprintf(pf,"%d %f %d ",gene-1,u_limit,histo[k][j]);
            fprintf(pf,"\n");
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
        fprintf(pf,"\n");
        fclose(pf);
    }
    
    for(j=target_number;j<target_number + Parasite_genome_size;j++)
    {
        FILE *pf = nullptr;
        char filename[50]={};
        sprintf(filename,"histo%d.txt",j);
        pf = fopen(filename,"a+");
        u_limit = (y_min + dw);
        l_limit =  y_min;
        
        fprintf(pf,"%d %f %d ",gene-1,l_limit,0);//1世代前の
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,l_limit,histogram[0][j]);//1世代前の
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,u_limit,histogram[0][j]);
        fprintf(pf,"\n");
        l_limit = u_limit;
        u_limit = u_limit +dw;
        
        for(k=1;k<width;k++)
        {
            fprintf(pf,"%d %f %d ",gene-1,l_limit,histogram[k][j]);//1世代前の
            fprintf(pf,"\n");
            fprintf(pf,"%d %f %d ",gene-1,u_limit,histogram[k][j]);
            fprintf(pf,"\n");
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
        fprintf(pf,"\n");
        u_limit = (y_min + dw);
        l_limit =  y_min;
        
        fprintf(pf,"%d %f %d ",gene-1,l_limit,0);
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,l_limit,histo[0][j]);
        fprintf(pf,"\n");
        fprintf(pf,"%d %f %d ",gene-1,u_limit,histo[0][j]);
        fprintf(pf,"\n");
        l_limit = u_limit;
        u_limit = u_limit +dw;
        for(k=1;k<width;k++)
        {
            fprintf(pf,"%d %f %d ",gene-1,l_limit,histo[k][j]);
            fprintf(pf,"\n");
            fprintf(pf,"%d %f %d ",gene-1,u_limit,histo[k][j]);
            fprintf(pf,"\n");
            l_limit = u_limit;
            u_limit = u_limit +dw;
        }
        fprintf(pf,"\n");
        fclose(pf);
    }
    
    for(k=0;k<width;k++)
    {
        for(j=0;j<target_number+Parasite_genome_size;j++)
        {
            histogram[k][j] = histo[k][j];
        }
    }
}

template <typename T> std::string tostr(const T& t)
{
    std::ostringstream os; os<<t; return os.str();
}

void CELL::Calculate_VipVg(int time)//ここでVipVgの計算をまとめて行う
{
    int cellnumber=rk_VipVg[0].GetTotalcellnumber();
    int divnumber= time;//timeは1秒の刻み数をかけたもの
    
    vector<vector<double> > g_rand;
    g_rand.resize(EM_h);
    for(int i=0;i<EM_h;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            g_rand[i].push_back(0.0);
        }
    }
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    
    Load_Jji_random();//いくつかの遺伝子集団を選ぶ
    SetFlag_of_VipVg(0);//Vip計算に戻す
    
    for(int t=1;t<(2*divnumber)+1;t++)
    {
        for(int i=0;i<cellnumber;i++)
        {
            for(int j=0;j<EM_h;j++)
            {
                for(int l=0;l<ProteinVAL;l++)
                {
                    g_rand[j][l] = Sigma*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));//乱数作成をここで行う
                    //g_rand[j][l] = cell.rk[i].GetSigma_Str()*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));//乱数作成をここで
                }
            }
            rk_VipVg[i].calculateNS_protein_g_EM_VipVg(g_rand);
            //cell.rk[i].calculateNS_protein_signal();
            //cell.rk[i].calculateNS_protein_signal_g_EM(g_rand);
        }
       
        if(t%(divnumber)==0)//2回行う
        {
            //cell.Calculate_Variance_random();//ここで分散を計算 Vip-->Vgの切り替えここでやる
            if(Getflag_of_VipVg()==0)
            {
                Calculate_Variance_random_Vip();
            }
            else if(Getflag_of_VipVg()==2)
            {
                Calculate_Variance_random_Vg();
            }
            if(Getflag_of_VipVg()==1)//Vip計算あと
            {
                //cell.Max_fitness_Jij();//世代ごとにfitnessが最大のJijを選ぶ
                //cell.Max_fitness_Jij_host(7);
                //cell.division_fitness_mutation(sigma,cell.GetGeneration());
                //cell.archive_sim_sta(cell.GetGeneration());//division_fitness_mutationで行う
                //archive_sim_sta_Vip(GetGeneration());//division_fitness_mutationで行う
                //archive_sim(GetGeneration());//division_fitness_mutationで行う sortしない
                //cell.archive_var(cell.GetGeneration());//division_fitness_mutationで行う
                //cell.archive_var_sta(cell.GetGeneration());//division_fitness_mutationで行う
                archive_var_Vip(GetGeneration());//division_fitness_mutationで行う
                //cell.archive_var_Vip_average(cell.GetGeneration());
                //cout<<"generarion="<<50*t/(20*divinumber)<<endl;
                //cout<<"generarion="<<cell.GetGeneration()<<endl;
                //medium.reset();
                //cell.SetGeneration(cell.GetGeneration()+25);//Generationラベルを次の25世代に更新
                SetFlag_of_VipVg(2);//Vip計算に戻す
            }
            if(Getflag_of_VipVg()==3)//Vg計算あと
            {
                //cell.Max_fitness_Jij();//世代ごとにfitnessが最大のJijを選ぶ
                //cell.Max_fitness_Jij_host(7);
                //cell.division_fitness_mutation(sigma,cell.GetGeneration());
                archive_sim_sta(GetGeneration());//division_fitness_mutationで行う
                //cell.archive_sim(cell.GetGeneration());//division_fitness_mutationで行う sortしない
                //cell.archive_var(cell.GetGeneration());//division_fitness_mutationで行う
                archive_var_sta(GetGeneration());//division_fitness_mutationで行う
                //cout<<"generarion="<<50*t/(20*divinumber)<<endl;
                cout<<"generarion="<<GetGeneration()<<endl;
                //medium.reset();
                SetGeneration(GetGeneration()+10);//Generationラベルを次の10世代に更新
                SetFlag_of_VipVg(0);//Vip計算に戻す
            }
            reset();
            if(Getflag_of_VipVg()==0)//Vip
            {
                //cell.Load_Jij();
                Load_Jji_random();//いくつかの遺伝子集団を選ぶ
            }
            if(Getflag_of_VipVg()==2)//Vg計算
            {
                //cell.Vg_cell_make();//0番目のcellに一つmutationを入れた集団を作る
                Vg_cell_make_random();//0番目のcellに一つmutationを入れた集団を作る
            }
        }
    }
}

void CELL::Calculate_Variance_random_Vip()//複数の遺伝子型から分散を求める時
{
    int i,j,l,totalcellnum,OneGlpNum=0;;
    double sum,sum_sta,var[ProteinVAL]={},sum_var,var_sta[ProteinVAL]={},sum_var_sta;
    
    totalcellnum = rk[0].GetTotalcellnumber();
    OneGlpNum = totalcellnum/Glp_Num;
    
    for(l=0;l<ProteinVAL;l++)
    {
        sum_var=0.0;
        sum_var_sta=0.0;
        for(i=0;i<Glp_Num;i++)
        {
            sum =0.0;
            sum_sta=0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                sum = sum +rk_VipVg[j].GetResult(l);
                sum_sta = sum_sta +rk_VipVg[j].GetResult_average(l);
            }
            sum = (double)(sum/OneGlpNum);//平均
            sum_sta = (double)(sum_sta/OneGlpNum);//平均
            //cout<<"sum="<<sum<<endl;
            var[l]=0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                var_sta[l] =  var_sta[l] + (sum_sta -rk_VipVg[j].GetResult_average(l))*(sum_sta -rk[j].GetResult_average(l));
                var[l] =  var[l] + (sum -rk_VipVg[j].GetResult(l))*(sum -rk_VipVg[j].GetResult(l));
            }
            var[l] = (double)(var[l]/(OneGlpNum-1));//同一遺伝子集団の総数
            var_sta[l] = (double)(var_sta[l]/(OneGlpNum-1));//同一遺伝子集団の総数
            SetVariance_Gloup_Vip(i, l, var_sta[l]);//ここでグループごとの値を記録
            //cout<<"var="<<var<<endl;
            sum_var = sum_var+var[l];//gloupで平均化
            sum_var_sta = sum_var_sta+var_sta[l];//gloupで平均化
        }
        sum_var=(double)(sum_var/Glp_Num);//選んだ遺伝子の数
        //cout<<"var["<<l<<"]="<<sum_var<<endl;
        SetVariance(sum_var,l);
        sum_var_sta=(double)(sum_var_sta/Glp_Num);//選んだ遺伝子の数
        //cout<<"var_sta["<<l<<"]="<<sum_var_sta<<endl;
        SetVariance_Sta(sum_var_sta,l);
    }
    SetFlag_of_VipVg(1);//Vip計算後,記録用に切り替える
}
void CELL::Calculate_Variance_random_Vg()//複数の遺伝子型から分散を求める時
{
    int i,j,l,totalcellnum,OneGlpNum=0;;
    double sum,var[ProteinVAL]={},sum_var;
    
    totalcellnum = rk_VipVg[0].GetTotalcellnumber();
    OneGlpNum = totalcellnum/Glp_Num;
    
    for(l=0;l<ProteinVAL;l++)
    {
        sum_var=0.0;
        for(i=0;i<Glp_Num;i++)
        {
            sum =0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                //sum = sum +rk[j].GetResult(l);
                sum = sum +rk_VipVg[j].GetResult_average(l);
            }
            sum = (double)(sum/OneGlpNum);//平均
            //cout<<"sum="<<sum<<endl;
            var[l]=0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                var[l] =  var[l] + (sum -rk_VipVg[j].GetResult_average(l))*(sum -rk_VipVg[j].GetResult_average(l));
                //var[l] =  var[l] + (sum -rk[j].GetResult(l))*(sum -rk[j].GetResult(l));
            }
            var[l] = (double)(var[l]/(OneGlpNum-1));//グループでの分散
            SetVariance_Gloup_Vg(i, l, var[l]);//ここでグループごとの値を記録
            //cout<<"var="<<var<<endl;
            var[l]=GetVariance_Gloup_Vg(i, l)-GetVariance_Gloup_Vip(i, l);//それぞれの時間平均を取ったもので引いている
            if(var[l]>0.0)
            {
                sum_var = sum_var+var[l];//gloupで平均化するために和をとる
            }
        }
        sum_var=(double)(sum_var/Glp_Num);
        SetVariance(sum_var,l);
        sum_var=0.0;
        for(i=0;i<Glp_Num;i++)
        {
            sum_var=sum_var+=GetVariance_Gloup_Vg(i, l);
        }
        sum_var=(double)(sum_var/Glp_Num);
        SetVarianceVg(sum_var,l);//時間平均を取る前のVgをvariance_vgに記録
    }
    SetFlag_of_VipVg(3);//Vg計算後に切り替える
}

void CELL::Calculate_Variance_random_Vip_average(int time)
{
    if(time%(int)(rk[0].GetDiviTime())==0)//1000秒の時
    {
        double sum_var[ProteinVAL]={};
        for(int l=0;l<ProteinVAL;l++)
        {
            sum_var[l] = 0.0;
        }
        for(int ith=0;ith<Glp_Num;ith++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                SetVariance_gloup_Vip_average(ith,l,GetVariance_gloup_Vip_average(ith,l)/(average_time/rk[0].GetDelt()));
            }
            for(int l=0;l<ProteinVAL;l++)
            {
                sum_var[l] = sum_var[l]+GetVariance_gloup_Vip_average_stack(ith,l);//gloupごとを足し合わす
            }
            cout<<"Variance_gloup_Vip_average["<<ith<<"]="<<GetVariance_gloup_Vip_average(ith,4)<<endl;
        }
        
        for(int l=0;l<ProteinVAL;l++)
        {
            sum_var[l] = sum_var[l]/Glp_Num;
            SetVarianceAverage(sum_var[l],l);
        }
        //cout<<GetVarianceAverage(4)<<endl;
    }
    else
    {
        for(int ith=0;ith<Glp_Num;ith++)
        {
            for(int l=0;l<ProteinVAL;l++)
            {
                SetVariance_gloup_Vip_average(ith,l,GetVariance_gloup_Vip_average(ith,l)+GetVariance_gloup_Vip_average_stack(ith,l));
            }
            //cout<<"Variance_gloup_Vip_average["<<ith<<"]="<<GetVariance_gloup_Vip_average(ith,4)<<endl;
        }
    }
}

void CELL::Calculate_Variance_random_Vip_average_stack()//複数の遺伝子型から分散を求める時
{
    int i,j,l,totalcellnum,OneGlpNum=0;;
    double sum,var[ProteinVAL]={},sum_var;
    
    totalcellnum = rk[0].GetTotalcellnumber();
    OneGlpNum = totalcellnum/Glp_Num;
    
    for(l=0;l<ProteinVAL;l++)
    {
        sum_var=0.0;
        for(i=0;i<Glp_Num;i++)
        {
            sum =0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                sum = sum +rk_VipVg[j].GetResult(l);
            }
            sum = (double)(sum/OneGlpNum);//平均
            //cout<<"sum="<<sum<<endl;
            var[l]=0.0;
            for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                var[l] =  var[l] + (sum -rk_VipVg[j].GetResult(l))*(sum -rk_VipVg[j].GetResult(l));
            }
            var[l] = (double)(var[l]/(OneGlpNum-1));//同一遺伝子集団の総数
            SetVariance_gloup_Vip_average_stack(i, l, var[l]);//ここでグループごとの値を記録
            //cout<<"var["<<l<<"]="<<var[l]<<endl;
        }
    }
}

void CELL::SetVariance_Sta(double var,int i)
{
    if(var<0)
    {
        cout <<"var="<<var;
        exit(-1);
    }
    else
    {
        variance_sta[i] = var;
    }
}

void CELL::SetVarianceVg(double var,int i)
{
    if(var<0)
    {
        cout <<"var="<<var;
        exit(-1);
    }
    else
    {
        variance_sta_Vg[i] = var;
    }
    
}

double CELL::GetVariance_Sta(int i)
{
    return variance_sta[i];
}

double CELL::GetVarianceVg(int i)
{
    return variance_sta_Vg[i];
}

void CELL::SetVarianceAverage(double var,int i)
{
    variance_average[i]=var;
}
double CELL::GetVarianceAverage(int i)
{
    return variance_average[i];
}
void CELL::SetVariance_Gloup_Vip(int ith,int Protein_jth,double var)
{
    variance_gloup_Vip[ith][Protein_jth] = var;
}

void CELL::SetVariance_Gloup_Vg(int ith,int Protein_jth,double var)
{
    variance_gloup_Vg[ith][Protein_jth] = var;
}

void CELL::SetVariance_gloup_Vip_average(int ith,int Protein_jth,double var)
{
    variance_gloup_Vip_average[ith][Protein_jth] = var;
}
void CELL::SetVariance_gloup_Vip_average_stack(int ith,int Protein_jth,double var)
{
    variance_gloup_Vip_average_stack[ith][Protein_jth] = var;
}

double CELL::GetVariance_Gloup_Vip(int ith,int Protein_jth)
{
    return variance_gloup_Vip[ith][Protein_jth];
}

double CELL::GetVariance_Gloup_Vg(int ith,int Protein_jth)
{
    return variance_gloup_Vg[ith][Protein_jth];
}

double CELL::GetVariance_gloup_Vip_average(int ith,int Protein_jth)
{
    return variance_gloup_Vip_average[ith][Protein_jth];
}
double CELL::GetVariance_gloup_Vip_average_stack(int ith,int Protein_jth)
{
    return variance_gloup_Vip_average_stack[ith][Protein_jth];
}

void CELL::Vg_cell_make()
{
    int i=0,totalcellnum=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0;
    totalcellnum = rk[0].GetTotalcellnumber();
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);//target以外
    
    for(i=1;i<totalcellnum;i++)//0番目の細胞のJijにmutationを加える
    {
        child1 = randProtein(mt);
        child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
        child3 = randProtein(mt);
        child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
        
        while(rk[0].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
        {
            child1 = randProtein(mt);//inputは変化しない
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        while(rk[0].GetJ(child3,child4) != 0 ||
              rk[0].GetSEL(child3)>=MaxPathNum ||
              rk[0].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
        {
            child3 = randProtein(mt);//inputは変化しない
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
        }
        
        //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
        //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
        
        for(int l=0;l<input_number;l++)
        {
            if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
            {
                
                child3 = randProtein_non_T(mt);
                while(rk[0].GetSEL(child3)>=MaxPathNum ||
                      rk[0].GetJ(child3,child4) != 0 )//別の0のpathを探す
                {
                    child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                }
                
                break;
            }
        }
        
        for(int t=0;t<target_number;t++)
        {
            if(child3 == t)//繋ぎ変える先がtargetの時
            {
                count1=0;count2=0;
                for(int j=0;j<target_number;j++)
                {
                    if(rk[0].GetJ(j,child4) == 1)
                        //すでに2からtargetがactivateされている時
                    {
                        count1 = count1 + 1;
                    }
                    else if(rk[0].GetJ(j,child4) == -1)
                    {
                        count2 = count2 +1;
                    }
                }
                if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                {
                    child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                    while(rk[0].GetSEL(child3)>=MaxPathNum ||
                          rk[0].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                    }
                }
                break;
            }
        }
        
        j12 = rk[0].GetJ(child1,child2);//-1 or 1
        j34 = rk[0].GetJ(child3,child4);//0
        
        if(j12 != j34)
        {
            //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[i].GetJ(child1, child2)<<endl;
            //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[i].GetJ(child3, child4)<<endl;
            
            rk[i].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
            rk[i].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
            
            //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[i].GetJ(child1, child2)<<endl;
            //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[i].GetJ(child3, child4)<<endl;
            rk[i].converterJ();//隣接リストにmutationを反映
            
        }
    }
    SetFlag_of_VipVg(1);
}

void CELL::Vg_cell_make_random()
{
    int i=0,j=0,totalcellnum=0;
    int child1=0,child2=0,child3=0,child4=0,count1=0,count2=0,j12=0,j34=0,OneGlpNum;
    totalcellnum = rk[0].GetTotalcellnumber();
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProtein_non_T(target_number,ProteinVAL-1);//target以外
    
    totalcellnum = rk[0].GetTotalcellnumber();
    OneGlpNum = totalcellnum/Glp_Num;
    
    for(i=0;i<Glp_Num;i++)
    {
        for(j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
        {
            child1 = randProtein(mt);
            child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがある方
            child3 = randProtein(mt);
            child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止 pathがない方
            
            while(rk[i*OneGlpNum].GetJ(child1,child2) == 0 )//-1,+1　のpathを探す
            {
                child1 = randProtein(mt);//inputは変化しない
                child2 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            while(rk[i*OneGlpNum].GetJ(child3,child4) != 0 ||
                  rk[i*OneGlpNum].GetSEL(child3)>=MaxPathNum ||
                  rk[i*OneGlpNum].GetRSEL(child4)>=MaxPathNum)//0のpathを探す。ただし総本数がMaxPathNum本を超えないように
            {
                child3 = randProtein(mt);//inputは変化しない
                child4 = randProtein_non_T(mt);//TargetのFBや自己複製を禁止
            }
            
            //cout<<"SEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetSEL(child3)<<endl;
            //cout<<"RSEL="<<rk[data_array_fitness[eliminate_cellnumber].num].GetRSEL(child4)<<endl;
            
            for(int l=0;l<input_number;l++)
            {
                if(child4 == ProteinVAL-1-l)//input-geneが選ばれた時
                {
                    child3 = randProtein_non_T(mt);
                    while(rk[0].GetSEL(child3)>=MaxPathNum ||
                          rk[0].GetJ(child3,child4) != 0 )//別の0のpathを探す
                    {
                        child3 = randProtein_non_T(mt);//inputは変化しない,input-targetの直結は禁止
                    }
                    
                    break;
                }
            }
            
            for(int t=0;t<target_number;t++)
            {
                if(child3 == t)//繋ぎ変える先がtargetの時
                {
                    count1=0;count2=0;
                    for(int j=0;j<target_number;j++)
                    {
                        if(rk[i*OneGlpNum].GetJ(j,child4) == 1)
                            //すでに2からtargetがactivateされている時
                        {
                            count1 = count1 + 1;
                        }
                        else if(rk[i*OneGlpNum].GetJ(j,child4) == -1)
                        {
                            count2 = count2 +1;
                        }
                    }
                    if(count1>=MaxTargetPath||count2>=MaxTargetPath)
                    {
                        child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                        while(rk[i*OneGlpNum].GetSEL(child3)>=MaxPathNum ||
                              rk[i*OneGlpNum].GetJ(child3,child4) != 0 )//別の0のpathを探す
                        {
                            child3 = randProtein_non_T(mt);//targetに繋ぎすぎるのは禁止
                        }
                    }
                    break;
                }
            }
            
            j12 = rk[i*OneGlpNum].GetJ(child1,child2);//-1 or 1
            j34 = rk[i*OneGlpNum].GetJ(child3,child4);//0
            
            if(j12 != j34)
            {
                //cout<<"oldJ["<<child1<<"]["<<child2<<"]="<<rk[i].GetJ(child1, child2)<<endl;
                //cout<<"oldJ["<<child3<<"]["<<child4<<"]="<<rk[i].GetJ(child3, child4)<<endl;
                
                rk[j].SetJ(j34,child1,child2);//path切り替え 2--->1を2--->3に
                rk[j].SetJ(j12,child3,child4);//path切り替え 2--->3を2--->1に
                
                //cout<<"newJ["<<child1<<"]["<<child2<<"]="<<rk[i].GetJ(child1, child2)<<endl;
                //cout<<"newJ["<<child3<<"]["<<child4<<"]="<<rk[i].GetJ(child3, child4)<<endl;
                rk[j].converterJ();//隣接リストにmutationを反映
                
            }
        }
    }
}

void CELL::Load_Jij()
{
    int Totalcellnum=0,l,t,gene50=0;
    char filename[50]={};
    ifstream ifs1;
    
    vector<vector<int> >  GRN_2jigen;
    
    GRN_2jigen.resize(ProteinVAL);
    for(int i=0;i<ProteinVAL;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            GRN_2jigen[i].push_back(0);
        }
    }
    
    Totalcellnum = rk[0].GetTotalcellnumber();
    gene50 =GetGeneration()*500;//1世代時間の500秒をかける
    
    //sprintf(filename,"GRN_gene%d_001.txt",gene50);
    sprintf(filename,"First_GRN_gene%d.txt",gene50);
    //sprintf(filename,"Second_GRN_gene%d.txt",gene50);
    
    for(int i=0;i<Totalcellnum;i++)
    {
        //sprintf(filename,"GRN_gene%d_%03d.txt",gene50,i+1);//First,Secondを作るとき
        //cout<<"読み取り:"<<filename<<endl;
        ifs1.open(filename);
        
        if (!ifs1)
        {
            cout << "ファイルオープンに失敗" <<endl;
            cout<<"読み取り:"<<filename<<endl;
            exit(1);
        }
        else{
            for(l=0;l<ProteinVAL;l++)
            {
                for(t=0;t<ProteinVAL;t++)
                {
                    ifs1 >>GRN_2jigen[t][l];//"GRN%d.txt"　を読み込む
                    rk[i].SetJ(GRN_2jigen[t][l],t,l);//初期のネットワーク構造 GRN form t to l
                }
            }
            rk[i].converterJ();//隣接リストを作成
            ifs1.close();
        }
    }
}

void CELL::Load_Jji_random()//複数の遺伝子型を選んでcopyを作成
{
    int Totalcellnum=0,l,t,gene50=0,rn = 0,OneGlpNum=0;
    char filename[50]={};
    
    ifstream ifs1;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> randCELL(1,MAXFILENUMBER);
    vector<vector<int> >  GRN_2jigen;
    
    GRN_2jigen.resize(ProteinVAL);
    for(int i=0;i<ProteinVAL;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            GRN_2jigen[i].push_back(0);
        }
    }
    
    Totalcellnum = rk[0].GetTotalcellnumber();
    OneGlpNum = Totalcellnum/Glp_Num;
    gene50 =GetGeneration()*500;//1世代時間の500秒をかける
    
    for(int i=0;i<Glp_Num;i++)
    {
        rn =randCELL(mt);
        if(GetGeneration()==0)
        {
            sprintf(filename,"J_3jigen1.txt");
        }
        else{
            sprintf(filename,"GRN_gene%d_%03d.txt",gene50,rn);
        }
        ifs1.open(filename);
        if (!ifs1)
        {
            cout << "ファイルオープンに失敗" <<endl;
            cout<<"読み取り:"<<filename<<endl;
            exit(1);
        }
        else{
            cout<<"読み取り:"<<filename<<endl;
            for(int j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
            {
                for(l=0;l<ProteinVAL;l++)
                {
                    for(t=0;t<ProteinVAL;t++)
                    {
                        ifs1 >>GRN_2jigen[t][l];//"GRN%d.txt"　を読み込む
                        rk_VipVg[j].SetJ(GRN_2jigen[t][l],t,l);//初期のネットワーク構造 GRN form t to l
                    }
                }
                rk_VipVg[j].converterJ();//隣接リストを作成
                
            }
        }
        ifs1.close();
    }
}

void CELL::ReLoad_Jji()//複数の遺伝子型を選んでcopyを作成
{
    int Totalcellnum=0,l=0,t=0,OneGlpNum=0;
    
    ifstream ifs1;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> randCELL(1,MAXCELLNUMBER);
    vector<vector<int> >  GRN_2jigen;
    
    Totalcellnum = rk[0].GetTotalcellnumber();
    OneGlpNum = Totalcellnum/Glp_Num;
    
    for(int i=0;i<Glp_Num;i++)
    {
        for(int j=i*OneGlpNum;j<(i+1)*OneGlpNum;j++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                for(t=0;t<ProteinVAL;t++)
                {
                    rk[j].SetJ(rk[i*OneGlpNum].GetJ(t, l),t,l);//初期のネットワーク構造 GRN form t to l
                }
            }
            rk[j].converterJ();//隣接リストを作成
        }
    }
    SetFlag_of_VipVg(1);
}
void Load_Jij_typical();//集団の典型的な遺伝子型

void CELL::SetGeneration(int gene)
{
    if(gene<0)
    {
        cout<<"gene is "<<gene;
        exit(-1);
    }
    else{
        generation = gene;
    }
}
int  CELL::GetGeneration()
{
    return generation;
}

void CELL::SetFlag_of_VipVg(int flag)
{
    if(flag==0||flag==1||flag==2||flag==3)
    {
        flag_of_VipVg = flag;
    }
    else{
        cout<<"flag_of_VipVg ="<<flag<<endl;
        exit(-1);
    }
    
}
int CELL::Getflag_of_VipVg()
{
    return flag_of_VipVg;
}

void CELL::archive_var_Vip(int time)
{
    FILE *pf = nullptr;
    FILE *pg = nullptr;
    FILE *pr = nullptr;
    
    char filename[50]={};
    char filenameg[50]={};
    char filenamer[50]={};
    double sum=0.0;
    
    sprintf(filenameg,"var0001_sta.txt");
    sprintf(filename,"var0001_First_Vip.txt");
    //sprintf(filename,"var0001_sta_Second.txt");
    //sprintf(filename,"var0001_sta_First_Vg.txt");
    pf = fopen(filename,"a+");
    //sprintf(filenameg,"genome_var0001.txt");
    pg = fopen(filenameg,"a+");
    //sprintf(filenamer,"fitness_var0001.txt");
    sprintf(filenamer,"fitness_var0001_First_Vip.txt");
    //sprintf(filenamer,"fitness_var0001_First_Vg.txt");
    //sprintf(filenamer,"fitness_var0001_Second.txt");
    pr = fopen(filenamer,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pf,"%.8f ",GetVariance(i));//"var0001_First_Vip.txt"
            sum = sum + GetVariance(i);
        }
        sum =(double)(sum/ProteinVAL);
        fprintf(pf,"%.8f ",sum);//分散の平均値
        fprintf(pf,"\n");
        fprintf(pg," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pg,"%.8f ",GetVariance_Sta(i));
            sum = sum + GetVariance_Sta(i);
        }
        sum =(double)(sum/ProteinVAL);
        fprintf(pg,"%.8f ",sum);//分散の平均値
        fprintf(pg,"\n");
        fprintf(pr," %d ",GetGeneration());
        fprintf(pr,"%f ",GetMean_fitness());//平均fitness
        fprintf(pr,"%.8f ",GetFitness_variance());
        fprintf(pr,"\n");
    }
    fclose(pf);
    fclose(pg);
    fclose(pr);
}

void CELL::archive_var_Vip_average(int time)
{
    FILE *pf = nullptr;
    char filename[50]={};
    double sum=0.0;
    
    sprintf(filename,"var0001_First_Vip_average.txt");
    pf = fopen(filename,"a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",time);
        for(int i=0;i<ProteinVAL;i++)
        {
            fprintf(pf,"%.8f ",GetVarianceAverage(i));
            sum = sum + GetVarianceAverage(i);
        }
        sum =(double)(sum/ProteinVAL);
        fprintf(pf,"%.8f ",sum);//分散の平均値
        fprintf(pf,"\n");
    }
    fclose(pf);
    
}
//Parasite用
double  CELL::GetParasite_population(int ith) const//ithはparasiteの種番号
{
    return parasite_population[ith];
}
void CELL::SetParasite_population(int ith,double pp)
{
    parasite_population[ith] = pp;
}
double CELL::GetParasite_Genotype(int ith,int jth) const
{
    return parasite_genotype[ith][jth];
}
void CELL::SetParasite_Genotype(int ith,int jth,double pg)//ithはparasiteの種番号
{
    parasite_genotype[ith][jth] = pg;
}
double CELL::GetHosto_num(int ith)//ithはparasiteの種番号
{
    return hosto_num[ith];
}
void CELL::SetHosto_num(int ith,int hn)
{
    hosto_num[ith] = hn;
}
double CELL::GetParasite_virulence(int ith) const//ithはparasiteの種番号
{
    return parasite_virulence[ith];
}
void CELL::SetParasite_virulence(int ith,double vir)
{
    if(vir>0.000001)
    {
        parasite_virulence[ith] = vir;
    }
    else{
        parasite_virulence[ith] = -vir;
    }
}

void CELL::calculate_parasite_dynamics()//ここで少しpopulationの計算を行う
{
    int i,dynamics_cal_time=0;
    dynamics_cal_time = (int)(parasite_population_dynamics_cal_times/parasite_dt);//刻み時間はhostよりも細かくとる
    
    for(i=0;i<dynamics_cal_time;i++)
    {
        population_dynamics();
    }
    
}

void CELL::calculate_parasite_dynamics_evo()//ここで少しpopulationの計算を行う
{
    int i,dynamics_cal_time=0;
    double rara1,delta=0.01;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> rand01(0, 1);
    
    rara1 = rand01(mt);
    
    dynamics_cal_time = (int)(parasite_population_dynamics_cal_times/parasite_dt);//刻み時間はhostよりも細かくとる
    
    for(i=0;i<dynamics_cal_time;i++)
    {
        population_dynamics();
    }
    
    i=0;
    while(GetParasite_population(i)>rara1)
    {
        SetParasite_virulence(i,rk[0].gabs_gauss_rand(GetParasite_virulence(i),delta));
        i++;
    }
}

void CELL::calculate_parasite_dynamics_dist_fitness()
{
    int i,dynamics_cal_time=0;
    dynamics_cal_time = (int)(parasite_population_dynamics_cal_times/parasite_dt);//刻み時間はhostよりも細かくとる
    
    for(i=0;i<dynamics_cal_time;i++)
    {
        population_dynamics_dist_fitness();
    }
    calculate_parasite_dist_fitness();
}

void CELL::calculate_parasite_dist_fitness()
{
    double r_in=0.0,c=0.0,ramda=10.0;
    int cell_num=0;
    
    ramda =LAMBDA;
    cell_num = rk[0].GetTotalcellnumber();
    c=maximum_virulence;
    
    for(int j=0;j<cell_num;j++)
    {
        for(int i=0;i<Parasite_species_num;i++)
        {
            r_in=0.0;
            for(int ith=0;ith<Parasite_genome_size;ith++)
            {
                r_in = r_in+(rk[j].GetResult(Parasite_target_genome+ith)-GetParasite_Genotype(i,ith))*(rk[j].GetResult(Parasite_target_genome+ith)-GetParasite_Genotype(i,ith));
            }
            //cout <<"dist="<<r_in<<endl;
            //cout <<"fitness="<<exp(-ramda*r_in)<<endl;
            //cout <<"parasite_density="<<rk[j].GetParasite_population_rate(i)<<endl;
            rk[j].SetInfection_virulence(rk[j].GetParasite_population_rate(i)*c*exp(-ramda*r_in),i);
        }
    }
    for(int i=0;i<Parasite_species_num;i++)
    {
        cout<<"Infection dist["<<i<<"]="<<rk[0].GetInfection_virulence(i)<<endl;
    }
}

void CELL::population_dynamics()
{
    int i=0,j=0,totalcell=0;
    double sum=0.0,h=0.0;
    
    double ParaPopu[Parasite_species_num];
    double k1[Parasite_species_num];
    double k2[Parasite_species_num];
    double k3[Parasite_species_num];
    double k4[Parasite_species_num];
    double tmp[Parasite_species_num];
    
    h=parasite_dt;
    totalcell = rk[0].GetTotalcellnumber();
    
    for(i=0;i<Parasite_species_num;i++)
    {
        ParaPopu[i] =0.0;
        k1[i] = 0.0;
        k2[i] = 0.0;
        k3[i] = 0.0;
        k4[i] = 0.0;
        tmp[i] = 0.0;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        ParaPopu[i] = GetParasite_population(i);
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        k1[i] = func_population_dynamics(ParaPopu,i);
        
        tmp[i] = ParaPopu[i] + h*k1[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        k2[i]= func_population_dynamics(tmp,i);
        
        tmp[i] = ParaPopu[i] + h*k2[i]/2;
        //cout<<"k2["<<i<<"]="<<k2[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<Parasite_species_num;i++)
    {
        k3[i] = func_population_dynamics(tmp,i);
        
        tmp[i] = ParaPopu[i] + h*k3[i];
        //cout<<"k3["<<i<<"]="<<k3[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        k4[i] = func_population_dynamics(tmp,i);
        
        ParaPopu[i] = ParaPopu[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        //cout<<"k4["<<i<<"]="<<k4[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
        sum = sum + ParaPopu[i];
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //SetParasite_population(i,ParaPopu[i]);
        SetParasite_population(i,ParaPopu[i]/sum);
    }
    
    for(i=0;i<totalcell;i++)
    {
        rk[i].SetParasite_population(GetParasite_population(rk[i].GetParasite_attacked_number()));
        for(j=0;j<Parasite_species_num;j++)
        {
            rk[i].SetParasite_population_rate(GetParasite_population(j),j);
        }
    }
    
}
void CELL::population_dynamics_dist_fitness()
{
    int i=0,j=0,totalcell=0;
    double sum=0.0,h=0.0;
    
    double ParaPopu[Parasite_species_num];
    double k1[Parasite_species_num];
    double k2[Parasite_species_num];
    double k3[Parasite_species_num];
    double k4[Parasite_species_num];
    double tmp[Parasite_species_num];
    
    h=parasite_dt;
    totalcell = rk[0].GetTotalcellnumber();
    
    for(i=0;i<Parasite_species_num;i++)
    {
        ParaPopu[i] =0.0;
        k1[i] = 0.0;
        k2[i] = 0.0;
        k3[i] = 0.0;
        k4[i] = 0.0;
        tmp[i] = 0.0;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        ParaPopu[i] = GetParasite_population(i);
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //k1[i] = func_population_dynamics(ParaPopu,i);
        k1[i] = func_population_dynamics_dist_fitness(ParaPopu,i);
        tmp[i] = ParaPopu[i] + h*k1[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //k2[i]= func_population_dynamics(tmp,i);
        k2[i] = func_population_dynamics_dist_fitness(tmp,i);
        tmp[i] = ParaPopu[i] + h*k2[i]/2;
        //cout<<"k2["<<i<<"]="<<k2[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<Parasite_species_num;i++)
    {
        //k3[i] = func_population_dynamics(tmp,i);
        k3[i] = func_population_dynamics_dist_fitness(tmp,i);
        tmp[i] = ParaPopu[i] + h*k3[i];
        //cout<<"k3["<<i<<"]="<<k3[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //k4[i] = func_population_dynamics(tmp,i);
        k4[i] = func_population_dynamics_dist_fitness(tmp,i);
        ParaPopu[i] = ParaPopu[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        //cout<<"k4["<<i<<"]="<<k4[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<Parasite_species_num;i++)
    {
        //SetParasite_population(i,ParaPopu[i]);
        SetParasite_population(i,ParaPopu[i]);
    }
    
    for(i=0;i<totalcell;i++)
    {
        rk[i].SetParasite_population(GetParasite_population(rk[i].GetParasite_attacked_number()));
        for(j=0;j<Parasite_species_num;j++)
        {
            rk[i].SetParasite_population_rate(GetParasite_population(j),j);
        }
    }
    
}

double CELL::func_population_dynamics(double *y,int ith)
{
    double D_P=0.001,r=0.0,r_out=0.0,r_Yin=0.0,c=0.0;
    
    //c = maximum_virulence;
    c = GetParasite_virulence(ith);
    //r_Yin = GetHosto_num(ith)/(double)(MAXCELLNUMBER)*y[ith];//生成項
    r_Yin = c*GetHosto_num(ith)*y[ith];//生成項
    
    if(Parasite_genome_size==3)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4]);
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7]);
                break;
        }
    }
    else if(Parasite_genome_size==1)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1]);
                break;
        }
    }
    else if(Parasite_genome_size==2)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3]);
                break;
        }
    }
    else if(Parasite_genome_size==4)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0])+D_P*(y[8]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1])+D_P*(y[9]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2])+D_P*(y[10]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3])+D_P*(y[14]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4])+D_P*(y[11]-y[4]);
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5])+D_P*(y[12]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6])+D_P*(y[13]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7])+D_P*(y[15]-y[7]);
                break;
            case 8:
                r_out = D_P*(y[0]-y[8])+D_P*(y[9]-y[8])+D_P*(y[10]-y[8])+D_P*(y[11]-y[8]);
                break;
            case 9:
                r_out = D_P*(y[1]-y[9])+D_P*(y[8]-y[9])+D_P*(y[12]-y[9])+D_P*(y[14]-y[9]);
                break;
            case 10:
                r_out = D_P*(y[2]-y[10])+D_P*(y[8]-y[10])+D_P*(y[13]-y[10])+D_P*(y[14]-y[10]);
                break;
            case 11:
                r_out = D_P*(y[4]-y[11])+D_P*(y[8]-y[11])+D_P*(y[12]-y[11])+D_P*(y[13]-y[11]);
                break;
            case 12:
                r_out = D_P*(y[5]-y[12])+D_P*(y[9]-y[12])+D_P*(y[11]-y[12])+D_P*(y[15]-y[12]);
                break;
            case 13:
                r_out = D_P*(y[6]-y[13])+D_P*(y[10]-y[13])+D_P*(y[11]-y[13])+D_P*(y[15]-y[13]);
                break;
            case 14:
                r_out = D_P*(y[3]-y[14])+D_P*(y[9]-y[14])+D_P*(y[10]-y[14])+D_P*(y[15]-y[14]);
                break;
            case 15:
                r_out = D_P*(y[7]-y[15])+D_P*(y[12]-y[15])+D_P*(y[13]-y[15])+D_P*(y[14]-y[15]);
                break;
        }
    }
    
    else if(Parasite_genome_size==5)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0])+D_P*(y[8]-y[0])+D_P*(y[16]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1])+D_P*(y[9]-y[1])+D_P*(y[17]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2])+D_P*(y[10]-y[2])+D_P*(y[18]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3])+D_P*(y[11]-y[3])+D_P*(y[19]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4])+D_P*(y[12]-y[4]+D_P*(y[20]-y[4]));
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5])+D_P*(y[13]-y[5])+D_P*(y[21]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6])+D_P*(y[14]-y[6])+D_P*(y[22]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7])+D_P*(y[15]-y[7])+D_P*(y[23]-y[7]);
                break;
            case 8:
                r_out = D_P*(y[0]-y[8])+D_P*(y[9]-y[8])+D_P*(y[10]-y[8])+D_P*(y[12]-y[8])+D_P*(y[24]-y[8]);
                break;
            case 9:
                r_out = D_P*(y[1]-y[9])+D_P*(y[8]-y[9])+D_P*(y[11]-y[9])+D_P*(y[13]-y[9])+D_P*(y[25]-y[9]);
                break;
            case 10:
                r_out = D_P*(y[2]-y[10])+D_P*(y[8]-y[10])+D_P*(y[11]-y[10])+D_P*(y[14]-y[10])+D_P*(y[26]-y[10]);
                break;
            case 11:
                r_out = D_P*(y[3]-y[11])+D_P*(y[9]-y[11])+D_P*(y[10]-y[11])+D_P*(y[15]-y[11])+D_P*(y[27]-y[11]);
                break;
            case 12:
                r_out = D_P*(y[4]-y[12])+D_P*(y[8]-y[12])+D_P*(y[13]-y[12])+D_P*(y[14]-y[12])+D_P*(y[28]-y[12]);
                break;
            case 13:
                r_out = D_P*(y[5]-y[13])+D_P*(y[9]-y[13])+D_P*(y[12]-y[13])+D_P*(y[15]-y[13])+D_P*(y[29]-y[13]);
                break;
            case 14:
                r_out = D_P*(y[6]-y[14])+D_P*(y[10]-y[14])+D_P*(y[12]-y[14])+D_P*(y[15]-y[14])+D_P*(y[30]-y[14]);
                break;
            case 15:
                r_out = D_P*(y[7]-y[15])+D_P*(y[11]-y[15])+D_P*(y[13]-y[15])+D_P*(y[14]-y[15])+D_P*(y[31]-y[15]);
                break;
            case 16:
                r_out = D_P*(y[0]-y[16])+D_P*(y[17]-y[16])+D_P*(y[18]-y[16])+D_P*(y[20]-y[16])+D_P*(y[24]-y[16]);
                break;
            case 17:
                r_out = D_P*(y[1]-y[17])+D_P*(y[16]-y[17])+D_P*(y[19]-y[17])+D_P*(y[21]-y[17])+D_P*(y[25]-y[17]);
                break;
            case 18:
                r_out = D_P*(y[2]-y[18])+D_P*(y[16]-y[18])+D_P*(y[19]-y[18])+D_P*(y[22]-y[18])+D_P*(y[26]-y[18]);
                break;
            case 19:
                r_out = D_P*(y[3]-y[19])+D_P*(y[17]-y[19])+D_P*(y[18]-y[19])+D_P*(y[23]-y[19])+D_P*(y[27]-y[19]);
                break;
            case 20:
                r_out = D_P*(y[4]-y[20])+D_P*(y[16]-y[20])+D_P*(y[21]-y[20])+D_P*(y[22]-y[20])+D_P*(y[28]-y[20]);
                break;
            case 21:
                r_out = D_P*(y[5]-y[21])+D_P*(y[17]-y[21])+D_P*(y[20]-y[21])+D_P*(y[23]-y[21])+D_P*(y[29]-y[21]);
                break;
            case 22:
                r_out = D_P*(y[6]-y[22])+D_P*(y[18]-y[22])+D_P*(y[20]-y[22])+D_P*(y[23]-y[22])+D_P*(y[30]-y[22]);
                break;
            case 23:
                r_out = D_P*(y[7]-y[23])+D_P*(y[19]-y[23])+D_P*(y[21]-y[23])+D_P*(y[22]-y[23])+D_P*(y[31]-y[23]);
                break;
            case 24:
                r_out = D_P*(y[8]-y[24])+D_P*(y[16]-y[24])+D_P*(y[25]-y[24])+D_P*(y[26]-y[24])+D_P*(y[28]-y[24]);
                break;
            case 25:
                r_out = D_P*(y[9]-y[25])+D_P*(y[17]-y[25])+D_P*(y[24]-y[25])+D_P*(y[27]-y[25])+D_P*(y[29]-y[25]);
                break;//ここまで
            case 26:
                r_out = D_P*(y[10]-y[26])+D_P*(y[18]-y[26])+D_P*(y[24]-y[26])+D_P*(y[27]-y[26])+D_P*(y[30]-y[26]);
                break;
            case 27:
                r_out = D_P*(y[11]-y[27])+D_P*(y[19]-y[27])+D_P*(y[25]-y[27])+D_P*(y[26]-y[27])+D_P*(y[31]-y[27]);
                break;
            case 28:
                r_out = D_P*(y[12]-y[28])+D_P*(y[20]-y[28])+D_P*(y[24]-y[28])+D_P*(y[29]-y[28])+D_P*(y[30]-y[28]);
                break;
            case 29:
                r_out = D_P*(y[13]-y[29])+D_P*(y[21]-y[29])+D_P*(y[25]-y[29])+D_P*(y[28]-y[29])+D_P*(y[31]-y[29]);
                break;
            case 30:
                r_out = D_P*(y[14]-y[30])+D_P*(y[22]-y[30])+D_P*(y[26]-y[30])+D_P*(y[28]-y[30])+D_P*(y[31]-y[30]);
                break;
            case 31:
                r_out = D_P*(y[15]-y[31])+D_P*(y[23]-y[31])+D_P*(y[27]-y[31])+D_P*(y[29]-y[31])+D_P*(y[30]-y[31]);
                break;
                
        }
    }
    
    //cout<<"r_out["<<ith<<"]="<<r_out<<endl;
    //cout<<"r_in["<<ith<<"]="<<r_Yin<<endl;
    r = r_out + r_Yin;
    
    return r;
}

double CELL::func_population_dynamics_dist_fitness(double *y,int ith)
{
    double D_P=0.01,r=0.0,r_out=0.0,sum=0.0,r_in=0.0,r_death=0.0,c=0.0,omega=0.0,ramda=10.0;
    int cell_num=0;
    c = maximum_virulence;
    omega = c*MAXCELLNUMBER*0.5;//ホストが100体を超えないと死ぬほうが早い
    cell_num =rk[0].GetTotalcellnumber();
    ramda=LAMBDA;
    
    for(int i=0;i<cell_num;i++)
    {
        r_in=0.0;
        for(int j=0;j<Parasite_genome_size;j++)
        {
            r_in = r_in+(rk[i].GetResult(Parasite_target_genome+j)-GetParasite_Genotype(ith,j))*(rk[i].GetResult(Parasite_target_genome+j)-GetParasite_Genotype(ith,j));//生成項
        }
        sum = sum +c*exp(-ramda*r_in);
    }
    
    r_in =sum*y[ith];
    r_death = -omega*y[ith];//生成項
    
    cout<<"r_death["<<ith<<"]="<<r_death<<endl;
    cout<<"r_in["<<ith<<"]="<<r_in<<endl;
    
    if(Parasite_genome_size==3)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4]);
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7]);
                break;
        }
    }
    else if(Parasite_genome_size==1)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1]);
                break;
        }
    }
    else if(Parasite_genome_size==2)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3]);
                break;
        }
    }
    else if(Parasite_genome_size==4)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0])+D_P*(y[8]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1])+D_P*(y[9]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2])+D_P*(y[10]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3])+D_P*(y[14]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4])+D_P*(y[11]-y[4]);
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5])+D_P*(y[12]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6])+D_P*(y[13]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7])+D_P*(y[15]-y[7]);
                break;
            case 8:
                r_out = D_P*(y[0]-y[8])+D_P*(y[9]-y[8])+D_P*(y[10]-y[8])+D_P*(y[11]-y[8]);
                break;
            case 9:
                r_out = D_P*(y[1]-y[9])+D_P*(y[8]-y[9])+D_P*(y[12]-y[9])+D_P*(y[14]-y[9]);
                break;
            case 10:
                r_out = D_P*(y[2]-y[10])+D_P*(y[8]-y[10])+D_P*(y[13]-y[10])+D_P*(y[14]-y[10]);
                break;
            case 11:
                r_out = D_P*(y[4]-y[11])+D_P*(y[8]-y[11])+D_P*(y[12]-y[11])+D_P*(y[13]-y[11]);
                break;
            case 12:
                r_out = D_P*(y[5]-y[12])+D_P*(y[9]-y[12])+D_P*(y[11]-y[12])+D_P*(y[15]-y[12]);
                break;
            case 13:
                r_out = D_P*(y[6]-y[13])+D_P*(y[10]-y[13])+D_P*(y[11]-y[13])+D_P*(y[15]-y[13]);
                break;
            case 14:
                r_out = D_P*(y[3]-y[14])+D_P*(y[9]-y[14])+D_P*(y[10]-y[14])+D_P*(y[15]-y[14]);
                break;
            case 15:
                r_out = D_P*(y[7]-y[15])+D_P*(y[12]-y[15])+D_P*(y[13]-y[15])+D_P*(y[14]-y[15]);
                break;
        }
    }
    
    else if(Parasite_genome_size==5)
    {
        switch(ith)
        {
            case 0:
                r_out = D_P*(y[1]-y[0])+D_P*(y[2]-y[0])+D_P*(y[4]-y[0])+D_P*(y[8]-y[0])+D_P*(y[16]-y[0]);
                break;
            case 1:
                r_out = D_P*(y[0]-y[1])+D_P*(y[3]-y[1])+D_P*(y[5]-y[1])+D_P*(y[9]-y[1])+D_P*(y[17]-y[1]);
                break;
            case 2:
                r_out = D_P*(y[0]-y[2])+D_P*(y[3]-y[2])+D_P*(y[6]-y[2])+D_P*(y[10]-y[2])+D_P*(y[18]-y[2]);
                break;
            case 3:
                r_out = D_P*(y[1]-y[3])+D_P*(y[2]-y[3])+D_P*(y[7]-y[3])+D_P*(y[11]-y[3])+D_P*(y[19]-y[3]);
                break;
            case 4:
                r_out = D_P*(y[0]-y[4])+D_P*(y[5]-y[4])+D_P*(y[6]-y[4])+D_P*(y[12]-y[4]+D_P*(y[20]-y[4]));
                break;
            case 5:
                r_out = D_P*(y[1]-y[5])+D_P*(y[4]-y[5])+D_P*(y[7]-y[5])+D_P*(y[13]-y[5])+D_P*(y[21]-y[5]);
                break;
            case 6:
                r_out = D_P*(y[2]-y[6])+D_P*(y[4]-y[6])+D_P*(y[7]-y[6])+D_P*(y[14]-y[6])+D_P*(y[22]-y[6]);
                break;
            case 7:
                r_out = D_P*(y[3]-y[7])+D_P*(y[5]-y[7])+D_P*(y[6]-y[7])+D_P*(y[15]-y[7])+D_P*(y[23]-y[7]);
                break;
            case 8:
                r_out = D_P*(y[0]-y[8])+D_P*(y[9]-y[8])+D_P*(y[10]-y[8])+D_P*(y[12]-y[8])+D_P*(y[24]-y[8]);
                break;
            case 9:
                r_out = D_P*(y[1]-y[9])+D_P*(y[8]-y[9])+D_P*(y[11]-y[9])+D_P*(y[13]-y[9])+D_P*(y[25]-y[9]);
                break;
            case 10:
                r_out = D_P*(y[2]-y[10])+D_P*(y[8]-y[10])+D_P*(y[11]-y[10])+D_P*(y[14]-y[10])+D_P*(y[26]-y[10]);
                break;
            case 11:
                r_out = D_P*(y[3]-y[11])+D_P*(y[9]-y[11])+D_P*(y[10]-y[11])+D_P*(y[15]-y[11])+D_P*(y[27]-y[11]);
                break;
            case 12:
                r_out = D_P*(y[4]-y[12])+D_P*(y[8]-y[12])+D_P*(y[13]-y[12])+D_P*(y[14]-y[12])+D_P*(y[28]-y[12]);
                break;
            case 13:
                r_out = D_P*(y[5]-y[13])+D_P*(y[9]-y[13])+D_P*(y[12]-y[13])+D_P*(y[15]-y[13])+D_P*(y[29]-y[13]);
                break;
            case 14:
                r_out = D_P*(y[6]-y[14])+D_P*(y[10]-y[14])+D_P*(y[12]-y[14])+D_P*(y[15]-y[14])+D_P*(y[30]-y[14]);
                break;
            case 15:
                r_out = D_P*(y[7]-y[15])+D_P*(y[11]-y[15])+D_P*(y[13]-y[15])+D_P*(y[14]-y[15])+D_P*(y[31]-y[15]);
                break;
            case 16:
                r_out = D_P*(y[0]-y[16])+D_P*(y[17]-y[16])+D_P*(y[18]-y[16])+D_P*(y[20]-y[16])+D_P*(y[24]-y[16]);
                break;
            case 17:
                r_out = D_P*(y[1]-y[17])+D_P*(y[16]-y[17])+D_P*(y[19]-y[17])+D_P*(y[21]-y[17])+D_P*(y[25]-y[17]);
                break;
            case 18:
                r_out = D_P*(y[2]-y[18])+D_P*(y[16]-y[18])+D_P*(y[19]-y[18])+D_P*(y[22]-y[18])+D_P*(y[26]-y[18]);
                break;
            case 19:
                r_out = D_P*(y[3]-y[19])+D_P*(y[17]-y[19])+D_P*(y[18]-y[19])+D_P*(y[23]-y[19])+D_P*(y[27]-y[19]);
                break;
            case 20:
                r_out = D_P*(y[4]-y[20])+D_P*(y[16]-y[20])+D_P*(y[21]-y[20])+D_P*(y[22]-y[20])+D_P*(y[28]-y[20]);
                break;
            case 21:
                r_out = D_P*(y[5]-y[21])+D_P*(y[17]-y[21])+D_P*(y[20]-y[21])+D_P*(y[23]-y[21])+D_P*(y[29]-y[21]);
                break;
            case 22:
                r_out = D_P*(y[6]-y[22])+D_P*(y[18]-y[22])+D_P*(y[20]-y[22])+D_P*(y[23]-y[22])+D_P*(y[30]-y[22]);
                break;
            case 23:
                r_out = D_P*(y[7]-y[23])+D_P*(y[19]-y[23])+D_P*(y[21]-y[23])+D_P*(y[22]-y[23])+D_P*(y[31]-y[23]);
                break;
            case 24:
                r_out = D_P*(y[8]-y[24])+D_P*(y[16]-y[24])+D_P*(y[25]-y[24])+D_P*(y[26]-y[24])+D_P*(y[28]-y[24]);
                break;
            case 25:
                r_out = D_P*(y[9]-y[25])+D_P*(y[17]-y[25])+D_P*(y[24]-y[25])+D_P*(y[27]-y[25])+D_P*(y[29]-y[25]);
                break;//ここまで
            case 26:
                r_out = D_P*(y[10]-y[26])+D_P*(y[18]-y[26])+D_P*(y[24]-y[26])+D_P*(y[27]-y[26])+D_P*(y[30]-y[26]);
                break;
            case 27:
                r_out = D_P*(y[11]-y[27])+D_P*(y[19]-y[27])+D_P*(y[25]-y[27])+D_P*(y[26]-y[27])+D_P*(y[31]-y[27]);
                break;
            case 28:
                r_out = D_P*(y[12]-y[28])+D_P*(y[20]-y[28])+D_P*(y[24]-y[28])+D_P*(y[29]-y[28])+D_P*(y[30]-y[28]);
                break;
            case 29:
                r_out = D_P*(y[13]-y[29])+D_P*(y[21]-y[29])+D_P*(y[25]-y[29])+D_P*(y[28]-y[29])+D_P*(y[31]-y[29]);
                break;
            case 30:
                r_out = D_P*(y[14]-y[30])+D_P*(y[22]-y[30])+D_P*(y[26]-y[30])+D_P*(y[28]-y[30])+D_P*(y[31]-y[30]);
                break;
            case 31:
                r_out = D_P*(y[15]-y[31])+D_P*(y[23]-y[31])+D_P*(y[27]-y[31])+D_P*(y[29]-y[31])+D_P*(y[30]-y[31]);
                break;
                
        }
    }
    
    cout<<"r_out["<<ith<<"]="<<r_out<<endl;
    r = r_death+r_out + r_in;
    
    return r;
}

int binToUInt(const std::string &str)
{
    unsigned int val = 0;
    for(int i = 0; i < (int)str.size(); ++i) {
        switch (str[i]) {
            case '0':
                val *= 2;
                break;
            case '1':
                val = val * 2 + 1;
                break;
        }
    }
    return val;
}

void CELL::caluculate_hosto_num()
{
    double r=0.0;
    int totalcellnum=0,i,j,type=0;
    int host_hako[Parasite_species_num]={};
    string b(Parasite_genome_size,'0');
    totalcellnum=rk[0].GetTotalcellnumber();
    
    for(j=0;j<Parasite_species_num;j++)
    {
        host_hako[j] =0;
    }
    
    //ここから各cellのtype決め
    for(i=0;i<totalcellnum;i++)
    {
        for(j=0;j<Parasite_genome_size;j++)
        {
            r = rk[i].GetResult(j+Parasite_target_genome);
            if(r>=0.5){
                b.push_back('1');
            }
            else{
                b.push_back('0');
            }
        }
        //cout<<"b="<<b<<endl;
        type = binToUInt(b);//2進数を10進数に変換して種のlabelにする
        b.clear();
        rk[i].SetParasite_attacked_number(type);//type:0~2**Parasite_genome_size-1
        host_hako[type] = host_hako[type]+1;//host_numのカウントをここで行う
        //cout<<"host_hako["<<type<<"]="<<host_hako[type]<<endl;
        //cout<<"rk["<<i<<"].Parasite_attacked_number="<<rk[i].GetParasite_attacked_number()<<endl;
    }
    
    for(j=0;j<Parasite_species_num;j++)
    {
        //cout<<"Host_num["<<j<<"]="<<host_hako[j]<<endl;
    }
    for(j=0;j<Parasite_genome_size;j++)
    {
        //cout<<"attacked_gene["<<j<<"]="<<rk[0].GetResult(j+Parasite_target_genome)<<endl;
    }
    //cout<<"Parasite_attacked_number="<<rk[0].GetParasite_attacked_number()<<endl;
    for(j=0;j<Parasite_species_num;j++)
    {
        SetHosto_num(j,host_hako[j]);
        //cout<<"GetHost_num["<<j<<"]="<<GetHosto_num(j)<<endl;
    }
}

void CELL::archive_parasite_population(int time)
{
    FILE *pf = nullptr;
    FILE *pg = nullptr;
    
    int l;
    pf = fopen("para_population.txt","a+");
    pg = fopen("host_num.txt","a+");
    
    if(pf==NULL||pg==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",time);
        fprintf(pf," %d ",rk[0].GetDivision_number());
        for(l=0;l<Parasite_species_num;l++)
        {
            fprintf(pf," %f ",GetParasite_population(l));
        }
        fprintf(pf,"\n");
        fprintf(pg," %d ",time);
        fprintf(pg," %d ",rk[0].GetDivision_number());
        for(l=0;l<Parasite_species_num;l++)
        {
            fprintf(pg," %f ",GetHosto_num(l));
        }
        fprintf(pg,"\n");
    }
    
    fclose(pf);
    fclose(pg);
}

void CELL::SetParasite_list(int list,int ith,int jth)
{
    parasite_list[ith][jth] = list;
}

int CELL::GetParasite_list(int ith,int jth)
{
    return parasite_list[ith][jth];
}

void CELL::parasite_list_maker()
{
    
}

void CELL::Initial_Jij_maker()
{
    
}


//以下Meduim Class
Medium::Medium(double t0,int div,double tn0,double medium_D,int _cell_number,double Volume,double S):time(t0),startTime(t0),LOOP(div),finishTime(tn0),Medium_D(medium_D),cell_number(_cell_number),V(Volume),S_out(S)
{
    Y.clear();
    for (int i=0;i<input_number; i++)
    {
        Y.push_back(S_resource);//resorverと濃度を合わせておく
    }
    X.clear();
    for (int i=0;i<ProteinVAL; i++)
    {
        X.push_back(0.0);//初期値はゼロ
    }
    delt = (double)((tn0-t0)/div);
}

void Medium::SetMedium_D(double medium_D)
{
    if(medium_D<0||medium_D>100)
    {
        cout<<"medium_D is inappropriate"<<endl;
        cout<<"medium_D = "<<medium_D<<endl;
    }
    else
    {
        Medium_D = medium_D;
    }
}
double Medium::GetMedium_D() const
{
    return Medium_D;
}

void Medium::SetY(double Y_med,int ith)
{
    if(Y_med<0||Y_med>100)
    {
        cout<<"medium_Y["<<ith<<"]="<<Y_med<<endl;
        exit(1);
    }
    else
    {
        Y[ith] = Y_med;
    }
}

double Medium::GetY(int ith) const
{
    return Y[ith];
}

void Medium::SetX(double X_med,int ith)
{
    if(X_med<0||X_med>100)
    {
        cout<<"medium_X["<<ith<<"]="<<X_med<<endl;
        exit(1);
    }
    else
    {
        X[ith] = X_med;
    }
}

double Medium::GetX(int ith) const
{
    return X[ith];
}

void Medium::calculateNS_med(RK *rk)
{
    double *k1,*k2,*k3,*k4,*tmp;
    double t=0,h=0;
    double *Y_medium;
    int i;
    
    Y_medium = new double[input_number];
    k1 = new double[input_number];
    k2 = new double[input_number];
    k3 = new double[input_number];
    k4 = new double[input_number];
    tmp = new double[input_number];
    
    t=rk[0].GetTime();
    h=rk[0].GetDelt();
    
    for(i=0;i<input_number;i++)
    {
        Y_medium[i] =0.0;
    }
    
    for(i=0;i<input_number;i++)
    {
        Y_medium[i] = GetY(i);
        //cout<<"Y["<<i<<"]="<<Y_medium[i]<<endl;
    }
    
    for(i=0;i<input_number;i++)
    {
        k1[i] = func_medium(t,Y_medium,i,rk);
        tmp[i] = Y_medium[i] + h*k1[i]/2;
        //cout<<"k1["<<i<<"]="<<k1[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    for(i=0;i<input_number;i++)
    {
        k2[i]= func_medium(t+h/2,tmp,i,rk);
        tmp[i] = Y_medium[i] + h*k2[i]/2;
        //cout<<"k2["<<i<<"]="<<k2[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<input_number;i++)
    {
        k3[i] = func_medium(t+h/2,tmp,i,rk);
        tmp[i] = Y_medium[i] + h*k3[i];
        //cout<<"k3["<<i<<"]="<<k3[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    for(i=0;i<input_number;i++)
    {
        k4[i] = func_medium(t+h,tmp,i,rk);
        Y_medium[i] = Y_medium[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
        SetY(Y_medium[i], i);
        //cout<<"k4["<<i<<"]="<<k4[i]<<endl;
        //cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;
    }
    
    delete [] Y_medium ;
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] tmp;
    
    Y_medium = 0;
    k1 = 0;
    k2 = 0;
    k3 = 0;
    k4 = 0;
    tmp = 0;
    
}

double Medium::func_medium(double t, double *y, int i,RK *rk)// input物質の ith 成分
{
    double D_Yin=0.0,r=0.0,r_out=0.0,r_Yin=0.0,x=0.0,totalcellnum=0.0;
    //double D_med=0.0;
    int k;
    //if(i==input_number-1){ D_med = S_D;}//mediumがResourceを受け取る拡散係数
    totalcellnum =rk[0].GetTotalcellnumber();
    for(k=0;k<totalcellnum;k++)
    {
        D_Yin = trspD;//inputの拡散係数を使う
        x = rk[k].GetResult(ProteinVAL-2*input_number +i);//transporter とりあえずラストのタンパクから割り当てておく
        r_Yin= r_Yin - D_Yin*x*y[i]*(double)(rk[k].GetVolume()/V);//能動拡散のみ
        //r_Yin= r_Yin - D_Yin*x*y[i]*(double)(rk[k].GetVolume()/V) - D_med*(y[i]-rk[k].GetYin(i))*(double)(rk[k].GetVolume()/V);
        //受動拡散も入れる
    }
    
    r_out = S_D*(S_out-y[i]);
    r = r_out + r_Yin;
    
    return r;
}

double Medium::func_medium_X(double t, double *x, int i,RK *rk)// input物質の ith 成分
{
    double D_pro=0.0,r_Xin=0.0,totalcellnum=0.0;
    int k;
    
    D_pro = rk[0].GetProD(i);//mediumが拡散物質を受け取る拡散係数
    totalcellnum =rk[0].GetTotalcellnumber();
    
    for(k=0;k<totalcellnum;k++)
    {
        r_Xin = r_Xin - D_pro*(x[i]-rk[k].GetResult(i))*(double)(rk[k].GetVolume()/V);
    }
    
    return r_Xin;
}

void Medium::archive_medium(int time)
{
    
    FILE *pf = nullptr;
    
    pf = fopen("Medium_input.txt","a+");
    
    if(pf==NULL)
    {
        cout << "Could not open file"<<endl;
        exit(1);
    }
    else
    {
        fprintf(pf," %d ",time);
        for(int i=0;i<input_number;i++)
        {
            fprintf(pf," %f",GetY(i));
        }
        fprintf(pf,"\n");
        
    }
    fclose(pf);
}

void Medium::reset()
{
    Y.clear();
    for (int i=0;i<input_number; i++)
    {
        Y[i] = S_resource;//resorverと濃度を合わせておく
    }
}

double func_log10(double sum)
{
    double x=0.0;
    
    if(sum<=0)
    {
        x=0.0;
    }
    else
    {
        x = log10(sum);
    }
    
    return x;
}




void shuffle(int array[], int size)
{
    int i = size;
    while (i > 1) {
        int j = rand() % i;
        i--;
        int t = array[i];
        array[i] = array[j];
        array[j] = t;
    }
}

int main()
{
    const int divnumber=100000000;//evolutionするまでの時間を決める divnumber-fitness_timeから
    const int Totalcellnum =100;//全typeの細胞数
    const int Totalparasitenum = 100;
    const int onecellnum = 100;//一つのtypeの細胞数
    const int oneparasitenum = 100;
    const double startTime=0.0,endTime=5000000.0,beta=25.0,mutnumber=25.0,divinumber=500.0,eps=0.01,medium_D=0.0,threshold = 0.05,sigma = 0.005;
    double SumValues[ProteinVAL] = {},p_threshold[ProteinVAL]={},target_p[target_number]={};
    //double SumTergetValues[MAXCELLNUMBER] = {};
    double q[Totalcellnum] = {},parent_result[ProteinVAL]={};
    int i,j,l,t,k,check_time;//,child1,child2,j_parent=0
    int terget[target_number] = {},jx[ProteinVAL]={};
    double protein_D[ProteinVAL],jpt[ProteinVAL]={},PThreshold[Parasite_ProteinVAL]={};
    double Yin_D[input_number]= {};
    //char filename[50]={};
    int cellnumber=0,now_gene=0,archive_time=0,memo=0,count1=0,count2=0,count3=0;
    int memo_tar[target_number]={};//不正な制御のgeneを記録するよう
    int eliminate_cell[MAXCELLNUMBER]={};//消去する細胞の入れ物
    double time=0.0,dt=0.0;
    int a=15,b=70,c=15;
    int d=0,e=100,f=0;
    int rnt,rnf;
    
    string sStr1,sStr2,sStr3,sStr4,sStr5,sStr6,sStr7,sStr8;
    
    vector<vector<int> > Y_2jigen;
    Y_2jigen.resize(input_number);
    for(int i=0;i<input_number;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            Y_2jigen[i].push_back(1);
        }
    }
   
    vector<vector<double> > pg_2jigen;
    pg_2jigen.resize(Species_numb);
    for(int i=0;i<Species_numb;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            pg_2jigen[i].push_back(0.0);
        }
    }
    
    vector<vector<double> > PJij_2jigen;
    PJij_2jigen.resize(Parasite_ProteinVAL);
    for(int i=0;i<Parasite_ProteinVAL;i++)
    {
        for(int l=0;l<Parasite_ProteinVAL;l++)
        {
            PJij_2jigen[i].push_back(0);
        }
    }
    
    vector<vector<double> > g_rand;
    g_rand.resize(EM_h);
    for(i=0;i<EM_h;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            g_rand[i].push_back(0.0);
        }
    }
    
    vector<vector<double> > g_rand_p;
    g_rand_p.resize(EM_h);
    for(i=0;i<EM_h;i++)
    {
        for(l=0;l<Parasite_ProteinVAL;l++)
        {
            g_rand_p[i].push_back(0.0);
        }
    }
    
    vector<vector<vector< double> > > g_rand_3jigen;
    g_rand_3jigen
    = vector<vector<vector<double> > >(Totalcellnum,vector<vector<double> >(EM_h,vector<double>(ProteinVAL,0)));
    
    for(i=0;i<Totalcellnum;i++)
    {
        for(j=0;j<EM_h;j++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                g_rand_3jigen[i][j][l] = 0.0;
            }
        }
    }
    
    vector<vector<vector<double> > > GRN_3jigen;
    GRN_3jigen
    = vector<vector<vector<double> > >(Totalcellnum,vector<vector<double> >(ProteinVAL,vector<double>(ProteinVAL,0)));
    
    for(i=0;i<Totalcellnum;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                GRN_3jigen[i][t][l]=0;//"GRN%d.txt"　を読み込む
                
            }
        }
    }
    
    vector<vector<vector<double> > > J_3jigen;
    //J_3jigen[i][j][k] i番目の種のGRNの入れ物 J[j][k] from k to j
    J_3jigen
    = vector<vector<vector<double> > >(Totalcellnum,vector<vector<double> >(ProteinVAL,vector<double>(ProteinVAL,0)));
    
    //boost::multi_array<double,2> fitness_2jigen(boost::extents[MAXCELLNUMBER][fitness_time]);
    int sstoi(std::string str);
    
    for(i=0;i<MAXCELLNUMBER;i++)
    {
        eliminate_cell[i]=i;
    }
    
    shuffle(eliminate_cell,MAXCELLNUMBER);
    
    for(i=0;i<MAXCELLNUMBER;i++)
    {
        //cout<<eliminate_cell[i]<<endl;
    }
    
    remove("condition.txt");//パラメタ設定削除
    remove("population.txt");
    remove("Meanfield.txt");
    remove("target.txt");
    
    ofstream ofs1,ofs2;
    ofs1.open("condition.txt");//パラメタ設定記録
    ofs2.open("exit.txt");//全滅終了時の記録
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> randProtein(0, ProteinVAL-1);
    uniform_int_distribution<> randProteinTarget(0, target_number-1);
    uniform_int_distribution<> randProtein_inter(target_number, ProteinVAL-input_number-1);//中間層
    uniform_int_distribution<> randYinJ(target_number+1, ProteinVAL-input_number-1);
    uniform_int_distribution<> randcell(0, MAXCELLNUMBER-1);
    uniform_real_distribution<> rand01(0, 1);
    uniform_real_distribution<> rand_double(-1, 1);
    uniform_int_distribution<> rand_int_01(0, 1);
    uniform_int_distribution<> rand_select(1,100);
    uniform_real_distribution<> randThreshold(threshold,threshold*3.0);
    uniform_real_distribution<> randPThreshold(-threshold,threshold);

    //normalize = double(terget_number);//growthの規格化
    
    Medium medium(startTime,divnumber,endTime,S_D,Totalcellnum,medium_V,S_resource);
    
    CELL cell(startTime,divnumber,endTime,ProteinVAL,target_number,onecellnum,mutnumber,divinumber,Totalcellnum,beta,eps,medium_D,initial_volume,threshold);
    
    cellnumber = cell.rk[0].GetTotalcellnumber();
    
    cell.SetData_array();
    
    for(l=0;l<ProteinVAL;l++)
    {
        SumValues[l] = 0.0;
        parent_result[l] =0.0;
    }
    for(i=0;i<Species_numb;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                //初期のネットワーク構造 作成 連続変化の時
                J_3jigen[i][l][t] = rand_double(mt);
            }
        }
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                //初期のネットワーク構造 作成　離散値
                rnt = rand_select(mt);
                if(rnt>=1 && rnt<=a)
                {
                    rnf = 1;
                }
                else if(rnt>=a && rnt<=a+b)
                {
                    rnf = 0;
                }
                else
                {
                    rnf = -1;
                }
                J_3jigen[i][l][t] = rnf;
                
            }
        }
        
        for(l=0;l<Parasite_ProteinVAL;l++)
        {
            for(t=0;t<Parasite_ProteinVAL;t++)
            {
                //初期のネットワーク構造 作成
                rnt = rand_select(mt);
                if(rnt>=1 && rnt<=a)
                {
                    rnf = 1;
                }
                else if(rnt>=a && rnt<=a+b)
                {
                    rnf = 0;
                }
                else
                {
                    rnf = -1;
                }
                PJij_2jigen[l][t] = rnf;
            }
            PThreshold[l] =randPThreshold(mt);
        }
    }
    
    
    for(i=0;i<Species_numb;i++)
    {
        for(l=0;l<target_number;l++)
        {
            J_3jigen[i][l][l] = 0;//target geneは自己複製できない
            for(t=0;t<ProteinVAL;t++)
            {
                J_3jigen[i][t][l] = 0;//target geneからのFBはない
            }
        }
        
        for(l=0;l<input_number;l++)
        {
            for(t=0;t<target_number;t++)
            {
                J_3jigen[i][t][ProteinVAL-1-l] = 0;
                //input geneからoutput geneへの直結はない
            }
        }
    }
    
    
    for(i=0;i<Species_numb;i++)
    {
        for(j=onecellnum*i;j<onecellnum*(i+1);j++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                for(t=0;t<ProteinVAL;t++)
                {
                    cell.rk[j].SetJ(J_3jigen[i][l][t],l,t);//初期のネットワーク構造 GRN form t to l
                }
            }
        }
    }
    
    for(i=0;i<Species_numb;i++)
    {
        cell.rk[onecellnum*i].converterJ();//隣接リストを作成
        //cell.rk[onecellnum*i].calculate_Input_list();//inputからの距離
    }
    
    for(i=0;i<Species_numb;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            if(cell.rk[onecellnum*i].GetRSEL(l)>=MaxPathNum+1)//MaxPathNum本以上出て入れば
            {
                while(cell.rk[onecellnum*i].GetRSEL(l)!=MaxPathNum)//MaxPathNum本にする
                {
                    k=randProtein(mt);
                    while(cell.rk[onecellnum*i].GetJ(k,l)==0)
                    {
                        k=randProtein(mt);
                    }
                    
                    cell.rk[onecellnum*i].SetJ(0,k,l);
                    cell.rk[onecellnum*i].converterJ();//隣接リストを作成
                    //cell.rk[onecellnum*i].calculate_Input_list();
                }
            }
            if(cell.rk[onecellnum*i].GetSEL(l)>=MaxPathNum+1)//MaxPathNum本以上入れば
            {
                while(cell.rk[onecellnum*i].GetSEL(l)!=MaxPathNum)//MaxPathNum本にする
                {
                    k=randProtein(mt);
                    while(cell.rk[onecellnum*i].GetJ(l,k)==0)
                    {
                        k=randProtein(mt);
                    }
                    cell.rk[onecellnum*i].SetJ(0,l,k);
                    cell.rk[onecellnum*i].converterJ();//隣接リストを作成
                    //cell.rk[onecellnum*i].calculate_Input_list();
                }
            }
        }
    }
    
    for(i=0;i<Species_numb;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                J_3jigen[i][l][t] = cell.rk[onecellnum*i].GetJ(l,t);//初期のネットワーク構造 GRN form t to l
            }
        }
    }
    
    for(i=0;i<Species_numb;i++)
    {
        for(t=0;t<ProteinVAL;t++)
        {
            count1=0;count2=0;//カウンターをリセット
            for(l=0;l<target_number;l++)
            {
                if(J_3jigen[i][l][t] == 1)
                {
                    count1 = count1 + 1;
                    memo_tar[l] = l;
                }
                else{
                    memo_tar[l] = 0;
                }
            }
            if(count1>MaxTargetPath)
            {
                for(l=0;l<count1-MaxTargetPath;l++)
                {
                    j=randProteinTarget(mt);
                    while(memo_tar[j] != 0)
                    {
                        J_3jigen[i][j][t] = 0;
                        j=randProteinTarget(mt);
                    }
                }
            }
        }
    }
    
    for(i=0;i<Species_numb;i++)
    {
        for(j=onecellnum*i;j<onecellnum*(i+1);j++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                for(t=0;t<ProteinVAL;t++)
                {
                    cell.rk[j].SetJ(J_3jigen[i][l][t],l,t);//初期のネットワーク構造 GRN form t to l
                    cell.rk[j].SetStationaryJ(cell.rk[j].GetJ(l,t), l, t);//分裂時のネットワークを更新
                    cell.SetJij_initila(J_3jigen[i][l][t],l,t);
                    //cout<<"J_3jigen["<<i<<"]["<<l<<"]["<<t<<"]="<<cell.rk[0].GetJ(l,t)<<endl;
                }
            }
        }
    }
    
    for(j=0;j<Totalcellnum;j++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                cell.rk[j].SetJ(J_3jigen[0][l][t],l,t);//初期のネットワーク構造 GRN form t to l
                //cell.rk[j].SetStationaryJ(cell.rk[j].GetJ(t,l), t, l);//分裂時のネットワークを更新
                //cout<<"J_3jigen["<<i<<"]["<<l<<"]["<<t<<"]="<<cell.rk[0].GetJ(l,t)<<endl;
            }
        }
    }
    
    
    for(i=0;i<Totalcellnum;i++)
    {
        cell.rk[i].converterJ();//隣接リストを作成
    }
    
    for(i=0;i<Totalcellnum;i++)
    {
        //cell.rk[i].calculate_Input_list();
    }
    
    for(l=0;l<input_number;l++)
    {
        for(t=0;t<ProteinVAL;t++)
        {
            Y_2jigen[l][t]=0;
        }
    }
    
    for(l=0;l<input_number;l++)
    {
        Y_2jigen[l][ProteinVAL-l-1] = 1;//最後から順にinputとして選ぶ
    }
    
    for(i=0;i<cellnumber;i++)
    {
        for(l=0;l<input_number;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                cell.rk[i].SetYin_J(Y_2jigen[l][t],l,t);
            }
        }
    }
    
    for(l=0;l<input_number;l++)
    {
        for(t=0;t<ProteinVAL;t++)
        {
            cout<<cell.rk[0].GetYin_J(l,t);
        }
        cout<<endl;
    }
    
    for(l=0;l<target_number;l++)
    {
        memo_tar[l] = 0;
    }
    
    count1=0;count2=0;count3=0;//カウンターをリセット
    for(i=0;i<1;i++)
    {
        for(t=0;t<ProteinVAL;t++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                if(J_3jigen[i][l][t] == 1)
                {
                    count1 = count1 + 1;
                }
                else if(J_3jigen[i][l][t] == -1)
                {
                    count2 = count2 + 1;
                }
                else
                {
                    count3 = count3 + 1;
                }
            }
        }
    }
    
    FILE *pf = nullptr;
    
    char filename[50]={};
    
    for(i=0;i<Species_numb;i++)
    {
        sprintf(filename,"J_3jigen%d.txt",i+1);
        pf = fopen(filename,"a+");
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                fprintf(pf," %f ",cell.rk[i].GetJ(t,l));
            }
            fprintf(pf,"\n");
        }
        fclose(pf);
    }
    
    for(i=0;i<Species_numb;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                jx[t] =J_3jigen[i][t][l];
                //archive_J_3jigen_graph(J_3jigen[i][t][l],i);//graph用
            }
        }
    }
    
    for(i=0;i<input_number;i++)
    {
        //cell.rk[0].archive_YinJ(i);//input物質の活性対象を表示
    }
    
    ifstream ifs;
    for(int species=0;species<Species_numb;species++)
    {
        //sprintf(filename,"J_3jigen%d.txt",species+1);
        //ifs.open(filename);
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                //ifs >>J_3jigen[species][t][l];//"J_3jigen%d.txt"　を読み込む
            }
        }
        //ifs.close();
    }
    
    ifstream ifs1,ifs2,ifs3,ifs4,ifs5,ifs6,ifs7,ifs8;
    
    //ここから
    
    for(l=0;l<ProteinVAL;l++)
    {
        //protein_D[l] = medium_D*0.5+rand01(mt)*medium_D*0.5;
        //protein_D[l] = medium_D;//tergetのみ拡散しない
        //protein_D[l] = medium_D*rand01(mt)*cell.rk[0].rand_select(d, e, f);
        protein_D[l] = medium_D*cell.rk[0].rand_select(d, e, f);
        //protein_D[l] = 0;//tergetのみ拡散する
        for(i=0;i<Totalcellnum;i++)
        {
            cell.rk[i].SetProteinD(protein_D[l],l);
        }
    }
    
    for(l=ProteinVAL-2*input_number;l<ProteinVAL-input_number;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            //protein_D[l] = 0;
            //cell.rk[i].SetProteinD(protein_D[l],l);//inputも拡散しない
        }
    }
    
    for(l=ProteinVAL-input_number;l<ProteinVAL;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            // protein_D[l] = medium_D;//transporterは拡散する
            //cell.rk[i].SetProteinD(protein_D[l],l);
        }
    }
    
    for(l=0;l<input_number;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            //protein_D[l] = medium_D;//transporterは拡散する
            //cell.rk[i].SetProteinD(protein_D[l],l);
        }
    }
    
    for(l=0;l<target_number;l++)
    {
        terget[l] = l;
        for(i=0;i<Totalcellnum;i++)
        {
            //protein_D[terget[l]] =medium_D;//tergetのみ拡散する
            //protein_D[terget[l]] =0;//tergetのみ拡散しない
            //cell.rk[i].SetProteinD(protein_D[terget[l]],terget[l]);
            //protein_D[ProteinVAL-1] =0;
        }
    }
    
    //protein_D[ProteinVAL-input_number+2]=0.0;
    //protein_D[ProteinVAL-input_number+3]=0.2;
    
    for(l=0;l<ProteinVAL;l++)
    {
        archive_ProteinD(protein_D[l],0);//拡散係数を記録
    }
    
    for(l=0;l<input_number;l++)
    {
        Yin_D[l]=trspD;//ここでinput物質の拡散係数を決める
    }
    //Yin_D[input_number-1]=0.0;//受動拡散のみ
    
    for(l=0;l<input_number;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            cell.rk[i].SetMedium_D(Yin_D[l], l);//ここでinput物質の拡散係数を与える
        }
    }
    
    for(int i=0;i<Species_numb;i++)
    {
        for(int l=0;l<ProteinVAL;l++)
        {
            pg_2jigen[i][l] = randThreshold(mt);
        }
        
        for(int l=0;l<input_number;l++)
        {
            pg_2jigen[i][ProteinVAL-1-l] = input_strength;//つねに入力は1になるように十分とっておく
        }
        
        for(int l=0;l<ProteinVAL;l++)
        {
            jpt[l] = pg_2jigen[i][l];
        }
        archive_threshold(jpt,i);
    }
    
    for(l=0;l<ProteinVAL;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            cell.rk[i].SetThreshold_g(pg_2jigen[0][l],l);//0番目のcellの閾値を使用
        }
        p_threshold[l] =pg_2jigen[0][l];
    }
    
    for(i=0;i<Totalcellnum;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                memo =cell.rk[i].GetJ(l,t);
                if(memo == -1)
                {
                    cell.rk[i].ori_Sstr_J = cell.rk[i].ori_Sstr_J + 'a' ;//初期のネットワーク構造 GRN -1を文字aで記録
                }
                else if(memo == 0)
                {
                    cell.rk[i].ori_Sstr_J = cell.rk[i].ori_Sstr_J + 'b' ;//初期のネットワーク構造 GRN 0を文字bで記録
                }
                else
                {
                    cell.rk[i].ori_Sstr_J = cell.rk[i].ori_Sstr_J + 'c' ;//初期のネットワーク構造 GRN 1を文字cで記録
                }
            }
        }
    }
    
    for(i=0;i<Totalcellnum;i++)
    {
        cell.rk[i].Sstr_J = cell.rk[i].ori_Sstr_J;
    }
    
    for(t=0;t<target_number;t++)
    {
        //target_p[t] = (double)rand_int_01(mt);//01パターン
        target_p[t] =1;//全てつける必要
    }
    
    for(i=0;i<Totalcellnum;i++)
    {
        for(t=0;t<target_number;t++)
        {
            cell.rk[i].SetTargetPattern(target_p[t],t);
        }
    }
    
    //ここまで terget setup
    for(i=0;i<Totalcellnum;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            for(t=0;t<ProteinVAL;t++)
            {
                //cell.rk[i].SetK_delta(1,ProteinVAL-1);//最後のgeneをinputにする
            }
        }
    }
    
    for(l=0;l<ProteinVAL;l++)
    {
        SumValues[l] = 0;
        for(i=0;i<Totalcellnum;i++)
        {
            SumValues[l] = SumValues[l] + cell.rk[i].GetResult(l);
        }
        for(i=0;i<Totalcellnum;i++)
        {
            //cell.rk[i].SetMeanfieldval(SumValues[l]/cellnumber,l);
        }
        //cout<<"SumValues["<<l<<"]="<<SumValues[l]<<endl;
    }
    
    for(i=0;i<Totalcellnum;i++)
    {
        for(l=0;l<ProteinVAL;l++)
        {
            //cell.rk[i].SetProteinD(protein_D[l],l);
        }
    }
    
    for(l=0;l<ProteinVAL;l++)
    {
        for(i=0;i<Totalcellnum;i++)
        {
            //cell.rk[i].SetThreshold_g(p_threshold[l],l);
        }
    }
    
    ofs1 << "#MAXCELLNUMBER = "<<MAXCELLNUMBER<<endl;
    ofs1 << "#MAXPARASITENUMBER = "<<MAXPARASITENUMBER<<endl;
    ofs1 << "#Species_numb = "<<Species_numb<<endl;
    ofs1 << "#ProteinVAL = "<<ProteinVAL<<endl;
    ofs1 << "#MaxPathNum = "<<MaxPathNum<<endl;
    ofs1 << "#MaxTargetPath = "<<MaxTargetPath<<endl;
    ofs1 << "#target_number = "<<target_number<<endl;
    ofs1 << "#input_number = " <<input_number<<endl;
    ofs1 << "#input_strength = "<<input_strength<<endl;
    //ofs1 << "#S_resource　=　"<<S_resource<<endl;
    ofs1 << "#death_time = "<<death_time<<endl;
    ofs1 << "#death_rate = "<<death_rate<<endl;
    ofs1 << "#a = "<<a<<endl;
    ofs1 << "#b = "<<b<<endl;
    ofs1 << "#c = "<<c<<endl;
    ofs1 << "#d = "<<d<<endl;
    ofs1 << "#e = "<<e<<endl;
    ofs1 << "#f = "<<f<<endl;
    ofs1 << "#activate rate= "<<(double)count1/(ProteinVAL*ProteinVAL)<<endl;
    ofs1 << "#inhibit rate = "<<(double)count2/(ProteinVAL*ProteinVAL)<<endl;
    ofs1 << "#no-interaction rate = "<<(double)count3/(ProteinVAL*ProteinVAL)<<endl;
    
    for(l=0;l<ProteinVAL;l++)
    {
        ofs1 <<"#proteinD["<<l+1<<"] = "<<protein_D[l]<<endl;
    }
    for(l=0;l<input_number;l++)
    {
        ofs1 <<"#Yin_D["<<l<<"]="<<Yin_D[l]<<endl;
    }
    for(l=0;l<cell.rk[0].GetInput_listSize(0);l++)
    {
        //ofs1 <<"#Input_list["<<l<<"]="<<cell.rk[0].GetInput_list(0,l)<<endl;
    }
    for(l=0;l<target_number;l++)
    {
        ofs1 <<"#target_pattern["<<l<<"]="<<target_p[l]<<endl;
    }
    ofs1 << "#Route = "<<Route<<endl;
    ofs1 << "#S_D="<<S_D<<endl;
    ofs1 << "#medium_V="<<medium_V<<endl;
    ofs1 << "#initial_volume = "<<initial_volume<<endl;
    for(l=0;l<ProteinVAL;l++)
    {
        ofs1 <<"#threshold["<<l<<"]="<<p_threshold[l]<<endl;
    }
    ofs1 << "#initial_volume = "<<initial_volume<<endl;
    ofs1 << "#beta = "<<beta<<endl;
    ofs1 << "#growth_const="<<growth_const<<endl;
    ofs1 << "#Sigma = "<<Sigma<<endl;
    ofs1 << "EM_h ="<<EM_h<<endl;//EMでの分割数
    ofs1 << "#dt = "<<cell.rk[0].GetDelt()<<endl;
    ofs1 << "#histdiv = "<<histdiv<<endl;//ヒストグラムの分割数
    ofs1 << "#elite_num = "<<elite_num<<endl;//次世代に残す個体数
    ofs1 << "#mutation = "<<mutation_num<<endl;//変異を与える個体数
    ofs1 << "#selection_num ="<<selection_num<<endl;//下位の集団は取り除く
    ofs1 << "#Host_mutaion_rate = "<<Host_mutaion_rate<<endl;
    ofs1 << "#Host_Sigma_mutaion_rate = "<<Host_Sigma_mutaion_rate<<endl;
    ofs1 << "#parasite_generation_times_over = "<<parasite_generation_times_over<<endl;
    ofs1 << "#parasite_generation_times = "<<parasite_generation_times<<endl;
    ofs1 << "#parasite_population_dynamics_cal_times = "<<parasite_population_dynamics_cal_times<<endl;
    
    ofs1 << "#Parasite_genome_size = "<<Parasite_genome_size<<endl;//parasiteの遺伝子型のサイズ
    ofs1 << "#Parasite_target_genome = "<<Parasite_target_genome<<endl;//cellのnontarget_geneを選ぶとき
    ofs1 << "#maximum_virulence = "<<maximum_virulence<<endl;//parasiteの毒性
    
    //printf("現在のディスク空きは%f%%です．\n", checkDiskSize(1024)*100);
    
    //第０秒目だけ先に計算
    cellnumber=cell.rk[0].GetTotalcellnumber();
    for(i=0;i<cellnumber;i++)
    {
        //cell.rk[i].ShowDate();
        //cell.rk[i].archive();
        //cell.rk[i].calculateNS_protein();
        //cell.rk[i].calculateNS_protein_MF();
        //cell.rk[i].calculateNS_protein_g();
        for(j=0;j<EM_h;j++)
        {
            for(l=0;l<ProteinVAL;l++)
            {
                g_rand[j][l] = Sigma*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));
            }
        }
        cell.rk[i].calculateNS_protein_g_EM(g_rand);
        for(l=0;l<input_number;l++)
        {
            //cell.rk[i].SetY(medium.GetY(l),l);//medium classのYをrk classのYに置き換え
        }
        for(l=0;l<ProteinVAL;l++)
        {
            //cell.rk[i].SetX(medium.GetX(l),l);//medium classのXをrk classのXに置き換え
        }
        //cell.rk[i].calculateNS_Yin();
        cell.rk[i].calculateNS_volume();
        //cell.rk[i].calculate_Mu();
    }
    //cell.calculate_meanfield();
    //cell.genosim_label_set();
    cell.caluculate_hosto_num();
    cell.population_dynamics();
    
    //cell.archive_parasite_population(0);
    //cell.rk[0].archive_GRN_every_gene(0,0);//初期のネットワークを記録
    
    time = cell.rk[0].GetTime();
    dt = cell.rk[0].GetDelt();
    for(i=0;i<cellnumber;i++)
    {
        //time = cell.rk[i].GetTime();
        cell.rk[i].SetTime(time+dt);//時間更新
        cell.rk[i].SetTime1000();
    }
    //ここまで
    
    //以下、計算
    
    for(t=1;t<divnumber;t++)
    {
        cellnumber=cell.rk[0].GetTotalcellnumber();
        for(i=0;i<cellnumber;i++)
        {
            //cell.rk[i].Show();
            //cell.rk[i].archive();
            //cell.rk[i].calculateNS_protein_MF();
            //cell.rk[i].calculateNS_protein();
            
            for(j=0;j<EM_h;j++)
            {
                for(l=0;l<ProteinVAL;l++)
                {
                    g_rand[j][l] = Sigma*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));//乱数作成をここで行う
                    //g_rand[j][l] =cell.rk[i].GetSigma_Str()*sqrt(-2.0*log(rand01(mt)))*sin(2.0*PI*rand01(mt));//sigmaが進化する場合
                }
            }
            
            cell.rk[i].calculateNS_protein_g_EM(g_rand);
            //parasite.para[i].calculate_Result(g_rand_p);
            //cell.rk[i].calculateNS_protein_delution_g_EM(g_rand);
            //cell.rk[i].calculateNS_protein_signal();
            //cell.rk[i].calculateNS_protein_signal_g_EM(g_rand);
            
            for(l=0;l<input_number;l++)
            {
                //cell.rk[i].SetY(medium.GetY(l),l);//medium classのYをrk classのYに置き換え
            }
            
            //cell.rk[i].calculateNS_Yin();
            cell.rk[i].calculate_growth();
            //cell.rk[i].calculate_growth_dist();
            cell.rk[i].calculateNS_volume();
            //cell.rk[i].calculate_Mu();
        }
        
        //medium.calculateNS_med(cell.rk);//mediumの計算
        //cell.calculate_meanfield();//拡散平均場の計算
        //cell.calculate_NFDS_growth();
        
        cellnumber=cell.rk[0].GetTotalcellnumber();
        now_gene = cell.rk[0].GetGene();
        
        if(t% ( 20*host_generation_times)==0)
        {
            cell.caluculate_hosto_num();
        }//ここで少しparasiteのpopulationを計算する
        
        if(t% ( 20*parasite_generation_times)==0)
        {
            cell.calculate_parasite_dynamics();
            //cell.calculate_parasite_dynamics_dist_fitness();
            //cell.calculate_parasite_dynamics_evo();
        }
        
        for(i=0;i<cellnumber;i++)
        {
            //cell.rk[i].Calculate_Stack_fitness_Target_parasite_average();//時間平均のfitnessを計算して値を確保しておく
            //cell.rk[i].Calculate_Stack_fitness_Target_noparasite_average();//nonattackedのfitness平均
        }
        
        //cell.calculate_Result_average();
        //cell.calculate_Fitness_Target_parasite_average();//Stackを足し合わせるFitness_averageをここで規格化
        //cell.calculate_Fitness_Target_noparasite_average();
        
        if(t % 200 == 0)//10秒ごと
        {
            cell.archive_parasite_population(t/20);
            //cell.archive_sim(t/20);//division_fitness_mutationで行う
            //cell.archive_growth_volume(t/20);
            for(j=0;j<cellnumber;j++)
            {
                //cout<<cell.rk[j].GetFitness_Target_parasite_average()<<endl;
            }
        }
        
        if(t%100000==0||t%250000==0)//5000秒 -->10世代
        {
            for(i=0;i<cellnumber;i++)
            {
                cell.rk[i].archive_GRN_every_gene(i,t/20);
            }
            //cell.rk[0].archive_GRN_every_gene_one(0,t/20);
            cout<<"generarion(every 10)="<<t/20<<endl;
            if(cellnumber<10)
            {
                ofs2 << "extinction"<<endl;
                exit(-1);
            }
        }
        
        cell.death_volume_check();
        
        for(i=0;i<cellnumber;i++)
        {
            //cell.division_volume_mutation_network_sigma(i,sigma,t/20);//sigmaが進化する
            cell.division_volume_mutation_network(i,sigma,t/20);//volume cell_sizeは固定
            //cell.division_volume_mutation_network_continuas(i,sigma,t/20);//volume cell_sizeは固定
            //cell.division_volume_mutation_network_couple(i,sigma,t/20);//volume
        }
        
        if(t%50000==0)//2500秒 -->5世代
        {
            cell.archive_Jij_Variance(t/20);
            //cell.division_fitness_mutation(sigma,t/20);
            //cell.archive_sim_sta(t/20);//division_fitness_mutationで行う
            //cell.archive_sim(t/20);//division_fitness_mutationで行う
            cell.archive_var(t/20);//division_fitness_mutationで行う
            //cell.archive_var_threshold(t/20);//division_fitness_mutationで行う
            //cell.archive_sim(t/20);//division_fitness_mutationで行う
            cell.archive_var_genome(t/20);
            //cell.archive_Sigma_Str_all(t/20);
            //cell.archive_var_sta(t/20);//division_fitness_mutationで行う
        }
        
        if(t%100000==0)//5000秒 -->10世代
        {
            //cell.archive_Sigma_Str(t/20);
            //cell.archive_threshold(t/20);
        }
        
        if(t%500000==0)//25000秒 -->50世代
        {
            //cell.archive_Jij_Variance(t/20);
            //cell.division_fitness_mutation(sigma,t/20);
            //cell.archive_sim_sta(t/20);//division_fitness_mutationで行う
            cell.archive_sim(t/20);//division_fitness_mutationで行う
            //cell.archive_var(t/20);//division_fitness_mutationで行う
            //cell.archive_var_sta(t/20);//division_fitness_mutationで行う
        }
        
        if(t%10000==0)//500秒-->1世代
        {
            //cell.archive_Jij_Variance(t/20);
            //cell.division_fitness_mutation(sigma,t/20);
            //cell.archive_sim_sta(t/20);//division_fitness_mutationで行う
            //cell.archive_sim(t/20);//division_fitness_mutationで行う
            //cell.archive_var(t/20);//division_fitness_mutationで行う
            //cell.archive_var_sta(t/20);//division_fitness_mutationで行う
            
            //cell.archive_Parasite_Virulence(t/20);
            //cell.archive_growth_volume(t/20);
            cout<<"generarion="<<t/(20*divinumber)<<endl;
            
            //medium.reset();
            //cell.reset();
            for(i=0;i<cellnumber;i++)
            {
                cell.rk[i].SetGene(cell.rk[i].GetGene()+1);
            }
        }
        
        dt = cell.rk[0].GetDelt();
        for(i=0;i<cellnumber;i++)
        {
            time = cell.rk[i].GetTime();
            cell.rk[i].SetTime(time+dt);//時間更新
            cell.rk[i].SetTime1000();
        }
    }
    
    ofs1.close();
    return 0;
}
