#include "Raw_Event.h"
//#include "AidaMover.h" 
#include <iostream>

using namespace std;

//---------------------------------------------------------------

//Raw_Event::Raw_Event(bool PADI_OR_PADIWA) : PLASTIC_Data(),PLASTIC_VME_Data(){
    Raw_Event::Raw_Event(){


    Event_Type = -1;

//     PLASTIC_Data.PADI_OR_PADIWA = PADI_OR_PADIWA;

}

//---------------------------------------------------------------

Raw_Event::~Raw_Event(){}

//---------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- FRS  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


void Raw_Event::set_DATA_MUSIC(Float_t* FRS_dE,Float_t* FRS_dE_cor, Int_t* MusicE1,Int_t* MusicE2,Int_t* MusicT1, Int_t* MusicT2 ){
  for(int i=0; i<2; ++i){

    MUSIC_dE[i] = 0;
    MUSIC_dE_cor[i] = 0;

    }
    for(int i=0; i<8;i++){
        MUSIC_e1[i] = 0;
        MUSIC_e2[i] = 0;
        MUSIC_t1[i] = 0;
        MUSIC_t2[i] =0;
    }
    
    
    for(int i=0; i<2; ++i){

    MUSIC_dE[i] = FRS_dE[i];
    MUSIC_dE_cor[i] = FRS_dE_cor[i];

    }
    for(int i=0; i<8;i++){
        MUSIC_e1[i] = MusicE1[i];
        MUSIC_e2[i] = MusicE2[i];
        MUSIC_t1[i] = MusicT1[i];
        MUSIC_t2[i] = MusicT2[i];
    }

}
void Raw_Event::set_DATA_SCI(Float_t* FRS_sci_l,Float_t* FRS_sci_r,Float_t* FRS_sci_e,Float_t* FRS_sci_tx,Float_t* FRS_sci_x){
for(int i=0; i <6; i++){
    sci_l[i] = 0;
    sci_r[i] = 0;
    sci_e[i] = 0;
    sci_tx[i] = 0;
    sci_x[i] = 0;
}
    
    for (int cnt=0;cnt<6;cnt++) // 
     {
       int idx = 0 ;
       //int mw_idx = 0;
       //Float_t mwx = 0;
       switch(cnt)
     {
     case 0:        /* SC21 */
       idx = 2; 
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;    
     case 1:        /* SC21 delayed */
       idx = 3; 
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;    
     case 2:        /* SC41 */
       idx = 5; 
       //mw_idx = 5;
       //mwx = clb.tpc_sc41_x;
       break;    
     case 3:        /* SC42 */
           idx = 6;
       break;
     case 4:
       idx = 7;     /* SC43 */
       break;
     case 5:
       idx = 10;    /* SC81 */
       break;
     default: idx = 2;
     }   

    sci_l[idx] = FRS_sci_l[idx];
    //cout<<"RAW sci_l[2]" <<  sci_l[2]  <<endl;
    sci_r[idx] = FRS_sci_r[idx];
    sci_e[idx] = FRS_sci_e[idx];
    sci_tx[idx] = FRS_sci_tx[idx];
    sci_x[idx] = FRS_sci_x[idx];

    }

}

void Raw_Event::set_DATA_SCI_dT(Int_t FRS_dt_21l_21r, Int_t FRS_dt_41l_41r,
                Int_t FRS_dt_21l_41l, Int_t FRS_dt_21r_41r,
                Int_t FRS_dt_42l_42r, Int_t FRS_dt_43l_43r,
                Int_t FRS_dt_21l_42l, Int_t FRS_dt_21r_42r,
                Int_t FRS_dt_81l_81r, Int_t FRS_dt_21l_81l, Int_t FRS_dt_21r_81r){


     dt_21l_21r = 0;
    dt_41l_41r = 0;
    dt_21l_41l = 0;
    dt_21r_41r = 0;
    dt_42l_42r = 0;
    dt_43l_43r = 0;
    dt_21l_42l = 0;
    dt_21r_42r = 0;
    dt_81l_81r = 0;
    dt_21l_81l = 0;
    dt_21r_81r = 0;
    
    
    dt_21l_21r = FRS_dt_21l_21r;
    dt_41l_41r = FRS_dt_41l_41r;
    dt_21l_41l = FRS_dt_21l_41l;
    dt_21r_41r = FRS_dt_21r_41r;
    dt_42l_42r = FRS_dt_42l_42r;
    dt_43l_43r = FRS_dt_43l_43r;
    dt_21l_42l = FRS_dt_21l_42l;
    dt_21r_42r = FRS_dt_21r_42r;
    dt_81l_81r = FRS_dt_81l_81r;
    dt_21l_81l = FRS_dt_21l_81l;
    dt_21r_81r = FRS_dt_21r_81r;




}
void Raw_Event::set_DATA_SCI_ToF(Float_t FRS_sci_tofll2,Float_t FRS_sci_tofll3,Float_t FRS_sci_tof2,Float_t FRS_sci_tofrr2,Float_t FRS_sci_tofrr3,Float_t FRS_sci_tof3){
    sci_tofll2 = 0;
    sci_tofll3 = 0;
    sci_tof2   = 0;
    sci_tofrr2 = 0;
    sci_tofrr3 = 0;
    sci_tof3   = 0;
    
    
    sci_tofll2 = FRS_sci_tofll2;
    sci_tofll3 = FRS_sci_tofll3;
    sci_tof2   = FRS_sci_tof2;
    sci_tofrr2 = FRS_sci_tofrr2;
    sci_tofrr3 = FRS_sci_tofrr3;
    sci_tof3   = FRS_sci_tof3;

}
void Raw_Event::set_DATA_TPC(Int_t** FRS_TPC_lt,Int_t** FRS_TPC_rt,Float_t* FRS_TPC_x, Float_t* FRS_TPC_y ,Float_t FRS_TPC_x0,Float_t FRS_TPC_x1){
 
    TPC_x0 = 0;
    TPC_x1 = 0;
for(int i=0; i<7; i++){
        TPC_x[i] = 0;
        TPC_y[i] = 0;
         for(int j=0; j<2; j++){
            TPC_lt[i][j] = 0;
            TPC_rt[i][j] = 0;
        }
    }
     for(int i=0; i<7; i++){
        TPC_x[i] = FRS_TPC_x[i];
        TPC_y[i] = FRS_TPC_y[i];
        for(int j=0; j<2; j++){
         TPC_lt[i][j] =  FRS_TPC_lt[i][j];
         TPC_rt[i][j] =  FRS_TPC_rt[i][j];
        }
    }
    TPC_x0 = FRS_TPC_x0;
    TPC_x1 = FRS_TPC_x1;
    
}
 
void Raw_Event::set_DATA_ID_2_4(Float_t FRS_ID_x2,Float_t FRS_ID_y2,Float_t FRS_ID_a2,Float_t FRS_ID_b2,Float_t FRS_ID_x4,Float_t FRS_ID_y4,Float_t FRS_ID_a4,Float_t FRS_ID_b4){
    ID_x2 = 0;
    ID_y2 = 0;
    ID_a2 = 0;
    ID_b2 = 0;
    ID_x4 = 0;
    ID_y4 = 0;
    ID_a4 = 0;
    ID_b4 = 0;
   
    
    ID_x2 = FRS_ID_x2;
    ID_y2 = FRS_ID_y2;
    ID_a2 = FRS_ID_a2;
    ID_b2 = FRS_ID_b2;
    ID_x4 = FRS_ID_x4;
    ID_y4 = FRS_ID_y4;
    ID_a4 = FRS_ID_a4;
    ID_b4 = FRS_ID_b4;

}
void Raw_Event::set_DATA_ID_Beta_Rho(Float_t* FRS_ID_brho,Float_t* FRS_ID_rho,Float_t FRS_beta,Float_t FRS_beta3,Float_t FRS_gamma){

    for(int i=0; i<2; ++i){

    ID_brho[i] = 0;
    ID_rho[i] = 0;

    }
    beta = 0;
    beta3 = 0;
    gamma = 0;
    
    for(int i=0; i<2; ++i){
   // cout<<"ID_brho[i] " << ID_brho[i]<<" i " << i << endl;
      
    ID_brho[i] = FRS_ID_brho[i];
    ID_rho[i] = FRS_ID_rho[i];

    }

    beta = FRS_beta;
    beta3 = FRS_beta3;
    gamma = FRS_gamma;

}
void Raw_Event::set_DATA_ID_Z_AoQ(Float_t FRS_AoQ,Float_t FRS_AoQ_corr,Float_t FRS_z,Float_t FRS_z2,Float_t FRS_z3,Float_t FRS_dEdeg,Float_t FRS_dEdegoQ){
    AoQ = 0;
    AoQ_corr = 0;
    z=0;
    z2=0;
    z3=0;
    dEdeg=0;
    dEdegoQ=0;
    AoQ = FRS_AoQ;
    AoQ_corr = FRS_AoQ_corr;
    z = FRS_z;
    z2 = FRS_z2;
    z3 = FRS_z3;
    dEdeg = FRS_dEdeg;
    dEdegoQ = FRS_dEdegoQ;

}
void Raw_Event::set_DATA_ID_Timestamp(Float_t FRS_timestamp,Float_t FRS_ts,Float_t FRS_ts2){

    timestamp = FRS_timestamp;
    ts = FRS_ts;
    ts2 = FRS_ts2;

}

void Raw_Event::set_DATA_FRS_SCALERS(Int_t FRS_time_in_ms, Int_t FRS_spill_count,Int_t FRS_ibin_for_s, Int_t FRS_ibin_for_100ms, Int_t FRS_ibin_for_spill,Int_t FRS_ibin_clean_for_s, Int_t FRS_ibin_clean_for_100ms,Int_t FRS_ibin_clean_for_spill, UInt_t* FRS_increase_scaler_temp){
    
    time_in_ms = FRS_time_in_ms;
    spill_count = FRS_spill_count;
    ibin_for_s = FRS_ibin_for_s;
    ibin_for_100ms = FRS_ibin_for_100ms;
    ibin_for_spill = FRS_ibin_for_spill;
    ibin_clean_for_s = FRS_ibin_clean_for_s;
    ibin_clean_for_100ms = FRS_ibin_clean_for_100ms;
    ibin_clean_for_spill = FRS_ibin_clean_for_spill;
       for (int i = 0; i < 64; i++) increase_scaler_temp[i] = FRS_increase_scaler_temp[i];
  
}

void Raw_Event::set_DATA_VFTX(Double_t vftx_21l, Double_t vftx_21r,Double_t vftx_22l, Double_t vftx_22r,Double_t vftx_41l, Double_t vftx_41r, Double_t vftx_42l, Double_t vftx_42r, Double_t* vftx_t){
    TRaw_vftx_21L = vftx_21l;
    TRaw_vftx_21R = vftx_21r;
    TRaw_vftx_22L = vftx_22l;
    TRaw_vftx_22R = vftx_22r;
    TRaw_vftx_41L = vftx_41l;
    TRaw_vftx_41R = vftx_41r;
    TRaw_vftx_42L = vftx_42l;
    TRaw_vftx_42R = vftx_42r;
    for(int i=0; i<32; i++){
     TRaw_vftx[i] = vftx_t[i];   
    }
}

// #################################################################

// int Raw_Event::get_Event_type(){
//     return Event_Type;
// }
//
//
// void Raw_Event::set_DATA_FATIMA(int QDC_FIRED,int TDC_FIRED,std::vector<double> &Ql_Raw,std::vector<double> &Qs_Raw,
//                                 std::vector<double> &Ql,std::vector<ULong64_t> &TDC,std::vector<double> &TDC_ns,
//                                 std::vector<ULong64_t> &QDC_c,std::vector<double> &QDC_f,std::vector<int> &det_ids_QDC,
//                                 std::vector<int> &det_ids_TDC){
//
//     FATIMA_Data.SetDATA(QDC_FIRED,TDC_FIRED,Ql_Raw,Qs_Raw,Ql,TDC,TDC_ns,QDC_c,QDC_f,det_ids_QDC,det_ids_TDC);
//  cout << FATIMA_Data.SetDATA(&Ql) << endl;
//     Event_Type = 3;
//
// }



int Raw_Event::get_Event_type(){
    return Event_Type;

}


// void Raw_Event::set_AIDA_Event(){}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------AIDA  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
 void Raw_Event::set_DATA_AIDA(double* Aida_E, int* Aida_feeID, int* Aida_ChID, ULong64_t* Aida_WR, int Aida_hits, bool* Aida_veto, int* Aida_side, int* Aida_strip, int* Aida_EvtID, ULong64_t* Aida_Fast, int* Aida_ADC){

    this->  AIDA_Hits = Aida_hits;
    for (int i =0; i< Aida_hits; i++){
        this->  AIDA_Energy[i] = Aida_E[i];
        this->  AIDA_FEE[i] = Aida_feeID[i];
        this->  AIDA_CHA_ID[i] = Aida_ChID[i];
        this->  AIDA_WR[i] = Aida_WR[i];

        this->  AIDA_HighE_VETO[i] = Aida_veto[i];
        this->  AIDA_SIDE[i] = Aida_side[i];
        this->  AIDA_STRIP[i] =  Aida_strip[i];
        this->  AIDA_EVT_ID[i] = Aida_EvtID[i];

        this-> AIDA_FastTime[i] = Aida_Fast[i];

        this->AIDA_ADC[i] = Aida_ADC[i];

       }
   Event_Type = 1;
 }



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- FINGER TAMEX  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Raw_Event::set_DATA_FINGER(int* it,double** Edge_Coarse,double** Edge_fine,UInt** ch_ed,double* Coarse_Trigger,double* Fine_Trigger,int amount_hit_tamex, int** Lead_Arr){

    this->amount_hit_tamex = amount_hit_tamex;
    //reset lead and trail hits
    for(int i = 0;i < amount_hit_tamex;i++){
        for(int j = 0;j < 32;j++){
            leading_hits_ch[i][j] = 0;
            trailing_hits_ch[i][j] = 0;
            leading_array[i][j] = 0;
        }
    }

    //loop over all 4 tamex modules
    for(int i = 0;i < amount_hit_tamex;i++){
        iterator[i] = it[i];
        trigger_coarse[i] = Coarse_Trigger[i];
        trigger_fine[i] = Fine_Trigger[i];
        //cout << "trigger_coarse[i] " << trigger_coarse[i]<< " trigger_fine[i] "<<trigger_fine[i]<< endl;
        fired_tamex[i] = (iterator[i] > 0);
        leading_hits[i] = 0;
        trailing_hits[i] = 0;
          
        for(int j = 0;j < iterator[i];j++){
            ch_ID[i][j] = ch_ed[i][j];
            leading_array[i][j] = Lead_Arr[i][j];
             if(ch_ID[i][j] % 2 == 1){
               
                coarse_T_edge_lead[i][j] = (double) Edge_Coarse[i][j];
                fine_T_edge_lead[i][j] = (double) Edge_fine[i][j];
                phys_channel[i][j] = (ch_ID[i][j]+1)/2;
                leading_hits[i]++;
                leading_hits_ch[i][phys_channel[i][j]]++;
               // cout<<"RAW " << " phys_channel[i][j] " << phys_channel[i][j] << " coarse_T_edge_lead[i][j] " <<coarse_T_edge_lead[i][j] << " fine_T_edge_lead[i][j] " <<fine_T_edge_lead[i][j] << " leading_hits[i] " <<leading_hits[i] << " i " << i << " j " << j <<endl;
      
            }
            else{
                coarse_T_edge_trail[i][j] = (double)  Edge_Coarse[i][j];
                fine_T_edge_trail[i][j] =(double)  Edge_fine[i][j];
            
                trailing_hits[i]++;
                phys_channel[i][j] = (ch_ID[i][j])/2;
                trailing_hits_ch[i][phys_channel[i][j]]++;
             //   cout<<"RAW " << " phys_channel[i][j] " << phys_channel[i][j] << " coarse_T_edge_trail[i][j] " <<coarse_T_edge_trail[i][j] << " fine_T_edge_trail[i][j] " <<fine_T_edge_trail[i][j] << " trailing_hits[i] " <<trailing_hits[i]  << " i " << i << " j " << j << endl;

            }
        }
    }

    Event_Type = 2;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//---------------------------------------PLASTIC VME  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Raw_Event::set_DATA_PLASTIC_VME(int QDC_ite, int TDC_ite, double* VME_QDC_Dat1, double* VME_QDC_Dat2, int* VME_QDC_Chan, double* VME_TDC_Dat, int* VME_TDC_Chan){
     QDC_IT=0;
    this-> QDC_IT = QDC_ite;
    this->TDC_IT = TDC_ite;
    for (int i = 0; i<QDC_IT; i++){
    this->VME_QDC_DAT1[i] = VME_QDC_Dat1[i];
    this->VME_QDC_DAT2[i] = VME_QDC_Dat2[i];
    this->VME_QDC_CHA[i] = VME_QDC_Chan[i];
   
      }

    for(int j=0; j<TDC_IT; j++) {

    this->VME_TDC_CHA[j] = VME_TDC_Chan[j];
    this->VME_TDC_DAT[j] = VME_TDC_Dat[j];


         }
       Event_Type=2;
        }
/*
    void Raw_Event::set_DATA_SCALAR(int Scalar_ite,  double* Scalar_Dat){

    this->SCALAR_ITERATOR = Scalar_ite;

    for(int i=1; i<SCALAR_ITERATOR; i++) {
    this->SCALAR_DATA[i] = Scalar_Dat[i];
    

         }*/
/*
    Event_Type=2;
        }*/
//---------------------------------------------------------------

bool Raw_Event::PLASTIC_CheckVME(){
    return VME_Event;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------PLASTIC TAMEX  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Raw_Event::set_DATA_PLASTIC_TAMEX(int* it_bPlas,double** Edge_Coarse_bPlas,double** Edge_fine_bPlas,UInt** ch_ed_bPlas,double* Coarse_Trigger_bPlas,double* Fine_Trigger_bPlas,int amount_hit_tamex_bPlas,int** Lead_Arr_bPlas, int* tamex_id){

    this->amount_hit_tamex_bPlas = amount_hit_tamex_bPlas;
    //reset lead and trail hits
    for(int i = 0;i < amount_hit_tamex_bPlas;++i){
        for(int j = 0;j < 48;j++){
            leading_hits_ch_bPlas[i][j] = 0;
            trailing_hits_ch_bPlas[i][j] = 0;
        }
    }

    //loop over all 4 tamex modules
    for(int i = 0;i < amount_hit_tamex_bPlas;i++){
        iterator_bPlas[i] = it_bPlas[i];
        tamex_id_bPlas[i] = tamex_id[i];
       
        trigger_coarse_bPlas[i] = Coarse_Trigger_bPlas[i];
        trigger_fine_bPlas[i] = Fine_Trigger_bPlas[i];
        fired_tamex_bPlas[i] = (iterator_bPlas[i] > 0);
        leading_hits_bPlas[i] = 0;
        trailing_hits_bPlas[i] = 0;

        for(int j = 0;j < iterator_bPlas[i];j++){
            ch_ID_bPlas[i][j] = ch_ed_bPlas[i][j];
            if(ch_ID_bPlas[i][j] % 2 == 1){

                coarse_T_edge_lead_bPlas[i][j] = (double) Edge_Coarse_bPlas[i][j];
                fine_T_edge_lead_bPlas[i][j] = (double) Edge_fine_bPlas[i][j];

                leading_array_bPlas[i][j] = Lead_Arr_bPlas[i][j];
                phys_channel_bPlas[i][j] =(tamex_id_bPlas[i]*16)+((ch_ID_bPlas[i][j]+1)/2-1);
                
                lead_T_bPlas[i][j] = (coarse_T_edge_lead_bPlas[i][j]-fine_T_edge_lead_bPlas[i][j]);
                
             
                leading_hits_bPlas[i]++;
                leading_hits_ch_bPlas[i][phys_channel_bPlas[i][j]]++;
                 }
            else{

                coarse_T_edge_trail_bPlas[i][j] = (double)  Edge_Coarse_bPlas[i][j];
                fine_T_edge_trail_bPlas[i][j] =(double)  Edge_fine_bPlas[i][j];

                trailing_hits_bPlas[i]++;
                phys_channel_bPlas[i][j] = (tamex_id_bPlas[i]*16)+(ch_ID_bPlas[i][j])/2-1;
                
                trail_T_bPlas[i][j] = (coarse_T_edge_trail_bPlas[i][j]-fine_T_edge_trail_bPlas[i][j]);
                
                trailing_hits_ch_bPlas[i][phys_channel_bPlas[i][j]]++;
          }
    
        }
    }

    Event_Type = 2;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------FATIMA VME  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Raw_Event::set_DATA_FATIMA(int QDC_FIRED,int TDC_FIRED,
                                double* Ql_Raw,double* Qs_Raw,
                                double* Ql,
                                ULong64_t* TDC, double* TDC_ns,
                                ULong64_t* QDC_c, double* QDC_f,
                                int* det_ids_QDC,int* det_ids_TDC){

    this->FAT_QDCs_FIRED = QDC_FIRED;
    this->FAT_TDCs_FIRED = TDC_FIRED;
    
    int dets_fired = 0;
    int tdc_dets_fired =0;
    for (int i=0; i<QDC_FIRED; i++) {
        this->FAT_QDC_id[i] = det_ids_QDC[i];
        //if (det_ids_QDC[i] == 35) cout<<"I am in QDC"<<endl;
        this->FAT_QLong[i]  = Ql[i];
        this->FAT_QLong_Raw[i]  = Ql_Raw[i];
        this->FAT_QShort_Raw[i] = Qs_Raw[i];
        this->FAT_QDC_t_coarse[i] = QDC_c[i];
        this->FAT_QDC_t_fine[i] = QDC_f[i];
      
        for (int j=0; j<TDC_FIRED; j++) {
            if (det_ids_QDC[i] == det_ids_TDC[j]) {
                this->FAT_id[dets_fired] = det_ids_QDC[i];
            //    if (det_ids_TDC[j] == 35) cout<<"I am in TDC"<<endl;

                this->FAT_E[dets_fired] = Ql[i];
                this->FAT_ratio[dets_fired] = (double) Qs_Raw[i]/Ql_Raw[i];
                this->FAT_t[dets_fired] = TDC_ns[j];
                this->FAT_t_qdc[dets_fired] = QDC_c[i];
                dets_fired++;
            }
        }
    }
    this->FAT_DET_FIRED = dets_fired;

    for (int i=0; i<TDC_FIRED; i++) {
        this->FAT_TDC_id[i]        = det_ids_TDC[i];
        this->FAT_TDC_timestamp[i] = TDC[i];
        this->FAT_TDC_timestamp_raw[i] = TDC[i];
     tdc_dets_fired++;

    }
       Event_Type = 3;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- SCALERS  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Raw_Event::set_DATA_SCALAR(int Scalar_ite,  double* Scalar_Dat,  bool Scaler_Used){
// 
     
      if(Scaler_Used==true){
          
        this->SCALAR_ITERATOR = Scalar_ite;
        
        for(int i=1; i<17; i++) {
        this->SCALAR_DATA[i] = Scalar_Dat[i];
         }

    Event_Type=2;
        }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- FATIMA TAMEX  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Raw_Event::set_DATA_FATIMA_TAMEX(int* it_fat,double** Edge_Coarse_fat,double** Edge_fine_fat,UInt** ch_ed_fat,double* Coarse_Trigger_fat,double* Fine_Trigger_fat,int amount_hit_tamex_fat, int** Lead_Arr_fat){

    this->amount_hit_tamex_fat = amount_hit_tamex_fat;
    //reset lead and trail hits
    for(int i = 0;i < amount_hit_tamex;i++){
        iterator_fat[i] = 0;
        for(int j = 0;j < 32;j++){
            leading_hits_ch_fat[i][j] = 0;
            trailing_hits_ch_fat[i][j] = 0;
            leading_array_fat[i][j] = 0;
//             coarse_T_edge_lead_fat[i][j] = 0;
//             fine_T_edge_lead_fat[i][j] = 0;
        }
    }

    //loop over all 4 tamex modules
    for(int i = 0;i < amount_hit_tamex_fat;i++){
        iterator_fat[i] = it_fat[i];
        trigger_coarse_fat[i] = Coarse_Trigger_fat[i];
        trigger_fine_fat[i] = Fine_Trigger_fat[i];
        fired_tamex_fat[i] = (iterator_fat[i] > 0);
        leading_hits_fat[i] = 0;
        trailing_hits_fat[i] = 0;

        for(int j = 0;j < iterator_fat[i];++j){
            ch_ID_fat[i][j] = ch_ed_fat[i][j];
            leading_array_fat[i][j] = Lead_Arr_fat[i][j];
            if(ch_ID_fat[i][j] % 2 == 1){

                coarse_T_edge_lead_fat[i][j] = (double) Edge_Coarse_fat[i][j];
                fine_T_edge_lead_fat[i][j] = (double) Edge_fine_fat[i][j];

                phys_channel_fat[i][j] = (ch_ID_fat[i][j]+1)/2 -1;
                //cout <<"RAW Ch Lead " << ch_ID_fat[i][j] <<" phys_channel_fat[i][j] " << phys_channel_fat[i][j]<<endl;
                leading_hits_fat[i]++;
                leading_hits_ch_fat[i][phys_channel_fat[i][j]]++;
                  // cout <<"RAW Lead" << "coarse_T_edge_lead_fat[i][j] " << coarse_T_edge_lead_fat[i][j] << " fine_T_edge_lead_fat[i][j] " <<fine_T_edge_lead_fat[i][j] <<" phys_channel_fat[i][j] " <<phys_channel_fat[i][j]<< " i " << i << " j " << j << endl;
                 }
            else{

                coarse_T_edge_trail_fat[i][j] = (double)  Edge_Coarse_fat[i][j];
                fine_T_edge_trail_fat[i][j] =(double)  Edge_fine_fat[i][j];

                trailing_hits_fat[i]++;
                phys_channel_fat[i][j] = (ch_ID_fat[i][j])/2-1;
                trailing_hits_ch_fat[i][phys_channel_fat[i][j]]++;
                //cout <<"RAW Ch Trail " << ch_ID_fat[i][j] <<" phys_channel_fat[i][j] " << phys_channel_fat[i][j] << endl;
               //cout <<"RAW Trail" << "coarse_T_edge_trail_fat[i][j] " << coarse_T_edge_trail_fat[i][j] << " fine_T_edge_trail_fat[i][j] " <<fine_T_edge_trail_fat[i][j] << " i " << i << " j " << j << endl;
          }
        }
    }

    Event_Type = 3;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- GALILEO FEBEX  ------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Raw_Event::set_DATA_GALILEO(int GAL_FIRED,ULong64_t* sum_time,int* hit_pattern,ULong64_t* chan_time,double* chan_en, int* det_ids, int* crystal_ids, bool* pileup, bool* overflow){
        this->GAL_FIRED = GAL_FIRED;
    
        for(int i = 0;i < GAL_FIRED;++i){
    
        GALILEO_Det_Nums[i] = det_ids[i];
        GALILEO_Crystal_Nums[i] = crystal_ids[i];
        GALILEO_sum_time[i] = sum_time[i];
        GALILEO_hit_pattern[i] = hit_pattern[i];
        GALILEO_chan_time[i] = chan_time[i];
        GALILEO_chan_energy[i] = chan_en[i];
        GALILEO_Pileup[i] = pileup[i];
        GALILEO_Overflow[i] = overflow[i];

    }

    Event_Type = 4;
}
//TEMPORARY GETTERS FOR FRS, FATIMA, PLASTIC, and GALILEO

// #############################################################

//------------------------------------------ FRS ------------------------------------------//

Float_t Raw_Event::get_FRS_MusicdE(int i){return MUSIC_dE[i];}
Float_t Raw_Event::get_FRS_MusicdE_corr(int i){return MUSIC_dE_cor[i];}
Int_t Raw_Event::get_FRS_MusicE1(int i){return MUSIC_e1[i];}
Int_t Raw_Event::get_FRS_MusicE2(int i){return MUSIC_e2[i];}
Int_t Raw_Event::get_FRS_MusicT1(int i){return MUSIC_t1[i];}
Int_t Raw_Event::get_FRS_MusicT2(int i){return MUSIC_t2[i];}

Float_t Raw_Event::get_FRS_sci_l(int i){return sci_l[i];}
Float_t Raw_Event::get_FRS_sci_r(int i){return sci_r[i];}
Float_t Raw_Event::get_FRS_sci_e(int i){return sci_e[i];}
Float_t Raw_Event::get_FRS_sci_tx(int i){return sci_tx[i];}
Float_t Raw_Event::get_FRS_sci_x(int i){return sci_x[i];}


Int_t Raw_Event::get_FRS_dt_21l_21r(){return dt_21l_21r;}
Int_t Raw_Event::get_FRS_dt_41l_41r(){return dt_41l_41r;}
Int_t Raw_Event::get_FRS_dt_21l_41l(){return dt_21l_41l;}
Int_t Raw_Event::get_FRS_dt_21r_41r(){return dt_21r_41r;}
Int_t Raw_Event::get_FRS_dt_42l_42r(){return dt_42l_42r;}
Int_t Raw_Event::get_FRS_dt_43l_43r(){return dt_43l_43r;}
Int_t Raw_Event::get_FRS_dt_21l_42l(){return dt_21l_42l;}
Int_t Raw_Event::get_FRS_dt_21r_42r(){return dt_21r_42r;}
Int_t Raw_Event::get_FRS_dt_81l_81r(){return dt_81l_81r;}
Int_t Raw_Event::get_FRS_dt_21l_81l(){return dt_21l_81l;}
Int_t Raw_Event::get_FRS_dt_21r_81r(){return dt_21r_81r;}

Float_t Raw_Event::get_FRS_tofll2(){return sci_tofll2;}
Float_t Raw_Event::get_FRS_tofll3(){return sci_tofll3;}
Float_t Raw_Event::get_FRS_tof2(){return sci_tof2;}
Float_t Raw_Event::get_FRS_tofrr2(){return sci_tofrr2;}
Float_t Raw_Event::get_FRS_tofrr3(){return sci_tofrr3;}
Float_t Raw_Event::get_FRS_tof3(){return sci_tof3;}

Float_t Raw_Event::get_FRS_x2(){return ID_x2;}
Float_t Raw_Event::get_FRS_y2(){return ID_y2;}
Float_t Raw_Event::get_FRS_a2(){return ID_a2;}
Float_t Raw_Event::get_FRS_b2(){return ID_b2;}

Float_t Raw_Event::get_FRS_x4(){return ID_x4;}
Float_t Raw_Event::get_FRS_y4(){return ID_y4;}
Float_t Raw_Event::get_FRS_a4(){return ID_a4;}
Float_t Raw_Event::get_FRS_b4(){return ID_b4;}

Float_t Raw_Event::get_FRS_tpcX(int i){return TPC_x[i];}
Float_t Raw_Event::get_FRS_tpcY(int i){return TPC_y[i];}
Int_t   Raw_Event::get_FRS_tpclt(int i,int j){return TPC_lt[i][j];}
Int_t   Raw_Event::get_FRS_tpcrt(int i,int j){return TPC_rt[i][j];}
Float_t Raw_Event::get_FRS_tpcx0(){return TPC_x0;}
Float_t Raw_Event::get_FRS_tpcx1(){return TPC_x1;}

Float_t Raw_Event::get_FRS_brho(int i){return ID_brho[i];}
Float_t Raw_Event::get_FRS_rho(int i){return ID_rho[i];}

Float_t Raw_Event::get_FRS_beta(){return beta;}
Float_t Raw_Event::get_FRS_beta3(){return beta3;}
Float_t Raw_Event::get_FRS_gamma(){return gamma;}

Float_t Raw_Event::get_FRS_AoQ(){return AoQ;}
Float_t Raw_Event::get_FRS_AoQ_corr(){return AoQ_corr;}
Float_t Raw_Event::get_FRS_z(){return z;}
Float_t Raw_Event::get_FRS_z2(){return z2;}
Float_t Raw_Event::get_FRS_z3(){return z3;}
Float_t Raw_Event::get_FRS_dEdeg(){return dEdeg;}
Float_t Raw_Event::get_FRS_dEdegoQ(){return dEdegoQ;}

Float_t Raw_Event::get_FRS_timestamp(){return timestamp;}
Float_t Raw_Event::get_FRS_ts(){return ts;}
Float_t Raw_Event::get_FRS_ts2(){return ts2;}

Int_t Raw_Event::get_FRS_time_in_ms(){return time_in_ms;}
Int_t Raw_Event::get_FRS_spill_count(){return spill_count;}
Int_t Raw_Event::get_FRS_ibin_for_s(){return ibin_for_s;}
Int_t Raw_Event::get_FRS_ibin_for_100ms(){return ibin_for_100ms;}
Int_t Raw_Event::get_FRS_ibin_for_spill(){return ibin_for_spill;}
Int_t Raw_Event::get_FRS_ibin_clean_for_s(){return ibin_clean_for_s;}
Int_t Raw_Event::get_FRS_ibin_clean_for_100ms(){return ibin_clean_for_100ms;}
Int_t Raw_Event::get_FRS_ibin_clean_for_spill(){return ibin_clean_for_spill;}
UInt_t* Raw_Event::get_FRS_increase_scaler_temp(){return increase_scaler_temp;}

Double_t Raw_Event::get_FRS_TRaw_vftx_21l(){return TRaw_vftx_21L;}
Double_t Raw_Event::get_FRS_TRaw_vftx_21r(){return TRaw_vftx_21R;}
Double_t Raw_Event::get_FRS_TRaw_vftx_22l(){return TRaw_vftx_22L;}
Double_t Raw_Event::get_FRS_TRaw_vftx_22r(){return TRaw_vftx_22R;}
Double_t Raw_Event::get_FRS_TRaw_vftx_41l(){return TRaw_vftx_41L;}
Double_t Raw_Event::get_FRS_TRaw_vftx_41r(){return TRaw_vftx_41R;}
Double_t Raw_Event::get_FRS_TRaw_vftx_42l(){return TRaw_vftx_42L;}
Double_t Raw_Event::get_FRS_TRaw_vftx_42r(){return TRaw_vftx_42R;}
Double_t Raw_Event::get_FRS_TRaw_vftx(int i){return TRaw_vftx[i];}


// #######################################################

//------------------------------------------ WHITE RABBIT ------------------------------------------//

void Raw_Event::set_WR(ULong64_t WR){this->WR = WR;}

ULong64_t Raw_Event::get_WR(){return WR;}


//------------------------------------------ AIDA ------------------------------------------//

     double Raw_Event::get_AIDA_Energy(int i){return AIDA_Energy[i];}
     int Raw_Event::get_AIDA_FEE_ID(int i){return AIDA_FEE[i];}
     int Raw_Event::get_AIDA_CHA_ID(int i){return AIDA_CHA_ID[i];}
     int Raw_Event::get_AIDA_HITS(){return AIDA_Hits;}
     ULong64_t Raw_Event::get_AIDA_WR(int i){return AIDA_WR[i];}
     bool Raw_Event::get_AIDA_HighE_VETO(int i){return AIDA_HighE_VETO[i];}
     int Raw_Event::get_AIDA_SIDE(int i){return AIDA_SIDE[i];}
     int Raw_Event::get_AIDA_STRIP(int i){return AIDA_STRIP[i];}
     int Raw_Event::get_AIDA_EVTID(int i){return AIDA_EVT_ID[i];}
     ULong64_t Raw_Event::get_AIDA_FastTime(int i){return AIDA_FastTime[i];}
     int Raw_Event::get_AIDA_ADC(int i){return AIDA_ADC[i];}



//------------------------------------ FATIMA VME ------------------------------------------//

      int Raw_Event::get_FAT_det_fired(){return FAT_DET_FIRED;}
      int Raw_Event::get_FAT_id(int i){return FAT_id[i];}
   double Raw_Event::get_FAT_E(int i){return FAT_E[i];}
   double Raw_Event::get_FAT_ratio(int i){return FAT_ratio[i];}
   double Raw_Event::get_FAT_t(int i){return FAT_t[i];}
   double Raw_Event::get_FAT_t_qdc(int i){return FAT_t_qdc[i];}

      int Raw_Event::get_FAT_QDCs_fired(){return FAT_QDCs_FIRED;}
      int Raw_Event::get_FAT_QDC_id(int i){return FAT_QDC_id[i];}
   double Raw_Event::get_FAT_QLong(int i){return FAT_QLong[i];}
   double Raw_Event::get_FAT_QShort_Raw(int i){return FAT_QShort_Raw[i];}
   double Raw_Event::get_FAT_QLong_Raw(int i){return FAT_QLong_Raw[i];}
   ULong64_t Raw_Event::get_FAT_QDC_t_Coarse(int i){return FAT_QDC_t_coarse[i];}
   double Raw_Event::get_FAT_QDC_t_Fine(int i){return FAT_QDC_t_fine[i];}

      int Raw_Event::get_FAT_TDCs_fired(){return FAT_TDCs_FIRED;}
      int Raw_Event::get_FAT_TDC_id(int i){return FAT_TDC_id[i];}
   double Raw_Event::get_FAT_TDC_timestamp(int i){return FAT_TDC_timestamp[i];}
   double Raw_Event::get_FAT_TDC_timestamp_raw(int i){return FAT_TDC_timestamp_raw[i];}

//FATIMA_DataStruct* Raw_Event::PassFATIMA(){return &FATIMA_Data;}

//------------------------------------- FATIMA TAMEX ------------------------------------------//

    int Raw_Event::get_FATIMA_tamex_hits(){return amount_hit_tamex_fat;}

    int Raw_Event::get_FATIMA_am_Fired(int i){return iterator_fat[i];}

    double Raw_Event::get_FATIMA_trigger_T(int i){return (trigger_coarse_fat[i] - trigger_fine_fat[i])*5000;}

    int Raw_Event::get_FATIMA_CH_ID(int i,int j){return ch_ID_fat[i][j];}

    double Raw_Event::get_FATIMA_lead_T(int i,int j){
        //cout << "SEND l" << coarse_T_edge_lead[i][j] << " " << fine_T_edge_lead[i][j]  << " " <<  coarse_T_edge_lead[i][j]*5 - fine_T_edge_lead[i][j] << endl;
  
        return (coarse_T_edge_lead_fat[i][j] - fine_T_edge_lead_fat[i][j]);
        
                    }

    double Raw_Event::get_FATIMA_coarse_lead(int i,int j){
       // cout <<"RAW getter coarse_T_edge_lead_fat[i][0]" <<coarse_T_edge_lead_fat[i][0]<<" i " <<i << " j " << j <<endl; 
        return coarse_T_edge_lead_fat[i][j];
   }
                    
   double Raw_Event::get_FATIMA_fine_lead(int i,int j){
        return fine_T_edge_lead_fat[i][j];
   }
   double Raw_Event::get_FATIMA_coarse_trail(int i,int j){
        return coarse_T_edge_trail_fat[i][j];
   }
   double Raw_Event::get_FATIMA_fine_trail(int i,int j){
        return fine_T_edge_trail_fat[i][j];
   }
   
    double Raw_Event::get_FATIMA_trail_T(int i,int j){
        //cout << "SEND t" << coarse_T_edge_trail[i][j] << " " << fine_T_edge_trail[i][j] << endl;
        return (coarse_T_edge_trail_fat[i][j] - fine_T_edge_trail_fat[i][j]); 
                    }
                    
    double Raw_Event::get_FATIMA_Lead_Lead(int i,int j){
        double T_lead_fat1 = (coarse_T_edge_lead_fat[i][j] - fine_T_edge_lead_fat[i][j]);
        double T_lead_fat2 = (coarse_T_edge_lead_fat[i][j+2] - fine_T_edge_lead_fat[i][j+2]);
               return T_lead_fat1 - T_lead_fat2;
                   }

    double Raw_Event::get_FATIMA_TOT(int i,int j){
        double T_lead_fat = (coarse_T_edge_lead_fat[i][j] - fine_T_edge_lead_fat[i][j]);
        double T_trail_fat = (coarse_T_edge_trail_fat[i][j+1] - fine_T_edge_trail_fat[i][j+1]);
               return T_trail_fat - T_lead_fat;
                   }
                   
//          double Raw_Event::get_FATIMA_TOT_added(int i,int j){
//         double T_lead1 = (coarse_T_edge_lead_fat[i][j] - fine_T_edge_lead_fat[i][j])*5000;
//         double T_trail1 = (coarse_T_edge_trail_fat[i][j+1] - fine_T_edge_trail_fat[i][j+1])*5000;
//         double T_lead2 = (coarse_T_edge_lead_fat[i][j-2] - fine_T_edge_lead_fat[i][j-2])*5000;
//         double T_trail2 = (coarse_T_edge_trail_fat[i][j-1] - fine_T_edge_trail_fat[i][j-1])*5000;
//                return (T_trail1 - T_lead1) + (T_trail2 - T_lead2);
//                    }

    int Raw_Event::get_FATIMA_trail_hits(int i){return trailing_hits_fat[i];}

    int Raw_Event::get_FATIMA_lead_hits(int i){return leading_hits_fat[i];}

    int Raw_Event::get_FATIMA_physical_channel(int i,int j){return phys_channel_fat[i][j];}

    int Raw_Event::get_FATIMA_physical_lead_hits(int i,int j){return leading_hits_ch_fat[i][j];}

    int Raw_Event::get_FATIMA_physical_trail_hits(int i,int j){return trailing_hits_ch_fat[i][j];}


//------------------------------------------ FINGER TAMEX ------------------------------------------//

    int Raw_Event::get_FINGER_tamex_hits(){return amount_hit_tamex;}

    int Raw_Event::get_FINGER_am_Fired(int i){return iterator[i];}

    double Raw_Event::get_FINGER_trigger_T(int i){return (trigger_coarse[i] - trigger_fine[i])*5000;}

    int Raw_Event::get_FINGER_CH_ID(int i,int j){return ch_ID[i][j];}

    double Raw_Event::get_FINGER_lead_T(int i,int j){
        //cout << "SEND l" << coarse_T_edge_lead[i][j] << " " << fine_T_edge_lead[i][j]  << " " <<  coarse_T_edge_lead[i][j]*5 - fine_T_edge_lead[i][j] << endl;
        return (coarse_T_edge_lead[i][j] - fine_T_edge_lead[i][j])*5000;
                    }

    double Raw_Event::get_FINGER_coarse_lead(int i,int j){

        return coarse_T_edge_lead[i][j]*5000;
   }
                    
   double Raw_Event::get_FINGER_fine_lead(int i,int j){
        return fine_T_edge_lead[i][j]*5000;
   }
   double Raw_Event::get_FINGER_coarse_trail(int i,int j){
        return coarse_T_edge_trail[i][j]*5000;
   }
   double Raw_Event::get_FINGER_fine_trail(int i,int j){
        return fine_T_edge_trail[i][j]*5000;
   }
   
    double Raw_Event::get_FINGER_trail_T(int i,int j){
        //cout << "SEND t" << coarse_T_edge_trail[i][j] << " " << fine_T_edge_trail[i][j] << endl;
        return (coarse_T_edge_trail[i][j] - fine_T_edge_trail[i][j])*5000; 
                    }

    double Raw_Event::get_FINGER_TOT(int i,int j){
        double T_lead = (coarse_T_edge_lead[i][j] - fine_T_edge_lead[i][j])*5000;
        double T_trail = (coarse_T_edge_trail[i][j+1] - fine_T_edge_trail[i][j+1])*5000;
               return T_trail - T_lead;
                   }
                   
         double Raw_Event::get_FINGER_TOT_added(int i,int j){
        double T_lead1 = (coarse_T_edge_lead[i][j] - fine_T_edge_lead[i][j])*5000;
        double T_trail1 = (coarse_T_edge_trail[i][j+1] - fine_T_edge_trail[i][j+1])*5000;
        double T_lead2 = (coarse_T_edge_lead[i][j-2] - fine_T_edge_lead[i][j-2])*5000;
        double T_trail2 = (coarse_T_edge_trail[i][j-1] - fine_T_edge_trail[i][j-1])*5000;
               return (T_trail1 - T_lead1) + (T_trail2 - T_lead2);
                   }

    int Raw_Event::get_FINGER_trail_hits(int i){return trailing_hits[i];}

    int Raw_Event::get_FINGER_lead_hits(int i){return leading_hits[i];}

    int Raw_Event::get_FINGER_physical_channel(int i,int j){return phys_channel[i][j];}

    int Raw_Event::get_FINGER_physical_lead_hits(int i,int j){return leading_hits_ch[i][j];}

    int Raw_Event::get_FINGER_physical_trail_hits(int i,int j){return trailing_hits_ch[i][j];}

//------------------------------------- bPlastic TAMEX ------------------------------------------//

    int Raw_Event::get_PLASTIC_tamex_hits(){return amount_hit_tamex_bPlas;}
    
    int Raw_Event::get_PLASTIC_TAMEX_ID(int i){return tamex_id_bPlas[i];}
    
    int Raw_Event::get_PLASTIC_am_Fired(int i){return iterator_bPlas[i];}

    double Raw_Event::get_PLASTIC_trigger_T(int i){return (trigger_coarse_bPlas[i] - trigger_fine_bPlas[i])*5000;}

    int Raw_Event::get_PLASTIC_CH_ID(int i,int j){return ch_ID_bPlas[i][j];}

    double Raw_Event::get_PLASTIC_lead_T(int i,int j){
       
        return (lead_T_bPlas[i][j]);
        
                    }

    double Raw_Event::get_PLASTIC_coarse_lead(int i,int j){

        return coarse_T_edge_lead_bPlas[i][j];
   }
                    
   double Raw_Event::get_PLASTIC_fine_lead(int i,int j){
        return fine_T_edge_lead_bPlas[i][j];
   }
   double Raw_Event::get_PLASTIC_coarse_trail(int i,int j){
        return coarse_T_edge_trail_bPlas[i][j];
   }
   double Raw_Event::get_PLASTIC_fine_trail(int i,int j){
        return fine_T_edge_trail_bPlas[i][j];
   }
   
    double Raw_Event::get_PLASTIC_trail_T(int i,int j){
        
        return (trail_T_bPlas[i][j]); 
                    }

    double Raw_Event::get_PLASTIC_TOT(int i,int j){
        double T_lead_bPlas = (coarse_T_edge_lead_bPlas[i][j] - fine_T_edge_lead_bPlas[i][j]);
        double T_trail_bPlas = (coarse_T_edge_trail_bPlas[i][j+1] - fine_T_edge_trail_bPlas[i][j+1]);
               return T_trail_bPlas - T_lead_bPlas;
                   }
                   
//          double Raw_Event::get_PLASTIC_TOT_added(int i,int j){
//         double T_lead1 = (coarse_T_edge_lead_bPlas[i][j] - fine_T_edge_lead_bPlas[i][j])*5000;
//         double T_trail1 = (coarse_T_edge_trail_bPlas[i][j+1] - fine_T_edge_trail_bPlas[i][j+1])*5000;
//         double T_lead2 = (coarse_T_edge_lead_bPlas[i][j-2] - fine_T_edge_lead_bPlas[i][j-2])*5000;
//         double T_trail2 = (coarse_T_edge_trail_bPlas[i][j-1] - fine_T_edge_trail_bPlas[i][j-1])*5000;
//                return (T_trail1 - T_lead1) + (T_trail2 - T_lead2);
//                    }

    int Raw_Event::get_PLASTIC_trail_hits(int i){return trailing_hits_bPlas[i];}

    int Raw_Event::get_PLASTIC_lead_hits(int i){return leading_hits_bPlas[i];}

    int Raw_Event::get_PLASTIC_physical_channel(int i,int j){return phys_channel_bPlas[i][j];}

    int Raw_Event::get_PLASTIC_physical_lead_hits(int i,int j){return leading_hits_ch_bPlas[i][j];}

    int Raw_Event::get_PLASTIC_physical_trail_hits(int i,int j){return trailing_hits_ch_bPlas[i][j];}
    
    
//------------------------------------------ bPlastic VME----------------------------------------//
 int Raw_Event::get_plastic_VME_QDC_fired(){return QDC_IT;}
 int Raw_Event::get_plastic_VME_TDC_fired(){return TDC_IT;}
 double Raw_Event::get_plastic_VME_QDC_dat1(int i){return VME_QDC_DAT1[i];}
 double Raw_Event::get_plastic_VME_QDC_dat2(int i){return VME_QDC_DAT2[i];}
 int Raw_Event::get_plastic_VME_QDC_cha(int i){return VME_QDC_CHA[i];}
 double Raw_Event::get_plastic_VME_TDC_dat(int i){return VME_TDC_DAT[i];}
 int Raw_Event::get_plastic_VME_TDC_cha(int i){return VME_TDC_CHA[i];}

 int Raw_Event::get_scalar_iterator(){return SCALAR_ITERATOR;}
 double Raw_Event::get_scalar_data(int i){return SCALAR_DATA[i];}
/*
PLASTIC_VME_DataStruct* Raw_Event::PassPLASTIC_VME(){ return &PLASTIC_VME_Data;}

//---------------------------------------------------------------

PLASTIC_DataStruct* Raw_Event::PassPLASTIC(){ return &PLASTIC_Data;}*/

//------------------------------------------ GALILEO FEBEX ------------------------------------------//


    int Raw_Event::get_GALILEO_am_Fired(){return GAL_FIRED;}
    
    ULong64_t Raw_Event::get_GALILEO_Event_T(int i){return GALILEO_sum_time[i];}
      
    
    
    bool Raw_Event::get_GALILEO_Pileup(int i){return GALILEO_Pileup[i];}
    
    bool Raw_Event::get_GALILEO_Overflow(int i){return GALILEO_Overflow[i];}
    
    int Raw_Event::get_GALILEO_Hit_Pattern(int i){return GALILEO_hit_pattern[i];}
    
    ULong64_t Raw_Event::get_GALILEO_Chan_T(int i){return GALILEO_chan_time[i];}
    
    double Raw_Event::get_GALILEO_Chan_E(int i){return GALILEO_chan_energy[i];}
    
    int Raw_Event::get_GALILEO_Det_id(int i){return GALILEO_Det_Nums[i];}   
    
    int Raw_Event::get_GALILEO_Crystal_id(int i){return GALILEO_Crystal_Nums[i];}
