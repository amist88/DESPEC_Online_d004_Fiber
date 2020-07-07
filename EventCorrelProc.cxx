// $Id: EventCorrelProc.cxx 754 2011-05-18 11:04:52Z adamczew $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "EventCorrelProc.h"

#include <cstdlib>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "TGo4WinCond.h"
#include "TGo4CondArray.h"
#include "TGo4Analysis.h"

#include "TGo4Picture.h"

#include "EventCorrelStore.h"
#include "EventUnpackStore.h"
#include "EventAnlStore.h"


#include "TAidaConfiguration.h"



//-----------------------------------------------------------
EventCorrelProc::EventCorrelProc() :
  TGo4EventProcessor()/*,
  fParam1(0),
  fTimeDiff(0),
  fGatedHist(0),
  fCoincQ1A1(0),
  fCoincQ1T1(0),
  fconHis1(0)*/
{
}
//-----------------------------------------------------------
EventCorrelProc::EventCorrelProc(const char* name) :
   TGo4EventProcessor(name)
{
    cout << "**** EventCorrelProc: Create" << endl;
    
    fCorrel = new CorrelParameter("CorrelPar");
  AddParameter(fCorrel);
  if (fCorrel) fCorrel->PrintParameter(0,0);
  else cout << "**** ERRR - CorrelPar doesn't exist - program will crash.\n";
  
    get_used_systems();

 }
//-----------------------------------------------------------
EventCorrelProc::~EventCorrelProc()
{
  cout << "**** EventCorrelProc: Delete" << endl;

}
//-----------------------------------------------------------
Bool_t EventCorrelProc::BuildEvent(TGo4EventElement* dest)
{
  Bool_t isValid=kFALSE; // validity of output event
  
  EventAnlStore* cInput  = (EventAnlStore*) GetInputEvent();
  EventCorrelStore* cOutput = (EventCorrelStore*) dest;

  if((cInput==0) || !cInput->IsValid()){ /// input invalid
    cOutput->SetValid(isValid); /// invalid
    return isValid; /// must be same is for SetValid
  }
    isValid=kTRUE;
    static bool create =false;
    if(!create){

        Make_FRS_AIDA_Histos();
        Make_FRS_Prompt_AIDA_FATIMA_Ge_Histos();
        Make_FRS_AIDA_bPlast_Histos();
    //  if(cInput->pUsed_Systems[5])
        Make_FRS_Delayed_AIDA_Gamma_Histos();
       // Make_Timemachine_Histos();
    
    
    create=true;
    }
    //cout<<"cInput->pUsed_Systems[1] " <<cInput->pUsed_Systems[1] <<" cInput->pPrcID_Conv[1] " << cInput->pPrcID_Conv[1] <<endl; 
    
    ///Demand at least FRS to be activated to do correlations
      if(Used_Systems[0]==1) {
  Process_FRS_AIDA(cInput, cOutput); 
  Process_FRS_AIDA_bPlast(cInput, cOutput); 
  Process_FRS_Prompt_AIDA_FATIMA_Ge(cInput, cOutput); ///Prompt gammas
  Process_FRS_Delayed_AIDA_Gamma(cInput, cOutput); ///Beta delayed gammas
      }
  //Process_Timemachine_Histos(cInput, cOutput);
  event_number = cInput->pEvent_Number;
  cOutput->cEvent_number = event_number;
  //cout<<"event " << event_number<<endl;
  ///White Rabbit inputs
    AIDA_WR = cInput->pAIDA_WR;   
    FRS_WR = cInput->pFRS_WR;
    bPLAS_WR = cInput->pbPLAS_WR;
    FAT_WR = cInput->pFAT_WR;
    GAL_WR = cInput->pGAL_WR;
    
    cOutput->cAIDA_WR = cInput->pAIDA_WR;   
    cOutput->cFRS_WR = cInput->pFRS_WR;
    cOutput->cbPlas_WR = cInput->pbPLAS_WR;
    cOutput->cFAT_WR = cInput->pFAT_WR;
    cOutput->cGAL_WR = cInput->pGAL_WR;
    
      ///FRS Outputs    
    cOutput->cFRS_AoQ = cInput->pFRS_AoQ;   
    cOutput->cFRS_ID_x2 = cInput-> pFRS_ID_x2;  
    cOutput->cFRS_ID_x4 = cInput-> pFRS_ID_x4;  
    cOutput->cFRS_z = cInput-> pFRS_z;  
    cOutput->cFRS_z2 = cInput-> pFRS_z2;
//    
    
   //  cout<<"AIDA_WR " << AIDA_WR <<" FAT_WR " << FAT_WR <<" AIDA_WR - FAT_WR "<< AIDA_WR - FAT_WR<< endl;
  ///Gates input TESTING!!!
//     cInput->pFRS_Z_Z2_pass=true;
//     cInput->pFRS_x2AoQ_pass=true;
 

 cOutput->SetValid(isValid);
  return isValid;
}
    //void Make_Timemachine_Histos(){}
    //void Process_Timemachine_Histos(){}
 /**----------------------------------------------------------------------------------------------**/
 /**----------------------------------     FRS-AIDA (Implanted ion)   -------------------------**/
 /**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Make_FRS_AIDA_Histos(){
      TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
      hA_FRS_Z1Z2_implants_strip_xy.resize(conf->DSSDs());
      //hA_FRS_Z1Z2_implants_pos_xy.resize(conf->DSSDs());
      hA_FRS_Z1Z2_implants_e.resize(conf->DSSDs());
      hA_FRS_Z1Z2_implants_e_xy.resize(conf->DSSDs());
      
      hA_FRS_dT = MakeTH1('I',"Correlations/AIDA-FRS/WR_timediff","AIDA-FRS WR Time difference ",16000,-40000,40000);
     for (int i = 0; i < conf->DSSDs(); ++i)
  {    
        hA_FRS_Z1Z2_implants_strip_xy[i] = MakeTH2('I', Form("Correlations/AIDA-FRS/Implants/Z1Z2_Gate/DSSD%d_implants_strip_XY_Z1Z2g", i+1), Form("DSSD %d implant hit pattern, FRS Z1 Z2 gated", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
        //hA_FRS_Z1Z2_implants_pos_xy[i] = MakeTH2('D', Form("Correlations/AIDA-FRS/Implants/Z1Z2_Gate/DSSD%d_implants_pos_XY_Z1Z2g", i+1), Form("DSSD %d implant position FRS Z1 Z2 gated", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
        
        hA_FRS_Z1Z2_implants_e[i] = MakeTH1('F', Form("Correlations/AIDA-FRS/Implants/Z1Z2_Gate/DSSD%d_implants_energy_Z1Z2g", i+1), Form("DSSD %d implant energy FRS Z1 Z2 gated", i+1), 1000, 0, 10000, "Implant Energy/MeV");
        
      //  hA_FRS_Z1Z2_implants_e_xy[i] = MakeTH2('F', Form("Correlations/AIDA-FRS/Implants/Z1Z2_Gate/DSSD%d_implants_energy_XY_Z1Z2g", i+1), Form("DSSD %d implant front energy vs back energy FRS Z1 Z2 gated", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
    
//         hA_FRS_ZAoQ_implants_strip_xy[i] = MakeTH2('I', Form("Correlations/AIDA-FRS/Implants/ZvsAoQ_Gate/DSSD%d_implants_strip_XY_ZvsAoQ_Gate", i+1), Form("DSSD %d implant hit pattern FRS ZvsAoQ Gate", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
//            
//         hA_FRS_ZAoQ_implants_e[i] = MakeTH1('F', Form("Correlations/AIDA-FRS/Implants/ZvsAoQ_Gate/DSSD%d_implants_energy_ZvsAoQ_Gate", i+1), Form("DSSD %d implant energy FRS AoQ vs Z Gate", i+1), 1000, 0, 10000, "Implant Energy/MeV");
        
        
        ///Z vs AoQ
        for(int g=0; g<8; g++){
        hA_FRS_ZAoQ_implants_strip_xy[g][i] = MakeTH2('I', Form("Correlations/AIDA-FRS/Implants/ZvsAoQ_Gated/DSSD%d_implants_strip_XY_ZvsAoQ_Gate%d", i+1,g), Form("DSSD %d implant hit pattern FRS ZvsAoQ Gate%d", i+1,g), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
           
        hA_FRS_ZAoQ_implants_e[g][i] = MakeTH1('F', Form("Correlations/AIDA-FRS/Implants/ZvsAoQ_Gated/DSSD%d_implants_energy_ZvsAoQ_Gate%d", i+1,g), Form("DSSD %d implant energy FRS AoQ vs Z Gate%d", i+1,g), 1000, 0, 10000, "Implant Energy/MeV");
        
        
         hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy[g][i] = MakeTH2('I', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_strip_Z1Z2_x2x4AoQ_Gate%d", i+1,g), Form("DSSD %d implant hit pattern FRS Z1Z2_x2x4AoQ_Gate Gate%d", i+1,g), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
           
        hA_FRS_Z1Z2_x2x4AoQ_implants_e[g][i] = MakeTH1('F', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_energy_Z1Z2_x2x4AoQ_Gate%d", i+1,g), Form("DSSD %d implant energy FRS Z1Z2_x2x4AoQ_Gate%d", i+1,g), 1000, 0, 10000, "Implant Energy/MeV");
        
        }
//         
       
       
       
      /* // hA_FRS_ZAoQ_implants_e_xy[i] = MakeTH2('F', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_energy_XY_Z1Z2IDx2x4AoQg", i+1), Form("DSSD %d implant front energy vs back energy FRS Z1 Z2, and X2AoQ or X4AoQ And ZAoQ gated", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
        hA_FRS_ZAoQ_implants_time_delta[i] = MakeTH1('F', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_time_delta_Z1Z2IDx2x4AoQg", i+1), Form("DSSD %d implant front vs back time FRS Z1 Z2, and X2AoQ And X4AoQ gated And ZAoQ", i+1), 1000, -10000, 10000, "Time Difference/ns");
        
        hA_FRS_ZAoQ_implants_strip_1d[i] = MakeTH1('I', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_strip_1d_Z1Z2IDx2x4AoQg", i+1), Form("DSSD %d implant 1D hit pattern FRS Z1 Z2, and X2AoQ or X4AoQ gated And ZAoQ", i+1), 256, 0, 256, "Strip number");
        
        hA_FRS_ZAoQ_implants_per_event[i] = MakeTH1('I', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_per_event_Z1Z2IDx2x4AoQg", i+1), Form("DSSD %d implants per event FRS Z1 Z2, and X2AoQ or X4AoQ And ZAoQ gated", i+1), 100, 0, 100, "Number of implants");  
        *///hA_FRS_ZAoQ_implants_strip_xy_dssdg[i] = MakeTH2('I', Form("Correlations/AIDA-FRS/Implants/Z1Z2_x2x4AoQ_Gate/DSSD%d_implants_strip_XY_Z1Z2IDx2x4AoQg_DSSDGate", i+1), Form("DSSD %d implant hit pattern FRS Z1 Z2, and X2AoQ or X4AoQ gated, DSSD Ion gate", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip"); 
    }
     ///2D AIDA Ion position Gates DSSD 1-3
  
    Float_t init_ID_AIDA_ION_DSSD1[7][2] =
    
    {{34, 55.45},{37,50.5},{40, 50.5}, {41, 59}, {37,60}, {33,60}, {34, 55}}; 
   
    Float_t init_ID_AIDA_ION_DSSD2[7][2] =
     {{0.0, 0.0},{40.0, 0.0},{80.0, 0.0}, {128.0, 0.0}, {128.0, 126.0}, {0.0, 128.0}, {0.0, 0.0}}; 
    
    Float_t init_ID_AIDA_ION_DSSD3[7][2] =
      {{0.0, 0.0},{40.0, 0.0},{80.0, 0.0}, {128.0, 0.0}, {128.0, 126.0}, {0.0, 128.0}, {0.0, 0.0}};  
    
      cAIDA_IMPgate_DSSD1 = MakePolyCond("cID_AIDA_IMP_DSSD1","FRS Gated AIDA ion pos DSSD1",7,init_ID_AIDA_ION_DSSD1, "AIDA Implantation DSSD1");
      
      cAIDA_IMPgate_DSSD2 = MakePolyCond("cID_AIDA_IMP_DSSD2","FRS Gated AIDA ion pos DSSD1",7,init_ID_AIDA_ION_DSSD2, "AIDA Implantation DSSD2");
      
      cAIDA_IMPgate_DSSD3 = MakePolyCond("cID_AIDA_IMP_DSSD3","FRS Gated AIDA ion pos DSSD1",7,init_ID_AIDA_ION_DSSD3, "AIDA Implantation DSSD3");
     
}
/**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Process_FRS_AIDA(EventAnlStore* cInputMain, EventCorrelStore* cOutput){
    
     TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
     //Branches from AnlProc     
      for (auto& cInputRef : cInputMain->pAida)
     {
       auto* cInput = &cInputRef;
       ///AIDA Implants
         std::vector<AidaHit> hits = cInput->Implants;     
        if(AIDA_WR>0 && FRS_WR>0){
         dT_AIDA_FRS = AIDA_WR - FRS_WR;
	 
	 //cout<<"AIDA_WR " << AIDA_WR << " FRS_WR  " << FRS_WR << " dT_AIDA_FRS " << dT_AIDA_FRS <<endl;
         cOutput->cdT_AIDA_FRS = dT_AIDA_FRS;
       
        hA_FRS_dT -> Fill(dT_AIDA_FRS); 
        }
        for (auto& i : hits)
      {
         AidaHit hit = i;
        //cout<<" hit.Energy " << hit.Energy << endl;
         if(dT_AIDA_FRS > fCorrel->GFRS_AIDA_TLow && dT_AIDA_FRS < fCorrel->GFRS_AIDA_THigh){
         
          ///Gate on FRS Z1_Z2 AND (X2vsAoQ OR X4vsAoQ)-> AIDA Implantation
          if(cInputMain->pFRS_Z_Z2_pass[0]==true ){ 
          
            hA_FRS_Z1Z2_implants_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
            //hA_FRS_Z1Z2_implants_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
            hA_FRS_Z1Z2_implants_e[hit.DSSD - 1]->Fill(hit.Energy);           
          //  hA_FRS_Z1Z2_implants_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
    
          }
         }
          
         // cout<<"cInputMain->pFRS_Z_Z2_pass  " << cInputMain->pFRS_Z_Z2_pass << " cInputMain->pFRS_x2AoQ_pass  " << cInputMain->pFRS_x2AoQ_pass << " cInputMain->pFRS_x4AoQ_pass " << cInputMain->pFRS_x4AoQ_pass << endl;
         
    /*if(cInputMain->pFRS_Z_Z2_pass==true && cInputMain->pFRS_ZAoQ_pass==true && (cInputMain->pFRS_x2AoQ_pass==true||cInputMain->pFRS_x4AoQ_pass==true)){  */ 
        ///Gate on FRS Z vs AoQ -> AIDA Implantation
        ///Gate1
        for(int g=0; g<8;g++){
        if( cInputMain->pFRS_ZAoQ_pass[g]==true ){     
          
         hA_FRS_ZAoQ_implants_strip_xy[g][hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
         hA_FRS_ZAoQ_implants_e[g][hit.DSSD - 1]->Fill(hit.Energy);
               }
               
         if( cInputMain->pFRS_Z_Z2_pass[g]==true && (cInputMain->pFRS_x2AoQ_pass[g]==true||cInputMain->pFRS_x4AoQ_pass[g]==true)){     
          
         hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy[g][hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
         hA_FRS_Z1Z2_x2x4AoQ_implants_e[g][hit.DSSD - 1]->Fill(hit.Energy);
                    }      
                }
             }  
          }
       }
    
 
  /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------  FRS-AIDA-FATIMA/GALILEO  (Isomers)  -------------**/
 /**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Make_FRS_Prompt_AIDA_FATIMA_Ge_Histos(){
     //hA_implant_FatdT =  MakeTH1('I',"Correlations/AIDA-Prompt-Fatima/Fat-Implant_dT","T Diff Fatima - Implant",90000,-90000,90000,"Time[ns]", "Counts");
     
      hA_FRSWR_FatWR =  MakeTH1('I',"Correlations/FRS-Prompt_Fatima/FRS-FATWR_dT","T Diff FRS WR -Fatima QR ",10000,-10000,10000,"Time[ns]", "Counts");
     
       
      hA_FRS_FatE = MakeTH1('D', "Correlations/FRS-Prompt_Fatima/Fat_EnergySum_allFRS", "Galileo Energy FRS (all) gated",500,0,4000);
      
      
       
       for(int i=0; i<8; i++){
        hA_FRS_PID_FatE[i]  = MakeTH1('F', Form("Correlations/FRS-Prompt_Fatima/AoQZ_Gated_Energy/Fat_EnergySum_PIDGate%d", i), Form("Fatima Energy FRS AoQvsZ gated %d", i), 500, 0, 4000, "Energy/keV");
      
        hA_FRS_Z1Z2_X2AoQX4AoQ_FatE[i]  = MakeTH1('F', Form("Correlations/FRS-Prompt_Fatima/Fat_EnergySum_Z1Z2_X2AoQX4AoQGate%d", i), Form("Fatima Energy FRS Z1Z2_X2AoQX4AoQ gated %d", i), 500, 0, 4000, "Energy/keV");
       
        hA_FRS_FatEvsT[i] = MakeTH2('D',Form("Correlations/FRS-Prompt_Fatima/AoQZ_Gated_FatE_vs_T/FRS_WR-FAT_WR_vs_FatE_gate%2d",i), Form("T Diff FRS SC41 - Fatima T vs Fatima Energy AoQZ Gate%2d",i),500, 0, 4000, 100,-10000,10000);
     
       }
 //     hA_FRS_FatEvsT[i]  = MakeTH2('D',"Correlations/FRS-Prompt_Fatima/FRS_WR-FAT_WR_vs_FatE","T Diff FRS WR - Fatima WR vs Fatima Energy", 1250, 0, 5000, 100,-10000,10000,"Fat Energy (keV)", "FRS WR - Fat WR time (ns)");
      
      
     
      hA_FRSWR_GalWR =  MakeTH1('I',"Correlations/FRS-Prompt_Ge/FRS-GALWR_dT","T Diff FRS WR -Galileo WR ",2000,-10000,10000,"Time[ns]", "Counts");
      
      hA_FRS_GalE = MakeTH1('D', "Correlations/FRS-Prompt_Ge/Gal_Energy_allFRS", "Galileo Energy FRS (all) gated",2000, 0, 4000);
      
      for(int i=0; i<8; i++){
      hA_FRS_PID_GalE[i]  = MakeTH1('F', Form("Correlations/FRS-Prompt_Ge/AoQZ_Gated_Energy/Gal_EnergySum_PIDGated%d", i), Form("Galileo Energy FRS PID gated %d", i), 4000, 0, 4000, "Energy/keV");
      
      hA_FRS_GalEvsT[i]  = MakeTH2('D',Form("Correlations/FRS-Prompt_Ge/AoQZ_Gated_E_vs_T/FRS_WR-Gal_WR_vs_GalE_gate%2d",i), Form("T Diff FRS - Ge WR vs Ge Energy AoQZ Gate%2d",i),4000, 0, 4000, 1500,-10000,5000);
      
      hA_SC41_GalEvsT[i]  = MakeTH2('D',Form("Correlations/FRS-Prompt_Ge/AoQZ_Gated_E_vs_SC41dT/Sc41-GaldT_vs_GalE_gate%2d",i), Form("T Diff Sc41 - Ge T vs Ge Energy AoQZ Gate%2d",i),4000, 0, 4000, 1500,-10000,5000);
      
      hA_FRS_Z1Z2_X2AoQX4AoQ_GalE[i]  = MakeTH1('F', Form("Correlations/FRS-Prompt_Ge/Z1Z2_X2AoQX4AoQGated/Gal_EnergySum_Z1Z2_X2AoQX4AoQGated%d", i), Form("Galileo Energy FRS Z1Z2_X2AoQX4AoQ gated %d", i), 2000, 0, 4000, "Energy/keV");
     
      }
//       hA_FRS_GalEvsT  = MakeTH2('D',"Correlations/FRS-Prompt_Ge/FRS_WR-GAL_WR_vs_GeE","T Diff FRS WR - Galileo WR vs Galileo Energy", 1250, 0, 5000, 100,-10000,10000,"Ge Energy (keV)", "FRS WR - Ge WR time (ns)");
      
     hA_bPlasWR_FatWR =  MakeTH1('I',"Correlations/bPlast-Prompt_Fatima/bPlasWR-FATWR_dT","T Diff bPlast WR -Fatima WR ",10000,-10000,10000,"Time[ns]", "Counts");
      
     hA_BplasDet1_FatE = MakeTH1('D', "Correlations/bPlast-Prompt_Fatima/Fat_EnergySum_bPlasGated_Det1", "Fatima Energy bPlas detector 1 gated",2000,0,2000);
       
     hA_BplasDet2_FatE = MakeTH1('D', "Correlations/bPlast-Prompt_Fatima/Fat_EnergySum_bPlasGated_Det1", "Fatima Energy bPlas detector 1 gated",2000,0,2000);
      
      hA_bPlasWR_GalWR =  MakeTH1('I',"Correlations/bPlast-Prompt_Ge/bPlast-GALWR_dT","T Diff bPlast WR -Galileo WR ",10000,-10000,10000,"Time[ns]", "Counts");
      
      hA_BplasDet1_GalE = MakeTH1('D', "Correlations/bPlast-Prompt_Ge/Gal_EnergySum_bPlasGated_Det1", "Galileo Energy bPlas detector 1 gated",2000,0,2000);
      
      hA_BplasDet2_GalE = MakeTH1('D', "Correlations/bPlast-Prompt_Ge/Gal_EnergySum_bPlasGated_Det2", "Galileo Energy bPlas detector 2 gated",2000,0,2000);
       
      hA_bPlasWR_FatWR =  MakeTH1('I',"Correlations/bPlast-Prompt_Fatima/bPlast-FATWR_dT","T Diff bPlast WR -Fatima WR ",10000,-10000,10000,"Time[ns]", "Counts");
    
     
}
                    
 void EventCorrelProc::Process_FRS_Prompt_AIDA_FATIMA_Ge(EventAnlStore* cInputMain, EventCorrelStore* cOutput){
   
     Long64_t aida_imptime_fatgal,aida_imptime_fatgal_noFRS;
     Long64_t fat_aida_imptime_dT, gal_aida_imptime_dT, bplas_aida_imptimedT ;
     Long64_t dT_FRS_Fatima = 0, dT_FRS_Galileo = 0;
     Long64_t dT_bPlas_Fatima = 0, dT_bPlas_Gal = 0;
//      Long64_t dT_AIDA_Fatima;
//      dT_AIDA_Fatima = AIDA_WR - FAT_WR;
  
 //  hA_Fatima_dT->Fill(dT_AIDA_Fatima);
    if(FRS_WR>0 && FAT_WR>0) dT_FRS_Fatima = FRS_WR - FAT_WR;
    if(dT_FRS_Fatima>0) hA_FRSWR_FatWR->Fill(dT_FRS_Fatima);

    if(FRS_WR>0 && GAL_WR>0)  dT_FRS_Galileo = FRS_WR - GAL_WR;
    if(dT_FRS_Galileo>0) hA_FRSWR_GalWR->Fill(dT_FRS_Galileo);
    
    if(bPLAS_WR>0 && FAT_WR>0) dT_bPlas_Fatima = bPLAS_WR - FAT_WR;
    if(dT_bPlas_Fatima>0) hA_bPlasWR_FatWR->Fill(dT_bPlas_Fatima);
    
    if(bPLAS_WR>0 && GAL_WR>0) dT_bPlas_Gal = bPLAS_WR - GAL_WR;
    if(dT_bPlas_Gal>0) hA_bPlasWR_GalWR->Fill(dT_bPlas_Gal);
     
     
 ///----------------------------- GALILEO FRS Prompt correlations----------------------------------------------- ///   
//     for(int gate=0;gate<8;gate++) if(cInputMain->pFRS_ZAoQ_pass[gate]==1)cout<<"cInputMain->pFRS_ZAoQ_pass[gate] "<< cInputMain->pFRS_ZAoQ_pass[gate] << " i " << gate << endl;
     for(int g=0; g<GALILEO_MAX_DETS; g++){
            for (int h=0; h<GALILEO_CRYSTALS; h++){
          
                 ///WR Time gate FRS-Ge
                if(dT_FRS_Galileo>fCorrel->GFRS_GAL_TLow && dT_FRS_Galileo < fCorrel->GFRS_GAL_THigh){
		
                    if(cInputMain->pGal_EAddback[g][h]>0){
                        hA_FRS_GalE->Fill(cInputMain->pGal_EAddback[g][h]);           
                        
                        for(int gate=0;gate<8;gate++){
			  //
                        if(cInputMain->pFRS_ZAoQ_pass[gate]==true){
//  cout<<"cInputMain->pFRS_ZAoQ_pass[gate] "<<cInputMain->pFRS_ZAoQ_pass[gate] << " gate " << gate <<" cInputMain->pGal_EAddback[g][h] " << cInputMain->pGal_EAddback[g][h]<< endl;
			 if(cInputMain->pGal_EAddback[g][h]<200&&(cInputMain->pGal_T[8][0]-cInputMain->pGal_T[g][h])>500){
                         hA_FRS_PID_GalE[gate]->Fill(cInputMain->pGal_EAddback[g][h]);
			 }
			  if(cInputMain->pGal_EAddback[g][h]>199&&cInputMain->pGal_T[8][0]-cInputMain->pGal_T[g][h]>50){
                         hA_FRS_PID_GalE[gate]->Fill(cInputMain->pGal_EAddback[g][h]);
			 }
			 
                         hA_FRS_GalEvsT[gate]->Fill(cInputMain->pGal_EAddback[g][h],dT_FRS_Galileo);
                         
                         hA_SC41_GalEvsT[gate]->Fill(cInputMain->pGal_EAddback[g][h],cInputMain->pGal_T[8][0]-cInputMain->pGal_T[g][h]);
                        }
                       if( cInputMain->pFRS_Z_Z2_pass[gate]==true && (cInputMain->pFRS_x2AoQ_pass[gate]==true||cInputMain->pFRS_x4AoQ_pass[gate]==true)){   
                        hA_FRS_Z1Z2_X2AoQX4AoQ_GalE[gate]->Fill(cInputMain->pGal_EAddback[g][h]);
                         }
                     }
                   }
                }
            }
     }
     
///----------------------------- FATIMA FRS Prompt correlations ---------------------------------------------------///
     for(int k=0; k<cInputMain->pFatmult; k++){
         
          ///Energy WR dT
                 
                 ///WR Time gate FRS-LaBr3
                if(dT_FRS_Fatima>fCorrel->GFRS_FAT_TLow && dT_FRS_Fatima < fCorrel->GFRS_FAT_THigh){
                    if(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]>0){
                        hA_FRS_FatE->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);           
                        
                         for(int gate=0;gate<8;gate++){
                        if(cInputMain->pFRS_ZAoQ_pass[gate]==true){
                         hA_FRS_PID_FatE[gate]->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);
                         hA_FRS_FatEvsT[gate]->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]],dT_FRS_Fatima);
                         
                   if( cInputMain->pFRS_Z_Z2_pass[gate]==true && (cInputMain->pFRS_x2AoQ_pass[gate]==true||cInputMain->pFRS_x4AoQ_pass[gate]==true)){  
                            hA_FRS_Z1Z2_X2AoQX4AoQ_FatE[gate]->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);
       
                         }
                        }
                      }
                   }
                }
            }
     ///-----------------------------  Beta E prompt gated Fatima-----------------------------------///

   
            // for(int j=1; j<3;j++){
              for(int k=0; k<16; k++) { ///loop on bPlast
                   for(int l=0; l<10; l++) {
                       ///WR Time gate bPlas-LaBr3
        if(dT_bPlas_Fatima>fCorrel->GbPlas_FAT_TLow && dT_bPlas_Fatima<fCorrel->GbPlas_FAT_THigh){
                     
            for(int k=0; k<cInputMain->pFatmult; k++){
                if(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]>0){
            if(cInputMain-> pbPlas_ToTCalib[1][k][l]>fCorrel->GbPlas_Egate_low && cInputMain-> pbPlas_ToTCalib[1][k][l]<fCorrel->GbPlas_Egate_high){               
                hA_BplasDet1_FatE->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);   
                                }
                                
               if(cInputMain-> pbPlas_ToTCalib[2][k][l]>fCorrel->GbPlas_Egate_low && cInputMain-> pbPlas_ToTCalib[2][k][l]<fCorrel->GbPlas_Egate_high){           
                hA_BplasDet2_FatE->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);   
                            }                 
                        }
                   }
                }
           
        
                       
       ///-----------------------------  Beta E prompt gated Galileo-----------------------------------///                
          if(dT_bPlas_Gal>fCorrel->GbPlas_GAL_TLow && dT_bPlas_Gal<fCorrel->GbPlas_GAL_THigh){
               for(int g=0; g<GALILEO_MAX_DETS; g++){
                    for (int h=0; h<GALILEO_CRYSTALS; h++){
             if(cInputMain-> pbPlas_ToTCalib[1][k][l]>fCorrel->GbPlas_Egate_low && cInputMain-> pbPlas_ToTCalib[1][k][l]<fCorrel->GbPlas_Egate_high && cInputMain->pGal_EAddback[g][h]>0){
                 
                hA_BplasDet1_GalE->Fill(cInputMain->pGal_EAddback[g][h]);   
                                }
               if(cInputMain-> pbPlas_ToTCalib[2][k][l]>fCorrel->GbPlas_Egate_low && cInputMain-> pbPlas_ToTCalib[2][k][l]<fCorrel->GbPlas_Egate_high){           
                hA_BplasDet2_GalE->Fill(cInputMain->pGal_EAddback[g][h]);   
                                    }
                                }
                            }
                        }
                    }
              }
  
//    cOutput->cFatE[cInputMain->pFat_QDC_ID[k]] = cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]];
//       if(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]>0 && FAT_WR>0){
//          //for(int j =0; j<50; j++){
//           
//           fat_aida_imptime_dT =  (FAT_WR - aida_imptime_fatgal)*0.001; /// in mus
//           hA_implant_FatdT->Fill(fat_aida_imptime_dT);
//           if(fat_aida_imptime_dT>10 &&  cOutput->cImplantIterator>-1){
//     
//                 hA_implant_FatE->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);
// //                hA_implant_FatdT_FatE->Fill( cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]],fat_aida_imptime_dT);
//                //     }              
//                 }
//             }
        
           
     
     
     
//     if(event_number >0){
//         // cout<<"event " << event_number <<" FAT_WR " << FAT_WR <<endl;
//     AidaAnlData* cInput = &cInputMain->pAida;
//     
//   
//     std::vector<AidaHit> imphits = cInput->Implants;
//      for(auto& i : imphits){
//       AidaHit imphit = i;
//       ///FRS -AIDA gate (inc. WR)
//       aida_imptime_fatgal_noFRS=  imphit.Time; 
// //       if(cInputMain->pFRS_Z_Z2_pass==true && (cInputMain->pFRS_x2AoQ_pass==true||cInputMain->pFRS_x4AoQ_pass==true) && cInputMain->pFRS_ZAoQ_pass==true  ){  
//          if(cInputMain->pFRS_ZAoQ_pass==true  ){ 
//       aida_imptime_fatgal=  imphit.Time; 
//     
//       }
//       
//       // for(int j =0; j<10000; j++){
//        ///FATIMA AIDA
//        ///For Combined system below
// //       for(int k=0; k<cInputMain->pFat_firedQDC_Comb; k++){
// //       if(cInputMain->pFat_QDC_E_Comb[cInputMain->pFat_QDC_ID_Comb[k]]>0 && FAT_WR>0){
// //          //for(int j =0; j<50; j++){
// //           
//            fat_aida_imptime_dT =  (FAT_WR - aida_imptime_fatgal); /// in mus
//       
//                  hA_implant_FatdT->Fill(fat_aida_imptime_dT);
//         
//     }
// ///End of combined system
// 
// ///Just VME Fatima
//         for(int k=0; k<cInputMain->pFatmult; k++){
// 
// //    cOutput->cFatE[cInputMain->pFat_QDC_ID[k]] = cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]];
//       if(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]>0 && FAT_WR>0){
//          //for(int j =0; j<50; j++){
//           
//           fat_aida_imptime_dT =  (FAT_WR - aida_imptime_fatgal)*0.001; /// in mus
//           hA_implant_FatdT->Fill(fat_aida_imptime_dT);
//           if(fat_aida_imptime_dT>10 &&  cOutput->cImplantIterator>-1){
//     
//                 hA_implant_FatE->Fill(cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]]);
// //                hA_implant_FatdT_FatE->Fill( cInputMain->pFat_QDC_E[cInputMain->pFat_QDC_ID[k]],fat_aida_imptime_dT);
//                //     }              
//                 }
//             }
//         }
//         
// 
//         ///GALILEO AIDA
//         //for(int j =0; j<cInputMain->pGalFired; j++){
// //        cOutput-> cGalE = cInputMain->pGal_EAddback;
//         
//         if(bPLAS_WR>0 && aida_imptime_fatgal_noFRS>0){
//         bplas_aida_imptimedT = (bPLAS_WR-aida_imptime_fatgal);
//         }
//         
//         hA_implant_bPlasdT->Fill(bplas_aida_imptimedT);
//         ///Loop over galilieo
//         for(int g=0; g<GALILEO_MAX_DETS; g++){
//             for (int h=0; h<GALILEO_CRYSTALS; h++){
//                 
//                 if(cInputMain->pGal_EAddback[g][h]>0 && GAL_WR>0){
//                     gal_aida_imptime_dT =  (GAL_WR - aida_imptime_fatgal); 
//              
//        //  cout<< " GAL_WR " << GAL_WR <<" aida_imptime_fatgal  " << aida_imptime_fatgal <<  " gal_aida_imptime_dT " << gal_aida_imptime_dT << endl;
//          
//              hA_implant_GaldT->Fill(gal_aida_imptime_dT);
//              if(gal_aida_imptime_dT>10 &&  cOutput->cImplantIterator>-1){
//              hA_implant_GalE->Fill(cInputMain->pGal_EAddback[g][h]);
//              }
//              
//              ///AIDA Implant - bPlas && Beta E gated Galileo
//              for(int j=1; j<3;j++){
//                    for(int k=0; k<16; k++) {
//                         for(int l=0; l<10; l++) {
//              if(bplas_aida_imptimedT*1e-6>0.5 && bplas_aida_imptimedT*1E-6<25 &&  cInputMain-> pbPlas_ToTCalib[j][k][l]>17000&& cInputMain-> pbPlas_ToTCalib[j][k][l]<50000){
//                  
//               hA_implantBplas_GalE->Fill(cInputMain->pGal_EAddback[g][h]);   
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//            //  cout<<"event " <<event_number << " cInputMain->pGal_EAddback " <<cInputMain->pGal_EAddback<< " gal_aida_imptime_dT " << gal_aida_imptime_dT<<endl;
//             // hA_implant_GaldT_GalE->Fill( cInputMain->pGal_EAddback,gal_aida_imptime_dT);
//             
//             
//           }
//    
//        aida_imptime_fatgal = 0;   
//        aida_imptime_fatgal_noFRS = 0;
//        // }
//     }        
}

 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------- (FRS)-AIDA-bPlastic  (Beta Decay)  -------------------**/
 /**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Make_FRS_AIDA_bPlast_Histos(){
     TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
     hA_bPlas_decays_strip_xy.resize(conf->DSSDs());
     hA_bPlas_decays_e.resize(conf->DSSDs());
     
     hA_bPlas_dT = MakeTH1('I',"Correlations/AIDA-bPlast/Aida_bPlas_WR_timediff","AIDA-bPlas WR Time difference ",20000,-200000,200000);
     
     hFRS_bPlas_dT = MakeTH1('I',"Correlations/bPlast-FRS/bPlas_FRS_WR_timediff","bPlas FRS WR Time difference ",20000,-200000,200000);
     
     hADec_bPlas_dT = MakeTH1('I',"Correlations/AIDA-bPlast/AidaDecay_bPlas_WR_timediff","AIDA Decay-bPlas WR Time difference ",20000,-200000,200000);
      
     hA_impdec_dT = MakeTH1('I',"Correlations/AIDA-bPlast/Implant_Decay_dT","T Diff Decay - Implant ",100000,-2000000,2000000,"Time[ms]", "Counts");
      
     hA_impdec_dT_FRS_gated = MakeTH1('I',"Correlations/AIDA-bPlast/Implant_Decay_dT_FRS_Gated","T Diff Decay - Implant FRS Gated",90000,-90000,90000,"Time[ms]", "Counts");
      
     hADecay_bPlas_E = MakeTH1('I',"Correlations/AIDA-bPlast/AidaDecay_bPlas_CoincEnergy","AIDA-bPlas Decay Energy correlated bPlast Energy ",25000, 0., 150000.);
     
    // hA_FRSgated_bPlas_E = MakeTH1('I',"Correlations/AIDA-bPlast/Aida_FRSGated_bPlas_CoincEnergy","AIDA-bPlas Implant-Decay Time Coincident FRS Gated bPlast Energy ",25000, 0., 150000.);
     
     hA_gatedE_bPlas_E = MakeTH1('I',"Correlations/AIDA-bPlast/AidaImplant_bPlas_CoincEnergy","AIDA-bPlas Implant gated coincident bPlast Energy ",25000, 0., 1500000.);
     //AIDA
     for (int i = 0; i < conf->DSSDs(); ++i)
  {      
     
     hA_bPlas_decays_e[i] = MakeTH1('F', Form("Correlations/AIDA-bPlast/AIDA Decay/DSSD%d_decay_E_Implant+bPlastFires", i+1), Form("DSSD %d decay energy: Implant and bPlast fires", i+1), 1000, 0, 20000, "Decay Energy/keV");
     
      hA_bPlas_decays_strip_xy[i] = MakeTH2('D', Form("Correlations/AIDA-bPlast/AIDA Decay/DSSD%d_decays_strip_XY_Implant+bPlastFires", i+1), Form("DSSD %d decay hit pattern:Implant and bPlast fires", i+1),128, 0, 128, 128, 0, 128, "X strip", "Y strip");
  }
  
        for (int i =1; i<4; i++)
        {
            for(int j =0; j<16; j++){
//        sprintf(chis,"Correlation/AIDA-bPlast/Aida_bPlas_CoincEnergy_ch/ToT Ch.%2d", i);
//        sprintf(chead,"AIDA-bPlas Coincident bPlast Energy Ch. %2d", i);
//        hA_bPlas_E_Ch[i] = MakeTH1('I', chis,chead, 25000, 0., 150000.);
       
	hA_bPlas_E_Ch[i][j] = MakeTH1('D', Form("Correlations/AIDA-bPlast/Aida_bPlas_CoincEnergy_ch/ToT Det.%2d Ch.%2d", i,j), Form("AIDA-bPlas Coincident bPlast Energy Det.%2d Ch.%2d", i,j),25000, 0., 1500000.);
         
	  
       hbPlas_FRS[i][j] = MakeTH1('D', Form("Correlations/bPlast-FRS/bPlas_FRS_CoincEnergy_ch/ToT_Plas Det.%2d Ch.%2d", i,j), Form("bPlas FRS AoQ Coincident bPlast Energy Det.%2d Ch.%2d", i,j),25000, 0., 150000.);
 
       hA_gatedE_bPlas_E_Ch[i][j] = MakeTH1('D', Form("Correlations/AIDA-bPlast/Aida_ImpEGate_bPlas_CoincEnergy_ch/ToT_AIDAEGate Det.%2d Ch.%2d", i,j), Form("AIDA-bPlas Implant gated Coincident bPlast Energy Det.%2d Ch.%2d", i,j),25000, 0., 1500000.);
        }
 }
  for(int k =0; k<16; k++){
    hFibre_FRS_AoQ[k] = MakeTH1('D', Form("Correlations/Fibre-FRS/Fibre_FRS_AoQGate_Energy_ch/ToT_Fibre_AoQGated_Ch.%2d", k), Form("Fibre FRS AoQ Coincident Energy Ch.%2d", k),25000, 0., 150000.);
    
    hFibre_FRS_all[k] = MakeTH1('D', Form("Correlations/Fibre-FRS/Fibre_FRS_CoincEnergy_ch/ToT_Fibre_FRSALL_Ch.%2d", k), Form("Fibre FRS ALL Coincident Energy Ch.%2d", k),25000, 0., 150000.);
  }
     //bPlastic
    /* for(int i =1; i<3; i++){    
        for (int j =0; j<16; j++){
      hA_bPlas_Energy_AidaImpGated[i][j] = MakeTH1('D', Form("Correlations/AIDA-bPlast/bPlastEnergy_AIDAImp_gated/bPlastEnergy_AIDA_Implant_Gated_Det. %2d Ch. %2d",  i,j), Form("Energy Calib Det. %2d Ch. %2d", i,j),6500, 0., 65000.);
      
      hA_bPlas_Energy_AidaDecayGated[i][j] = MakeTH1('D', Form("Correlations/AIDA-bPlast/bPlastEnergy_AIDADecay_gated/bPlastEnergy_AIDA_Implant_Gated_Det. %2d Ch. %2d",  i,j), Form("Energy Calib Det. %2d Ch. %2d", i,j),6500, 0., 65000.);

        }
     }*/     
}
 /**------------------------------------(FRS)-AIDA-bPlastic  (Beta Decay)----------------------------**/
 void EventCorrelProc::Process_FRS_AIDA_bPlast(EventAnlStore* cInputMain, EventCorrelStore* cOutput){    
    
     TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
      
        Long64_t AIDA_dT_imp_decay;
        Long64_t AIDA_dT_imp_decay_FRS_gated = 0;
        Long64_t AIDA_dT_decay_bPlas = 0;
        bool bPlasfired =false;
        bool SC41_fired = false;
        int ImpIterator =0, DecayIterator=0;
         //for(int i=0; i<32; i++)  cInputMain-> pbPlas_ToTCalib[i]=100;
         ///Get Plastic fired
        for(int i=0; i<48; i++){  

       if ( cInputMain-> pbPlas_ToTCalib[i]>0 && i>15) bPlasfired=true; ///VME-TAMEX 
             }
        
       
    
       
       dT_FRS_bPlas = FRS_WR-bPLAS_WR;
       hFRS_bPlas_dT->Fill(dT_FRS_bPlas);
             
             
             ///FRS bPlas 
             for(int j=1; j<4; j++){
                 for(int k=0; k<16; k++) {
                        for(int l=0; l<10; l++) {
	       if( cInputMain-> pbPlas_ToTCalib[j][k][l]>0 &&  cInputMain->pFRS_ZAoQ_pass[0]==true ) {
		if(j<3) hbPlas_FRS[j][k]->Fill(cInputMain-> pbPlas_ToTCalib[j][k][l]);
         
                        }
         
                    }
                 }
             }
             ///Fibre detector
              for(int k=0; k<16; k++) {
                        for(int l=0; l<10; l++) {
        
	       if(cInputMain-> pbPlas_ToTCalib[FIBRE_BOARD][k][l]>0){ 
               
       if(cInputMain->pFRS_ZAoQ_pass[0]==true)   hFibre_FRS_AoQ[k]->Fill(cInputMain-> pbPlas_ToTCalib[FIBRE_BOARD][k][l]);
    
            
        if(FRS_WR>0) hFibre_FRS_all[k]->Fill(cInputMain-> pbPlas_ToTCalib[FIBRE_BOARD][k][l]);
                        }
                    }
              }
              
              
             for (auto& cInputRef : cInputMain->pAida)
        {
          auto* cInput = &cInputRef;  
               
          if(AIDA_WR>0 && bPLAS_WR>0){
       dT_AIDA_bPlas = AIDA_WR - bPLAS_WR;
       hA_bPlas_dT->Fill(dT_AIDA_bPlas);
                 }
     std::vector<AidaHit> imphits = cInput->Implants;
     std::vector<AidaHit> dechits = cInput->Decays; 
    // if(cInputMain-> pbPlas_ToTCalib[9]>4000) SC41_fired=true;
   // if(event_number >0){
     ///IMPLANTS
     for(auto& i : imphits){
      AidaHit imphit = i;
    
      int ImpstripX = imphit.StripX;
      int ImpstripY = imphit.StripY;  
      int ImpDSSD = imphit.DSSD-1;
      //cout<<"CORR EVENT " << event_number << " ImpstripX " << ImpstripX <<" imphit.Energy " <<imphit.Energy <<endl;
                ///Define the implantation pixel time for all AIDA (non FRS gated)
                    aida_imptime[ImpDSSD][ImpstripX][ImpstripY] =  imphit.Time;  
                     
               ///Define the implantation pixel time for FRS Gated AIDA
  /*if(cInputMain->pFRS_Z_Z2_pass==true && cInputMain->pFRS_ZAoQ_pass==true && (cInputMain->pFRS_x2AoQ_pass==true||cInputMain->pFRS_x4AoQ_pass==true)&& cOutput->cdT_AIDA_FRS> fCorrel->GFRS_AIDA_TLow && dT_AIDA_FRS < fCorrel->GFRS_AIDA_THigh ){       */  
  ///CHECK LATER NOTE
  if(cInputMain->pFRS_ZAoQ_pass[0]==true && cOutput->cdT_AIDA_FRS> fCorrel->GFRS_AIDA_TLow && dT_AIDA_FRS < fCorrel->GFRS_AIDA_THigh ){ 

            aida_imptime_FRS_gated[ImpDSSD][ImpstripX][ImpstripY] =  imphit.Time; 
            }
            
          //  cout<<"ImpIterator " << ImpIterator << " hit.Energy "<< imphit.Energy<<endl;
           /// AIDA implant energy gated bPlastic
       if(dT_AIDA_bPlas>12000 && dT_AIDA_bPlas < 15000){
           
           for(int j=1; j<4; j++){
               for(int k=0; k<16; k++) {
                        for(int l=0; l<10; l++) {
              cOutput-> cbPlasE[j] = cInputMain-> pbPlas_ToTCalib[j][k][l];
                   
                    
           
       //    cout<<"event " << event_number << "cOutput-> cbPlasE[j] " << cOutput-> cbPlasE[j] << endl;
            if(imphit.Energy>0 && cInputMain-> pbPlas_ToTCalib[j][k][l]>0  && j>15 && ImpIterator==0){
                 //cout<<"CORR PLAS " << event_number << " cInputMain-> pbPlas_ToTCalib[j][k][l] " <<cInputMain-> pbPlas_ToTCalib[j][k][l]<<" j " << j <<" imphit.Energy " << imphit.Energy<<  endl;
                 hA_gatedE_bPlas_E->Fill(cInputMain-> pbPlas_ToTCalib[j][k][l]);
                 hA_gatedE_bPlas_E_Ch[j][k]->Fill(cInputMain-> pbPlas_ToTCalib[j][k][l]);
         
                    }
                }
               }
            }
           }
             cOutput->cAIDAImplantE[ImpIterator] = imphit.Energy;

    
    
     for (auto& j : dechits)
      {
         AidaHit dechit = j;
      int DecstripX = dechit.StripX;
      int DecstripY = dechit.StripY;    
      int DecDSSD = dechit.DSSD-1;
     //cout<<"1 " << endl;
      //   if(DecstripX==ImpstripX && DecstripY==ImpstripY && DecDSSD==ImpDSSD){
             ///Prompt WR AIDA-bPlast gate 
      // if (dT_AIDA_bPlas>fCorrel->GAIDA_bPlas_TLow && dT_AIDA_bPlas<fCorrel->GAIDA_bPlas_THigh){
               
     ///AIDA Decay bPlas coincident Energy
    // cout<<"DecayIterator "<<DecayIterator << " Energy " << dechit.Energy << endl;
        if(dechit.Energy>0&& DecayIterator==0 && (AIDA_WR-bPLAS_WR)> 12000&&  (AIDA_WR-bPLAS_WR)<15000){
      
                   ///
                   for(int j=1; j<4; j++){
                        for(int k=0; k<16; k++) {
                         for(int l=0; l<10; l++) {
                if(j>15)  hADecay_bPlas_E->Fill(cInputMain-> pbPlas_ToTCalib[j][k][l]);
                            }
                        }
                   }
         }
         // cout<<"aida_imptime[DecDSSD][DecstripX][DecstripY] " << aida_imptime[DecDSSD][DecstripX][DecstripY] << " DecDSSD " << DecDSSD << " DecstripX "<< DecstripX <<" DecstripY " << DecstripY<<endl;
     if( dechit.Time>0 && aida_imptime[DecDSSD][DecstripX][DecstripY]>0 && bPlasfired==true && (AIDA_WR-bPLAS_WR)> 12000&&  (AIDA_WR-bPLAS_WR)<15000){
 
             AIDA_dT_imp_decay =  dechit.Time - aida_imptime[DecDSSD][DecstripX][DecstripY] ;///dT decay - implant time
            
           // cout<<"1 AIDA_dT_imp_decay " << AIDA_dT_imp_decay << endl;
            AIDA_dT_imp_decay_FRS_gated =  dechit.Time - aida_imptime_FRS_gated[DecDSSD][DecstripX][DecstripY];
           // cout<<"1 event_number " << event_number << " AIDA_dT_imp_decay " << AIDA_dT_imp_decay <<endl;
         
          
            cOutput->cAIDA_dT_imp_decay[cOutput->cAIDA_dT_imp_decay_hits] = AIDA_dT_imp_decay;
            cOutput->cAIDA_dT_imp_decay_FRS_gated[cOutput->cAIDA_dT_imp_decay_hits] = AIDA_dT_imp_decay_FRS_gated;
                   
            AIDA_dT_decay_bPlas = dechit.Time-bPLAS_WR;
            hADec_bPlas_dT->Fill(AIDA_dT_decay_bPlas); ///Decay time - bPlas Time
              
  
	 //   cout<<"  AIDA_dT_imp_decay*1E-6 "<<AIDA_dT_imp_decay*1E-6<<endl;
            hA_impdec_dT->Fill(AIDA_dT_imp_decay*1E-6); ///In ms
            hA_impdec_dT_FRS_gated->Fill(AIDA_dT_imp_decay_FRS_gated*1E-6); ///In ms
         //  }
            if(bPlasfired==true){
            hA_bPlas_decays_e[DecDSSD]->Fill(dechit.Energy); ///AIDA Decay Energy when bPlast Fires
            hA_bPlas_decays_strip_xy[DecDSSD]->Fill(DecstripX,DecstripY);
            }
   //  cout<<"event_number " << event_number <<" dT_AIDA_bPlas*ms " << dT_AIDA_bPlas*1E-6 <<" aida_imptime " << aida_imptime[DecDSSD][DecstripX][DecstripY] <<" DecstripX " <<DecstripX << " DecstripX " << DecstripX <<" dechit.Energy " << dechit.Energy << " AIDA_dT_imp_decay*1E-6 " <<AIDA_dT_imp_decay*1E-6 <<  " AIDA_dT_decay_bPlas ms " << AIDA_dT_decay_bPlas*1E-6 <<  endl;
             cOutput->cAIDA_dT_imp_decay_hits++;
         }
                    

            ///Reset the implantation pixel
              aida_imptime[DecDSSD][DecstripX][DecstripY] = 0;
              aida_imptime_FRS_gated[DecDSSD][DecstripX][DecstripY] = 0; 
       
              
            } ///END OF DECAY HIT EVENTS

                ImpIterator++;
              //cout<<"2ImpIterator " <<  ImpIterator << " imphit.Energy " << imphit.Energy << endl;
                cOutput->cImplantIterator =ImpIterator;
       } ///END OF IMPLANT HIT EVENTS
        
     //  cout<<"2event " << event_number <<  " ImpIterator " << ImpIterator <<endl;
//      if(cOutput->cImplantIterator>-1){
//         cout<<"cOutput->cImplantIterator " << cOutput->cImplantIterator <<  " ImpIterator " << ImpIterator <<endl;}
        ///For the decay events iterator
         for (auto& j : dechits)
      {
         
         AidaHit dechit = j;
        
        cOutput->cDecayIterator = DecayIterator;
    if(DecayIterator<1000)cOutput->cAIDADecayE[DecayIterator] = dechit.Energy;  
             DecayIterator++;               
            }
         
     
     AIDA_dT_imp_decay=0;
     AIDA_dT_imp_decay_FRS_gated=0;
    }
 }
 
  /**----------------------------------------------------------------------------------------------**/
  /**--------------- FRS-AIDA-bPlastic-Delayed FATIMA/GALILEO  (Isomers/Beta decay spectroscopy) ----------**/
  /**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Make_FRS_Delayed_AIDA_Gamma_Histos(){
     ///Germanium
     hA_Ge_WRdT = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Ge/Aida_Gal_WR_timediff","AIDA-Germanium WR Time difference ",4000,-100000,100000);
          
   //  hA_dT_imp_decay_vs_GeE  = MakeTH2('D',"Correlations/AIDA-BetaDelayed-Ge/Aida_Ge_E_vs.impdecdT","Ge Addback Energy vs. dT implantation-decay", 5000, 0, 5000, 4500,-90000,90000,"Ge Energy (keV)", "AIDA Implantation-Decay time (ms)");

     hA_dT_GeE = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Ge/Aida_Ge_E","AIDA Implant-Decay Time Ge Energy Addback Sum",5000,0,5000);
     
   //  hA_dT_FRS_Gated_imp_decay_vs_GeE  = MakeTH2('D',"Correlations/AIDA-BetaDelayed-Ge/Aida_FRS_Gated_Ge_E_vs.impdecdT","Ge Addback Energy vs. dT implantation-decay", 5000, 0, 5000, 4500,-90000,90000,"Ge Energy (keV)", "AIDA Implantation-Decay time (ms)");

     hA_dT_FRS_Gated_GeE = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Ge/Aida_FRS_Gated_Ge_E","AIDA Implant-Decay Time Ge Energy Addback Sum",5000,0,5000);
     
    ///Fatima
     hA_Fat_WRdT = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Fatima/Aida_Fat_WR_timediff","AIDA-Germanium WR Time difference ",4000,-100000,100000);
    
     hA_dT_imp_decay_vs_FatE  = MakeTH2('D',"Correlations/AIDA-BetaDelayed-Fatima/Aida_Fatima_E_vs.impdecdT","Fatima Energy vs. dT implantation-decay", 5000, 0, 5000, 4500,-90000,90000,"Ge Energy (keV)", "AIDA Implantation-Decay time (ms)");

     hA_dT_FatE = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Fatima/Aida_Fatima_E","AIDA Implant-Decay Time Fatima Energy Addback Sum",5000,0,5000);
     
     hA_dT_FRS_Gated_imp_decay_vs_FatE  = MakeTH2('D',"Correlations/AIDA-BetaDelayed-Ge/Aida_FRS_Gated_Fatima_E_vs.impdecdT","Fatima Addback Energy vs. dT implantation-decay", 5000, 0, 5000, 4500,-90000,90000,"Fatima Energy (keV)", "AIDA Implantation-Decay time (ms)");

     hA_dT_FRS_Gated_FatE = MakeTH1('I',"Correlations/AIDA-BetaDelayed-Fatima/Aida_FRS_Gated_Fatima_E","AIDA Implant-Decay Time Fatima Energy Addback Sum",5000,0,5000);
}
 /**-------------------FRS-AIDA-bPlastic-Delayed FATIMA/GALILEO  (Isomers/Beta decay spectroscopy)------*/

 void EventCorrelProc::Process_FRS_Delayed_AIDA_Gamma(EventAnlStore* cInput, EventCorrelStore* cOutput){
     if(AIDA_WR>0 && GAL_WR>0){
     dT_AIDA_Gal = AIDA_WR - GAL_WR;
     hA_Ge_WRdT->Fill(dT_AIDA_Gal);
    // cout<<"event_num " << event_number << " dT_AIDA_Gal " << dT_AIDA_Gal << " AIDA_WR " << AIDA_WR << " GAL_WR " << GAL_WR << endl;
     }
     
     if(AIDA_WR>0 && FAT_WR>0){ 
     dT_AIDA_Fat = AIDA_WR - FAT_WR;
     hA_Fat_WRdT->Fill(dT_AIDA_Fat);
     }
          
     
      
      for(int j =0; j<cOutput->cAIDA_dT_imp_decay_hits; j++){
           /// AIDA Fat 
            if(dT_AIDA_Fat>13000 && dT_AIDA_Fat<15000){
                   for(int k=0; k<cInput->pFatmult; k++){
                       // hA_dT_imp_decay_vs_FatE->Fill(cInput->pFat_QDC_E[cInput->pFat_QDC_ID[k]],cOutput->cAIDA_dT_imp_decay[j]*1E-6);
                        if(cInput->pFat_QDC_E[cInput->pFat_QDC_ID[k]]>0){
                        hA_dT_FatE->Fill(cInput->pFat_QDC_E[cInput->pFat_QDC_ID[k]]);       
                        }
                   }   
                }
   
                ///Germanium
                    for(int g=0; g<GALILEO_MAX_DETS; g++){
                        for (int h=0; h<GALILEO_CRYSTALS; h++){
              if( cOutput->cAIDA_dT_imp_decay[j]*1E-6>500 && cOutput->cAIDA_dT_imp_decay[j]*1E-6<24000 && j==1){
                  if(cInput->pGal_EAddback[g][h]>0  &&dT_AIDA_Gal>-1000 && dT_AIDA_Gal<500 ){
                        hA_dT_GeE->Fill(cInput->pGal_EAddback[g][h]);
                                /// AIDA Gal
                    ///disabled for now 02.12.19 AKM
                  // hA_dT_imp_decay_vs_GeE->Fill(cInput->pGal_EAddback,cOutput->cAIDA_dT_imp_decay[j]*1E-6);
                    } 
              }
                       ///FRS Gated AIDA-Gal
    if(cOutput->cAIDA_dT_imp_decay_FRS_gated[j]*1E-6>500 && cOutput->cAIDA_dT_imp_decay_FRS_gated[j]*1E-6<24000 && j==1&&dT_AIDA_Gal>-2000 && dT_AIDA_Gal<2000){
        ///disabled for now 02.12.19 AKM
      //   hA_dT_FRS_Gated_imp_decay_vs_GeE->Fill(cInput->pGal_EAddback,cOutput->cAIDA_dT_imp_decay_FRS_gated[j]*1E-6);
                        
                      if(cInput->pGal_EAddback[g][h]>0){
                         hA_dT_FRS_Gated_GeE->Fill(cInput->pGal_EAddback[g][h]);
                                }
                            }
                        }
                    }
                        ///FRS Gated AIDA-Fatima
                        for(int k=0; k<cInput->pFatmult; k++){
                      if(cInput->pFat_QDC_E[cInput->pFat_QDC_ID[k]]>0){
                        hA_dT_FRS_Gated_FatE->Fill(cInput->pFat_QDC_E[cInput->pFat_QDC_ID[k]]);
                            }
                        }
                   
                } ///End of FRS condition              
             }
          
        
 
    
 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------  FRS-AIDA-bPlastic-FATIMA  (Beta Decay Timing)  -------------**/
 /**----------------------------------------------------------------------------------------------**/
 void EventCorrelProc::Make_FRS_AIDA_bPlas_FATIMA_Histos(){}
// void EventCorrelProc::Process_FRS_AIDA_FATIMA(EventAnlStore* cInput, EventCorrelStore* cOutput){}
 
 /**----------------------------------------------------------------------------------------------**/
 /**------------------------------------- End of Correlations ------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
  TGo4WinCond* EventCorrelProc::MakeWindowCond(const char* fname,
                                           const char* cname,
                                           float left,
                                           float right,
                                           const char* HistoName) {
  // TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4WinCond*>(res);
   
   TGo4WinCond* cond = new TGo4WinCond((Text_t*)cname);
   cond->SetValues(left, right);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}
/**----------------------------------------------------------------------------------------------**/
 TGo4PolyCond* EventCorrelProc::MakePolyCond(const char* fname,
                                          const char* cname,
                                          Int_t size,
                                          Float_t (*points)[2],
                                          const char* HistoName) {
   //TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4PolyCond*>(res);
   
   Float_t fullx[size+1], fully[size+1];
   int numpoints = size;
   
   for (int i=0;i<numpoints;i++) {
     fullx[i] = points[i][0];
     fully[i] = points[i][1];
   }
   
   // connect first and last points
   if ((fullx[0]!=fullx[numpoints-1]) || (fully[0]!=fully[numpoints-1])) {
      fullx[numpoints] = fullx[0];
      fully[numpoints] = fully[0];
      numpoints++;
   }
 
   TCutG mycat("initialcut", numpoints, fullx, fully);
   TGo4PolyCond* cond = new TGo4PolyCond((Text_t*)cname);
   cond->SetValues(&mycat);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}

void EventCorrelProc::get_used_systems(){
    for(int i = 0;i < 6;i++) Used_Systems[i] = false;

  ifstream data("Configuration_Files/Used_Systems.txt");
  if(data.fail()){
    cerr << "Could not find Used_Systems config file!" << endl;
    exit(0);
  }
  int i = 0;
  int id = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d",s_tmp,&id);
    Used_Systems[i] = (id == 1);
    i++;
  }
  string DET_NAME[7] = {"FRS","AIDA","PLASTIC","FATIMA_VME","FATIMA_TAMEX","GALILEO","FINGER"};

    cout << "\n=====================================================" << endl;
    cout << "USED SYSTEMS" << endl;
    cout << "-----------------------------------------------------" << endl;
    for(int j = 0;j < 6;++j){
        if(Used_Systems[j]) cout << DET_NAME[j] << endl;
    }
    cout << "=====================================================" << endl;


}
/**----------------------------------------------------------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
