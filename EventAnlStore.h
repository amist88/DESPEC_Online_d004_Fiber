// $Id: EventAnlStore.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef TSCNCALEVENT_H
#define TSCNCALEVENT_H

#include "TGo4EventElement.h"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "AIDA_Event.h"

  
struct AidaAnlData {
      std::vector<AidaHit> Implants;
      std::vector<AidaHit> Decays;
      };
class EventAnlStore : public TGo4EventElement {
   public:
      EventAnlStore() : TGo4EventElement() {}
      EventAnlStore(const char* name) : TGo4EventElement(name) {}
      virtual ~EventAnlStore() {}

      virtual void  Clear(Option_t *t="");

      ///General Outputs 
      Int_t pTrigger;
      Int_t pEvent_Number;
      Int_t pPrcID_Conv[7];
      bool pOnSpill;
      
      ///White Rabbit Outputs
      Long64_t pFRS_WR;
      Long64_t pbPLAS_WR;
      Long64_t pFAT_WR;
      Long64_t pFAT_Tamex_WR;
      Long64_t pAIDA_WR;
      Long64_t pGAL_WR;
      ///FRS Outputs
    
      Float_t pFRS_AoQ;
      Float_t pFRS_ID_x2, pFRS_ID_y2, pFRS_ID_a2, pFRS_ID_b2;
      Float_t pFRS_ID_x4, pFRS_ID_y4, pFRS_ID_a4, pFRS_ID_b4;
      Float_t pFRS_z;
      Float_t pFRS_z2;
      Float_t pFRS_dEdeg;
      Float_t pFRS_dEdegoQ;
      Int_t   pSci_num;
      Float_t pFRS_sci_l[12];
      Float_t pFRS_sci_r[12];
      Float_t pFRS_sci_e[12];
      Float_t pFRS_sci_tx[12];
      Float_t pFRS_sci_x[12];
      Float_t pFRS_tpc_x[7];
      Float_t pFRS_tpc_y[7];
      UInt_t pFRS_scaler[64];
      UInt_t pFRS_scaler_delta[64];
      
      Double_t pTRaw_vftx_21l;
      Double_t pTRaw_vftx_21r;
      Double_t pTRaw_vftx_22l;
      Double_t pTRaw_vftx_22r;
      Double_t pTRaw_vftx_41l;
      Double_t pTRaw_vftx_41r;
      Double_t pTRaw_vftx_42l;
      Double_t pTRaw_vftx_42r;
      
      Bool_t pFRS_ZAoQ_pass[8];
      Bool_t pFRS_x2AoQ_pass[8];
      Bool_t pFRS_x4AoQ_pass[8];
      Bool_t pFRS_Z_Z2_pass[8];
      
     ///AIDA output 
       //AidaAnlData pAida;
      std::vector<AidaAnlData> pAida;
      
      Int_t         pFat_QDC_ID[50];
      Double_t      pFat_QDC_E[50];
      Double_t      pFat_QDC_E_Raw[50];
      Long64_t      pFat_QDC_T_coarse[50];
      Double_t      pFat_QDC_T_fine[50];
      Int_t         pFat_TDC_ID[50];
      Long64_t      pFat_TDC_T[50];
      Long64_t      pFat_TDC_T_Raw[50];
    
      Int_t         pFatmult;
      Int_t         pSC40mult;
      Int_t         pSC41mult;
            
      Int_t pFat_TDC_Singles_ID[50];
      Long64_t pFat_TDC_Singles_t[50]; 
      Long64_t pFat_TDC_Singles_t_Raw[50]; 
      Long64_t pFat_QDC_Singles_t_coarse[50]; 
      double_t pFat_QDC_Singles_t_fine[50]; 

      Int_t pFat_QDC_Singles_ID[50];
      Double_t pFat_QDC_Singles_E[50];
      Double_t pFat_QDC_Singles_E_Raw[50];
      UInt_t pstdcmult;
      UInt_t psqdcmult; 
      
      //SC 40 and 41 arrays
      Long64_t pSC40[10];
      Long64_t pSC41[10];

      Double_t pFat_ToTCalib[50][10];
      Double_t pFat_LeadT[50][10];
      Double_t pFat_TrailT[50][10];
      Int_t    pFat_LeadHits;
     // Int_t    pFat_TrailHits;
      Int_t    pFat_Tamex_chan;
     
      Double_t pbPlas_ToTCalib[4][16][10];
      Int_t    pbPlas_PMT_Lead_N[4][16];
      Double_t pbPlas_LeadT[4][16][10];
      Double_t pbPlas_LeadT_Avg;
      Int_t    pbPlas_PMT_Trail_N[4][16];
      Double_t pbPlas_TrailT[4][16][10];
      Int_t    pbPlas_LeadHits;
      Int_t    pbPlas_TrailHits;
      
          ULong64_t pGal_Event_T[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
          ULong64_t pGal_T[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
          double   pGal_E[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
          double   pGal_E_Raw[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
          double   pGal_EAddback[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
          bool     pGalPileUp[GALILEO_MAX_HITS];
          bool     pGalOverFlow[GALILEO_MAX_HITS];
         
      

//       Int_t pFing_tot;
//     //  Double_t pFing_lead_lead[52];
//       Double_t pFing_Lead_Up[52][100];
//       Double_t pFing_Lead_Down[52][100];
//       Double_t pFing_Trail_Up[52][100];
//       Double_t pFing_Trail_Down[52][100];   
//       Int_t pFing_stripID;
//       Double_t pFing_maxtot;
//       Int_t pFing_maxtotchan;

   ClassDef(EventAnlStore,1)
};
#endif //TSCNCALEVENT_H



