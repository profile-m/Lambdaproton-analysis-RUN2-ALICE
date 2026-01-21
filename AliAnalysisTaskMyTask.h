#ifndef AliAnalysisTaskMyTask_H
#define AliAnalysisTaskMyTask_H

#include "AliAODVertex.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"   
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TF1.h"


class AliAODEvent;

    class AliAnalysisTaskMyTask : public AliAnalysisTaskSE  
{	
	
	public:
                                AliAnalysisTaskMyTask();
                                AliAnalysisTaskMyTask(const char *name);
	
        virtual                 ~AliAnalysisTaskMyTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
	
	AliAODEvent* 		fAOD;
    	TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram
        TH1F*                   fHistNEvents;        //! dummy histogram
        TH1F*                   fHistEta;        //! dummy histogram
        TH1F*                   fHistPhi;        //! dummy histogram

        TH1F*                   fHistMult;        //! dummy histogram
     
        
    	TH2F*                   fHistNsigma;       
        TH2F*                   fHistdEdx; 

        TH2F* 			fHistDecayLengthVsPt; 

	   TH2F* fHistdEdxVsP;
	   TH1F* fHistCheck;    
	   TH1F* fHistDCAv0;
	   TH1F* fHistOpenAngle;
	   TH1F*  fHistPtP;
	   TH1F*  fHistEtaP;
	   TH1F*  fHistPhiP;
	   TH1F*  fHistPx;
	   TH1F*  fHistPy;
	   TH1F*  fHistPz;
	   TH2F* fHistPtPvsLength;
	   TH1F* fHistInvMassLambda; 
	   TH1F* fHistPtLambda;
	   TH1F* fHistEtaLambda; 
	   TH1F* fHistPhiLambda; 
	   TH1F* fHistTrackOne;//mass of 1st daughter
	   TH1F* fHistCosThetaLambda; 
	   TH1F* fHistCosThetaLambda_with_cut;
	   TH1F* fHistLambdasNbeforecutCosTheta;
	   TH1F* fHistInvMassLambdap;
	   TH1F* fHistPtlambdaproton;
	   TH1F* fHistEtalambdaproton; 
	   TH1F* fHistPhilambdaproton; 
	   TH1F* fHistCosThetalambdaproton; 
	   TH1F* fHistCosThetalambdaproton_with_cut;
	   TH1F* fHistDecayLength;
	
	   TH1F* fHistMassXiPlus;
	   TH1F*fHistMassXiMinus;
	   TH2F* fHistMassXiPlusVsPt;
	   TH2F*fHistMassXiMinusVsPt;
	 
	   TH2F*  fHistptvsalpha;
  
		AliPIDResponse *        fPIDResponse;
	
	     AliAnalysisTaskMyTask(const AliAnalysisTaskMyTask&); // not implemented
	        AliAnalysisTaskMyTask& operator=(const AliAnalysisTaskMyTask&); // not implemented
	        ClassDef(AliAnalysisTaskMyTask, 1);


};


#endif






