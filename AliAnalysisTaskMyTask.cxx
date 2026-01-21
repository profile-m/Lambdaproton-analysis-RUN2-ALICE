
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
#include "AliPIDResponse.h"   
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "TGraph.h"

#include "AliAnalysisTaskMyTask.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"   

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskMyTask.h"
#include "AliVEvent.h"

#include "AliPIDResponse.h"
#include "TH2F.h"
#include "AliPID.h"
#include "AliDetectorPID.h"

#include "TLorentzVector.h"

#include "TF1.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliESDpid.h"

#include "AliAODcascade.h"

#include "vector"
#include "TVector3.h"

//Constructor OR .h (Header File)


typedef std::vector<double> v;

Int_t NLambdas =0;
const double LambdaMass = 1.115; // GeV/c^2
const double ProtonMass = 0.938; // GeV/c^2
const double  PionMass = 0.135;
const double XiZeroMass = 1.3217; // GeV/c^2
double c = 2.99792457999999984e-02;


// NEW: Event buffer for mixing
static const int fMixDepth = 10; // Mix with 10 previous events
std::deque<std::vector<TLorentzVector>> fLambdaBuffer; //!
std::deque<std::vector<TLorentzVector>>  pp; //!fProtonBuffer;


class AliAnalysisTaskMyTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(),
    fAOD(0),fHistPtP(0),fHistEtaP(0),fHistPhiP(0),fHistPtPvsLength(0),
	fOutputList(0), fHistPt(0), fHistNEvents(0), fHistEta(0), fHistPhi(0), fHistMult(0),fHistNsigma(0), fHistdEdx(0), fHistDecayLengthVsPt (0), fHistDecayLength(0),  fHistInvMassLambda(0), fHistInvMassLambdap(0),
	fHistptvsalpha(0),// ArmenPodol
       fHistptvsalphawithcut(0),
	myfuncptvsalpha(0),

	fHistMassXiMinus(0),
	fHistMassXiPlus(0),
        fHistMassXiPlusVsPt(0),
	fHistMassXiMinusVsPt(0),

	fHistTrackOne(0), fHistDCAv0(0),
	fHistLambdasNbeforecutCosTheta(0),
	fHistCosThetaLambda(0),
	fHistCosThetaLambda_with_cut(0),
	fHistCosThetalambdaproton(0),
	fHistCosThetalambdaproton_with_cut(0),
        fHistPtLambda(0),
        fHistPtlambdaproton(0),
	fHistEtaLambda(0),
	fHistPhiLambda(0),
	fHistEtalambdaproton(0),
	fHistPhilambdaproton(0),
	fHistOpenAngle(0),
       cos0(0),cos1(0),cos10(0),cos50(0),

	fPIDResponse(0)
{
   // default constructor, don't allocate memory here!
   // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0),fHistPtP(0),fHistEtaP(0),fHistPhiP(0),fHistPtPvsLength(0),
	fOutputList(0), fHistPt(0), fHistNEvents(0), fHistEta(0), fHistPhi(0), fHistMult(0), fHistNsigma(0), fHistdEdx(0),  fHistDecayLengthVsPt (0), fHistDecayLength(0), fHistInvMassLambda(0), fHistInvMassLambdap(0),
    fHistptvsalpha(0),
	fHistptvsalphawithcut(0),

	fHistMassXiPlus(0),
	fHistMassXiMinus(0),
	fHistMassXiPlusVsPt(0),
	fHistMassXiMinusVsPt(0),

	fHistOpenAngle(0), fHistDCAv0(0),
	fHistLambdasNbeforecutCosTheta(0),
	fHistCosThetaLambda(0),
	fHistCosThetaLambda_with_cut(0),
	fHistCosThetalambdaproton(0),
	fHistCosThetalambdaproton_with_cut(0),
        fHistPtLambda(0),
	fHistPtlambdaproton(0),
	fHistEtaLambda(0),
	fHistPhiLambda(0),
	fHistEtalambdaproton(0),
	fHistPhilambdaproton(0),
    cos0(0),cos1(0),cos10(0),cos50(0),
       
	fPIDResponse(0)
{

    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{

    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}

//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

    // example of a histogram

    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogram
    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 1, 0, 1);       // create your histogram
    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -10, 10);       // create your histogram
    fHistPhi = new TH1F("fHistPhi", "fHistPhi", 100, 0, 10);       // create your histogram
    fHistMult = new TH1F("fHistMult", "fHistMult", 100, 0, 1000);       // create your histogram, was 10000

    fHistdEdx      = new TH2F("fHistdEdx","Momentum vs dEdx;p (GeV/#it{c});dE/dx (arb. units)",500,0,10,1000,0,5000);
    fHistNsigma      = new TH2F("fHistNsigma","Momentum vs N sigma TPC Proton;p (GeV/#it{c});n_{#sigma,TPC}",500,-15,15,500,-20,20);

    fHistLambdasNbeforecutCosTheta= new TH1F("fHistLambdasNbeforecutCosTheta", "N LAmbdas",1,0,1);
    fHistInvMassLambdap = new TH1F("fHistInvMassLambdap", "Invariant Mass of Lambdaproton; M_lambda{p} (GeV/c^2); Counts", 100, 2.0, 3.0); // 1000
    fHistInvMassLambda = new TH1F("fHistInvMassLambda", "Invariant Mass of Lambda; M_{Lambda} (GeV/c^2); Counts", 50, 1.0, 1.4); //10000

    fHistDecayLengthVsPt = new TH2F("fHistDecayLengthVsPt", "Decay Length vs pT; decay cm (pT (GeV/c); pT GeV/c ( Decay Length (cm)", 500, 0, 10, 500, 0, 16); //500

    fHistDecayLength = new TH1F ("fHistDecayLength", "Decay of Lambda(cm)", 500, 0, 5);

    fHistTrackOne = new TH1F ("fHistTrackPne", "Daugth one mass; Gev/c^2", 500, -3, 3);
    fHistDCAv0 =  new TH1F ("fHistDCAv0", "dca lambda daughters/cm", 500, -2,5);


     fHistPtP=new TH1F("fHistPtP","Momentum of proton",100, 0, 20);
     fHistEtaP = new TH1F("fHistEtaP", "Eta proton", 100, -5, 5.0);
     fHistPhiP = new TH1F("fHistPhiP", "Phi P", 100, -5, 10.0);

     fHistPtPvsLength = new TH2F ("fHistPtPvsLength", "Length vs Pt of proton, Gev/c*{#cm}", 500, 0 ,10,200, 0,16);

    fHistptvsalpha =  new TH2F("fHistptvsalpha", "ArmenterosPodolansi diagram, alpha vs pt", 200, -1.0,1.0,500, 0,1);

    fHistptvsalphawithcut =  new TH2F("fHistptvsalphawithcut", "ArmenterosPodolansi diagram with all cuts+ AP cuts, alpha vs pt", 200, -1.0,1.0,500, 0,1);


    fHistPtLambda = new TH1F("fHistPtLambda", "Pt Momentum lambda, Gev/c^2", 100, 0, 16);
	fHistPtlambdaproton =  new TH1F ("fHistPtlambdaproton", "Pt lambda proton, Gev/c^2", 100, 0, 16);
	fHistCosThetaLambda = new TH1F ("fHistCosThetaLambda", "CosTheta lambda",100, -1, 1);
	fHistCosThetaLambda_with_cut = new TH1F ("fHistCosThetaLambda_with_cut", "Cos Theta with cut", 100, -1, 1);
	fHistCosThetalambdaproton = new TH1F( "fHistCosThetalambdaproton", "CosTheta lambda proton",100, -1, 1);
	fHistCosThetalambdaproton_with_cut = new TH1F ("fHistCosThetalambdaproton_with_cut", "Cos Theta lambda-proton with cut", 100, -1,1 );
	fHistEtaLambda = new TH1F("fHistEtaLambda", "Eta lambda", 100, -5, 5);
	fHistPhiLambda = new TH1F("fHistPhiLambda", "phi lambda", 100, 0, 10);

	fHistEtalambdaproton = new TH1F ("fHistEtalambdaproton", "eta lambda proton", 100, -5, 5);
	fHistPhilambdaproton = new TH1F ("fHistPhilambdaproton", "phi lambda proton", 100, 0, 6);
	       // cascades
	    fHistMassXiPlusVsPt   = new TH2F  ("fHistMassXiMinusVsPt"," mass vs pt ", 200, 0, 3, 500, 0, 10 );
	    fHistMassXiMinusVsPt  = new TH2F  ("fHistMassXiPlusVsPt", "mass vs pt", 500,0 , 4, 500, 0 , 10 );

	   fHistMassXiMinus  = new TH1F ("fHistMassXiMinus", "plus", 500, 0, 4) ;
	   fHistMassXiPlus = new TH1F( "fHistMassXiPlus", "xi massi", 500, 0, 4 );

	   //
	fHistOpenAngle = new TH1F("fHistOpenAngle","Open Angle v0's", 100, 0,2.0);


     cos0 = new TH1F("cos0", "Invariant Mass of Lambdaproton; M_lambda{p} (GeV/c^2); Counts", 100, 2.0, 4.0);
    cos1 = new TH1F("cos1", "Invariant Mass of Lambda; M_{Lambdap} (GeV/c^2); Counts", 50, 2.0, 5);
    cos10 = new TH1F("cos10", "Invariant Mass of Lambda; M_{Lambdap} (GeV/c^2); Counts", 50, 2, 5); //10000
    cos50 = new TH1F("cos50", "Invariant Mass of Lambda; M_{Lambdap} (GeV/c^2); Counts", 50, 2.0, 5);

cos0_mixed = new TH1F("cos0_mixed", "Mixed M(#Lambdap) cos<-0.5; M (GeV/c^{2}); Counts", 100, 2.0, 4.0);
cos1_mixed = new TH1F("cos1_mixed", "Mixed M(#Lambdap) -0.5<cos<0; M (GeV/c^{2}); Counts", 100, 2.0, 4.0);
cos10_mixed = new TH1F("cos10_mixed", "Mixed M(#Lambdap) 0<cos<0.5; M (GeV/c^{2}); Counts", 100, 2.0, 4.0);
cos50_mixed = new TH1F("cos50_mixed", "Mixed M(#Lambdap) cos>0.5; M (GeV/c^{2}); Counts", 100, 2.0, 4.0);
fHistMass_Mixed = new TH1F("fHistMass_Mixed", "Invariant Mass of mixed for Lambda ; M_{Lambda} (GeV/c^2); Counts", 100, 0.3, 2.0); //10000


	//set TStyle
	fHistptvsalpha->GetXaxis()->SetTitle("#alpha");
	fHistptvsalpha->GetYaxis()->SetTitle("q_{t}");
	fHistptvsalpha->SetMarkerStyle(kFullCircle);


              // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    fOutputList->Add(fHistNEvents);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistPhi);
    fOutputList->Add(fHistMult);

    fOutputList->Add(fHistdEdx);
    fOutputList->Add(fHistNsigma);

    fOutputList-> Add(fHistInvMassLambdap);// Lambda = Xi0
    fOutputList-> Add(fHistInvMassLambda);

    fOutputList-> Add(fHistDecayLengthVsPt);
    fOutputList->Add(fHistDecayLength);

      fOutputList -> Add (fHistMassXiMinus);
       fOutputList -> Add (fHistMassXiPlus);
	 fOutputList -> Add (fHistMassXiPlusVsPt);
	 fOutputList-> Add( fHistMassXiMinusVsPt);


    fOutputList->Add(fHistPtP);
    fOutputList->Add(fHistDCAv0);
    fOutputList->Add(fHistPhiP);
    fOutputList->Add(fHistEtaP);
    fOutputList->Add(fHistPtPvsLength);

    fOutputList->Add(fHistptvsalpha);
    //get pt and costheta lambda & lambda proton
    fOutputList->Add( fHistLambdasNbeforecutCosTheta);
    fOutputList-> Add(fHistPtLambda);
    fOutputList->Add( fHistPtlambdaproton);

    fOutputList->Add(fHistCosThetaLambda);
    fOutputList->Add(fHistCosThetaLambda_with_cut);
    fOutputList->Add(fHistCosThetalambdaproton);
    fOutputList->Add(fHistCosThetalambdaproton_with_cut);
    fOutputList->Add(fHistEtaLambda);
    fOutputList->Add(fHistPhiLambda);

    fOutputList->Add(fHistEtalambdaproton);
    fOutputList->Add(fHistPhilambdaproton);

    fOutputList->Add(fHistOpenAngle);

    fOutputList->Add(cos0);
    fOutputList->Add(cos1);
    fOutputList->Add(cos50);
    fOutputList->Add(cos10);
    fOutputList->Add(fHistMass_Mixed);
     fOutputList->Add(fHistptvsalphawithcuts);
                  
    // fOutputList->Add(fHistLambdaNBeforecut);
     PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________


void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();//Get Analysis Manager
    	if(man)  		{
      	 	  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());  //Get Input Handler
      		  if(inputHandler) fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());//inputHandler->GetPIDResponse();    //Retrieve the AliPIDResponse object from the analysis manager
      		 else AliWarning("No Input Handler!");
	}
    	else AliWarning("No Analysis Manager!");
	//Get Primary Vertex

 // here mulpiplicity 

	AliAODVertex *PrimaryVertex = fAOD->GetPrimaryVertex();//cut error
  if(!PrimaryVertex)::Warning("AliAnalsisTaskMyTask::UserExec","No AliAODVertex object found!");


  Double_t  fMag =fAOD->GetMagneticField();
 //get Primery vertex of track and put the cut
 //
 //
 Int_t  nTracks   =fAOD->GetNumberOfTracks();
 Int_t  nV0s      =fAOD->GetNumberOfV0s();
 Int_t  nCascades =fAOD->GetNumberOfCascades();

 //get initial primary vertex for secondary vertex using KF
 /* vertex       = fAOD->GetPrimaryVertexSPD();
  primVertex   = AliAnalysisTaskpOmegaDibaryon::AODToESDVertex(*vertex);
  primKFVertex  = CreateKFVertex(*primVertex);
  PrimVertexKF[0]   = primKFVertex.GetX();
  PrimVertexKF[1]   = primKFVertex.GetY();
  PrimVertexKF[2]   = primKFVertex.GetZ();
*/


  Double_t PrimaryVertexPos[3] = {-999.0,-999.0,-999.0};// array ?
 PrimaryVertex->GetXYZ(PrimaryVertexPos);
  double PrimaryVertexZ = PrimaryVertexPos[2]; // 3rd el
  //cm
  // error check
  if(TMath::IsNaN(PrimaryVertexZ)) return;
 // if(TMath::Abs(PrimaryVertexZ) > 10.) return;


  Int_t iTracks(fAOD->GetNumberOfTracks());  // see how many tracks there are in the event



    for(Int_t i(0); i < iTracks; i++) {   // loop ove rall these tracks


	AliAODTrack* Track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i)); // get a track (type AliAODTrack) from the event
      	if(!Track) continue;                  // if we failed, skip this track

	 float xv[2];
   	 float yv[3];

   	 Track->GetImpactParameters(xv,yv);


    	float DCAxy = xv[0];
   	double TOF_m2_nSigma    = -999.0;

    	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,Track);

    	bool TOFisOK = false;
    	if(statusTOF == AliPIDResponse::kDetPidOk) TOFisOK = true;

    	double ITS_dEdx	    = -999.0;
    	double ITS_dEdx_nSigma  = -999.0;

    	int ITS_nCluster	    = 0;
   	 AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,Track);

   	 bool ITSisOK = false;
    if(statusITS == AliPIDResponse::kDetPidOk) ITSisOK = true;

   	 Double_t dEdx   = Track->GetTPCsignal();
   	 Double_t p    = Track->GetTPCmomentum();

    Double_t len = Track->GetIntegratedLength();
   	Double_t tim = 0;
    Double_t beta = -1;

	//get its_momentum
	//fHistITSSIGNAL->Fill(p,p_its);
	//
    	/*if(TOFisOK){ if(Track->GetTPCmomentum()!=0 && len>=350.0) tim = Track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(Track->GetTPCmomentum());
 		else  tim = Track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(TMath::Sqrt(Track->Px()*Track->Px()+Track->Py()*Track->Py()+Track->Pz()*Track->Pz()));
    	if(tim != 0.) beta = len / (tim * c);}*/

	 Double_t nSigmaTPCProt = -99999.0;
     	 Double_t nSigmaTOFProt = -99999.0;

//    	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,Track)==AliPIDResponse::kDetPidOk)nSigmaTPCProt = fPIDResponse->NumberOfSigmasTPC(Track,AliPID::kProton);
  // 	if (p!=0 && tim!=0 && len>=350.) nSigmaTOFProt = fPIDResponse->NumberOfSigmasTOF(Track,AliPID::kProton);

	fHistdEdx->Fill(p,dEdx);
	fHistNsigma->Fill(p,nSigmaTPCProt);

        fHistPt->Fill(Track->Pt());    // plot the pt value of the track in a histogram
        fHistEta->Fill(Track->Eta());
        fHistPhi->Fill(Track->Phi());


	//General cuts chi2 tpc eta //
	if (Track->Chi2perNDF()>=4)continue;
	if( Track->GetTPCNcls()<=70)continue; // its~ like length of flight parameter. //    if (Track->GetITS()<30)continue
	if ((TMath::Abs(Track->Eta() >= 0.8)))continue;


	//get the Proton vector for Xip decay
        // proton selection
	AliAODTrack* protonTrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
	// cut on sigma TPC
	if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(protonTrack,AliPID::kProton))>3)continue;


       //get Cascades
		int nCAsc = fAOD -> GetNumberOfCascades();

		Double_t arrayvertex [3]= {0.};

   	arrayvertex[0] = fAOD -> GetPrimaryVertex()->GetX();
		arrayvertex[1] = fAOD -> GetPrimaryVertex()->GetY();
		arrayvertex[2] = fAOD -> GetPrimaryVertex()->GetZ();

	for ( Int_t n = 0; n < nCAsc ; n++) {

		AliAODcascade* cascade = fAOD->GetCascade(n);
		if(!cascade) continue;
		Double_t Casc_Phi = cascade -> Phi();
		Double_t Casc_Eta = cascade -> Eta();
		Double_t Casc_Pt = cascade -> Pt();
		Double_t MassXiMinus = cascade -> MassXi();
		Double_t MassXiPlus = cascade -> MassXi();
		Short_t charge = cascade->ChargeXi();
		    if ( charge <0 ){

			 fHistMassXiMinus -> Fill (MassXiMinus ) ;
	     fHistMassXiMinusVsPt -> Fill (MassXiMinus, Casc_Pt);
		    } if(charge >0) {
			fHistMassXiPlus ->Fill (MassXiPlus);
	    fHistMassXiPlusVsPt ->Fill (MassXiPlus, Casc_Pt);

		}
	}

	// create TLorentz Vector
   TLorentzVector ProtonVector;
	ProtonVector.SetPxPyPzE( protonTrack->Px(), protonTrack->Py(), protonTrack->Pz(),TMath::Sqrt(ProtonMass*ProtonMass+protonTrack->P()*protonTrack->P()));

// get parameters from the Track Proton (lambda-proton)
	//
        Double_t Ptproton = protonTrack->Pt();
        //if(Ptproton<1.0)continue;
        Double_t Lenproton = protonTrack->GetIntegratedLength();

	//dynamic_cast<TH1F*>(fOutputList -> FindObject("fHistPt"))->Fill(Pt);

	fHistPtP->Fill(Ptproton);
    fHistPtPvsLength->Fill( Lenproton, Ptproton);
	fHistEtaP->Fill(protonTrack->Eta());
    fHistPhiP->Fill(protonTrack->Phi());

      //buffer proton

     std::vector <TLorentzVector> p;

     p.pushback(ProtonVector);



	//check v0 methods and AODTracks;
	
//for(Int_t
for (Int_t j = 0; j < fAOD->GetNumberOfV0s(); j++){

	       	AliAODv0* v0 = fAOD->GetV0(j);
	    	if(v0->CosPointingAngle(PrimaryVertex)<=0.999)continue;
		 //
	   	 Double_t OpenAnglev0;
	   	 OpenAnglev0 = v0->OpenAngleV0();
	// check the full range of the OpenAngle
	         fHistOpenAngle->Fill(OpenAnglev0);
	         if(OpenAnglev0 >0.8) continue;// check 0.9 , was >0.8 cos of angle < ~1->0 grad

	         //check OnFlyStatus -> checke redused background nicely
		 ///
		 Bool_t FLy = v0->GetOnFlyStatus();
	   	 if (!FLy)continue;
		// get the Lambda Daughters : p an n
	 	 AliAODTrack* Trackone = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
      	 	 AliAODTrack* Tracktwo = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));

 	 	// set the DCA daughters p and n and lambda cut
		 if(v0->DcaV0Daughters()>0.1)continue;
 		 if(v0->DcaV0ToPrimVertex()>0.05)continue;
  		//DCA histograms
		 Double_t DcaV0 = v0->DcaV0ToPrimVertex();
	 	 fHistDCAv0->Fill(DcaV0);
  		//get TLorentz Vectors for Lambda, daughters, lambdaproton
  		//
		 TLorentzVector lambdaVector,  LambdaProtonVector;
		 TLorentzVector posProtonVector;
	 	 TLorentzVector negPionVector;

         if (( Trackone->Charge())< 0) continue ;
		//set cuts Sigmas kProton and kPion and their Charges
        if((( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Trackone,AliPID::kProton))) > 3)&&(Trackone->Charge()< 0))continue;// pro =true;
	  	 
        if (( Tracktwo->Charge())> 0)continue;
         if((( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Tracktwo, AliPID::kPion))) > 3)&&(Tracktwo->Charge()>0))continue;//  pio = true;
	 	
         // get the components of P of Proton and Pions track
	 	posProtonVector.SetPxPyPzE( Trackone->Px(),Trackone->Py(), Trackone->Pz(),TMath::Sqrt(ProtonMass*ProtonMass+Trackone->P()*Trackone->P()));
	  	negPionVector.SetPxPyPzE(Tracktwo->Px(),Tracktwo->Py(),Tracktwo->Pz(),TMath::Sqrt(PionMass*PionMass+Tracktwo->P()*Tracktwo->P()));



	  //  if ((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton)) > 3)||(TMath::Abs(fPIDResponse->NumberofSigmasTPC(posTrack, AliPID::kPion)))>3) continue;
	  //  if (negTrack->Charge()>0)continue;
	  //  if (posTrack->Charge()<0)continue;
	  //  Double_t TrackoneMass = Trackone ->GetMass();
	  //  Double_t TracktwoMass = Tracktwo ->GetMass();
	  //  if(((TMath::Abs(fPIDResponce->NumberOfSigmasTPC(Trackone, AliPID::kPion)))||(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Tracktwo, AliPID::kPion)))){


	    //set lambda vector
        lambdaVector = posProtonVector + negPionVector; // combine inv lambda mass
	     fHistLambdasNbeforecutCosTheta->Fill(0);

	    //Cos Theta angle for Lambda's
	     fHistCosThetaLambda->Fill(lambdaVector.CosTheta());

	     // CosTHeta Lambda cut
	     if(lambdaVector.CosTheta()>0.5)continue;
	     fHistCosThetaLambda_with_cut->Fill(lambdaVector.CosTheta());

	     fHistPtLambda->Fill(lambdaVector.Pt());
	     //fHistInvMassLambda->Fill(lambdaVector.M());


	  	  //get impulse from TLoretz vector
	  	  //

	  	 //GEneral variables for LAmbdas'
	    fHistEtaLambda->Fill(lambdaVector.Eta());

	   // check if it depandce of varibles or not + Compare NEnties for MLambdas amd Eta+check the lamda and lambdaproton distributions
	    //
	    //

	    fHistInvMassLambda->Fill(lambdaVector.M());


           //set Eta Phi Lambda in a vector just in case and read Phi
	    std::vector<double> letaphi= {lambdaVector.Eta(),lambdaVector.Phi()};
            fHistPhiLambda->Fill(letaphi[1]);

	    // ArmentosPodolanski
	    v l = {lambdaVector.Px(),lambdaVector.Py(),lambdaVector.Pz()};
	    v t1 = {posProtonVector.Px(),posProtonVector.Py(), posProtonVector.Pz()};
	    v t2 = {negPionVector.Px(), negPionVector.Py(), negPionVector.Pz()};

	    TVector3 vecN(t1[0],t1[1],t1[2]);
 	    TVector3 vecP(t2[0],t2[1],t2[2]);
    	   TVector3 vecM(l[0],l[1],l[2]);

   	 Double_t thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
   	 Double_t thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
   	 Double_t alfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/
     	 ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) ;
    	Double_t qt = vecP.Mag()*sin(thetaP);
    	fHistptvsalpha->Fill(alfa,qt); // Armennteros-Podolanski
    	//
	//
        //
        // K0's cut
         if (qt>0.2*alfa)continue;// smash K0s
	// lambda proton vector
       fHistptvsalphawithcuts->Fill(alfa,qt); // Armennteros-Podolansk


         std::vector <TLorentzVector> vv;
         vv.pushback(lambdaVector);


    	LambdaProtonVector = lambdaVector + ProtonVector; // combine lambda p for ksip look
	    //get eta phi distributions for lambda-proton
	    //
	    fHistEtalambdaproton->Fill(LambdaProtonVector.Eta());
	    std::vector<double> lpphi = {LambdaProtonVector.Eta(),LambdaProtonVector.Phi()};
	    // fHist Eta Phi lambda proton
	    fHistPhilambdaproton->Fill(lpphi[1]);

	    //get Pt for lambda-proton distridution
	    fHistPtlambdaproton->Fill(LambdaProtonVector.Pt());
	    // get the CosTheta distribution for lambda-proton (after lambda cos theta cut)
	    fHistCosThetalambdaproton->Fill(LambdaProtonVector.CosTheta());// cheched

	    //cut Costheta lambdaproton


	   // if (LambdaProtonVector.CosTheta()>0.8)continue;
	    // Draw CosTheta lambda proton with cut
	    fHistCosThetalambdaproton_with_cut->Fill(LambdaProtonVector.CosTheta());

	    // check DCA-> DCA small in ranges
	    fHistDCAv0->Fill(DcaV0);
	    // lambda mass
	    fHistInvMassLambda->Fill(lambdaVector.M());
	    // lambda p
            fHistInvMassLambdap->Fill(LambdaProtonVector.M());
	    // proton mass
	    fHistTrackOne->Fill(posProtonVector.M());
            Double_t decayLength = 0 ;
	    decayLength = (lambdaVector.Vect().Mag() / lambdaVector.P()) * lambdaVector.M();
	    // fill decay legth and decayLength versus transverse momentum  histogram
	    fHistDecayLength->Fill(decayLength);
        fHistDecayLengthVsPt->Fill(decayLength,lambdaVector.Pt());
   

   Double_t cosTheta = LambdaProtonVector.CosTheta();
     if(cosTheta < -0.5) {
    cos0->Fill(LambdaProtonVector.M());}

if((cosTheta > -0.5) && (cosTheta < 0)) {
    cos1->Fill(LambdaProtonVector.M());
}

// Bin 3: 0 < cosθ < 0.5 (slightly forward)
        if((cosTheta >0) && (cosTheta < 0.5)) {
        cos10->Fill(LambdaProtonVector.M());}

// Bin 4: cosθ > 0.5 (forward)
        if(cosTheta > 0.5) {
        cos50->Fill(LambdaProtonVector.M());}
        


}
}

     // mix buffer proton lambda
       for (auto& vi: vv ){  //
        for (auto& ppi: pp) {   //
            for(auto& pi: p){

                TLorentzVector mixedPair = vi + pi; //lambda+proton
                Double_t mixedMass = mixedPair.M();
                Double_t mixedCosTheta = mixedPair.CosTheta();


if(mixedCosTheta < -0.5) {
 cos0_mixed->Fill(mixedMass); }
if(mixedCosTheta > -0.5 && mixedCosTheta < 0) {
 cos1_mixed->Fill(mixedMass); }
if(mixedCosTheta > 0 && mixedCosTheta < 0.5) {
cos10_mixed->Fill(mixedMass);
 }
if(mixedCosTheta > 0.5) {
cos50_mixed->Fill(mixedMass);
        }

       }
   }
}
    pp.pushback(p);
    if (pp.size()>10) pp.pop_front();

    //std::cout<<NLambdas<<std::endl;
    fHistNEvents->Fill(0);
    fHistMult->Fill(iTracks);
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}

//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
