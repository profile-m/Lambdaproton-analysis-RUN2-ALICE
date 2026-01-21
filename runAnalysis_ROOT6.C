#include "AliAnalysisTaskMyTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"   

#include "AliAnalysisAlien.h"   
#include "AliAODCluster.h"
#include "AliAODBaseMult.h"   

void runAnalysis_ROOT6()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
   // Bool_t local = kTRUE;
    Bool_t local = kFALSE;   
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kTRUE;//TRUE;
    //Bool_t gridTest = k;//TRUE;
    Bool_t isMC = kFALSE;   

gSystem->Load("libCore");
	  gSystem->Load("libGeom");
	  gSystem->Load("libVMC");
	  gSystem->Load("libPhysics");
	  gSystem->Load("libTree");
	  gSystem->Load("libSTEERBase");
	  gSystem->Load("libESD");
	  gSystem->Load("libAOD");
	  gSystem->Load("libANALYSIS");
	  gSystem->Load("libOADB");
	  gSystem->Load("libANALYSISalice");
	  gSystem->Load("libpythia6.so");
	  gSystem->Load("libAliPythia6.so");
      gSystem->Load("libANALYSISalice");   

    // since we will compile a class, tell root where to look for headers  
    //gInterpreter->ProcessLine(".include $ROOTSYS/include");
    //gInterpreter->ProcessLine(".include $ALICE_ROOT/include");


    #ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"
#include "ANALYSIS/macros/AddTaskPIDResponse.C"
#endif   

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#endif

AliPhysicsSelectionTask *PhysSelTask = AddTaskPhysicsSelection(isMC, true);   
    #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);


 //__________________________________________________
 /*  char PhysSelTaskInput[200];
    char MultTaskInput[200];
    char PIDResponseInput[200];
*/
    bool TuneOnData = true;
    char RecoPass[3] = "1";
//it was commented ---->
    //strcpy(PhysSelTaskInput,Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%i,true)",isMC));
    //strcpy(MultTaskInput,"$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    //strcpy(PIDResponseInput,Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%i,true,%i,\"%s\")",isMC,TuneOnData,RecoPass));

  //AliPhysicsSelectionTask *PhysSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(PhysSelTaskInput));
  //AliMultSelectionTask *MultTask = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro(MultTaskInput));
  //AliAnalysisTaskPIDResponse *PIDTask=reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro(PIDResponseInput));

//---instead i add --->
AliPhysicsSelectionTask *PhysSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(
    gInterpreter->ProcessLine(Form(".x $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d,true)", isMC))
);

AliMultSelectionTask *MultTask = reinterpret_cast<AliMultSelectionTask*>(
    gInterpreter->ProcessLine(".x $ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));

AliAnalysisTaskPIDResponse *PIDTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(
    gInterpreter->ProcessLine(Form(".x $ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d,true,%d,\"%s\")", isMC, TuneOnData, RecoPass)));   
//<-----


#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
    AliAnalysisTaskMyTask *task = reinterpret_cast<AliAnalysisTaskMyTask*>(gInterpreter->ExecuteMacro("AddMyTask.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
    gROOT->LoadMacro("AddMyTask.C");
    AliAnalysisTaskMyTask *task = AddMyTask();
#endif
    // compile the class (locally)
//    gROOT->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
    // load the addtask macro
//     gInterpreter->LoadMacro("AddMyTask.C");
    // create an instance of your analysis task
  //  AliAnalysisTaskMyTask *task=reinterpret_cast<AliAnalysisTaskMyTask*>(
   //     gInterpreter->ExecuteMacro("AddMyTask.C()"));

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("AliAOD.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);

    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskMyTask.cxx AliAnalysisTaskMyTask.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskMyTask.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20210808_ROOT6-1"); //vAN-20210808_ROOT6)
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18b"); //alice/data/2011/LHC11h_2

        //commented aber egal ich glaube alienHandler->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root"); 
        //COMMENTED alienHandler->SetGridDataDir("/alice/data/2018/LHC18b");
        alienHandler->SetDataPattern("*AOD/*/AliAOD.root");   
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber
        alienHandler->AddRunNumber(296623);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir("myWorkingDir"); // what i should write here ?
        alienHandler->SetGridOutputDir("myOutputDir"); // and here also

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");  // and here 
            mgr->StartAnalysis("grid");

        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }

gROOT->ProcessLine(".q");