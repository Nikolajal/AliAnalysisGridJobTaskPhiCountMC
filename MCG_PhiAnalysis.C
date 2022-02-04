
// Pythia
#include "Pythia8/Pythia.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1F.h"

using namespace std;
using namespace ROOT;
using namespace Pythia8;

TBenchmark *fBenchmark      =   new TBenchmark();

//-// InvMass range Pythia MC
const   float   fMinIMMC  =   0.75;
const   float   fMaxIMMC  =   1.25;


//-// pT bins 2D
const   Int_t     nBinPT2D  =   10;
const   Float_t   fMinPT2D  =   0.4;
const   Float_t   fMaxPT2D  =   10.;
        Float_t  *fArrPT2D  =   new Float_t [nBinPT2D+1];
        Float_t  *fArrPT2D_Comp  =   new Float_t [nBinPT2D+2];
void
fSetBinPT2D
()      {
    fArrPT2D[0]     =   0.40; //0.3
    fArrPT2D[1]     =   0.80; //0.2
    fArrPT2D[2]     =   1.00; //0.1
    fArrPT2D[3]     =   1.10; //0.1
    fArrPT2D[4]     =   1.20; //0.2
    fArrPT2D[5]     =   1.40; //0.2
    fArrPT2D[6]     =   1.60; //0.4
    fArrPT2D[7]     =   2.00; //0.8
    fArrPT2D[8]     =   2.80; //1.2
    fArrPT2D[9]     =   4.00; //6.0
    fArrPT2D[10]    =   10.0;
}

TString     fMSG_PrintTimer =   "[INFO] Event # %4.f %s | %02.0f %% | %2.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
//
//_____________________________________________________________________________
//
void    fStartTimer             ( TString fTimerName )    {
    fBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fStopTimer              ( TString fTimerName )     {
    fBenchmark->Stop(fTimerName.Data());
    printf("[INFO] Stopping %s \n", fTimerName.Data());
    Float_t fElapsedS   = (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   = (Int_t)(fElapsedS/60.);
    printf("[INFO] It took %02.0f:%02.0f \n",   fElapsedM,  fElapsedS - 60.*fElapsedM);
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fPrintLoopTimer         ( TString fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
    if ( iEvent%iPrintInterval != 0 || iEvent == 0 ) return;
    
    // Suffix for events
    TString     fSuffix =   "";
    Int_t       fSfxCor =   iPrintInterval;
    if ( iPrintInterval/1e3 < 1e3         && iPrintInterval/1e3 >= 1 )       {
        fSuffix =   "k";
        fSfxCor =   (int)(iPrintInterval/1e3) + iPrintInterval%(int)1e3;
    }
    if ( iPrintInterval/1e6 < 1e6         && iPrintInterval/1e6 >= 1 )       {
        fSuffix =   "mln";
        fSfxCor =   (int)(iPrintInterval/1e6) + iPrintInterval%(int)1e6;
        cout << fSfxCor << endl;
    }
    if ( iPrintInterval/1e9 < 1e9         && iPrintInterval/1e9 >= 1 )       {
        fSuffix =   "mld";
        fSfxCor =   (int)(iPrintInterval/1e9) + iPrintInterval%(int)1e9;
        cout << fSfxCor << endl;
    }
    
    // Stopping timer
    fBenchmark->Stop(fTimerName.Data());
    
    // Evaluating informations
    Float_t fFraction   =   (float)iEvent/((float)nEntries);
    Float_t fElapsedS   =   (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   =   (Int_t)(fElapsedS/60.);
    Float_t fPrintEvt   =   (float)iEvent*(float)fSfxCor/((float)iPrintInterval);
    Float_t fSpeedvsS   =   fPrintEvt/fElapsedS;
    Float_t fEta____S   =   (Int_t)(fElapsedS*((float)nEntries/((float)iEvent) -1));
    Float_t fEta____M   =   (Int_t)(fEta____S/60.);
    
    // Printing
    //"[INFO] Event # %4.f %s | %02.0f %% | %1.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
    printf(fMSG_PrintTimer.Data(),  fPrintEvt,  fSuffix.Data(), 100.*fFraction, fSpeedvsS,  fSuffix.Data(), fElapsedM,  fElapsedS -60.*fElapsedM,  fEta____M,  fEta____S -60.*fEta____M);
    fflush(stdout);
    
    // Resuming timer
    fBenchmark->Start(fTimerName.Data());
}
//
//_____________________________________________________________________________
//
int
main
 ( int argc, char *argv[] ) {
    // Check everything is good
    if (argc < 4 )  {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents] [seed] [option]" << endl;
        return -1;
    }
    
    // Setting bins
    fSetBinPT2D();
    
    // Definition of number of events
    int   nEvents = atoi(argv[2]);
    
    // Definition of option
    int   fOption = argc < 4 ? 0 : atoi(argv[4]);
    
    // Input File
    TFile * inHIFile    = new   TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/DATA/COMMON_sysErrRelInleasticXSec.root");
    TFile * inMBFile    = new   TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/DATA/COMMON_sysErrRelMaterial.root");
    
    // Input File
    TH1F*   inHIRefKplus        =   (TH1F*)(inHIFile->Get("hKaonPlusMinBiasSysErrRel"));
    TH1F*   inHIRefKminus       =   (TH1F*)(inHIFile->Get("hKaonMinusMinBiasSysErrRel"));
    TH1F*   inMBRefKplus        =   (TH1F*)(inMBFile->Get("hKaonMinusMinBiasSysErrRel"));
    TH1F*   inMBRefKminus       =   (TH1F*)(inMBFile->Get("hKaonMinusMinBiasSysErrRel"));
    
    // Output File
    TFile * outFile     = new   TFile   (Form("%s.root",argv[1]),   "recreate", "", 101);
    
    // Output - Standard
    TProfile*   hMeanPT     =   new TProfile("hMeanPT","hMeanPT",nBinPT2D,fArrPT2D);
    TProfile*   hHISysErr   =   new TProfile("hHISysErr","hHISysErr",100,0,20);
    TProfile*   hMBSysErr   =   new TProfile("hMBSysErr","hMBSysErr",100,0,20);
    TProfile2D* hHISysEr2   =   new TProfile2D("hHISysEr2","hHISysEr2",100,0,10,100,0,20);
    TProfile2D* hMBSysEr2   =   new TProfile2D("hMBSysEr2","hMBSysEr2",100,0,10,100,0,20);
    TProfile*   hYPhiProd   =   new TProfile("hYPhiProd","hYPhiProd",10,-0.5,9.5);
    TProfile*   hYPhiPro1   =   new TProfile("hYPhiPro1","hYPhiPro1",100,0.,100.);
    TProfile*   hYPhiPro2   =   new TProfile("hYPhiPro2","hYPhiPro2",100,0.,100.);
    TH1F*       hPhiProd    =   new TH1F("hPhiProd","hPhiProd",10,-0.5,9.5);
    //
    //  --- Output Correlation search
    TH2F*       hRapCorrXiPM            =   new TH2F("hRapCorrXiPM",        "hRapCorrXiPM",         500,-10,10,     500,-10,10);
    TH2F*       hRapCorrXiPP            =   new TH2F("hRapCorrXiPP",        "hRapCorrXiPP",         500,-10,10,     500,-10,10);
    TH2F*       hRapCorrXiMM            =   new TH2F("hRapCorrXiMM",        "hRapCorrXiMM",         500,-10,10,     500,-10,10);
    TH1F*       hRapCorrXiPM1D          =   new TH1F("hRapCorrXiPM1D",      "hRapCorrXiPM1D",       500,-10,10);
    TH1F*       hPhiCorrXiPM1D          =   new TH1F("hPhiCorrXiPM1D",      "hPhiCorrXiPM1D",       360,-180,180);
    //
    TH2F*       hRapCorrPhi             =   new TH2F("hRapCorrPhi",         "hRapCorrPhi",          500,-10,10,     500,-10,10);
    TH1F*       hRapCorrPhi1D           =   new TH1F("hRapCorrPhi1D",       "hRapCorrPhi1D",        500,-10,10);
    TH1F*       hPhiCorrPhi1D           =   new TH1F("hPhiCorrPhi1D",       "hPhiCorrPhi1D",        360,-180,180);
    TH1F*       hPhiCorrPhi1D_BKG       =   new TH1F("hPhiCorrPhi1D_BKG",   "hPhiCorrPhi1D_BKG",    360,-180,180);
    //
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    // Settings
    pythia.readString("SoftQCD:all = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5000");
    pythia.readString(Form("333:mMin = %f",fMinIMMC));
    pythia.readString(Form("333:mMax = %f",fMaxIMMC));
    //
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %i",atoi(argv[3])));
    pythia.readString("Print:quiet = on");
    pythia.readString(Form("Random:seed = %i",atoi(argv[3])));
    switch (fOption)
    {
        case 1:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("BeamRemnants:remnantMode = 1");
            cout << "[INFO] Using Colour reconnection Mode 1" << endl;
            cout << "The new more QCD based scheme." << endl;
            break;
        case 2:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 2");
            cout << "[INFO] Using Colour reconnection Mode 2" << endl;
            cout << "The new gluon-move model." << endl;
            break;
        case 3:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 3");
            pythia.readString("ColourReconnection:forceResonance = on");
            pythia.readString("PartonLevel:earlyResDec = off");
            cout << "[INFO] Using Colour reconnection Mode 3" << endl;
            cout << "The SK I e^+ e^- CR model." << endl;
            break;
        case 4:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 4");
            pythia.readString("ColourReconnection:forceResonance = on");
            pythia.readString("PartonLevel:earlyResDec = off");
            cout << "[INFO] Using Colour reconnection Mode 4" << endl;
            cout << "The SK II e^+ e^- CR model." << endl;
            break;
        case 5:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 2.9");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.43");
            pythia.readString("ColourReconnection:timeDilationMode = 0");
            cout << "[INFO] Using Colour reconnection Mode 5" << endl;
            cout << "Mode 0 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 6:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.20");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            cout << "[INFO] Using Colour reconnection Mode 6" << endl;
            cout << "Mode 2 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 7:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.15");
            pythia.readString("ColourReconnection:timeDilationMode = 3");
            pythia.readString("ColourReconnection:timeDilationPar = 0.073");
            cout << "[INFO] Using Colour reconnection Mode 7" << endl;
            cout << "Mode 3 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 8:
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.2");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            pythia.readString("Ropewalk:RopeHadronization = on");
            pythia.readString("Ropewalk:doShoving = on");
            pythia.readString("Ropewalk:tInit = 1.5");
            pythia.readString("Ropewalk:deltat = 0.05");
            pythia.readString("Ropewalk:tShove = 0.1");
            pythia.readString("Ropewalk:gAmplitude = 0");
            pythia.readString("Ropewalk:doFlavour = on");
            pythia.readString("Ropewalk:r0 = 0.5");
            pythia.readString("Ropewalk:m0 = 0.2");
            pythia.readString("Ropewalk:beta = 0.1");
            pythia.readString("PartonVertex:setVertex = on");
            pythia.readString("PartonVertex:protonRadius = 0.7");
            pythia.readString("PartonVertex:emissionWidth = 0.1");
            cout << "[INFO] Using Colour reconnection Mode 8" << endl;
            cout << "Custom mode 1: QCD based CR w/ ropes" << endl;
            break;
        case 9:
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:junctionCorrection = 1.2");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            pythia.readString("Ropewalk:RopeHadronization = off");
            cout << "[INFO] Using Colour reconnection Mode 9" << endl;
            cout << "Custom mode 2: QCD based CR w/o ropes" << endl;
            break;
        default:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 0");
            cout << "[INFO] Using default settings" << endl;
            cout << "To change settings choose from below:" << endl;
            cout << "D. The MPI-based original Pythia 8 scheme." << endl;
            cout << "1. The new more QCD based scheme." << endl;
            cout << "2. The new gluon-move model." << endl;
            cout << "3. The SK I e^+ e^- CR model." << endl;
            cout << "4. The SK II e^+ e^- CR model." << endl;
            cout << "5. Mode 0 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "6. Mode 2 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "7. Mode 3 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "8. Custom mode 1: QCD based CR w/ ropes" << endl;
            cout << "9. Custom mode 2: QCD based CR w/o ropes" << endl;
            break;
    }
    pythia.init();
    
    // Start Timer
    fStartTimer("Production");
    
    // Cycling through events
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )  {
        // Next event
        pythia.next();
        fPrintLoopTimer("Production",iEvent,nEvents,10000);
        //
        std::vector<Float_t> vMeanPT;
        //
        std::vector<Int_t>  kLabel_P3312_XiMinus;
        std::vector<Int_t>  kLabel_M3312_XiPlus;
        std::vector<Int_t>  kLabel_P0333_Phi;
        std::vector<Int_t>  kLabel_P0321_KaonPlus;
        std::vector<Int_t>  kLabel_M0321_KaonMinus;
        //
        // Starting cycling through event particles
        auto nPhi   =   0;
        auto nMultiplicity  =   0;
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( Current_Particle.id() == +3312 && (fabs(Current_Particle.p().rap()) < 0.5) )   kLabel_P3312_XiMinus    .push_back(iParticle);
            if ( Current_Particle.id() == -3312 && (fabs(Current_Particle.p().rap()) < 0.5) )   kLabel_M3312_XiPlus     .push_back(iParticle);
            if ( Current_Particle.id() == +333  && (fabs(Current_Particle.p().rap()) < 0.5) )   kLabel_P0333_Phi        .push_back(iParticle);
            if ( Current_Particle.id() == +321  && (fabs(Current_Particle.p().rap()) < 0.5) )   kLabel_P0321_KaonPlus   .push_back(iParticle);
            if ( Current_Particle.id() == -321  && (fabs(Current_Particle.p().rap()) < 0.5) )   kLabel_M0321_KaonMinus  .push_back(iParticle);
            
            //  Calculating Multiplicity
            
            //Skipping non-final particles
            if ( !Current_Particle.isFinal() )      continue;
            
            //Skipping non-charged particles
            if ( !Current_Particle.isCharged() )    continue;
            
            //Skipping particles outside eta range
            if ( !(fabs(Current_Particle.eta()) < 0.5) )  continue;
            
            // Mutliplicity evaluation
            nMultiplicity++;
            
        }
        //
        //  Correlation Search
        //
        for ( auto iXiPlus : kLabel_M3312_XiPlus )   {
            const auto  iParticle_XiPlus    =   pythia.event[iXiPlus];
            const auto  iRapidity_XiPlus    =   iParticle_XiPlus.p().rap();
            for ( auto jXiPlus : kLabel_M3312_XiPlus )   {
                if ( iXiPlus == jXiPlus ) continue;
                const auto  jParticle_XiPlus    =   pythia.event[jXiPlus];
                const auto  jRapidity_XiPlus    =   jParticle_XiPlus.p().rap();
                hRapCorrXiPP    ->  Fill    ( iRapidity_XiPlus, jRapidity_XiPlus );
            }
            for ( auto jXiMinus : kLabel_P3312_XiMinus )   {
                const auto  jParticle_XiMinus    =   pythia.event[jXiMinus];
                const auto  jRapidity_XiMinus    =   jParticle_XiMinus.p().rap();
                hRapCorrXiPM        ->  Fill    ( iRapidity_XiPlus, jRapidity_XiMinus );
                hRapCorrXiPM1D      ->  Fill    ( iRapidity_XiPlus - jRapidity_XiMinus );
                auto    kDeltaPhi   =   ( iParticle_XiPlus.phi() - jParticle_XiMinus.phi()  )*360/( TMath::Pi()*2 );
                kDeltaPhi = kDeltaPhi < -180 ? kDeltaPhi + 360 : kDeltaPhi > 180 ? kDeltaPhi -360 : kDeltaPhi;
                hPhiCorrXiPM1D      ->  Fill    ( kDeltaPhi );
            }
        }
        for ( auto iPhi : kLabel_P0333_Phi )   {
            const auto  iParticle_Phi       =   pythia.event[iPhi];
            const auto  iRapidity_Phi       =   iParticle_Phi.p().rap();
            for ( auto jPhi : kLabel_P0333_Phi )   {
                if ( iPhi >= jPhi ) continue;
                const auto  jParticle_Phi       =   pythia.event[jPhi];
                const auto  jRapidity_Phi       =   jParticle_Phi.p().rap();
                hRapCorrPhi         ->  Fill    ( iRapidity_Phi, jRapidity_Phi );
                hRapCorrPhi1D       ->  Fill    ( iRapidity_Phi - jRapidity_Phi );
                auto    kDeltaPhi   =   ( iParticle_Phi.phi() - jParticle_Phi.phi()  )*360/( TMath::Pi()*2 );
                kDeltaPhi = kDeltaPhi < -180 ? kDeltaPhi + 360 : kDeltaPhi > 180 ? kDeltaPhi -360 : kDeltaPhi;
                hPhiCorrPhi1D       ->  Fill    ( kDeltaPhi );
            }
        }
        //
        //  Standard
        //
        for ( auto kiPhi : kLabel_P0333_Phi ) {
            // Saving particles
            const auto  Current_Particle    =   pythia.event[kiPhi];
            const auto  uUtilityPhiRap      =   Current_Particle.p().rap();
            //
            nPhi++;
            vMeanPT.push_back(Current_Particle.pT());
            if ( Current_Particle.daughterList().size() != 2 ) continue;
            auto    nKaon1  =   Current_Particle.daughter1();
            auto    fKaon1  =   pythia.event[nKaon1];
            if ( fabs( fKaon1.id() ) != 321 ) continue;
            auto    nKaon2  =   Current_Particle.daughter2();
            auto    fKaon2  =   pythia.event[nKaon2];
            if ( fabs( fKaon2.id() ) != 321 ) continue;
            for ( auto iKmnus = kLabel_M0321_KaonMinus.begin(); iKmnus != kLabel_M0321_KaonMinus.end(); ++iKmnus ) { if ( (*iKmnus == Current_Particle.daughter1()) || (*iKmnus == Current_Particle.daughter2()) ){ kLabel_M0321_KaonMinus.erase(iKmnus); break; }  }
            for ( auto iKplus = kLabel_P0321_KaonPlus.begin(); iKplus != kLabel_P0321_KaonPlus.end(); ++iKplus ) { if ( (*iKplus == Current_Particle.daughter1()) || (*iKplus == Current_Particle.daughter2()) ){ kLabel_P0321_KaonPlus.erase(iKplus); break; } }
            auto    fUncerHI1   =   fKaon1.charge() > 0 ? inHIRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKaon1.pT(),1.299) ) ) : inHIRefKminus->GetBinContent( inHIRefKminus->FindBin( fKaon1.pT() ) ) ;
            auto    fUncerMB1   =   fKaon1.charge() > 0 ? inMBRefKplus->GetBinContent( inMBRefKplus->FindBin( min(fKaon1.pT(),1.299) ) ) : inMBRefKminus->GetBinContent( inMBRefKminus->FindBin( fKaon1.pT() ) ) ;
            auto    fUncerHI2   =   fKaon2.charge() > 0 ? inHIRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKaon2.pT(),1.299) ) ) : inHIRefKminus->GetBinContent( inHIRefKminus->FindBin( fKaon2.pT() ) ) ;
            auto    fUncerMB2   =   fKaon2.charge() > 0 ? inHIRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKaon2.pT(),1.299) ) ) : inMBRefKminus->GetBinContent( inHIRefKminus->FindBin( fKaon2.pT() ) ) ;
            hMBSysErr->Fill(Current_Particle.pT(), fUncerMB1 + fUncerMB2 );
            hHISysErr->Fill(Current_Particle.pT(), fUncerHI1 + fUncerHI2 );
            for ( auto kjPhi : kLabel_P0333_Phi ) {
                // Saving particles
                if ( kiPhi == kjPhi ) continue;
                const auto Current_Particl2 = pythia.event[kjPhi];
                if ( Current_Particl2.daughterList().size() != 2 ) continue;
                auto    nKao21  =   Current_Particl2.daughter1();
                auto    fKao21  =   pythia.event[nKao21];
                if ( fabs( fKao21.id() ) != 321 ) continue;
                auto    nKao22  =   Current_Particl2.daughter2();
                auto    fKao22  =   pythia.event[nKao22];
                if ( fabs( fKao22.id() ) != 321 ) continue;
                auto    fUnce2HI1   =   fKao21.charge() > 0 ? inHIRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKao21.pT(),1.299) ) ) : inHIRefKminus->GetBinContent( inHIRefKminus->FindBin( fKao21.pT() ) ) ;
                auto    fUnce2MB1   =   fKao21.charge() > 0 ? inMBRefKplus->GetBinContent( inMBRefKplus->FindBin( min(fKao21.pT(),1.299) ) ) : inMBRefKminus->GetBinContent( inMBRefKminus->FindBin( fKao21.pT() ) ) ;
                auto    fUnce2HI2   =   fKao22.charge() > 0 ? inHIRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKao22.pT(),1.299) ) ) : inHIRefKminus->GetBinContent( inHIRefKminus->FindBin( fKao22.pT() ) ) ;
                auto    fUnce2MB2   =   fKao22.charge() > 0 ? inMBRefKplus->GetBinContent( inHIRefKplus->FindBin( min(fKao22.pT(),1.299) ) ) : inMBRefKminus->GetBinContent( inHIRefKminus->FindBin( fKao22.pT() ) ) ;
                hMBSysEr2->Fill(Current_Particle.pT(), Current_Particl2.pT(), fUncerMB1 + fUncerMB2 + fUnce2MB1 + fUnce2MB2, 0.5 );
                hHISysEr2->Fill(Current_Particle.pT(), Current_Particl2.pT(), fUncerHI1 + fUncerHI2 + fUnce2HI1 + fUnce2HI2, 0.5 );
            }
        }
        hYPhiProd->Fill(1,nPhi);
        hYPhiProd->Fill(2,nPhi*(nPhi-1)*0.5);
        hYPhiPro1->Fill(nMultiplicity,nPhi);
        hYPhiPro2->Fill(nMultiplicity,nPhi*(nPhi-1)*0.5);
        hPhiProd->Fill(nPhi);
        if ( nPhi < 2 ) continue;
        for ( auto iPT : vMeanPT )   {
            for ( auto jPT : vMeanPT )   {
                if ( iPT == jPT ) continue;
                hMeanPT->Fill(iPT,jPT);
            }
        }
    }
    
    fStopTimer("Production");
    //
    hRapCorrXiPP->Write();
    hRapCorrXiMM->Write();
    hRapCorrXiPM->Write();
    hRapCorrXiPM1D->Write();
    hPhiCorrXiPM1D->Write();
    hRapCorrPhi ->Write();
    hRapCorrPhi1D->Write();
    hPhiCorrPhi1D->Write();
    hPhiCorrPhi1D_BKG->Write();
    //
    hPhiProd->Scale(1./nEvents);
    hHISysErr->Write();
    hMBSysErr->Write();
    hHISysEr2->Write();
    hMBSysEr2->Write();
    hYPhiProd->Write();
    hYPhiPro1->Write();
    hYPhiPro2->Write();
    hMeanPT->Write();
    hPhiProd->Write();
    //
    
    //
    outFile     ->Close();
    return 0;
}
