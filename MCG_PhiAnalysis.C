
// Pythia
#include "Pythia8/Pythia.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TTree.h"

using namespace std;
using namespace ROOT;
using namespace Pythia8;

TBenchmark *fBenchmark      =   new TBenchmark();

const auto fMinIMMC = 0.75;
const auto fMaxIMMC = 1.25;

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
    if (argc < 5 )  {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents] [seed] [option] [Energy]" << endl;
        return -1;
    }
    
    // Definition of number of events
    int   nEvents = atoi(argv[2]);
    
    // Definition of option
    int   fOption = argc < 4 ? 0 : atoi(argv[4]);
    
    // Definition of energy
    int   kEnergy = argc < 5 ? 7000 : atoi(argv[5]);

    // Output File
    TFile * outFile     =   new   TFile   (Form("%s.root",argv[1]),   "recreate", "", 101);
    
    //Output Tree
    Int_t   nPhi;
    Float_t fPx[1024];
    Float_t fPy[1024];
    Float_t fPz[1024];
    TTree*  outTree     =   new TTree(Form("Phi_S%I_E%i_M%i",atoi(argv[3]),kEnergy,fOption),Form("Phi_S%I_E%i_M%i",atoi(argv[3]),kEnergy,fOption));
    outTree->Branch       ("nPhi",  &nPhi,  "nPhi/I");
    outTree->Branch       ("fPx",   &fPx,   "fPx[nPhi]/F");
    outTree->Branch       ("fPy",   &fPy,   "fPy[nPhi]/F");
    outTree->Branch       ("fPz",   &fPz,   "fPz[nPhi]/F");
    
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    // Settings
    pythia.readString("SoftQCD:all = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString(Form("Beams:eCM = %i",kEnergy));
    pythia.readString(Form("333:mMin = %f",fMinIMMC));
    pythia.readString(Form("333:mMax = %f",fMaxIMMC));
    //
    pythia.readString("Random:setSeed = on");
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
        std::vector<Int_t>  kLabel_P3312_XiMinus;
        std::vector<Int_t>  kLabel_M3312_XiPlus;
        std::vector<Int_t>  kLabel_P0333_Phi;
        std::vector<Int_t>  kLabel_P0321_KaonPlus;
        std::vector<Int_t>  kLabel_M0321_KaonMinus;
        //
        // Starting cycling through event particles
        nPhi   =   0;
        auto nMultiplicity  =   0;
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Particles
            if ( Current_Particle.id() == +3312 && (fabs(Current_Particle.p().rap()) < 1.5) )   kLabel_P3312_XiMinus    .push_back(iParticle);
            if ( Current_Particle.id() == -3312 && (fabs(Current_Particle.p().rap()) < 1.5) )   kLabel_M3312_XiPlus     .push_back(iParticle);
            if ( Current_Particle.id() == +333  && (fabs(Current_Particle.p().rap()) < 1.5) )   kLabel_P0333_Phi        .push_back(iParticle);
            if ( Current_Particle.id() == +321  && (fabs(Current_Particle.p().rap()) < 1.5) )   kLabel_P0321_KaonPlus   .push_back(iParticle);
            if ( Current_Particle.id() == -321  && (fabs(Current_Particle.p().rap()) < 1.5) )   kLabel_M0321_KaonMinus  .push_back(iParticle);
            
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
        //  Save for later
        for ( auto kCurrent_Phi : kLabel_P0333_Phi )    {
            // Saving particles
            const auto Current_Particle = pythia.event[kCurrent_Phi];
            //
            fPx[nPhi]   =   Current_Particle.px();
            fPy[nPhi]   =   Current_Particle.py();
            fPz[nPhi]   =   Current_Particle.pz();
            nPhi++;
        }
        outTree->Fill();
    }
    fStopTimer("Production");
    //
    outTree->Write();
    //
    outFile     ->Close();
    return 0;
}
