
#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "CommandLineParser.hh"

#include "Geant/RunManager.h"
#include "Geant/PhysicsProcessHandler.h"
#include "Geant/PhysicsListManager.h"
#include "Geant/UserFieldConstruction.h"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#endif

// FULL-CMS application
#include "CMSFullApp.h"
#include "CMSDetectorConstruction.h"
#include "CMSFieldConstruction.h"
#include "CMSParticleGun.h"
#include "CMSPhysicsList.h"


using namespace G4DNAPARSER;
CommandLineParser* parser(0);

void Parse(int& argc, char** argv);


int main(int argc,char** argv){
printf("inicia \n");
Parse(argc, argv);

CLHEP::RanecuEngine defaultEngine(1234567);
G4Random::setTheEngine(&defaultEngine);
CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

Command* commandLine(0);

  #ifdef G4MULTITHREADED
  G4MTRunManager* runManager= new G4MTRunManager;
  if ((commandLine = parser->GetCommandIfActive("-mt")))
  {
    int nThreads = 2;
    if(commandLine->GetOption() == "NMAX")
    {
     nThreads = G4Threading::G4GetNumberOfCores();
    }
    else
    {
     nThreads = G4UIcommand::ConvertToInt(commandLine->GetOption());
    }
    G4cout << "===== PDB4DNA is started with "
       << runManager->GetNumberOfThreads()
       << " threads =====" << G4endl;

    runManager->SetNumberOfThreads(nThreads);
  }
#else
  G4RunManager* runManager = new G4RunManager();
#endif

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive* ui(0);


  if ((commandLine = parser->GetCommandIfActive("-gui")))  {
    ui = new G4UIExecutive(argc, argv,
                           commandLine->GetOption());

    if(parser->GetCommandIfActive("-novis") == 0)
    {
      if ((commandLine = parser->GetCommandIfActive("-vis")))

      {
        UImanager->ApplyCommand(G4String("/vis/open ")+
                                commandLine->GetOption());
      }
      else

      {
        UImanager->ApplyCommand("/vis/open OGL 600x600-0+0");
      }
      UImanager->ApplyCommand("/control/execute vis.mac");
    }
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
  }
  else

  {
    if ((commandLine = parser->GetCommandIfActive("-vis")))
    {
      UImanager->ApplyCommand(G4String("/vis/open ")+
                              commandLine->GetOption());
      UImanager->ApplyCommand("/control/execute vis.mac");
    }
  }

  if ((commandLine = parser->GetCommandIfActive("-mac")))
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + commandLine->GetOption());
  }
  else
  {
    UImanager->ApplyCommand("/control/execute init.mac");
  }

  if ((commandLine = parser->GetCommandIfActive("-gui")))
  {
#ifdef G4UI_USE_QT
    G4UIQt* UIQt = static_cast<G4UIQt*> (UImanager->GetG4UIWindow());
    if ( UIQt) {
      UIQt->AddViewerTabFromFile("README", "README from "+ G4String(argv[0]));
    }
#endif
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  delete visManager;
  delete runManager;
return 0;
}


std::string parDetGDMFile = ""; 

std::string parGunPrimaryParticleName = "";           
int parGunPrimaryPerEvent             = 0;        
double parGunPrimaryKinEnergy         = 0.;     
double parGunPrimaryDir[3]            = {0., 0., 0.}; 

int parConfigNumBufferedEvt = 4;    
int parConfigNumRunEvt      = 10;   
int parConfigNumThreads     = 4;     
int parConfigNumPropagators = 1;   
bool parConfigIsPerformance = false; 


static struct option options[] = {{"gun-set-primary-energy", required_argument, 0, 'a'},
                                  {"gun-set-primary-type", required_argument, 0, 'b'},
                                  {"gun-set-primary-per-event", required_argument, 0, 'c'},
                                  {"gun-set-primary-direction", required_argument, 0, 'd'},

                                  {"det-set-gdml", required_argument, 0, 'e'},

                                  {"config-number-of-buffered-events", required_argument, 0, 'm'},
                                  {"config-total-number-of-events", required_argument, 0, 'n'},
                                  {"config-number-of-threads", required_argument, 0, 'p'},
                                  {"config-number-of-propagators", required_argument, 0, 'q'},
                                  {"config-run-performance", no_argument, 0, 'r'},

                                  {"help", no_argument, 0, 'h'},
                                  {0, 0, 0, 0}};

enum PRIMDIR_OPTIONS { PRIMDIR_X_OPT = 0, PRIMDIR_Y_OPT, PRIMDIR_Z_OPT };
char *const primdir_token[] = {[PRIMDIR_OPTIONS::PRIMDIR_X_OPT] = (char *const) "x",
                               [PRIMDIR_OPTIONS::PRIMDIR_Y_OPT] = (char *const) "y",
                               [PRIMDIR_OPTIONS::PRIMDIR_Z_OPT] = (char *const) "z",
                               NULL};


void help(){
  printf("\nUsage: fullLHCbApp [OPTIONS] INPUT_FILE\n\n");
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}


void GetArguments(int argc, char *argv[]) {

  int errfnd = 0;
  char *subopts;
  char *value;

  while (true) {
    int c, optidx = 0;
    //
    c = getopt_long(argc, argv, "", options, &optidx);
    //
    if (c == -1) break;
    //
    switch (c) {
    case 0:
      c = options[optidx].val; 
      
    case 'a':
      parGunPrimaryKinEnergy = strtod(optarg, NULL);
      if (parGunPrimaryKinEnergy <= 0) errx(1, "primary particle energy must be positive");
      break;
    case 'b':
      parGunPrimaryParticleName = optarg;
      break;
    case 'c':
      parGunPrimaryPerEvent = (int)strtol(optarg, NULL, 10);
      break;
    case 'd': // primary direction sub-optarg
      subopts = optarg;
      while (*subopts != '\0' && !errfnd) {
        switch (getsubopt(&subopts, primdir_token, &value)) {
        case PRIMDIR_OPTIONS::PRIMDIR_X_OPT:
          parGunPrimaryDir[0] = strtod(value, NULL);
          break;
        case PRIMDIR_OPTIONS::PRIMDIR_Y_OPT:
          parGunPrimaryDir[1] = strtod(value, NULL);
          break;
        case PRIMDIR_OPTIONS::PRIMDIR_Z_OPT:
          parGunPrimaryDir[2] = strtod(value, NULL);
          break;
        default:
          fprintf(stderr, "No match found for token: [%s] among PRIMDIR_OPTIONS", value);
          errfnd = 1;
          exit(0);
          break;
        }
      }

      break;
    //---- Detector
    case 'e':
      parDetGDMFile = optarg;
      break;
    //---- Run configuration
    case 'm':
      parConfigNumBufferedEvt = (int)strtol(optarg, NULL, 10);
      break;
    case 'n':
      parConfigNumRunEvt = (int)strtol(optarg, NULL, 10);
      break;
    case 'p':
      parConfigNumThreads = (int)strtol(optarg, NULL, 10);
      break;
    case 'q':
      parConfigNumPropagators = (int)strtol(optarg, NULL, 10);
      break;
    case 'r':
      parConfigIsPerformance = true;
      break;
    //---- Help
    case 'h':
      help();
      exit(0);
      break;
    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }
}

geant::RunManager *RunManager() {

  geant::GeantConfig *runConfig = new geant::GeantConfig();
  geant::RunManager *runManager = new geant::RunManager(parConfigNumPropagators, parConfigNumThreads, runConfig);

  runConfig->fNtotal = parConfigNumRunEvt;
  runConfig->fNbuff  = parConfigNumBufferedEvt;
  runConfig->fUseV3         = true;
  runConfig->fNminThreshold = 5 * parConfigNumThreads;
  // Set threshold for tracks to be reused in the same volume
  runConfig->fNminReuse = 100000;
  //
  // Activate standard scoring
  runConfig->fUseStdScoring = true;
  runManager->SetPhysicsInterface(new geantphysics::PhysicsProcessHandler(*runConfig));

  return runManager;
}


void Parse(int& argc, char** argv){

 parser = CommandLineParser::GetParser();

  parser->AddCommand("-gui",
                     Command::OptionNotCompulsory,
                    "Select geant4 UI or just launch a geant4 terminal session",
                    "qt");

  parser->AddCommand("-mac",
                     Command::WithOption,
                     "Give a mac file to execute",
                     "pdb4dna.in");

  parser->AddCommand("-mt", Command::WithOption,
                     "Launch in MT mode (events computed in parallel)",
                     "2");

  parser->AddCommand("-vis",
                     Command::WithOption,
                     "Select a visualization driver",
                     "OGL 600x600-0+0");

  parser->AddCommand("-novis",
                     Command::WithoutOption,
                     "Deactivate visualization when using GUI");


  if (parser->Parse(argc, argv) != 0) {

    CommandLineParser::DeleteInstance();
    std::exit(0);
  }


  if (parser->CheckIfNotHandledOptionsExists(argc, argv)){

   abort();
  }
}




