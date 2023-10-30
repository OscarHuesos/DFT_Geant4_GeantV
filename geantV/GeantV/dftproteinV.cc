

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

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
#include "CommandLineParser.hh"

using namespace G4DNAPARSER;
CommandLineParser* parser(0);

void Parse(int& argc, char** argv);

int main(int argc,char** argv){


Parse(argc, argv);
printf("se inicia \n");

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

  if ((commandLine = parser->GetCommandIfActive("-gui")))
  {
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

  delete visManager;
  delete runManager;

  return 0;
}



void Parse(int& argc, char** argv) {

  parser = CommandLineParser::GetParser();

  parser->AddCommand("-gui",
                     Command::OptionNotCompulsory,
                    "Select geant4 UI or just launch a geant4 terminal session",
                    "qt");

  parser->AddCommand("-mac",
                     Command::WithOption,
                     "Give a mac file to execute",
                     "pdb4dna.in");


#ifdef G4MULTITHREADED
  parser->AddCommand("-mt", Command::WithOption,
                     "Launch in MT mode (events computed in parallel)",
                     "2");
#endif

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


  if (parser->CheckIfNotHandledOptionsExists(argc, argv)) {
    abort();
  }
}
