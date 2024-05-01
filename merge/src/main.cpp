#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

#define CHANNEL_NUMBER=160

// Author: Jiaoyang Li (jiaoyang.li@ed.ac.uk)

// reads grid output and creates slimmed file for semi-analytic parameterisation generation
int main( int argc, char * argv[] )
{
  if ( argc != 2 )
  {
    std::cerr << "Please provide list of files to merge (as txt file)" << std::endl;
    return 1;
  }

  //std::string lista_files = "./List_of_Grid_files.txt";
  std::string lista_files = argv[1];
  int numberDevices = CHANNEL_NUMBER; // this is the number of lines in the optical map
  int genPhotons = 100000;

  // Structures to collect data
  struct PointData
  {
    int VUV_hits[CHANNEL_NUMBER];
    int Vis_hits[CHANNEL_NUMBER];
    double Wavelength;
    int GenPhotons;
    int NumberDevices;
  };
  std::map< std::vector<double>, int > originCounter;
  std::map< std::vector<double>, PointData > originData;

  // getting the file names to merge
  std::ifstream Traks_file2( lista_files );
  if( !Traks_file2 ) std::cerr << "WARNING:  Failed to open file with Input file names"<< std::endl;
  Traks_file2.seekg( 0 );
  std::vector< std::string > names;
  std::string nombre;
  while( std::getline( Traks_file2, nombre ) )
  {
    names.push_back(nombre);
  }
  const int n_files = names.size();
  std::cout << "----> number of files: " << n_files << std::endl;

  // create output structures
  const TString root_file_name = "output.root";
  TFile *hfile = new TFile(root_file_name,"RECREATE");
  TTree *myTree = new TTree("myTree","A ROOT tree");
  int VUV_hits[numberDevices];
  int Vis_hits[numberDevices];
  double posX, posY, posZ;
  double wavelength_new;
  myTree->Branch("numberDevices",&numberDevices,"numberDevices/I");
  myTree->Branch("X", &posX, "X/D");
  myTree->Branch("Y", &posY, "Y/D");
  myTree->Branch("Z", &posZ, "Z/D");
  myTree->Branch("VUV_hits",VUV_hits,"VUV_hits[numberDevices]/I");
  myTree->Branch("Vis_hits",Vis_hits,"Vis_hits[numberDevices]/I");
  myTree->Branch("genPhotons", &genPhotons, "genPhotons/I");
  myTree->Branch("Wavelength", &wavelength_new, "Wavelength/D");

  //loop over files
  for(int n=0; n<n_files; n++) {

    // read data file
    std::string input_file = names.at(n);
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"File " << n << ": " << input_file << std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;

    TFile *infile = TFile::Open(input_file.c_str());

    if (!infile || infile->IsZombie()) {
      std::cerr << "Error opening file" << std::endl;
      continue;
    }

    TTree *AllPhotons = nullptr;
    infile->GetObject("pmtresponse/AllPhotons", AllPhotons);

    int EventID, OpChannel; float Wavelength, Time;
    Float_t X, Y, Z;
    if (AllPhotons) {
      AllPhotons->SetBranchAddress("EventID", &EventID);
      AllPhotons->SetBranchAddress("OpChannel", &OpChannel);
      AllPhotons->SetBranchAddress("Wavelength", &Wavelength);
      AllPhotons->SetBranchAddress("Time", &Time);
      AllPhotons->SetBranchAddress("originX", &X);
      AllPhotons->SetBranchAddress("originY", &Y);
      AllPhotons->SetBranchAddress("originZ", &Z);
    }else{
      std::cout<<"corrupted file!"<<std::endl<<std::endl<<std::endl;
      continue;
    }

    // loop through photons
    int num_entries = AllPhotons->GetEntries();
    for (int i = 0; i < num_entries; i++){

      AllPhotons->GetEntry(i);

      // Make a new entry for a new origin
      auto thisOrigin = {X*0.1, Y*0.1, Z*0.1};
      auto searchResult = originData.find( thisOrigin );
      if ( searchResult == originData.end() ) {

        originData[ thisOrigin ] = PointData();

        // Reset counters
        for(int i=0; i<numberDevices; i++) {
          originData[ thisOrigin ].VUV_hits[i] = 0;
          originData[ thisOrigin ].Vis_hits[i] = 0;
        }

        originData[ thisOrigin ].GenPhotons = 100*100000;
        originData[ thisOrigin ].NumberDevices = numberDevices;
        originData[ thisOrigin ].Wavelength = Wavelength; // This doesn't make sense, and I'm pretty sure it's not used
      }
      originCounter[ thisOrigin ]++;

      // apply wavelength cut to separate vuv and visible
      if (Wavelength < 300) {
        originData[ thisOrigin ].VUV_hits[OpChannel]++;
      }
      else {
        originData[ thisOrigin ].Vis_hits[OpChannel]++;
      }
    } // end photon loop

    delete infile;
  } // end file loop

  // Write stored data to tree
  for ( auto const& entry : originData )
  {
    auto point = entry.first;
    auto data = entry.second;

    posX = point[0];
    posY = point[1];
    posZ = point[2];

    for(int i=0; i<numberDevices; i++) {
      VUV_hits[i] = data.VUV_hits[i];
      Vis_hits[i] = data.Vis_hits[i];
    }

    wavelength_new = data.Wavelength;
    genPhotons = data.GenPhotons;

    myTree->Fill();
  }

  // write output file
  hfile->Write();
  hfile->Close();
}
