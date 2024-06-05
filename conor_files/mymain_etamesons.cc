// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.


//CONOR: copied form mymain91
//Now rewritten to use MyGenParticle as simple struct, not custom class
// CONOR: now does specifically eta/eta' meson studies
// pass it: etameson.card


// Stdlib header file for input and output.
#include <iostream>
// for sort
#include <algorithm>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming and trees
#include "TH1.h"
#include "TTree.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

// my own GenParticle class
//#include "MyGenParticle.h"
//CONOR: struct version
#include "MyGenParticleStruct.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  //Conor: seemed like the line Tapp below was eating up the argv args!
  // and it looks like things run okay wihout it ...
  //TApplication theApp("hist", &argc, argv);

  // Create Pythia instance
  
  Pythia pythia;

  // read setup from external card file passed as arg
  cout << "argc = " <<argc <<endl;
  cout << argv[0] <<endl;
  cout << argv[1] <<endl;

  
  if (argc != 2) {
    cout << "Setup Pythia from external card file, to be passed as arg" <<endl;
    return -1;
  }
  //    pythia.readFile("main04.cmnd");
  else {
    // Check that the provided input name corresponds to an existing file.
    ifstream is(argv[1]);
    if (!is) {
      cerr << " Command-line file " << argv[1] << " was not found. \n"
           << " Program stopped! " << endl;
      return 1;
    }
    pythia.readFile(argv[1]);
  }

  
  // original main91:
  //and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 20.");
  //pythia.readString("Beams:eCM = 14000.");

  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");

  pythia.init();

  // print some info:
  // see: https://pythia.org/latest-manual/ParticleDataScheme.html
  pythia.particleData.listChanged();
  const int eta_code = 221;
  const int etap_code = 331;
  pythia.particleData.list(eta_code);
  pythia.particleData.list(etap_code);


  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("hist.root", "RECREATE");

  // Book histogram.
  TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
  TH1F *nPhotons_h = new TH1F("nPhotons_h","N Photons",10,0,10); 

  // Eta hists - count all etas per event
  // and count N in LHCB acc too
  TH1F *nEta_h = new TH1F("nEta_h","",100,0,50);
  TH1F *nEta_lhcbacc_h = new TH1F("nEta_lhcbacc_h","",100,0,50);

  // pT, ptot of etas in LHCB acc
  TH1F *etameson_pt_lhcbacc_h = new TH1F("etameson_pt_lhcbacc_h","",100,0,100);
  TH1F *etameson_ptot_lhcbacc_h = new TH1F("etameson_ptot_lhcbacc_h","",100,0,100);

  // spread of decay prods when eta meson is in LHCb acc
  TH1F *etameson_decay_spread_h = new TH1F("etameson_decay_spread_h","",100,0,10);

  // pT of all muon daughters from 4mu decay, when all 4 in LHCb acc
  TH1F *muon_pt_all4lhcbacc_h = new TH1F("muon_pt_all4lhcbacc_h","",100,0,10);
  
  //Tree and assoc objects/structs to store info
  TTree *outTree = new TTree("outTree","Pythia Tree");
  CH::MyGenParticleStruct quark_in;
  CH::MyGenParticleStruct qbar_in;
  CH::MyGenParticleStruct W;
  CH::MyGenParticleStruct lepton;
  CH::MyGenParticleStruct neutrino;
  


  CH::InitMyGenParticleStruct(quark_in);
  CH::InitMyGenParticleStruct(qbar_in);
  CH::InitMyGenParticleStruct(W);
  CH::InitMyGenParticleStruct(lepton);
  CH::InitMyGenParticleStruct(neutrino);
  
  
  //  CH::MyGenParticleStruct selObj2;
  // I can adjust the split level to have different ways of storing the class objects in the branch
  // default level is 1
  // let me also try 0 ? This stores is *Object*
  /// Note: split level is 4ht arg, so put buf size as 3rd
  //  outTree->Branch("Obj1",&selObj1,32000,0);
  //  outTree->Branch("Obj1",&selObj1);
  // outTree->Branch("Obj2",&selObj2);

  // struct version of Tree:
  outTree->Branch("Q",&quark_in,CH::MyGenParticleStructBranchDefString.c_str());
  outTree->Branch("Qbar",&qbar_in,CH::MyGenParticleStructBranchDefString.c_str());
  outTree->Branch("W",&W,CH::MyGenParticleStructBranchDefString.c_str());
  outTree->Branch("Lepton",&lepton,CH::MyGenParticleStructBranchDefString.c_str());
  outTree->Branch("Neutrino",&neutrino,CH::MyGenParticleStructBranchDefString.c_str());
  
  //  outTree->Branch("Obj2",&selObj2,CH::MyGenParticleStructBranchDefString.c_str());
  
  outTree->Print();
  
  /// COunters over all events
  int nEtasWithAllFourDaughtersInAcc = 0;
  int nEtasWithAllFourDaughtersInAccAndTwoMuons500 = 0;
  int nEtasWithAllFourDaughtersInAccAndFourMuons250 = 0;
  int nEtasWithAllFourDaughtersInAccAndEitherTrig = 0;
  int nEtasWithAllFourDaughtersInAccAndBothTrigs = 0;

  // FLAG to set eta vs eta' meson mode
  const int eta_meson_pid = 221;
  const int etap_meson_pid = 331;
  //const int meson_selected_pid = eta_meson_pid;
    const int meson_selected_pid = etap_meson_pid;
  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    
    //now loop over particles in this event
    // event is basically a vector, and event[i] is Particle
    
    // Find number of all final charged particles.
    int nCharged = 0;
    int nSelectedObj = 0;

    //use the Pythia own particle class, and only copy to my struct for Tree at very end
    std::vector<Particle> selectedObjects;


    // we will count all objects with some min pt cut as we loop through,
    //then sort at end to get top two


    CH::InitMyGenParticleStruct(quark_in);
    CH::InitMyGenParticleStruct(qbar_in);
    CH::InitMyGenParticleStruct(W);
    CH::InitMyGenParticleStruct(lepton);
    CH::InitMyGenParticleStruct(neutrino);
    
    // now step through all particles in this event

    int nEtasThisEvent = 0;
    int nEtasInLHCbAccThisEvent = 0;
    
    for (int i = 0; i < pythia.event.size(); ++i) {

      Particle &thisParticle = pythia.event[i];
      // particle.id() is PDG code
      // isFinal()

      // select e/mu if FInal and some min pt (20 geV)
      // if(pythia.event[i].isFinal() && pythia.event[i].pT()>0 && (fabs(pythia.event[i].id())==24 || fabs(pythia.event[i].id())==13)  ) {

      //      if()


      if(pythia.event[i].id()==meson_selected_pid) {
	// then its an eta meson (eta', depending on selection)
	nEtasThisEvent++;

	// count inLHCb acc spefici
	if(pythia.event[i].eta()>2.0 && pythia.event[i].eta()<5.0) {
	  nEtasInLHCbAccThisEvent++;

	  // fill pt for those in LHCb acc
	  etameson_pt_lhcbacc_h->Fill(pythia.event[i].pT());
	  // using E as short for ptot here
	  etameson_ptot_lhcbacc_h->Fill(pythia.event[i].e());


	  // check spread in angle between decay products

	  double min_eta_for_decay_prod = +9999;
	  double max_eta_for_decay_prod = -9999;

	  // can also get vector of the indices
	  vector<int> dlist = pythia.event[i].daughterList();
	  // set the min, max to the first one, then see if any is different
	  
	  // loop over the daughters
	  for(int id=0;id<dlist.size();id++) {
	    //get each daughter by index
	    //	  pythia.event[dlist[id]];
	    //	  cout << pythia.event[dlist[id]].id()<<endl;
	    if(pythia.event[dlist[id]].id()==13 || pythia.event[dlist[id]].id()==-13) {
	      // then its a muon daughter
	      if(pythia.event[dlist[id]].eta()<min_eta_for_decay_prod) {
		min_eta_for_decay_prod = pythia.event[dlist[id]].eta();	
	      }
	      if(pythia.event[dlist[id]].eta()>max_eta_for_decay_prod) {
		max_eta_for_decay_prod = pythia.event[dlist[id]].eta();
	      }		      
	    }
	    
	  } //end daughter loop

	  // calc spread in decay products
	  double spread_decay_prod = max_eta_for_decay_prod-min_eta_for_decay_prod;
	  etameson_decay_spread_h->Fill(spread_decay_prod);
	} //end of lHCb acc check for eta meson

	// for ANY eta, now check its daughters
	// duaghter list stored as 'daughter1 = index of first entry in list'
	// and daughter2= index of last
	//	int dghtr1 = pythia.event[i].daughter1();
	//	int dghtr2 = pythia.event[i].daughter2();
	// can also get vector of the indices
	vector<int> dlist = pythia.event[i].daughterList();
	//cout<< dghtr1 <<", "<<dghtr2 <<endl;
	//cout << dlist.size()<<endl;  
	// loop over the daughters
	// for this eta, count number of daughters in Acc
	int nDaughtersInAcc = 0;
	int nDaughtersInAccAndPT250 = 0;
	int nDaughtersInAccAndPT500 = 0;
	
	for(int id=0;id<dlist.size();id++) {
	  //get each daughter by index
	  //	  pythia.event[dlist[id]];
	  //cout << pythia.event[dlist[id]].id()<<endl;
	  if(pythia.event[dlist[id]].id()==13 || pythia.event[dlist[id]].id()==-13) {
	      // then its a muon daughter
	    // is this daughter in the LHCb acc?
	    if(pythia.event[dlist[id]].eta()>2.0 && pythia.event[dlist[id]].eta()<5.0) {
	      // tehn yes, its in LHCb acc
	      nDaughtersInAcc++;
	      // now check pT - first 0.25 GeV
	      if(pythia.event[dlist[id]].pT()>0.25)  {
		nDaughtersInAccAndPT250++;
		// then 0.5 GeV
		if(pythia.event[dlist[id]].pT()>0.5)  {
		  nDaughtersInAccAndPT500++;
		}
	      }
	      
	    }
	  }
	
	} //end daughter loop
	// now check: were all 4 daughters in acc?
	if(nDaughtersInAcc==4) {
	  //	  cout << "4 daughters in acc"<<endl;
	  nEtasWithAllFourDaughtersInAcc++;
	  if(nDaughtersInAccAndPT250==4) {
	    // then it passes 4 muon trig
	    nEtasWithAllFourDaughtersInAccAndFourMuons250++;
	  }
	  if(nDaughtersInAccAndPT500>=2) {
	    // then it passes dimuon trig
	    nEtasWithAllFourDaughtersInAccAndTwoMuons500++;
	  }
	  if(nDaughtersInAccAndPT250==4 || nDaughtersInAccAndPT500>=2) {
	    // then it passes either trigger!
	    nEtasWithAllFourDaughtersInAccAndEitherTrig++;
	  }

	  if(nDaughtersInAccAndPT250==4 && nDaughtersInAccAndPT500>=2) {
	    // then it passes BOTH triggers
	    nEtasWithAllFourDaughtersInAccAndBothTrigs++;
	  }
	}

      } //ends if this particle is an eta meson
	
	  // fill struct for quark
      //	  quark_in.fPt = pythia.event[i].pT();
      //	  quark_in.fEta = pythia.event[i].eta();
      //  quark_in.fPhi = pythia.event[i].phi();
      //  quark_in.fRapidity = pythia.event[i].y();
      //  quark_in.fPz = pythia.event[i].pz();
      //  quark_in.fPDGId = pythia.event[i].id();
      //  quark_in.fStatus = pythia.event[i].status();

      
      // get a final lepton whose mother is a W?
      if(abs(pythia.event[i].id())==13 && pythia.event[i].isFinal()) {

	// now find the Mother from mother1()
	int mother1_index = pythia.event[i].mother1();
	//cout << "Found final particle ID = " << pythia.event[i].id() <<", pt = " << pythia.event[i].pT()<<"eta = "<< pythia.event[i].eta()<<"; spin = "<<pythia.event[i].pol()<<"status = "<< pythia.event[i].status()<< ", Mother ID and status = "<< pythia.event[mother1_index].id() <<", "<< pythia.event[mother1_index].status()<<endl;

	//	if()
	// but some leptons have mothers that are not the W, eg due to some perturb
	// just take the highest pt final lepton?

	//so, if theres a lepton, push it onto array, then later sort and take highest pt
	selectedObjects.push_back(pythia.event[i]);
 	nSelectedObj++;

	
      }

	


      

      // try to find the outgoing particles from hard scatt (status 23 according to Pythia 8 manual)
      // but maybe this is not robust enough?
      // Actually, its probably -23, with convention that + means particles remaining, while - status means its not in final state
      // if(pythia.event[i].status()==23) {
      //	cout << "Outgoing hard scatt particle: ID = " << pythia.event[i].id() <<", pt = "<< pythia.event[i].pT() <<endl;
        
      //     if(pythia.event[i].isFinal() && pythia.event[i].pT()>10) {
      //	nSelectedObj++;

      //	CH::MyGenParticle thisObj(pythia.event[i]);
	//	cout << "Photon pt = " <<thisPhoton.GetPt() <<"id = " <<thisPhoton.GetPDGId() <<endl;
      // 	selectedObjects.push_back(thisObj);
	// get mother - this returns index in event record of mother
	//int motherId = 	pythia.event[i].mother1()
	//int mother2Id =  pythia.event[i].mother2()
	  // 'normal' decay case: mother1>0 and mother2==0	  	
      //      }

      
      // count Charged Particles in event (just because it was already there in example)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;


    } //end Particle loop within this event
    //cout << "N sel obj = " << nSelectedObj <<endl;
    // fill per-event hists
    //    nPhotons_h->Fill(nPhotons);


    //sort objects by pt, to select top 
    //    sort(selectedObjects.begin(),selectedObjects.end(),CH::OrderPythiaParticlesByPt);
    // for (int j=0;j<photons.size();j++) {
    //  cout << "Sorted photon pt = " <<photons[j].GetPt() <<endl;
    //}
    
    //set placeholders for tree filling

    // first default to unphysical
    //    InitMyGenParticleStruct(selObj1);
    //   InitMyGenParticleStruct(selObj2);
    
    // if(selectedObjects.size()>0) {
    //   lepton.fPt = selectedObjects[0].pT();
    //   lepton.fEta = selectedObjects[0].eta();
    //   lepton.fPhi = selectedObjects[0].phi();
    //   lepton.fRapidity = selectedObjects[0].y();
    //   lepton.fPz = selectedObjects[0].pz();
    //   lepton.fPDGId = selectedObjects[0].id();
    //   lepton.fStatus = selectedObjects[0].status();
    //   lepton.fSpin = selectedObjects[0].pol();
      
    // }
    
    /*if (selectedObjects.size()>1) {
      selObj2.fPt = selectedObjects[1].pT();
      selObj2.fEta = selectedObjects[1].eta();
      selObj2.fPhi = selectedObjects[1].phi();
      selObj2.fPDGId = selectedObjects[1].id();
      selObj2.fStatus = selectedObjects[1].status();
      } */
    //then fill tree
    outTree->Fill();

    
    // Fill charged multiplicity in histogram. 
    mult->Fill( nCharged );

    // Fill per-event histos
    nEta_h->Fill(nEtasThisEvent);
    nEta_lhcbacc_h->Fill(nEtasInLHCbAccThisEvent);
    
  } //end EVENT loop

  // Statistics on event generation.
  pythia.stat();

  // Show histogram. Possibility to close it.
  //  mult->Draw();
  //  nPhotons_h->Draw();
  // std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();

  //  cout <<"Out tree entries =" <<outTree->GetEntries()<<endl;
  cout << "N etas with all 4 daughters in LHCb acc = "<< nEtasWithAllFourDaughtersInAcc<<endl;
  cout << "N etas with all 4 daughters in LHCb acc and 2 muons pT>0.5 GeV= "<< nEtasWithAllFourDaughtersInAccAndTwoMuons500<<endl;
  cout << "N etas with all 4 daughters in LHCb acc and 4 muons pT>0.25 GeV= "<< nEtasWithAllFourDaughtersInAccAndFourMuons250<<endl;

  cout << "N etas with all 4 daughters in LHCb acc and EITHER triggers= "<<   nEtasWithAllFourDaughtersInAccAndEitherTrig++ <<endl;

  cout << "N etas with all 4 daughters in LHCb acc and both triggers= "<< nEtasWithAllFourDaughtersInAccAndBothTrigs<<endl;
  
  
  // Save histogram on file and close file.
  mult->Write();
  nPhotons_h->Write();
  nEta_h->Write();
  nEta_lhcbacc_h->Write();
  etameson_pt_lhcbacc_h->Write();
  etameson_ptot_lhcbacc_h->Write();
  etameson_decay_spread_h->Write();
  
  outTree->Write();
  outFile->Close();
  delete outFile;

  // Done.
  return 0;
}
