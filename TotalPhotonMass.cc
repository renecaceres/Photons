{
    gROOT->Reset();
    gSystem->Load("/home/rene/SiLab/ClasTool/slib/Linux/libClasTool.so");
    gSystem->Load("/home/rene/SiLab/Analyser/analysis_lib/slib/libTIdentificator.so");
    TClasTool *ct=new TClasTool();
    TIdentificator *partident=new TIdentificator(ct);
    ct->InitDSTReader("ROOTDSTR"); //Inicializa el DSTReader

    //Lista de Archivos
  
    ct->AddFile("/home/rene/SiLab/Photons/out.root");  //Agrega el archivo a la lista por analizar del DSTReader
  
    Int_t nentries = (Int_t)ct->GetEntries();
    cout<<nentries<<endl;

    //Archivos y Arboles

    TFile *plots = new TFile("tree.root","RECREATE");
    TTree *tree_photon = new TTree("tree_photon","Photons");
    TH1F *h1 = new TH1F("Mass (Momentum)", "Total Photons Mass (Momentum) [GeV]", 100,0,1.2);
    TH1F *h2 = new TH1F("Mass (Energy)", "Total Photons Mass (Energy) [GeV]", 100,0,1.2);
  
    Int_t nmb_photon=20;
    Float_t etot[nmb_photon],px[nmb_photon],py[nmb_photon],pz[nmb_photon],moment[nmb_photon],etot[nmb_photon],mass;

    tree_photon->Branch("nmb_photon",&nmb_photon,"nmb_photon/I");
    tree_photon->Branch("etot",&etot,"etot[nmb_photon]/F");
    tree_photon->Branch("px",&px,"px[nmb_photon]/F");
    tree_photon->Branch("py",&py,"py[nmb_photon]/F");
    tree_photon->Branch("pz",&pz,"pz[nmb_photon]/F");
    tree_photon->Branch("mass",&mass,"mass/F");

    Int_t number;
    for(Int_t i=0;i<ct->GetEntriesFast();i++){ //se usa GetEntriesFast pq GetEntries ya fue llamado, más rápido
        ct->Next();  //Importante!!!
        number=ct->GetNRows("EVNT"); //Número de filas en EVNT 
        if(number>0){
            nmb_photon=0;
            for(int k=1;k<number;k++){
                if(partident->GetCategorization(k).CompareTo("photon")==0) {
                    px[nmb_photon]=partident.Px(k,0);
                    py[nmb_photon]=partident.Py(k,0);
                    pz[nmb_photon]=partident.Pz(k,0);
                    etot[nmb_photon]=partident.Etot(k);
                    nmb_photon++;
                }
            }
            if(nmb_photon==2) {
                mass=sqrt(pow(sqrt(px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1])+sqrt(px[2]*px[2]+py[2]*py[2]+pz[2]*pz[2]),2)-pow(px[1]+px[2],2)-pow(py[1]+py[2],2)-pow(pz[1]+pz[2],2));
                tree_photon->Fill();
                h1->Fill(mass);
                h2->Fill(sqrt(pow(etot[1],2)+pow(etot[2],2)));
            }
        }
    }
    plots->Write();
    return 0;
}
