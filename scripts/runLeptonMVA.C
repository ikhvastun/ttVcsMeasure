{
    TString outfileName( "TMVA.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

    
    dataloader->AddVariable( "pt", 'F' );
    dataloader->AddVariable( "eta", 'F' );
    dataloader->AddVariable( "trackMult", 'I' );
    dataloader->AddVariable( "miniIsoCharged", 'F' );
    dataloader->AddVariable( "miniIsoNeutral", 'F' );
    dataloader->AddVariable( "ptrel", 'F' );
    dataloader->AddVariable( "min(ptratio,1.5)", 'F' );
    dataloader->AddVariable( "max(jetbtagCSV,0)", 'F' );
    dataloader->AddVariable( "sip3d", 'F' );
    dataloader->AddVariable( "log(abs(dxy))", 'F' );
    dataloader->AddVariable( "log(abs(dz))", 'F' );
    dataloader->AddVariable( "eleMVA", 'F');

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;

    TString fname = "fileDummy.root";
    TFile *input = TFile::Open( fname );

    TTree *signalTree     = (TTree*)input->Get("signalTree");
    TTree *bkgTree = (TTree*)input->Get("bkgTree");

    dataloader->AddSignalTree    ( signalTree,     signalWeight     );
    dataloader->AddBackgroundTree( bkgTree, backgroundWeight );

    dataloader->SetSignalWeightExpression( "weight" );
    dataloader->SetBackgroundWeightExpression( "weight" );

    // _weight > 0 && nJLoc > 1 && _met > 50
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=100000:nTrain_Background=100000:nTest_Signal=100000:nTest_Background=100000:NormMode=None:SplitMode=Random:!V" );


    //factory->BookMethod( TMVA::Types::kCuts, "Cuts",
    //                       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
    //factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
    //                       "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=5000" );
    
    //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
    //                       "!H:!V:NTrees=850:MinNodeSize=1%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1:IgnoreNegWeightsInTraining=True" );
    
    // Boosted Decision Trees
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.01:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=200:MaxDepth=4:IgnoreNegWeightsInTraining=True" );
    // Bagging
    //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
    //                       "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining=True" );
                           

    // Decorrelation + Adaptive Boost
    //factory->BookMethod( TMVA::Types::kBDT, "BDTD",
    //                       "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate:IgnoreNegWeightsInTraining:Pray" );

    // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    //factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
    //                       "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining:Pray" );

    //factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------
    // Save the output
    outputFile->Close();
 
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
 
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
 
}
