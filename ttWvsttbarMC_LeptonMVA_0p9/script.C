{
    TString outfileName( "TMVA.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    factory->AddVariable( "HTLoc", 'F' );
    factory->AddVariable( "nJLoc", 'I' );
    factory->AddVariable( "nBLoc", 'I' );
    factory->AddVariable( "_met", 'F' );
    //factory->AddVariable( "mWjetjet", 'F' );
    //factory->AddVariable( "mtopWbjet", 'F' );
    factory->AddVariable( "minDeltaR", 'F' );
    //factory->AddVariable( "ele_mll", 'F' );
    factory->AddVariable( "mt", 'F' );
    factory->AddVariable( "mtlow", 'F' );
    factory->AddVariable( "leadpt", 'F' );
    factory->AddVariable( "trailpt", 'F' );
    factory->AddVariable( "leadingJetPt", 'F' );
    factory->AddVariable( "trailJetPt", 'F' );
    //factory->AddVariable("leptonsCharge",'F');
    //factory->AddVariable( "leadingJetCSV", 'F' );
    //factory->AddVariable( "trailJetCSV", 'F' );  

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;

    TString fname = "fileDummy.root";
    TFile *input = TFile::Open( fname );

    TTree *signalTree     = (TTree*)input->Get("signalTree");
    TTree *bkgTree = (TTree*)input->Get("bkgTree");

    factory->AddSignalTree    ( signalTree,     signalWeight     );
    factory->AddBackgroundTree( bkgTree, backgroundWeight );

    factory->SetSignalWeightExpression( "_weight" );
    factory->SetBackgroundWeightExpression( "_weight" );

    // _weight > 0 && nJLoc > 1 && _met > 50
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:NormMode=None:SplitMode=Random:!V" );


    //factory->BookMethod( TMVA::Types::kCuts, "Cuts",
    //                       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
    //factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
    //                       "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=5000" );
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=1%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1:IgnoreNegWeightsInTraining" );
    // Boosted Decision Trees
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
    // Bagging
    factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

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

    delete factory;
}
