void Run(char* inputFile){

    gSystem->Load("/MG5/Delphes/libDelphes.so");
    gSystem->AddIncludePath("-I/MG5/Delphes/");
    //gSystem->CompileMacro("HelperClasses.cxx","k");
    gSystem->CompileMacro("cutflow.c","k");
    //gSystem->AddIncludePath("-I/home/ytchou/MadGraph5_v1_5_14/Delphes/");^M
    //std::cout << gSystem->GetIncludePath();
    cutflow(inputFile);


}
