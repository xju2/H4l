#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <fstream>

#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"

#include "H4l/H4lAnalysis.h"
#include <TChain.h>

using namespace std;

bool cmdline(int argc, char** argv, map<string,string> &opts);
void usage();

int main( int argc, char* argv[] ) 
{
    map<string, string> opts;
    if (!cmdline(argc, argv, opts)) return 0;
    if(opts["in"] == ""){
        cout << "ERROR: please indicate input by \'-in toberun.txt\'" << endl;
        return 1;
    }

    bool inclusive = (bool)atoi(opts["inclusive"].c_str());
    bool smearing = (bool)atoi(opts["smearing"].c_str()); 
    bool weights = (bool)atoi(opts["weights"].c_str()); 
    bool doCR = (bool)atoi(opts["doCR"].c_str()); 
    bool debug = (bool)atoi(opts["debug"].c_str()); 
    bool is_data = (bool)atoi(opts["isData"].c_str());
    bool is_atlfast = (bool) atoi(opts["atlFast"].c_str());

    int numEvents = atoi(opts["numEvts"].c_str());

    std::string submitDir = "out";
    if (opts["out"] != ""){
        submitDir = opts["out"];
    }
    xAOD::Init().ignore();
    SH::SampleHandler sh;
    cout << "--------------------" << endl;
    cout << "       summary      " << endl;
    cout << "isData: " << is_data << endl;
    cout << "numEvents: " << numEvents << endl;
    cout << "do CR: " << doCR << endl;
    cout << "is AtlFast: " << is_atlfast << endl;
    cout << "--------------------" << endl;
    
    TChain chain("CollectionTree","CollectionTree");
    ifstream input_file(opts["in"].c_str());
    string name_file;
    while (input_file>>name_file){
        chain.Add(name_file.c_str());
    }

    SH::Sample* testsample = SH::makeFromTChain("CollectionTree", chain);
    sh.add(testsample);

    sh.print();
    EL::Job job;
    job.sampleHandler( sh );
    if(debug){
        cout << "In DEBUG mode!" << endl;
    }
    H4lAnalysis *analysis = new H4lAnalysis(inclusive, smearing, weights, 
            doCR, is_data, is_atlfast, debug);
    job.algsAdd(analysis);

    // set job
    if (numEvents != 0) job.options()->setDouble(EL::Job::optMaxEvents, numEvents);

    EL::DirectDriver driver;
    // try {
    driver.submit( job, submitDir );
    // } catch (...) {
    //   cout <<"out directory \'" << submitDir << "\' exists, use another one"<<endl;
    // }

    return 0;
}

bool cmdline(int argc, char** argv, map<string,string> &opts) {
    opts.clear();

    // defaults
    opts["out"] = "";
    opts["in"] = "";
    opts["inclusive"] = "0";
    opts["smearing"] = "1";
    opts["weights"] = "0";
    opts["numEvts"] = "0";
    opts["doCR"] = "0";
    opts["debug"] = "0";
    opts["isData"] = "0";
    opts["atlFast"] = "0";

    for (int i=1;i<argc;i++) {

        string opt=argv[i];

        if (opt.find("help") != string::npos) {
            usage(); return false;
        }

        if (0 != opt.find("-")) {
            cout<<"ERROR: options start with '-'!"<<endl;
            return false;
        }
        opt.erase(0,1);
        if (opts.find(opt) == opts.end()) {
            cout<<"ERROR: invalid option '"<<opt<<"'!"<<endl;
            return false;
        }
        string nxtopt=argv[i+1];
        if (0 == nxtopt.find("-") || i+1 >= argc) {
            cout<<"ERROR: option '"<<opt<<"' requires value!"<<endl;
            return false;
        }

        opts[opt] = nxtopt;
        i++;
    }

    return true;
}

void usage()
{
    cout<<"USAGE: run [-option value]\n\n"
        <<"options [default]:\n\n"
        <<"-out (required!)\n"
        <<"-in (required!)\n"
        <<"-inclusive <0/1> [0]\n"
        <<"-smearing <0/1> [1]\n"
        <<"-weight <0/1> [0]\n"
        <<"-numEvts [all]\n"
        <<"-doCR <0/1> [0]\n"
        <<"-isData <0/1> [0]\n"
        <<"-atlFast <0/1> [0]\n"
        <<"-debug <0/1> [0]\n"
        <<endl;

    return;
}
