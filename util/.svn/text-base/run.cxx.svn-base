#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/DiskListEOS.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/GridDriver.h"

#include "TChain.h"

#include "H4l/H4lAnalysis.h"

using namespace std;

bool cmdline(int argc, char** argv, map<string,string> &opts);
void usage();

int main (int argc, char **argv) {
    map<string, string> opts;
    if (!cmdline(argc, argv, opts)) return 0;

    bool grid = false, wisc = false, inclusive = false, smearing = false, weights = false, doCR = false, debug = false;
    stringstream ss_grid; ss_grid << opts["grid"] << " "; ss_grid >> grid;
    stringstream ss_wisc; ss_wisc << opts["wisc"] << " "; ss_wisc >> wisc;
    stringstream ss_inc; ss_inc << opts["inclusive"] << " "; ss_inc >> inclusive;
    stringstream ss_sm; ss_sm << opts["smearing"] << " "; ss_sm >> smearing;
    stringstream ss_w; ss_w << opts["weights"] << " "; ss_w >> weights;
    stringstream ss_cr; ss_cr << opts["doCR"] << " "; ss_cr >> doCR;
    stringstream ss_d; ss_d << opts["debug"] << " "; ss_d >> debug;

    int numEvents = atoi(opts["numEvents"].c_str());

    string outDir_in = opts["outdir"];
    if (outDir_in == "") {
        cout << "Please enter output directory!  All output saved to $WORK/workarea/outData." << endl;
        return 1;
    }
    string outDir = outDir_in;
    if (!wisc) outDir = string(getenv("WORK")) + "/workarea/outData/" + outDir_in;

    xAOD::Init().ignore();

    SH::SampleHandler sh;

    if (opts["in"] == "") {
        cout << "Name of input file(s) or dataset(s) required!" << endl;
        return 2;
    }

    vector<string> ins;
    istringstream instream(opts["in"]);
    string in;
    while (std::getline(instream, in, ',')) {
        ins.push_back(in);
    }
    
    
    if (grid) {
        for (auto &d : ins) {
            SH::scanDQ2(sh, d);
        }
    }
    else {
        TChain chain("CollectionTree");
        for (auto &f : ins) {
            chain.Add(f.c_str());
        }
        sh.add(SH::makeFromTChain("xAOD", chain));
    }

    sh.setMetaString("nc_tree", "CollectionTree");

    sh.print();

    EL::Job job;
    job.sampleHandler(sh);

    H4lAnalysis *analysis = new H4lAnalysis(inclusive, smearing, weights, doCR, debug);
    job.algsAdd(analysis);

    if (numEvents != 0) job.options()->setDouble(EL::Job::optMaxEvents, numEvents);

    if (grid) {
        EL::GridDriver driver;
        driver.outputSampleName = "user.xju." + outDir_in + ".%in:name[2]%";
        driver.nFilesPerJob = 1;
        driver.submitOnly(job, outDir);
    }
    else {
        EL::DirectDriver driver;
        driver.submit(job, outDir);
    }

    return 0;
}

bool cmdline(int argc, char** argv, map<string,string> &opts) {
    opts.clear();

    // defaults
    opts["outdir"] = "";
    opts["in"] = "";
    opts["grid"] = "0";
    opts["wisc"] = "0";
    opts["inclusive"] = "0";
    opts["smearing"] = "0";
    opts["weights"] = "0";
    opts["numEvents"] = "0";
    opts["doCR"] = "0";
    opts["debug"] = "0";

    for (int i=1;i<argc;i++) {

        string opt=argv[i];

        if (opt=="--help") {usage(); return false;}

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
        <<"-outdir (required!)\n"
        <<"-in (required!)\n"
        <<"-grid <0/1> [0]\n"
        <<"-wisc <0/1> [0]\n"
        <<"-inclusive <0/1> [0]\n"
        <<"-smearing <0/1> [0]\n"
        <<"-weight <0/1> [0]\n"
        <<"-numEvents [all]\n"
        <<"-doCR <0/1> [0]\n"
        <<"-debug <0/1> [0]\n"
        <<endl;

    return;
}
