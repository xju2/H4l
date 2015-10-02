#include <iostream>

#include "EventLoopGrid/GridDriver.h"

using namespace std;

void usage();

int main (int argc, char **argv) {
    if (argc < 2) {
        cout << "Please specify action!" << endl;
        usage();
        return 0;
    }
    if (argc > 3) {
        cout << "Too much input!" << endl;
        usage();
        return 0;
    }

    EL::GridDriver driver;
    string outDir = string(getenv("WORK")) + "/workarea/outData/";

    if (argc == 2) {
        string option = argv[1];
        if (option == "killAll") {
            driver.killAll();
            return 0;
        }
        else if (option == "listActive") {
            driver.listActive();
            return 0;
        }
        else {
            cout << "Invalid input!" << endl;
            usage();
            return 0;
        }
    }
    else if (argc == 3) {
        string option = argv[1];
        string name = argv[2];
        if (option == "status") {
            driver.status(outDir + name);
            return 0;
        }
        else if (option == "kill") {
            driver.kill(outDir + name);
            return 0;
        }
        else {
            cout << "Invalid input!" << endl;
            usage();
            return 0;
        }
    }

    return 0;
}

void usage() {
    cout << "gangaControl usage:" << endl;
    cout << "gangaControl [option] [name]" << endl;
    cout << "Options that do not require a name: killAll, listActive" << endl;
    cout << "Options that do require a name: kill, status" << endl;
    return;
}
