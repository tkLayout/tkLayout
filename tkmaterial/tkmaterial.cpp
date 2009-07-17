//
// File:   tkmaterial.cpp
// Author: ndemaio
//
// Created on June 8, 2009, 10:51 AM
//

#include <stdlib.h>
#include <Squid.h>
#include <boost/algorithm/string/trim.hpp>
namespace balgo = boost::algorithm;
int main(int argc, char** argv) {
    std::string geomfile, settingsfile, matfile, htmlout, rootout, graphout;
    bool switch_processed = false, files_processed = false, u = false, m = false, d = false, h = false, r = false, g = false, t = false;
    int cfiles = 0, tracks = 0, pos = 0, i = 0;
    
    // idiot usage check
    if (argc < 2) {
        // print usage message
        std::cout << "Error: tkmaterial needs at least a geometry configuration file to run." << std::endl;
        std::cout << "The full call syntax is as follows:" << std::endl << std::endl;
        std::cout << "./tkmaterial [-umd] geomfile [settingsfile] [materialfile] [-t n_of_tracks] [-h htmlfile] [-r rootfile] [-g graphfile]";
        std::cout << std::endl << std::endl << "u   print geometry summary after volume creation" << std::endl;
        std::cout << "m   print material budget summary after material assignment" << std::endl;
        std::cout << "d   write detailed geometry to root file - WARNING: may cause root to crash later!" << std::endl;
        std::cout << "geomfile   geometry configuration file" << std::endl << "settingsfile   module settings configuration file" << std::endl;
        std::cout << "materialfile   material configuration file" << std::endl << "n_of_tracks   number of tracks used for analysis" << std::endl;
        std::cout << "htmlfile   output file for material budget analysis" << std::endl;
        std::cout << "rootfile   output file for ROOT (geometry)" << std::endl << "graphfile   output file for feeder/neighbour graph" << std::endl;
        return (EXIT_FAILURE);
    }
    
    // argument processing
    i++;
    // argument can be option switch or geomfile
    if (argv[i][0] == '-') {
        std::string tmp = argv[i];
        tmp = tmp.substr(1);
        balgo::trim(tmp);
        for (std::string::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
            switch (*iter) {
                case 'u': u = true;
                break;
                case 'm': m = true;
                break;
                case 'd': d = true;
                break;
                default: std::cerr << "Error: unknown parameter " << *iter << ". Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
        }
        switch_processed = true;
    }
    else {
        geomfile = argv[i];
        cfiles++;
    }
    i++;
    if (i < argc) {
        // next argument can be geomfile, settingsfile, tracks or outfileswitch
        if (switch_processed) {
            if (argv[i][0] == '-') {
                std::cerr << "Error: tkmaterial needs at least a geometry configuration file to run." << std::endl;
                std::cout << "Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            else {
                geomfile = argv[i];
                cfiles++;
            }
        }
        else {
            if (argv[i][0] == '-') {
                if (argv[i][1]){
                    if (argv[i][1] == 'h') {
                        if (h) {
                            std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                        else {
                            if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                                htmlout = argv[i+1];
                                i++;
                            }
                            h = true;
                        }
                    }
                    else if (argv[i][1] == 'r') {
                        if (r) {
                            std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                        else {
                            if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                                rootout = argv[i+1];
                                i++;
                            }
                            r = true;
                        }
                    }
                    else if (argv[i][1] == 'g') {
                        if (g) {
                            std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                        else {
                            if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                                graphout = argv[i+1];
                                i++;
                            }
                            g = true;
                        }
                    }
                    else if (argv[i][1] == 't') {
                        if (t) {
                            std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                        else {
                            if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                                tracks = atoi(argv[i+1]);
                                i++;
                                t = true;
                            }
                            else {
                                std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                                return (EXIT_FAILURE);
                            }
                        }
                    }
                    else {
                        std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
                files_processed = true;
            }
            else {
                settingsfile = argv[i];
                cfiles++;
            }
        }
    }
    i++;
    if (i < argc) {
        // next argument can be settingsfile, matfile, outfileswitch or tracks
        if (argv[i][0] != '-') {
            if (cfiles == 1) {
                settingsfile = argv[i];
                cfiles++;
            }
            else {
                matfile = argv[i];
                cfiles++;
                files_processed = true;
            }
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        // next argument can be outfileswitch, tracks or materialfile
        if (argv[i][0] != '-') {
            if (files_processed) {
                std::cerr << "Error: unexpected config file found. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            else {
                matfile = argv[i];
                cfiles++;
                files_processed = true;
            }
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        // next argument can be outfileswitch or tracks (up to four to go)
        if (argv[i][0] != '-') {
            std::cerr << "Error: too many config files. Aborting tkmaterial." << std::endl;
            return (EXIT_FAILURE);
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            std::cout << "tracks = " << tracks << std::endl;
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        // next argument can be outfileswitch or tracks (up to three to go)
        if (argv[i][0] != '-') {
            std::cerr << "Error: too many config files. Aborting tkmaterial." << std::endl;
            return (EXIT_FAILURE);
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        // next argument can be outfileswitch or tracks (up to two to go)
        if (argv[i][0] != '-') {
            std::cerr << "Error: too many config files. Aborting tkmaterial." << std::endl;
            return (EXIT_FAILURE);
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        // next argument can be outfileswitch or tracks (last one)
        if (argv[i][0] != '-') {
            std::cerr << "Error: too many config files. Aborting tkmaterial." << std::endl;
            return (EXIT_FAILURE);
        }
        else {
            if (argv[i][1]){
                if (argv[i][1] == 'h') {
                    if (h) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            htmlout = argv[i+1];
                            i++;
                        }
                        h = true;
                    }
                }
                else if (argv[i][1] == 'r') {
                    if (r) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            rootout = argv[i+1];
                            i++;
                        }
                        r = true;
                    }
                }
                else if (argv[i][1] == 'g') {
                    if (g) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-')) {
                            graphout = argv[i+1];
                            i++;
                        }
                        g = true;
                    }
                }
                else if (argv[i][1] == 't') {
                    if (t) {
                        std::cerr << "Error: redefinition of parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                        return (EXIT_FAILURE);
                    }
                    else {
                        if (((i + 1) < argc) && (argv[i+1][0] != '-') && (atoi(argv[i+1]) != 0)) {
                            tracks = atoi(argv[i+1]);
                            i++;
                            t = true;
                        }
                        else {
                            std::cerr << "Error parsing numeric parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                            return (EXIT_FAILURE);
                        }
                    }
                }
                else {
                    std::cerr << "Error: unknown parameter " << argv[i] << ". Aborting tkmaterial." << std::endl;
                    return (EXIT_FAILURE);
                }
            }
            else {
                std::cerr << "Error: unknown parameter. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
            files_processed = true;
        }
    }
    i++;
    if (i < argc) {
        std::cerr << "Error: too many arguments. The tkmaterial call syntax is as follows:" << std::endl;
        std::cout << "./tkmaterial [-umd] geomfile [settingsfile [materialfile]] [-t n_of_tracks] [-h [htmlfile]] [-r [rootfile]] [-g [graphfile]]";
        std::cout << std::endl << std::endl << "u   print geometry summary after volume creation" << std::endl;
        std::cout << "m   print material budget summary after material assignment" << std::endl;
        std::cout << "d   write detailed geometry to root file - WARNING: may cause root to crash later!" << std::endl;
        std::cout << "geomfile   geometry configuration file" << std::endl << "settingsfile   module settings configuration file" << std::endl;
        std::cout << "materialfile   material configuration file" << std::endl << "n_of_tracks   number of tracks used for analysis" << std::endl;
        std::cout << "htmlfile   output file for material budget analysis" << std::endl;
        std::cout << "rootfile   output file for ROOT (geometry)" << std::endl << "graphfile   output file for feeder/neighbour graph" << std::endl;
        return (EXIT_FAILURE);
    }
    
    //DEBUG: print internal status
    /*std::cout << std::endl << "Internal status after parsing " << argc << " arguments (i = " << i << "):" << std::endl;
     * std::cout << "switch_processed = " << (switch_processed ? "true" : "false") << ", files_processed = ";
     * std::cout << (files_processed ? "true" : "false") << ", cfiles = " << cfiles << std::endl;
     * std::cout << "usher_verbose = " << (u ? "true" : "false") << ", mat_verbose = " << (m ? "true" : "false");
     * std::cout << ", detailed = " << (d ? "true" : "false") << std::endl;
     * std::cout << "geomfile = " << geomfile << ", settingsfile = " << settingsfile << ", matfile = " << matfile << "." <<std::endl;
     * std::cout << "htmlflag = " << (h ? "true" : "false") << ", htmlout = " << htmlout << "." << std::endl;
     * std::cout << "rootflag = " << (r ? "true" : "false") << ", rootout = " << rootout << "." << std::endl;
     * std::cout << "graphflag = " << (g ? "true" : "false") << ", graphout = " << graphout << "." << std::endl;
     * std::cout << "trackflag = " << (t ? "true" : "false") << ", tracks = " << tracks << std::endl << std::endl;*/
    
    // here comes the heavy lifting...
    insur::Squid s;
    if (h) {
        if (s.buildFullSystem(geomfile, settingsfile, matfile, u, m)) {
            if (htmlout.empty()) {
                htmlout = geomfile;
                pos = htmlout.find_last_of('/');
                if (pos != (int)htmlout.npos) {
                    pos++;
                    htmlout = htmlout.substr(pos);
                }
                pos = htmlout.find('.');
                if (pos != (int)htmlout.npos) {
                    if (pos > (int)htmlout.npos - 1) htmlout.erase(pos + 1);
                }
                else htmlout.push_back('.');
                htmlout = htmlout + "html";
                std::cout << "HTML file will be written to " << htmlout << std::endl;
            }
            if (tracks != 0) {
                std::cout << "Calling analyzer with " << tracks << " tracks." << std::endl;
                if (!s.analyzeMaterialBudget(htmlout, tracks)) return (EXIT_FAILURE);
            }
            else {
                std::cout << "Calling analyzer with the default number of tracks." << std::endl;
                if (!s.analyzeMaterialBudget(htmlout)) return (EXIT_FAILURE);
            }
            if (r) {
                if (rootout.empty()) {
                    rootout = geomfile;
                    pos = rootout.find_last_of('/');
                    if (pos != (int)rootout.npos) {
                        pos++;
                        rootout = rootout.substr(pos);
                    }
                    pos = rootout.find('.');
                    if (pos != (int)rootout.npos) {
                        if (pos > (int)rootout.npos - 1) rootout.erase(pos + 1);
                    }
                    else rootout.push_back('.');
                    rootout = rootout + "root";
                    std::cout << "ROOT file will be written to " << rootout << std::endl;
                }
                if (!s.analyzeGeometry(rootout, !d)) return (EXIT_FAILURE);
                if (g) {
                    if (graphout.empty()) {
                        graphout = geomfile;
                        pos = graphout.find_last_of('/');
                        if (pos != (int)graphout.npos) {
                            pos++;
                            graphout = graphout.substr(pos);
                        }
                        pos = graphout.find('.');
                        if (pos != (int)graphout.npos) {
                            if (pos > (int)graphout.npos - 1) graphout.erase(pos + 1);
                        }
                        else graphout.push_back('.');
                        graphout = graphout + "graph";
                        std::cout << "Graph file will be written to " << graphout << std::endl;
                    }
                    if (!s.analyzeNeighbours(graphout)) return (EXIT_FAILURE);
                }
            }
        }
        else return (EXIT_FAILURE);
    }
    else {
        if (r) {
            if (!s.buildTracker(geomfile)) return (EXIT_FAILURE);
            if (settingsfile.empty()) std::cout << "Warning: using tracker geometry without a settings file to dress it from." << std::endl;
            else {
                if (!s.dressTracker(settingsfile)) return (EXIT_FAILURE);
            }
            if (s.buildInactiveSurfaces(u)) {
                if (rootout.empty()) {
                    rootout = geomfile;
                    pos = rootout.find_last_of('/');
                    if (pos != (int)rootout.npos) {
                        pos++;
                        rootout = rootout.substr(pos);
                    }
                    pos = rootout.find('.');
                    if (pos != (int)rootout.npos) {
                        if (pos > (int)rootout.npos - 1) rootout.erase(pos + 1);
                    }
                    else rootout.push_back('.');
                    rootout = rootout + "root";
                    std::cout << "ROOT file will be written to " << rootout << std::endl;
                }
                if (!s.analyzeGeometry(rootout, !d)) return (EXIT_FAILURE);
                if (g) {
                    if (graphout.empty()) {
                        graphout = geomfile;
                        pos = graphout.find_last_of('/');
                        if (pos != (int)graphout.npos) {
                            pos++;
                            graphout = graphout.substr(pos);
                        }
                        pos = graphout.find('.');
                        if (pos != (int)graphout.npos) {
                            if (pos > (int)graphout.npos - 1) graphout.erase(pos + 1);
                        }
                        else graphout.push_back('.');
                        graphout = graphout + "graph";
                        std::cout << "Graph file will be written to " << graphout << std::endl;
                    }
                    if (!s.analyzeNeighbours(graphout)) return (EXIT_FAILURE);
                }
            }
            else return (EXIT_FAILURE);
            if (m && !matfile.empty()) s.createMaterialBudget(matfile, m);
        }
        else if (g) {
            if (!s.buildTracker(geomfile)) return (EXIT_FAILURE);
            if (settingsfile.empty()) std::cout << "Warning: using tracker geometry without a settings file to dress it from." << std::endl;
            else {
                if (!s.dressTracker(settingsfile)) return (EXIT_FAILURE);
            }
            if (s.buildInactiveSurfaces(u)) {
                if (graphout.empty()) {
                    graphout = geomfile;
                    pos = graphout.find_last_of('/');
                    if (pos != (int)graphout.npos) {
                        pos++;
                        graphout = graphout.substr(pos);
                    }
                    pos = graphout.find('.');
                    if (pos != (int)graphout.npos) {
                        if (pos > (int)graphout.npos - 1) graphout.erase(pos + 1);
                    }
                    else graphout.push_back('.');
                    graphout = graphout + "graph";
                    std::cout << "Graph file will be written to " << graphout << std::endl;
                }
                if (!s.analyzeNeighbours(graphout)) return (EXIT_FAILURE);
            }
            if (m && !matfile.empty()) s.createMaterialBudget(matfile, m);
        }
        else {
            //no output files required
            switch (cfiles) {
                case 1 : if (!s.buildInactiveSurfaces(geomfile, u)) return (EXIT_FAILURE);
                break;
                case 2 : if (!s.buildInactiveSurfaces(geomfile, settingsfile, u)) return (EXIT_FAILURE);
                break;
                case 3 : if (!s.buildFullSystem(geomfile, settingsfile, matfile, u, m)) return (EXIT_FAILURE);
                break;
                default: std::cerr << "Something truly strange happened during processing: the command line passed the parsing stage ";
                std::cerr << "but the number of input files is not between 1 and 3. Aborting tkmaterial." << std::endl;
                return (EXIT_FAILURE);
            }
        }
    }
    std::cout << "Done." << std::endl;
    return (EXIT_SUCCESS);
}

