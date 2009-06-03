//
// File:   demo.cpp.cc
// Author: ndemaio
//
// Created on October 13, 2008, 11:55 AM
//

#include <stdlib.h>
#include <iostream>
#include <tracker.hh>
#include <MatCalc.h>
#include <MatCalcDummy.h>
#include <MatParser.h>
#include <MaterialTable.h>
#include <TrackerActions.h>
#include <MaterialBudget.h>
#include <InactiveSurfaces.h>
#include <Usher.h>
#include <Vizard.h>
int main(int argc, char** argv) {
    /* GEOMETRY TEST PROGRAM */
    insur::MatParser p;
    Tracker* tr = NULL;
    insur::InactiveSurfaces is;
    insur::TrackerActions ta;
    insur::Usher u;
    insur::MatCalc c;
    insur::Vizard v;
    insur::Analyzer a;
    // tracker instance builds up active volumes
    if ((argc == 3) || (argc == 4)) tr = ta.createActiveSurfaces(argv[1], argv[2]);
    else if (argc == 2) tr = ta.createActiveSurfaces(argv[1]);
    else {
        std::cout << "The demo program needs at least a geometry configuration file." << std::endl;
        std::cout << "Either call './demo <geometry_config_file>', for geometry only" << std::endl;
        std::cout << "or './demo <geometry_config_file> <settings_file>' for fully dressed modules." << std::endl;
        std::cout << "./demo <geometry_config_file> <settings_file> <materials_list> creates the" << std::endl;
        std::cout << "complete material budget." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (tr) {
        // Usher instance builds inactive geometry around tracker elements
        is = u.arrange(*tr, is, argv[1], true);
        if (argc == 4) {
            // MatParser instance parses material file and initialises MatCalc instance with its contents
            if (p.initMatCalc(argv[3], c)) {
                // MaterialBudget instance builds up internal structures in constructor
                insur::MaterialBudget m(*tr, is);
                std::cout << std::endl << "Assigning materials...";
                // MaterialBudget instance assigns materials from material file to geometry using MatCalc instance
                m.materialsAll(c);
                std::cout << "done." << std::endl << std::endl;
                m.print();
                // Analyse the material budget over eta
                std::cout << "Analyzing material budget...";
                a.analyzeMaterialBudget(m);
                std::cout << "done." << std::endl << std::endl;
            }
            else std::cout << "main(): Reading of material parameter file failed." << std::endl;
        }
        // display result: write summary to .html file in /matsum, write (simplified) geometry to .root file in rootfiles/,
        //                        write neighbour graph to .graph file in graphs/
        std::string tmp(argv[1]);
        int pos = tmp.find('/') + 1;
        if (argc == 4) {
            tmp = tmp.substr(pos);
            pos = tmp.find('.') + 1;
            tmp.erase(pos);
            tmp = tmp + "html";
            v.histogramSummary(a, tmp);
            pos = tmp.find('/') + 1;
        }
        tmp = tmp.substr(pos);
        pos = tmp.find('.') + 1;
        tmp.erase(pos);
        tmp = tmp + "root";
        v.display(*tr, is, tmp, true);
        pos = tmp.find('/') + 1;
        tmp = tmp.substr(pos);
        pos = tmp.find('.') + 1;
        tmp.erase(pos);
        tmp = tmp + "graph";
        v.writeNeighbourGraph(is, tmp);
    }
    return (EXIT_SUCCESS);
}

