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
    /* MATPARSER TEST PROGRAM *
    insur::MatParser pars;
    insur::MatCalc calc;

    if (argc == 2) {
        if (pars.initMatCalc(argv[1], calc)) {
            calc.printInternals();
        }
        else std::cout << "main(): Reading of parameter file failed." << std::endl;
    }
    else std::cout << "Expecting a single argument: a material config file." << std::endl;*/
    /* GEOMETRY TEST PROGRAM */
    Tracker* tr = NULL;
    insur::InactiveSurfaces is;
    insur::TrackerActions ta;
    insur::Usher u;
    //insur::Vizard v;
    // tracker instance builds up active volumes
    if (argc == 3) tr = ta.createActiveSurfaces(argv[1], argv[2]);
    else if (argc == 2) tr = ta.createActiveSurfaces(argv[1]);
    else {
        std::cout << "The demo program needs at least a geometry configuration file." << std::endl;
        std::cout << "Either call './demo <geometry_config_file>', for geometry only" << std::endl;
        std::cout << "or './demo <geometry_config_file> <settings_file>' for fully dressed modules." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (tr) {
        // Usher instance builds inactive geometry around tracker elements
        is = u.arrange(*tr, is, argv[1]);
        // display result
        /*v.buildVisualization(*tr, is, true);
        v.display();
        std::string tmp(argv[1]);
        int pos = tmp.find('/') + 1;
        tmp = tmp.substr(pos);
        pos = tmp.find('.') + 1;
        tmp.erase(pos);
        tmp = tmp + "graph";
        v.writeNeighbourGraph(is, tmp);*/
        insur::MatCalcDummy mc;
        std::cout << std::endl << "Creating MaterialBudget instance...";
        insur::MaterialBudget mb(*tr, is);
        std::cout << "done." << std::endl << "Assigning dummy materials...";
        mb.materialsAll(mc);
        std::cout << "done." << std::endl << std::endl;
        mc.getMaterialTable().print();
        mb.print();
    }
    return (EXIT_SUCCESS);
}

