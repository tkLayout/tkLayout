// 
// File:   xmltest.cpp
// Author: ndemaio
//
// Created on September 11, 2009, 11:23 AM
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <tk2CMSSW.h>

//
// 
//
int main(int argc, char** argv) {
    insur::tk2CMSSW tc;
    std::ostringstream stream;
    insur::tk2CMSSW::LogicalInfo li;
    std::vector<LogicalInfo> lv;
    /*std::string name, solid, material;
    double d1, d2, d3, d4;
    int i1;
    //elementaryMaterial()
    name = "Mat";
    d1 = 123.4;
    i1 = 10;
    d2 = 56.7;
    tc.elementaryMaterial(name, d1, i1, d2, stream);
    //logicalPart()
    name = "Logic";
    solid = "Vol";
    material = "Mat";
    tc.logicalPart(name, solid, material, stream);
    //box()
    name = "Boxx";
    d1 = 15;
    d2 = 20;
    d3 = 25;
    tc.box(name, d1, d2, d3, stream);
    //tubs()
    name = "Tube";
    d1 = 30;
    d2 = 35;
    d3 = 40;
    tc.tubs(name, d1, d2, d3, stream);
    //trapezoid()
    name = "Trap";
    d1 = 12;
    d2 = 23;
    d3 = 4.5;
    d4 = 67;
    tc.trapezoid(name, d1, d2, d3, d4, stream);
    //rotation
    insur::tk2CMSSW::Rotation rot;
    rot.name = "Rot";
    rot.phix = 1.1;
    rot.phiy = 1.2;
    rot.phiz = 1.3;
    rot.thetax = 2.1;
    rot.thetay = 2.2;
    rot.thetaz = 2.3;
    //translation
    insur::tk2CMSSW::Translation trans;
    trans.dx = 5.0;
    trans.dy = 6.0;
    trans.dz = 7.0;
    //posPart()
    i1 = 0;
    solid = "Parent";
    material = "Child";
    tc.posPart(i1, solid, material, rot, trans, stream);*/
    li.name_tag = "Top";
    li.shape_tag = "Tube";
    li.material_tag = "Air";
    lv.push_back(li);
    li.name_tag = "Middle";
    li.shape_tag = "Box";
    li.material_tag = "CF";
    lv.push_back(li);
    li.name_tag = "Leaf";
    li.shape_tag = "Box";
    li.material_tag = "Si";
    tc.logicalPartSection(lv, "Logic", stream);
    std::cout << stream.str() << std::endl;
    return (EXIT_SUCCESS);
}

