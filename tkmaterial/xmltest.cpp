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
    insur::tk2CMSSW::Element el;
    insur::tk2CMSSW::Composite co;
    insur::tk2CMSSW::LogicalInfo li;
    insur::tk2CMSSW::ShapeInfo si;
    std::vector<insur::tk2CMSSW::Element> ev;
    std::vector<insur::tk2CMSSW::Composite> cv;
    std::vector<insur::tk2CMSSW::LogicalInfo> lv;
    std::vector<insur::tk2CMSSW::ShapeInfo> sv;
    std::pair<std::string, double> p;
    /*//rotation
     * insur::tk2CMSSW::Rotation rot;
     * rot.name = "Rot";
     * rot.phix = 1.1;
     * rot.phiy = 1.2;
     * rot.phiz = 1.3;
     * rot.thetax = 2.1;
     * rot.thetay = 2.2;
     * rot.thetaz = 2.3;
     * //translation
     * insur::tk2CMSSW::Translation trans;
     * trans.dx = 5.0;
     * trans.dy = 6.0;
     * trans.dz = 7.0;
     * //posPart()
     * i1 = 0;
     * solid = "Parent";
     * material = "Child";
     * tc.posPart(i1, solid, material, rot, trans, stream);*/
    //materialSection
    el.tag = "Air";
    el.density = 123.4;
    el.atomic_number = 10;
    el.atomic_weight = 56.7;
    ev.push_back(el);
    p.first = el.tag;
    p.second = .99;
    co.elements.push_back(p);
    el.tag = "CF";
    el.density = 789.0;
    el.atomic_number = 88;
    el.atomic_weight = 98.7;
    ev.push_back(el);
    p.first = el.tag;
    p.second = .099;
    co.elements.push_back(p);
    el.tag = "Si";
    el.density = 333.3;
    el.atomic_number = 27;
    el.atomic_weight = 44.4;
    ev.push_back(el);
    p.first = el.tag;
    p.second = .001;
    co.elements.push_back(p);
    co.name = "Evilmix";
    co.density = 666.;
    co.method = insur::tk2CMSSW::wt;
    cv.push_back(co);
    co.name = "Quickmix";
    co.density = 55.55;
    co.method = insur::tk2CMSSW::vl;
    co.elements.clear();
    p.first = cv.at(0).elements.at(1).first;
    p.second = 0.55;
    co.elements.push_back(p);
    p.first = cv.at(0).elements.at(2).first;
    p.second = 0.45;
    co.elements.push_back(p);
    cv.push_back(co);
    co.name = "Mutantmix";
    co.density = 5432.1;
    co.method = insur::tk2CMSSW::ap;
    co.elements.clear();
    p.first = cv.at(0).elements.at(0).first;
    p.second = 1.;
    co.elements.push_back(p);
    cv.push_back(co);
    tc.materialSection("Materials", ev, cv, stream);
    //logicalPartSection
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
    lv.push_back(li);
    tc.logicalPartSection(lv, "Logic", stream);
    //solidSection
    si.name_tag = "Boxx";
    si.type = insur::tk2CMSSW::bx;
    si.dx = 15;
    si.dy = 20;
    si.dz = 25;
    sv.push_back(si);
    si.name_tag = "Trap";
    si.type = insur::tk2CMSSW::tp;
    si.dx = 12;
    si.dxx = 23;
    si.dy = 4.5;
    si.dz = 67;
    sv.push_back(si);
    si.name_tag = "Tube";
    si.type = insur::tk2CMSSW::tb;
    si.rmin = 30;
    si.rmax = 35;
    si.dz = 40;
    sv.push_back(si);
    tc.solidSection(sv, "Shapes", stream);
    std::cout << stream.str() << std::endl;
    return (EXIT_SUCCESS);
}

