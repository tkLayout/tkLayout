//
// File:   MaterialTable.h
// Author: ndemaio
//
// Created on October 14, 2008, 10:53 AM
// MaterialTable is nearly not used and was replaced by MaterialTab. This should be reworked.

/**
 * @file MaterialTable.h
 * @brief This is the header file for the internal data structure for material properties.
 */

#ifndef _MATERIALTABLE_H
#define	_MATERIALTABLE_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <math.h>

#include <global_funcs.hh>

namespace insur {
    /**
     * Error messages that may be included in exceptions
     */
    static const std::string errEntryNotFound = "Entry does not exist.";
    static const std::string errIndexOutOfRange = "Index out of range.";
    
    /**
     * @struct MaterialRow
     * @brief One instance of the MaterialRow struct contains the relevant information about a single basic, i.e. non-mixture, material.
     * @param tag The name of the basic material; used in search functions
     * @param density The density of the basic material
     * @param rlength The radiation length of the basic material
     * @param ilength The interaction length of the basic material
     */
    struct MaterialRow {
        std::string tag;
        double density;
        double rlength;
        double ilength;
    };



    /**
     * @class MaterialTable
     * @brief Essentially, this is a collection class for <i>MaterialRow</i> instances.
     *
     * It provides access functions for individual entries and a couple of bookkeeping calls.
     * The access functions may trow exceptions if an entry doesn't exist or if an index is out of range.
     */
    class MaterialTable {
    public:
        MaterialTable() {}
        virtual ~MaterialTable() {}
        void addMaterial(MaterialRow mat);
        void addMaterial(std::string tag, double density, double rlength, double ilength);
        MaterialRow& getMaterial(std::string tag); // throws exception
        MaterialRow& getMaterial(int index); // throws exception
        bool replaceMaterial(std::string oldtag, MaterialRow newmat);
        bool replaceMaterial(int index, MaterialRow newmat);
        unsigned int rowCount();
        bool empty();
        void print();
    protected:
        std::vector<MaterialRow> materials;
    private:
        int findIndex(std::string tag);
    };



    class BaseMaterial {  // this wasn't made abstract because it's used as key class in a set, therefore it is necessary to create an instance of it to perform key-based lookups
    protected:
      std::string tag_;
      std::set<std::string> aliases_;
    public:
      virtual ~BaseMaterial() {}
      std::string getTag() const { return tag_; }
      const std::set<std::string>& getAliases() const { return aliases_; }
      virtual double getDensity() = 0;
      virtual double getRadiationLength() = 0;
      virtual double getInteractionLength() = 0;
      virtual std::string toXML() = 0;
      virtual void fromCfg(std::string) = 0;
    };

    // TODO: restart the development of this "new" material table so that supports composite materials and syncs with CMSSW
    class MaterialTable2 {
      std::map<std::string, BaseMaterial*> materials_;
    public:
      void parseMaterial(std::string str);
      BaseMaterial* getMaterial(std::string tags) const { return materials_.at(tags); }
      const std::map<std::string, BaseMaterial*>& getMaterials() const { return materials_; }
      ~MaterialTable2() { for (std::map<std::string, BaseMaterial*>::iterator it = materials_.begin(); it != materials_.end(); ++it) delete it->second; }
    }; 




    class ElementaryMaterial : public BaseMaterial {
      double density_;
      double z_, a_;
      double rlength_, ilength_;
    public:
      ElementaryMaterial() : density_(0.), z_(0.), a_(0.), rlength_(0.), ilength_(0.) {}
      double getDensity() { return density_; }
      double getAtomicNumber() { return z_; } 
      double getAtomicWeight() { return a_; } 
      double getRadiationLength();
      double getInteractionLength();
      std::string toXML() { return std::string(); }
      void fromCfg(std::string);
    };

    class CompositeMaterial : public BaseMaterial {
      const MaterialTable2* const mTable_;
      struct MaterialFraction { BaseMaterial* material; double fraction; };
      std::vector<MaterialFraction> materialFractions_;
      double density_;
      double rlength_, ilength_;
    public:
      CompositeMaterial(const MaterialTable2* const mTable) : mTable_(mTable), density_(0.), rlength_(0.), ilength_(0.) {} 
      double getDensity() { return density_; }
      double getRadiationLength();
      double getInteractionLength();
      std::string toXML() { return std::string(); }
      void fromCfg(std::string);
    };
}
#endif	/* _MATERIALTABLE_H */

