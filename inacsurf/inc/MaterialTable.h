//
// File:   MaterialTable.h
// Author: ndemaio
//
// Created on October 14, 2008, 10:53 AM
//

/**
 * @file MaterialTable.h
 * @brief This is the header file for the internal data structure for material properties.
 */

#ifndef _MATERIALTABLE_H
#define	_MATERIALTABLE_H

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
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
        uint rowCount();
        bool empty();
        void print();
    protected:
        std::vector<MaterialRow> materials;
    private:
        int findIndex(std::string tag);
    };
}
#endif	/* _MATERIALTABLE_H */

