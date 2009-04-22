// 
// File:   MatCalcDummy.h
// Author: ndemaio
//
// Created on April 21, 2009, 2:43 PM
//

#ifndef _MATCALCDUMMY_H
#define	_MATCALCDUMMY_H

#include <string>
#include <MatCalc.h>
namespace insur {
    class MatCalcDummy : public MatCalc {
    public:
        MatCalcDummy();
        virtual ~MatCalcDummy() {}
        virtual bool calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps);
        virtual bool calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps);
        virtual bool calculateBarrelServiceMaterials(
            std::vector<std::vector<ModuleCap> >& barrelcaps,
                std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
        virtual bool calculateEndcapServiceMaterials(
            std::vector<std::vector<ModuleCap> >& endcapcaps,
                std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
        virtual bool calculateSupportMaterials(std::vector<InactiveElement>& supports);
    private:
        std::string ta, ts, tl;
        double ab, ae, sb, se, lb, le;
        void init();
    };
}
#endif	/* _MATCALCDUMMY_H */

