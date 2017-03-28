/**
 * @file MaterialBudget.cc
 * @brief
 */

#include <MaterialBudget.hh>
namespace insur {
  /**
   * The constructor registers a tracker and a collection of inactive surfaces with the material budget before 
   * creating a module cap object for each active module found in the tracker.
   */
  MaterialBudget::MaterialBudget(Tracker& tr, InactiveSurfaces& is) {
    tracker = &tr;
    inactive = &is;

    class CapsVisitor : public GeometryVisitor {
      typedef std::vector<std::vector<ModuleCap>> Caps;
      Caps &capsbarrelmods_, &capsendmods_;
    public:
      CapsVisitor(Caps& capsbarrelmods, Caps& capsendmods) : capsbarrelmods_(capsbarrelmods), capsendmods_(capsendmods) {}
      void visit(Layer&) { capsbarrelmods_.push_back(std::vector<ModuleCap>()); }
      void visit(BarrelModule& m) {
		/*
        ModuleCap* cap = new ModuleCap(m);
        cap->setCategory(MaterialProperties::b_mod);
        capsbarrelmods_.back().push_back(*cap);
	m.setModuleCap(& (capsbarrelmods_.back().back()));
		*/
        capsbarrelmods_.back().push_back(*m.getModuleCap());
      }
      void visit(Disk&) { capsendmods_.push_back(std::vector<ModuleCap>()); }
      void visit(EndcapModule& m) {
		/*
        ModuleCap* cap = new ModuleCap(m);
        cap->setCategory(MaterialProperties::e_mod);
        capsendmods_.back().push_back(*cap);
        m.setModuleCap(& (capsendmods_.back().back()));
		*/
        capsendmods_.back().push_back(*m.getModuleCap());
      }
    };

    CapsVisitor v(capsbarrelmods, capsendmods);
    tr.accept(v);
  }

  /**
   * Nothing to do for the destructor...
   */
    MaterialBudget::~MaterialBudget() {}
    
    /**
     * Get the tracker object that belongs to this material budget.
     * @return A reference to the registered tracker object
     */
    Tracker& MaterialBudget::getTracker() { return *tracker; }
    
    /**
     * Get the collection of inactive surfaces that belongs to this material budget.
     * @return A reference to the registered inactive surfaces collection
     */
    InactiveSurfaces& MaterialBudget::getInactiveSurfaces() { return *inactive; }
    
    /**
     * Get the collection of barrel module caps.
     * @return A reference to the vector of vectors listing the module caps that are associated with the barrel modules
     */
    std::vector<std::vector<ModuleCap> >& MaterialBudget::getBarrelModuleCaps() { return capsbarrelmods; }
    
    /**
     * Get the collection of endcap module caps.
     * @return A reference to the vector of vectors listing the module caps that are associated with the endcap modules
     */
    std::vector<std::vector<ModuleCap> >& MaterialBudget::getEndcapModuleCaps() { return capsendmods; }
    
    /**
     * Print a summary of the material budget to <i>cout</i>.
     */
    void MaterialBudget::print() {
        std::cout << std::endl << "----------Material Budget Internal State----------" << std::endl;
        int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, o = 0, x = 0, s = 0, t = 0;
        for (unsigned int i = 0; i < inactive->getSupports().size(); i++) {
            if (inactive->getSupportPart(i).getCategory() == MaterialProperties::b_sup) a++;
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::e_sup) b++;
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::t_sup) g++;
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::u_sup) h++;
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::o_sup) o++;
            else x++;
        }
        for (unsigned int i = 0; i < capsbarrelmods.size(); i++) d = d + capsbarrelmods.at(i).size();
        for (unsigned int i = 0; i < capsendmods.size(); i++) f = f + capsendmods.at(i).size();
        std::cout << "MaterialBudget: " << d << " barrel modcaps in " << capsbarrelmods.size() << " vectors." << std::endl;
        std::cout << "MaterialBudget: " << f << " endcap modcaps in " << capsendmods.size() << " vectors." << std::endl;
        std::cout << "InactiveSurfaces: " << inactive->getBarrelServices().size() << " barrel services and ";
        std::cout << inactive->getEndcapServices().size() << " endcap services." << std::endl;
        std::cout << "InactiveSurfaces: " << inactive->getSupports().size() << " supports in total, "  << a << " in barrels, ";
        std::cout << b << " in endcaps, " << g << " in barrel support tubes, " << h << " user defined, " << o;
        std::cout << " in outer support tube and " << x << " unclassified.";
        std::cout << std::endl << std::endl;
        double am = 0, bm = 0, cm = 0, dm = 0, em = 0, fm = 0, gm = 0, hm = 0, lm = 0;
        double ar = 0, br = 0, cr = 0, dr = 0, er = 0, fr = 0, gr = 0, hr = 0, lr = 0;
        double ai = 0, bi = 0, ci = 0, di = 0, ei = 0, fi = 0, gi = 0, hi = 0, li = 0;
        for (unsigned int i = 0; i < capsbarrelmods.size(); i++) {
            for (unsigned int j = 0; j < capsbarrelmods.at(i).size(); j++) {
                if (capsbarrelmods.at(i).at(j).getTotalMass() > 0) {
                am = am + capsbarrelmods.at(i).at(j).getTotalMass();
                ar = ar + capsbarrelmods.at(i).at(j).getRadiationLength();
                ai = ai + capsbarrelmods.at(i).at(j).getInteractionLength();
                }
                else d--;
            }
        }
        std::cout << "-----Barrel Modules-----" << std::endl;
        std::cout << "Total valid modules in barrels: " << d << std::endl;
        std::cout << "Total module mass in barrels: " << am << std::endl;
        if (d > 0) std::cout << "Average module mass in barrels: " << (am / (double)d) << std::endl;
        if (d > 0) std::cout << "Average module radiation length in barrels: " << (ar / (double)d) << std::endl;
        if (d > 0) std::cout << "Average module interaction length in barrels: " << (ai / (double)d) << std::endl;
        for (unsigned int i = 0; i < capsbarrelmods.size(); i++) {
            std::cout << "Sample boundary module in layer " << (i + 1) << ":" << std::endl;
            try {
                capsbarrelmods.at(i).at(onBoundary(capsbarrelmods, i)).print();
            }
            catch (std::out_of_range& oor) {
                std::cerr << "Error: unable to find a module index for layer " << i << "." << std::endl;
                std::cerr << oor.what() << std::endl;
            }
        }
        std::cout << std::endl;
        for (unsigned int i = 0; i < capsendmods.size(); i++) {
            for (unsigned int j = 0; j < capsendmods.at(i).size(); j++) {
                if (capsendmods.at(i).at(j).getTotalMass() > 0) {
                bm = bm + capsendmods.at(i).at(j).getTotalMass();
                br = br + capsendmods.at(i).at(j).getRadiationLength();
                bi = bi + capsendmods.at(i).at(j).getInteractionLength();
                }
                else f--;
            }
        }
        std::cout << "-----Endcap Modules-----" << std::endl;
        std::cout << "Total valid modules in endcaps: " << f << std::endl;
        std::cout << "Total module mass in endcaps: " << bm << std::endl;
        if (f > 0) std::cout << "Average module mass in endcaps: " << (bm / (double)f) << std::endl;
        if (f > 0) std::cout << "Average module radiation length in endcaps: " << (br / (double)f) << std::endl;
        if (f > 0) std::cout << "Average module interaction length in endcaps: " << (bi / (double)f) << std::endl;
        for (unsigned int i = 0; i < capsendmods.size(); i++) {
            std::cout << "Sample boundary module on disc " << i << ":" << std::endl;
            try {
                capsendmods.at(i).at(onBoundary(capsendmods, i)).print();
            }
            catch (std::out_of_range& oor) {
                std::cerr << "Error: unable to find a module index for disc " << i << "." << std::endl;
                std::cerr << oor.what() << std::endl;
            }
        }
        std::cout << std::endl;
        for (unsigned int i = 0; i < inactive->getBarrelServices().size(); i++) {
            if (inactive->getBarrelServicePart(i).getTotalMass() > 0) {
            cm = cm + inactive->getBarrelServicePart(i).getTotalMass();
            cr = cr + inactive->getBarrelServicePart(i).getRadiationLength();
            ci = ci + inactive->getBarrelServicePart(i).getInteractionLength();
            }
            else s++;
        }
        std::cout << "-----Barrel Services-----" << std::endl;
        std::cout << "Total valid services in barrels: " << (inactive->getBarrelServices().size() - s) << std::endl;
        std::cout << "Total service mass in barrels: " << cm << std::endl;
        if (inactive->getBarrelServices().size() > 0)
            std::cout << "Average service mass in barrels: "<< (cm / (double)(inactive->getBarrelServices().size())) << std::endl;
        if (inactive->getBarrelServices().size() > 0)
            std::cout << "Average service radiation length in barrels: " << (cr / (double)(inactive->getBarrelServices().size())) << std::endl;
        if (inactive->getBarrelServices().size() > 0)
            std::cout << "Average service interaction length in barrels: " << (ci / (double)(inactive->getBarrelServices().size())) << std::endl;
        std::cout << std::endl;
        for (unsigned int i = 0; i < inactive->getEndcapServices().size(); i++) {
            if (inactive->getEndcapServicePart(i).getTotalMass() > 0) {
            dm = dm + inactive->getEndcapServicePart(i).getTotalMass();
            dr = dr + inactive->getEndcapServicePart(i).getRadiationLength();
            di = di + inactive->getEndcapServicePart(i).getInteractionLength();
            }
            else t++;
        }
        std::cout << "-----Endcap Services-----" << std::endl;
        std::cout << "Total valid services in endcaps: " << (inactive->getEndcapServices().size() - t) << std::endl;
        std::cout << "Total service mass in endcaps: " << dm << std::endl;
        if (inactive->getEndcapServices().size() > 0)
            std::cout << "Average service mass in endcaps: " << (dm / (double)(inactive->getEndcapServices().size())) << std::endl;
        if (inactive->getEndcapServices().size() > 0)
            std::cout << "Average service radiation length in endcaps: " << (dr / (double)(inactive->getEndcapServices().size())) << std::endl;
        if (inactive->getEndcapServices().size() > 0)
            std::cout << "Average service interaction length in endcaps: " << (di / (double)(inactive->getEndcapServices().size())) << std::endl;
        std::cout << std::endl;
        for (unsigned int i = 0; i < inactive->getSupports().size(); i++) {
            if (inactive->getSupportPart(i).getCategory() == MaterialProperties::b_sup) {
                if (inactive->getSupportPart(i).getTotalMass() > 0) {
                em = em + inactive->getSupportPart(i).getTotalMass();
                er = er + inactive->getSupportPart(i).getRadiationLength();
                ei = ei + inactive->getSupportPart(i).getInteractionLength();
                }
                else a--;
            }
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::e_sup){
                if (inactive->getSupportPart(i).getTotalMass() > 0) {
                fm = fm + inactive->getSupportPart(i).getTotalMass();
                fr = fr + inactive->getSupportPart(i).getRadiationLength();
                fi = fi + inactive->getSupportPart(i).getInteractionLength();
                }
                else b--;
            }
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::t_sup) {
                if (inactive->getSupportPart(i).getTotalMass() > 0) {
                gm = gm + inactive->getSupportPart(i).getTotalMass();
                gr = gr + inactive->getSupportPart(i).getRadiationLength();
                gi = gi + inactive->getSupportPart(i).getInteractionLength();
                }
                else g--;
            }
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::u_sup) {
                if (inactive->getSupportPart(i).getTotalMass() > 0) {
                hm = hm + inactive->getSupportPart(i).getTotalMass();
                hr = hr + inactive->getSupportPart(i).getRadiationLength();
                hi = hi + inactive->getSupportPart(i).getInteractionLength();
                }
                else h--;
            }
            else if (inactive->getSupportPart(i).getCategory() == MaterialProperties::o_sup) {
                if (inactive->getSupportPart(i).getTotalMass() > 0) {
                lm = lm + inactive->getSupportPart(i).getTotalMass();
                lr = lr + inactive->getSupportPart(i).getRadiationLength();
                li = li + inactive->getSupportPart(i).getInteractionLength();
                }
                else o--;
            }
        }
        std::cout << "-----Supports-----" << std::endl;
        std::cout << "Total valid support parts: " << a << " in barrels, " << b << " in endcaps, ";
        std::cout << g << " in barrel support tubes, " << o << " in outer support tube and " << h;
        std::cout << " user defined." << std::endl << std::endl;
        std::cout << "Total support mass in barrels: " << (em + gm) << std::endl;
        if ((a + g) > 0) std::cout << "Average support mass in barrels: " << ((em + gm)/ (double)(a + g)) << std::endl;
        if ((a + g) > 0) std::cout << "Average support radiation length in barrels: " << ((er + gr) / (double)(a + g)) << std::endl;
        if ((a + g) > 0) std::cout << "Average support interaction length in barrels: " << ((ei + gi) / (double)(a + g)) << std::endl;
        std::cout << std::endl;
        std::cout << "Total support mass in endcaps: " << fm << std::endl;
        if (b > 0) std::cout << "Average support mass in endcaps: " << (fm / (double)b) << std::endl;
        if (b > 0) std::cout << "Average support radiation length in endcaps: " << (fr / (double)b) << std::endl;
        if (b > 0) std::cout << "Average support interaction length in endcaps: " << (fi / (double)b) << std::endl;
        std::cout << std::endl;
        std::cout << "Total mass in outer support tube: " << lm << std::endl;
        if (o > 0) std::cout << "Average mass in outer support tube: " << (lm / (double)o) << std::endl;
        if (o > 0) std::cout << "Average radiation length in outer support tube: " << (lr / (double)o) << std::endl;
        if (o > 0) std::cout << "Average interaction length in outer support tube: " << (li / (double)o) << std::endl;
        std::cout << std::endl;
        std::cout << "Total mass in user defined supports: " << hm << std::endl;
        if (h > 0) std::cout << "Average mass in user defined supports: " << (hm / (double)h) << std::endl;
        if (h > 0) std::cout << "Average radiation length in user defined supports: " << (hr / (double)h) << std::endl;
        if (h > 0) std::cout << "Average interaction length in user defined supports: " << (hi / (double)h) << std::endl;
        std::cout << std::endl;
        std::cout << "Total mass in barrels: " << (am + cm + em + gm) << std::endl;
        std::cout << "Total mass in endcaps: " << (bm + dm + fm) << std::endl;
        std::cout << std::endl;
        double arb = 0, aib = 0, are = 0, aie = 0, cnt = 0;
        if (d > 0) {
            arb = arb + ar / (double)d;
            aib = aib + ai / (double)d;
            cnt++;
        }
        if (inactive->getBarrelServices().size() > 0) {
            arb = arb + cr / (double)(inactive->getBarrelServices().size());
            aib = aib + ci / (double)(inactive->getBarrelServices().size());
            cnt++;
        }
        if (a > 0) {
            arb = arb + er / (double)a;
            aib = aib + ei / (double)a;
            cnt++;
        }
        if (g > 0) {
            arb = arb + gr / (double)g;
            aib = aib + gi / (double)g;
            cnt++;
        }
        if (cnt > 0) {
            arb = arb / cnt;
            aib = aib / cnt;
        }
        cnt = 0;
        if (f > 0) {
            are = are + br / (double)f;
            aie = aie + bi / (double)f;
            cnt++;
        }
        if (inactive->getEndcapServices().size()) {
            are = are + dr / (double)(inactive->getEndcapServices().size());
            aie = aie + di / (double)(inactive->getEndcapServices().size());
            cnt++;
        }
        if (b > 0) {
            are = are + fr / (double) b;
            aie = aie + fi / (double)b;
            cnt++;
        }
        if (cnt > 0) {
            are = are / cnt;
            aie = aie / cnt;
        }
        std::cout << "Average radiation length in barrels: " << arb << std::endl;
        std::cout << "Average radiation length in endcaps: " << are << std::endl;
        std::cout << std::endl;
        std::cout << "Average interaction length in barrels: " << aib << std::endl;
        std::cout << "Average interaction length in endcaps: " << aie << std::endl;
        std::cout << std::endl;
        std::cout << "Total mass in tracker: " << ((am + bm + cm + dm + em + fm + gm +hm + lm) / 1000.0) << "kg." << std::endl;
        std::cout << "----------Material Budget Internal State----------" << std::endl << std::endl;
    }
    
    // protected
    /**
     * This finds a sample module in the given layer that sits at the end of its rod.
     * @param source The complete collection of barrel module caps
     * @param layer The layer within that the query applies to
     * @return The index of the sample module within the vector of module caps of the given layer
     */
    int MaterialBudget::onBoundary(std::vector<std::vector<ModuleCap> >& source, int layer) { //throws exception
        int ring = 0, index = 0;
        if ((layer >= 0) && (layer < (int)source.size())) {
            for (unsigned int mod = 0; mod < source.at(layer).size(); mod++) {
                int myring = source.at(layer).at(mod).getModule().as<BarrelModule>()->ring();
                if (myring > ring) {
                    ring = myring; 
                    index = mod;
                }
            }
        }
        else throw std::out_of_range("Layer index is out of range: " + layer);
        return index;
    }

  std::vector<InactiveElement> MaterialBudget::getAllServices() {
    std::vector<InactiveElement> allServices;
    auto& barrelServices = getInactiveSurfaces().getBarrelServices();
    auto& endcapServices = getInactiveSurfaces().getEndcapServices();
    auto& supports = getInactiveSurfaces().getSupports();

    // We put all services inside the same container
    allServices.reserve( barrelServices.size() + endcapServices.size() + supports.size() ); // preallocate memory
    allServices.insert( allServices.end(), barrelServices.begin(), barrelServices.end() );
    allServices.insert( allServices.end(), endcapServices.begin(), endcapServices.end() );
    allServices.insert( allServices.end(), supports.begin(), supports.end() );

    return allServices;
  }

}
