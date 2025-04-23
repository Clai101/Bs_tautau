#include <framework/core/Module.h>
#include <analysis/VariableManager/Manager.h>
#include <analysis/dataobjects/Particle.h>
#include <analysis/dataobjects/ParticleList.h>
#include <framework/datastore/StoreArray.h>
#include <framework/datastore/StoreObjPtr.h>
#include <mdst/dataobjects/MCParticle.h>
#include <analysis/utility/PCmsLabTransform.h>
#include <analysis/variables/PIDVariables.h>
#include <analysis/utility/ReferenceFrame.h>
#include <cmath>

#include <analysis/dataobjects/Particle.h>
#include <framework/gearbox/Const.h>
#include <framework/datastore/StoreObjPtr.h>
#include <framework/dataobjects/Helix.h>
#include <fstream>
// dataobjects from the MDST
#include <mdst/dataobjects/Track.h>
#include <mdst/dataobjects/MCParticle.h>
#include <mdst/dataobjects/TrackFitResult.h>
#include <mdst/dataobjects/EventLevelTrackingInfo.h>
#include <mdst/dataobjects/HitPatternVXD.h>
#include <mdst/dataobjects/ECLCluster.h>
#include <Math/Vector4D.h>
#include "TVectorF.h"


namespace Belle2 {
  namespace Variable {
    void getDaughterTracks(std::vector<const Particle*>& daughters,
                         const Particle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDGCode());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        daughters.push_back(particle);
      } else {
        n = particle->getNDaughters();
        for (i = 0; i < n; i++)
          getDaughterTracks(daughters, particle->getDaughter(i));
      }
    }

    void getDaughterGamma(std::vector<const Particle*>& daughters,
                         const Particle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDGCode());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        if (l == 22)
           daughters.push_back(particle);
      } else {
        n = particle->getNDaughters();
        for (i = 0; i < n; i++)
          getDaughterGamma(daughters, particle->getDaughter(i));
      }
    }
  int N_tracks_in_ROE(const Particle* part) {
      std::vector<const Particle*> daughters;
      getDaughterTracks(daughters, part);
      std::string listName = "pi+:alle";

      StoreObjPtr<ParticleList> listPi(listName);

      int N_tracks = 0;
      for (int j = 0; j < listPi->getListSize(); j++) {
        bool inROE = true;
        const Particle* part1 = listPi->getParticle(j);
        for (std::vector<const Particle*>::iterator it = daughters.begin(); it != daughters.end(); ++it) {
        if (part1->getMdstArrayIndex() == (*it)->getMdstArrayIndex()) {
          inROE = false;
          break;
        }
      }
      if (inROE) N_tracks+=1;
      }
  return N_tracks;
  }   

double pseudo_skim_y4s_MC_charm(const Particle*) {
        StoreObjPtr<EventMetaData> evtMetaData;
        int expNo = evtMetaData->getExperiment();
        int runNo = evtMetaData->getRun();
        int evtNo = evtMetaData->getEvent();
        std::ostringstream file_name;
        file_name << "/home/belle2/matrk/basf2/Labda_c_skim/"<<expNo<<"_"<<runNo<<".txt";
        std::ifstream file_in;
        file_in.open(file_name.str().c_str());
        std::string s;
        int index = 0;
        while (getline(file_in,s)) {
          if (evtNo == atof(s.c_str())) {
            index = 1;
            break;
          }
        }
        file_in.close();
        if (index == 1) return 1.0;
        return 0;
    }

 double pseudo_skim_data_hadron_y4s(const Particle*) {
        StoreObjPtr<EventMetaData> evtMetaData;
        int expNo = evtMetaData->getExperiment();
        int runNo = evtMetaData->getRun();
        int evtNo = evtMetaData->getEvent();
        std::ostringstream file_name;
        file_name << "/home/belle2/matrk/basf2/data_Lc_skim/"<<expNo<<"_"<<runNo<<".txt";
        std::ifstream file_in;
        file_in.open(file_name.str().c_str());
        std::string s;
        int index = 0;
        while (getline(file_in,s)) {
          if (evtNo == atof(s.c_str())) {
            index = 1;
            break;
          }
        }
        file_in.close();
        if (index == 1) return 1.0;
        return 0;
    }

     double pseudo_skim_data_exl_hadron_y4s(const Particle*) {
        StoreObjPtr<EventMetaData> evtMetaData;
        int expNo = evtMetaData->getExperiment();
        int runNo = evtMetaData->getRun();
        int evtNo = evtMetaData->getEvent();
        std::ostringstream file_name;
        file_name << "/home/belle2/matrk/basf2/data_exl_skim/"<<expNo<<"_"<<runNo<<".txt";
        std::ifstream file_in;
        file_in.open(file_name.str().c_str());
        std::string s;
        int index = 0;
        while (getline(file_in,s)) {
          if (evtNo == atof(s.c_str())) {
            index = 1;
            break;
          }
        }
        file_in.close();
        if (index == 1) return 1.0;
        return 0;
    }

    double pseudo_skim_data_exl_all_events(const Particle*) {
      StoreObjPtr<EventMetaData> evtMetaData;
      int expNo = evtMetaData->getExperiment();
      int runNo = evtMetaData->getRun();
      int evtNo = evtMetaData->getEvent();
      std::ostringstream file_name;
      file_name << "/home/belle2/matrk/basf2/all_events_Xc_skim/"<<expNo<<"_"<<runNo<<".txt";
      std::ifstream file_in;
      file_in.open(file_name.str().c_str());
      std::string s;
      int index = 0;
      while (getline(file_in,s)) {
        if (evtNo == atof(s.c_str())) {
          index = 1;
          break;
        }
      }
      file_in.close();
      if (index == 1) return 1.0;
      return 0;
  }

  void getGammaInROE(std::vector<const Particle*>& daughters, const Particle* part) {
    std::string listName = "gamma:alle";
    StoreObjPtr<ParticleList> listGamma(listName);
    std::vector<const Particle*> daughters1;
    getDaughterGamma(daughters1, part);
    for (unsigned int j = 0; j < listGamma->getListSize(); j++) {
      bool inROE = true;
      const Particle* part1 = listGamma->getParticle(j);
      for (std::vector<const Particle*>::iterator it = daughters1.begin(); it != daughters1.end(); ++it) {
      if (part1->getMdstArrayIndex() == (*it)->getMdstArrayIndex()) {
        inROE = false;
        break;
      }
    }
    if (inROE) daughters.push_back(part1);
    }
  }

  double E_gamma_in_ROE(const Particle* part) {
    std::vector<const Particle*> daughters;
    getGammaInROE(daughters, part);
    ROOT::Math::PxPyPzEVector p4;
    if (daughters.size() != 0) {
       int N_photons = daughters.size();
       for (int i = 0; i < N_photons; i++)
               p4 += daughters[i]->get4Vector();
       return p4.E();
    }
    else return 0;
  }
  
  double mcPx(const Particle* part)
    {
      StoreArray<MCParticle> mcparticles;
      if (mcparticles.getEntries() < 1)
        return -999;

      // MC truth momentum
      const MCParticle* mcpart = part->getRelatedTo<MCParticle>();
      if (mcpart == nullptr)
        return -999;
      ROOT::Math::PxPyPzEVector pMomentum = mcpart->get4Vector();

      // MC truth boost
      ROOT::Math::PxPyPzEVector cmsMomentum = mcparticles[0]->get4Vector();
      B2Vector3D cmsBoost    = (cmsMomentum.BoostToCM());
      pMomentum=ROOT::Math::Boost(cmsBoost)*pMomentum;

      B2Vector3D pMomentum3D = pMomentum.Vect();
      double p = pMomentum3D.Px();
      return p;
    }

  double mcPy(const Particle* part)
    {
      StoreArray<MCParticle> mcparticles;
      if (mcparticles.getEntries() < 1)
        return -999;

      // MC truth momentum
      const MCParticle* mcpart = part->getRelatedTo<MCParticle>();
      if (mcpart == nullptr)
        return -999;
      ROOT::Math::PxPyPzEVector pMomentum = mcpart->get4Vector();

      // MC truth boost
      ROOT::Math::PxPyPzEVector cmsMomentum = mcparticles[0]->get4Vector();
      B2Vector3D cmsBoost    = (cmsMomentum.BoostToCM());
      pMomentum=ROOT::Math::Boost(cmsBoost)*pMomentum;

      B2Vector3D pMomentum3D = pMomentum.Vect();
      double p = pMomentum3D.Py();
      return p;
    }

  double mcPz(const Particle* part)
    {
      StoreArray<MCParticle> mcparticles;
      if (mcparticles.getEntries() < 1)
        return -999;

      // MC truth momentum
      const MCParticle* mcpart = part->getRelatedTo<MCParticle>();
      if (mcpart == nullptr)
        return -999;
      ROOT::Math::PxPyPzEVector pMomentum = mcpart->get4Vector();

      // MC truth boost
      ROOT::Math::PxPyPzEVector cmsMomentum = mcparticles[0]->get4Vector();
      B2Vector3D cmsBoost    = (cmsMomentum.BoostToCM());
      pMomentum=ROOT::Math::Boost(cmsBoost)*pMomentum;

      B2Vector3D pMomentum3D = pMomentum.Vect();
      double p = pMomentum3D.Pz();
      return p;
    }

    B2Vector3D vector_prod(B2Vector3D p1, B2Vector3D p2)
    {
      return B2Vector3D(-p1.Pz()*p2.Py() + p1.Py()*p2.Pz(), -p1.Px()*p2.Pz() + p1.Pz()*p2.Px(), -p1.Py()*p2.Px() + p1.Px()*p2.Py());
    }


  double cos_theta_L(const Particle* particle) 
  {
    const Particle* part = particle;

    PCmsLabTransform T;

    ROOT::Math::PxPyPzEVector pSum = T.rotateLabToCms()*part->get4Vector();
    ROOT::Math::PxPyPzEVector pdaughter1 = T.rotateLabToCms()*part->getDaughter(0)->get4Vector();
    ROOT::Math::PxPyPzEVector pbeam(0, 0, 0, T.getCMSEnergy());

    B2Vector3D boost_vect = pSum.BoostToCM();

    pbeam = ROOT::Math::Boost(boost_vect)*pbeam;
    pdaughter1 = ROOT::Math::Boost(boost_vect)*pdaughter1;

    B2Vector3D p1 = -pbeam.Vect();
    B2Vector3D p2 = pdaughter1.Vect();

    return cos(p1.Angle(p2));
  }
  

double cos_theta_p(const Particle* particle) 
  {
    const Particle* part = particle;
    PCmsLabTransform T;
    ROOT::Math::PxPyPzEVector pbeam = T.rotateLabToCms()*part->get4Vector();
    const Particle* parts = part->getDaughter(0);

    ROOT::Math::PxPyPzEVector pSum = T.rotateLabToCms()*parts->get4Vector();
    ROOT::Math::PxPyPzEVector pdaughter1 = T.rotateLabToCms()*parts->getDaughter(0)->get4Vector();

    B2Vector3D boost_vect = pSum.BoostToCM();

    pbeam = ROOT::Math::Boost(boost_vect)*pbeam;
    pdaughter1 = ROOT::Math::Boost(boost_vect)*pdaughter1;

    B2Vector3D p1 = -pbeam.Vect();
    B2Vector3D p2 = pdaughter1.Vect();

    return cos(p1.Angle(p2));
  }

  double det(B2Vector3D v1, B2Vector3D v2, B2Vector3D v3) 
  {
    return v1.Px()*v2.Py()*v3.Pz() - v1.Px()*v2.Pz()*v3.Py() + v1.Pz()*v2.Px()*v3.Py() - v1.Pz()*v2.Py()*v3.Px() + v1.Py()*v2.Pz()*v3.Px() - v1.Py()*v2.Px()*v3.Pz();
  }


template <typename T> 
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}



double Delta(const Particle* particle) 
  {
    PCmsLabTransform T;
    const Particle* Lc = particle;
    const Particle* L = Lc->getDaughter(0);
    const Particle* prot = L->getDaughter(0);
    const Particle* pi = L->getDaughter(1);
    ROOT::Math::PxPyPzEVector pbeam4D(0, 0, 0, T.getCMSEnergy());

    ROOT::Math::PxPyPzEVector pSum = T.rotateLabToCms()*Lc->get4Vector();
    B2Vector3D boost_vect = pSum.BoostToCM();

    ROOT::Math::PxPyPzEVector pbeam4d = pbeam4D;
    B2Vector3D pbeam = (ROOT::Math::Boost(boost_vect)*pbeam4d).Vect();
    ROOT::Math::PxPyPzEVector pprot4d = T.rotateLabToCms()*prot->get4Vector();
    B2Vector3D pprot = (ROOT::Math::Boost(boost_vect)*pprot4d).Vect();
    ROOT::Math::PxPyPzEVector ppi4d = T.rotateLabToCms()*pi->get4Vector();
    B2Vector3D ppi = (ROOT::Math::Boost(boost_vect)*ppi4d).Vect();
    ROOT::Math::PxPyPzEVector pL4d = T.rotateLabToCms()*L->get4Vector();
    B2Vector3D pL = (ROOT::Math::Boost(boost_vect)*pL4d).Vect();

    B2Vector3D v1 = vector_prod(pprot, ppi);
    B2Vector3D v2 = vector_prod(pL, pbeam);

    return v1.Angle(v2)*sign(det(v1, v2, pL));
  }


    VARIABLE_GROUP("CUSTOM_VARIABLES");
    REGISTER_VARIABLE("pseudo_skim_y4s_MC_charm",pseudo_skim_y4s_MC_charm,"");
    REGISTER_VARIABLE("pseudo_skim_data_hadron_y4s",pseudo_skim_data_hadron_y4s,"");
    REGISTER_VARIABLE("pseudo_skim_data_exl_hadron_y4s",pseudo_skim_data_exl_hadron_y4s,"");
    REGISTER_VARIABLE("pseudo_skim_data_exl_all_events",pseudo_skim_data_exl_all_events,"");
    REGISTER_VARIABLE("mcPx",mcPx,"");
    REGISTER_VARIABLE("mcPy",mcPy,"");
    REGISTER_VARIABLE("mcPz",mcPz,"");
    REGISTER_VARIABLE("cos_theta_L",cos_theta_L,"");
    REGISTER_VARIABLE("cos_theta_p",cos_theta_p,"");
    REGISTER_VARIABLE("Delta",Delta,"");
    REGISTER_VARIABLE("N_tracks_in_ROE",N_tracks_in_ROE,"");
    REGISTER_VARIABLE("E_gamma_in_ROE",E_gamma_in_ROE,"");
  }
  // Create an empty module which does nothing at all. What it does is allowing
  // basf2 to easily find the library and load it from the steering file
  class EnableMyVariableModule: public Module {}; // And register this module to create a .map lookup file.
  REG_MODULE(EnableMyVariable);
}
