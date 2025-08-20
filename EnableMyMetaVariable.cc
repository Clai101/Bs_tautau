#include <framework/core/Module.h>
#include <analysis/VariableManager/Manager.h>
#include <analysis/dataobjects/Particle.h>
#include <analysis/dataobjects/ParticleList.h>
#include <framework/datastore/StoreArray.h>
#include <framework/datastore/StoreObjPtr.h>
#include <mdst/dataobjects/MCParticle.h>
#include <analysis/utility/PCmsLabTransform.h>
#include <framework/utilities/Conversion.h>
#include <analysis/utility/ReferenceFrame.h>
#include <Math/Vector4D.h>
#include <cmath>
#include <Math/Boost.h>
#include <fstream>
#include "TVectorF.h"

namespace Belle2 {
  namespace Variable {

      
    void getDaughterTracks(std::vector<const Particle*>& daughters,
                         const Particle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDGCode());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321) daughters.push_back(particle);
      } else {
        n = particle->getNDaughters();
        for (i = 0; i < n; i++)
          getDaughterTracks(daughters, particle->getDaughter(i));
      }
    }
   
    void getDaughterMCTracks(std::vector<const MCParticle*>& daughters,
                         const MCParticle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDG());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321) daughters.push_back(particle);
      } else 
          for (auto daug : particle->getDaughters())
		  getDaughterMCTracks(daughters, daug);
      
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

    void getMCDaughterGamma(std::vector<const MCParticle*>& daughters,
                         const MCParticle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDG());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        if (l == 22)
           daughters.push_back(particle);
      } else
           for (auto daug : particle->getDaughters())
             getMCDaughterGamma(daughters, daug);
    }

    void getDaughterMCGamma(std::vector<const MCParticle*>& daughters,
                         const MCParticle* particle)
    {
      int i, n, l;
      l = abs(particle->getPDG());
      if (l == 11 || l == 13 || l == 211 || l == 2212 || l == 321 || l == 22) {
        if (l == 22) daughters.push_back(particle);
      } else 
          for (auto daug : particle->getDaughters())
                  getDaughterMCTracks(daughters, daug);

    }


    Manager::FunctionPtr varValue(const std::vector<std::string>& arguments)
    {
      if (arguments.size() != 3) {
        B2FATAL("Wrong number of arguments for varValue");
      }

      auto listName = arguments[0];
      const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[1]);
      int index = Belle2::convertString<int>(arguments[2]);

      auto func = [listName, var, index](const Particle*)-> double {

        StoreObjPtr<ParticleList> list(listName);

        if (!(list.isValid())) {
          B2FATAL("Invalid list name " << listName << " is given to varValue");
        }

        if (list->getListSize() == 0) return -999.;
        if (index >= int(list->getListSize())) {
                return -999;
              }

        const Particle* particle = list->getParticle(index);
        return std::get<double>(var->function(particle));
      };
      return func;
    }

   Manager::FunctionPtr InROE(const std::vector<std::string>& arguments)
   {
     if (arguments.size() != 2) {
        std::cout << arguments.size() << "\n";
        //return -999.;
        B2FATAL("Wrong number of arguments for inROE");
      }
      auto listName = arguments[0];
      int index = Belle2::convertString<int>(arguments[1]);

      auto func = [listName, index](const Particle* part)-> double {

        std::vector<const Particle*> daughter_tracks;
        std::vector<const Particle*> daughter_gamma;
        unsigned int particleArrayIndex;
        bool addParticle;

        addParticle = true;
        StoreObjPtr<ParticleList> list(listName);

        if (!(list.isValid())) {
          B2FATAL("Invalid list name " << listName << " is given to varValue");
        }

        if (list->getListSize() == 0) return 1.;
        if (index >= int(list->getListSize())) {
          return -999;
        }
        
        getDaughterTracks(daughter_tracks, list->getParticle(index));
        getDaughterGamma(daughter_gamma, list->getParticle(index));

        if (part->getPDGCode() == 22)
          for (auto* daughter : daughter_gamma) {
            particleArrayIndex = daughter->getMdstArrayIndex();
            if (part->getMdstArrayIndex() == particleArrayIndex) {
              addParticle = false;
              break;
            }
          }
        else
          for (auto* daughter : daughter_tracks) {
            particleArrayIndex = daughter->getMdstArrayIndex();
            if (part->getMdstArrayIndex() == particleArrayIndex) {
              addParticle = false;
              break;
            }
          }

        if (addParticle)  return 1.0;
        else return 0.0;
      };
      return func;
    }


   Manager::FunctionPtr InROEMC(const std::vector<std::string>& arguments)
   {
     if (arguments.size() != 2) {
        std::cout << arguments.size() << "\n";
        //return -999.;
        B2FATAL("Wrong number of arguments for inROEMC");
      }
      auto listName = arguments[0];
      int index = Belle2::convertString<int>(arguments[1]);

      auto func = [listName, index](const Particle* part)-> double {

        std::vector<const MCParticle*> daughter_tracks;
        std::vector<const MCParticle*> daughter_gamma;
        bool addParticle;
	const MCParticle* mcpart = part->getMCParticle();
        
	if (!mcpart)
           addParticle = true;
        else {
        StoreObjPtr<ParticleList> list(listName);

        if (!(list.isValid())) {
          B2FATAL("Invalid list name " << listName << " is given to varValue");
        }

        if (list->getListSize() == 0) return 1.;
        if (index >= int(list->getListSize())) {
          return -999;
        }


        getDaughterMCTracks(daughter_tracks, list->getParticle(index)->getMCParticle());
        getDaughterMCGamma(daughter_gamma, list->getParticle(index)->getMCParticle());

        if (part->getPDGCode() == 22)
          for (auto* daughter : daughter_gamma) {
            if (mcpart->getIndex() == daughter->getIndex()) {
              addParticle = false;
              break;
            }
          }
        else
          for (auto* daughter : daughter_tracks) {
            if (mcpart->getIndex() == daughter->getIndex()) {
              addParticle = false;
              break;
            }
          }
        }
        if (addParticle)  return 1.0;
        else return 0.0;
      };
      return func;
    }



    Manager::FunctionPtr FSPProductOf(const std::vector<std::string>& arguments)
    {
      if (arguments.size() == 1) {
        const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[0]);
        auto func = [var](const Particle * particle) -> double {
          std::vector<const Particle*> daughters;
	  getDaughterTracks(daughters, particle);
	  double product = 1.0;
          for (std::vector<const Particle*>::iterator it = daughters.begin(); it != daughters.end(); ++it)
	  {
            product *= std::get<double>(var->function(*it));
          }
          return product;
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function FSPProductOf");
      }
    }

    Manager::FunctionPtr FSPSumOf(const std::vector<std::string>& arguments)
    {
      if (arguments.size() == 1) {
        const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[0]);
        auto func = [var](const Particle * particle) -> double {
          std::vector<const Particle*> daughters;
          getDaughterTracks(daughters, particle);
          double sum = 0.0;
          for (std::vector<const Particle*>::iterator it = daughters.begin(); it != daughters.end(); ++it)
          { 
            sum += std::get<double>(var->function(*it));
          }
          return sum;
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function FSPProductOf");
      }
    }


    Manager::FunctionPtr daughterHelAngle(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 1) {
        std::vector<int> daughterIndices;
        try { 
          for (auto& argument : arguments) daughterIndices.push_back(Belle2::convertString<int>(argument));
        } catch (boost::bad_lexical_cast&) {
          B2WARNING("The arguments of daughterHelAngle meta function must be integers!");
          return nullptr;
        }
        auto func = [daughterIndices](const Particle * particle) -> double {
          if (particle == nullptr)
            return -999;
          else {
            const auto& frame = ReferenceFrame::GetCurrent();
            ROOT::Math::PxPyPzEVector pSum;

            for (auto& index : daughterIndices)
            {
              if (index >= int(particle->getNDaughters())) {
                return -999; 
              } else pSum += frame.getMomentum(particle->getDaughter(index));
            }

            B2Vector3D daughtersBoost = pSum.BoostToCM();

            ROOT::Math::PxPyPzEVector pMother;
            ROOT::Math::PxPyPzEVector pDaughter;
            pMother = frame.getMomentum(particle);

            pMother = ROOT::Math::Boost(daughtersBoost)*pMother;
            pDaughter = frame.getMomentum(particle->getDaughter(daughterIndices[0]));
            pDaughter = ROOT::Math::Boost(daughtersBoost)*pDaughter;

            B2Vector3D p1(pMother.Px(), pMother.Py(), pMother.Pz());
	    B2Vector3D p2(pDaughter.Px(), pDaughter.Py(), pDaughter.Pz());
            //ROOT::Math::XYZVector p2 = pDaughter.Vect();

            return cos(p1.Angle(p2));
            }
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function daughterHelAngle. At least two integers are needed.");
      }
    }


    Manager::FunctionPtr minimumValueInList(const std::vector<std::string>& arguments)
    {
      if (arguments.size() == 2) {
        std::string listName = arguments[0];
        const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[1]);

        auto func = [listName, var](const Particle*) -> double {
          StoreObjPtr<ParticleList> listOfParticles(listName);

          if (!(listOfParticles.isValid())) B2FATAL("Invalid list name " << listName << " given to averageValueInList");
          int nParticles = listOfParticles->getListSize();
          if (nParticles == 0)
             return 999;
          else {
          double minimum = std::get<double>(var->function(listOfParticles->getParticle(0))); 
          for (int i = 0; i < nParticles; i++)
          {
            const Particle* part = listOfParticles->getParticle(i);
            if (minimum > std::get<double>(var->function(part)))
               minimum = std::get<double>(var->function(part));
          }
          return minimum;
          }
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function minimumValueInList");
      }
    }

       Manager::FunctionPtr maximumValueInList(const std::vector<std::string>& arguments)
    {
      if (arguments.size() == 2) {
        std::string listName = arguments[0];
        const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[1]);

        auto func = [listName, var](const Particle*) -> double {
          StoreObjPtr<ParticleList> listOfParticles(listName);

          if (!(listOfParticles.isValid())) B2FATAL("Invalid list name " << listName << " given to averageValueInList");
          int nParticles = listOfParticles->getListSize();
          if (nParticles == 0)
             return -999;
          else {
          double maximum = std::get<double>(var->function(listOfParticles->getParticle(0)));
          for (int i = 0; i < nParticles; i++)
          {
            if (maximum < std::get<double>(var->function(listOfParticles->getParticle(i))))
                maximum = std::get<double>(var->function(listOfParticles->getParticle(i)));
          }
          return maximum;
          }
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function maximumValueInList");
      }
    }

    Manager::FunctionPtr ChargeInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          double chrg = 0;

          std::vector<Particle*> particlePool;

          (void) particle;
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to ChrgInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                chrg += part->getCharge();
                particlePool.push_back(part);
              }
            }
          }

          return chrg;

        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function ChrgInLists");
      }
    }

    Manager::FunctionPtr PxInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          ROOT::Math::PxPyPzEVector total4VectorCMS;

          std::vector<Particle*> particlePool;

          (void) particle; 
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to PxInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                total4VectorCMS += PCmsLabTransform::labToCms(part->get4Vector());
                particlePool.push_back(part);
              }
            }
          }

          B2Vector3D total3VectorCMS = total4VectorCMS.Vect();
          double Px = total3VectorCMS.Px();
          return Px;
        
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function PxInLists");
      }
    }

    Manager::FunctionPtr PyInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          ROOT::Math::PxPyPzEVector total4VectorCMS;

          std::vector<Particle*> particlePool;

          (void) particle; 
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to PyInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                total4VectorCMS += PCmsLabTransform::labToCms(part->get4Vector());
                particlePool.push_back(part);
              }
            }
          }

	  B2Vector3D total3VectorCMS = total4VectorCMS.Vect();
          double Py = total3VectorCMS.Py();
          return Py;
        
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function PyInLists");
      }
    }

    Manager::FunctionPtr PzInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          ROOT::Math::PxPyPzEVector total4VectorCMS;

          std::vector<Particle*> particlePool;

          (void) particle;
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to PzInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                total4VectorCMS += PCmsLabTransform::labToCms(part->get4Vector());
                particlePool.push_back(part);
              }
            }
          }

	  B2Vector3D total3VectorCMS = total4VectorCMS.Vect();
          double Pz = total3VectorCMS.Pz();
          return Pz;

        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function PzInLists");
      }
    }

    Manager::FunctionPtr EInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          ROOT::Math::PxPyPzEVector total4VectorCMS;

          std::vector<Particle*> particlePool;

          (void) particle;
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to EInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                total4VectorCMS += PCmsLabTransform::labToCms(part->get4Vector());
                particlePool.push_back(part);
              }
            }
          }


          double E = total4VectorCMS.E();
          return E;

        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function EInLists");
      }
    }

    Manager::FunctionPtr ElabInLists(const std::vector<std::string>& arguments)
    {
      if (arguments.size() > 0) {

        auto func = [arguments](const Particle * particle) -> double {

          ROOT::Math::PxPyPzEVector total4VectorCMS;

          std::vector<Particle*> particlePool;

          (void) particle;
          for (const auto& argument : arguments)
          {
            StoreObjPtr <ParticleList> listOfParticles(argument);

            if (!(listOfParticles.isValid())) B2FATAL("Invalid Listname " << argument << " given to EInLists");
            int nParticles = listOfParticles->getListSize();
            for (int i = 0; i < nParticles; i++) {
              bool overlaps = false;
              Particle* part = listOfParticles->getParticle(i);
              for (auto poolPart : particlePool) {
                if (part->overlapsWith(poolPart)) {
                  overlaps = true;
                  break;
                }
              }
              if (!overlaps) {
                total4VectorCMS += (part->get4Vector());
                particlePool.push_back(part);
              }
            }
          }


          double E = total4VectorCMS.E();
          return E;

        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function EInLists");
      }
    }

  void getGammaInROE(std::vector<const Particle*>& daughters, const Particle* part) {
      std::vector<const Particle*> daughters1;
      getDaughterGamma(daughters1, part);
      std::string listName = "gamma:good";
      StoreObjPtr<ParticleList> listGamma(listName);

      for (int j = 0; j < listGamma->getListSize(); j++) {
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

  void getGammaInMCROE(std::vector<const Particle*>& daughters, const Particle* part) {
      std::string listName = "gamma:good";
      StoreObjPtr<ParticleList> listGamma(listName);
      std::vector<const MCParticle*> daughters1;
      getMCDaughterGamma(daughters1, part->getMCParticle());
      for (unsigned int j = 0; j < listGamma->getListSize(); j++) {
        bool inROE = true;
        const Particle* part1 = listGamma->getParticle(j);
        if (part1->getMCParticle())
        for (std::vector<const MCParticle*>::iterator it = daughters1.begin(); it != daughters1.end(); ++it) {
          if (part1->getMCParticle()->getArrayIndex() == (*it)->getArrayIndex()) {
            inROE = false;
            break;
          }
        }
      if (inROE) daughters.push_back(part1);
      }
  }

   double HardestGammaInROEMomentum(const Particle* part, const std::vector<double>& index) {
    if (index.size() != 1) {
        B2FATAL("Wrong number of arguments for cosHelicityAngleIfCMSIsTheMother");
      }
      int my_arg = std::lround(index[0]);

      std::vector<const Particle*> daughters;
      getGammaInROE(daughters, part);
      if (daughters.size() == 0) return -999;
      else {
          ROOT::Math::PxPyPzEVector pSum;
          ROOT::Math::PxPyPzEVector p4 = daughters[0] -> get4Vector();

          if (my_arg == 0) return p4.Px();
          else if (my_arg == 1) return p4.Py();
          else if (my_arg == 2) return p4.Pz();
          else return -999;
     }
}

double HardestGammaInROECMSMomentum(const Particle* part, const std::vector<double>& index) {
    if (index.size() != 1) {
        B2FATAL("Wrong number of arguments for cosHelicityAngleIfCMSIsTheMother");
      }
      int my_arg = std::lround(index[0]);

      std::vector<const Particle*> daughters;
      getGammaInROE(daughters, part);
      if (daughters.size() == 0) return -999;
      else {
          ROOT::Math::PxPyPzEVector pSum;
          PCmsLabTransform T;
          ROOT::Math::PxPyPzEVector p4 = daughters[0] -> get4Vector();
          ROOT::Math::PxPyPzEVector p4CMS = T.rotateLabToCms() * p4;

          if (my_arg == 0) return p4CMS.Px();
          else if (my_arg == 1) return p4CMS.Py();
          else if (my_arg == 2) return p4CMS.Pz();
          else return -999;
     }
}

   double HardestGammaInMCROEMomentum(const Particle* part, const std::vector<double>& index) {
    if (index.size() != 1) {
        B2FATAL("Wrong number of arguments for cosHelicityAngleIfCMSIsTheMother");
      }
      int my_arg = std::lround(index[0]);

      std::vector<const Particle*> daughters;
      getGammaInMCROE(daughters, part);
      if (daughters.size() == 0) return -999;
      else {
          ROOT::Math::PxPyPzEVector pSum;
          ROOT::Math::PxPyPzEVector p4 = daughters[0] -> get4Vector();

          if (my_arg == 0) return p4.Px();
          else if (my_arg == 1) return p4.Py();
          else if (my_arg == 2) return p4.Pz();
          else return -999;
     }
}

   double HardestGammaInMCROECMSMomentum(const Particle* part, const std::vector<double>& index) {
    if (index.size() != 1) {
        B2FATAL("Wrong number of arguments for cosHelicityAngleIfCMSIsTheMother");
      }
      int my_arg = std::lround(index[0]);

      std::vector<const Particle*> daughters;
      getGammaInMCROE(daughters, part);
      if (daughters.size() == 0) return -999;
      else {
          ROOT::Math::PxPyPzEVector pSum;
          PCmsLabTransform T;
          ROOT::Math::PxPyPzEVector p4 = daughters[0] -> get4Vector();
          ROOT::Math::PxPyPzEVector p4CMS = T.rotateLabToCms() * p4;

          if (my_arg == 0) return p4CMS.Px();
          else if (my_arg == 1) return p4CMS.Py();
          else if (my_arg == 2) return p4CMS.Pz();
          else return -999;
     }
}


 double pseudo_skim_y5s_MC(const Particle* part, const std::vector<std::string>& arguments) {
        if (arguments.size() != 2) {
        B2FATAL("Wrong number of arguments for pseudo_skim_y5s_MC");
      }

	StoreObjPtr<EventMetaData> evtMetaData;
        int expNo = evtMetaData->getExperiment();
        int runNo = evtMetaData->getRun();
        int evtNo = evtMetaData->getEvent();
        std::ostringstream file_name;
        file_name << "/home/belle/yasaveev/skim_files/Y5S_skim/" << arguments[0] << "/" << arguments[1] << "/"<<expNo<<"_"<<runNo<<".txt";
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
    Manager::FunctionPtr useRecoilPiRestFrame(const std::vector<std::string>& arguments)
    {
      if (arguments.size() == 1) {
        const Variable::Manager::Var* var = Manager::Instance().getVariable(arguments[0]);
        auto func = [var](const Particle * particle) -> double {
          const Particle* p = particle->getDaughter(1);
          PCmsLabTransform T;
          ROOT::Math::PxPyPzEVector recoil = T.getBeamFourMomentum() - p->get4Vector();
          Particle pRecoil(recoil, 0);
          pRecoil.setVertex(particle->getVertex());
          UseReferenceFrame<RestFrame> frame(&pRecoil);
          double result = std::get<double>(var->function(particle));
          return result;
        };
        return func;
      } else {
        B2FATAL("Wrong number of arguments for meta function useRestFrame");
      }
    }

    VARIABLE_GROUP("CUSTOM_METAVARIABLES");
    REGISTER_METAVARIABLE("FSPProductOf(variable)", FSPProductOf,
                      "Returns product of variables.",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("FSPSumOf(variable)", FSPSumOf,
                      "Returns sum of variables.",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("varValue(listName,variable,index)", varValue,
                      "Returns variable value for the particle in listName with given index.",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("InROE(listName,index)", InROE,
                      "Returns 1 if the particle in the ROE of the list and 0 otherwise.",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("InROEMC(listName,index)", InROEMC,
                      "Returns 1 if the particle in the ROE of the list and 0 otherwise.",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("daughterHelAngle(i,j)", daughterHelAngle, "Returns helicity angle",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("minimumValueInList(particleListName, variable)", minimumValueInList, "Returns minimum value",Manager::VariableDataType::c_double); 
    REGISTER_METAVARIABLE("maximumValueInList(particleListName, variable)", maximumValueInList, "Returns maximum value",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("ChargeInLists", ChargeInLists, "Returns charge of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("PxInLists", PxInLists, "Returns the x component of momentum of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("PyInLists", PyInLists, "Returns the y component of momentum of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("PzInLists", PzInLists, "Returns the z component of momentum of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("EInLists", EInLists, "Returns E component of momentum of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("ElabInLists", ElabInLists, "Returns E component of momentum of particles in the given particle list",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("HardestGammaInROEMomentum", HardestGammaInROEMomentum, "Returns component of 4-momentum of the hardest gamma in roe",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("HardestGammaInROECMSMomentum", HardestGammaInROECMSMomentum, "Returns component of 4-momentum of the hardest gamma in roe",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("HardestGammaInMCROEMomentum", HardestGammaInMCROEMomentum, "Returns component of 4-momentum of the hardest gamma in roe",Manager::VariableDataType::c_double);
    REGISTER_METAVARIABLE("HardestGammaInMCROECMSMomentum", HardestGammaInMCROECMSMomentum, "Returns component of 4-momentum of the hardest gamma in roe",Manager::VariableDataType::c_double);
    }
  // basf2 to easily find the library and load it from the steering file
  class EnableMyMetaVariableModule: public Module {}; // And register this module to create a .map lookup file.
  REG_MODULE(EnableMyMetaVariable);
}
