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

namespace Belle2 {
  namespace Variable {

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
    VARIABLE_GROUP("CUSTOM_METAVARIABLES");
    REGISTER_METAVARIABLE("varValue(listName,variable,index)", varValue,
                      "Returns variable value for the first particle in listName.", Manager::VariableDataType::c_double);

   }
  // basf2 to easily find the library and load it from the steering file
  class EnableMyMetaVariableModule: public Module {}; // And register this module to create a .map lookup file.
  REG_MODULE(EnableMyMetaVariable);
}
