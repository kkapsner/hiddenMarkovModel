#include "HMMConfiguration.h"
#include <iostream>
#ifndef JSON_IS_AMALGAMATION
# define JSON_IS_AMALGAMATION
#endif
#include "include/json/json.h"

#define READ_CONFIG(name, as){\
	configuration.name = root.get(#name, configuration.name).as();\
}

using namespace hiddenMarkovModel;

HMMConfiguration::HMMConfiguration():
	verbose(false),
	verboseOutputTransition(true),
	verboseOutputEmission(false),
	
	pauseAfterIteration(false),

	minSelfTransition(0),
	minEmission(1e-6),

	doEmissionUpdate(true),
	doTransitionUpdate(true),
	binningCount(300),
	maxIterations(100),
	abortStateChanges(5)
	{};

HMMConfiguration HMMConfiguration::fromFile(std::istream &file){
	Json::Value root;
	Json::Reader reader;

	HMMConfiguration configuration;
	if (reader.parse(file, root)){
		if (root.isMember("verbose") && root["verbose"].isObject()){
			configuration.verbose = root["verbose"].get("enabled", configuration.verbose).asBool();
			configuration.verboseOutputTransition = root["verbose"].get("outputTransition", configuration.verboseOutputTransition).asBool();
			configuration.verboseOutputEmission = root["verbose"].get("outputEmission", configuration.verboseOutputEmission).asBool();
		}
		READ_CONFIG(pauseAfterIteration, asBool);
		READ_CONFIG(minSelfTransition, asDouble);
		READ_CONFIG(minEmission, asDouble);
		READ_CONFIG(doEmissionUpdate, asBool);
		READ_CONFIG(doTransitionUpdate, asBool);
		READ_CONFIG(binningCount, asUInt);
		READ_CONFIG(maxIterations, asUInt);
		READ_CONFIG(abortStateChanges, asUInt);
	}
	else {
		std::cerr << "Failed to parse configuration" << std::endl
			<< reader.getFormattedErrorMessages();
	}

	return configuration;
}