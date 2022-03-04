#include "../../lib/tinyxml/tinyxml.h"
# include <vector>
# include <string>

#ifndef SEMESTRAL_WORK_DATASETXMLREADER_H
#define SEMESTRAL_WORK_DATASETXMLREADER_H

/**
 Class facilitating of reading the .xml datasets in "VRP-REP" format.
 */
class VrpRepXmlReader {
    TiXmlDocument xmlFile;
    TiXmlHandle xmlFileHandle;

public:
    std::string datasetName;

    VrpRepXmlReader(const char *filePath);

    void getDatasetName();

    int getNumberOfNodes();

    int getFleetSize();

    std::vector<std::pair<float, float>> getNodesCoordinates(int numberOfNodes);

    std::vector<std::pair<float, float>> getTimeWindows(int numberOfNodes);

    std::vector<float> getServiceTimes(int numberOfNodes);

    std::vector<float> getNodeDemands(int numberOfNodes);

    float getVehicleCapacity();
};


#endif //SEMESTRAL_WORK_DATASETXMLREADER_H
