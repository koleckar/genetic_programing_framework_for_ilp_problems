# include "VrpRepXmlReader.h"
# include <iostream>
# include <cstring>

// todo: TinyXml bug:  Child(char*, int) doesn't work, only Child(int)

VrpRepXmlReader::VrpRepXmlReader(const char *filePath) : xmlFile(filePath), xmlFileHandle(&xmlFile) {
    if (xmlFile.LoadFile() == false) {
        throw std::exception("XML file not loaded. File path: " + *filePath);
    }

    getDatasetName();
    std::cout << "Dataset \"" << datasetName << "\" loaded." << std::endl;
}

void VrpRepXmlReader::getDatasetName() {
    TiXmlElement *datasetElem = xmlFileHandle
            .FirstChild("instance")
            .FirstChild("info")
            .FirstChild("dataset").ToElement();

    if (datasetElem) {
        datasetName = datasetElem->FirstChild()->Value();
    }
}

int VrpRepXmlReader::getNumberOfNodes() {
    TiXmlElement *node = xmlFileHandle
            .FirstChild("instance")
            .ChildElement(1)
            .FirstChild("nodes")
            .FirstChild().ToElement();

    int numberOfNodes = 0;
    for (; node; node = node->NextSiblingElement()) {
        numberOfNodes++;
    }

    return numberOfNodes;
}

// todo: not using the read node type ( node 0 = depot, ...)
std::vector<std::pair<float, float>> VrpRepXmlReader::getNodesCoordinates(int numberOfNodes) {
    std::vector<std::pair<float, float>> nodesCoordinates;
    nodesCoordinates.reserve(numberOfNodes);

    TiXmlElement *node = xmlFileHandle
            .FirstChild("instance")
            .ChildElement(1)
            .FirstChild("nodes")
            .FirstChild().ToElement();

    char *p;
    int ID, nodeType;
    float cx, cy;
    for (node; node; node = node->NextSiblingElement()) {
        ID = std::strtol(node->FirstAttribute()->Value(), &p, 10);
        nodeType = std::strtol(node->LastAttribute()->Value(), &p, 10);
        cx = std::strtof(node->FirstChild()->FirstChild()->Value(), &p);
        cy = std::strtof(node->LastChild()->FirstChild()->Value(), &p);
        nodesCoordinates.emplace_back(cx, cy);
    }

    return nodesCoordinates;
}

int VrpRepXmlReader::getFleetSize() {

    const char *fleetSizeSpecification = xmlFileHandle
            .FirstChild("instance")
            .Child(2) // <fleet>
            .FirstChild("vehicle_profile").ToElement()->Attribute("number");

    char *p;
    return std::strtol(fleetSizeSpecification, &p, 10);
}

std::vector<std::pair<float, float>> VrpRepXmlReader::getTimeWindows(int numberOfNodes) {
    std::vector<std::pair<float, float>> timeWindows;
    timeWindows.reserve(numberOfNodes);
    timeWindows.emplace_back(-1, -1);// depot doesn't have time window, invalid for depot
    TiXmlElement *request = xmlFileHandle
            .FirstChild("instance")
            .Child(3)
            .FirstChild("request").ToElement();

    char *p;
    int ID;
    float start, end;
    for (request; request; request = request->NextSiblingElement()) {
        ID = std::strtol(request->FirstAttribute()->Value(), &p, 10);

        TiXmlElement *twElem = request->FirstChildElement();
        TiXmlNode *startNode = twElem->FirstChildElement()->FirstChild();
        TiXmlNode *endNode = twElem->FirstChildElement()->NextSiblingElement()->FirstChild();

        start = std::strtof(startNode->Value(), &p);
        end = std::strtof(endNode->Value(), &p);
        timeWindows.emplace_back(start, end);
    }

    return timeWindows;
}

std::vector<float> VrpRepXmlReader::getServiceTimes(int numberOfNodes) {
    std::vector<float> serviceTimes;
    serviceTimes.reserve(numberOfNodes);
    serviceTimes.push_back(-1); // invalid for depot.

    TiXmlElement *request = xmlFileHandle
            .FirstChild("instance")
            .Child(3)
            .FirstChild("request").ToElement();

    char *p;
    int ID;
    float service_time;
    for (request; request; request = request->NextSiblingElement()) {

        // todo: really tinyxml doesnt work by tag names? only indices?
        TiXmlNode *service_timeNode
                = request->FirstChildElement()->NextSiblingElement()->NextSiblingElement()->FirstChild();

        service_time = std::strtof(service_timeNode->Value(), &p);
        serviceTimes.push_back(service_time);
    }

    return serviceTimes;
}

std::vector<float> VrpRepXmlReader::getNodeDemands(int numberOfNodes) {
    std::vector<float> demands;
    demands.reserve(numberOfNodes);
    demands.push_back(-1); // invalid for depot.

    TiXmlElement *request = xmlFileHandle
            .FirstChild("instance")
            .Child(3) // <requests>
            .FirstChild("request").ToElement();

    char *p;
    int ID;
    float demand;
    for (request; request; request = request->NextSiblingElement()) {
        // <tw>, <quantity>, <service_time>
        // <quantity>
        //------- TODO ------------
        TiXmlNode *service_timeNode
                = request->FirstChildElement()->NextSiblingElement()->FirstChild();
              //= request->FirstChildElement()->FirstChild();

        demand = std::strtof(service_timeNode->Value(), &p);
        demands.push_back(demand);
    }

    return demands;
}

float VrpRepXmlReader::getVehicleCapacity() {
    const char *vehicleCapacity = xmlFileHandle
            .FirstChild("instance")
            .Child(2) // <fleet>
            .FirstChild("vehicle_profile")
            .Child(2) // <capacity>
            .ToElement()->FirstChild()->Value();

    char *p;
    return std::strtof(vehicleCapacity, &p);
}
