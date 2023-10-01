#include"Master.hpp"
#include<vector>
#include<cassert>

FiniteDifferenceGrid::FiniteDifferenceGrid(std::vector<Node>::size_type numNodes, double xMin, double xMax)
{
    int znumNodes = static_cast<int>(numNodes);
    double stepsize= (xMax-xMin)/((double) numNodes-1);
    for (int i=0; i < znumNodes; i++)
    {
        Node node;
        node.coordinate = xMin+i*stepsize;
        mNodes.push_back(node);
    }
    assert(mNodes.size() == numNodes);
}