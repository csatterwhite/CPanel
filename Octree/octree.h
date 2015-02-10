//
//  octree.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__octree__
#define __CPanel__octree__

#include <iostream>
#include <Eigen/Dense>
#include "node.h"
#include "member.h"

template<typename type>
class octree
{
    node<type> *root_node;
    short maxMembersPerNode;
    std::vector<member<type>> members;
    
    void boundingBox(Eigen::Vector3d &boxMin, Eigen::Vector3d &boxMax)
    {
        // Initialize min and max with first data point
        boxMin = members[0].getRefPoint();
        boxMax = members[0].getRefPoint();
        
        for (int i=0; i<members.size(); i++)
        {
            Eigen::Vector3d pnt = members[i].getRefPoint();;
            for (int j=0; j<3; j++)
            {
                if (pnt[j]<boxMin[j])
                {
                    boxMin[j] = pnt[j];
                }
                if (pnt[j]>boxMax[j])
                {
                    boxMax[j] = pnt[j];
                }
            }
        }

    }
    
    void setDimensions(Eigen::Vector3d &center, Eigen::Vector3d &halfDimension)
    {
        halfDimension << 1,1,1;
        Eigen::Vector3d boxMin;
        Eigen::Vector3d boxMax;
        boundingBox(boxMin,boxMax);
        
        center = 0.5*(boxMin+boxMax);
        double maxD = (boxMax-boxMin).maxCoeff();
        halfDimension *= 0.51*maxD; // Added 1% to half dimension to handle unique cases where floating point error leaves a node outside the octree.
    }
    
    int findCorner(const member<type> &member)
    {
        // Corner   0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
        Eigen::Vector3d origin = root_node->getOrigin();
        int corner = 0;
        if (member.getRefPoint()[0] > origin[0])
        {
            corner |= 4;
        }
        if (member.getRefPoint()[1] > origin[1])
        {
            corner |= 2;
        }
        if (member.getRefPoint()[2] > origin[2])
        {
            corner |= 1;
        }
        return corner;
    }
    
    member<type> createMember(type* obj)
    {
        Eigen::Vector3d refPoint = findRefPoint(*obj);
        member<type> newMember(obj,refPoint);
        return newMember;
    }
    
public:
    
    octree() : maxMembersPerNode(10), root_node(NULL) {}
    
    ~octree()
    {
        delete root_node;
    }
    
    octree(const octree& copy) : maxMembersPerNode(copy.maxMembersPerNode), members(copy.members)
    {
        root_node = new node<type>(*copy.root_node);
    }
    
    octree<type>& operator=(const octree<type> &rhs)
    {
        if (this == &rhs)
        {
            return (*this);
        }
        maxMembersPerNode = rhs.maxMembersPerNode;
        members = rhs.members;
        root_node = new node<type>(*rhs.root_node);
        return *this;
    }
    
    void setMaxMembers(const int &maxMembers)
    {
        maxMembersPerNode = maxMembers;
    }

    
    void addData(const std::vector<type*> &newData)
    {
        size_t iter = 0;
        if (root_node != NULL)
        {
            iter = members.size();
        }
        
        for (int i=0; i<newData.size(); i++)
        {
            members.push_back(createMember(newData[i]));
        }
        
        if (root_node == NULL)  // If octree hasn't been created, add data and create octree containing data
        {
            Eigen::Vector3d center;
            Eigen::Vector3d halfDimension;
            setDimensions(center,halfDimension);
            
            root_node = new node<type>(NULL,center,halfDimension,0,maxMembersPerNode);
            
            for (int i=0; i<members.size(); i++)
            {
                root_node->addMember(members[i]);
            }
        }
        else  // If octree exists, check points to see if they are inside octree.  If not, expand octree to contain point by giving root_node a parent.
        {
            for (size_t i=iter; i<members.size(); i++)
            {
                if (!isInsideOctree(members[i]))
                {
                    root_node->createParent(findCorner(members[i]));
                    root_node = root_node->getParent(); //Update root_node to be parent just created

                    root_node->addMember(members[i]);
                }
                else
                {
                    root_node->addMember(members[i]);
                }
            }
        }
    }
    
    void addData(type* newData)
    {
        member<type> newMember = createMember(newData);
        members.push_back(newMember);

        if (root_node == NULL)  // If octree hasn't been created, add data and create octree containing data
        {
            Eigen::Vector3d center;
            Eigen::Vector3d halfDimension;
            setDimensions(center,halfDimension);
            
            root_node = new node<type>(NULL,center,halfDimension,0,maxMembersPerNode);
            
            root_node->addMember(newMember);
        }
        else  // If octree exists, check point to see if it is inside octree.  If not, expand octree to contain point by giving root_node a parent.
        {
            if (!isInsideOctree(newMember))
            {
                root_node->createParent(findCorner(newMember));
                root_node = root_node->getParent(); //Update root_node to be parent just created
                
                root_node->addMember(newMember);
            }
            else
            {
                root_node->addMember(newMember);
            }
        }
    }
    
    node<type>* findNodeContainingPnt(const Eigen::Vector3d pnt)
    {
        node<type>* current_node = root_node;
        while (!current_node->isLeafNode())
        {
            current_node = current_node->getChild(current_node->getChildContainingPnt(pnt));
        }
        
        return current_node;
    }
    
    node<type>* findNodeContainingMember(type* obj)
    {
        member<type> temp = createMember(obj);
        Eigen::Vector3d pnt = temp.getRefPoint();
        return findNodeContainingPnt(pnt);
    }
                                                                                          
    bool isInsideOctree(const member<type> &member)
    {
        return isInsideOctree(member.getRefPoint());
    }
    
    bool isInsideOctree(const Eigen::Vector3d &pnt)
    {
        Eigen::Vector3d center = root_node->getOrigin();
        Eigen::Vector3d halfDimension = root_node->getHalfDimension();
        for (int i=0; i<3; i++)
        {
            if (pnt(i) < (center[i]-halfDimension[i]) || pnt(i) > (center[i]+halfDimension[i]))
            {
                return false;
            }
        }
        return true;
    }
    
    std::vector<node<type>*> getNodes()
    {
        return root_node->getSubNodes();
    }
    
    virtual Eigen::Vector3d findRefPoint(const type &obj) = 0;
    // Returns 3 element array of X,Y,Z locations of point used to determine which node the member belongs to. i.e. (return center of triangle for unstructured grid)
    
    int getMaxMembersPerNode() {return maxMembersPerNode;}
    std::vector<member<type>> getMembers() {return members;}
    node<type> *getRootNode() {return root_node;}
};

#endif /* defined(__CPanel__octree__) */
