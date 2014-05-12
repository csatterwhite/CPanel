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
#include "node.h"
#include "member.h"

template<typename type>
class octree
{
    typedef std::array<double,3>    point;
    
    node<type> *root_node;
    short maxMembersPerNode;
    std::vector<member<type>> members;
    
    void boundingBox(point &boxMin, point &boxMax)
    {
        // Initialize min and max with first data point
        boxMin = members[0].getRefPoint();
        boxMax = members[0].getRefPoint();
        
        for (int i=0; i<members.size(); i++)
        {
            point pnt = members[i].getRefPoint();;
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
    
    void setDimensions(point &center, point &halfDimension)
    {
        point boxMin;
        point boxMax;
        boundingBox(boxMin,boxMax);
        for (int i=0; i<3; i++)
        {
            center[i] = 0.5*(boxMax[i]+boxMin[i]);
            halfDimension[i] = 0.5*(boxMax[i]-boxMin[i]);
        }
    }
    
    int findCorner(const member<type> &member)
    {
        // Corner   0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
        point origin = root_node->getOrigin();
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
        point refPoint = findRefPoint(obj);
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
    
    octree<type> operator=(octree<type> other)
    {
        maxMembersPerNode = other.maxMembersPerNode;
        members = other.members;
        root_node = new node<type>(*other.root_node);
        return this;
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
            point center;
            point halfDimension;
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
            point center;
            point halfDimension;
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
    
    node<type> *findNodeContainingMember(type* obj)
    {
        member<type> temp = createMember(obj);
        assert(isInsideOctree(temp));
        node<type>* current_node = root_node;
        while (!current_node->isLeafNode())
        {
            current_node = current_node->getChild(current_node->getChildContainingMember(temp));
        }
        
        return current_node;
    }
                                                                                          
    bool isInsideOctree(const member<type> &member)
    {
        point center = root_node->getOrigin();
        point halfDimension = root_node->getHalfDimension();
        for (int i=0; i<3; i++)
        {
            if (member.getRefPoint()[i] < (center[i]-halfDimension[i]) || member.getRefPoint()[i] > (center[i]+halfDimension[i]))
            {
                return false;
            }
        }
        return true;
    }
    
    virtual std::array<double,3> findRefPoint(const type* obj) = 0;
    // Returns 3 element array of X,Y,Z locations of point used to determine which node the member belongs to. i.e. (return center of triangle for unstructured grid)
    
    int getMaxMembersPerNode() {return maxMembersPerNode;}
    std::vector<member<type>> getMembers() {return members;}
    node<type> *getRootNode() {return root_node;}
};

#endif /* defined(__CPanel__octree__) */
