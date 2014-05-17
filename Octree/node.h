//
//  node.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__node__
#define __CPanel__node__

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <assert.h>
#include "member.h"


template<typename type>
class node
{
    node<type>* parent;
    node<type>* children[8];
    std::array<double,3> origin;
    std::array<double,3> halfDimension;
    std::vector<member<type>> members;
    short level;
    short maxMembers;
    
    void createChild(int childNumber)
    {
        // Child    0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
            std::array<double,3> tempOrigin = origin;
            std::array<double,3> tempHalfDimension;
            tempOrigin[0] += halfDimension[0] * (childNumber&4 ? 0.5 : -0.5);
            tempOrigin[1] += halfDimension[1] * (childNumber&2 ? 0.5 : -0.5);
            tempOrigin[2] += halfDimension[2] * (childNumber&1 ? 0.5 : -0.5);
            for (int i=0; i<3; i++)
            {
                tempHalfDimension[i] = 0.5*halfDimension[i];
            }
            children[childNumber] = new node<type>(this,tempOrigin,tempHalfDimension,level,maxMembers);
    }
    
    void pushMember(const member<type> &member)
    {
        // Checks if child is already created. If not, creates the right child before adding member.
        int child = getChildContainingMember(member);
        if (children[child] != NULL)
        {
            children[child]->addMember(member);
        }
        else
        {
            createChild(child);
            children[child]->addMember(member);
        }
    }
    
    void addChild(node<type>* child, int childNumber)
    {
        children[childNumber] = child;
    }
    
    void addLevel()
    {
        //Recursively adds one to the level of each node if parent is added
        level++;
        if (!isLeafNode())
        {
            for (int i=0; i<8; i++)
            {
                children[i]->addLevel();
            }
        }
    }
    
    std::vector<type*> membersToObjects()
    {
        std::vector<type*> objects;
        for (int i=0; i<members.size(); i++)
        {
            type* obj = members[i].getObj();
            objects.push_back(obj);
        }
        return objects;
    }
    
public:
    node(node<type>* parent_ptr,std::array<double,3> origin,std::array<double,3> halfDimension, short parent_level,  short maxMembers) : parent(parent_ptr),origin(origin),halfDimension(halfDimension),maxMembers(maxMembers)
    {
        for (int i=0; i<8; i++)
        {
            children[i] = NULL;
        }
        
        level = parent_level+1;
    }
    
    ~node()
    {
        if (!isLeafNode())
        {
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL)
                {
                    delete children[i];
                }
            }
        }
    }
    
    node(const node<type>& copy) : parent(copy.parent), origin(copy.origin), halfDimension(copy.halfDimension), members(copy.members), level(copy.level), maxMembers(copy.maxMembers)
    {
        for (int i=0; i<8; i++)
        {
            children[i] = new node<type>(*copy.children[i]);
        }
    }
    
    node<type> operator=(node<type> other)
    {
        parent = other.parent;
        origin = other.origin;
        halfDimension = other.halfDimension;
        members = other.members;
        level = other.level;
        maxMembers = other.maxMembers;
        for (int i=0; i<8; i++)
        {
            children[i] = new node<type>(*other.children[i]);
        }
        return this;
    }
    
    
    void setMaxMembers(const int &max)
    {
        maxMembers = max;
    }
    
    void createParent(int corner)
    {
        // Corner   0 1 2 3 4 5 6 7
        //   x      - - - - + + + +
        //   y      - - + + - - + +
        //   z      - + - + - + - +
        std::array<double,3> tempOrigin(origin);
        std::array<double,3> tempHalfDimension;
        tempOrigin[0] = origin[0]+halfDimension[0] * (corner&4 ? 1 : -1);
        tempOrigin[1] = origin[1]+halfDimension[1] * (corner&2 ? 1 : -1);
        tempOrigin[2] = origin[2]+halfDimension[2] * (corner&1 ? 1 : -1);
        for (int i=0; i<3; i++)
        {
            tempHalfDimension[i] = 2*halfDimension[i];
        }
        parent = new node<type>(NULL,tempOrigin,tempHalfDimension,level-1,maxMembers);
        
        int child = 0;
        for (int i=0; i<3; i++)
        {
            child = corner ^ 1<<i; //Flips the bits set for the corner with bitwise XOR
        }
        parent->addChild(this,child); //Adds current node to parent
        addLevel();
    }
    
    void addMember(const member<type> &member)
    {
        if (isLeafNode())  //If no children exist, add member to current node.  If they do, add member to proper child.
        {
            members.push_back(member);
            if (members.size() > maxMembers)
            {
                for (int i=0; i<members.size(); i++)
                {
                    pushMember(members[i]);
                }
                members.clear();
            }
        }
        else
        {
            pushMember(member);
        }
    }
    
    bool isLeafNode()
    {
        bool flag = true;
        for (int i=0; i<8; i++)
        {
            if (children[i] != NULL)
            {
                flag = false;
                break;
            }
        }
        return flag;
    }
    
    int getChildContainingMember(const member<type> &member)
    {
        int child = 0;
        if (member.getRefPoint()[0] > origin[0])
        {
            child |= 4;
        }
        if (member.getRefPoint()[1] > origin[1])
        {
            child |= 2;
        }
        if (member.getRefPoint()[2] > origin[2])
        {
            child |= 1;
        }
        
        return child;
    }
    
    std::vector<type*> getMembers()
    {
        if (isLeafNode())
        {
            return membersToObjects();
        }
        else
        {
            std::vector<type*> recursiveMembers;
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL)
                {
                    std::vector<type*> temp = children[i]->getMembers();
                    for (int j=0; j<temp.size(); j++)
                    {
                        recursiveMembers.push_back(temp[j]);
                    }
                }
            }
            return recursiveMembers;
        }
        
    }
    
    std::vector<type*> getMembers(node<type>* exception)
    {
        // Used if searching tree starting at bottom to avoid searching the same node twice.
        if (isLeafNode() && this != exception)
        {
            return membersToObjects();
        }
        else if (isLeafNode() && this == exception)
        {
            std::vector<type*> empty;
            return empty;
        }
        else
        {
            std::vector<type*> recursiveMembers;
            for (int i=0; i<8; i++)
            {
                if (children[i] != NULL && children[i] != exception)
                {
                    std::vector<type*> temp = children[i]->getMembers(exception);
                    for (int j=0; j<temp.size(); j++)
                    {
                        recursiveMembers.push_back(temp[j]);
                    }
                }
            }
            return recursiveMembers;
        }
    }

    
    std::array<double,3> getOrigin() {return origin;}
    std::array<double,3> getHalfDimension() {return halfDimension;}
    node<type>* getChild(int childNumber) {return children[childNumber];}
    short getLevel() {return level;}
    node<type>* getParent() {return parent;}
};

#endif /* defined(__CPanel__node__) */
