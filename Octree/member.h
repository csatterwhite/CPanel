//
//  member.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef CPanel_member_h
#define CPanel_member_h

template<typename type>
class member
{
    type* pObj;
    Eigen::Vector3d refPoint;
    
public:
    member(type* object, Eigen::Vector3d referencePoint) : pObj(object), refPoint(referencePoint) {}
    
    Eigen::Vector3d getRefPoint() const {return refPoint;}
    
    type *getObj() {return pObj;}
};


#endif
